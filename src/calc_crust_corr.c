//#include "iscloc.h"
#include "GridLib.h"
#include "calc_crust_corr.h"
#include "crust_corr_model.h"
#include "crust_type.h"
#include "crust_type_key.h"

#define KM_PER_DEG (10000.0/90.0)

/*#DOC  Title:																*/
/*#DOC    calc_crust_corr													*/
/*#DOC  Desc:																*/
/*#DOC    Returns travel time correction for a source or station.			*/
/*#DOC  Input Arguments:													*/
/*#DOC    ps    - either 'P' or 'S'.										*/
/*#DOC    lat - lattitude of either source or station.						*/
/*#DOC    lon - longitude of either source or station.						*/
/*#DOC    depth - depth of source.  0 for station.							*/
/*#DOC    elev  - elevation for station. NULLVAL for source.				*/
/*#DOC    dtdd  - dt/dd for phase.											*/
/*#DOC  Return:																*/
/*#DOC    Travel time correction.											*/

/* Isostatic correction assume total density*height constant - this will	*/
/* include an average elevation for the tile.								*/

/* At each end:																*/
/* Real tt = tt + crust time + time in extra mantle - JB crust time.		*/
/* Source more complicated if it's half way down crust.						*/

/* When considering effects of take off angle use p = sinQ/vn = dtdd		*/
/* where vn is velocity at top of mantle and dtdd in s/km					*/
/* in combination with  Snells law to get:									*/
/* t = h/v * 1 / ( 1 - (v * dtdd)2 )1/2										*/

/* Ignore water and ice layers for travel times but include in isostasy.	*/

double calc_crust_corr (char ps, double lat, double lon, double depth, double elev, double dtdd )
{
	double g_vel,b_vel,n_vel,vel[8];
	double crust_time,jb_crust_time;
	double layer_time,cumulative_thick;
	double iso_height,extra_mantle,iso_corr;
	double uplift = 0.0, elev_diff = 0.0, elev_corr = 0.0;
	double total_corr;
	int c,col,row;
	int i, diagnostic;

	/* Switch debug messages from this function on/off (1/0). */
	diagnostic = message_flag >= 5;
//diagnostic = 1;

	/* Look up crust-type number 2x2 deg tile. */
	col =  (int)((90 - lat)/2);
	row =  (int)((180 + lon)/2);
	c = crust_type[col][row];

	if (ps == 'P'){
		g_vel = PGVEL;
		b_vel = PBVEL;
		n_vel = PNVEL;
		for(i=0;i<8;i++)
			vel[i]=c_type[c].p_vel[i];
	}
	else if (ps == 'S'){
		g_vel = SGVEL;
		b_vel = SBVEL;
		n_vel = SNVEL;
		for(i=0;i<8;i++)
			vel[i]=c_type[c].s_vel[i];
	}
	else {
		sprintf(MsgStr,"calc_crust_corr: wrong ps value %c",ps);
		nll_putmsg(1, MsgStr);
		return 0.0;
	}

	/* Only for sources above Moho (both Mohos). */
	if (depth > MOHO || depth > c_type[c].thick[8])
		return 0.0;

	/* Isostatic correction. Include water/ice. */
	iso_height=0;
	for (i=0; i<7; i++)
		iso_height += c_type[c].thick[i]*c_type[c].density[i];

	/* Difference in height*density is upper mantle height*density. */
	extra_mantle = (JB_ISO_HEIGHT - iso_height)/BELOW_RHO;

	/* Isostatic correction is time taken to traverse extra mantle. */
	iso_corr = n_vel*(dtdd/KM_PER_DEG);
	iso_corr *= iso_corr;
	iso_corr  = 1.0 / sqrt( 1.0 - iso_corr);
	iso_corr *= extra_mantle/n_vel;


	/* Time to traverse J-B crust. */
	if (depth < CONRAD){

		/* Whole of lower crust. */
		layer_time = b_vel*(dtdd/KM_PER_DEG);
		layer_time *= layer_time;
		layer_time = 1.0 / sqrt( 1.0 - layer_time);
		layer_time *= (MOHO-CONRAD)/b_vel;

		jb_crust_time = layer_time;

		/* Part of upper crust. */
		layer_time = g_vel*(dtdd/KM_PER_DEG);
		layer_time *= layer_time;
		layer_time = 1.0 / sqrt( 1.0 - layer_time);
		layer_time *= (CONRAD-depth)/g_vel;

		jb_crust_time += layer_time;
	}
	else{
		/* Part of lower crust. */
		layer_time = b_vel*(dtdd/KM_PER_DEG);
		layer_time *= layer_time;
		layer_time = 1.0 / sqrt( 1.0 - layer_time);
		layer_time *= (MOHO-depth)/b_vel;

		jb_crust_time = layer_time;
	}

	/* Time to traverse model crust. */
	crust_time = 0;
	cumulative_thick=0;
	for (i=2; i<7; i++){
		cumulative_thick+=c_type[c].thick[i];
		if (depth > cumulative_thick)
			continue;
		/* Whole layer. */
		if (crust_time){

			layer_time  = vel[i]*(dtdd/KM_PER_DEG);
			layer_time *= layer_time;
			if (layer_time >= 1){  /* Goes through layer horizontaly. */
				sprintf(MsgStr,"WARNING: calc_crust_corr: layer_time = %f",layer_time);
				nll_putmsg(1, MsgStr);
				continue;
			}
			layer_time  = 1.0 / sqrt( 1.0 - layer_time);
			layer_time *= c_type[c].thick[i]/vel[i];

			crust_time += layer_time;
		}
		/* Part of layer the source is in. */
		else {
			layer_time  = vel[i]*(dtdd/KM_PER_DEG);
			layer_time *= layer_time;
			if (layer_time >= 1){  /* Goes through layer horizontaly. */
				sprintf(MsgStr,"WARNING: calc_crust_corr: layer_time = %f",layer_time);
				continue;
			}
			layer_time  = 1.0 / sqrt( 1.0 - layer_time);
			layer_time *= cumulative_thick - depth;
			layer_time /= vel[i];

			crust_time += layer_time;
		}
	}


	/* Station elevation correction.					*/
	/* thick[8] from crust_type is total thickness.		*/
	if (depth == 0 && elev < VERY_LARGE_DOUBLE){
		uplift = c_type[c].thick[8] + extra_mantle - MOHO;
		elev_diff = elev/1000 - uplift;
		if (c_type[c].thick[2] + elev_diff < 0){
			sprintf(MsgStr,"WARNING: calc_crust_corr: elev_diff=%f",elev_diff);
		}
		else {
			layer_time  = vel[2]*(dtdd/KM_PER_DEG);
			layer_time *= layer_time;
			layer_time  = 1.0 / sqrt( 1.0 - layer_time);
			layer_time *= elev_diff/vel[2];

			elev_corr = layer_time;
		}
	}
	else
		elev_corr = 0;


	total_corr = crust_time + iso_corr - jb_crust_time + elev_corr;

	if (diagnostic /*|| fabs(total_corr) > 2.0*/){
		sprintf(MsgStr,"\ncalc_crust_corr: lat=%.3f lon=%.3f depth=%.3f elev=%.3f", lat,lon, depth, elev);
		nll_putmsg(1, MsgStr);
		sprintf(MsgStr,"c[%d][%d]=%d",col,row,c);
		nll_putmsg(1, MsgStr);
		sprintf(MsgStr,"calc_crust_corr: extra_mantle=%.3f ",extra_mantle);
		nll_putmsg(1, MsgStr);
		sprintf(MsgStr,"iso_corr=%.3f",iso_corr);
		nll_putmsg(1, MsgStr);
		sprintf(MsgStr,"calc_crust_corr:dtdd=%.3f ",dtdd);
		nll_putmsg(1, MsgStr);
		sprintf(MsgStr,"time %.3f cf jb %.3f",crust_time,jb_crust_time);
		nll_putmsg(1, MsgStr);
		sprintf(MsgStr,"calc_crust_corr: uplift=%.3f ",uplift);
		nll_putmsg(1, MsgStr);
		sprintf(MsgStr,"elev_diff=%.3f ",elev_diff);
		nll_putmsg(1, MsgStr);
		sprintf(MsgStr,"elev_corr=%.3f",elev_corr);
		nll_putmsg(1, MsgStr);
		sprintf(MsgStr,"calc_crust_corr: lat=%.3f lon=%.3f depth=%.3f elev=%.3f", lat,lon, depth, elev);
		nll_putmsg(1, MsgStr);
		sprintf(MsgStr,"calc_crust_corr: tot_corr=%.3f \n", total_corr);
		nll_putmsg(1, MsgStr);
	}

	return (total_corr);
}



