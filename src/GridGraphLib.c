/*
 * Copyright (C) 1999 Anthony Lomax <lomax@faille.unice.fr>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */


/* Grid graphics functions

	07OCT1997	AJL	modified from wls graphics functions

*/

#define EXTERN_MODE 1
#define INPUT_DEV 1
#include "GridLib.h"
#include "GridGraphLib.h"
#include "vector/vector.h"



/*** function to set constants */

void SetGraphicConstants(void)
{
	gr_num_sta = 0;			/* number of stations to draw */
	gr_num_txt = 0;			/* number of text strings to draw */
	gr_num_line = 0;		/* number of lines to draw */
	gr_num_map_files = 0;		/* number of map lines to draw */
}


/* function to read input file */

int ReadGraphicsInput(FILE* fp_input, char* param, char* line, int istat)
{


	istat = -1;

	/* read parameter line */


	/* skip comment line or white space */

	if (strncmp(param, "#", 1) == 0 || isspace(param[0]))
		istat = 0;


#if 0
	/* read station params */
	if (strncmp(param, "STATION", 7) == 0)
		if ((istat = get_stations(strchr(line, ' ') + 1, fp_input)) < 0)
			fprintf(stderr, "Error reading station parameters.\n");

	/* read graph text params */
	if (strncmp(param, "GTEXT", 5) == 0)
		if ((istat = get_gtext(strchr(line, ' ') + 1)) < 0)
			fprintf(stderr, "Error reading graphics text parameters.\n");

	/* read graph line params */
	if (strncmp(param, "GLINE", 5) == 0)
		if ((istat = get_gline(strchr(line, ' ') + 1, fp_input)) < 0)
			fprintf(stderr, "Error reading graphics line parameters.\n");

#endif

	/* read mapline params */
	if (strncmp(param, "MAPLINE", 7) == 0)
		if ((istat = get_maplines(strchr(line, ' ') + 1)) < 0)
			fprintf(stderr, "Error reading mapline parameters.\n");

	return(istat);


}




/*** function to read mapline parameters from input line ***/

int get_maplines(char* input_line)
{

		/* add map file to mapfile array */

	if (gr_num_map_files == MAX_NUM_MAP) {
		fprintf(stderr, "Error - Too many map line files.\n");
		return (-1);
	}

		/* read mapline input line */

	sscanf(input_line, "%s %s %lf %lf %lf %s",
		mapfile[gr_num_map_files].format,
		mapfile[gr_num_map_files].name,
		&(mapfile[gr_num_map_files].rgb.r),
		&(mapfile[gr_num_map_files].rgb.g),
		&(mapfile[gr_num_map_files].rgb.b),
		mapfile[gr_num_map_files].line_style);

	sprintf(MsgStr, "MAPLINE  fmt %s  `%s' ",
      		mapfile[gr_num_map_files].format,
		mapfile[gr_num_map_files].name);
 	nll_putmsg(2, MsgStr);
	sprintf(MsgStr, "  col r %f g %f b %f   line_style %s",
       		mapfile[gr_num_map_files].rgb.r,
		mapfile[gr_num_map_files].rgb.g,
		mapfile[gr_num_map_files].rgb.b,
		mapfile[gr_num_map_files].line_style);
 	nll_putmsg(2, MsgStr);

	gr_num_map_files++;

    return (0);
}




/*** function to draw map lines */

#if 0 /*!!!!!!!!!!*/
#define MAXLINEPTS 10000
int PlotMapLines()
{
	int  nmapfile, npts;
	char chardummy;
	double maplong, maplat, xtemp, ytemp, xloc, yloc;
	double c111, angle, cosang, sinang;
	double xpt[MAXLINEPTS], ypt[MAXLINEPTS];
	FILE *fp_map;

	if (!prog_mode_3d)
		return;

	c111 = 10000.0 / 90.0;

	      /* draw map lines */

 	set_drawing_window(norm_trans_map);

	for (nmapfile = 0; nmapfile < gr_num_map_files; nmapfile++) {

    		set_line_color(mapfile[nmapfile].rgb, );
		set_line_style(mapfile[nmapfile].line_style);

			/* open input file */

		if ((fp_map = fopen(mapfile[nmapfile].name, "r")) == NULL) {
			fprintf(stderr, "Error - Cannot open map file `%s'.\n",
				mapfile[nmapfile].name);
    	}

		angle = -cRPD * mapfile[nmapfile].map_rot;
		cosang = cos(angle);
		sinang = sin(angle);

			/* draw each segment */

    	while (fgets(line, MAXLINE, fp_map) != NULL) {

			npts = 0;
			for ( ; ; ) {
				fscanf(fp_map, "%lf %lf \n", &maplong, &maplat);
				if (maplong > 200.0 || maplong < -200.0) {
					if (npts > 1)
						drawpolyline(npts, xpt, ypt);
					break;
				}
				if (npts < MAXLINEPTS) {
					xtemp = (maplong - mapfile[nmapfile].map_orig_long)
						* c111 * cos(cRPD * maplat);
					ytemp = (maplat - mapfile[nmapfile].map_orig_lat) * c111;
					xpt[npts] = xtemp * cosang - ytemp * sinang;
					ypt[npts++] = ytemp * cosang + xtemp * sinang;
				}
			}
		}
	}

	set_drawing_window(norm_trans_all);
    setline(SOLID);

}
#endif  /*!!!!!!!!!!*/



/*** function to draw stations on model */

#if 0  /*!!!!!!!!!!*/
int Stations2GMT()
{
	int  nsta;
	int slen;
	static double xpos[8] = {-1.0, -1.0, 1.0, 1.0, -1.0, 0.0, 1.0, 0.0};
	static double ypos[8] = {-1.0, 1.0, 1.0, -1.0, 0.0, 1.5, 0.0, -1.5};
	static double len[8] = {-1.0, -1.0, 0.0, 0.0, -1.0, -0.5, 0.0, -0.5};
	double side;

	/* draw and label each station */

	side = (xmodmax - xmodmin) / 80.0;

	for (nsta = 0; nsta < gr_num_sta; nsta++) {
		set_text_color(sta_array[nsta].col);
		reset_char_hgt(sta_array[nsta].txt_size);
		drawtriangle(sta_array[nsta].xloc, sta_array[nsta].zloc + 1.0,
			side, 0, get_col(sta_array[nsta].col), 1,
			get_col(sta_array[nsta].col));
		if (prog_mode_3d) {
			set_drawing_window(norm_trans_map);
			drawtriangle(sta_array[nsta].xloc, sta_array[nsta].yloc,
				side, 0, get_col(sta_array[nsta].col), 1,
				get_col(sta_array[nsta].col));
				gtextw(sta_array[nsta].xloc - chrwthw / 4.0 +
					xpos[sta_array[nsta].txt_pos] * chrwthw +
					len[sta_array[nsta].txt_pos] * chrwthw *
				(double) (slen = strlen(sta_array[nsta].name)),
					sta_array[nsta].yloc - chrhtw / 2.0 +
					ypos[sta_array[nsta].txt_pos] * chrhtw,
					sta_array[nsta].name);
			set_drawing_window(norm_trans_all);
		}
		reset_char_hgt(1.0);
	}
}
#endif  /*!!!!!!!!!!*/


/*** function to draw user text on model */

#if 0  /*!!!!!!!!!!*/
int Text2GMT()
{
	int  ntxt;

	for (ntxt = 0; ntxt < gr_num_txt; ntxt++) {
		set_text_color(txt_array[ntxt].col);
		reset_char_hgt(txt_array[ntxt].size);
		gcharup(-sin(cRPD * txt_array[ntxt].rot),
			cos(cRPD * txt_array[ntxt].rot));
		gchar_rot(txt_array[ntxt].rot);
		if (prog_mode_3d && strncmp(txt_array[ntxt].plot, "MAP", 3) == 0) {
			set_drawing_window(norm_trans_map);
			gtextw(txt_array[ntxt].xloc, txt_array[ntxt].yloc,
				txt_array[ntxt].text);
			set_drawing_window(norm_trans_all);
		} else if (strncmp(txt_array[ntxt].plot, "SECT", 4) == 0) {
			gtextw(txt_array[ntxt].xloc, -txt_array[ntxt].zloc,
				txt_array[ntxt].text);
		} else if (strncmp(txt_array[ntxt].plot, "TIME", 4) == 0) {
			set_clip(0);
			gtextw(txt_array[ntxt].xloc,
				ttowind(0.0, (double) txt_array[ntxt].yloc),
				txt_array[ntxt].text);
			set_clip(1);
		}
		gcharup(0.0, 1.0);
		gchar_rot(0.0);
		reset_char_hgt(1.0);
	}

}
#endif  /*!!!!!!!!!!*/


/*** function to draw user lines on model */

#if 0  /*!!!!!!!!!!*/
int Line2GMT()
{
	int  nline;

	for (nline = 0; nline < gr_num_line; nline++) {
		set_line_color(line_array[nline].col);
		set_line_style(line_array[nline].style);
		if (prog_mode_3d && strncmp(line_array[nline].plot, "MAP", 3) == 0) {
			set_drawing_window(norm_trans_map);
   			 drawline(line_array[nline].x1, line_array[nline].y1,
				line_array[nline].x2, line_array[nline].y2);
			set_drawing_window(norm_trans_all);
		} else if (strncmp(line_array[nline].plot, "SECT", 4) == 0) {
   			 drawline(line_array[nline].x1, line_array[nline].z1,
				line_array[nline].x2, line_array[nline].z2);
		}
	}

	setline(SOLID);
}

#endif  /*!!!!!!!!!!*/



/*** function to get current date and time in character format ***/

char *get_timestamp()
{
	static char timestr[21];
	time_t tp;
	struct tm *ltime;

	(void) time(&tp);

	ltime = localtime(&tp);

	strftime(timestr, (size_t) 20, "%d %b %y  %H:%M", ltime);

	return(timestr);

}



/** function to write GMT psxy data to draw an error cross */

int Err2GMT(FILE* fp_io, double xcent, double ycent,
	double xerrlen, double xbarlen, double yerrlen, double ybarlen)
{
	double xloc, yloc;


	fprintf(fp_io, ">\n");

	/* cross */

	xloc = xcent;
	yloc = ycent - yerrlen / 2.0;
	fprintf(fp_io, "%lf %lf\n", xloc, yloc);
	yloc = ycent + yerrlen / 2.0;
	fprintf(fp_io, "%lf %lf\n>\n", xloc, yloc);

	xloc = xcent - xerrlen / 2.0;
	yloc = ycent;
	fprintf(fp_io, "%lf %lf\n", xloc, yloc);
	xloc = xcent + xerrlen / 2.0;
	fprintf(fp_io, "%lf %lf\n>\n", xloc, yloc);

	/* bars */

	xloc = xcent - xerrlen / 2.0;
	yloc = ycent - ybarlen / 2.0;
	fprintf(fp_io, "%lf %lf\n", xloc, yloc);
	yloc = ycent + ybarlen / 2.0;
	fprintf(fp_io, "%lf %lf\n>\n", xloc, yloc);

	xloc = xcent + xerrlen / 2.0;
	yloc = ycent - ybarlen / 2.0;
	fprintf(fp_io, "%lf %lf\n", xloc, yloc);
	yloc = ycent + ybarlen / 2.0;
	fprintf(fp_io, "%lf %lf\n>\n", xloc, yloc);

	xloc = xcent - xbarlen / 2.0;
	yloc = ycent + yerrlen / 2.0;
	fprintf(fp_io, "%lf %lf\n", xloc, yloc);
	xloc = xcent + xbarlen / 2.0;
	fprintf(fp_io, "%lf %lf\n>\n", xloc, yloc);

	xloc = xcent - xbarlen / 2.0;
	yloc = ycent - yerrlen / 2.0;
	fprintf(fp_io, "%lf %lf\n", xloc, yloc);
	xloc = xcent + xbarlen / 2.0;
	fprintf(fp_io, "%lf %lf\n>\n", xloc, yloc);

	return(0);

}



/** function to convert covariance matrix to an array of points defining  an
		error ellipse (covariance ellipse) */

/* from "Bivariate Normal Distribution Web page
	http://www.cern.ch/Physics/DataAnalysis/BriefBook/AN10pp/node17.html */

Vect2D* Cov2Ellipse(double covxx,  double covyy, double covxy,
		double xcent, double ycent, Vect2D* xyarray, int npts)
{
	int n;

	double sigma1_2, sigma1, sigma2_2, sigma2, rho, rho_2;
	double angle, cosa, cosa_2, sina, sina_2, p1, p2;
	double theta, dtheta, xloc, yloc;


	/* sigma**2, sigma */
	sigma1_2 = covxx;
	sigma1 = sqrt(sigma1_2);
	sigma2_2 = covyy;
	sigma2 = sqrt(sigma2_2);

	/* corellation coefficient */
	rho = covxy / (sigma1 * sigma2);
	rho_2 = rho * rho;

	/* angle between x-axis and p1 */
	angle = 0.5 * atan2(2.0 * rho * sigma1 * sigma2, sigma1_2 - sigma2_2);
	cosa = cos(angle);
	cosa_2 = cosa * cosa;
	sina = sin(angle);
	sina_2 = sina * sina;

	/* semi-diameters */
	p1 = sqrt( sigma1_2 * sigma2_2 * (1.0 - rho_2) /
		(sigma2_2 * cosa_2 - 2.0 * rho * sigma1 * sigma1 * sina * cosa
			+ sigma1_2 * sina_2) );
	p2 = sqrt( sigma1_2 * sigma2_2 * (1.0 - rho_2) /
		(sigma2_2 * sina_2 + 2.0 * rho * sigma1 * sigma1 * sina * cosa
			+ sigma1_2 * cosa_2) );

	/* prob inside ellipse is ~39%, multiply by 2 to get ~??% */

	p1 *= 2.0;
	p2 *= 2.0;

	/* generate points defining ellipse */

	dtheta = 2.0 * cPI / (double) (npts - 1);
	for (n = 0, theta = 0.0; n < npts; n++, theta += dtheta)
	{
		/* unrotated ellipse */
		xloc = p1 * cos(theta);
		yloc = p2 * sin(theta);
		/* rotate by angle */
		xyarray[n].x = xcent + xloc * cosa - yloc * sina;
		xyarray[n].y = ycent + yloc * cosa + xloc * sina;
	}

	return(xyarray);

}



/** function to write array of Vect2D to GMT psxy data */

int Vect2DArray2GMT(FILE* fp_io, Vect2D* array, int npts)
{
	int n;


	for (n = 0; n < npts; n++)
		fprintf(fp_io, "%lf %lf\n", array[n].x, array[n].y);

	fprintf(fp_io, ">\n");

	return(0);

}



/** method to convert 2 semi-axis in 3D to an array of points
		defining  an error ellipse (covariance ellipse) */

/* converted from method toEllipsoid3D in Java class Ellipsoid3D (04DEC1998) */

/* !!! remember to free ellArray */

Vect3D* toEllipsoid3D(Vect3D ax1, Vect3D ax2, Vect3D center, int npts)
{
	Vect3D *ellArray;
	double d_angle;
	double angle;
	int n;
	double cosang, sinang;


	ellArray = malloc((size_t) npts * sizeof(Vect3D));

	d_angle = 2.0 * cPI / (double) (npts - 1);
	angle = 0.0;

	for (n = 0; n < npts; n++) {

		cosang = cos(angle);
		sinang = sin(angle);
		ellArray[n].x = center.x
			+ ax1.x * cosang + ax2.x * sinang;
		ellArray[n].y = center.y
			+ ax1.y * cosang + ax2.y * sinang;
		ellArray[n].z = center.z
			+ ax1.z * cosang + ax2.z * sinang;
			angle += d_angle;
	}

	return(ellArray);

}



/** function to project 3D array of points to 2D */

Vect2D* Vect3D2To2D(Vect3D *array3d, Vect2D* array2d, int npts, int iComps)
{
	int n;

	if (iComps == 12) {
		for (n = 0; n < npts; n++)
		{
			array2d[n].x = array3d[n].x;
			array2d[n].y = array3d[n].y;
		}
	} else if (iComps == 13) {
		for (n = 0; n < npts; n++)
		{
			array2d[n].x = array3d[n].x;
			array2d[n].y = array3d[n].z;
		}
	} else if (iComps == 32) {
		for (n = 0; n < npts; n++)
		{
			array2d[n].x = array3d[n].z;
			array2d[n].y = array3d[n].y;
		}
	}

	return(array2d);

}


/*------------------------------------------------------------/ */
/* Anthony Lomax           | email: lomax@faille.unice.fr     / */
/* UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax / */
/* 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        / */
/* 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        / */
/*------------------------------------------------------------/ */

