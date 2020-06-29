/*
   calculate magnitudes
 */
#include <stdio.h>
#include <math.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/types.h>
#include <stdlib.h>
#include <string.h>


#include "complex.h"
#include "grid_search.h"
#include "new_sedlib.h"

#define DEBUG 0

float c_magnitude();
float median_select();
void ctransfer2();
void WoodAnderson();

// AJL 20040527 added args: o_event_type, g13_depth
void magnitude(int jmin, float* mag, float* magerr, char* cmaxmag, int* num_mag, char o_event_type, double g13_depth)
{
   int nm,nmedian,nclip, got_WA;
   int i,j,k,l;
   float avmag,sqmag,maxmag,tmag,clmag;
   float edist,cormag=0.0,median[300],clmedian[300];
   char tchan[6];
   
   int cindx[100], Ncindx;

   if (DEBUG) printf("TP 1\n");
   /** strcpy(tchan,"SHZ"); **/
   tchan[0]= 0;
   got_WA= 0; /* if WA amplitudes available then ignore clipped data */
   nclip= 0;
   for (i=0;i<nphas;i++) {
      phases[i].mag= 0.0;
      phases[i].magrmk = '?';
// AJL 20040527 added args: o_event_type
//      if (location_st[0].o_event_type[0] == 'T' )  edist= -phases[i].delta;
      if (o_event_type == 'T' )  edist= -phases[i].delta;
//
      else edist= phases[i].delta*111.17;

// AJL 20040527 k here not used!
//      k= phases[i].stntab_index;

   if (DEBUG) printf("TP 2\n");
      GetCal4StnTime(phases[i].name,jmin,&cindx,&Ncindx);
   if (DEBUG) printf("TP 3 Ncindx=%d\n", Ncindx);

      phases[i].mag= -13.0;
      if ( Ncindx <= 0 ) continue;
      l= -1;
      if ( phases[i].rec_sys[0] != 'N' ) { /* any recording system other than Nanometrics */
   if (DEBUG) printf("TP 4\n");
         for (j=0;j<Ncindx;j++) {
	    if (calib[cindx[j]].comp[2] == 'Z' &&
	       calib[cindx[j]].auxid[0] == phases[i].rec_sys[0]) {
	       l= cindx[j];
	       break;
	    }
	 }
      }
      else {
	 if ( phases[i].rec_sys[1] == 'W' ) {   /*   Wood-Anderson   */
   if (DEBUG) printf("TP 5\n");
            for (j=0;j<Ncindx;j++) {
   if (DEBUG) printf("TP 5.1 %c\n", calib[cindx[j]].auxid[0]);
	       if ( calib[cindx[j]].auxid[0] == 'W') {   /*   Wood-Anderson   */
	          l= cindx[j];
	          break;
	       }
	    }
	    got_WA= 1;
	 }
	 else {
   if (DEBUG) printf("TP 6\n");
	    /* search Z-comp */
            for (j=0;j<Ncindx;j++) {
	       if (calib[cindx[j]].comp[2] == 'Z' &&
	          calib[cindx[j]].auxid[0] == phases[i].rec_sys[0]) {
	          l= cindx[j];
	          break;
	       }
	    }	       

	 }
      }
      if ( l != -1 ) {
   if (DEBUG) printf("TP 7\n");
// AJL 20040528 added args: g13_depth
         phases[i].mag= c_magnitude(calib[l], edist, (float) g13_depth,
				    phases[i].amplitude,
				    phases[i].period,cormag);
      }
   if (DEBUG) printf("TP 7.1 %lf\n",phases[i].mag);
      if ( phases[i].mag > -13.0 ) {
   if (DEBUG) printf("TP 7.2 %d\n",phases[i].clipped);
	if(phases[i].clipped < 0) {
	   phases[i].magrmk = '<';
	}
	else if ( phases[i].clipped == 0 ){
	   phases[i].magrmk = ' ';
	}
        else phases[i].magrmk = '?';
      }
      else phases[i].magrmk = '?';
   }
// AJL 20040527 added args: o_event_type
//   switch (location_st[0].o_event_type[0]) {
   switch (o_event_type) {
//
      case 'T': 	/* teleseismic events */
   if (DEBUG) printf("TP 8\n");
         nm= 0;
         avmag= 0.0;
         sqmag= 0.0;
         maxmag= 0.0;
         nmedian= 0;
	 for (i=0;i<nphas;i++) {
	    if ( strcmp(phases[i].phase,"P") == 0 ) {
	       if ( phases[i].magrmk == ' ') {
	          nm++;
      	          avmag+= phases[i].mag;
      	          sqmag+= phases[i].mag*phases[i].mag;
	          median[nmedian]= phases[i].mag;
	          nmedian++;
	       }
	       else if ( phases[i].magrmk == '<') {
		  clmedian[nclip]= phases[i].mag;
                  nclip++;
	       }
	    }
	    else {
	       phases[i].magrmk = '-';
	    }
	 }
	 break;
      default:		/* local and region events */
   if (DEBUG) printf("TP 9\n");
	 for (i=0;i<nphas;i++) {
	    tmag= phases[i].mag;
   if (DEBUG) printf("TP 9.1\n");
	    for(j=i+1;j<nphas;j++) {
   if (DEBUG) printf("TP 9.11 %s %s\n",phases[i].name,phases[j].name);
	       if( strcmp(phases[i].name,phases[j].name) == 0) {
		  if ( phases[j].mag > tmag ) {
	             if ( phases[i].magrmk == ' ') phases[i].magrmk = '-';
		  }
		  else  {
	             if ( phases[j].magrmk == ' ') phases[j].magrmk= '-';
		  }
	       }
   if (DEBUG) printf("TP 9.12 %c %c\n",phases[i].magrmk, phases[j].magrmk);
	    }
	 }
   if (DEBUG) printf("TP 9.2\n");
         nm= 0;
         avmag= 0.0;
         sqmag= 0.0;
         maxmag= 0.0;
         nmedian= 0;
	 for (i=0;i<nphas;i++) {
   if (DEBUG) printf("TP 9.3 %c\n", phases[i].magrmk);
	    if ( phases[i].magrmk == ' ') {
	       nm++;
      	       avmag+= phases[i].mag;
      	       sqmag+= phases[i].mag*phases[i].mag;
	       median[nmedian]= phases[i].mag;
	       nmedian++;
	    }
	    else if ( phases[i].magrmk == '<') {
	       clmedian[nclip]= phases[i].mag;
               nclip++;
	    }
	 }
	 break;
   }
   if (DEBUG) printf("TP 10\n");
   tmag= 0;
   if ( nm > 2 ) tmag= sqrt(((float)nm*sqmag-avmag*avmag)/((float)nm*(float)(nm-1)));
   if ( nm > 0 ) avmag/=(float)nm;
   if( nmedian > 1) {
   if (DEBUG) printf("TP 11\n");
      k= (nmedian-1)/2;
      avmag= median_select(k,nmedian,median);
      if( (k+1)*2 == nmedian) {
          k++;
          avmag+= median_select(k,nmedian,median);
	  avmag*=0.5;
      }
   }
   cmaxmag[0]= ' ';
   if( nclip > 1 && got_WA == 0 ) {
   if (DEBUG) printf("TP 12\n");
      k= (nclip-1)/2;
      clmag= median_select(k,nclip,clmedian);
      if( (k+1)*2 == nclip) {
          k++;
          clmag+= median_select(k,nclip,clmedian);
	  clmag*=0.5;
      }
      if ( clmag > avmag ) {
	  avmag= clmag;
	  if ( avmag > 2.5 ) cmaxmag[0]= '<';
      }
   }
   if (DEBUG) printf("TP 13\n");
   *mag= avmag;
   *magerr=tmag;
   *num_mag=nm;
}
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;

float median_select(k,n,arr)
float arr[];
unsigned int k,n;
{
	unsigned int i,ir,j,l,mid;
	float a,temp;

	l=0;
	ir=n-1;
	for (;;) {
		if (ir <= l+1) {
			if (ir == l+1 && arr[ir] < arr[l]) {
				SWAP(arr[l],arr[ir])
			}
			return arr[k];
		} else {
			mid=(l+ir) >> 1;
			SWAP(arr[mid],arr[l+1])
			if (arr[l+1] > arr[ir]) {
				SWAP(arr[l+1],arr[ir])
			}
			if (arr[l] > arr[ir]) {
				SWAP(arr[l],arr[ir])
			}
			if (arr[l+1] > arr[l]) {
				SWAP(arr[l+1],arr[l])
			}
			i=l+1;
			j=ir;
			a=arr[l];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				SWAP(arr[i],arr[j])
			}
			arr[l]=arr[j];
			arr[j]=a;
			if (j >= k) ir=j-1;
			if (j <= k) l=i;
		}
	}
}
#undef SWAP
/************ C-Version for magnitude calculation ***********************/
float c_magnitude(cal, epdist,depth, ampl,eperiod,stacor)
   CalEntry cal;
   float epdist,depth,ampl,eperiod,stacor;
{
   /* C-Adaption of muk2 (fortran) */

   float epdistkm,disp,veloc,accel,perctg,waampl;
   float a,b,delta,sigma,mag;
   int i;
   //int resultok;
   //char comment[80];
      
   static float rdelt[] = {4.0,6.,13.,16.,28.,87.,114.,120.,134.,141.,146.,180.};
      /* c SIGAR values scaled for amplitude in nano-meters */
   static float sigar[] = {3.2,4.,4.0,2.8,3.6,4.0,4.5, 4.5, 3.9, 4.1, 3.9, 3.9 };

if (DEBUG) printf("TP c_m 1 %s %f %f %f %f %f\n", cal.name, epdist, depth, ampl, eperiod, stacor);
   ampl/= 2.0;		/* peak-peak --> 0-peak */
   mag= -13.;			/* falls xmag nicht anders berechnet wird !! */
   if(eperiod == 0.0 || ampl == 0.0) return(mag);
   if(epdist < 0.) {
         epdistkm= -epdist/360.*40030;    /* epdist < 0. --> [epdist] = Grad */
   }
   else {
         epdistkm=sqrt(epdist*epdist+depth*depth);
   }

if (DEBUG) printf("TP c_m 2a %s %f %f %f %f %f %f %f\n", cal.name, ampl,eperiod,disp,veloc,accel,perctg,waampl);
   ctransfer2(cal,ampl,eperiod,&disp,&veloc,&accel,&perctg,&waampl);
if (DEBUG) printf("TP c_m 2b %s %f %f %f %f %f %f %f\n", cal.name, ampl,eperiod,disp,veloc,accel,perctg,waampl);

   if(epdistkm > 2200.) {
      delta= -epdist;    /* delta:=epdistkm in Grad */
      if(delta > 180.) {
            delta=180.;
            sigma=sigar[11];
      }  
      else {
	 i=0;
         while ( delta > rdelt[i]) i++;
         /*****************************************
            jetzt ist delta < rdelt(i)     --> Linearisieren
            y=a*x+b
         *****************************************/
         a=(sigar[i]-sigar[i-1])/(rdelt[i]-rdelt[i-1]);
         b=sigar[i]-a*rdelt[i];
         sigma=a*delta+b;
      }
      mag=log10(disp/eperiod) + sigma + stacor ;
   }
   else {
      if( (epdistkm > 0) && ( epdistkm <= 60.) ) {
            mag=log10(waampl) + 0.018 *epdistkm+1.77 + 0.40;
      }
      if( (epdistkm > 60.) && (epdistkm <= 700.) ) {
            mag=log10(waampl) + 0.0038*epdistkm+2.62 + 0.40;
      }
      if( (epdistkm > 1100.) && (epdistkm <= 1700.) ){
            mag=log10(waampl) + 0.0029*epdistkm+3.40 + 0.40 - 2.;  /* EMPIRISCH */
      }
      if( mag > -13.)  {
	 mag+= stacor;  /* Stationskorrektur Cs */
         if ( cal.auxid[0] == 'W' ) mag-= 0.3;  /* vertical to horizontal factor: log(2) */
      }
   }
   return( mag );
}
void ctransfer2(cal,ampl,period,disp,veloc,accel,perctg,waampl)
   CalEntry cal;
   float ampl,period;
   float *disp,*veloc,*accel,*perctg,*waampl;
{
   float om, g=9.81e+9;
   float rzwa[2],izwa[2],rpwa[2],ipwa[2],fac,scalefactor, *x;
   int nz,np,i,k;
   dcomplex jom, j, allpoles, allzeros, t;
   
   j= Complex((float) 0, (float) 1);
   om= twopi/period;
   jom= RCmul(om,j);

   allzeros= Complex((float) 1,(float) 0);
   allpoles= Complex((float) 1,(float) 0);
   scalefactor= 1.0;
   
   for ( k=0; k<cal.Nstage; k++) {
      if ( cal.stage[k].digit != NULL ) {
         scalefactor*= cal.stage[k].digit->counts_volt;
         continue;
      }
      if ( cal.stage[k].paz == NULL ) continue;
      x= cal.stage[k].paz->zeros;
      for(i=0;i<cal.stage[k].paz->nzeros*2;i+=2) {
         allzeros= Cmul(allzeros,
		     Csub(jom,Complex(x[i],x[i+1])));
      }
      x= cal.stage[k].paz->poles;
      for(i=0;i<cal.stage[k].paz->npoles*2;i+=2) {
         allpoles= Cmul(allpoles,
		    Csub(jom,Complex(x[i],x[i+1])));
      }
      scalefactor*= cal.stage[k].paz->factor;
   }

   t= RCmul((float)(scalefactor*cal.gain),Cdiv(allzeros,allpoles));
      
   *disp= ampl/Cabs(t);
   *veloc=(*disp)*twopi/period;
   *accel=(*veloc)*twopi/period;
   *perctg=(*accel)*twopi*100.0/(g*period);
   
   if ( cal.auxid[0] != 'W' ) {
if (DEBUG) printf("TP ct2 2a %d %f %f %f %d %f %f\n", nz,fac,rzwa[0],izwa[0],np,rpwa[0],ipwa[0]);
      WoodAnderson(&nz,rzwa,izwa,&np,rpwa,ipwa,&fac);
if (DEBUG) printf("TP ct2 2b %d %f %f %f %d %f %f\n", nz,fac,rzwa[0],izwa[0],np,rpwa[0],ipwa[0]);

      for(i=0;i<nz;i++) {
         allzeros= Cmul(allzeros, Csub(jom,Complex(rzwa[i],izwa[i])));
      }
      for(i=0;i<np;i++) {
         allpoles= Cmul(allpoles, Csub(jom,Complex(rpwa[i],ipwa[i])));
      }
if (DEBUG) printf("TP ct2 3a1 %f\n", t.r);
if (DEBUG) printf("TP ct2 3a2 %f %f %f %f %f %f %f\n", fac, scalefactor, cal.gain, allzeros.r,allzeros.i,allpoles.r,allpoles.i);
      t= RCmul((float)(fac*scalefactor*cal.gain),Cdiv(allzeros,allpoles));
if (DEBUG) printf("TP ct2 3b2 %f %f %f %f %f %f %f\n", fac, scalefactor, cal.gain, allzeros.r,allzeros.i,allpoles.r,allpoles.i);
if (DEBUG) printf("TP ct2 3b1 %f\n", t.r);
   }
   *waampl= ampl/Cabs(t);
   
}
void WoodAnderson(nz,rzwa,izwa,np,rpwa,ipwa,fac)
   int *nz, *np;
   float *rzwa,*izwa,*rpwa,*ipwa,*fac;
{
   float Twa= 0.8, hw=0.78;
   float om0,sqrthw;
   //int i;
  
   om0= twopi/Twa;
   sqrthw=sqrt(1.0-hw*hw);
   rzwa[0]= -hw;
   izwa[0]= sqrthw;
   rzwa[1]= -hw;
   izwa[1]= -sqrthw;
   *nz= 2;
   rpwa[0]= 0.0;
   ipwa[0]= 0.0;
   rpwa[1]= 0.0;
   ipwa[1]= 0.0;
   *np= 2;
   *fac= 1.0e+6/WA_gain;
   
}
