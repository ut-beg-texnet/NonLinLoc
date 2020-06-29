/* SH 02/27/2004
    routines to convert lat/lon into Swiss coordinates and to find region name and number
    of local earthquakes in Switzerland
	original routine from M. Baer  */
	
#include <stdio.h>
#include <math.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/types.h>
#include <stdlib.h>

#define EXTERN_MODE 1
#include "../GridLib.h"
#include "eth_functions.h"
//#include "get_region_name_nr.h"
#include "new_sedlib.h"
//#include "gsetype.h"
//#include "gselib.h"
#include "get_region_name_nr.h"


int get_region_name_nr(lat,lon,type,reg_nr,name,lx,ly)
   float lat, lon;
   int type;
   int   *reg_nr;
   char  *name;
   float *lx,*ly;
{
   //char reg_type[4];
   float x,y,xx,yy;
   int ityp, i;
   Aname region;
   
   if(lat < 45 || lat > 49.0) ityp= 1; /* F-E */
   else {
      x=lon;
      y=lat;
      celleb(x,y,&xx,&yy);
      if(yy < 62.0 || yy > 302.0 || xx < 480.0 || xx > 847.5) {
	 *lx=999;
	 *ly=999;
      }
      else {
         *lx=xx;
         *ly=yy;
	 lat=yy;
	 lon=xx;
      }
   }

   i= set_up_names(snap_file_PATH1,
		   snap_file_PATH2,
		   snap_file_PATH3);
// AJL 20040527
   if (i > 0)
      return(i);

   region= region_attributes(lat,lon);
   *reg_nr= region.reg_number;
   strcpy(name,region.name);
   return(i);
} /* end of get_region_name_nr */


Aname region_attributes(lat,lon)
   float lat,lon;
{
   int func_value;
   int invert= 1;
   int i,k,hit,nv;
   Aname error;
   float y,x;
   /*******************************************/
   hit=-1;
   y=lat;
   x=lon;
   if ( invert ) {      
      for(i=0;i<nseg;i++) {
	     k= a_name[i].start;
            nv= a_name[i].npairs;
	     func_value=inside(y,x,&la_tab[k],&lo_tab[k],nv);
         if ( func_value ) { 
	        return(a_name[i]); 
	     }
      }
   }
   else {
      for(i=0;i<nseg;i++) {
	 k= a_name[i].start;
         nv= a_name[i].npairs;
      
	 func_value= inside(x,y,&lo_tab[k],&la_tab[k],nv);
         if ( func_value ) { 
	    return(a_name[i]); 
	 }
      }
   }
   error.reg_number=-1;
   return(error);
}
int set_up_names(fe,ch,bi)
   char *fe,*ch,*bi;
{
   FILE *fp_in;	/* Pointer to input  file   */
   char line[100],string[100];
   //int i;
   //Aname prev_name, *B_name;

   /*** input of area file ********************/
   if ((fp_in = fopen(bi, "rb")) == NULL) {
      fprintf(stderr,"input file %s not found --> try to create it\n",bi);
      nseg=-1;
      np= 0;
      a_name=(Aname *) malloc(sizeof(Aname));
      
      /* open and read fe96 file */
      if ((fp_in = fopen(fe, "r")) == NULL) {
         fprintf(stderr,"input file %s not found --> terminate\n",fe);
	 return(1);
      }

      lo_tab=(float *) malloc(sizeof(float));
      la_tab=(float *) malloc(sizeof(float));
      fgets(line,100,fp_in);
      while ( feof(fp_in) == 0  ) {
         if ( line[0] == '#' ) {
            fgets(line,100,fp_in);
	    continue;
         }
         if ( line[0] == '>' ) {
	    line[strlen(line)-1]= 0;
	    if( nseg > -1 ) {
	       a_name[nseg].npairs=np-a_name[nseg].start;
	    }
	    nseg++;
	    a_name=(Aname *) realloc(a_name,sizeof(Aname)*(nseg+1));
	    csplitstring(line,string,line,':');
	    csplitstring(line,string,line,'.');
	    sscanf(string,"%d",&a_name[nseg].reg_number);
	    csplitstring(line,string,a_name[nseg].name,' ');
	    a_name[nseg].attribute[0]= 0;
	    a_name[nseg].start= np;
	    reg_index[a_name[nseg].reg_number]= nseg;
            fgets(line,100,fp_in);
	    continue;
         }
         lo_tab=(float *) realloc(lo_tab,sizeof(float)*(np+1));
         la_tab=(float *) realloc(la_tab,sizeof(float)*(np+1));
         sscanf(line,"%f %f",&la_tab[np],&lo_tab[np]);
         np++;
         fgets(line,100,fp_in);
      }
      a_name[nseg].npairs=np-a_name[nseg].start;
      // AJL 20040510 corection
      //if ( fclose(fp_in) != NULL) fprintf(stderr,"Warning: unable to close file.\n");
      if ( fclose(fp_in) != 0 ) fprintf(stderr,"Warning: unable to close file.\n");

      /* open and read ch25 file */
      if ((fp_in = fopen(ch, "r")) == NULL) {
         fprintf(stderr,"input file %s not found --> terminate\n",ch);
	 return(1);
      }

      fgets(line,100,fp_in);
      while ( feof(fp_in) == 0  ) {
         if ( line[0] == '#' ) {
            fgets(line,100,fp_in);
	    continue;
         }
         if ( line[0] == '>' ) {
	    line[strlen(line)-1]= 0;
	    if( nseg > -1 ) {
	       a_name[nseg].npairs=np-a_name[nseg].start;
	    }
	    nseg++;
	    a_name=(Aname *) realloc(a_name,sizeof(Aname)*(nseg+1));
	    csplitstring(line,string,line,':');
	    sscanf(string,"%*s %d",&a_name[nseg].reg_number);
	    csplitstring(line,a_name[nseg].attribute,a_name[nseg].name,' ');
	    reg_index[a_name[nseg].reg_number]= nseg;
	    a_name[nseg].start= np;
            fgets(line,100,fp_in);
	    continue;
         }
         lo_tab=(float *) realloc(lo_tab,sizeof(float)*(np+1));
         la_tab=(float *) realloc(la_tab,sizeof(float)*(np+1));
         sscanf(line,"%f %f",&la_tab[np],&lo_tab[np]);
         np++;
         fgets(line,100,fp_in);
      }
      a_name[nseg].npairs=np-a_name[nseg].start;
      nseg++;
      // AJL 20040510 corection
      //if ( fclose(fp_in) != NULL) fprintf(stderr,"Warning: unable to close file.\n");
      if ( fclose(fp_in) != 0 ) fprintf(stderr,"Warning: unable to close file.\n");
      if ((fp_in = fopen(bi, "wb")) == NULL) {
         fprintf(stderr,"input file %s not found --> terminate\n",bi);
	 return(1);
      }

      fwrite(&nseg,sizeof(int),1,fp_in);
      fwrite(&np,sizeof(int),1,fp_in);
      fwrite(la_tab,sizeof(float),np,fp_in);
      fwrite(lo_tab,sizeof(float),np,fp_in);
      fwrite(reg_index,sizeof(int),1500,fp_in);
      fwrite(a_name,sizeof(Aname),nseg,fp_in);

      // AJL 20040510 corection
      //if ( fclose(fp_in) != NULL) fprintf(stderr,"Warning: unable to close file.\n");
      if ( fclose(fp_in) != 0 ) fprintf(stderr,"Warning: unable to close file.\n");
   }
   else {
      fread(&nseg,sizeof(int),1,fp_in);
      fread(&np,sizeof(int),1,fp_in);
      
      la_tab=(float *) malloc(sizeof(float)*np);
      lo_tab=(float *) malloc(sizeof(float)*np);      
      fread(la_tab,sizeof(float),np,fp_in);
      fread(lo_tab,sizeof(float),np,fp_in);
      fread(reg_index,sizeof(int),1500,fp_in);	  
      
      a_name=(Aname *) malloc(sizeof(Aname)*nseg);
      fread(a_name,sizeof(Aname),nseg,fp_in);
      fclose(fp_in);
   }
   return(0);
}
void celleb(l,b,y,x)
/**********************
C
C  SCHWEIZ. PROJEKTIONSSYSTEM  FORMELN VON H. ODERMATT
C  TRANSFORMATION ELLIPSOID - EBENE
C  L,B  laenge und breite in grad
C  Y,X LANDESKOORDINATEN IN KILO-METER y= e-w; x=n-s
C  MY  MERIDIANKONVERGENZ ( SEXAG. SEK.)
C
*************************/

float l;  /* [deg] */
float b;  /* [deg] */
float *y; /* [km] */
float *x; /* [km] */
{
#define TOP 8
   double bb = 169028.66;
   double bl = 26782.5;
   short i;
   double d[9],e[9],f[9],rw[9],iw[9],p,q,a,c,my,dx,dy;

   a= l;
   a = a*3600.0 - bl;
   c = b;
   c = c * 3600.0 - bb;
   d[1]= 0.68382546262761;
   d[2]=-3.91798328045E-8;
   d[3]= 1.4965410352E-15;
   d[4]=-8.039471422E-23;
   d[5]= 7.0021390E-30;
   d[6]=-5.586904E-37;
   d[7]= 4.0402E-44;
   d[8]=-3.06E-51;
   e[1]= 2.3635916074715E-2;
   e[2]= 0.;
   e[3]= 4.527219881E-17;
   e[4]=-3.89081120E-24;
   e[5]= 2.3407700E-31;
   e[6]=-1.59674E-38;
   e[7]= 1.287E-45;
   e[8]= 0.0;
   f[1]= 4.515344386039E1;
   f[2]= 1.17912305209E-4;
   f[3]= 5.8474201864E-10;
   f[4]= 2.73386187E-15;
   f[5]= 1.4308547E-20;
   f[6]= 7.66562E-26;
   f[7]= 4.2445E-31;
   f[8]= 2.40E-36;

   p = 30.91849390613 * a;
   q = c * f[8];
   i = TOP;
   label1:
   i = i - 1;
   q = c * (q + f[i]);
   if (i > 1) goto label1;

   rw[1] = q;
   iw[1] = p;
   for(i=2;i<=TOP;i++)
   {
      rw[i]=q*rw[i-1]-p*iw[i-1];
      iw[i]=p*rw[i-1]+q*iw[i-1];
   }

   dx=d[TOP]*rw[TOP];
   dy=d[TOP]*iw[TOP];
   my= 0.0;
   i=TOP;
   label3:
   i=i-1;
   dx=dx+d[i]*rw[i];
   dy=dy+d[i]*iw[i];
   my=my+e[i]*iw[i];
   if (i > 1) goto label3;
   dx=dx+200000.;
   dy=dy+600000.;
   dy= dy/1000.;
   dx= dx/1000.;
   *y =  dy;
   *x =  dx;
}

void cebell(YL,XB,L,B)
float YL; /* [km] */
float XB; /* [km] */
float *L; /* [deg] */
float *B; /* [deg] */
{
/*********************************
C
C  SCHWEIZ. PROJEKTIONSSYSTEM FORMELN VON H. ODERMATT
C  TRANSFORMATION EBENE ELLIPSOID
C  YL,XB LANDESKOORDINATEN IN KILOMETER (BERN 600/200)
C  L,B LAENGE UND BREITE  (GRAD)
C  XY MERIDIEANKONVERGENZ (SEXAG.SEK.)
C
***********************************/

#define TOP 8
   double A[9],BB[9],C[9],RZ[9],QZ[9],P,Q,X,Y,DB,DL,MY;
   short I;
   XB=XB*1000;
   YL=YL*1000;
   X=XB-200000.;
   Y=YL-600000.;
   A[1] =1.4623614572021;
   A[2] =1.225255821052E-7;
   A[3] =1.3687923002E-14;
   A[4] =1.971224191E-21;
   A[5] =2.97898051E-28;
   A[6] =4.650273E-35;
   A[7] =7.48203E-42;
   A[8] =1.229E-48;
   BB[1]=3.4564252673326E-2;
   BB[2]=2.89600437564E-9;
   BB[3]=4.651046030E-16;
   BB[4]=6.43850954E-23;
   BB[5]=9.600412E-30;
   BB[6]=1.50512E-36;
   BB[7]=2.422E-43;
   C[1] = 2.2146704979846E-2;
   C[2] =-1.280815253730E-9;
   C[3] = 7.4775676024E-18;
   C[4] = 4.691943327E-24;
   C[5] =-3.6550101E-31;
   C[6] = 3.71615E-39;
   C[7] = 1.6901E-45;
   C[8] = 1.96E-52;
   RZ[1]= X;
   QZ[1]= Y;
   for (I=2;I<=TOP;I++)
   {
      RZ[I]= X*RZ[I-1]-Y*QZ[I-1];
      QZ[I]= Y*RZ[I-1]+X*QZ[I-1];
   }
   I=TOP;
   Q=A[I]*RZ[I];
   P=A[I]*QZ[I];
   MY=0.;
   label2:
   I=I-1;
   Q=Q+A[I]*RZ[I];
   P=P+A[I]*QZ[I];
   MY=MY+BB[I]*QZ[I];
   if (I > 1) goto label2;
   DL=3.2343101932327E-2 * P;
   DB=Q*C[TOP];
   I=8;
   label3:
   I=I-1;
   DB=Q*(DB+C[I]);
   if(I > 1) goto label3;
   *L = ((DL + 26782.5E00)/3600.);
   *B = ((DB + 169028.66E00)/3600.);
}
void cgrkr(lat1,lon1,lat2,lon2,dist,az)
double lat1, lon1, lat2, lon2, *dist, *az;
{

   double bog, fq, th1, th2, ph2, x1, z1, x2, y2, z2;
   double p12, xa, za, xb, zb, pab, dab;
  
   bog=atan(1.0)/45.0;
   fq=.9933;
   th1=atan(fq*tan(bog*lat1));
   th2=atan(fq*tan(bog*lat2));
   ph2=(lon2-lon1)*bog;
   x1=cos(th1);
   z1=sin(th1);
   x2=cos(th2)*cos(ph2);
   y2=cos(th2)*sin(ph2);
   z2=sin(th2);
   p12=x1*x2+z1*z2;
   *dist=acos(p12)/bog;
   xa=x2-x1*p12;
   za=z2-z1*p12;
   xb=-x1*z1;
   zb=1.0-z1*z1;
   pab=xa*xb+za*zb;
/**  dab=x1*y2*zb+y1*za*xb-xb*y2*z1-zb*xa*y1; **/
   dab=x1*y2*zb         -xb*y2*z1;
   *az=atan2(dab,pab)/bog;
}
int ccoord(lat1,lon1,dist,az,lat2,lon2)
float lat1,lon1,dist,az,*lat2,*lon2;
{
   float galor,galar,xazr,xdr,gelac,xlatr,gelor,ck1,ck2;
   float pi,pi2,pirad,q,as,ac,alon;
      
   pi= atan(1.0)*4.0;
   pirad= atan(1.0)/45.0;
   pi2=pi*2.; 
   /**   CONVERT ARRAY COORDINATES TO RADIANS **/
   galor=lon1*pirad;
   galar=lat1*pirad;

   xazr=  az * pirad;
   xdr=  dist*pirad;

   q=sin(galar)*cos(xdr)+cos(galar)*sin(xdr)*cos(xazr);
   if(fabs(q) > 1.0) {
      *lon2=0.0;
      *lat2=0.0;
      return(1);
   }
   as=asin(q);
   gelac=as;
   xlatr=gelac; 
   
   gelor= 0.0;
   if ( cos(gelac) && cos(galar) ) {
      ck1=(cos(xdr)-sin(galar)*sin(xlatr))/(cos(galar)*cos(xlatr)) ;
      ck2=(sin(xazr)*sin(xdr))/cos(xlatr) ;
      if ( (ck1 >= 0.0) && (ck2 >= 0.0)) {
	 if ( fabs(ck1) <= 1.0 ) {
	    ac= acos(ck1);
	    gelor=galor+ac;
	 }
      }
      else if ( (ck1 < 0.0) && (ck2 < 0.0)) {
	 if ( fabs(ck2) <= 1.0 ) {
	    as= asin(ck2);
	    gelor=galor-pi-as;
	 }
      }
      else if (ck1 < 0.0 ) {
	 if ( fabs(ck1) <= 1.0 ) {
	    ac= acos(ck1);
	    gelor=galor+ac;
	 }
      }
      else if (ck2 < 0.0 ) {
	 if ( fabs(ck2) <= 1.0 ) {
	    as= asin(ck2);
	    gelor=galor+as;
	 }
      }
   }
/***********************************************************
      CONVERT LONGITUDE AND LATITUDE FROM RADIANS TO DEGREES
      AND REDUCE LONGITUDE TO A FORM LESS THAN 360 DEGREES
***********************************************************/
   alon= fabs(gelor);
   while ( alon >= pi2 ) {
      alon-= pi2;
   }

   if( gelor < 0.0 ) gelor=-fabs(alon);
   else              gelor= fabs(alon);
   if ( fabs(gelor) > pi) {
      if( gelor < 0.0 ) gelor+= pi2;
      else              gelor-= pi2;
   }
/***********************************************************
   LONGITUDE IN DEGREES 
***********************************************************/
   *lon2= gelor/pirad;
/***********************************************************
   LATITUDE IS CONVERTED FROM GEOCENTRIC TO GEODETIC COORDINATES 
***********************************************************/
   *lat2= gelac/pirad;
   return(0);
}

/*
int csplitstring(s1,result,s2,c)
   char *s1, *result, *s2, c;
 */
   /******************************************************* 
      split a string 's1' into result (first part) and
                               's2' (remaining part)
      where the terminating character for is 'c'
   
   example: strcpy(s1,"beispiel a-b;c");
            csplitstring("beispiel a-b;c",result,s2,' ');
         --> result: "beispiel"
             s2:     "a-b;c"
            csplitstring("beispiel a-b;c",result,s2,'-');
         --> result: "beispiel a"
             s2:     "b;c"
   *******************************************************/
   
/*{
   int i,j,l,k,m;
   
   l= strlen(s1);
   m=0;
   for (i=0;i<l;i++) {
      if( s1[i] == c ) {
	    result[i]= 0;
	    s2[0]= 0;
	    k=i+1;
	    for (j=k;j<l;j++) {
	       if ( s1[j] == c ) k++;
	       else              break;
	    }
	    strcpy(s2,&s1[k]);
	    return (0);
      }
      else result[i]= s1[i];
      m++;
   }
   s2[0]= 0;
   result[m]=0;
   return(0);
}*/

int inside (float x0,float y0,float* px,float* py,float n)
	
/*	GODKIN,//.B. AND J.J.PULLI: APPLI//ATION OF THE "WINDING-NUMBER
	ALGORITHM" TO THE SPATIAL SORTING OF //ATALOGED EARTHQUAKE DATA.
	BSSA 74, 5, PP. 1845-1848, O//TOBER 1984

	CHECK IF POINT X0,Y0 LIES INSIDE POLYGONE

	PARAMETERS:
	X0,Y0    POINT TO TEST
	PX       ARRAY OF X-//OORDINATES OF POLYGON
	PY       ARRAY OF Y-//OORDINATES OF POLYGON
       N        NUMBER OF ELEMENTS OF PX AND PY OR OF SIDES OF POLYGON

	RETURN VALUE:    0  IF POINT OUTSIDE
					+/-1  IF POINT INSIDE 
					2  IF POINT IS ON AN EDGE OR VERTEX  
					
	routine was converted into C code  from original Fortan code SH 03/02/2004  */

{
	int i;
	int ins=0;
	int isicr=0;
	float x1,x2,y1,y2;

/*	ACCUMULATE SIGNED CROSSING NUMBERS WITH INSIDE */
	
	for (i=0; i < n-1; i++) {
	
/*	   PROCEED AROUND POLYGON CHECKING EACH SEGMENT TO SEE IF 
	   NEGATIVE X-AXIS WAS CROSSED
       TRANSLATE COORDINATES OF POLYGON TO PUT TEST POINT AT ORIGIN */
	  x1 = px[i] - x0;
	  x2 = px[i + 1] - x0;
	  y1 = py[i] - y0;
	  y2 = py[i + 1] - y0;
	  isicr = KSICR(x1, y1, x2, y2);
	  
	  
/*	 STOP COUNTING IF POINT IS ON EDGE  */
	  if (isicr == 4) {
	    return(2);
	   }
	  ins = ins + isicr;
       } /* for (i=0... */

/*       CHECK SEGMENT FROM LAST VERTEX TO FIRST VERTEX
    mb: not necessary for closed area. next 2 statements commented
 	    ISICR=KSICR(PX(N)-X0,PY(N)-Y0,PX(1)-X0,PY(1)-Y0)
	    IF(ISICR .EQ. 4) GO TO 600
	    INSIDE=(INSIDE+ISICR)/2  */
	ins = (int)ins/2;
/*	type *, ' N= ',N,' INSIDE=',INSIDE,' ISICR=',ISICR */

	return(ins);
}

int  KSICR(float x1, float y1, float x2, float y2)
/*	COMPUTE SIGNED CROSSING NUMBER

	A "SIGNED CROSSING NUMBER" TELLS WETHER A SEGMENT 
	(IN THIS CASE THE SEGMENT FROM (X1,Y1) TO (X2,Y2))
	CROSSES THE NEGATIVE X-AXIS OR GOES THROUGH THE ORIGIN

	THE RETURN VALUES ARE:
		  +2 IF SEGMENT CROSSES FROM BELOW
	      +1 IF SEGMENT EITHER ENDS ON -X-AXIS FROM BELOW OR STARTS
                           UPWARDS FROM -X-AXIS ("HALF CROSSING")
	       0 IF THERE IS NO CROSSING
	      -1 IF SEGMENT EITHER ENDS ON -X-AXIS FROM ABOVE OR STARTS
                           DOWNWARDS FROM -X-AXIS ("HALF CROSSING")
		  -2 IF SEGMENT CROSSES FROM ABOVE
		  +4 IF SEGMENT CROSSES THROUGH THE ORIGIN  
		  
       routine was converted into C code  from original Fortan code SH 03/02/2004  */

{
	
// IF BOTH POINTS ARE ON THE SAME SIDE OF X-AXIS, RETURN 0
	if (y1*y2 > 0) return(0);

// CHECK IF SEGMENT CROSSES THROUGH THE ORIGIN
	if (x1*y2 != x2*y1 || x1*x2 > 0) goto label100;
	return(4);
	
// check for complete crossing
 label100 : if (y1*y2 < 0) goto label300;
// half crossing ?
// one end of segment touches x-axis! which end?
       if (y2 == 0) goto label200;
// here y1=0  check if segment touches +x-axis
       if (x1 > 0) return(0); // no crossing
// upward or downward ?
       if (y2 > 0) return(1); // upward half crossing
       return(-1); // downward half crossing
       
// here y2=0  check if segment touches +x-axis
label200 :  if (y1 == 0 || x2 > 0) return(0);  // no crossing
// upward or downward ?
       if (y1 > 0) return(-1); // downward half crossing
       return(1); // upward half crossing
       
// complete crossing of -x-axis ?
// break into cases according to crossing direction
 label300 :  if (y1 > 0) goto label400;
// case y1<y2
       if (x1*y2 >= y1*x2) return(0); // no crossing
       return(2); // upward crossing
       
// case y1 > 0 > y2
 label400 :  if (y1*x2 >= x1*y2) return(0); // no crossing
       return(-2);  // downward crossing
       

} /* end of KSCIR */
