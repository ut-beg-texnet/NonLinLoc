/*****************************************************************************
* NAME       gselib.c                                                        *
* VERSION    1.0                                                             *
* DATE       25.11.1991                                                      *
* AUTHOR     H. Maurer                                                       *
******************************************************************************
*  The GSELIB library is a collection of voids and functions, which          *
*  deals with i/o's of GSE data files. It includes the headerfile            *
*  ngsetype.h which contains the neccessary type definitions.                *
*                                                                            *
*  The library contains the following voids and functions:                   *
*                                                                            *
*  - juliam     (Converts a date of the form " year, month, day, hour, min " *
*	 	into a minute expression. This expression gives the number   *
*	 	of minutes after B.C. (julian minute)).                      *
*                                                                            *
*  - datum      (This is the opposite of juliam. A julian minute is converted*
*	 	into a date).                                                *
*                                                                            *
*  - juliam_KP   Converts an sed filename (e.g. KP198412112345) into a julian*
*                minute                                                      *
*                                                                            *
*  - datum_KP    Creates an sed filename from a julian minute                *
*                                                                            *
*  - juliad  	(Converts a date of the form  "month, day" to the julian day)*
*                                                                            *
*  - juliadtodate  (Opposite of juliad).                                     *
*                                                                            *
*  - ReadGseHeader (Reads all header informations of a GSE data file into    *
*	    	   a structure of the type masterheader (see gsetype.h)      *
*                                                                            *
*  - Comp6toBuf    (This void performs 6 bit data compression of an 4 byte   *
*	    integer array into a global defined character buffer).           *
*                                                                            *
*  - Decomp6toArray (Decompression alorithmus for a character buffer created *
*	     with Comp6toBuf).                                               *
*                                                                            *
*  - MakeDiff     (Calculates the n'th differences of a 4 byte integer       *
*	     	   array).                                                   *
*                                                                            *
*  - RemoveDiff     (Opposite of MakeDiff).                                  *
*****************************************************************************/

#ifdef StrC
#     define strcmp StrCmp
#     define strncmp StrNCmp
#endif

#ifndef MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "gsetype.h" 
#include "gselib.h"
// AJL 20040518 change
//#include "sedlib.h"
#include "new_sedlib.h"
#endif
long int juliam_KP(filename)
char *filename;
/* juliam_KP returns the julian minute calculated from an sed filename */
   {
	static long kmo[13]={0,0,31,59,90,120,151,181,212,243,273,304,334};
	long leap=1;
	long ky;
	long km;
	long kd;
	long temp;
	long l=0;
	long kl;
        long ky4;
        long ky1;
        long ky0;
        long khour;
        long kmin;
	char temp_str[5];
	int pos;

        if (filename[0] == '.')
	{
	   pos = 1;
        }
	else
	{
	   pos = 0;
        }
	strncpy(temp_str,filename+2+pos,4);
	temp_str[4]=0;
	ky =  atol(temp_str);
	strncpy(temp_str,filename+6+pos,2);
	temp_str[2]=0;
	km =  atol(temp_str);
	strncpy(temp_str,filename+8+pos,2);
	temp_str[2]=0;
	kd =  atol(temp_str);
	strncpy(temp_str,filename+10+pos,2);
	temp_str[2]=0;
	khour =  atol(temp_str);
	strncpy(temp_str,filename+12+pos,2);
	temp_str[2]=0;
	kmin =  atol(temp_str);

	if (km <= 0)
		km = 1;
        if (kd <= 0)
                kd = 1;
	temp = 365 * ky;
	kd = kmo[(short)km] + kd;
	ky4  = ky / 4;
	ky1  = ky / 100;
	ky0  = ky / 1000;
	kl   = leap * (ky4 - ky1 + ky0);
	if (
             (  (ky4 *4) == ky
             )
             &&
             (  (  (ky1 * 100 ) != ky ) ||
                (  (ky0 * 1000) == ky)
             )
           )
		l = leap;
	if ((l != 0) && (km < 3))
		kl = kl - leap;
	temp = temp + kd + kl;
	temp = temp * 24 + khour;
	temp = temp * 60 + kmin;
	return(temp);
	}

void datum_KP(itf,filename)
long *itf;
char *filename;
   {
	static long kmo[13]={0,31,28,31,30,31,30,31,31,30,31,30,31};
	long ky;
	long km;
	long kd;
	long khour;
	long kmin;

	long temp;
        long k;
        long kh;
	long  ky4;
	long  ky1;
	long  ky0;
	long  l;
	long  ld;
	short  a;
	short  syear;
	short  sday;
	short  shour;
	short  smin;

	temp = *itf;
	k    = temp/60;
	kmin = temp - k * 60;
	kh   = k/24;
	khour	= k - kh * 24;
	ky    = kh/365;
label5:
	kd    = kh - ky * 365;
	l     = 0;
	ky4   = ky/4;
	ky1   = ky/100;
	ky0   = ky/1000;
	ld    = ky4 - ky1 + ky0;
	if   (
	     (  (ky4 *4) == ky
	     )
	     &&
	     (  (  (ky1 * 100 ) != ky ) ||
		(  (ky0 * 1000) == ky)
	     )
	     )
		    l = 1;
	kd = kd - ld + l;
	if (kd > 0) goto label10;
	if ((kd == 0) && ((kmin == 0) && (khour == 0)))
		{
		/*************
		*day = 0;
		*month = 0;
		**************/
		}
	ky--;
	goto label5;
label10:
	kmo[2] = 28 + l;
	for (a=1;a<=12;a++)
		{
		kd = kd - kmo[a];
		if (kd <= 0) goto label30;
		}
	a = 12;
label30:
	kd = kd + kmo[a];
	syear  = (short)ky;
	sday   = (short)kd;
	shour  = (short)khour;
	smin   = (short)kmin;
	sprintf(filename,"KP%04hd%02hd%02hd%02hd%02hd",syear,a,sday,shour,smin);
	}   

   
short int juliad(year,mon,day)
short year;
short mon;
short day;

{
static short nrofdays[13] = {0,0,31,59,90,120,151,181,212,243,273,304,334};

if ((mon >= 3) && (!(year % 4)))
    return(day + nrofdays[mon] + 1);
else
    return(day + nrofdays[mon]);
}

void juliadtodate(jday,year,tmon,tday)
short jday;
short year;
short *tmon;
short *tday;

{
static short nrofdays[14] = {0,0,31,59,90,120,151,181,212,243,273,304,334,365};
short int a;

if ((year % 4) == 0)
   {
   for(a=3;a<=13;a++)
     nrofdays[a]++;
   }
a = 0;
while(nrofdays[a] < jday)
    a++;
(*tmon) = a -1;
(*tday) = jday - nrofdays[a-1];
if ((year % 4) == 0)
   {
   for(a=3;a<=13;a++)
     nrofdays[a]--;
   }
}


int GseBrowse(infile)
char *infile;
{
FILE           *gsef;
char           line[82];
int            nrwf;


/* Trying to open the GSE data file                                       */
gsef = fopen(infile,"rt");
if (gsef==NULL)
   {
      fprintf(stderr,"\nUnable to open file %s \n\n",infile);
   exit(1);
   }
else
   {
   nrwf = 0;
   do
	{
	fgets(line,82,gsef);
	if (!strncmp(line,"WID1 ",5))
	   nrwf++;
	}     /* do while strncmp ...*/
   while  (!feof(gsef));
   } /* else */
fclose(gsef);
return(nrwf);
}     /* GseBrowse */





void ReadGseHeader(infile,mhp)
char *infile;
masterheader *mhp;
{
int            a;
int            b;
int            c;
FILE           *gsef;
char           line[82];
long           date;
short int      flag;
char           temp_str[20];
char           tmp3[3];
float          no_sec_f;
int            old_station;


/* Trying to open the GSE data file                                       */
gsef = fopen(infile,"rt");
if (gsef==NULL)
   {
      fprintf(stderr,"\nUnable to open file %s \n\n",infile);
   exit(1);
   }
else
   /* Init some variables                                                 */
   {
   a = 0;
   mhp->nst = 0;
   mhp->nwf = GseBrowse(infile);
   flag = 1;

   /* Allocating memory for the badsamp array */
   mhp->bads = (short **)malloc(sizeof(short *));
   if (mhp->bads == NULL)
       {
	  fprintf(stderr,"\nMemory allocation error ReadGseHeader (1) !\n");
	  exit(1);
       }

/* Allocating memomry for the gseinfo array */
   mhp->gseinfo = (gseh *)malloc(sizeof(gseh) * (size_t)((*mhp).nwf + 1));
   if (mhp->gseinfo == NULL)
       {
	  fprintf(stderr,"\nMemory allocation error ReadGseHeader (2) !\n");
	  exit(1);
       }
   memset(mhp->gseinfo,0,sizeof(gseh) * (size_t)((*mhp).nwf + 1));
/* Allocating memomry for the wfindex array */
   mhp->wfindex = (str9 *)malloc(sizeof(str9) * (size_t)((*mhp).nwf + 1));
   if (mhp->wfindex == NULL)
       {
	  fprintf(stderr,"\nMemory allocation error ReadGseHeader (3) !\n");
	  exit(1);
       }
   memset(mhp->wfindex,0,sizeof(str9) * (size_t)((*mhp).nwf + 1));

   /* Searching for WID1                                                  */
   while (flag)
     {
     do
	{
	fgets(line,82,gsef);

	/* Setting the flag, if the STOP id is reached                    */
/* New Version 18.11.1991: The STOP marker in the GSE file is now ignored */
/*	if (!(strncmp("STOP",line,4)))                                   */
	if (feof(gsef))
	   flag = 0;


/********************************************************************************************
 *	 Error message, if an EOF is detected and the flag is on
 *	if (feof(gsef) && flag)
 *	   {
 *	   fprintf(stderr,"\n\nEnd of GSE-file %s reached without finding \"STOP\"\n\n",infile);
 *	   exit(1);
 *	   }
 **********************************************************************************************/

	}
     while  ((strncmp("WID1",line,4)) && (flag));

     /* If WID1 is found  and the flag is on then the header lines of
	the waveform datas are read  to a temporary structure gh          */
     if (flag)
       {
       a++;
       mhp->gseinfo[a].use_flag = 1;
       strcpy(mhp->gseinfo[a].name,"");
       strcpy(mhp->gseinfo[a].chan,"");
       strcpy(mhp->gseinfo[a].chid,"");
       strcpy(mhp->gseinfo[a].systype,"");
       strcpy(mhp->gseinfo[a].form,"");
       strcpy(temp_str,"");
       sscanf(line+5,"%8ld",&date);
       sscanf(line+14,"%2hd",&(mhp->gseinfo[a].tfsamp.hour));
       sscanf(line+17,"%2hd",&(mhp->gseinfo[a].tfsamp.min));
       sscanf(line+20,"%2hd",&(mhp->gseinfo[a].tfsamp.sec));
       sscanf(line+23,"%3hd",&(mhp->gseinfo[a].tfsamp.msec));
       sscanf(line+27,"%8ld",&(mhp->gseinfo[a].nosamp));
       sscanf(line+36,"%6s",mhp->gseinfo[a].name);
       sscanf(line+43,"%8s",mhp->gseinfo[a].chid);
       sscanf(line+52,"%2s",mhp->gseinfo[a].chan);
       sscanf(line+55,"%11s",temp_str);
       mhp->gseinfo[a].sra = (float) atof(temp_str);
       sscanf(line+67,"%6s",mhp->gseinfo[a].systype);
       sscanf(line+74,"%4s",mhp->gseinfo[a].form);
       sscanf(line+79,"%1hd",&(mhp->gseinfo[a].dflag));
       while (strlen(mhp->gseinfo[a].name) < 6)
       {
	  strncat(mhp->gseinfo[a].name," ",1);
       }
       mhp->gseinfo[a].name[6] = '\0';
       while (strlen(mhp->gseinfo[a].chan) < 2)
       {
	  strncat(mhp->gseinfo[a].chan," ",1);
       }
       mhp->gseinfo[a].chan[2] = '\0';
       mhp->gseinfo[a].chid[8] = '\0';
       mhp->gseinfo[a].form[4] = '\0';
       mhp->gseinfo[a].systype[6] = '\0';

       /* Determing the julian day and the year                           */
       mhp->gseinfo[a].tfsamp.year = (short int)(date / 1000);
       mhp->gseinfo[a].tfsamp.jday = (short int)(date % 1000);
       juliadtodate(mhp->gseinfo[a].tfsamp.jday,
		    mhp->gseinfo[a].tfsamp.year,
		    &(mhp->gseinfo[a].tfsamp.month),
		    &(mhp->gseinfo[a].tfsamp.day));
       mhp->gseinfo[a].tfsamp.jmin = juliam(&(mhp->gseinfo[a].tfsamp.year),
					    &(mhp->gseinfo[a].tfsamp.month),
					    &(mhp->gseinfo[a].tfsamp.day),
					    &(mhp->gseinfo[a].tfsamp.hour),
					    &(mhp->gseinfo[a].tfsamp.min));

       /* The station name and the channel identifier are formatted and
	  then concatenated to the waveform index                         */
       strcpy(mhp->wfindex[a],"");
       strcpy((*mhp).wfindex[a],mhp->gseinfo[a].name);
       if ((!strncmp(mhp->gseinfo[a].chid,"MLR2",4))  &&
	   (strncmp(mhp->gseinfo[a].chid+6,"NG",2)))
	   {
	      strncpy(mhp->wfindex[a]+4,mhp->gseinfo[a].chid+6,2);
	      mhp->wfindex[a][6] = 0;
	   }  
       if (strncmp(mhp->gseinfo[a].chid+6,"LG",2))
       {
	      strncpy(mhp->wfindex[a]+4,mhp->gseinfo[a].chid+6,2);
	      mhp->wfindex[a][6] = 0;
       }
	  
       while(strlen((*mhp).wfindex[a]) < 6)
	    strcat((*mhp).wfindex[a]," ");
       while(strlen(mhp->gseinfo[a].chan) < 2)
	    strcat(mhp->gseinfo[a].chan," ");
       strncat((*mhp).wfindex[a],mhp->gseinfo[a].chan,2);
       
       (*mhp).wfindex[a][8] = '\0';
       while(strlen(mhp->gseinfo[a].systype) < 6)
	     strcat(mhp->gseinfo[a].systype," ");
       /* Reading the second header line of the waveform id               */
       fgets(line,82,gsef);

       strcpy(temp_str,"");
       strncpy(temp_str,line,9);
       mhp->gseinfo[a].gain = (float) atof(temp_str);
       strcpy(temp_str,"");
       strncpy(temp_str,line+9,1);
       temp_str[1] = '\0';
       mhp->gseinfo[a].units = (short) atoi(temp_str);
       strcpy(temp_str,"");
       strncpy(temp_str,line+10,7);
       temp_str[7] = '\0';
       mhp->gseinfo[a].cperiod = (float) atof(temp_str);
       strcpy(temp_str,"");
       strncpy(temp_str,line+18,9);
       temp_str[9] = '\0';
       mhp->gseinfo[a].lat = (float) atof(temp_str);
       strcpy(temp_str,"");
       strncpy(temp_str,line+28,9);
       temp_str[9] = '\0';
       mhp->gseinfo[a].lon = (float) atof(temp_str);
       strcpy(temp_str,"");
       strncpy(temp_str,line+38,9);
       temp_str[9] = '\0';
       mhp->gseinfo[a].alt = (float) atof(temp_str);
       strcpy(temp_str,"");
       strncpy(temp_str,line+48,9);
       temp_str[9] = '\0';
       mhp->gseinfo[a].depsens = (float) atof(temp_str);

       /* Check, if the station name is already presen */
       if (a == 1)
       {
	  mhp->nst ++;
       }
       else
       {
	  old_station = 0;
          for (b=1;b<=(a-1);b++)
          {
             if (!strcmp(mhp->gseinfo[a].name,(*mhp).gseinfo[b].name))
             {
	         old_station = 1;  
             }
          }
	  if (old_station)
	  {
	     mhp->nst ++;
	  }
       }
       mhp->gseinfo[a].stat_nr = mhp->nst;
       if ((!strncmp(mhp->gseinfo[a].chan,"TC",2)) &&
	   (!strncmp(mhp->gseinfo[a].chid,"MLR2",4)))
       {
	  mhp->bads = (short **)realloc(mhp->bads,sizeof(short *) * (size_t)(mhp->nst + 1));
          if (mhp->bads == NULL)
          {
             fprintf(stderr,"\nMemory allocation error ReadGseHeader (4) !\n");
             exit(1);
          }
	  no_sec_f = (float)mhp->gseinfo[a].nosamp/mhp->gseinfo[a].sra;
	  mhp->gseinfo[a].no_sec = (size_t)floor(no_sec_f + 0.5);
	  c = a-1;
	  while(!strncmp(mhp->gseinfo[c].name,mhp->gseinfo[a].name,6))
	  {
	     mhp->gseinfo[c].no_sec = (size_t)floor(no_sec_f + 0.5);
	     c--;
	  }  
	  mhp->bads[mhp->nst] = (short *)malloc(mhp->gseinfo[a].no_sec * sizeof(short));
	  if (mhp->bads[mhp->nst] == NULL)
	  {
	     fprintf(stderr,"\nMemory allocation error ReadGseHeader (5) !\n");
	     exit(1);
	  }
	  for (b=0;b<mhp->gseinfo[a].no_sec;b++)
	  {
	     fscanf(gsef,"%hd",&(mhp->bads[mhp->nst][b]));
	  }
       }
       do
       {
	  fgets(line,82,gsef);
       }
          while(strncmp(line,"DAT1",4));
       /* Determing the number of 512 byte blocks of the datas and
	  calculating the amount of dummy values to fill the last
	  block up to 512 bytes                                           */
       if (strncmp(mhp->gseinfo[a].form,"CMP6",4))  /* if not compressed  */
	   {
	   mhp->gseinfo[a].nblock = mhp->gseinfo[a].nosamp/128 + 1;     /* (+ 1) are added for:        */
					      /*   - 1 for the remaining     */
					      /*     datas and               */
	   mhp->gseinfo[a].dummyval = 128 - (mhp->gseinfo[a].nosamp % 128);
	   if (mhp->gseinfo[a].dummyval == 0)
	      mhp->gseinfo[a].nblock--;
	   mhp->gseinfo[a].cmpbytes = 0;
	   while ((strncmp("CHK1 ",line,5)))
		 fgets(line,82,gsef);
	   sscanf(line,"%*s %15ld",&mhp->gseinfo[a].checksum);
	   }
       else                 /* if compressed  */
	   {
	   b = 0;
	   do
		{
		fgets(line,82,gsef);
		b++;
       /* New Version 18.11.1991                                      */
		if (feof(gsef))
       /* 		if (feof(gsef) || (!strncmp(line,"STOP",4)))   */
		   {
		      fprintf(stderr,"End of file %s reached without finding CHK1",infile);
		   exit(1);
		   }
		}
	   while ((strncmp("CHK1 ",line,5)));
	   sscanf(line,"%*s %15ld",&mhp->gseinfo[a].checksum);
	   mhp->gseinfo[a].cmpbytes = 80 * (--b);
	   mhp->gseinfo[a].dummyval = 512 - (mhp->gseinfo[a].cmpbytes % 512);
	   mhp->gseinfo[a].nblock = mhp->gseinfo[a].cmpbytes / 512 + 1; /* (+ 2) are added for:        */
					      /*   - 1 for the remaining     */
					      /*    datas and               */
					      /*  - 1 for the signal header */
	   if (mhp->gseinfo[a].dummyval == 0)
	      mhp->gseinfo[a].nblock--;
	   }
       }  /* if flag */
     }    /* while flag */
   if ((*mhp).nwf != a)
      {
	 fprintf(stderr,"\nError while reading the GSE data headers\n");
      exit(1);
      }
   }      /* else */
fclose(gsef);
}     /* ReadGseHeader */





char *Comp6toBuf(in4,nrofint,nrofchar,cbuf)
long *in4;
long nrofint;
long *nrofchar;
char *cbuf;
{
   /* Adjust this statement for your memory model or operating system */

#ifdef UNIX
#    define MBS 128*1024*1024
#endif
#ifdef MSDOS
#    define MBS 65530
#endif

   static char    LookUp[65]=
	       {' ','+','-','0','1','2','3','4','5','6','7','8','9',
		'A','B','C','D','E','F','G','H','I','J','K','L','M',
		'N','O','P','Q','R','S','T','U','V','W','X','Y','Z',
		'a','b','c','d','e','f','g','h','i','j','k','l','m',
		'n','o','p','q','r','s','t','u','v','w','x','y','z'};
   long n4 = 16;
   long n5 = 32;
   long n5m1 = 31;
   long n6  = 64;
   long n9  = 512;
   long n10 = 1024;
   long n10m1 = 1023;
   long n14  = 16384;
   long n15  = 32768;
   long n15m1 = 32767;
   long n19  = 524288;
   long n20  = 1048576;
   long n20m1 = 1048575;
   long n24  = 16777216;
   long n25  = 33554432;
   long n25m1 = 33554431;
   long n28  = 134217728;
   long n28m1 = 134217727;
   long mflag = 32;
   long nflag;
   long a;
   long count = 0;
   long jn;
   long j;
   long nblank;

   cbuf = malloc(1);
   if (cbuf == NULL)
   {
      fprintf(stderr,"Memory allocation error in Comp6toBuf (1) !\n");
      exit(1);
   }
   for (a=0;a<nrofint;a++)
   {
      nflag = 1;
      jn = in4[a];
      if (jn < 0)
      {
          nflag = nflag + n4;
          jn = -jn;
      }
      if (jn >= n4 )
      {
         /* more than 1 byte */
         if (jn >= n9)
	      {
	         /* more than 2 bytes */
	         if (jn >= n14)
	         {
	            /* more than 3 bytes */
	            if (jn >= n19)
		         {
		            /* more than 4 bytes */
		            if (jn >= n24)
		            {
		               /* more than 5 bytes */
		               if (jn >= n28m1)
		               jn = n28m1;
		               j = jn/n25 + nflag + mflag;
		               if ((count + 1) > MBS)
		               {
		                  fprintf(stderr,"\nCompression buffer became too large!!\n");
		                  exit(1);
                     }
		               if ((cbuf = realloc(cbuf,(size_t)(count + 1))) == NULL)
		               {
		                  fprintf(stderr,"\nMemory allocation error  in Comp6toBuf (2) !\n");
		                  exit(1);
		               }
		               cbuf[count] = LookUp[j];
		               count++;
		               jn = jn & n25m1;
		               nflag = 1;
		            }
		            j = (jn/n20 + nflag + mflag);
		            if ((count + 1) > MBS)
		            {
		               fprintf(stderr,"\nCompression buffer became too large!!\n");
		               exit(1);
                  }
		            if ((cbuf = realloc(cbuf,(size_t)(count + 1))) == NULL)
		            {
		               fprintf(stderr,"\nMemory allocation error in Comp6toBuf (3) !\n");
		               exit(1);
		            }
		            cbuf[count] = LookUp[j];
		            count++;
		            jn = jn & n20m1;
		            nflag = 1;
		         }
	            j = (jn/n15 + nflag + mflag);
	            if ((count + 1) > MBS)
		         {
		            fprintf(stderr,"\nCompression buffer became too large!!\n");
		            exit(1);
               }
	            if ((cbuf = realloc(cbuf,(size_t)(count + 1))) == NULL)
		         {
		            fprintf(stderr,"\nMemory allocation error in Comp6toBuf (4) !\n");
		            exit(1);
		         }
	            cbuf[count] = LookUp[j];
	            count++;
	            jn = jn & n15m1;
	            nflag = 1;
	         }
	         j =(jn/n10 + nflag + mflag);
	         if ((count + 1) > MBS)
		      {
		         fprintf(stderr,"\nCompression buffer became too large!!\n");
		         exit(1);
            }
	         if ((cbuf = realloc(cbuf,(size_t)(count + 1))) == NULL)
	         {
	            fprintf(stderr,"\nMemory allocation error in Comp6toBuf (5) !\n");
	            exit(1);
	         }
	         cbuf[count] = LookUp[j];
	         count++;
	         jn = jn & n10m1;
	         nflag = 1;
	      }
         j = (jn/n5 + nflag + mflag);
         if ((count + 1) > MBS)
		   {
		      fprintf(stderr,"\nCompression buffer became too large!!\n");
		      exit(1);
         }
         if ((cbuf = realloc(cbuf,(size_t)(count + 1))) == NULL)
	      {
	         fprintf(stderr,"\nMemory allocation error in Comp6toBuf (6) !!\n");
	         exit(1);
	      }
         cbuf[count] = LookUp[j];
         count++;
         jn = jn & n5m1;
         nflag = 1;
      }
      j = (jn + nflag);
      if ((count + 1) > MBS)
		{
		   fprintf(stderr,"\nCompression buffer became too large!!\n");
		   exit(1);
      }
      if ((cbuf = realloc(cbuf,(size_t)(count + 1))) == NULL)
      {
         fprintf(stderr,"\nMemory allocation error in Comp6toBuf (7) !\n");
         exit(1);
      }
      cbuf[count] = LookUp[j];
      count++;
   }
   nblank = 80 - (count % 80);
   if (nblank < 2) nblank += 80;
   for (a=0;a<nblank;a++)
   {
      if ((count + 1) > MBS)
		{
		   fprintf(stderr,"\nCompression buffer became too large!!\n");
		   exit(1);
      }
      if ((cbuf = realloc(cbuf,(size_t)(count + 1))) == NULL)
      {
         fprintf(stderr,"\nMemory allocation error in Comp6toBuf (8) !\n");
         exit(1);
      }
      cbuf[count] = ' ';
      count++;
   }
   (*nrofchar) = count;
   return(cbuf);
}

void Decomp6toArray(lb,lout,iout,cbuf)
long lb;
long *lout;
long *iout;
char *cbuf;
{
static long ichar[129] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		   0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,2,
		   3,4,5,6,7,8,9,10,11,0,0,0,0,0,0,0,
		   12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,
		   28,29,30,31,32,33,34,35,36,37,0,0,0,0,0,0,
		   38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,
		   54,55,56,57,58,59,60,61,62,63,0,0,0,0,0,0};
static char test[5] = {'1','2','3','4','\0'};
char achar[5];
char aspace = 32;
char lfeed = 10;
char cretn = 13;
long *itest;
long *inn;
long imax;
long isign;
long jsign;
long ioflow;
long joflow;
long mask1;
long mask2;
long ibyte;
long icount;
long i;
long j;
long k;
long icheck;
long itemp;

imax = (*lout);
isign   = 16;
ioflow = 32;
mask1   = 15;
mask2   = 31;
ibyte   =  3;
itest = (long *)test;
inn   = (long *)achar;
/* Work out which way bytes are stored in computer */
icheck = (long)256 * (long)256 * (long)256 * (long)52;
icheck += (long)256 * (long)256 * (long)51;
icheck += (long)256 * (long)50;
icheck += (long)49;
if ((*itest) == icheck)
    {
    ibyte = 0;
    }
icount = 0;
i = 0;
j = 0;
achar[4] = '\0';
label1:
i++;
achar[(unsigned int)ibyte] = cbuf[(unsigned int)i - 1];
/* If a carriage return or line feed ignore, if a space then end of data */
if (achar[(unsigned int)ibyte] == cretn) goto label1;
if (achar[(unsigned int)ibyte] == lfeed) goto label1;
if (achar[(unsigned int)ibyte] == aspace) goto label5;
icount++;
/* Strip off any higher order byte   */
k = (*inn) & (long)127;
/* Get number representation of input character  */
(*inn) = ichar[(unsigned int)k];
/* Get sign bit */
jsign = (*inn) & isign;
/* get continuation bit (if any)   */
joflow = (*inn) & ioflow;
/* Remove bits we dont want    */
itemp = (*inn) & mask1;
label2:
if (joflow == 0) goto label4;
/* There is another byte in this sample   */
itemp *= 32;
label3:
i++;
achar[(unsigned int)ibyte] = cbuf[(unsigned int)i - 1];
if (achar[(unsigned int)ibyte] == cretn) goto label3;
if (achar[(unsigned int)ibyte] == lfeed) goto label3;
icount++;
/* Strip off any higher order bits   */
k = (*inn) & (long)127;
(*inn) = ichar[(unsigned int)k];
/* get continuation bit (if any)   */
joflow = (*inn) & ioflow;
k = (*inn) & mask2;
itemp += k;
goto label2;
label4:
if (jsign != 0) itemp *= -1;
j++;
if (j > imax) goto label5;
iout[(unsigned int)j-1] = itemp;


if (icount < lb) goto label1;
label5:
if (j > imax)
   {
   (*lout) = imax;
   }
}

void MakeDiff(in4,nelements,ndiff)
long *in4;
long nelements;
short ndiff;
{
   int a,b;

   b = 0;
   while (b < ndiff)
   {
      b++;
      for (a=nelements-1;a>=1;a--)
	   {
	      in4[a] = in4[a] - in4[a-1];
	   }
   }
}

void RemoveDiff(in4,nelements,ndiff)
long *in4;
long nelements;
short ndiff;
{
   int a,b;

   b = 0;
   while (b < ndiff)
   {
      b++;
      for (a=1;a<nelements;a++)
	   {
	      in4[a] += in4[a-1];
	   }
   }
}

int ReadGseFromMemory(g_ptr,gsize,mhp,stnfile,ncal)
char *g_ptr;
int gsize;
masterheader *mhp;
char *stnfile;
int *ncal;
{
   int  a,b,c;
   long date;
   char temp_str[20];
   float  no_sec_f;
   int  old_station;
   int  pos;
   char line[128];
   int fok;
   FILE *inf;
   int nstat;
   station_struct *stn;
   int ngse2;
   int ncal_local;


   /* Init */
   a = 0;
   mhp->nst = 0;
   mhp->nwf = 0;
   ngse2 = 0;
   *ncal = 0;
   ncal_local= 0;

   
   /* Allocating memomry for the gseinfo array */
   mhp->gseinfo = (gseh *)malloc(sizeof(gseh));
   if (mhp->gseinfo == NULL)
   {
      fprintf(stderr,"\nMemory allocation error in ReadGseFromMemory (1) !\n");
      exit(1);
   }
   memset(mhp->gseinfo,0,sizeof(gseh));
   
   /* Allocating memomry for the wfindex array */
   mhp->wfindex = (str9 *)malloc(sizeof(str9));
   if (mhp->wfindex == NULL)
   {
      fprintf(stderr,"\nMemory allocation error in ReadGseFromMemory (2) !\n");
      exit(1);
   }
   memset(mhp->wfindex,0,sizeof(str9));

   /* Main loop of reading procedure */
   pos = 0;
   a = 0;
   while (1)
   {
      while((memcmp(&g_ptr[pos],"WID1 ",5)) && (memcmp(&g_ptr[pos],"WID2 ",5)))
      {
	 if (!(strncmp(&g_ptr[pos],"CAL1 ",5)) || (!strncmp(&g_ptr[pos],"CAL2 ",5)))
	 {
	    ncal_local++;
            *ncal= ncal_local;
         }
	 pos++;
	 if ((pos+5) >= gsize)
	 {
	    if (nstat > 0) 
	    {
	       free(stn);
	       stn = NULL;
	    }
	    if (mhp->nwf > 0)
	    {
	       return(1);
	    }
	    else
	    {
	       fprintf(stderr,"Warning: No waveforms found in gse file%s\n");
	       return(0);
	    }
	 }
      }
      
      /* if here, a WID1 or WID2 identifier is found */
      a++;
      mhp->gseinfo = (gseh *)realloc(mhp->gseinfo,sizeof(gseh)*(a+1));
      if (mhp->gseinfo == NULL)
      {
	 fprintf(stderr,"Memory allocation error in ReadGseFromMemory (3) \n");
	 exit(1);
      }
      memset(&mhp->gseinfo[a],0,sizeof(gseh));
      mhp->wfindex = (str9 *)realloc(mhp->wfindex,sizeof(str9)*(a+1));
      if (mhp->wfindex == NULL)
      {
         fprintf(stderr,"\nMemory allocation error in ReadGseFromMemory (4) !\n");
         exit(1);
      }
      memset(&mhp->wfindex[a],0,sizeof(str9));
      mhp->gseinfo[a].use_flag = 1;
      if (g_ptr[pos+3] == '1')
      {
	 /* if here, gse type WID1 section is found */
	 if ((pos+80) >= gsize)
	 {
	    fprintf(stderr,"Error in module ReadGseToMemory: Inappropriate input file\n");
	    if (nstat > 0) 
	    {
	       free(stn);
	       stn = NULL;
	    }
	    return(0);
	 }
	 strncpy(line,g_ptr+pos,80);
	 line[80] = 0;
         sscanf(line+5,"%8ld",&date);
         sscanf(line+14,"%2hd",&(mhp->gseinfo[a].tfsamp.hour));
         sscanf(line+17,"%2hd",&(mhp->gseinfo[a].tfsamp.min));
         sscanf(line+20,"%2hd",&(mhp->gseinfo[a].tfsamp.sec));
         sscanf(line+23,"%3hd",&(mhp->gseinfo[a].tfsamp.msec));
         mhp->gseinfo[a].tfsamp.year = (short int)(date / 1000);
         mhp->gseinfo[a].tfsamp.jday = (short int)(date % 1000);
         juliadtodate(mhp->gseinfo[a].tfsamp.jday,
	              mhp->gseinfo[a].tfsamp.year,
		    &(mhp->gseinfo[a].tfsamp.month),
		    &(mhp->gseinfo[a].tfsamp.day));
	 mhp->gseinfo[a].tfsamp.fsec = (float)mhp->gseinfo[a].tfsamp.sec +
					(float)mhp->gseinfo[a].tfsamp.msec/(float)1000.0;
         sscanf(line+27,"%8ld",&(mhp->gseinfo[a].nosamp));
         sscanf(line+36,"%6s",mhp->gseinfo[a].name);
         sscanf(line+43,"%8s",mhp->gseinfo[a].chid);
         sscanf(line+52,"%2s",mhp->gseinfo[a].chan);
         sscanf(line+55,"%f",&mhp->gseinfo[a].sra);
         sscanf(line+67,"%6s",mhp->gseinfo[a].systype);
         sscanf(line+74,"%4s",mhp->gseinfo[a].form);
         sscanf(line+79,"%1hd",&(mhp->gseinfo[a].dflag));
         while (strlen(mhp->gseinfo[a].name) < 6)
         {
	    strcat(mhp->gseinfo[a].name," ");
         }
         while (strlen(mhp->gseinfo[a].chan) < 2)
         {
	    strcat(mhp->gseinfo[a].chan," ");
         }       
	 mhp->gseinfo[a].name[6] = '\0';
	 mhp->gseinfo[a].chan[2] = '\0';
         mhp->gseinfo[a].chid[8] = '\0';
         mhp->gseinfo[a].form[4] = '\0';
         mhp->gseinfo[a].systype[6] = '\0';
	 
         /* The station name and the channel identifier are formatted and
	    then concatenated to the waveform index */
         strcpy(mhp->wfindex[a],mhp->gseinfo[a].name);
         if ((!strncmp(mhp->gseinfo[a].chid,"MLR2",4))  &&
	      (strncmp(mhp->gseinfo[a].chid+6,"NG",2)))
         {
            strncpy(mhp->wfindex[a]+4,mhp->gseinfo[a].chid+6,2);
	    mhp->wfindex[a][6] = 0;
         }  
         if (strncmp(mhp->gseinfo[a].chid+6,"LG",2))
         {
            strncpy(mhp->wfindex[a]+4,mhp->gseinfo[a].chid+6,2);
	    mhp->wfindex[a][6] = 0;
         }
	  
         while(strlen(mhp->wfindex[a]) < 6) strcat((*mhp).wfindex[a]," ");
         strncat(mhp->wfindex[a],mhp->gseinfo[a].chan,2);
       
         mhp->wfindex[a][8] = '\0';
         while(strlen(mhp->gseinfo[a].systype) < 6) strcat(mhp->gseinfo[a].systype," ");
	 
	 strcpy(temp_str,mhp->gseinfo[a].name);
	 while (strlen(temp_str) < 6) strcat(temp_str," ");
	 
	 /* MLR stuff */	       
         if ((!strncmp(mhp->gseinfo[a].chid,"MLR2",4))  &&
              (strncmp(mhp->gseinfo[a].chid+6,"NG",2)))
         {
            strncpy(temp_str+6,mhp->gseinfo[a].chid+6,2);
            if (!strncmp(temp_str+6,"TC",2))
            {
               strncpy(mhp->gseinfo[a].name+4,"TC",2);
               temp_str[6] = 0;
               temp_str[7] = 0;
	    }
            if (!strncmp(temp_str+6,"HG",2))
            {
               strncpy(temp_str+6,"SZ",2);
	    }
            if (!strncmp(temp_str+6,"MG",2))
            {
               strncpy(temp_str+6,"SM",2);
	    }
            if (!strncmp(temp_str+6,"LG",2))
            {
               strncpy(temp_str+6,"SS",2);
	    }
            if (!strncmp(temp_str+6,"GR",2))
            {
               strncpy(temp_str+6,"SG",2);
	    }
	    temp_str[8] = 0;
	 }
 		
         /* All strong motion channels */
         if (!strncmp(mhp->wfindex[a],"LG",2))
         {
            strncpy(temp_str+6,"SS",2);
	    temp_str[8] = 0;
	 }
	 if (strlen(temp_str) < 8) strcat(temp_str,mhp->gseinfo[a].chan);
	 while (strlen(temp_str) < 8) strcat(temp_str," ");
	 build_name_chanid(temp_str,mhp->gseinfo[a].name,mhp->gseinfo[a].chan);
	 while (strlen(mhp->gseinfo[a].name) < 5) strcat(mhp->gseinfo[a].name," ");
	 while (strlen(mhp->gseinfo[a].chan) < 3) strcat(mhp->gseinfo[a].chan," ");

         /* Reading the second header line of the waveform id   */
	 pos += 79;
	 while(g_ptr[pos] != '\n') 
	 {
	    pos++;
	    if ((pos+81) >= gsize)
	    {
	       fprintf(stderr,"Error in module ReadGseToMemory: Inappropriate input file\n");
	       if (nstat > 0) 
	       {
	          free(stn);
	          stn = NULL;
	       }
	       return(0);
	    }
	 }
	 pos++;

         strncpy(temp_str,g_ptr+pos,9);
	 temp_str[9] = 0;
         mhp->gseinfo[a].gain = (float) atof(temp_str);
         strncpy(temp_str,g_ptr+pos+9,1);
         temp_str[1] = '\0';
         mhp->gseinfo[a].units = (short) atoi(temp_str);
         strncpy(temp_str,g_ptr+pos+10,7);
         temp_str[7] = '\0';
         mhp->gseinfo[a].cperiod = (float) atof(temp_str);
         strncpy(temp_str,g_ptr+pos+18,9);
         temp_str[9] = '\0';
         mhp->gseinfo[a].lat = (float) atof(temp_str);
         strncpy(temp_str,g_ptr+pos+28,9);
         temp_str[9] = '\0';
         mhp->gseinfo[a].lon = (float) atof(temp_str);
         strncpy(temp_str,g_ptr+pos+38,9);
         temp_str[9] = '\0';
         mhp->gseinfo[a].alt = (float) atof(temp_str);
         strncpy(temp_str,g_ptr+pos+48,9);
         temp_str[9] = '\0';
         mhp->gseinfo[a].depsens = (float) atof(temp_str);
      }
      else
      {
	 /* if here, gse type WID2 section is found */
	 ngse2++;
	 if ((pos+106) >= gsize)
	 {
	    fprintf(stderr,"Error in module ReadGseToMemory: Inappropriate input file\n");
	    if (nstat > 0) 
	    {
	       free(stn);
	       stn = NULL;
	    }
	    return(0);
	 }
	 strncpy(line,g_ptr+pos,105);
	 line[105] = 0;
         sscanf(line+5,"%hd/%hd/%hd",&mhp->gseinfo[a].tfsamp.year,
		                     &mhp->gseinfo[a].tfsamp.month,
				     &mhp->gseinfo[a].tfsamp.day);
         sscanf(line+16,"%hd:%hd:%f",&mhp->gseinfo[a].tfsamp.hour,
		                     &mhp->gseinfo[a].tfsamp.min,
				     &mhp->gseinfo[a].tfsamp.fsec);
         sscanf(line+29,"%5s",mhp->gseinfo[a].name);
         sscanf(line+35,"%3s",mhp->gseinfo[a].chan);
         sscanf(line+39,"%4s",mhp->gseinfo[a].chid);
         sscanf(line+44,"%3s",mhp->gseinfo[a].form);
	 fok = 0;
	 if (!strncmp(mhp->gseinfo[a].form,"CM6",3))
	 {
	    fok = 1;
	    strcpy(mhp->gseinfo[a].form,"CMP6");
	    mhp->gseinfo[a].dflag = 2;
	 }
	 if (!strncmp(mhp->gseinfo[a].form,"INT",3))
	 {
	    fok = 1;
	    strcpy(mhp->gseinfo[a].form,"INTV");
	    mhp->gseinfo[a].dflag = 0;
	 }
	 if (!fok)
	 {
	    fprintf(stderr,"Unsupported data type %s\n",mhp->gseinfo[a].form);
	    if (nstat > 0) 
	    {
	       free(stn);
	       stn = NULL;
	    }
	    return(0);
	 }
         sscanf(line+48,"%8ld",&(mhp->gseinfo[a].nosamp));
         sscanf(line+57,"%f",&mhp->gseinfo[a].sra);
         sscanf(line+69,"%f",&mhp->gseinfo[a].gain);
         sscanf(line+80,"%f",&mhp->gseinfo[a].cperiod);
         sscanf(line+88,"%6s",mhp->gseinfo[a].systype);
         while (strlen(mhp->gseinfo[a].name) < 5)
         {
	    strcat(mhp->gseinfo[a].name," ");
         }
         while (strlen(mhp->gseinfo[a].chan) < 3)
         {
	    strcat(mhp->gseinfo[a].chan," ");
         }       
	 mhp->gseinfo[a].name[5] = '\0';
	 mhp->gseinfo[a].chan[3] = '\0';
         mhp->gseinfo[a].chid[4] = '\0';
         mhp->gseinfo[a].form[4] = '\0';
         mhp->gseinfo[a].systype[6] = '\0';
	 
         /* The station name and the channel identifier are formatted and
	    then concatenated to the waveform index */
         strcpy(mhp->wfindex[a],mhp->gseinfo[a].name);
	 
         while(strlen(mhp->wfindex[a]) < 5) strcat((*mhp).wfindex[a]," ");
         strncat(mhp->wfindex[a],mhp->gseinfo[a].chan,2);
       
         mhp->wfindex[a][8] = '\0';
         while(strlen(mhp->gseinfo[a].systype) < 6) strcat(mhp->gseinfo[a].systype," ");
	 
	 /* Reading coordinates from station list */
	 if (ngse2 == 1)
	 {
	    nstat = 0;
	    stn = NULL;
            stn = GetStationList(stn,stnfile,&nstat);
	    if (stn == NULL) return(0);
	 }
	 for(b=0;b<nstat;b++)
	 {
	    if (!strncmp(mhp->gseinfo[a].chan,"A",1)) 
	    {
	       b = nstat;
	       break;
	    }
	    if (!strncmp(mhp->gseinfo[a].name,stn[b].name,5)) break;
	 }
	 if (b < nstat)
	 {
	    mhp->gseinfo[a].lat = stn[b].lat;
	    mhp->gseinfo[a].lon = stn[b].lon;
	    mhp->gseinfo[a].alt = stn[b].elev;
	 }
	 else
	 {
	    if (strncmp(mhp->gseinfo[a].chan,"A",1))
	    {
	       fprintf(stderr,"Warning! no coordinates available for station %s\n",mhp->gseinfo[a].name);
	    }
	 }
      } /* else GSE2 */
      
      /* Converting name and chan to capital letters */
      c = strlen(mhp->gseinfo[a].name);
      for (b=0;b<c;b++) mhp->gseinfo[a].name[b] = toupper(mhp->gseinfo[a].name[b]);
      c = strlen(mhp->gseinfo[a].chan);
      for (b=0;b<c;b++) mhp->gseinfo[a].chan[b] = toupper(mhp->gseinfo[a].chan[b]);

      /* Determing the julian minute    */
      mhp->gseinfo[a].tfsamp.jmin = juliam(&(mhp->gseinfo[a].tfsamp.year),
					   &(mhp->gseinfo[a].tfsamp.month),
					   &(mhp->gseinfo[a].tfsamp.day),
					   &(mhp->gseinfo[a].tfsamp.hour),
					   &(mhp->gseinfo[a].tfsamp.min));

      /* Check, if the station name is already present */
      if (a == 1)
      {
	 mhp->nst ++;
      }
      else
      {
	 old_station = 0;
         for (b=1;b<a;b++)
         {
            if (!strcmp(mhp->gseinfo[a].name,mhp->gseinfo[b].name))
            {
	        old_station = 1;  
		break;
            }
         }
	 if (!old_station)
	 {
	    mhp->nst ++;
	 }
      }
      mhp->gseinfo[a].stat_nr = mhp->nst;
       
      /* Moving to the beginning of the data section */
      while((memcmp(&g_ptr[pos],"DAT1",4)) && (memcmp(&g_ptr[pos],"DAT2",4)))
      {
	 pos++;
	 if ((pos+5) >= gsize)
	 {
	    fprintf(stderr,"Error in module ReadGseToMemory: Inappropriate input file\n");
	    if (nstat > 0) 
	    {
	       free(stn);
	       stn = NULL;
	    }
	    return(0);
	 }
      }
      while(g_ptr[pos] != '\n')
      {
	 pos++;
	 if ((pos+1) >= gsize)
	 {
	    fprintf(stderr,"Error in module ReadGseToMemory: Inappropriate input file\n");
	    if (nstat > 0) 
	    {
	       free(stn);
	       stn = NULL;
	    }
	    return(0);
	 }
      }
      pos++;
       
      /* Determing the number of 512 byte blocks of the datas and
	  calculating the amount of dummy values to fill the last
	  block up to 512 bytes and reading data */
      if (strncmp(mhp->gseinfo[a].form,"CMP6",4))  /* if not compressed */
      {
         mhp->gseinfo[a].nblock = mhp->gseinfo[a].nosamp/128 + 1;     
	 mhp->gseinfo[a].dummyval = 128 - (mhp->gseinfo[a].nosamp % 128);
	 if (mhp->gseinfo[a].dummyval == 0) mhp->gseinfo[a].nblock--;
	 mhp->gseinfo[a].cmpbytes = 0;
	 
	 mhp->gseinfo[a].trace_i4 = (int *)malloc(sizeof(int)*mhp->gseinfo[a].nosamp);
	 if (mhp->gseinfo[a].trace_i4 == NULL)
	 {
	    fprintf(stderr,"Memory allocation error in ReadGseFromMemory (5)\n");
	    exit(1);
	 }
	 for (b=0;b<mhp->gseinfo[a].nosamp;b++)
	 {
	    pos = i4gets(&mhp->gseinfo[a].trace_i4[b],g_ptr,pos);
            if (pos == -1) pos = gsize;
	    if ((pos+1) >= gsize)
	    {
	       if (nstat > 0) 
	       {
	          free(stn);
	          stn = NULL;
	       }
	       fprintf(stderr,"Error in module ReadGseToMemory: Inappropriate input file\n");
	       return(0);
	    }
	 }
         while((memcmp(&g_ptr[pos],"CHK1 ",5)) && (memcmp(&g_ptr[pos],"CHK2 ",5)))
	 {
	    pos++;
	    if ((pos+5) >= gsize)
	    {
	       fprintf(stderr,"Error in module ReadGseToMemory: Inappropriate input file\n");
	       if (nstat > 0) 
	       {
	          free(stn);
	          stn = NULL;
	       }
	       return(0);
	    }
	 }
      }
      else   /* if compressed  */
      {
	 b = 0;
	 mhp->gseinfo[a].cbuf = malloc(sizeof(char));
	 if (mhp->gseinfo[a].cbuf == NULL)
	 {
	    fprintf(stderr,"Memory allocation error in ReadGseFromMemory (6)\n");
	    exit(1);
	 }
         while((memcmp(&g_ptr[pos],"CHK1 ",5)) && (memcmp(&g_ptr[pos],"CHK2 ",5)))
         {
	    c = 0;
	    while (g_ptr[pos+c] != '\n') c++;
	    mhp->gseinfo[a].cbuf = realloc(mhp->gseinfo[a].cbuf,sizeof(char)*(b+c));
	    if (mhp->gseinfo[a].cbuf == NULL)
	    {
	       fprintf(stderr,"Memory allocation error in ReadGseFromMemory (7)!\n");
	       exit(1);
	    }
	    memcpy(&mhp->gseinfo[a].cbuf[b],&g_ptr[pos],sizeof(char)*c);
	    b += c;
	    pos += c+1;
	    if ((pos+1) >= gsize)
	    {
	       fprintf(stderr,"Error in module ReadGseToMemory: Inappropriate input file\n");
	       if (nstat > 0) 
	       {
	          free(stn);
	          stn = NULL;
	       }
	       return(0);
	    }
	 }
	 mhp->gseinfo[a].cmpbytes = b;
	 mhp->gseinfo[a].nblock = mhp->gseinfo[a].cmpbytes / 512 + 1; 
	 mhp->gseinfo[a].dummyval = 512 - (mhp->gseinfo[a].cmpbytes % 512);
	 if (mhp->gseinfo[a].dummyval == 0) mhp->gseinfo[a].nblock--;
      }
      
      /* Reading checksum */
      pos += 5;
      pos = i4gets(&mhp->gseinfo[a].checksum,g_ptr,pos);
      if (pos == -1) 
      {
	 if (nstat > 0) 
	 {
	    free(stn);
	    stn = NULL;
	 }
	 return(0);
      }
      mhp->nwf = a;
   }    /* while(1) */
}     /* ReadGseToMemory */


int FreeMasterHeader(mhp)
masterheader *mhp;   
{
   int a,b;
   for(a=1;a<=mhp->nwf;a++)
   {
      if (mhp->gseinfo[a].cbuf != NULL) 
      {
	 free(mhp->gseinfo[a].cbuf);
	 mhp->gseinfo[a].cbuf = NULL;
      }
      if (mhp->gseinfo[a].trace_i4 != NULL) 
      {
	 free(mhp->gseinfo[a].trace_i4);
	 mhp->gseinfo[a].trace_i4 = NULL;
      }
      if (mhp->gseinfo[a].cal_data != NULL)
      {
         if (mhp->gseinfo[a].cal_data[0].realpole != NULL) 
         {
	    free(mhp->gseinfo[a].cal_data[0].realpole); 
	    mhp->gseinfo[a].cal_data[0].realpole = NULL;
         }
         if (mhp->gseinfo[a].cal_data[0].imagpole != NULL) 
         {
	    free(mhp->gseinfo[a].cal_data[0].imagpole); 
	    mhp->gseinfo[a].cal_data[0].imagpole = NULL;
         }
         if (mhp->gseinfo[a].cal_data[0].realzero != NULL) 
         {
	    free(mhp->gseinfo[a].cal_data[0].realzero); 
	    mhp->gseinfo[a].cal_data[0].realzero = NULL;
         }
         if (mhp->gseinfo[a].cal_data[0].imagzero != NULL) 
         {
	    free(mhp->gseinfo[a].cal_data[0].imagzero); 
	    mhp->gseinfo[a].cal_data[0].imagzero = NULL;
         }
	 free(mhp->gseinfo[a].cal_data);
	 mhp->gseinfo[a].cal_data = NULL;
      }
   }
   free(mhp->gseinfo); mhp->gseinfo = NULL;
   free(mhp->wfindex); mhp->wfindex = NULL;
   
   if (mhp->hydainfo.nrofreadphases > 0) 
   {
      free(mhp->hydainfo.phaselist);
      mhp->hydainfo.phaselist = NULL;
      if (mhp->hydainfo.new_stat != NULL)
      {
	 free(mhp->hydainfo.new_stat);
	 mhp->hydainfo.new_stat = NULL;
      }
   }
         
   return(1);
}



int i2gets(val,c_ptr,pos)
short int *val;
char *c_ptr;
int pos;
{
   char tmp[12];
   int a = 0;
   
   memset(tmp,0,12);
   while(isspace((int)c_ptr[pos])) pos++;
   while(!isspace((int)c_ptr[pos])) 
   {
      tmp[a] = c_ptr[pos];
      pos++;
      a++;
      if (a == 12)
      {
	 tmp[11] = 0;
	 fprintf(stderr,"Error in i2gets: Attempting to read a number that is too long\n");
	 fprintf(stderr,"Characters read: %s\n",tmp);
	 return(-1);
      }
   }
   tmp[a] = 0;
   *val = (short)atoi(tmp);
   return(pos);
}

int i4gets(val,c_ptr,pos)
int *val;
char *c_ptr;
int pos;
{
   char tmp[12];
   int a = 0;
   
   memset(tmp,0,12);
   while(isspace((int)c_ptr[pos])) pos++;
   while(!isspace((int)c_ptr[pos])) 
   {
      tmp[a] = c_ptr[pos];
      pos++;
      a++;
      if (a == 12)
      {
	 tmp[11] = 0;
	 fprintf(stderr,"Error in i4gets: Attempting to read a number that is too long\n");
	 fprintf(stderr,"Characters read: %s\n",tmp);
	 return(-1);
      }
   }
   tmp[a] = 0;
   *val = atoi(tmp);
   return(pos);
}


station_struct *ReadStnInfo(stn,nst,stnfile,sedpath,calfilename,jmin)
station_struct *stn;   
int *nst;
char *stnfile;   
char *sedpath;
char *calfilename;
int jmin;
{
   int a,b;
   cal_struct *cal_info;
   int ncal;
   char hcr;
   char tchan[6];
   
   if (stn == NULL) stn = GetStationList(stn,stnfile,nst);
   
   cal_info = NULL;
   ncal = 0;
   if (strlen(sedpath))
   {
      cal_info = ReadSedCalPath(cal_info,sedpath,&ncal);
   }
   if (strlen(calfilename))
   {
      cal_info = ReadCalFromFile(cal_info,calfilename,&ncal);
   }
   
#if 0 /* Test output */   
   for(a=0;a<ncal;a++)
   {
      printf("name %s  chan %s\n",cal_info[a].name,cal_info[a].chan);
      printf("seismotype %s instrype %s  auxid %s\n",cal_info[a].seismotype,cal_info[a].instrtype,cal_info[a].auxid);
      printf("jminb %d  jmine %d\n",cal_info[a].jminb,cal_info[a].jmine);
      printf("resptype %s  npol %d  nzero %d\n",cal_info[a].resptype,cal_info[a].npoles,cal_info[a].nzeros);
      /*
      printf("Poles\n");
      for (b=0;b<cal_info[a].npoles;b++) printf("%e  %e\n",cal_info[a].realpole[b],cal_info[a].imagpole[b]);
      printf("Zeros\n");
      for (b=0;b<cal_info[a].nzeros;b++) printf("%e  %e\n",cal_info[a].realzero[b],cal_info[a].imagzero[b]);
      */
      printf("scalef %e  gain %6.0f  cpv %6.0f  sra %6.0f  units %c\n",cal_info[a].scalefactor,cal_info[a].gain,cal_info[a].cpv,cal_info[a].sra,cal_info[a].units);
      printf("\ncoil %f  damp %f  period %f\n",cal_info[a].seismocoil,cal_info[a].seismodamp,cal_info[a].seismoperiod);
      fflush(stdin);
      scanf("%c",&hcr);
   }
#endif
   
   for (a=0;a<(*nst);a++)
   {
      stn[a].chan = NULL;
      strcpy(tchan,"SRZ ");
      stn[a].chan = AssignCalInfo(stn[a].name,tchan,&stn[a].nchan,cal_info,ncal,jmin,stn[a].chan);
      if (stn[a].nchan > 0)
      {
	 for(b=0;b<stn[a].nchan;b++) printf("%s %s\n",stn[a].chan[b].name,stn[a].chan[b].chan);
      }
   }
   cal_info = FreeCalInfo(cal_info,ncal);
   
   return(stn);
}   

station_struct *GetStationList(stn,stnfile,nst)
station_struct *stn;   
char *stnfile;
int *nst;
{
   FILE *inf;
   int a;
   
   *nst = 0;
   inf = fopen(stnfile,"r");
   if (inf == NULL)
   {
      strcpy(stnfile,"/signals/public/trt_tab/stn_struct.bin");
      inf = fopen(stnfile,"r");
   }
   if (inf != NULL)
   {
      fread(nst,sizeof(int),1,inf);
      stn = (station_struct *)malloc(sizeof(station_struct)*(*nst));
      if (stn == NULL)
      {
         printf("Memory allocation error in GetStationList!\n");
         exit(1);
      }
      a = fread(stn,sizeof(station_struct),(*nst),inf);
      fclose(inf);
      if (a != (*nst))
      {
         fprintf(stderr,"Error while reading station file %s",stnfile);
         free(stn);
	 *nst = 0;
         stn = NULL;
      }
   }
   return(stn);
}


cal_struct *ReadSedCalPath(cal_info,sedpath,ncal)
cal_struct *cal_info;   
char *sedpath;
int *ncal;
{   
   /* Read stndef file */
   cal_info = ReadStnDef(cal_info,sedpath,ncal);
   if (cal_info == NULL) return(NULL);
   
   /* Read auxchanpaz */
   ReadSedPaz(sedpath,cal_info,(*ncal),"aux");
   
   /* Read instrpaz */
   ReadSedPaz(sedpath,cal_info,(*ncal),"instr");
      
   /* That's it */
   return(cal_info);
}

void ReadSedPaz(sedpath,cal_info,ncal,task)
char *sedpath;
cal_struct *cal_info;
int ncal;
char *task;
{
   int a,b;
   FILE *inf;
   char path[256];
   char filename[256];
   char line[128];
   int newpol,newzero;
   char tunits;
   float tscal;
   
   strcpy(path,sedpath);
   if (path[strlen(path)-1] != '/') strcat(path,"/");
   if (!strcmp(task,"aux")) strcat(path,"auxchanpaz/");
   if (!strcmp(task,"instr")) strcat(path,"instrpaz/");
   for (a=0;a<ncal;a++)
   {
      if (strlen(cal_info[a].auxid))
      {
         strcpy(filename,path);
         if (!strcmp(task,"aux")) strcat(filename,cal_info[a].auxid);
         if (!strcmp(task,"instr")) strcat(filename,cal_info[a].seismotype);
         inf = fopen(filename,"r");
         if (inf == NULL)
         {
	    fprintf(stderr,"Unable to open %s\n",filename);
         }
         else
         {
	    while(!feof(inf))
	    {
	       memset(line,0,128);
	       fgets(line,128,inf);
	       if (!strncmp(line,"( Sconst",8)) sscanf(line+1,"%*s %f",&cal_info[a].seismocoil);
	       if (!strncmp(line,"( Period",8)) sscanf(line+1,"%*s %f",&cal_info[a].seismoperiod);
	       if (!strncmp(line,"( Damping",9)) sscanf(line+1,"%*s %f",&cal_info[a].seismodamp);
	       if (!strncmp(line,"DIG2",4))
	       {
		  sscanf(line,"%*s %*d %f %f",&cal_info[a].cpv,&cal_info[a].sra);
	       }
	       if (!strncmp(line,"PAZ2",4))
	       {
		  /* Read header lines */
		  strcpy(cal_info[a].resptype,"PAZ");
                  sscanf(line,"%*s %*d %c %f",&tunits,&tscal);
                  /* sscanf(line,"%*s %*d %c %f",&cal_info[a].units,&cal_info[a].scalefactor); */
                  sscanf(line+40,"%d %d",&newpol,&newzero);
		  
		  if ((cal_info[a].npoles == 0) && (cal_info[a].nzeros == 0))
		  {
		     cal_info[a].units = tunits;
		     cal_info[a].scalefactor = tscal;
		     if (tunits == 'V') 
		     {
			if (cal_info[a].cpv != 0.0) 
			{
			   cal_info[a].scalefactor = tscal * cal_info[a].cpv;
			   cal_info[a].units = 'C';
			}
			else
			{
			   cal_info[a].scalefactor = 1.0;
			   cal_info[a].units = 'C';
			}
		     }	   
		  }
		  else
		  {
		     if (tunits == 'V') 
		     {
			cal_info[a].scalefactor *= tscal * cal_info[a].cpv;
			cal_info[a].units = 'C';
		     }
		     else
		     {
			cal_info[a].scalefactor *= tscal;
			cal_info[a].units = 'C';
		     }
		  }
		  
		  /* Allocate Memory */
		  if (cal_info[a].npoles == 0)
		  {
		     cal_info[a].realpole = (float *)malloc(sizeof(float)*newpol);
		     cal_info[a].imagpole = (float *)malloc(sizeof(float)*newpol);
		  }
		  else
		  {
		     cal_info[a].realpole = (float *)realloc(cal_info[a].realpole,
							     sizeof(float)*(newpol+cal_info[a].npoles));
		     cal_info[a].imagpole = (float *)realloc(cal_info[a].imagpole,
							     sizeof(float)*(newpol+cal_info[a].npoles));
		  }
		  if (cal_info[a].nzeros == 0)
		  {
		     cal_info[a].realzero = (float *)malloc(sizeof(float)*newzero);
		     cal_info[a].imagzero = (float *)malloc(sizeof(float)*newzero);
		  }
		  else
		  {
		     cal_info[a].realzero = (float *)realloc(cal_info[a].realzero,
							     sizeof(float)*(newzero+cal_info[a].nzeros));
		     cal_info[a].imagzero = (float *)realloc(cal_info[a].imagzero,
							     sizeof(float)*(newzero+cal_info[a].nzeros));
		  }
		     
		  if ((cal_info[a].realpole == NULL) || (cal_info[a].imagpole == NULL) ||
		      (cal_info[a].realzero == NULL) || (cal_info[a].imagzero == NULL))
		  {
		     fprintf(stderr,"\nMemory allocation error in ReadSedPaz\n");
		     exit(1);
		  }
		  
		  /* Read poles and zeros */
		  for (b=cal_info[a].npoles;b<cal_info[a].npoles+newpol;b++)
		  {
		       fgets(line,128,inf);
		       sscanf(line,"%f %f",&cal_info[a].realpole[b],&cal_info[a].imagpole[b]);
		  }
		  for (b=cal_info[a].nzeros;b<cal_info[a].nzeros+newzero;b++)
		  {
		       fgets(line,128,inf);
		       sscanf(line,"%f %f",&cal_info[a].realzero[b],&cal_info[a].imagzero[b]);
		  }
		  cal_info[a].npoles += newpol;
		  cal_info[a].nzeros += newzero;
	       }
	    }
	    fclose(inf);
         }
      }
   }
}

cal_struct *ReadStnDef(cal_info,sedpath,ncal)
cal_struct *cal_info;
char *sedpath;
int *ncal;
{
   char filename[256];
   int a;
   int nitems;
   short int year,month,day,hour,min;
   char t1[16],t2[16];
   char line[128];
   FILE *inf;
   
   /* Initial Allocation */
   cal_info = (cal_struct *)malloc(sizeof(cal_struct));
   if (cal_info == NULL)
   {
      fprintf(stderr,"Memory allocation error in ReadStnDef (1)!\n");
      exit(1);
   }
   
   strcpy(filename,sedpath);
   strcat(filename,"stndef/all.def");
   inf = fopen(filename,"r");
   if (inf == NULL)
   {
      fprintf(stderr,"Unable to open %s\n",filename);
      return(NULL);
   }
   a = 0;
   while(!feof(inf))
   {
      fgets(line,128,inf);
      if (line[0] != '#')
      {
         cal_info = (cal_struct *)realloc(cal_info,sizeof(cal_struct)*(a+1));
         if (cal_info == NULL)
         {
            fprintf(stderr,"Memory allocation error ReadStnDef (2)!\n");
            exit(1);
	 }
	 memset(&cal_info[a],0,sizeof(cal_struct));
	 nitems = sscanf(line,"%s %s %s %s %s %s %f %s",cal_info[a].name,
			                                cal_info[a].chan,
						        t1,t2,
						        cal_info[a].seismotype,
						        cal_info[a].instrtype,
						        &cal_info[a].gain,
						        cal_info[a].auxid);
	 if (nitems == 8)
	 {
	    while(strlen(cal_info[a].name) < 5) strcat(cal_info[a].name," ");
	    while(strlen(cal_info[a].chan) < 3) strcat(cal_info[a].name," ");
	    sscanf(t1,"%4hd%2hd%2hd%2hd%2hd",&year,&month,&day,&hour,&min);
	    cal_info[a].jminb = juliam(&year,&month,&day,&hour,&min);
	    if (!strncmp(t2,"----",4))
	    {
	       year = 2100;
	       month = 12;
	       day = 31;
	       hour = 23;
	       min = 59;
	    }
	    else
	    {
 	       sscanf(t2,"%4hd%2hd%2hd%2hd%2hd",&year,&month,&day,&hour,&min);
	    }
	    cal_info[a].jmine = juliam(&year,&month,&day,&hour,&min);
	    a++;
	 }
      }
   }
   fclose(inf);
   *ncal = a;
   return(cal_info);
}

cal_struct *ReadCalFromFile(cal_info,calfilename,ncal)
cal_struct *cal_info;   
char *calfilename;
int *ncal;
{
   
   int a,b,c,d;
   FILE *inf;
   char line[128],line_rest[128];
   short int byear,bmonth,bday,bhour,bmin;
   short int eyear,emonth,eday,ehour,emin;
   char descr[16],tmp_str[16];
   int npoles=0, nzeros=0 ;
   int np,nz;
   float scalefactor=1.0;
      
   /* Initial Allocation */
   if (cal_info == NULL)
   {
      cal_info = (cal_struct *)malloc(sizeof(cal_struct));
      if (cal_info == NULL)
      {
         fprintf(stderr,"Memory allocation error ReadCalFromFile (1) !\n");
         exit(1);
      }
   }
   
   inf = fopen(calfilename,"r");
   if (inf == NULL)
   {
      fprintf(stderr,"Unable to open %s\n",calfilename);
      return(NULL);
   }
   a = (*ncal);
   memset(line,0,sizeof(char)*128);
   while(!feof(inf))
   {
      if ((!strncmp(line,"CAL1 ",5)) || (!strncmp(line,"CAL2 ",5)))
      {
	 /* Reallocate memory */
         cal_info = (cal_struct *)realloc(cal_info,sizeof(cal_struct)*(a+1));
         if (cal_info == NULL)
         {
            fprintf(stderr,"Memory allocation error a is %d!\n",a);
            exit(1);
	 }
	 memset(&cal_info[a],0,sizeof(cal_struct));
	 
	 /* Read calibration section */
	 if (line[3] == '1')
	 {
	    eyear=300;
	    emonth= 0;
	    eday= 0;
	    ehour= 0;
	    emin= 0;
	    sscanf(line+4,"%s %s %s %s %s %2hd%2hd%2hd %2hd%2hd %2hd%2hd%2hd %2hd%2hd",
		   cal_info[a].name,cal_info[a].instrtype,cal_info[a].chan,cal_info[a].seismotype,
		   cal_info[a].resptype,&byear,&bmonth,&bday,&bhour,&bmin,
		   &eyear,&emonth,&eday,&ehour,&emin);	    
	    strcpy(tmp_str,cal_info[a].name);
	    while (strlen(tmp_str) < 6) strcat(tmp_str," ");
	    
	    /* MLR stuff */
            if ((!strncmp(cal_info[a].instrtype,"MLR2",4))  &&
                 (strncmp(cal_info[a].instrtype+6,"NG",2)))
            {
               strncpy(tmp_str+6,cal_info[a].instrtype+6,2);
               if (!strncmp(tmp_str+6,"TC",2))
               {
                  strncpy(cal_info[a].name+4,"TC",2);
                  tmp_str[6] = 0;
                  tmp_str[7] = 0;
	       }
               if (!strncmp(tmp_str+6,"HG",2))
               {
                  strncpy(tmp_str+6,"SZ",2);
	       }
               if (!strncmp(tmp_str+6,"MG",2))
               {
                  strncpy(tmp_str+6,"SM",2);
	       }
               if (!strncmp(tmp_str+6,"LG",2))
               {
                  strncpy(tmp_str+6,"SS",2);
	       }
               if (!strncmp(tmp_str+6,"GR",2))
               {
                  strncpy(tmp_str+6,"SG",2);
	       }
	       tmp_str[8] = 0;
	    }
 		
            /* All strong motion channels */
            if (!strncmp(cal_info[a].instrtype,"LG",2))
            {
               strncpy(tmp_str+6,"SS",2);
	       tmp_str[8] = 0;
	    }
	    
	    if (strlen(tmp_str) < 8) strcat(tmp_str,cal_info[a].chan);
	    while (strlen(tmp_str) < 8) strcat(tmp_str," ");
	    build_name_chanid(tmp_str,cal_info[a].name,cal_info[a].chan);
	    while (strlen(cal_info[a].name) < 5) strcat(cal_info[a].name," ");
	    while (strlen(cal_info[a].chan) < 3) strcat(cal_info[a].chan," ");
	    byear += 1900;
	    eyear += 1900;
	    cal_info[a].jminb = juliam(&byear,&bmonth,&bday,&bhour,&bmin);
	    cal_info[a].jmine = juliam(&eyear,&emonth,&eday,&ehour,&emin);
	    
	    /* Reading poles */
	    fgets(line,128,inf);
	    sscanf(line,"%d",&cal_info[a].npoles);
	    cal_info[a].realpole = (float *)malloc(sizeof(float)*cal_info[a].npoles);
            if (cal_info[a].realpole == NULL)
            {
	       fprintf(stderr,"Memory allocation error  ReadCalFromFile (3) \n");
               exit(1);
	    }
	    cal_info[a].imagpole = (float *)malloc(sizeof(float)*cal_info[a].npoles);
            if (cal_info[a].imagpole == NULL)
            {
	       fprintf(stderr,"\nMemory allocation error ReadCalFromFile (4) \n");
               exit(1);
	    }
	    
           for (b=0;b<cal_info[a].npoles;b++)
            {
               fgets(line,128,inf);
               sscanf(line,"%8f%8f",&(cal_info[a].realpole[b]),&(cal_info[a].imagpole[b]));
	    }
	    
	    /* Reading zeros */
 	    fgets(line,128,inf);
	    sscanf(line,"%d",&cal_info[a].nzeros);
	    cal_info[a].realzero = (float *)malloc(sizeof(float)*cal_info[a].nzeros);
            if (cal_info[a].realzero == NULL)
            {
	       fprintf(stderr,"\nMemory allocation error ReadCalFromFile (5) \n");
               exit(1);
	    }
	    cal_info[a].imagzero = (float *)malloc(sizeof(float)*cal_info[a].nzeros);
            if (cal_info[a].imagzero == NULL)
            {
	       fprintf(stderr,"\nMemory allocation error ReadCalFromFile (6) \n");
               exit(1);
	    }
            for (b=0;b<cal_info[a].nzeros;b++)
            {
               fgets(line,128,inf);
               sscanf(line,"%8f%8f",&(cal_info[a].realzero[b]),&(cal_info[a].imagzero[b]));
	    }
	    
            /* Reading the scale factor */
            fgets(line,128,inf);
            sscanf(line,"%f",&(cal_info[a].scalefactor));
            while((strncmp(line,"CAL",3)) && (!feof(inf)))
            {
               if (line[0] == '$')
               {
                  c = 1;
		  d = 0;
                  strcpy(descr,"");
                  while(line[c] != '&')
                  {
                     if (line[c] != ' ')
                     {
                        descr[d] = line[c];
                        d++;
                        descr[d] = 0;
		     }
                     c++;
		  }
                  c++;
                  if (!strncmp(descr,"RS",2)) sscanf(line+c,"%f",&(cal_info[a].resolution));
                  if (!strncmp(descr,"GA",2)) sscanf(line+c,"%f",&(cal_info[a].gain));
                  if (!strncmp(descr,"VO",2)) sscanf(line+c,"%f",&(cal_info[a].volts));
                  if (!strncmp(descr,"SD",2)) sscanf(line+c,"%f",&(cal_info[a].seismodamp));
                  if (!strncmp(descr,"SP",2)) sscanf(line+c,"%f",&(cal_info[a].seismoperiod));
                  if (!strncmp(descr,"SC",2)) 
                  {
                     sscanf(line+c,"%f",&(cal_info[a].seismocoil));
                     cal_info[a].seismocoil /= 100.0;
		  }
	       }   /* if line = $ */
               fgets(line,128,inf);
	    } /* while strncmp CAL */
	    
	    /* Increment cal_info counter */
	    a++;
	 }
	 else /* if here CAL2 is found */
	 {
	    /* Read CAL2 line */
	    eyear=2200;
	    emonth= 0;
	    eday= 0;
	    ehour= 0;
	    emin= 0;
	    sscanf(line+5,"%s %s %s %s %f %f %f %hd/%hd/%hd %hd:%hd %hd/%hd/%hd %hd:%hd",
		   cal_info[a].name,cal_info[a].chan,cal_info[a].auxid,cal_info[a].instrtype,
		   &cal_info[a].widg,&cal_info[a].calper,&cal_info[a].sra,
		   &byear,&bmonth,&bday,&bhour,&bmin,&eyear,&emonth,&eday,&ehour,&emin);
	    while (strlen(cal_info[a].name) < 5) strcat(cal_info[a].name," ");
	    while (strlen(cal_info[a].chan) < 3) strcat(cal_info[a].chan," ");
	    cal_info[a].jminb = juliam(&byear,&bmonth,&bday,&bhour,&bmin);
	    cal_info[a].jmine = juliam(&eyear,&emonth,&eday,&ehour,&emin);
            cal_info[a].scalefactor = 1.0;
	    line[0]= 0;
            while((strncmp(line,"CAL",3)) && (!feof(inf)))
	    {
	       fgets(line,128,inf);
	       if (!strncmp(line,"DIG2",4))
	       {
		  sscanf(line,"%*s %*d %f %f",&cal_info[a].cpv,&cal_info[a].sra);
	       }
	       if (!strncmp(line,"PAZ2",4))
	       {
		  /* Read header lines */
		  strcpy(cal_info[a].resptype,"PAZ");
		  /***********
                  sscanf(line,"%*s %*d %c %f",&cal_info[a].units,&cal_info[a].scalefactor);
                  sscanf(line+40,"%d %d",&cal_info[a].npoles,&cal_info[a].nzeros);
		  ***********/
                  sscanf(line,"%*s %*d %c %f",&cal_info[a].units,&scalefactor);
                  sscanf(line+40,"%d %d",&npoles,&nzeros);
		  cal_info[a].scalefactor *=scalefactor;
		  
		  /* Allocate Memory */
		  np= cal_info[a].npoles;
		  nz= cal_info[a].nzeros;
		  cal_info[a].npoles += npoles;
		  cal_info[a].nzeros += nzeros;
		  if ( cal_info[a].realpole ) {
		     cal_info[a].realpole = (float *)realloc(cal_info[a].realpole,sizeof(float)*cal_info[a].npoles);
		     cal_info[a].imagpole = (float *)realloc(cal_info[a].imagpole,sizeof(float)*cal_info[a].npoles);
		     cal_info[a].realzero = (float *)realloc(cal_info[a].realzero,sizeof(float)*cal_info[a].nzeros);
		     cal_info[a].imagzero = (float *)realloc(cal_info[a].imagzero,sizeof(float)*cal_info[a].nzeros);
		  }
		  else {
		     cal_info[a].realpole = (float *)malloc(sizeof(float)*cal_info[a].npoles);
		     cal_info[a].imagpole = (float *)malloc(sizeof(float)*cal_info[a].npoles);
		     cal_info[a].realzero = (float *)malloc(sizeof(float)*cal_info[a].nzeros);
		     cal_info[a].imagzero = (float *)malloc(sizeof(float)*cal_info[a].nzeros);
		  }
		     
		  if ((cal_info[a].realpole == NULL) || (cal_info[a].imagpole == NULL) ||
		      (cal_info[a].realzero == NULL) || (cal_info[a].imagzero == NULL))
		  {
		     fprintf(stderr,"\nMemory allocation error ReadCalFromFile (7) \n");
		     exit(1);
		  }
		  
		  /* Read poles and zeros */
		  for (b=np;b<cal_info[a].npoles;b++)
		  {
		       fgets(line,128,inf);
		       sscanf(line,"%f %f",&cal_info[a].realpole[b],&cal_info[a].imagpole[b]);
		  }
		  for (b=nz;b<cal_info[a].nzeros;b++)
		  {
		       fgets(line,128,inf);
		       sscanf(line,"%f %f",&cal_info[a].realzero[b],&cal_info[a].imagzero[b]);
		  }
	       } /* if PAZ */
	       if ( !strncmp(line," (",2) || line[0] == '(' ) { 
		  /*interpret comment parameters 
					( Sconst 5.7000   [V/cm/s])
 					( Period 2.00     [sec])
 					( Damping 0.60 )
					( Gain 12345.6 )
 					( CalVolt 0.70    [mV])
					( Resol 24        [bits])
					( MaxVolt 16.0	  [Volts])
 					( Lat 46.789      [deg])
 					( Lon 5.678       [deg])
					( Elev 1.2345     [km])
		  examples of valid parameters */
		  csplitstring(line,line_rest,line,'(');
		  remove_leading(line,' ');
		  csplitstring(line,descr,line,' ');
		  
                  if (!strncmp(descr,"Resol",4))   sscanf(line,"%f",&(cal_info[a].resolution));
                  if (!strncmp(descr,"Gain",4))    sscanf(line,"%f",&(cal_info[a].gain));
                  if (!strncmp(descr,"MaxVolt",7)) sscanf(line,"%f",&(cal_info[a].volts));
                  if (!strncmp(descr,"Damping",7)) sscanf(line,"%f",&(cal_info[a].seismodamp));
                  if (!strncmp(descr,"Period",6))  sscanf(line,"%f",&(cal_info[a].seismoperiod));
                  if (!strncmp(descr,"Lat",3))     sscanf(line,"%f",&(cal_info[a].lat));
                  if (!strncmp(descr,"Lon",3))     sscanf(line,"%f",&(cal_info[a].lon));
                  if (!strncmp(descr,"Elev",4))    sscanf(line,"%f",&(cal_info[a].elev));
                  if (!strncmp(descr,"Sconst",5))  sscanf(line,"%f",&(cal_info[a].seismocoil));

	       }
	    } /* while strncmp CAL */
	    cal_info[a].scalefactor *=cal_info[a].cpv;

	    /* Increment cal_info counter */
	    a++;
	 } /* else (if CAL2 ) */
      } /* if CAL1 || CAL2 */
      else
      {
         fgets(line,128,inf);
      }
   } /* while !feof(inf) */
   fclose(inf);
   *ncal = a;
   
   /* That's it */
   return(cal_info);
}
cal_struct *PointerCalInfo(sname,schan,rsys,nchan,cal_info,ncal,jmin,channels)
char *sname;
char *schan;
char *rsys;
int *nchan;
cal_struct *cal_info;
int ncal;
int jmin;
cal_struct *channels;
{
   int a,b,c;
   int found;
   b = 0;
   for(a=0;a<ncal;a++)
   {
      if (!strncmp(sname,cal_info[a].name,5))
      {
	 /* if here a channel section is found */
	 if ((jmin >= cal_info[a].jminb) && (jmin <= cal_info[a].jmine) 
	     && (rsys[0] == cal_info[a].auxid[0] ))
	 {
	       if (!strncmp(schan,cal_info[a].chan,3))
	       {
		  *nchan = 1;
		  return(&cal_info[a]);
	       }
	 } /* if jmin */
      } /* if strncmp sname */
   } /* for a */
   *nchan = b;
   return(NULL);
}
cal_struct *AssignCalInfo(sname,schan,nchan,cal_info,ncal,jmin,channels)
char *sname;
char *schan;
int *nchan;
cal_struct *cal_info;
int ncal;
int jmin;
cal_struct *channels;
{
   int a,b,c;
   int found;
   
   /* Initial allocation */
   channels = (cal_struct *)malloc(sizeof(cal_struct));
   if (channels == NULL)
   {
      fprintf(stderr,"Memory allocation error AssignCalInfo (1) !\n");
      exit(1);
   }
   memset(channels,0,sizeof(cal_struct));
   b = 0;
   for(a=0;a<ncal;a++)
   {
      if (!strncmp(sname,cal_info[a].name,5))
      {
	 /* if here a channel section is found */
	 if ((jmin >= cal_info[a].jminb) && (jmin <= cal_info[a].jmine))
	 {
	    if (strlen(schan))
	    {
	       if (!strncmp(schan,cal_info[a].chan,3))
	       {
	          memcpy(channels,&cal_info[a],sizeof(cal_struct));
		  channels[0].realpole = (float *)malloc(sizeof(float) * cal_info[a].npoles);
		  channels[0].imagpole = (float *)malloc(sizeof(float) * cal_info[a].npoles);
		  channels[0].realzero = (float *)malloc(sizeof(float) * cal_info[a].nzeros);
		  channels[0].imagzero = (float *)malloc(sizeof(float) * cal_info[a].nzeros);
		  if ((channels[0].realpole == NULL) || (channels[0].imagpole == NULL) || 
		      (channels[0].realzero == NULL) || (channels[0].imagzero == NULL))
		  {
		     fprintf(stderr,"Memory allocation error AssignCalInfo (2) !\n");
		     exit(1);
		  }
		  memcpy(channels[0].realpole,cal_info[a].realpole,cal_info[a].npoles*sizeof(float));
		  memcpy(channels[0].imagpole,cal_info[a].imagpole,cal_info[a].npoles*sizeof(float));
		  memcpy(channels[0].realzero,cal_info[a].realzero,cal_info[a].nzeros*sizeof(float));
		  memcpy(channels[0].imagzero,cal_info[a].imagzero,cal_info[a].nzeros*sizeof(float));
		  *nchan = 1;
		  return(channels);
	       }
	    }
	    else
	    {
	       /* if here a valid channel section is found */
	       /* Now check, if the information is already read */
	       found = 0;
	       for (c=0;c<b;c++)
	       {
	          if (!strcmp(channels[c].chan,cal_info[a].chan) &&
		      channels[c].auxid[0] == cal_info[a].auxid[0] )
	          {
		     found = 1;
		     break;
	          }
	       }
	       if (!found)
	       {
	          /* if here, cal_info section must be added to channels */
	          channels = (cal_struct *)realloc(channels,sizeof(cal_struct)*(b+1));
	          if (channels == NULL)
	          {
		     fprintf(stderr,"Memory allocatio error!\n");
		     exit(1);
	          }
	          memset(&channels[b],0,sizeof(cal_struct));
	          memcpy(&channels[b],&cal_info[a],sizeof(cal_struct));
		  channels[b].realpole = (float *)malloc(sizeof(float) * cal_info[a].npoles);
		  channels[b].imagpole = (float *)malloc(sizeof(float) * cal_info[a].npoles);
		  channels[b].realzero = (float *)malloc(sizeof(float) * cal_info[a].nzeros);
		  channels[b].imagzero = (float *)malloc(sizeof(float) * cal_info[a].nzeros);
		  if ((channels[b].realpole == NULL) || (channels[b].imagpole == NULL) || 
		      (channels[b].realzero == NULL) || (channels[b].imagzero == NULL))
		  {
		     fprintf(stderr,"Memory allocation error AssignCalInfo (3) !\n");
		     exit(1);
		  }
		  memcpy(channels[b].realpole,cal_info[a].realpole,cal_info[a].npoles*sizeof(float));
		  memcpy(channels[b].imagpole,cal_info[a].imagpole,cal_info[a].npoles*sizeof(float));
		  memcpy(channels[b].realzero,cal_info[a].realzero,cal_info[a].nzeros*sizeof(float));
		  memcpy(channels[b].imagzero,cal_info[a].imagzero,cal_info[a].nzeros*sizeof(float));
	          b++;
	       } /* if !found */
	    } /* else */
	 } /* if jmin */
      } /* if strncmp sname */
   } /* for a */
   *nchan = b;
   return(channels);
}

cal_struct *FreeCalInfo(cal_info,ncal)
cal_struct *cal_info;
int ncal;
{
   int a;
   
   if (cal_info != NULL)
   {
      for (a=0;a<ncal;a++)
      {
	 if (cal_info[a].npoles > 0)
	 {
	    free(cal_info[a].realpole); cal_info[a].realpole = NULL;
	    free(cal_info[a].imagpole); cal_info[a].imagpole = NULL;
	 }
	 if (cal_info[a].nzeros > 0)
	 {
	    free(cal_info[a].realzero); cal_info[a].realzero = NULL;
	    free(cal_info[a].imagzero); cal_info[a].imagzero = NULL;
	 }
      }
      free(cal_info); cal_info = NULL;
   }
		 
   return(cal_info);
}


