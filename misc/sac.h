/*
 * Copyright (C) 1999-2010 Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser Public License for more details.

 * You should have received a copy of the GNU Lesser Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */


/*  sac.h

	include file for SAC file functions


	Anthony Lomax
	Dept. of Theoretical Geophysics
	University of Utrecht
	PO Box 80.021
	3508 TA, Utrecht
	The Netherlands


	history:

	ver 1.00  27JUL1993  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

*/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <limits.h>
/*#include <float.h>*/

#ifdef EXTERN_MODE
#   define EXTERN_TXT extern
#else
#   define EXTERN_TXT
#endif


	/* structure defining SAC header
		(see SAC User's Manual, for description) */

#define LOGIC unsigned long

EXTERN_TXT struct sac_hdr {
	float delta; float depmin; float depmax; float scale; float odelta;
	float b; float e; float o; float a; float intr_00;
	float t0; float t1; float t2; float t3; float t4;
	float t5; float t6; float t7; float t8; float t9;
	float f; float resp0; float resp1; float resp2; float resp3;
	float resp4; float resp5; float resp6; float resp7; float resp8;
	float resp9; float stla; float stlo; float stel; float stdp;
	float evla; float evlo; float evel; float evdp; float unu_00;
	float user0; float user1; float user2; float user3; float user4;
	float user5; float user6; float user7; float user8; float user9;
	float dist; float az; float baz; float gcarc; float intr_01;
	float intr_02; float depmen; float cmpaz; float cmpinc; float unu_01;
	float unu_02; float unu_03; float unu_04; float unu_05; float unu_06;
	float unu_07; float unu_08; float unu_09; float unu_10; float unu_11; 
	long nzyear; long nzjday; long nzhour; long nzmin; long nzsec;
	long nzmsec; long nvhdr; long intr_03; long intr_04; long npts;
	long intr_05; long intr_06; long unu_12; long unu_13; long unu_14;
	long iftype; long idep; long iztype; long unu_15; long iinst;
	long istreg; long ievreg; long ievtyp; long iqual; long isynth;
	long unu_16; long unu_17; long unu_18; long unu_19; long unu_20;
	long unu_21; long unu_22; long unu_23; long unu_24; long unu_25; 
	LOGIC leven; LOGIC lpspol; LOGIC lovrok; LOGIC lcalda; LOGIC unu_26;
	char kstnm[8]; char kevnm[16]; 
	char kt0[8]; char kt1[8]; char kt2[8];
	char kt3[8]; char kt4[8]; char kt5[8];
	char kt6[8]; char kt7[8]; char kt8[8];
	char kt9[8]; char kf[8]; char kuser0[8];
	char kuser1[8]; char kuser2[8]; char kcmpnm[8];
	char knetwk[8]; char kdatrd[8]; char kinst[8];
};


/*** function to read sac header  */

EXTERN_TXT int read_sac_hdr(fp_sacfile, phdr)
FILE *fp_sacfile;
struct sac_hdr *phdr;
{

	rewind(fp_sacfile);

	if (fread(phdr, (size_t) sizeof(struct sac_hdr), (size_t) 1,
			fp_sacfile) != 1)
		return(-1);

		/* add null charcter to strings */

	phdr->kstnm[7] = '\0';
	phdr->kevnm[15] = '\0';
	phdr->kt0[7] = '\0';
	phdr->kt1[7] = '\0';
	phdr->kt2[7] = '\0';
	phdr->kt3[7] = '\0';
	phdr->kt4[7] = '\0';
	phdr->kt5[7] = '\0';
	phdr->kt6[7] = '\0';
	phdr->kt7[7] = '\0';
	phdr->kt8[7] = '\0';
	phdr->kt9[7] = '\0';
	phdr->kf[7] = '\0';
	phdr->kuser0[7] = '\0';
	phdr->kuser1[7] = '\0';
	phdr->kuser2[7] = '\0';
	phdr->kcmpnm[7] = '\0';
	phdr->knetwk[7] = '\0';
	phdr->kdatrd[7] = '\0';
	phdr->kinst[7] = '\0';

	return(0);

}


/*** function to dump selected variables from sac header  */

EXTERN_TXT int dump_sac_hdr(phdr)
struct sac_hdr *phdr;
{

	fprintf(OUT_LEVEL_1, "\nDump of SAC file header:\n\n");

	fprintf(OUT_LEVEL_1,
		"delta %lf  depmin %lf  depmax %lf  scale %lf  odelta %lf\n",
		phdr->delta, phdr->depmin, phdr->depmax, phdr->scale,
		phdr->odelta);

	fprintf(OUT_LEVEL_1,
		"stla %lf  stlo %lf  stel %lf  stdp %lf\n",
		phdr->stla, phdr->stlo, phdr->stel, phdr->stdp);

	fprintf(OUT_LEVEL_1,
		"evla %lf  evlo %lf  evel %lf  evdp %lf\n",
		phdr->evla, phdr->evlo, phdr->evel, phdr->evdp);

	fprintf(OUT_LEVEL_1,
		"az %lf  baz %lf  gcarc %lf\n",
		phdr->az, phdr->baz, phdr->gcarc);

	fprintf(OUT_LEVEL_1,
		"depmen %lf  cmpaz %lf  cmpinc %lf\n",
		phdr->depmen, phdr->cmpaz, phdr->cmpinc);

	fprintf(OUT_LEVEL_1,
		"nzyear %d  nzjday %d  nzhour %d  nzmin %d  nzsec %d\n",
		phdr->nzyear, phdr->nzjday, phdr->nzhour, phdr->nzmin,
		phdr->nzsec);

	fprintf(OUT_LEVEL_1,
		"nzmsec %d  nvhdr %d  npts %d\n",
		phdr->nzmsec, phdr->nvhdr, phdr->npts);

	fprintf(OUT_LEVEL_1,
		"iftype %d  idep %d  iztype %d  iinst %d\n",
		phdr->iftype, phdr->idep, phdr->iztype, phdr->iinst);

	fprintf(OUT_LEVEL_1,
		"istreg %d  ievreg %d  ievtyp %d  iqual %d  isynth %d\n",
		phdr->istreg, phdr->ievreg, phdr->ievtyp, phdr->iqual,
		phdr->isynth);

	fprintf(OUT_LEVEL_1,
		"leven %d  lpspol %d  lovrok %d  lcalda %d\n",
		phdr->leven, phdr->lpspol, phdr->lovrok, phdr->lcalda);

	fprintf(OUT_LEVEL_1,
		"kstnm %s  kevnm %s\n", phdr->kstnm, phdr->kevnm);

	fprintf(OUT_LEVEL_1,
		"kt0 %s  kt1 %s  kt2 %s  kt3 %s\n",
		phdr->kt0, phdr->kt1, phdr->kt2, phdr->kt3);

	fprintf(OUT_LEVEL_1,
		"kt4 %s  kt5 %s  kt6 %s  kt7 %s\n",
		phdr->kt4, phdr->kt5, phdr->kt6, phdr->kt7);

	fprintf(OUT_LEVEL_1,
		"kt8 %s  kt9 %s  kf %s  kuser0 %s\n",
		phdr->kt8, phdr->kt9, phdr->kf, phdr->kuser0);

	fprintf(OUT_LEVEL_1,
		"kuser1 %s  kuser2 %s  kuser3 %s\n",
		phdr->kuser1, phdr->kuser2, phdr->kcmpnm);

	fprintf(OUT_LEVEL_1,
		"knetwk %s  kdatrd %s  kinst %s\n",
		phdr->knetwk, phdr->kdatrd, phdr->kinst);

	fprintf(OUT_LEVEL_1,"\n");

	return(0);

}



/*** function to read sac data  */

EXTERN_TXT double *read_sac_data(fp_sacfile, phdr)
FILE *fp_sacfile;
struct sac_hdr *phdr;
{
	float *pdata, *ptmp;
	double *pdata2, *ptmp2;
	int npts_read, n;

	if ((pdata = (float *)
		    malloc((size_t) (sizeof(float) * phdr->npts))) == NULL) {
		fprintf(OUT_LEVEL_0, "ERROR allocating data array.\n");
		return(NULL);
	}

	if ((pdata2 = (double *)
		    malloc((size_t) (sizeof(double) * phdr->npts))) == NULL) {
		fprintf(OUT_LEVEL_0, "ERROR allocating data2 array.\n");
		return(NULL);
	}

	rewind(fp_sacfile);

	if (fseek(fp_sacfile, (long) sizeof(struct sac_hdr), SEEK_SET)) {
		fprintf(OUT_LEVEL_0, "ERROR seeking in sac file.\n");
		return(NULL);
	}

	if ((npts_read =
		    fread(pdata, (size_t) sizeof(float), (size_t) phdr->npts,
		    fp_sacfile)) != phdr->npts) {
		fprintf(OUT_LEVEL_0, "ERROR reading data in sac file.\n");
		fprintf(OUT_LEVEL_0, "  %d/%d points read.\n",
			npts_read, phdr->npts);
		return(NULL);
	}

	fprintf(OUT_LEVEL_0, "SAC file data (0-4): %e %e %e %e %e \n",
		pdata, pdata + 1, pdata + 2, pdata + 3, pdata + 4);
	fprintf(OUT_LEVEL_0, "  %d points read.\n", npts_read);

	ptmp = pdata;
	ptmp2 = pdata2;
	for (n = 0; n < npts_read; n++)
		*(ptmp2++) = (double) *(ptmp++);

	free(pdata);

	return(pdata2);
	
}

/*** function to dump selected variables from sac header  */

EXTERN_TXT int dump_sac_data(phdr, pdata)
struct sac_hdr *phdr;
double *pdata;
{

	int n;

	fprintf(OUT_LEVEL_1, "\nDump of SAC file data:\n\n");

	for (n = 0; n < 100; n++)
		fprintf(OUT_LEVEL_1, "%f  ", *pdata++);

	fprintf(OUT_LEVEL_1, "\n");

	return(0);

}


/*** function to write P file from seis data with sac header  */

EXTERN_TXT int write_p_seis(fn_pfile, phdr, pdata, chrcomm, chrproc)
char *fn_pfile;
struct sac_hdr *phdr;
double *pdata;
char *chrcomm, *chrproc;
{

	FILE *fp_pfile;
	int n;

		/* open P file */

	if ((fp_pfile = fopen(fn_pfile, "w")) == NULL) {
		fprintf(OUT_LEVEL_0, "ERROR opening P file.\n");
		exit(EXIT_ERROR_FILEIO);
	}


		/* write header */

	write_p_hdr(fp_pfile, phdr, chrcomm, chrproc, 0);

		/* write data */

	write_p_data(fp_pfile, phdr, pdata);

	fclose(fp_pfile);

	return(0);

}


/*** function to write seis data in P file  */

EXTERN_TXT int write_p_data(fp_pfile, phdr, pdata)
FILE *fp_pfile;
struct sac_hdr *phdr;
double *pdata;
{

	int n;

		/* write data */

	fprintf(fp_pfile, "1 %f %f %d\n", phdr->b, phdr->delta, phdr->npts);
	for (n = 0; n < phdr->npts; n++)
		fprintf(fp_pfile, "%.6e  ", *pdata++);

	fprintf(fp_pfile, "\n");

	return(0);

}


/*** function to write seis data in P file  */

EXTERN_TXT int write_p_data_binary(fp_pfile, phdr, pdata)
FILE *fp_pfile;
struct sac_hdr *phdr;
double *pdata;
{

	int n;
	float fdata;

/*  WARNING - THIS FUNCTION IS NOT WORKING */

		/* write data */

	fprintf(fp_pfile, "1.0 %.1f %.1f %.1f\n",
		phdr->b, phdr->delta, phdr->npts);
	for (n = 0; n < phdr->npts; n++) {
		fdata = *(pdata++);
		fwrite(&fdata, sizeof(float), 1, fp_pfile);
	}

	fprintf(fp_pfile, "\n");

	return(0);

}


/*** function to write seis header in P file  */

EXTERN_TXT int write_p_hdr(fp_pfile, phdr, chrcomm, chrproc, n_addl_hdr)
FILE *fp_pfile;
struct sac_hdr *phdr;
char *chrcomm, *chrproc;
int n_addl_hdr;
{

		/* write header stuff */

	fprintf(fp_pfile,
		"-997 NHDR: %d\n", 11 + n_addl_hdr);
	fprintf(fp_pfile, "2 0\n");

	fprintf(fp_pfile,
		"-997 ID: %s %s %s\n", phdr->kstnm, phdr->kinst, phdr->kevnm);
	fprintf(fp_pfile, "2 0\n");

	fprintf(fp_pfile,
		"-997 COMP: az %f  inc %f  name %s\n",
		phdr->cmpaz, phdr->cmpinc, phdr->kcmpnm);
	fprintf(fp_pfile, "2 0\n");

	fprintf(fp_pfile,
		"-997 REF_TIME: %4.4d %3.3d %2.2d:%2.2d:%2.2d.%3.3d\n",
		phdr->nzyear, phdr->nzjday, phdr->nzhour, phdr->nzmin,
		phdr->nzsec, phdr->nzmsec);
	fprintf(fp_pfile, "2 0\n");

	fprintf(fp_pfile, "-997 DATA: start %f  delta %f  npts %d\n",
		phdr->b, phdr->delta, phdr->npts);
	fprintf(fp_pfile, "2 0\n");

	fprintf(fp_pfile,
		 "-997 DEPVAR: min %.4e  max %.4e  mean %.4e\n",
		phdr->depmin, phdr->depmax, phdr->depmen);
	fprintf(fp_pfile, "2 0\n");

	fprintf(fp_pfile,
		"-997 STA: lat %.3f  long %.3f  elev %.2f  depth %.2f\n",
		phdr->stla, phdr->stlo, phdr->stel, phdr->stdp);
	fprintf(fp_pfile, "2 0\n");

	fprintf(fp_pfile,
	    "-997 EVENT: lat %.3f  long %.3f  elev %.2f  depth %.2f  time %f\n",
		phdr->evla, phdr->evlo, phdr->evel, phdr->evdp, phdr->o);
	fprintf(fp_pfile, "2 0\n");

	fprintf(fp_pfile,
	    "-997 DIST: dist %.2f  az %.2f  baz %.2f  gcarc %.4f\n",
		phdr->dist, phdr->az, phdr->baz, phdr->gcarc);
	fprintf(fp_pfile, "2 0\n");

	fprintf(fp_pfile,
	    "-997 COMM: %s\n", chrcomm);
	fprintf(fp_pfile, "2 0\n");

	fprintf(fp_pfile,
	    "-997 PROC: %s\n", chrproc);
	fprintf(fp_pfile, "2 0\n");

	return(0);

}
