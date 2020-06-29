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


/*   ingv_list2hyp.c

	Program to convert ISC ISF hypocenter lines to NLLoc .hyp format

*/


/*
	history:

	ver 01    29Dec2003  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

*/



#include "../src/GridLib.h"


/* defines */



/* globals */



/* functions */

int isc_isf2Hyp(int , char **);
int ReadIscIsfLines(FILE *fp_in, HypoDesc *phypo);


/*** program to sum event scatter files */

#define PNAME  "isc_isf2hyp"


int main(int argc, char *argv[])
{

	int istat, narg;


	/* set program name */

	strcpy(prog_name, PNAME);


	/* check command line for correct usage */

	fprintf(stdout, "\n%s Arguments: ", prog_name);
	for (narg = 0; narg < argc; narg++)
		fprintf(stdout, "<%s> ", argv[narg]);
	fprintf(stdout, "\n");

	if (argc < 3) {
		nll_puterr("ERROR wrong number of command line arguments.");
		disp_usage(PNAME,
"<isc_isf_file> <out_hyp_file> [EllLenMax [RMSMax [NRdgsMin [GapMax]]]]");
		exit(-1);
	}

	if ((istat = isc_isf2Hyp(argc, argv)) < 0) {
		nll_puterr("ERROR converting hypoellise summary file.");
		exit(-1);
	}



	exit(0);

}



int isc_isf2Hyp(int argc, char *argv[])
{
	int istat;
	int nLocWritten, nLocRead;
	char fn_hyp_out[FILENAME_MAX];
	char fn_hypoell_in[FILENAME_MAX];
	FILE *fp_hyp_out, *fp_hypoell_in;

	double EllLenMax, RMSMax;
	int NRdgsMin, GapMax;

	GridDesc locgrid;
	HypoDesc Hypo;



	/* get command line parameters */
	strcpy(fn_hypoell_in, argv[1]);
	strcpy(fn_hyp_out, argv[2]);


	EllLenMax = 1.0e6;
	if (argc > 4) {
		sscanf(argv[3], "%lf", &EllLenMax);
	}
	fprintf(stdout, "  Ellipsoid Len3 Maximum: %lf\n", EllLenMax);
	RMSMax = 1.0e6;
	if (argc > 5) {
		sscanf(argv[4], "%lf", &RMSMax);
	}
	fprintf(stdout, "  RMS Maximum: %lf\n", RMSMax);
	NRdgsMin = 0;
	if (argc > 6) {
		sscanf(argv[5], "%d", &NRdgsMin);
	}
	fprintf(stdout, "  Num Readings Minimum: %d\n", NRdgsMin);
	GapMax = 360;
	if (argc > 7) {
		sscanf(argv[6], "%d", &GapMax);
	}
	fprintf(stdout, "  Gap Maximum: %d\n", GapMax);



	/* open isc_isf file */

	if ((fp_hypoell_in = fopen(fn_hypoell_in, "r")) == NULL) {
		nll_puterr("ERROR: opening isc_isf hypocenter input file.");
		return(-1);
	}

	/* open ascii hypocenter output file */
	if ((fp_hyp_out = fopen(fn_hyp_out, "w")) == NULL) {
		nll_puterr("ERROR: opening hypocenter output file.");
		return(-1);
	}



	nLocWritten = 0;
	nLocRead = 0;
	while (1) {

		strcpy(Hypo.fileroot, fn_hypoell_in);
		strcpy(Hypo.locStat, "Converted from isc_isf List");
		Hypo.grid_misfit_max = -1.0;

		/* read next hypocenter */
		if ((istat = ReadIscIsfLines(fp_hypoell_in, &Hypo)) == EOF)
			break;
		else if (istat < 0) {
			nll_puterr2("ERROR: reading isc_isf file", fn_hypoell_in);
			break;
		}
		nLocRead++;

		if (strcmp(Hypo.locStat, "ABORTED") == 0) {
			nll_puterr("WARNING: location ABORTED, ignoring event");
			continue;
		} else if (strcmp(Hypo.locStat, "REJECTED") == 0) {
			nll_puterr("WARNING: location REJECTED, ignoring event");
			continue;
		} else if (Hypo.ellipsoid.len1 > EllLenMax
				|| Hypo.ellipsoid.len2 > EllLenMax
				|| Hypo.ellipsoid.len3 > EllLenMax) {
/*			nll_puterr(
"WARNING: location ellipsoid Len is greater than EllLenMax, ignoring event");*/
			continue;
		} else if (Hypo.rms > RMSMax) {
/*			nll_puterr(
"WARNING: location RMS is Greater than RMSMax, ignoring event");*/
			continue;
		} else if (Hypo.nreadings < NRdgsMin) {
/*			nll_puterr(
"WARNING: location num readings is less than NRdgsMin, ignoring event");*/
			continue;
		} else if (Hypo.gap > GapMax) {
/*			nll_puterr(
"WARNING: location gap is greater than GapMax, ignoring event");*/
			continue;
		} else {
			NumArrivals = 0;
			WriteLocation(fp_hyp_out, &Hypo,
				Arrival, NumArrivals, fn_hyp_out, 0, 0, 1,
				&locgrid, 0);
		}


		/* write ellipsoid */
/*
		phypo = &Hypo;
		fprintf(fp_hyp_out,
			"ELLIPSOID  Hyp  %lf %lf %lf",
			phypo->dlat, phypo->dlong, phypo->depth);
		fprintf(fp_hyp_out,
			" Ell1  %.1lf %.1lf %.2le",
			phypo->ellipsoid.az1, phypo->ellipsoid.dip1,
				phypo->ellipsoid.len1);
		fprintf(fp_hyp_out, " Ell2  %.1lf %.1lf %.2le",
			phypo->ellipsoid.az2, phypo->ellipsoid.dip2,
				phypo->ellipsoid.len2);
		fprintf(fp_hyp_out, " Ell3  %.2le\n", phypo->ellipsoid.len3);
*/

		/* write end line and blank line */
		fprintf(fp_hyp_out, "END_NLLOC\n\n");

		nLocWritten++;

	}


	fclose(fp_hyp_out);

	/* write message */
	fprintf(stdout,
"%d locations read, %d written to ascii hyp file <%s>\n",
		nLocRead, nLocWritten, fn_hyp_out);


	return(0);

}



/*** function to read isc_isf Bulletin List record to HypoDesc structure */

int ReadIscIsfLines(FILE *fp_in, HypoDesc *phypo)
{

	int istat;
	char *cstat;
	double mag;

	static char line[MAXLINE_LONG];


	/* read next line */
	while (1) {
		cstat = fgets(line, MAXLINE_LONG, fp_in);
		if (cstat == NULL)
			return(EOF);
		// check if magnitude line
		if (line[0] == 'm' || line[0] == 'M' )
			continue;
		break;
	}

	/* read hypocenter parameters */
/*
   Date       Time        Err   RMS Latitude Longitude  Smaj  Smin  Az Depth   Err Ndef Nsta Gap  mdist  Mdist Qual   Author      OrigID
1985/03/17 10:41:36.72              -32.7350  -71.7780                  23.2        204            0.70               IWREF      6353788
2001/02/19 08:24:20.16   0.27  1.10  23.5680   70.0720   8.0   4.8  16  10.0f       155       91  10.43 146.54     fe NEIC       4051004
 (Felt in Gujarat.)
2001/02/19 08:24:20.68   0.50  0.85  23.4849   70.0627  15.6  12.9  31   0.0f        30       96  10.51 142.83        IDC        4228590
2001/02/19 08:24:20.90   0.00  0.00  22.6400   69.3000   0.0   0.0   1         1.0   19   19      23.38  30.26     ke STR        4293675
2001/02/19 08:24:21.05               23.5960   70.1580                   8.6        264            0.60               IWREF      6353841
2001/02/19 08:24:21.40f        1.00  23.5790   70.1430f  5.9   5.4  -1   5.0f             42 148                  */
	istat = 0;
	istat += ReadFortranInt(line, 1, 4, &phypo->year);
	istat += ReadFortranInt(line, 6, 2, &phypo->month);
	istat += ReadFortranInt(line, 9, 2, &phypo->day);
	istat += ReadFortranInt(line, 12, 2, &phypo->hour);
	istat += ReadFortranInt(line, 15, 2, &phypo->min);
	istat += ReadFortranReal(line, 18, 5, &phypo->sec);

	istat += ReadFortranReal(line, 37, 8, &phypo->dlat);
	istat += ReadFortranReal(line, 46, 9, &phypo->dlong);

	istat += ReadFortranReal(line, 72, 5, &phypo->depth);

	istat += ReadFortranReal(line, 27, 3, &mag);
	istat += ReadFortranInt(line, 89, 4, &phypo->nreadings);
	istat += ReadFortranInt(line, 94, 3, &phypo->gap);
	istat += ReadFortranReal(line, 98, 6, &phypo->dist);
	istat += ReadFortranReal(line, 31, 5, &phypo->rms);

	istat += ReadFortranReal(line, 68, 3, &(phypo->ellipsoid.az1));
	phypo->ellipsoid.dip1 = 0.0;
	istat += ReadFortranReal(line, 62, 5, &(phypo->ellipsoid.len1));

	phypo->ellipsoid.az2 = phypo->ellipsoid.az1 + 90.0;
	phypo->ellipsoid.dip2 = 0.0;
	istat += ReadFortranReal(line, 56, 5, &(phypo->ellipsoid.len2));

	phypo->ellipsoid.len3 = 0.0;

	//istat += ReadFortranReal(line, 69, 2, &phypo->amp_mag);
	//phypo->amp_mag /= 10.0;
	//istat += ReadFortranReal(line, 71, 2, &phypo->dur_mag);
	//phypo->dur_mag /= 10.0;


	return(istat);

}



