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


/*   fpfit2hyp.c

	Program to convert fpfit summary files to NLLoc .hyp format

*/

/*------------------------------------------------------------/ */
/* Anthony Lomax           | email: lomax@faille.unice.fr     / */
/* UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax / */
/* 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        / */
/* 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        / */
/*------------------------------------------------------------/ */


/*
	history:

	ver 01    24Nov1998  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

*/



#include "GridLib.h"


/* defines */



/* globals */



/* functions */

int Fpfit2Hyp(int , char **);



/*** program to sum event scatter files */

#define PNAME  "fpfit2hyp"


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
"<fpfit_file> <out_hyp_file> [MechMisfitMax [RMSMax [NRdgsMin [GapMax]]]]");
		exit(-1);
	}

	if ((istat = Fpfit2Hyp(argc, argv)) < 0) {
		nll_puterr("ERROR converting fpfitise summary file.");
		exit(-1);
	}



	exit(0);

}



int Fpfit2Hyp(int argc, char *argv[])
{
	int istat;
	int nLocWritten, nLocRead;
	char fn_hyp_out[FILENAME_MAX];
	char fn_fpfit_in[FILENAME_MAX];
	FILE *fp_hyp_out, *fp_fpfit_in;


	double MechMisfitMax, RMSMax;
	int NRdgsMin, GapMax;

	GridDesc locgrid;
	HypoDesc Hypo, *phypo;



	/* get command line parameters */
	strcpy(fn_fpfit_in, argv[1]);
	strcpy(fn_hyp_out, argv[2]);

	MechMisfitMax = 1.0e6;
	if (argc > 3) {
		sscanf(argv[3], "%lf", &MechMisfitMax);
	}
	fprintf(stdout, "  Mechanism Misfit Maximum: %lf\n", MechMisfitMax);
	RMSMax = 1.0e6;
	if (argc > 4) {
		sscanf(argv[4], "%lf", &RMSMax);
	}
	fprintf(stdout, "  RMS Maximum: %lf\n", RMSMax);
	NRdgsMin = 0;
	if (argc > 5) {
		sscanf(argv[5], "%d", &NRdgsMin);
	}
	fprintf(stdout, "  Num Readings Minimum: %d\n", NRdgsMin);
	GapMax = 360;
	if (argc > 6) {
		sscanf(argv[6], "%d", &GapMax);
	}
	fprintf(stdout, "  Gap Maximum: %d\n", GapMax);



	/* open fpfit summary file */

	if ((fp_fpfit_in = fopen(fn_fpfit_in, "r")) == NULL) {
		nll_puterr("ERROR: opening fpfit summary file.");
		return(-1);
	}

	/* open ascii hypocenter output file */
	if ((fp_hyp_out = fopen(fn_hyp_out, "w")) == NULL) {
		nll_puterr("ERROR: opening NLL hyp output file.");
		return(-1);
	}



	nLocWritten = 0;
	nLocRead = 0;
	while (1) {


		strcpy(Hypo.fileroot, fn_fpfit_in);
		strcpy(Hypo.locStat, "Converted from FPFIT");
		Hypo.grid_misfit_max = -1.0;

		/* read next hypocenter */
	    	if ((istat = ReadFpfitSum(fp_fpfit_in, &Hypo)) == EOF)
			break;
		else if (istat < 0) {
			nll_puterr2(
"ERROR: reading fpfit summary file", fn_fpfit_in);
			break;
		}
		nLocRead++;

		if (Hypo.focMech.misfit > MechMisfitMax) {
/*			nll_puterr(
"WARNING: solution misfit is greater than MFMax, ignoring event");*/
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
		fprintf(fp_hyp_out, " Ell3  %.2le", phypo->ellipsoid.len3);
		fprintf(fp_hyp_out, "\n");


		/* write mechanism */
		phypo = &Hypo;
		fprintf(fp_hyp_out,
			"FOCALMECH  Hyp  %lf %lf %lf",
			phypo->dlat, phypo->dlong, phypo->depth);
		fprintf(fp_hyp_out,
			" Mech  %.1lf %.1lf %.1lf",
			phypo->focMech.dipDir, phypo->focMech.dipAng,
			phypo->focMech.rake);
		fprintf(fp_hyp_out, " mf  %.2lf nObs %d",
			phypo->focMech.misfit, phypo->focMech.nObs);
		fprintf(fp_hyp_out, "\n");


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




/*------------------------------------------------------------/ */
/* Anthony Lomax           | email: lomax@faille.unice.fr     / */
/* UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax / */
/* 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        / */
/* 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        / */
/*------------------------------------------------------------/ */

