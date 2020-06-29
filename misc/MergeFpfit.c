/*
 * Copyright (C) 1999-2003 Anthony Lomax <anthony@alomax.net>
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


/*   MergeFpfit.c

	Program to merge fpfit results into  NonLinLoc hyp files

*/



/*
	history:

	ver 01    01Jul2004  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

*/



#include "../src/GridLib.h"


/* defines */

#define MAX_NUM_INPUT_FILES 4096


/* globals */



/* functions */

int MergeFpfit(int argc, char** argv);
int compareTimes(int iy1, int iy2, int im1, int im2, int id1, int id2, int ih1, int ih2,
		int imn1, int imn2, double s1, double s2, double tolerance);


/*** program to merge fpfit results into NLL hyp file */

#define PNAME  "MergeFpfit"


int main(int argc, char** argv)
{

	int istat, narg;


	/* set program name */

	strcpy(prog_name, PNAME);


	/* check command line for correct usage */

	fprintf(stdout, "%s Arguments: ", prog_name);
	for (narg = 0; narg < argc; narg++)
		fprintf(stdout, "<%s> ", argv[narg]);
	fprintf(stdout, "\n");

	if (argc < 5) {
		nll_puterr("ERROR wrong number of command line arguments.");
		disp_usage(PNAME,
"<NLLhyp_file> <fpfit_file> <out_hyp_directory> time_tolerance [MechMisfitMax [RMSMax [NRdgsMin [GapMax]]]]");
		exit(-1);
	}

	SetConstants();
	message_flag = 1;
	DispProgInfo();
	message_flag = 0;

	if ((istat = MergeFpfit(argc, argv)) < 0) {
		nll_puterr("ERROR doing fpfit merge process.");
		exit(-1);
	}



	exit(0);

}



int MergeFpfit(int argc, char** argv)
{

	int istat1, istat2;
	int ok;
	int numNLLFiles;
	int nHypMerged, nFpfitRead, nNLLRead;

	char fn_hyp_in_list[MAX_NUM_INPUT_FILES][FILENAME_MAX];

	char fn_nllhyp_in[FILENAME_MAX];
	char fn_fpfit_in[FILENAME_MAX];
	char fdir_hyp_out[FILENAME_MAX], fn_out[FILENAME_MAX], filename[FILENAME_MAX];
	FILE *fp_fpfit_in;


	double time_tolerance;
	double idiff;

	double MechMisfitMax, RMSMax;
	int NRdgsMin, GapMax;

	GridDesc locgrid;
	HypoDesc HypoFpfit, HypoNLL;

	char test_str[10];




	strcpy(test_str, ".hyp");

	/* get command line parameters */
	strcpy(fn_nllhyp_in, argv[1]);
	strcpy(fn_fpfit_in, argv[2]);
	strcpy(fdir_hyp_out, argv[3]);
	sscanf(argv[4], "%lf", &time_tolerance);

	MechMisfitMax = 1.0e6;
	if (argc > 5) {
		sscanf(argv[5], "%lf", &MechMisfitMax);
	}
	fprintf(stdout, "  Mechanism Misfit Maximum: %lf\n", MechMisfitMax);
	RMSMax = 1.0e6;
	if (argc > 6) {
		sscanf(argv[6], "%lf", &RMSMax);
	}
	fprintf(stdout, "  RMS Maximum: %lf\n", RMSMax);
	NRdgsMin = 0;
	if (argc > 7) {
		sscanf(argv[7], "%d", &NRdgsMin);
	}
	fprintf(stdout, "  Num Readings Minimum: %d\n", NRdgsMin);
	GapMax = 360;
	if (argc > 8) {
		sscanf(argv[8], "%d", &GapMax);
	}
	fprintf(stdout, "  Gap Maximum: %d\n", GapMax);



	/* open fpfit summary input files */

	if ((fp_fpfit_in = fopen(fn_fpfit_in, "r")) == NULL) {
		nll_puterr2("ERROR: opening fpfit summary file: ", fn_fpfit_in);
		return(-1);
	}


	/* check for wildcards in input file name */
	strcat(fn_nllhyp_in, test_str);
	if ((numNLLFiles = ExpandWildCards(fn_nllhyp_in, fn_hyp_in_list, MAX_NUM_INPUT_FILES)) < 1) {
		nll_puterr("ERROR: no matching .hyp files found.");
		return(-1);
	}
	if (numNLLFiles >= MAX_NUM_INPUT_FILES) {
		sprintf(MsgStr, "WARNING: maximum number of event files exceeded, only first %d will be processed.", MAX_NUM_INPUT_FILES);
		nll_puterr(MsgStr);
	}


	/* merge cooresponding hypocenters */

	nHypMerged = 0;
	nFpfitRead = 0;
	nNLLRead = 0;

	while (1) {

		// read next fpfit hypocenter
	    	if ((istat2 = ReadFpfitSum(fp_fpfit_in, &HypoFpfit)) == EOF)
			break;
		else if (istat2 < 0) {
			nll_puterr2("ERROR: reading fpfit summary file", fn_fpfit_in);
			continue;
		}
		nFpfitRead++;

		while (istat2 >= 0 && istat2 != EOF && nNLLRead < numNLLFiles) {

			// read next NLL hypo
			sprintf(strstr(fn_hyp_in_list[nNLLRead], test_str), "\0");
			istat1 = GetHypLoc(NULL, fn_hyp_in_list[nNLLRead], &HypoNLL, Arrival,
					&NumArrivals, 1, &locgrid, 0);

			nNLLRead++;

			if (istat1 >= 0 && istat1 != EOF) {

				ok = 0;
				while (1) {

					idiff = compareTimes(HypoFpfit.year, HypoNLL.year, HypoFpfit.month, HypoNLL.month,
							HypoFpfit.day, HypoNLL.day, HypoFpfit.hour, HypoNLL.hour,
							HypoFpfit.min, HypoNLL.min, HypoFpfit.sec, HypoNLL.sec,
							time_tolerance);

/*
		printf("\nFPFIT:\n");
		WriteLocation(stdout, &HypoFpfit, Arrival, 0, fn_nllhyp_in, 0, 0, 1, &locgrid, 0);
		printf("NLL (%s):\n", fn_hyp_in_list[nNLLRead-1]);
		WriteLocation(stdout, &HypoNLL, Arrival, 0, fn_fpfit_in, 0, 0, 1, &locgrid, 0);
*/

					if (idiff == 0) {   	// match
						ok = 1;
						break;
					} else if (idiff < 0) {
						// fpfit hyp earlier, read next fpfit hyp
						if ((istat2 = ReadFpfitSum(fp_fpfit_in, &HypoFpfit)) == EOF)
							break;
						else if (istat2 < 0) {
							nll_puterr2("ERROR: reading fpfit summary file", fn_fpfit_in);
							continue;
						}
						nFpfitRead++;
					} else {
						// fpfit hyp later, read next NLL hyp
						break;
					}
				}
				if (!ok)
					continue;

				// check optional limits

				if (HypoFpfit.focMech.misfit > MechMisfitMax) {
		/*			nll_puterr(
		"WARNING: solution misfit is greater than MFMax, ignoring event");*/
					continue;
				} else if (HypoFpfit.rms > RMSMax) {
		/*			nll_puterr(
		"WARNING: location RMS is Greater than RMSMax, ignoring event");*/
					continue;
				} else if (HypoFpfit.nreadings < NRdgsMin) {
		/*			nll_puterr(
		"WARNING: location num readings is less than NRdgsMin, ignoring event");*/
					continue;
				} else if (HypoFpfit.gap > GapMax) {
		/*			nll_puterr(
		"WARNING: location gap is greater than GapMax, ignoring event");*/
					continue;
				}

				// merge fpfit info into NLL hyp

				HypoNLL.focMech = HypoFpfit.focMech;

				// output merged hyp

				strcpy(fn_out, fdir_hyp_out);
				strcat(fn_out, "/");
				if (strrchr(fn_hyp_in_list[nNLLRead - 1], '/') != NULL)
					strcpy(filename, strrchr(fn_hyp_in_list[nNLLRead - 1], '/'));
				else
					strcpy(filename, fn_hyp_in_list[nNLLRead - 1]);
				strcat(fn_out, filename);
				WriteGrid3dHdr(&locgrid, NULL, fn_out, NULL);
				strcat(fn_out, ".hyp");
				WriteLocation(NULL, &HypoNLL, Arrival, NumArrivals, fn_out, 1, 1, 0, &locgrid, 0);

				nHypMerged++;
//printf("MERGED!!!!!!!!!!! %d\n", nHypMerged);

			}

			break;  // next fpfit hypo

		}

	}

	printf("\n%d fpfit and %d NLL hypocenters read, %d merged.\n:", nFpfitRead, nNLLRead, nHypMerged);

	return(0);

}




int compareTimes(int iy1, int iy2, int im1, int im2, int id1, int id2, int ih1, int ih2,
		int imn1, int imn2, double s1, double s2, double tolerance)
{
	double time1, time2;

	if (iy1 != iy2)
		return(iy1 < iy2 ? -1 : 1);
	if (im1 != im2)
		return(im1 < im2 ? -1 : 1);
	if (id1 != id2)
		return(id1 < id2 ? -1 : 1);

	time1 = 60.0 * (60.0 * (double) ih1 + (double) imn1) + s1;
	time2 = 60.0 * (60.0 * (double) ih2 + (double) imn2) + s2;
	if (fabs(time1 - time2) > tolerance)
		return(time1 < time2 ? -1 : 1);

	return(0);

}


