/*
 * Copyright (C) 1999-2018 Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
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


/*   hyp2pretty_phases.c

        Program to convert NLLoc .hyp to a pretty phase format (originally for Alaska_historical_2017 project)

 */

/*-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
Mouans-Sartoux, France
e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/


/*
        history:

        ver 01    29Jul2008  AJL  Original version


 */



#include "../src/GridLib.h"


/* defines */

#define MAX_NUM_INPUT_FILES 5000


/* globals */



/* functions */

int hyp2pretty_phases_AlaskaHistorical(int, char **);



/*** program to sum event scatter files */

#define PNAME  "hyp2pretty_phases_AlaskaHistorical"

int main(int argc, char *argv[]) {

    int istat, narg;


    /* set program name */

    strcpy(prog_name, PNAME);


    /* check command line for correct usage */

    fprintf(stdout, "\n%s Arguments: ", prog_name);
    for (narg = 0; narg < argc; narg++)
        fprintf(stdout, "<%s> ", argv[narg]);
    fprintf(stdout, "\n");

    if (argc < 2) {
        nll_puterr("ERROR wrong number of command line arguments.");
        disp_usage(PNAME,
                "<hyp_files> ");
        exit(-1);
    }

    if ((istat = hyp2pretty_phases_AlaskaHistorical(argc, argv)) < 0) {
        nll_puterr("ERROR converting hyp file.");
        exit(-1);
    }



    exit(0);

}

int hyp2pretty_phases_AlaskaHistorical(int argc, char *argv[]) {

    int istat;
    int nLocWritten, nLocRead;
    int nFile;
    char fn_out[FILENAME_MAX];
    char fn_hyp_in[FILENAME_MAX];
    //char label_last[PHASE_LABEL_LEN];
    FILE *fp_out, *fp_hypo;

    GridDesc locgrid;
    HypoDesc Hypo;
    ArrivalDesc *parr;
    int narr;


    SetConstants();

    /* get command line parameters */
    strcpy(fn_hyp_in, argv[1]);


    for (nFile = 1; nFile < argc; nFile++) {

        nLocWritten = 0;
        nLocRead = 0;

        fprintf(OUT_LEVEL_1, "Adding location <%s>.\n", argv[nFile]);

        /* open hypocenter file */
        if ((fp_hypo = fopen(argv[nFile], "r")) == NULL) {
            nll_puterr2("ERROR: opening hypocenter file, ignoring event, file",
                    argv[nFile]);
            continue;
        }
        nLocRead++;

        /* open ascii output file */
        sprintf(fn_out, "%s.phs.txt", argv[nFile]);
        if ((fp_out = fopen(fn_out, "w")) == NULL) {
            nll_puterr("ERROR: opening output file.");
            return (-1);
        }
        fprintf(fp_out,
                "Station Phase OrigPhase Date HourMinSec TimeShift(sec)"
                " PriorWeight Residual(sec) PosteriorWeight"
                " Distance(deg) Azimuth(deg) Lat(deg) Lon(deg) Elev(m)\n");

        // loop over events in file

        while (1) {

            istat = GetHypLoc(fp_hypo, NULL, &Hypo, Arrival, &NumArrivals, 1, &locgrid, 0);
            if (istat == EOF) {
                break;
            }
            if (istat < 0) {
                nll_puterr2("ERROR: reading hypocenter file, ignoring event, file",
                        argv[nFile]);
                break;
            }

            if (strcmp(Hypo.locStat, "ABORTED") == 0) {
                //nll_puterr("WARNING: location ABORTED, ignoring event");
                continue;
            } else if (strcmp(Hypo.locStat, "REJECTED") == 0) {
                //nll_puterr("WARNING: location REJECTED, ignoring event");
                continue;
            }

            nLocWritten++;

            // write hypocenter
            /*
            fprintf(fp_out,
                    "evenement      date        heure           lat.         long.       prof.     rms\n");
            //    3       11/02/2006    20:26:30.16     43.62 N       5.60 E      7 km     0.10 sec

            fprintf(fp_out, "  %3d", nLocWritten);
            fprintf(fp_out, "       %2.2d/%2.2d/%4.4d", Hypo.day, Hypo.month, Hypo.year);
            fprintf(fp_out, "    %2.2d:%2.2d:%5.2f", Hypo.hour, Hypo.min, Hypo.sec);
            fprintf(fp_out, "    %6.2f N", Hypo.dlat);
            fprintf(fp_out, "    %7.2f E", Hypo.dlong);
            fprintf(fp_out, "    %3.0f km", Hypo.depth);
            fprintf(fp_out, "    %6.3f sec", Hypo.rms);

            fprintf(fp_out, "\n");

            for (narr = 0; narr < NumArrivals; narr++) {
                Arrival[narr].obs_time = (long double) Arrival[narr].sec
                        + 60.0L * ((long double) Arrival[narr].min
                        + 60.0L * (long double) Arrival[narr].hour);
            }
            if (SortArrivalsTime(Arrival, NumArrivals) < 0)
                    nll_puterr("ERROR: sorting arrivals by time.");
            if (SortArrivalsDist(Arrival, NumArrivals) < 0)
                    nll_puterr("ERROR: sorting arrivals by distance.");

            strcpy(label_last, "XXX");
            for (narr = 0; narr < NumArrivals; narr++) {
                    parr = Arrival + narr;
                    if (strcmp(parr->phase, "S") == 0) {	// S
                            if (strcmp(parr->label, label_last) == 0) {	// same station, S after P
                                    fprintf(fp_out, "      S:");
                                    strcpy(label_last, "$$$");
                            } else {	// lone S
                                    fprintf(fp_out, "\n");
                                    fprintf(fp_out, "    %s", parr->label);
                                    fprintf(fp_out, "                                            S:");
                                    strcpy(label_last, "$$$");
                            }
                            fprintf(fp_out, "  %2.2d:%2.2d:%5.2f", parr->hour, parr->min, parr->sec);
                    } else {	// P
                            fprintf(fp_out, "\n");
                            fprintf(fp_out, "    %s", parr->label);
                            fprintf(fp_out, "        %s%s:", parr->phase, strcmp(parr->first_mot, "?") == 0 ? " " : parr->first_mot);
                            strcpy(label_last, parr->label);
                            fprintf(fp_out, "  %2.2d:%2.2d:%6.3f", parr->hour, parr->min, parr->sec);
                    }
                    fprintf(fp_out, " +/-%5.2f sec", parr->error);
            }

            if (strcmp(label_last, "$$$") == 0)	// S after P
                    fprintf(fp_out, "\n");
            fprintf(fp_out, "\n");
             */

            // write observation part of FORMAT_PHASE_2 phase line

            // Station Phase OrigPhase Date HourMinSec TimeShift(sec) PriorWeight Residual(sec) PosteriorWeight Distance(deg) Azimuth(deg) Lat(deg) Lon(deg) Elev(m)
            //PWA	P	P	19750101	03:55:22.4	-	1.0	-1.4	1.02	0.1	211.5	61.6508	-149.8790	137

            for (narr = 0; narr < NumArrivals; narr++) {

                parr = Arrival + narr;

                istat = fprintf(fp_out,
                        "%-5s %-6s %-6s %8.8d %2.2d:%2.2d:%05.2lf %5s %7.2lf %7.2lf %7.2lf %6.1lf %6.1lf %9.4lf %9.4lf %5.0lf",
                        parr->label,
                        parr->phase,
                        parr->inst,
                        parr->year * 10000 + parr->month * 100 + parr->day,
                        parr->hour, parr->min, parr->sec,
                        strncmp(parr->comp, "TS", 2) == 0 ? parr->comp : "-",
                        parr->apriori_weight,
                        parr->residual,
                        parr->weight,
                        parr->dist,
                        parr->azim,
                        parr->station.y > -9999.9 ? parr->station.y : -999.0,
                        parr->station.x > -9999.9 ? parr->station.x : -999.0,
                        parr->station.z > -9999.9 ? -parr->station.z * 1000.0 : -999.0
                        );

                fprintf(fp_out, "\n");

            }

            fclose(fp_out);
            fclose(fp_hypo);

            /* write message */
            fprintf(stdout,
                    "%d locations read, %d written to ascii bulletin file <%s>\n",
                    nLocRead, nLocWritten, fn_out);

        }

    }

    return (0);

}



