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


/*   Time2EQ.c

        Program to calculate travel times for 3-D grids

 */

/*------------------------------------------------------------/ */
/* Anthony Lomax           | email: lomax@faille.unice.fr     / */
/* UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax / */
/* 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        / */
/* 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        / */
/*------------------------------------------------------------/ */


/*
        history:

        ver 01    25SEP1997  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

 */



#define PNAME  "Time2EQ"

#include "GridLib.h"
#include "ran1/ran1.h"


/* defines */

#define MODE_SRCE_TO_STA 0
#define MODE_STA_TO_SRCE 1
#define MODE_UNDEF  -1

/* globals  */
char EventStr[MAXLINE], fn_hypo_output[MAXLINE], fn_simulps_output[MAXLINE], fn_invc_output[MAXLINE];
char fn_gmt_output[MAXLINE];
char fn_eq_input[MAXLINE], fn_eq_output[MAXLINE];
int time2eq_mode;
double VpVsRatio;
int NumStationPhases;
//
// auto set origin time so multiple sources do not overlap in time and have the same hyp filename
// 20201013 AJL - added
double otime_EQSRCE = 0.0;
double otime_EQSRCE_step = 60.0;

/* mechanism */
#define  MECH_NONE   0
#define  MECH_DOUBLE   1
#define  MECH_ISOTROPIC  2
int imech;
char mech_type[11];
double mech_phi, mech_del, mech_lam;
char mech_str[MAXLINE];


/* event */
SourceDesc *Event;




/* function declarations */

int GetTime2EQ_Mode(char* line1);
int GetTime2EQ_Event(char* in_line);
int GetTime2EQ_Stations(char* in_line);
int GetTime2EQ_Source(char*);
int ReadTime2EQ_Input(FILE*);
int GetTime2EQ_Files(char*);
int get_vp_vs(char* line1);
double CalcArrivalTime(FILE*, GridDesc*, SourceDesc*, StationDesc*);
double AddNoise(double arrival_time, StationDesc* psta);
int CalcFirstMotion(char *, GridDesc*, SourceDesc*, StationDesc*, int);
int WritePhaseArrival(double, int, FILE*, FILE*, FILE*, FILE*, FILE*, FILE*, SourceDesc*,
        StationDesc*, int);
int get_mech(char*);
double calc_rad(double, double, char);




/*** program to generate synthetic travel times to stations */

#define NARGS 2

int main(int argc, char *argv[]) {

    int istat;
    int nsta, nsource;
    double arrival_time;
    int ipolarity;
    int nsta_written;

    char filename[MAXLINE], filename_angle[MAXLINE];
    FILE *fp_eq_output, *fp_hypo_output, *fp_simulps_output, *fp_invc_output;
    FILE *fp_gmt_output, *fp_gmt_output_az;
    FILE *fp_time_grid, *fp_time_hdr;

    GridDesc grid0;
    SourceDesc srce0;

    char MsgStr_sta[2 * MAXLINE], MsgStr_srce[2 * MAXLINE];

    char last_label[MAXLINE], last_phs_label[ARRIVAL_LABEL_LEN];
    double last_arrival_time = -1.0;


    /* set program name */
    strcpy(prog_name, PNAME);

    /* check command line for correct usage */

    if (argc != NARGS) {
        disp_usage(prog_name, "<control file>");
        exit(EXIT_ERROR_USAGE);
    }



    /* set constants */

    SetConstants();
    NumStationPhases = 0;
    NumSources = 0;
    imech = MECH_NONE;
    time2eq_mode = MODE_SRCE_TO_STA;
    VpVsRatio = -1.0;

    grid0.iSwapBytes = 0;


    /* read control file */

    strcpy(fn_control, argv[1]);
    if ((fp_control = fopen(fn_control, "r")) == NULL) {
        nll_puterr("ERROR: opening control file.");
        exit(EXIT_ERROR_FILEIO);
    }

    if ((istat = ReadTime2EQ_Input(fp_control)) < 0) {
        exit(EXIT_ERROR_FILEIO);
    }

    /* initialize random number generator */

    SRAND_FUNC(RandomNumSeed);
    if (message_flag >= 3)
        test_rand_int();
    if (message_flag >= 3)
        test_normal_dist_deviate();



    /* convert source location coordinates  */

    if ((istat = ConvertSourceLoc(0, Source, NumSources, 1, 1)) < 0)
        nll_puterr("ERROR: converting source location to x/y and lat/lon.");


    /* open arrivals output file */

    if ((fp_eq_output = fopen(fn_eq_output, "w")) == NULL) {
        nll_puterr("ERROR: opening eq times output file.");
        exit(EXIT_ERROR_FILEIO);
    }

    sprintf(fn_hypo_output, "%s.hypo", fn_eq_output);
    if ((fp_hypo_output = fopen(fn_hypo_output, "w")) == NULL) {
        nll_puterr("ERROR: opening eq times HYPO output file.");
        exit(EXIT_ERROR_FILEIO);
    }

    sprintf(fn_simulps_output, "%s.sim", fn_eq_output);
    if ((fp_simulps_output = fopen(fn_simulps_output, "w")) == NULL) {
        nll_puterr("ERROR: opening eq times SIMULPS (Alberto) output file.");
        exit(EXIT_ERROR_FILEIO);
    }


    sprintf(fn_invc_output, "%s.invc", fn_eq_output);
    if ((fp_invc_output = fopen(fn_invc_output, "w")) == NULL) {
        nll_puterr("ERROR: opening eq times INVCOUCHE output file.");
        exit(EXIT_ERROR_FILEIO);
    }


    sprintf(fn_gmt_output, "%s.gmtxy", fn_eq_output);
    if ((fp_gmt_output = fopen(fn_gmt_output, "w")) == NULL) {
        nll_puterr("ERROR: opening eq times GMT output file.");
        exit(EXIT_ERROR_FILEIO);
    }

    sprintf(fn_gmt_output, "%s_az.gmtxy", fn_eq_output);
    if ((fp_gmt_output_az = fopen(fn_gmt_output, "w")) == NULL) {
        nll_puterr("ERROR: opening eq times GMT azimuth output file.");
        exit(EXIT_ERROR_FILEIO);
    }



    /* generate arrival times for each source */

    for (nsource = 0; nsource < NumSources; nsource++) {

        Event = Source + nsource;

        sprintf(MsgStr_srce,
                "Calculating travel time for Source: %s  X %.2lf  Y %.2lf  Z %.2lf  OT %.2lf",
                Event->label, Event->x, Event->y, Event->z, Event->otime);
        if (time2eq_mode == MODE_STA_TO_SRCE)
            nll_putmsg(3, MsgStr_srce);
        else
            nll_putmsg(2, MsgStr_srce);

        sprintf(EventStr,
                "EQEVENT:  Label: %s  Loc:  X %.2lf  Y %.2lf  Z %.2lf  OT %.2lf",
                Event->label, Event->x, Event->y, Event->z, Event->otime);
        fprintf(fp_eq_output, "# %s\n", EventStr);
        fprintf(fp_eq_output, "PUBLIC_ID  %s\n", Event->label);

        // INVCOUCHE
        // T5  864979  177041  230   0.1490
        fprintf(fp_invc_output, "%s  %f %f elevation_toit-%f  %f", Event->label, 1000.0 * Event->x, 1000.0 * Event->y, 1000.0 * Event->z, Event->otime);

        /* generate arrival times for each station */

        nsta_written = 0;
        for (nsta = 0; nsta < NumStationPhases; nsta++) {

            // check for S and valid VpVs
            if (strcmp((Station + nsta)->label, last_label) == 0
                    && strcmp((Station + nsta)->phs[0].label, "S") == 0
                    && strcmp(last_phs_label, "P") == 0
                    && VpVsRatio > 0.0 && last_arrival_time > 0.0) {
                arrival_time = Event->otime + VpVsRatio * (last_arrival_time - Event->otime);

                sprintf(MsgStr_sta,
                        "Calculating travel time for station / S phase: %s  %s",
                        (Station + nsta)->label,
                        (Station + nsta)->phs[0].label);
            } else {
                // not using VpVs

                sprintf(filename, "%s.%s.%s.time", fn_eq_input,
                        (Station + nsta)->phs[0].label,
                        (Station + nsta)->label);
                if ((istat = OpenGrid3dFile(filename, &fp_time_grid,
                        &fp_time_hdr, &grid0, "time", &srce0, grid0.iSwapBytes))
                        < 0) {
                    CloseGrid3dFile(&grid0, &fp_time_grid, &fp_time_hdr);
                    nll_puterr2("ERROR: opening time grid files", filename);
                    continue;
                }
                sprintf(filename_angle, "%s.%s.%s.angle", fn_eq_input,
                        (Station + nsta)->phs[0].label,
                        (Station + nsta)->label);
                (Station + nsta)->x = srce0.x;
                (Station + nsta)->y = srce0.y;
                (Station + nsta)->z = srce0.z;

                arrival_time = CalcArrivalTime(fp_time_grid, &grid0, Event, Station + nsta);
                CloseGrid3dFile(&grid0, &fp_time_grid, &fp_time_hdr);
                sprintf(MsgStr_sta,
                        "Calculating travel time for station / phase: %s  %s  X %.2lf  Y %.2lf  Z %.2lf",
                        (Station + nsta)->label,
                        (Station + nsta)->phs[0].label,
                        (Station + nsta)->x,
                        (Station + nsta)->y,
                        (Station + nsta)->z);
            }

            strcpy(last_label, (Station + nsta)->label);
            strcpy(last_phs_label, (Station + nsta)->phs[0].label);
            last_arrival_time = arrival_time;

            if (arrival_time < 0.0) {
                nll_puterr(MsgStr_srce);
                nll_puterr(MsgStr_sta);
                sprintf(MsgStr,
                        "ERROR: calculating travel time: t = %lf", arrival_time);
                nll_puterr(MsgStr);
                continue;
            }

            if (time2eq_mode == MODE_SRCE_TO_STA)
                nll_putmsg(3, MsgStr_sta);
            else
                nll_putmsg(2, MsgStr_sta);
            ipolarity = 0;
            if (imech == MECH_DOUBLE) {
                int isP = strncmp((Station + nsta)->phs[0].label, "P", 1) == 0;
                ipolarity = CalcFirstMotion(filename_angle, &grid0, Event, Station + nsta, isP);
                if (!isP) {
                    // not P, no polarity possible
                    ipolarity = 0;
                }
            } else if (imech == MECH_ISOTROPIC) {
                ipolarity = 1;
            }
            // AJL 20070731 - bug fix: new function AddNoise treats S with VpVs correctly
            arrival_time = AddNoise(arrival_time, Station + nsta);

            // check if active
            if ((Station + nsta)->prob_active < 1.0) {
                if (get_rand_double(0.0, 1.0) > (Station + nsta)->prob_active)
                    continue;
            }


            if ((istat = WritePhaseArrival(arrival_time, ipolarity,
                    fp_eq_output, fp_hypo_output, fp_simulps_output, fp_invc_output,
                    fp_gmt_output, fp_gmt_output_az,
                    Event, Station + nsta, nsta_written++)) < 0)
                nll_puterr("ERROR: writing phase arrival.");

        }
        fprintf(fp_eq_output, "#\n");
        fprintf(fp_hypo_output, "\n                 10\n");
        fprintf(fp_invc_output, "\n  END_STA\n\n");
        fprintf(fp_simulps_output, "\n0\n");
    }

    fprintf(fp_gmt_output, ">\n");
    fprintf(fp_gmt_output_az, ">\n");

    fclose(fp_eq_output);

    fclose(fp_hypo_output);
    fclose(fp_invc_output);
    fclose(fp_simulps_output);
    fclose(fp_gmt_output);
    fclose(fp_gmt_output_az);

    exit(EXIT_NORMAL);

}

/*** function to calc arrival time from time grid file and add noise */

double CalcArrivalTime(FILE* fpgrid, GridDesc* ptgrid, SourceDesc* pevent, StationDesc* psta) {

    double arrival_time, yval_grid;


    /* get travel time */


    if (ptgrid->type == GRID_TIME) {
        /* 3D grid */
        arrival_time = pevent->otime +
                ReadAbsInterpGrid3d(fpgrid, ptgrid, pevent->x, pevent->y, pevent->z, 0);
    } else {
        /* 2D grid (1D model) */
        yval_grid = GetEpiDistSta(psta, pevent->x, pevent->y);
        //printf("CalcArrivalTime: xyz= %f %f %f  yval_grid=%f\n", pevent->x, pevent->y, pevent->z, yval_grid);
        arrival_time = pevent->otime + ReadAbsInterpGrid2d(fpgrid, ptgrid, yval_grid, pevent->z);
        //printf("CalcArrivalTime: arrival_time=%f  dist=%f %f\n", arrival_time, sqrt(yval_grid * yval_grid + (pevent->z - psta->z) * (pevent->z - psta->z)), sqrt((pevent->x - psta->x) * (pevent->x- psta->x)+ (pevent->y - psta->y) * (pevent->y - psta->y) + (pevent->z - psta->z) * (pevent->z - psta->z)));
    }

    return (arrival_time);

}

/*** function to calc arrival time from time grid file and add noise */

double AddNoise(double arrival_time, StationDesc* psta) {

    double noise = 0.0;


    /* add noise */

    double error = psta->phs[0].error;
    if (psta->phs[0].prob_outlier > 0.0) {
        if (get_rand_double(0.0, 1.0) < psta->phs[0].prob_outlier) {
            error *= psta->phs[0].outlier_err_factor;
        }
    }

    if (strcmp(psta->phs[0].error_type, "GAU") == 0)
        noise = error * normal_dist_deviate();
    else if (strcmp(psta->phs[0].error_type, "BOX") == 0)
        noise = get_rand_double(-error, error);
    else if (strcmp(psta->phs[0].error_type, "FIX") == 0)
        noise = error;
    else if (strcmp(psta->phs[0].error_type, "NONE") == 0)
        noise = 0.0;
    else
        nll_puterr2("ERROR: unrecognized error type:", psta->phs[0].error_type);

    sprintf(MsgStr, "Station %s  Phase %s  Error Type %s  Error %f  ->  Noise %f + ArrivalTime %f = Sum %f", psta->label, psta->phs[0].label, psta->phs[0].error_type, error, noise, arrival_time, noise + arrival_time
            );
    nll_putmsg(2, MsgStr);

    return (arrival_time + noise);

}

/*** function to calc first motion from angle grid file */

int CalcFirstMotion(char *filename, GridDesc* ptgrid, SourceDesc* pevent, StationDesc* psta, int isP) {

    double yval_grid, azim;
    double radamp = 0.0;
    int ipolarity = 0;
    double ray_azim, ray_dip;
    int ray_qual;

    /* get take-off angles */

    if (ptgrid->type == GRID_TIME) {
        /* 3D grid */
        if (ReadTakeOffAnglesFile(filename,
                pevent->x, pevent->y, pevent->z,
                &ray_azim, &ray_dip, &ray_qual, -1.0, ptgrid->iSwapBytes) < 0)
            return (0);
    } else {
        /* 2D grid (1D model) */
        yval_grid = GetEpiDistSta(psta, pevent->x, pevent->y);
        azim = GetEpiAzimSta(psta, pevent->x, pevent->y);
        if (ReadTakeOffAnglesFile(filename,
                0.0, yval_grid, pevent->z,
                &ray_azim, &ray_dip, &ray_qual, azim, ptgrid->iSwapBytes) < 0)
            return (0);
    }

    /* calc radiation amplitude and polarity */

    if (isP) {
        radamp = calc_rad(ray_dip * cRPD, ray_azim * cRPD, 'P');
        ipolarity = radamp < 0.0 ? -1 : 1;
    }

    sprintf(MsgStr, "FM: Sta: %s  Pha: %s  ray_dip %.1lf  ray_azim %.1lf  radamp %.2le  ipol %d\n",
            psta->label, psta->phs[0].label, ray_dip, ray_azim, radamp, ipolarity);
    nll_putmsg(2, MsgStr);

    return (ipolarity);

}

/*** function to calc arrival time from grid file and save to event file */

int WritePhaseArrival(double arrival_time, int ipolarity,
        FILE* fpout, FILE* fpout_hypo, FILE* fpout_simulps, FILE* fpout_invc,
        FILE* fpout_gmt, FILE* fpout_gmt_az,
        SourceDesc* pevent, StationDesc* psta, int num_sta) {

    int iEOL;
    int nhr = 0, nmin = 0;
    int nlat = 0, nlon = 0;
    double osec = 0.0;
    ArrivalDesc arr, *parr;
    double az;

    static ArrivalDesc lastarr;

    parr = &arr;

    /* determine label */

    if (time2eq_mode == MODE_SRCE_TO_STA)
        strcpy(parr->label, psta->label);
    else if (time2eq_mode == MODE_STA_TO_SRCE)
        strcpy(parr->label, pevent->label);

    strcpy(parr->inst, "?");
    strcpy(parr->comp, "?"),
            strcpy(parr->onset, "?"),
            strcpy(parr->phase, psta->phs[0].label);
    if (ipolarity == 1)
        strcpy(parr->first_mot, "U");
    else if (ipolarity == -1)
        strcpy(parr->first_mot, "D");
    else
        strcpy(parr->first_mot, "?");
    parr->quality = -1;
    parr->year = 1900;
    parr->month = 01;
    parr->day = 01;

    parr->hour = (int) (arrival_time / 3600.0);
    osec = arrival_time - 3600.0 * (double) parr->hour;
    parr->min = (int) (osec / 60.0);
    parr->sec = osec - 60.0 * (double) parr->min;

    strcpy(parr->error_type, psta->phs[0].error_report_type);
    parr->error = psta->phs[0].error_report;
    parr->coda_dur = 0.0;
    parr->amplitude = 0.0;
    parr->period = 0.0;

    parr->apriori_weight = 1.0; // 20110105 AJL - to support NLL_FORMAT_VER_2

    /* write arrival time to NLLoc output file */

    PhaseFormat = FORMAT_PHASE_2; // 20110105 AJL - to allow long station names
    if (WriteArrival(fpout, parr, IO_ARRIVAL_OBS) < 0) {
        nll_puterr("ERROR: writing arrival time to disk.\n");
        return (-1);
    }

    /* write arrival time to INVCOUCHE output file */

    if (fpout_invc != NULL) {
        fprintf(fpout_invc, "\n   %-20.20s %f", parr->label, arrival_time);
    }

    /* write arrival time to HYPO output file */

    if (fpout_hypo != NULL) {
        iEOL = 1;
        if (strcmp(lastarr.label, parr->label) == 0
                && strcmp(lastarr.phase, "P") == 0
                && strcmp(parr->phase, "S") == 0)
            iEOL = 0;
        if (WriteArrivalHypo(fpout_hypo, parr, iEOL) < 0) {
            nll_puterr("ERROR: writing hypo71 arrival time to disk.\n");
            return (-1);
        }
        lastarr = *parr;
    }

    /* write arrival time to sim (alberto) output file */

    if (fpout_simulps != NULL) {

        nhr = (int) (pevent->otime / 3600.0);
        osec = pevent->otime - 3600.0 * (double) nhr;
        nmin = (int) (osec / 60.0);
        osec = osec - 60.0 * (double) nmin;
        if (num_sta == 0) {
            //87 1 1  1 2  0.00 35N14.26 120W44.00   5.15   0.00
            nlat = (int) fabs(pevent->dlat);
            nlon = (int) fabs(pevent->dlong);
            fprintf(fpout_simulps,
                    "%2.2d%2.2d%2.2d %2.2d%2.2d%6.2f %2.2d%c%5.2f %3.3d%c%5.2f %6.2f %6.2f",
                    99, 01, 01, nhr, nmin, osec,
                    nlat, pevent->dlat > 0.0 ? 'N' : 'S',
                    60.0 * (fabs(pevent->dlat) - (double) nlat),
                    nlon, pevent->dlong > 0.0 ? 'E' : 'W',
                    60.0 * (fabs(pevent->dlong) - (double) nlon),
                    pevent->depth, 0.0
                    );
        }
        if (num_sta % 5 == 0)
            fprintf(fpout_simulps, "\n");
        //ST06iP 013.8274
        fprintf(fpout_simulps, "%4s%1s%1s%2.2d%7.4f",
                parr->label,
                strcmp(parr->onset, "?") == 0 ? "i" : parr->onset,
                parr->phase,
                parr->min, parr->sec
                );
    }

    /* write dist/az to GMT output file */

    if (fpout_gmt != NULL) {
        if ((az = GetEpiAzimSta(psta, pevent->x, pevent->y)) < 90.0
                || az > 270.0) {
            fprintf(fpout_gmt, "%lf %lf %s %s\n",
                    GetEpiDistSta(psta, pevent->x, pevent->y), arrival_time,
                    pevent->label, psta->label);
        } else {
            fprintf(fpout_gmt, "%lf %lf %s %s\n",
                    -GetEpiDistSta(psta, pevent->x, pevent->y), arrival_time,
                    pevent->label, psta->label);
        }
    }
    if (fpout_gmt_az != NULL) {
        fprintf(fpout_gmt_az, "%lf %lf %s %s\n",
                GetEpiAzim(pevent, psta->x, psta->y), arrival_time,
                pevent->label, psta->label);
    }


    return (0);

}

/*** function to read input file */

int ReadTime2EQ_Input(FILE* fp_input) {
    int istat, iscan;
    char param[MAXLINE];
    char line[4 * MAXLINE], *fgets_return;

    int flag_control = 0, flag_outfile = 0, flag_stations = 0,
            flag_event = 0, flag_source = 0, flag_mode = 0,
            flag_trans = 0, flag_qual2err = 0, flag_mech = 0, flag_vp_vs = 0;

    int flag_include = 1;


    /* read each input line */

    while ((fgets_return = fgets(line, 4 * MAXLINE, fp_input)) != NULL
            || fp_include != NULL) {


        /* check for end of include file */

        if (fgets_return == NULL && fp_include != NULL) {
            SwapBackIncludeFP(&fp_input);
            continue;
        }


        istat = -1;

        /*read parameter line */

        if ((iscan = sscanf(line, "%s", param)) < 0)
            continue;

        /* skip comment line or white space */

        if (strncmp(param, "#", 1) == 0 || isspace(param[0]))
            istat = 0;


        /* read include file params and set input to include file */

        if (strcmp(param, "INCLUDE") == 0)
            if ((istat = GetIncludeFile(strchr(line, ' '),
                    &fp_input)) < 0) {
                nll_puterr("ERROR: processing include file.");
                flag_include = 0;
            }


        /* read control params */

        if (strcmp(param, "CONTROL") == 0) {
            if ((istat = get_control(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading control params.");
            else
                flag_control = 1;
        }


        /*read transform params */

        if (strcmp(param, "TRANS") == 0) {
            if ((istat = get_transform(0, strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading transformation parameters.");
            else
                flag_trans = 1;
        }


        /* read file names */

        if (strncmp(param, "EQFILES", 7) == 0) {
            if ((istat = GetTime2EQ_Files(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading output file name.");
            else
                flag_outfile = 1;
        }


        /* read grid mode names */

        if (strcmp(param, "EQMODE") == 0) {
            if ((istat = GetTime2EQ_Mode(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading mode.");
            else
                flag_mode = 1;
        }


        /* read source params */

        if (strncmp(param, "EQEVENT", 7) == 0) {
            if ((istat = GetTime2EQ_Event(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading Time2EQ event params.");
            else
                flag_event = 1;
        }


        /* read mechanism params */

        if (strncmp(param, "EQMECH", 6) == 0) {
            if ((istat = get_mech(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading Time2EQ mechanism params.");
            else
                flag_mech = 1;
        }


        /* read VpVs params */

        if (strncmp(param, "EQVPVS", 6) == 0) {
            if ((istat = get_vp_vs(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading Time2EQ Vp/Vs params.");
            else
                flag_vp_vs = 1;
        }


        /* read source params */

        if (strcmp(param, "EQSRCE") == 0) {
            if ((istat = GetNextSource(strchr(line, ' '))) < 0) {
                nll_puterr("ERROR: reading source params.");
            } else {
                flag_source = 1;
                // auto set origin time so multiple sources do not overlap in time and have the same hyp filename
                // 20201013 AJL - added
                SourceDesc *srce = Source + NumSources - 1;
                srce->otime = otime_EQSRCE;
                otime_EQSRCE += otime_EQSRCE_step;
            }
        }


        /* read station params */

        if (strncmp(param, "EQSTA", 5) == 0) {
            if ((istat = GetTime2EQ_Stations(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading station params.");
            else
                flag_stations = 1;
        }


        if (strcmp(param, "EQQUAL2ERR") == 0) {
            if ((istat = GetQuality2Err(strchr(line, ' '))) < 0)
                nll_puterr("ERROR: reading quality2error values.");
            else
                flag_qual2err = 1;
        }


        /* unrecognized input */

        if (istat < 0) {
            sprintf(MsgStr, "Skipping input: %s", line);
            nll_putmsg(4, MsgStr);
        }

    }


    /* check for missing input */

    if (!flag_control)
        nll_puterr("ERROR: no control (CONTROL) params read.");
    if (!flag_trans)
        nll_puterr("INFO: no transformation (TRANS) params read.");
    if (!flag_mode) {
        sprintf(MsgStr, "INFO: no mode (EQMODE) params read.");
        nll_putmsg(1, MsgStr);
        sprintf(MsgStr, "INFO: DEFAULT: \"SRCE_TO_STA\"");
        nll_putmsg(1, MsgStr);
    }
    if (!flag_outfile)
        nll_puterr("ERROR: no i/o file (EQFILES) params read.");
    if (!flag_event && !flag_source)
        nll_puterr(
            "ERROR: no event (EQEVENT) or source (EQSRCE) params read.");
    if (!flag_stations)
        nll_puterr("ERROR: no station (EQSTA) params read.");
    if (!flag_qual2err)
        nll_puterr("ERROR: no quality to error (EQQUAL2ERR) params read.");
    if (!flag_mech)
        nll_puterr("INFO: no mechanism (EQMECH) params read.");
    if (!flag_vp_vs)
        nll_puterr("INFO: no Vp/Vs (EQVPVS) params read.");


    return (flag_include * flag_control * flag_outfile * flag_qual2err *
            (flag_event || flag_source) * flag_stations - 1);
}

/*** function to read output file name ***/

int GetTime2EQ_Files(char* line1) {

    sscanf(line1, "%s %s", fn_eq_input, fn_eq_output);

    sprintf(MsgStr, "Time2EQ FILES:  Input: %s.*  Output: %s",
            fn_eq_input, fn_eq_output);
    nll_putmsg(2, MsgStr);

    return (0);
}

/*** function to read mode params ***/

int GetTime2EQ_Mode(char* line1) {
    char str_mode[MAXLINE];


    sscanf(line1, "%s", str_mode);

    sprintf(MsgStr, "Time2EQ EQMODE:  %s", str_mode);
    nll_putmsg(2, MsgStr);

    if (strcmp(str_mode, "SRCE_TO_STA") == 0)
        time2eq_mode = MODE_SRCE_TO_STA;
    else if (strcmp(str_mode, "STA_TO_SRCE") == 0)
        time2eq_mode = MODE_STA_TO_SRCE;
    else {
        time2eq_mode = MODE_UNDEF;
        nll_puterr2("ERROR: unrecognized mode:", str_mode);
        return (-1);
    }

    return (0);

}

/** function to read event params fom input line */

int GetTime2EQ_Event(char* in_line) {
    int istat;


    Event = Source + NumSources;
    NumSources++;

    Event->is_coord_xyz = 1;
    Event->is_coord_latlon = 0;


    /* read event input line */

    istat = sscanf(in_line, "%s %lf %lf %lf %lf",
            Event->label, &(Event->x), &(Event->y), &(Event->z),
            &(Event->otime));

    sprintf(EventStr,
            "EQEVENT:  Label: %s  Loc:  X %.2lf  Y %.2lf  Z %.2lf  OT %.2lf",
            Event->label, Event->x, Event->y, Event->z, Event->otime);
    nll_putmsg(2, EventStr);

    if (istat != 5)
        return (-1);

    return (0);
}

/** function to read station params fom input line */

int GetTime2EQ_Stations(char* in_line) {
    int istat, ierr;
    StationDesc *sta_in;


    /* check number of stations */
    if (NumStationPhases >= MAX_NUM_STATIONS) {
        nll_puterr("ERROR: to many stations, ignoring station.");
        return (0);
    }
    sta_in = Station + NumStationPhases;
    NumStationPhases++;


    /* read source input line */

    sta_in->prob_active = 1.0;
    // 20201015 AJL - added probabilistic data outlier, with probability prob_outlier, phase error is increased by outlier_err_factor
    sta_in->phs[0].prob_outlier = 0.0;
    sta_in->phs[0].outlier_err_factor = 1.0;

    istat = sscanf(in_line, "%s %s %s %lf %s %lf %lf %lf %lf",
            sta_in->label, sta_in->phs[0].label,
            sta_in->phs[0].error_type, &(sta_in->phs[0].error),
            sta_in->phs[0].error_report_type,
            &(sta_in->phs[0].error_report), &(sta_in->prob_active), &(sta_in->phs[0].prob_outlier), &(sta_in->phs[0].outlier_err_factor));

    sprintf(MsgStr,
            "STATION:  %d  Name: %s  Phase: %s  ErrorCalc: %s +/- %.2le  ErrorReport: %s +/- %.2le  ProbActive: %.2lf"
            "  ProbOutlier: %.2lf  OutlierErrorFactor: %.2lf",
            NumStationPhases - 1, sta_in->label, sta_in->phs[0].label,
            sta_in->phs[0].error_type, sta_in->phs[0].error,
            sta_in->phs[0].error_report_type, sta_in->phs[0].error_report,
            sta_in->prob_active,
            sta_in->phs[0].prob_outlier, sta_in->phs[0].outlier_err_factor);
    nll_putmsg(2, MsgStr);

    ierr = 0;

    if (ierr < 0 || (istat != 6 && istat != 7 && istat != 9))
        return (-1);

    return (0);
}

/*** function to read Vp / Vs ratio ***/

int get_vp_vs(char* line1) {
    int istat;

    istat = sscanf(line1, "%lf", &VpVsRatio);

    fprintf(stdout, "EQVPVS: VpVsRatio=%lf\n", VpVsRatio);

    if (istat != 1)
        return (-1);

    return (0);

}

/*** function to read mechanism ***/

int get_mech(char* line1) {
    int istat, ierr;
    char mstr[11];

    istat = sscanf(line1, "%s %lf %lf %lf",
            mech_type, &mech_phi, &mech_del, &mech_lam);

    if (strncmp(mech_type, "double", 6) == 0 ||
            strncmp(mech_type, "DOUBLE", 6) == 0) {
        imech = MECH_DOUBLE;
        strcpy(mstr, "DC");
    } else if (strncmp(mech_type, "iso", 3) == 0 ||
            strncmp(mech_type, "ISO", 3) == 0) {
        imech = MECH_ISOTROPIC;
        strcpy(mstr, "Iso");
    } else if (strncmp(mech_type, "none", 4) == 0 ||
            strncmp(mech_type, "NONE", 4) == 0) {
        imech = MECH_NONE;
        strcpy(mstr, "None");
    } else {
        imech = MECH_NONE;
        strcpy(mstr, "None");
        nll_puterr2("ERROR: unrecognized mechanism:", mech_type);
        return (-1);
    }

    fprintf(stdout, "MECHANISM: %s  Stike=%5.1f Dip=%5.1f Rake=%5.1f\n",
            mech_type, mech_phi, mech_del, mech_lam);

    ierr = 0;
    if (checkRangeDouble("EQMECH", "Stike", mech_phi, 1, 0.0, 1, 360.0) != 0)
        ierr = -1;
    if (checkRangeDouble("EQMECH", "Dip", mech_del, 1, 0.0, 1, 90.0) != 0)
        ierr = -1;
    if (checkRangeDouble("EQMECH", "Rake", mech_lam, 1, -180.0, 1, 180.0) != 0)
        ierr = -1;


    if (ierr < 0 || istat != 4)
        return (-1);


    sprintf(mech_str, "%s(s%.0fd%.0fr%.0f)", mstr, mech_phi, mech_del, mech_lam);

    /* adjust azimuth  for projection azimuth */
    mech_phi = latlon2rectAngle(0, mech_phi);

    mech_phi *= cRPD;
    mech_del *= cRPD;
    mech_lam *= cRPD;

    return (0);

}


/*** function to calculate radiation amplitude */

/* A & R figs 4.20 & 5.5 */

double calc_rad(double ray_vert_ang, double ray_horiz_ang, char orig_wave) {
    double radamp = 1.0;

    /* cal radiation pattern (from Aki & Richards eqs. 4.84 - 4.86) */

    if (orig_wave == 'P') {
        radamp = cos(mech_lam) * sin(mech_del) * pow(sin(ray_vert_ang), 2.0)
                * sin(2.0 * (ray_horiz_ang - mech_phi))
                - cos(mech_lam) * cos(mech_del) * sin(2.0 * ray_vert_ang)
                * cos(ray_horiz_ang - mech_phi)
                + sin(mech_lam) * sin(2.0 * mech_del)
                * (pow(cos(ray_vert_ang), 2.0) - pow(sin(ray_vert_ang)
                * sin(ray_horiz_ang - mech_phi), 2.0))
                + sin(mech_lam) * cos(2.0 * mech_del) * sin(2.0 * ray_vert_ang)
                * sin(ray_horiz_ang - mech_phi)
                ;
    } else if (orig_wave == 'V') {
        radamp = sin(mech_lam) * cos(2.0 * mech_del) * cos(2.0 * ray_vert_ang)
                * sin(ray_horiz_ang - mech_phi)
                - cos(mech_lam) * cos(mech_del) * cos(2.0 * ray_vert_ang)
                * cos(ray_horiz_ang - mech_phi)
                + 0.5 * cos(mech_lam) * sin(mech_del)
                * sin(2.0 * ray_vert_ang) * sin(2.0 * (ray_horiz_ang - mech_phi))
                - 0.5 * sin(mech_lam) * sin(2.0 * mech_del) * sin(2.0 * ray_vert_ang)
                * (1.0 + pow(sin(ray_horiz_ang - mech_phi), 2.0))
                ;
        radamp *= -1.0;
        /* mult by -1 since A & R change dir of + in figs 4.20 & 5.5 */
    } else if (orig_wave == 'H') {
        radamp = cos(mech_lam) * cos(mech_del) * cos(ray_vert_ang)
                * sin(ray_horiz_ang - mech_phi)
                + cos(mech_lam) * sin(mech_del) * sin(ray_vert_ang)
                * cos(2.0 * (ray_horiz_ang - mech_phi))
                + sin(mech_lam) * cos(2.0 * mech_del)
                * cos(ray_vert_ang) * cos(ray_horiz_ang - mech_phi)
                - 0.5 * sin(mech_lam) * sin(2.0 * mech_del) * sin(ray_vert_ang)
                * sin(2.0 * (ray_horiz_ang - mech_phi))
                ;
    } else {
        sprintf(MsgStr,
                "WARNING: Mechanism: unrecognized original wave type: %c", orig_wave);
        nll_putmsg(2, MsgStr);
    }
    return (radamp);
}



/*------------------------------------------------------------/ */
/* Anthony Lomax           | email: lomax@faille.unice.fr     / */
/* UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax / */
/* 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        / */
/* 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        / */
/*------------------------------------------------------------/ */

