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



/*   GridLib.c

        grid library functions

 */

/*-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
161 Allee du Micocoulier, 06370 Mouans-Sartoux, France
tel: +33(0)493752502  e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/


/*
        history:

        ver 01    22SEP1997  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

 */



#define EXTERN_MODE 1

#include "GridLib.h"

// private functions
int _WriteLocation(FILE *fpio, HypoDesc* phypo, ArrivalDesc* parrivals,
        int narrivals, char* filename,
        int iWriteArrivals, int iWriteEndLoc, int iWriteMinimal,
        GridDesc* pgrid, int n_proj, int io_arrival_mode);

/** function to set constants */

void SetConstants(void) {

    MAX_NUM_STATIONS = X_MAX_NUM_STATIONS;
    MAX_NUM_ARRIVALS = X_MAX_NUM_ARRIVALS;
    // AEH/AJL 20080709
    //      if ((Arrival = (ArrivalDesc *) malloc(MAX_NUM_ARRIVALS * sizeof(ArrivalDesc))) == NULL) {
    if ((Arrival = (ArrivalDesc *) calloc(MAX_NUM_ARRIVALS, sizeof (ArrivalDesc))) == NULL) {
        nll_puterr("ERROR: re-allocating Arrival array.");
        exit(EXIT_ERROR_MEMORY);
    }


    strcpy(package_name, PACKAGE);
    strcpy(prog_ver, PVER);
    strcpy(prog_date, PDATE);
    strcpy(prog_copyright, PCOPYRIGHT);
    message_flag = 0;

    // 20171122 AJL  cPI = 4. * atan(1.); /* PI */
    cPI = M_PI; /* PI */
    cRPD = cPI / 180.; /* radians per degree */
    // 20171122 AJL - changed km/deg scaling to be based on sphere with radius 6371, average Earth radius.
    //c111 = 10000.0 / 90.0; /* kilometers per degree */
    c111 = DEG2KM;

    fp_include = NULL; /* set include file ptr */

    /* bookkeeping */
    NumFilesOpen = 0;
    NumGridBufFilesOpen = 0;
    NumGridHdrFilesOpen = 0;
    NumAllocations = 0;

    /* program variables */
    NumQuality2ErrorLevels = 0;
    PhaseFormat = FORMAT_PHASE_1;

    /* set null angles indicator */
    AnglesNULL = SetTakeOffAngles(400.0, 200.0, 0);

    /* set null elliposid */
    /*Ellipsoid3D Ell3NULL = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
    EllipsoidNULL = Ell3NULL;
    Ellipse2D Ell2NULL = {-1.0, -1.0, -1.0};
    EllipseNULL = Ell2NULL;
     */

}

/** function to read control params ***/

int get_control(char* line1) {
    int istat;

    istat = sscanf(line1, "%d %d", &message_flag, &RandomNumSeed);

    if (istat == 1)
        RandomNumSeed = 837465;

    /* display program information */
    DispProgInfo();

    sprintf(MsgStr, "CONTROL:  MessageFlag: %d  RandomNumSeed: %d", message_flag, RandomNumSeed);
    nll_putmsg(3, MsgStr);

    if (checkRangeInt("CONTROL", "MessageFlag", message_flag, 1, -1, 0, 0) != 0)
        return (-1);

    if (istat != 1 && istat != 2)
        return (-1);

    return (0);

}

/** function to read output file name ***/

int get_outfile(char* line1) {

    int istat;

    istat = sscanf(line1, "%s", fn_output);

    sprintf(MsgStr, "OUTPUT FILES: %s.*", fn_output);
    nll_putmsg(3, MsgStr);

    if (istat != 1)
        return (-1);

    return (0);
}

/** function to read include file name and reset input to include file ***/

int GetIncludeFile(char* line1, FILE **fp_io) {

    /* read include file name */

    sscanf(line1, "%s", fn_include);

    sprintf(MsgStr, "Reading from INCLUDE FILE: %s", fn_include);
    nll_putmsg(3, MsgStr);


    /* open include file */

    if ((fp_include = fopen(fn_include, "r")) == NULL) {
        nll_puterr2("ERROR: opening INCLUDE file", fn_include);
        return (-1);
    }
    NumFilesOpen++;


    /* swap file pointers */

    fp_input_save = *fp_io;
    *fp_io = fp_include;


    return (0);
}

/** function to reset input from include file to inlput file ***/

void SwapBackIncludeFP(FILE **fp_io) {

    /* swap file pointers */

    *fp_io = fp_input_save;

    // AJL - 20080709 (valgrind)
    //fp_include = NULL;
    if (fp_include != NULL) {
        fclose(fp_include);
        NumFilesOpen--;
        fp_include = NULL;
    }

    sprintf(MsgStr, "Returning from INCLUDE FILE: %s.*", fn_include);
    nll_putmsg(3, MsgStr);

}

/** function to read source params from input line */

int GetNextSource(char* in_line) {
    int istat;
    SourceDesc *srce_in = NULL;


    /* check number of sources */
    if (NumSources >= MAX_NUM_SOURCES) {
        nll_puterr2("ERROR: to many sources, ignoring source", srce_in->label);
        return (0);
    }

    srce_in = Source + NumSources;
    istat = GetSource(in_line, srce_in, NumSources);
    if (istat < 0)
        return (istat);

    // check if duplicate
    if (FindSource(srce_in->label) != NULL) {
        if (message_flag >= 2) {
            sprintf(MsgStr, "WARNING: duplicated source, ignoring source: %s", srce_in->label);
            nll_putmsg(2, MsgStr);
            return (istat);
        }
    }

    NumSources++;

    return (istat);

}

/** function to read source params fom input line */

int GetSource(char* in_line, SourceDesc *srce_in, int num_sources) {
    int istat, ierr;
    char chr1, chr2, coord_type[MAXLINE];
    double val1, val1a, val1b, val2, val2a, val2b, val3, val4;
    double sign;
    char label[10 * ARRIVAL_LABEL_LEN];


    /* initialize some source fields */
    srce_in->is_coord_xyz = 0;
    srce_in->is_coord_latlon = 0;
    srce_in->otime = 0.0;


    /* read coordinate type */

    istat = sscanf(in_line, "%*s %s", coord_type);

    /* read coordinate type and coordinates */

    if (strncmp(coord_type, "XYZ", 3) == 0) {
        istat = sscanf(in_line, "%s %s %lf %lf %lf %lf",
                label,
                coord_type, &val1, &val2, &val3, &val4);
        strncpy(srce_in->label, label, ARRIVAL_LABEL_LEN - 1);
        srce_in->x = val1;
        srce_in->y = val2;
        srce_in->z = val3 - val4;
        srce_in->is_coord_xyz = 1;
        if (message_flag >= 5) {
            sprintf(MsgStr,
                    "SOURCE: %3d  Name: %s  Loc:  type: %s  X(east) %lg  Y(north) %lg  Z(pos DOWN) %lg",
                    num_sources, srce_in->label,
                    coord_type, srce_in->x, srce_in->y, srce_in->z);
            nll_putmsg(5, MsgStr);
        }
        if (istat != 6)
            return (-1);
    } else if (strcmp(coord_type, "LATLON") == 0) {
        istat = sscanf(in_line, "%s %s %lf %lf %lf %lf", label,
                coord_type, &val1, &val2, &val3, &val4);
        strncpy(srce_in->label, label, ARRIVAL_LABEL_LEN - 1);
        srce_in->dlat = val1;
        srce_in->dlong = val2;
        srce_in->depth = val3 - val4;
        srce_in->is_coord_latlon = 1;
        ierr = 0;
        if (checkRangeDouble("SRCE",
                "Lat", srce_in->dlat, 1, -90.0, 1, 90.0) != 0)
            ierr = -1;
        if (checkRangeDouble("SRCE",
                "Long", srce_in->dlong, 1, -180.0, 1, 180.0) != 0)
            ierr = -1;
        if (message_flag >= 5) {
            sprintf(MsgStr,
                    "SOURCE:  %d  Name: %s  Loc:  type: %s  Lat %lg  Long %lg  Depth %lg",
                    num_sources, srce_in->label,
                    coord_type, srce_in->dlat, srce_in->dlong,
                    srce_in->depth);
            nll_putmsg(5, MsgStr);
        }
        if (ierr < 0 || istat != 6)
            return (-1);
    } else if (strcmp(coord_type, "LATLONDM") == 0) {
        istat = sscanf(in_line,
                "%s %s %lf %lf %c %lf %lf %c %lf %lf",
                label, coord_type, &val1, &val1a, &chr1,
                &val2, &val2a, &chr2, &val3, &val4);
        strncpy(srce_in->label, label, ARRIVAL_LABEL_LEN - 1);
        if ((toupper(chr1) != 'N' && toupper(chr1) != 'S')
                || (toupper(chr2) != 'E' && toupper(chr2) != 'W'))
            return (-1);
        sign = toupper(chr1) == 'N' ? 1.0 : -1.0;
        srce_in->dlat = sign * (val1 + val1a / 60.0);
        sign = toupper(chr2) == 'E' ? 1.0 : -1.0;
        srce_in->dlong = sign * (val2 + val2a / 60.0);
        srce_in->depth = val3 - val4;
        srce_in->is_coord_latlon = 1;
        ierr = 0;
        if (checkRangeDouble("SRCE",
                "Lat", srce_in->dlat, 1, -90.0, 1, 90.0) != 0)
            ierr = -1;
        if (checkRangeDouble("SRCE",
                "Long", srce_in->dlong, 1, -180.0, 1, 180.0) != 0)
            ierr = -1;
        if (message_flag >= 5) {
            sprintf(MsgStr,
                    "SOURCE:  %d  Name: %s  Loc:  type: %s  Lat %lg  Long %lg  Depth %lg",
                    num_sources, srce_in->label,
                    coord_type, srce_in->dlat, srce_in->dlong,
                    srce_in->depth);
            nll_putmsg(5, MsgStr);
        }
        if (ierr < 0 || istat != 10)
            return (-1);
    } else if (strcmp(coord_type, "LATLONDS") == 0) {
        istat = sscanf(in_line,
                "%s %s %lf %lf %lf %c %lf %lf %lf %c %lf %lf",
                label, coord_type,
                &val1, &val1a, &val1b, &chr1,
                &val2, &val2a, &val2b, &chr2, &val3, &val4);
        strncpy(srce_in->label, label, ARRIVAL_LABEL_LEN - 1);
        if ((toupper(chr1) != 'N' && toupper(chr1) != 'S')
                || (toupper(chr2) != 'E' && toupper(chr2) != 'W'))
            return (-1);
        sign = toupper(chr1) == 'N' ? 1.0 : -1.0;
        srce_in->dlat = sign * (val1 + (val1a + val1b / 60.0) / 60.0);
        sign = toupper(chr2) == 'E' ? 1.0 : -1.0;
        srce_in->dlong = sign * (val2 + (val2a + val2b / 60.0) / 60.0);
        srce_in->depth = val3 - val4;
        srce_in->is_coord_latlon = 1;
        ierr = 0;
        if (checkRangeDouble("SRCE",
                "Lat", srce_in->dlat, 1, -90.0, 1, 90.0) != 0)
            ierr = -1;
        if (checkRangeDouble("SRCE",
                "Long", srce_in->dlong, 1, -180.0, 1, 180.0) != 0)
            ierr = -1;
        if (message_flag >= 5) {
            sprintf(MsgStr,
                    "SOURCE:  %d  Name: %s  Loc:  type: %s  Lat %lg  Long %lg  Depth %lg",
                    num_sources, srce_in->label,
                    coord_type, srce_in->dlat, srce_in->dlong,
                    srce_in->depth);
            nll_putmsg(5, MsgStr);
        }
        if (ierr < 0 || istat != 12)
            return (-1);
    } else {
        nll_puterr2("ERROR: unrecognized coordinate type:", in_line);
        return (-1);
    }




    return (0);
}

/** function to find source from label */

SourceDesc* FindSource(char* label) {

    int len = strlen(label);

    int nsrce;
    int len2;
    for (nsrce = 0; nsrce < NumSources; nsrce++) {
        len2 = strlen((Source + nsrce)->label);
        // 20100512 AJL
        if (len2 != len)
            continue;
        // 20100409 AJL if(strncmp((Source + nsrce)->label, label, ARRIVAL_LABEL_LEN) == 0)
        if (strncmp((Source + nsrce)->label, label, len) == 0)
            return (Source + nsrce);
    }
    return (NULL);
}

/** function to convert map projection string to map transformation parameter string by removing parameter tags
 *
 * Example: converts
 * TRANSFORM  LAMBERT RefEllipsoid Clarke-1880  LatOrig 40.780400  LongOrig 15.415500  FirstStdParal 43.199300  SecondStdParal 44.996100  RotCW 0.000000
 * to
 * LAMBERT Clarke-1880  40.780400  15.415500  43.199300  44.996100  0.000000
 *
 ***/

char* projection_str2transform_str(char* trans_str, char* proj_str) {

    char *proj_ptr;
    char *trans_ptr;

    proj_ptr = proj_str;
    trans_ptr = trans_str;
    while (*proj_ptr != '\0') {
        // skip tag
        while (*proj_ptr != ' ' && *proj_ptr != '\0') {
            proj_ptr++;
        }
        // skip space
        while (*proj_ptr == ' ' && *proj_ptr != '\0') {
            proj_ptr++;
        }
        // copy parameter value
        while (*proj_ptr != ' ' && *proj_ptr != '\0') {
            *trans_ptr = *proj_ptr;
            trans_ptr++;
            proj_ptr++;
        }
        // copy space
        while (*proj_ptr == ' ' && *proj_ptr != '\0') {
            *trans_ptr = *proj_ptr;
            trans_ptr++;
            proj_ptr++;
        }
    }

    *trans_ptr = '\0';
    return (trans_str);

}

/** function to read map transformation parameters from input line ***/

int get_transform(int n_proj, char* in_line) {
    int istat, ierr;
    double angle;

    // SDC
    double dlt1, dlt2, del, r, bc;


    map_itype[n_proj] = MAP_TRANS_UNDEF;
    GeometryMode = MODE_RECT;

    /* read transform input line */

    sscanf(in_line, "%s", map_trans_type[n_proj]);

    if (strcmp(map_trans_type[n_proj], "GLOBAL") == 0) {

        // mode
        GeometryMode = MODE_GLOBAL;

        map_itype[n_proj] = MAP_TRANS_GLOBAL;
        istat = sscanf(in_line, "%s", map_trans_type[n_proj]);

        angle = 0.0;
        map_cosang[n_proj] = cos(angle);
        map_sinang[n_proj] = sin(angle);

        sprintf(MapProjStr[n_proj],
                "TRANSFORM  %s",
                map_trans_type[n_proj]);
        nll_putmsg(3, MapProjStr[n_proj]);

        ierr = 0;
        if (ierr < 0 || istat != 1) {
            nll_puterr("ERROR: reading GLOBAL transformation parameters");
            return (-1);
        }


    } else if (strcmp(map_trans_type[n_proj], "NONE") == 0) {

        map_itype[n_proj] = MAP_TRANS_NONE;
        istat = sscanf(in_line, "%s", map_trans_type[n_proj]);

        angle = 0.0;
        map_cosang[n_proj] = cos(angle);
        map_sinang[n_proj] = sin(angle);

        sprintf(MapProjStr[n_proj],
                "TRANSFORM  %s",
                map_trans_type[n_proj]);
        nll_putmsg(3, MapProjStr[n_proj]);

        ierr = 0;
        if (ierr < 0 || istat != 1) {
            nll_puterr("ERROR: reading NONE transformation parameters");
            return (-1);
        }


    } else if (strcmp(map_trans_type[n_proj], "SIMPLE") == 0) {

        map_itype[n_proj] = MAP_TRANS_SIMPLE;
        istat = sscanf(in_line, "%s %lf %lf %lf",
                map_trans_type[n_proj], &map_orig_lat[n_proj], &map_orig_long[n_proj],
                &map_rot[n_proj]);

        angle = -cRPD * map_rot[n_proj];
        map_cosang[n_proj] = cos(angle);
        map_sinang[n_proj] = sin(angle);

        sprintf(MapProjStr[n_proj],
                "TRANSFORM  %s LatOrig %lf  LongOrig %lf  RotCW %lf",
                map_trans_type[n_proj], map_orig_lat[n_proj], map_orig_long[n_proj], map_rot[n_proj]);
        nll_putmsg(3, MapProjStr[n_proj]);

        ierr = 0;
        if (checkRangeDouble("TRANS",
                "LatOrig", map_orig_lat[n_proj], 1, -90.0, 1, 90.0) != 0)
            ierr = -1;
        if (checkRangeDouble("TRANS",
                "LongOrig", map_orig_long[n_proj], 1, -180.0, 1, 180.0) != 0)
            ierr = -1;
        if (checkRangeDouble("TRANS",
                "RotCW", map_rot[n_proj], 1, -360.0, 1, 360.0) != 0)
            ierr = -1;
        if (ierr < 0 || istat != 4) {
            nll_puterr("ERROR: reading SIMPLE transformation parameters");
            return (-1);
        }


    } else if (strcmp(map_trans_type[n_proj], "SDC") == 0) {

        map_itype[n_proj] = MAP_TRANS_SDC;
        istat = sscanf(in_line, "%s %lf %lf %lf",
                map_trans_type[n_proj], &map_orig_lat[n_proj], &map_orig_long[n_proj],
                &map_rot[n_proj]);

        angle = -cRPD * map_rot[n_proj];
        map_cosang[n_proj] = cos(angle);
        map_sinang[n_proj] = sin(angle);


        sprintf(MapProjStr[n_proj],
                "TRANSFORM  %s LatOrig %lf  LongOrig %lf  RotCW %lf",
                map_trans_type[n_proj], map_orig_lat[n_proj], map_orig_long[n_proj], map_rot[n_proj]);
        nll_putmsg(3, MapProjStr[n_proj]);

        ierr = 0;
        if (checkRangeDouble("TRANS",
                "LatOrig", map_orig_lat[n_proj], 1, -90.0, 1, 90.0) != 0)
            ierr = -1;
        if (checkRangeDouble("TRANS",
                "LongOrig", map_orig_long[n_proj], 1, -180.0, 1, 180.0) != 0)
            ierr = -1;
        if (checkRangeDouble("TRANS",
                "RotCW", map_rot[n_proj], 1, -360.0, 1, 360.0) != 0)
            ierr = -1;
        if (ierr < 0 || istat != 4) {
            nll_puterr("ERROR: reading SDC transformation parameters");
            return (-1);
        }

        // SDC initialization
        //  conversion factor for latitude
        dlt1 = atan(MAP_TRANS_SDC_DRLT * tan(map_orig_lat[n_proj] * DE2RA));
        dlt2 = atan(MAP_TRANS_SDC_DRLT * tan((map_orig_lat[n_proj] + 1.0) * DE2RA));
        del = dlt2 - dlt1;
        r = ERAD * (1.0 - pow(sin(dlt1), 2) * FLATTENING);
        map_sdc_xltkm[n_proj] = r * del;
        //  conversion factor for longitude
        del = acos(1.0 - (1.0 - cos(DE2RA)) * pow(cos(dlt1), 2));
        bc = r * del;
        map_sdc_xlnkm[n_proj] = bc / cos(dlt1);



    } else if (strcmp(map_trans_type[n_proj], "LAMBERT") == 0) {

        map_itype[n_proj] = MAP_TRANS_LAMBERT;
        istat = sscanf(in_line, "%s %s %lf %lf %lf %lf %lf",
                map_trans_type[n_proj], map_ref_ellipsoid[n_proj],
                &map_orig_lat[n_proj], &map_orig_long[n_proj],
                &map_lambert_1st_std_paral[n_proj], &map_lambert_2nd_std_paral[n_proj],
                &map_rot[n_proj]);

        ierr = 0;
        if (checkRangeDouble("TRANS",
                "LatOrig", map_orig_lat[n_proj], 1, -90.0, 1, 90.0) != 0)
            ierr = -1;
        if (checkRangeDouble("TRANS",
                "LongOrig", map_orig_long[n_proj], 1, -180.0, 1, 180.0) != 0)
            ierr = -1;
        if (checkRangeDouble("TRANS",
                "FirstStdParal", map_lambert_1st_std_paral[n_proj], 1, -90.0, 1, 90.0) != 0)
            ierr = -1;
        if (checkRangeDouble("TRANS",
                "SecondStdParal", map_lambert_2nd_std_paral[n_proj], 1, -90.0, 1, 90.0) != 0)
            ierr = -1;
        if (checkRangeDouble("TRANS",
                "RotCW", map_rot[n_proj], 1, -360.0, 1, 360.0) != 0)
            ierr = -1;

        angle = -cRPD * map_rot[n_proj];
        map_cosang[n_proj] = cos(angle);
        map_sinang[n_proj] = sin(angle);

        /* initialize GMT projection values */
        if (map_setup_proxy(n_proj, map_ref_ellipsoid[n_proj]) < 0) {
            nll_puterr(
                    "ERROR: initializing general transformation parameters, RefEllipsoid may be invalid");
            return (-1);
        }

        /* initialize lambert projection */
        vlamb(n_proj, map_orig_long[n_proj], map_orig_lat[n_proj],
                map_lambert_1st_std_paral[n_proj], map_lambert_2nd_std_paral[n_proj]);

        sprintf(MapProjStr[n_proj],
                "TRANSFORM  %s RefEllipsoid %s  LatOrig %lf  LongOrig %lf  FirstStdParal %lf  SecondStdParal %lf  RotCW %lf",
                map_trans_type[n_proj], map_ref_ellipsoid[n_proj],
                map_orig_lat[n_proj], map_orig_long[n_proj],
                map_lambert_1st_std_paral[n_proj], map_lambert_2nd_std_paral[n_proj],
                map_rot[n_proj]);
        nll_putmsg(3, MapProjStr[n_proj]);

        if (ierr < 0 || istat != 7) {
            nll_puterr("ERROR: reading LAMBERT transformation parameters");
            return (-1);
        }

    } else if (strcmp(map_trans_type[n_proj], "TRANS_MERC") == 0) {

        // 20160725 AJL - added TRANSVERSE MERCATOR PROJECTION (TM)

        map_itype[n_proj] = MAP_TRANS_TM;
        int use_false_easting = 0;
        istat = sscanf(in_line, "%s %s %lf %lf %lf %d",
                map_trans_type[n_proj], map_ref_ellipsoid[n_proj],
                &map_orig_lat[n_proj], &map_orig_long[n_proj],
                &map_rot[n_proj], &use_false_easting);

        ierr = 0;
        if (checkRangeDouble("TRANS",
                "LatOrig", map_orig_lat[n_proj], 1, -90.0, 1, 90.0) != 0)
            ierr = -1;
        if (checkRangeDouble("TRANS",
                "LongOrig", map_orig_long[n_proj], 1, -180.0, 1, 180.0) != 0)
            ierr = -1;
        if (checkRangeDouble("TRANS",
                "RotCW", map_rot[n_proj], 1, -360.0, 1, 360.0) != 0)
            ierr = -1;

        angle = -cRPD * map_rot[n_proj];
        map_cosang[n_proj] = cos(angle);
        map_sinang[n_proj] = sin(angle);

        // initialize GMT projection values
        if (map_setup_proxy(n_proj, map_ref_ellipsoid[n_proj]) < 0) {
            nll_puterr(
                    "ERROR: initializing general transformation parameters, RefEllipsoid may be invalid");
            return (-1);
        }

        // initialize projection
        vtm(n_proj, map_orig_long[n_proj], map_orig_lat[n_proj], use_false_easting);

        sprintf(MapProjStr[n_proj],
                "TRANSFORM  %s RefEllipsoid %s  LatOrig %lf  LongOrig %lf  RotCW %lf  UseFalseEasting %d",
                map_trans_type[n_proj], map_ref_ellipsoid[n_proj],
                map_orig_lat[n_proj], map_orig_long[n_proj],
                map_rot[n_proj], use_false_easting);
        nll_putmsg(3, MapProjStr[n_proj]);

        if (ierr < 0 || (istat != 5 && istat != 6)) {
            nll_puterr("ERROR: reading TRANS_MERC transformation parameters");
            return (-1);
        }

    } else if (strcmp(map_trans_type[n_proj], "AZIMUTHAL_EQUIDIST") == 0) {

        // 20171120 AJL - added AZIMUTHAL EQUIDISTANT PROJECTION

        map_itype[n_proj] = MAP_TRANS_AZ_EQUID;
        istat = sscanf(in_line, "%s %s %lf %lf %lf",
                map_trans_type[n_proj], map_ref_ellipsoid[n_proj],
                &map_orig_lat[n_proj], &map_orig_long[n_proj],
                &map_rot[n_proj]);

        ierr = 0;
        if (checkRangeDouble("TRANS",
                "LatOrig", map_orig_lat[n_proj], 1, -90.0, 1, 90.0) != 0)
            ierr = -1;
        if (checkRangeDouble("TRANS",
                "LongOrig", map_orig_long[n_proj], 1, -180.0, 1, 180.0) != 0)
            ierr = -1;
        if (checkRangeDouble("TRANS",
                "RotCW", map_rot[n_proj], 1, -360.0, 1, 360.0) != 0)
            ierr = -1;

        angle = -cRPD * map_rot[n_proj];
        map_cosang[n_proj] = cos(angle);
        map_sinang[n_proj] = sin(angle);

        // initialize GMT projection values
        if (map_setup_proxy(n_proj, map_ref_ellipsoid[n_proj]) < 0) {
            nll_puterr(
                    "ERROR: initializing general transformation parameters, RefEllipsoid may be invalid");
            return (-1);
        }

        // initialize projection
        vazeqdist(n_proj, map_orig_long[n_proj], map_orig_lat[n_proj]);

        sprintf(MapProjStr[n_proj],
                "TRANSFORM  %s RefEllipsoid %s  LatOrig %lf  LongOrig %lf  RotCW %lf",
                map_trans_type[n_proj], map_ref_ellipsoid[n_proj],
                map_orig_lat[n_proj], map_orig_long[n_proj],
                map_rot[n_proj]);
        nll_putmsg(3, MapProjStr[n_proj]);

        if (ierr < 0 || istat != 5) {
            nll_puterr("ERROR: reading AZIMUTHAL_EQUIDIST transformation parameters");
            return (-1);
        }

    } else {

        nll_puterr("ERROR: unrecognized map transformation type");
        return (-1);
    }

    return (0);
}

/** function to read grid params */

int get_grid(char* input_line) {
    int istat, ierr = 0;

    istat = sscanf(input_line, "%d %d %d %lf %lf %lf %lf %lf %lf %s",
            &(grid_in.numx), &(grid_in.numy), &(grid_in.numz),
            &(grid_in.origx), &(grid_in.origy), &(grid_in.origz),
            &(grid_in.dx), &(grid_in.dy), &(grid_in.dz), grid_in.chr_type);

    grid_in.iSwapBytes = 0;

    convert_grid_type(&grid_in, 1);
    if (message_flag >= 2)
        display_grid_param(&grid_in);

    if (checkRangeInt("GRID", "xNum", grid_in.numx, 1, 2, 0, 0) != 0)
        ierr = -1;
    if (checkRangeInt("GRID", "yNum", grid_in.numy, 1, 2, 0, 0) != 0)
        ierr = -1;
    if (checkRangeInt("GRID", "zNum", grid_in.numz, 1, 2, 0, 0) != 0)
        ierr = -1;

    if (ierr < 0)
        return (-1);

    if (istat != 10)
        return (-1);

    return (0);
}

/** function to convert grid type from char to numeric or v.v. */

int convert_grid_type(GridDesc* pgrid, int char2numeric) {

    int i;


    int ntypes = 20;
    //
    static char char_types[][32] = {
        "VELOCITY", "VELOCITY_METERS",
        "SLOWNESS", "SLOW_LEN", "VEL2", "SLOW2", "SLOW2_METERS",
        "TIME", "TIME2D",
        "ANGLE", "ANGLE2D",
        "INCLINATION", "INCLINATION2D",
        "PROB_DENSITY", "MISFIT", "LIKELIHOOD",
        "DEPTH", "LENGTH", "COULOMB", "SSST_TIMECORR"
    };
    //
    static int types[] = {
        GRID_VELOCITY, GRID_VELOCITY_METERS,
        GRID_SLOWNESS, GRID_SLOW_LEN, GRID_VEL2, GRID_SLOW2, GRID_SLOW2_METERS,
        GRID_TIME, GRID_TIME_2D,
        GRID_ANGLE, GRID_ANGLE_2D,
        GRID_INCLINATION, GRID_INCLINATION_2D,
        GRID_PROB_DENSITY, GRID_MISFIT, GRID_LIKELIHOOD,
        GRID_DEPTH, GRID_LENGTH, GRID_COULOMB, GRID_SSST_TIMECORR
    };



    /* check grid type */

    if (char2numeric) {
        for (i = 0; i < ntypes; i++) {
            if (strcmp(pgrid->chr_type, char_types[i]) == 0) {
                pgrid->type = types[i];
                return (pgrid->type);
            }
        }
        pgrid->type = GRID_UNDEF;
        nll_puterr2("WARNING: unrecognized grid type", pgrid->chr_type);
        return (pgrid->type);

    } else {
        for (i = 0; i < ntypes; i++) {
            if (pgrid->type == types[i]) {
                strcpy(pgrid->chr_type, char_types[i]);
                return (pgrid->type);
            }
        }
        nll_puterr("WARNING: unrecognized grid type code");
        return (pgrid->type);
    }

}

/** function to display grid params */

int display_grid_param(GridDesc* pgrid) {

    fprintf(stdout,
            "GRID: {x, y, z}\n  Num: {%d, %d, %d}\n  Orig: {%lg, %lg, %lg}\n  LenSide: {%lg, %lg, %lg}\n",
            pgrid->numx, pgrid->numy, pgrid->numz,
            pgrid->origx, pgrid->origy, pgrid->origz,
            pgrid->dx, pgrid->dy, pgrid->dz);
    fprintf(stdout, "  Type: %s\n", pgrid->chr_type);

    return (0);
}

/** function to determine if a grid is type Cascading 3D grid
 *
 * 20161019 AJL - added
 */

int isCascadingGrid(GridDesc* pgrid) {

    return (pgrid->flagGridCascading == IS_CASCADING);

}

/** function set a grid as type Cascading 3D grid
 *
 * 20161020 AJL - added
 */

void setCascadingGrid(GridDesc* pgrid) {

    pgrid->flagGridCascading = IS_CASCADING;

    // flag dynamic arrays as uninitialized
    pgrid->gridDesc_Cascading.xyz_scale = NULL;
    pgrid->gridDesc_Cascading.zindex = NULL;

}

/** function to write grid buffer and header to disk ***/

int WriteGrid3dBuf(GridDesc* pgrid, SourceDesc* psrce, char* filename, char* file_type) {


    int istat;
    FILE *fpio;
    char fname[FILENAME_MAX];


    /* write buffer file */

    if (file_type != NULL) {
        sprintf(fname, "%s.%s.buf", filename, file_type);
    } else {
        sprintf(fname, "%s.buf", filename);
    }
    if ((fpio = fopen(fname, "w")) == NULL) {
        nll_puterr("ERROR: opening buffer output file.");
        return (-1);
    }
    NumFilesOpen++;

    if (fwrite((char *) pgrid->buffer, pgrid->buffer_size, 1, fpio) != 1) {
        nll_puterr("ERROR: writing grid buffer output file.");
        return (-1);
    }

    fclose(fpio);
    NumFilesOpen--;


    /* write header file */

    istat = WriteGrid3dHdr(pgrid, psrce, filename, file_type);

    return (istat);
}

/** function to write grid header to disk ***/

int WriteGrid3dHdr(GridDesc* pgrid, SourceDesc* psrce,
        char* filename, char* file_type) {

    FILE *fpio;
    char fname[FILENAME_MAX];


    /* write header file */

    if (file_type != NULL) {
        sprintf(fname, "%s.%s.hdr", filename, file_type);
    } else {
        sprintf(fname, "%s.hdr", filename);
    }

    if ((fpio = fopen(fname, "w")) == NULL) {
        nll_puterr("ERROR: opening grid output header file.");
        return (-1);
    }
    NumFilesOpen++;

    fprintf(fpio, "%d %d %d  %lf %lf %lf  %lf %lf %lf %s",
            pgrid->numx, pgrid->numy, pgrid->numz,
            pgrid->origx, pgrid->origy, pgrid->origz,
            pgrid->dx, pgrid->dy, pgrid->dz, pgrid->chr_type);

#ifdef GRID_FLOAT_TYPE_DOUBLE
    fprintf(fpio, " DOUBLE\n");
#else
    fprintf(fpio, " FLOAT\n");
#endif

    if (pgrid->type == GRID_TIME || pgrid->type == GRID_TIME_2D
            || pgrid->type == GRID_ANGLE || pgrid->type == GRID_ANGLE_2D
            || pgrid->type == GRID_INCLINATION || pgrid->type == GRID_INCLINATION_2D
            )
        fprintf(fpio, "%s %lf %lf %lf\n",
            psrce->label, psrce->x, psrce->y, psrce->z);

    fprintf(fpio, "%s\n", MapProjStr[0]);

    // write extra cascading grid header lines
    if (isCascadingGrid(pgrid)) {
        fprintf(fpio, "CASCADING_GRID %d ", pgrid->gridDesc_Cascading.num_z_merge_depths);
        for (int n = 0; n < pgrid->gridDesc_Cascading.num_z_merge_depths; n++) {
            fprintf(fpio, "%f,", pgrid->gridDesc_Cascading.z_merge_depths[n]);
        }
    }
    fprintf(fpio, "\n");


    fclose(fpio);
    NumFilesOpen--;

    return (0);
}

//#define DEBUG_CASC

/** function to allocate buffer for Cascading 3D grid
 *
 * 20161019 AJL - added
 */
void* AllocateGrid_Cascading(GridDesc* pgrid, int allocate_buffer) {

    // free any existing allocations
    if (allocate_buffer) {
        FreeGrid(pgrid);
    } else {
        FreeGrid_Cascading(pgrid);
    }

    // allocate cascading grid description index arrays
    pgrid->gridDesc_Cascading.zindex = (void *) malloc((size_t) (pgrid->numz * sizeof (double)));
    NumAllocations++;
    pgrid->gridDesc_Cascading.xyz_scale = (void *) malloc((size_t) (pgrid->numz * sizeof (int)));
    NumAllocations++;

    // set indices and determine grid memory size
    double reg_grid_depth = pgrid->origz;
    int reg_grid_z_index = 0;
    int reg_grid_z_index_below_last_merge_depth = 0;
    double casc_grid_z_index = 0.0;
    int merge_factor = 1;
    long size_casc_grid = 0;
    long size_reg_grid = 0;
    // loop over specified merge depths
    int below_deepest_merge_depth = 0;
    for (int n = 0; n <= pgrid->gridDesc_Cascading.num_z_merge_depths; n++) {
#ifdef DEBUG_CASC
        printf("DEBUG: ndep %d depth=%f ===========================================\n", n,
                (n < pgrid->gridDesc_Cascading.num_z_merge_depths) ? pgrid->gridDesc_Cascading.z_merge_depths[n] : 999);
#endif
        if (n == pgrid->gridDesc_Cascading.num_z_merge_depths) {
            below_deepest_merge_depth = 1;
        }
        // set zindex and xyz_scale while regular grid depth < merge depth
        while (reg_grid_z_index < pgrid->numz // above bottom of regular grid
                && (below_deepest_merge_depth // below deepest specified merge depth
                || (!below_deepest_merge_depth && reg_grid_depth < pgrid->gridDesc_Cascading.z_merge_depths[n]) // above current merge depth
                || reg_grid_z_index_below_last_merge_depth % merge_factor != 0 // not at bottom of new cascading grid cell
                )) {
            if (reg_grid_z_index_below_last_merge_depth % merge_factor == 0) { // at top of new cascading grid depth level, increment grid memory size
                //size_casc_grid += (pgrid->numx / merge_factor) * (pgrid->numy / merge_factor);  // truncate in x, y
                //size_casc_grid += (int) ((double) pgrid->numx / (double) merge_factor) * (int) ((double) pgrid->numy / (double) merge_factor);  // allow fraction of full cascading cell x, y
                int nx_use = 1 + (pgrid->numx - 1) / merge_factor;
                nx_use += (pgrid->numx - 1) % merge_factor > 0 ? 1 : 0;
                //nx_use += pgrid->numx % merge_factor != 1 ? 1 : 0;
                int ny_use = 1 + (pgrid->numy - 1) / merge_factor;
                ny_use += (pgrid->numy - 1) % merge_factor > 0 ? 1 : 0;
                //ny_use += pgrid->numy % merge_factor != 1 ? 1 : 0;
                size_casc_grid += nx_use * ny_use; // allow fraction of full cascading cell x, y
            }
            size_reg_grid += pgrid->numx * pgrid->numy;
#ifdef DEBUG_CASC
            printf("DEBUG: ireg %d, dep_reg %.1f, ds_casc %.1f, icasc %.2f, size_casc %ld/%ld [%.2f]\n",
                    reg_grid_z_index, reg_grid_depth, pgrid->dz * (double) merge_factor, casc_grid_z_index,
                    size_casc_grid, size_reg_grid, (double) size_casc_grid / (double) size_reg_grid);
#endif
            pgrid->gridDesc_Cascading.zindex[reg_grid_z_index] = casc_grid_z_index;
            casc_grid_z_index += (1.0 / (double) merge_factor);
            pgrid->gridDesc_Cascading.xyz_scale[reg_grid_z_index] = merge_factor;
            reg_grid_z_index++;
            reg_grid_z_index_below_last_merge_depth++;
            reg_grid_depth += pgrid->dz;
        }
        if (!below_deepest_merge_depth && reg_grid_z_index >= pgrid->numz) {
            sprintf(MsgStr, "WARNING: AllocateGrid_Cascading: z merge depth: %f below grid bottom: %f",
                    pgrid->gridDesc_Cascading.z_merge_depths[n], pgrid->origz + (double) (pgrid->numz - 1) * pgrid->dz);
            nll_puterr(MsgStr);
            //return (NULL);
        }
        reg_grid_z_index_below_last_merge_depth = 0;
        merge_factor *= 2; // creates a new cascade level
    }
    //pgrid->gridDesc_Cascading.num_z_cascading = pgrid->gridDesc_Cascading.zindex[pgrid->numz - 1];
    pgrid->buffer_size = (size_t) (size_casc_grid * sizeof (GRID_FLOAT_TYPE));
    // DEBUG pgrid->buffer_size -= 128  * sizeof (GRID_FLOAT_TYPE);
    if (allocate_buffer) {
        pgrid->buffer = (void *) malloc(pgrid->buffer_size);
        if (pgrid->buffer != NULL)
            NumAllocations++;
    }
#ifdef DEBUG_CASC
    printf("DEBUG: byte size_casc %ld/%ld [%.2f], val size_casc %ld/%ld\n",
            size_casc_grid * sizeof (GRID_FLOAT_TYPE), size_reg_grid * sizeof (GRID_FLOAT_TYPE), (double) size_casc_grid / (double) size_reg_grid,
            size_casc_grid, size_reg_grid);
#endif

    return (pgrid->buffer);
}

/** function to create array for accessing Cascading 3D grid
 *
 * 20161019 AJL - added
 */
void*** CreateGridArray_Cascading(GridDesc* pgrid) {

#ifdef DEBUG_CASC
    for (int i = 0; i < pgrid->numz; i++) { // DEBUG
        printf("iz_reg %d, zindex[i] %d, xyz_scale[i] %d\n",
                i, pgrid->gridDesc_Cascading.zindex[i], pgrid->gridDesc_Cascading.xyz_scale[i]);
    }
#endif

    GRID_FLOAT_TYPE ***garray;
    if ((garray = (GRID_FLOAT_TYPE ***) malloc((size_t) pgrid->numx * sizeof (GRID_FLOAT_TYPE **))) == NULL)
        return ((void***) garray);
    NumAllocations++;

    int ix_casc, iy_casc, iz_reg, numz;
    int current_numz;
    int xyz_scale;
    GRID_FLOAT_TYPE *buf_ptr = (GRID_FLOAT_TYPE *) pgrid->buffer;
    long bfactor = 1;
    if (buf_ptr == NULL) { // buffer not created, arrays used to fseek in direct read of buffer on disk in ReadGrid3dValue()
        bfactor = 2; // conversion from pointer addressing to absolute byte addressing
        buf_ptr = 0;
    }
#ifdef DEBUG_CASC
    printf("garray begin:  buf_ptr %ld  buf_ptr_end %ld\n", (long) buf_ptr, (buf_ptr - (GRID_FLOAT_TYPE *) pgrid->buffer));
#endif
    for (ix_casc = 0; ix_casc < pgrid->numx; ix_casc++) {
        if ((garray[ix_casc] = (GRID_FLOAT_TYPE **) malloc((size_t) pgrid->numy * sizeof (GRID_FLOAT_TYPE *))) == NULL)
            return (NULL);
        NumAllocations++;
        for (iy_casc = 0; iy_casc < pgrid->numy; iy_casc++) {
            current_numz = -1;
            // calculate num z cascading at this x,y
            numz = 0;
            for (iz_reg = 0; iz_reg < pgrid->numz; iz_reg++) { // loop over nominal regular grid z levels
                if (pgrid->gridDesc_Cascading.zindex[iz_reg] == current_numz) { // still in current cascading grid z level
                    continue;
                }
                // into new cascading grid z level
                current_numz = pgrid->gridDesc_Cascading.zindex[iz_reg];
                // check if xyz cascading cell exists at this ix_casc,iy_casc,numz
                xyz_scale = pgrid->gridDesc_Cascading.xyz_scale[iz_reg];
                if (ix_casc * xyz_scale < pgrid->numx + xyz_scale - 1 && iy_casc * xyz_scale < pgrid->numy + xyz_scale - 1) {
                    //if (ix_casc <= pgrid->numx / xyz_scale && iy_casc <= pgrid->numy / xyz_scale) {
                    numz++; // increment count of z levels for this x,y, cascading
                } else {
                    break; // reached bottom of z levels for this x,y cascading
                }
            }
            garray[ix_casc][iy_casc] = buf_ptr;
            //buf_ptr += numz * sizeof (GRID_FLOAT_TYPE);
            buf_ptr += numz * bfactor;
            if (buf_ptr - (GRID_FLOAT_TYPE *) pgrid->buffer > 2 * pgrid->buffer_size / sizeof (GRID_FLOAT_TYPE)) {
                sprintf(MsgStr, "ERROR: CreateGridArray_Cascading: buf_ptr > buffer_size: x%d y%d numz%d (offset %ld buf_size %ld diff %ld) in: %s",
                        ix_casc, iy_casc, numz, (buf_ptr - (GRID_FLOAT_TYPE *) pgrid->buffer), 2 * pgrid->buffer_size / sizeof (GRID_FLOAT_TYPE), (buf_ptr - (GRID_FLOAT_TYPE *) pgrid->buffer) - 2 * pgrid->buffer_size / sizeof (GRID_FLOAT_TYPE), pgrid->title);
                nll_puterr(MsgStr);
            }
#ifdef DEBUG_CASC
            //if ((ix_casc == 0 || ix_casc == pgrid->numx - 1) && (iy_casc == 0 || iy_casc == pgrid->numy - 1))
            if (ix_casc < 4 && iy_casc < 4 && iz_reg < 4)
                printf("garray [ix_casc %d, iy_casc %d]  buf_ptr %ld  buf_ptr_end %ld  numz %d\n",
                    ix_casc, iy_casc, (long) buf_ptr, (buf_ptr - (GRID_FLOAT_TYPE *) pgrid->buffer), numz);
#endif
        }
    }
#ifdef DEBUG_CASC
    printf("garray end:  buf_ptr %ld  buf_ptr_end %ld  buf_size %ld\n", (long) buf_ptr, (buf_ptr - (GRID_FLOAT_TYPE *) pgrid->buffer), pgrid->buffer_size / sizeof (GRID_FLOAT_TYPE));
#endif

    pgrid->array = (void***) garray;

    return ((void***) garray);
}

/** function to allocate buffer for 3D grid ***/

void* AllocateGrid(GridDesc * pgrid) {

    if (isCascadingGrid(pgrid)) {
        AllocateGrid_Cascading(pgrid, 1);
        return (pgrid->buffer);
    }

    pgrid->buffer_size = (size_t) (pgrid->numx * pgrid->numy * pgrid->numz * sizeof (GRID_FLOAT_TYPE));
    pgrid->buffer = (void *) malloc(pgrid->buffer_size);
    if (pgrid->buffer != NULL)
        NumAllocations++;

    return (pgrid->buffer);
}

/** function to free buffer for 3D grid ***/

void FreeGrid_Cascading(GridDesc * pgrid) {

    // 20170207 AJL - z_merge_depths moved to fixed array to ease memory management
    /*if (pgrid->gridDesc_Cascading.z_merge_depths != NULL) {
        free(pgrid->gridDesc_Cascading.z_merge_depths);
        pgrid->gridDesc_Cascading.z_merge_depths = NULL;
        NumAllocations--;
    }
    pgrid->gridDesc_Cascading.z_merge_depths = NULL;*/
    //printf("DEBUG: FreeGrid_Cascading: pgrid->gridDesc_Cascading.zindex %ld\n", (long) pgrid->gridDesc_Cascading.zindex);
    if (pgrid->gridDesc_Cascading.zindex != NULL) {
        free(pgrid->gridDesc_Cascading.zindex);
        pgrid->gridDesc_Cascading.zindex = NULL;
        NumAllocations--;
    }
    pgrid->gridDesc_Cascading.zindex = NULL;
    if (pgrid->gridDesc_Cascading.xyz_scale != NULL) {
        free(pgrid->gridDesc_Cascading.xyz_scale);
        pgrid->gridDesc_Cascading.xyz_scale = NULL;
        NumAllocations--;
    }
    pgrid->gridDesc_Cascading.xyz_scale = NULL;
}

/** function to free buffer for 3D grid ***/

void FreeGrid(GridDesc * pgrid) {

    if (isCascadingGrid(pgrid)) {
        FreeGrid_Cascading(pgrid);
    }

    if (pgrid->buffer != NULL) {
        free(pgrid->buffer);
        pgrid->buffer = NULL;
        NumAllocations--;
    }
}

/** function to initialize buffer for 3D grid ***/

int InitializeGrid(GridDesc* pgrid, GRID_FLOAT_TYPE init_value) {

    GRID_FLOAT_TYPE *gbuf;

    gbuf = (GRID_FLOAT_TYPE *) pgrid->buffer + pgrid->numx * pgrid->numy * pgrid->numz;

    while (gbuf-- > (GRID_FLOAT_TYPE *) pgrid->buffer)
        *gbuf = init_value;

    return (0);
}

/** function to create array for accessing 3D grid ***/

void*** CreateGridArray(GridDesc * pgrid) {

    if (isCascadingGrid(pgrid)) {
        return (CreateGridArray_Cascading(pgrid));
    }

    int ix, iy, numyz;
    GRID_FLOAT_TYPE ***garray;

    if ((garray = (GRID_FLOAT_TYPE ***) malloc((size_t) pgrid->numx * sizeof (GRID_FLOAT_TYPE **))) == NULL)
        return ((void***) garray);
    NumAllocations++;

    numyz = pgrid->numy * pgrid->numz;
    for (ix = 0; ix < pgrid->numx; ix++) {
        if ((garray[ix] = (GRID_FLOAT_TYPE **) malloc((size_t) pgrid->numy * sizeof (GRID_FLOAT_TYPE *))) == NULL)
            return (NULL);
        NumAllocations++;
        for (iy = 0; iy < pgrid->numy; iy++) {
            garray[ix][iy] = (GRID_FLOAT_TYPE *) pgrid->buffer + ix * numyz + iy * pgrid->numz;
        }
    }

    pgrid->array = (void***) garray;

    return ((void***) garray);
}

/** function to free array for accessing 3D grid ***/

void DestroyGridArray(GridDesc * pgrid) {

    int ix;

    if (pgrid->array != NULL) {

        for (ix = 0; ix < pgrid->numx; ix++) {
            //printf("ix %d/%d pgrid->array[ix] %ld\n", ix, pgrid->numx, pgrid->array[ix]);
            free(pgrid->array[ix]);
            pgrid->array[ix] = NULL;
            NumAllocations--;
        }

        free(pgrid->array);
        pgrid->array = NULL;
        NumAllocations--;


    }

}

/** function to duplicate a 3D grid description and allocate grid memory */

void DuplicateGrid(GridDesc* pnew_grid, GridDesc* pold_grid, char *new_chr_type) {

    /* copy grid description */
    *pnew_grid = *pold_grid;

    /* set grid type */
    strcpy(pnew_grid->chr_type, new_chr_type);
    convert_grid_type(pnew_grid, 1);


    /* allocate grid */
    pnew_grid->buffer = AllocateGrid(pnew_grid);
    if (pnew_grid->buffer == NULL) {
        nll_puterr(
                "ERROR: allocating memory for duplicate 3D grid buffer.");
        exit(EXIT_ERROR_MEMORY);
    }

    /* create grid array access pointers */
    pnew_grid->array = CreateGridArray(pnew_grid);
    if (pnew_grid->array == NULL) {
        nll_puterr(
                "ERROR: creating array for accessing duplicate 3D grid buffer.");
        exit(EXIT_ERROR_MEMORY);
    }

}

/** function to perform several checks on grid data ***/

int CheckGridArray(GridDesc* pgrid, double gridMax, double gridMaxReplace,
        double gridMin, double gridMinReplace) {

    int ix, iy, iz;
    int ierror = 0, inegative = 0, imax = 0, imin = 0;
    GRID_FLOAT_TYPE val;

    for (ix = 0; ix < pgrid->numx; ix++) {
        for (iy = 0; iy < pgrid->numy; iy++) {
            for (iz = 0; iz < pgrid->numz; iz++) {

                /* check for negative values */
                if ((val = ((GRID_FLOAT_TYPE ***) pgrid->array)[ix][iy][iz]) < 0.0)
                    inegative++;

                /* check for out of range values */
                if (val > gridMax) {
                    val = gridMaxReplace;
                    imax++;
                } else if (val < gridMin) {
                    val = gridMinReplace;
                    imin++;
                }

            }
        }
    }

    if (inegative) {
        sprintf(MsgStr,
                "WARNING: %d negative values in grid.", inegative);
        nll_putmsg(1, MsgStr);
        ierror = -1;
    }
    if (imax) {
        sprintf(MsgStr,
                "WARNING: %d values > %e in grid replaced with %e",
                imax, gridMax, gridMaxReplace);
        nll_putmsg(1, MsgStr);
        ierror = -1;
    }
    if (imin) {

        sprintf(MsgStr,
                "WARNING: %d values < %e in grid replaced with %e",
                imin, gridMin, gridMinReplace);
        nll_putmsg(1, MsgStr);
        ierror = -1;
    }

    return (ierror);
}




/** function to sum 2 grids ***/

/* reads from  array if fp_grid_new == NULL */

//#define DEBUG_SUMGRIDS

int SumGrids(GridDesc* pgrid_sum, GridDesc* pgrid_new, FILE* fp_grid_new, double factor) {

    int ix, iy, iz;
    GRID_FLOAT_TYPE xval, yval, zval, newval;

    int identical = pgrid_new->array != NULL && testIdentical(pgrid_sum, pgrid_new);

    xval = pgrid_sum->origx;
    for (ix = 0; ix < pgrid_sum->numx; ix++) {

        yval = pgrid_sum->origy;
        for (iy = 0; iy < pgrid_sum->numy; iy++) {

            zval = pgrid_sum->origz;
            for (iz = 0; iz < pgrid_sum->numz; iz++) {

                if (identical) {
                    newval = ((GRID_FLOAT_TYPE ***) pgrid_new->array)[ix][iy][iz];
                } else {
                    newval = ReadAbsInterpGrid3d(fp_grid_new, pgrid_new, xval, yval, zval, 0);
                }
                if (newval > -LARGE_FLOAT) {
#ifdef DEBUG_SUMGRIDS
                    float result = ((GRID_FLOAT_TYPE ***) pgrid_sum->array)[ix][iy][iz] + factor * newval;
                    if (factor < 0.0 && fabs(result) > 1.0)
                        printf("DEBUG: SumGrids: %d %d %d  %f + %f * %f = %f\n",
                            ix, iy, iz, ((GRID_FLOAT_TYPE ***) pgrid_sum->array)[ix][iy][iz], factor, newval, result);
#endif
                    ((GRID_FLOAT_TYPE ***) pgrid_sum->array)[ix][iy][iz] += factor * newval;
                }

                zval += pgrid_sum->dz;
            }

            yval += pgrid_sum->dy;
        }

        xval += pgrid_sum->dx;
    }


    return (0);
}


/** function to multiply grid by a constant factor ***/

/* reads from  array if fp_grid_in == NULL */

//#define DEBUG_MULCONSTGRID

int MulConstGrid(GridDesc* pgrid_out, GridDesc* pgrid_in, FILE* fp_grid_in, double factor) {

    int ix, iy, iz;
    GRID_FLOAT_TYPE xval, yval, zval, newval;

    int identical = pgrid_in->array != NULL && testIdentical(pgrid_out, pgrid_in);

    xval = pgrid_out->origx;
    for (ix = 0; ix < pgrid_out->numx; ix++) {

        yval = pgrid_out->origy;
        for (iy = 0; iy < pgrid_out->numy; iy++) {

            zval = pgrid_out->origz;
            for (iz = 0; iz < pgrid_out->numz; iz++) {

                if (identical) {
                    newval = ((GRID_FLOAT_TYPE ***) pgrid_in->array)[ix][iy][iz];
                } else {
                    newval = ReadAbsInterpGrid3d(fp_grid_in, pgrid_in, xval, yval, zval, 0);
                }
                if (newval > -LARGE_FLOAT) {
#ifdef DEBUG_MULCONSTGRID
                    float result = ((GRID_FLOAT_TYPE ***) pgrid_out->array)[ix][iy][iz] + factor * newval;
                    if (factor < 0.0 && fabs(result) > 1.0)
                        printf("DEBUG: SumGrids: %d %d %d  %f + %f * %f = %f\n",
                            ix, iy, iz, ((GRID_FLOAT_TYPE ***) pgrid_out->array)[ix][iy][iz], factor, newval, result);
#endif
                    ((GRID_FLOAT_TYPE ***) pgrid_out->array)[ix][iy][iz] = factor * newval;
                }

                zval += pgrid_out->dz;
            }

            yval += pgrid_out->dy;
        }

        xval += pgrid_out->dx;
    }


    return (0);
}

/** check if two grids origin and dimensions are identical
 *
 * requires that grid are identical in size,
 * returns 1 if identical, returns 0 if not
 *
 */

int testIdentical(GridDesc* pGrid1, GridDesc* pGrid2) {

    //printf("DEBUG: testIdentical: test %s / %s\n", pGrid2->title, pGrid1->title);

    // check all relevant grid parameters are identical
    if (pGrid1->origx != pGrid2->origx
            || pGrid1->origy != pGrid2->origy
            || pGrid1->origz != pGrid2->origz
            ) {
        //printf("return 0 %f %f %f : %f %f %f\n", pGrid1->origx, pGrid1->origy, pGrid1->origz, pGrid2->origx, pGrid2->origy, pGrid2->origz);
        return (0);
    }
    if (pGrid1->dx != pGrid2->dx
            || pGrid1->dy != pGrid2->dy
            || pGrid1->dz != pGrid2->dz
            ) {
        //printf("return 0 %f %f %f : %f %f %f\n", pGrid1->dx, pGrid1->dy, pGrid1->dz, pGrid2->dx, pGrid2->dy, pGrid2->dz);
        return (0);
    }
    if (pGrid1->numx != pGrid2->numx
            || pGrid1->numy != pGrid2->numy
            || pGrid1->numz != pGrid2->numz
            ) {
        //printf("return 1 %d %d %d : %d %d %d\n", pGrid1->numx, pGrid1->numy, pGrid1->numz, pGrid2->numx, pGrid2->numy, pGrid2->numz);
        return (0);
    }

    return (1);

}

/** function to check if a point is on grid boundary */

/* return 10,11,20,21,30,31 if within tolerance of x,y,z orig/end boundary, return 0 otherwise */

int isOnGridBoundary(double xloc, double yloc, double zloc, GridDesc* pgrid,
        double tolerance_xy, double tolerance_z, int i_check_top) {

    // 20151118 AJL - bug fix, should not perform these tests for GLOBAL mode, if whole earth TODO: add check if global search volume is whole earth?
    // if (GeometryMode == MODE_GLOBAL) {
    if (GeometryMode != MODE_GLOBAL) {
        if (fabs(xloc - pgrid->origx) <= tolerance_xy)
            return (10);
        if (fabs(xloc - (pgrid->origx + (double) (pgrid->numx - 1) * pgrid->dx))
                <= tolerance_xy)
            return (11);
        if (fabs(yloc - pgrid->origy) <= tolerance_xy)
            return (20);
        if (fabs(yloc - (pgrid->origy + (double) (pgrid->numy - 1) * pgrid->dy))
                <= tolerance_xy)
            return (21);
    }

    if (i_check_top && fabs(zloc - pgrid->origz) <= tolerance_z)
        return (30);
    if (fabs(zloc - (pgrid->origz + (double) (pgrid->numz - 1) * pgrid->dz))
            <= tolerance_z)

        return (31);

    return (0);

}

/** function to check if a point is contained within a grid */

int IsPointInsideGrid(GridDesc* pgrid, double xloc, double yloc, double zloc) {

    if (xloc < pgrid->origx ||
            xloc > pgrid->origx + (double) (pgrid->numx - 1) * pgrid->dx)
        return (0);

    if (yloc < pgrid->origy ||
            yloc > pgrid->origy + (double) (pgrid->numy - 1) * pgrid->dy)
        return (0);

    if (zloc < pgrid->origz ||
            zloc > pgrid->origz + (double) (pgrid->numz - 1) * pgrid->dz)

        return (0);

    return (1);

}




/** function to check if a grid is entirely contained within another grid */

/* return 1 if pgrid_inside is inside pgrid, return 0 otherwise */

/* if iShiftFlag==1 shifts grid to attempt to get it inside */

int IsGridInside(GridDesc* pgrid_inside, GridDesc* pgrid, int iShiftFlag) {

    double xmin_in, xmax_in, ymin_in, ymax_in, zmin_in, zmax_in;
    double xmin, xmax, ymin, ymax, zmin, zmax;

    if (pgrid_inside == pgrid)
        return (1);

    xmin_in = pgrid_inside->origx;
    xmax_in = pgrid_inside->origx +
            (double) (pgrid_inside->numx - 1) * pgrid_inside->dx;
    ymin_in = pgrid_inside->origy;
    ymax_in = pgrid_inside->origy +
            (double) (pgrid_inside->numy - 1) * pgrid_inside->dy;
    zmin_in = pgrid_inside->origz;
    zmax_in = pgrid_inside->origz +
            (double) (pgrid_inside->numz - 1) * pgrid_inside->dz;

    xmin = pgrid->origx;
    xmax = pgrid->origx + (double) (pgrid->numx - 1) * pgrid->dx;
    ymin = pgrid->origy;
    ymax = pgrid->origy + (double) (pgrid->numy - 1) * pgrid->dy;
    zmin = pgrid->origz;
    zmax = pgrid->origz + (double) (pgrid->numz - 1) * pgrid->dz;

    if (!iShiftFlag) {
        if (xmin_in < xmin || xmax_in > xmax ||
                ymin_in < ymin || ymax_in > ymax ||
                zmin_in < zmin || zmax_in > zmax)
            return (0);
        else
            return (1);

    } else {
        if (xmin_in < xmin)
            pgrid_inside->origx += xmin - xmin_in;
        else if (xmax_in > xmax)
            pgrid_inside->origx -= xmax_in - xmax;
        if (ymin_in < ymin)
            pgrid_inside->origy += ymin - ymin_in;
        else if (ymax_in > ymax)
            pgrid_inside->origy -= ymax_in - ymax;
        if (zmin_in < zmin)
            pgrid_inside->origz += zmin - zmin_in;
        else if (zmax_in > zmax)
            pgrid_inside->origz -= zmax_in - zmax;

        return (IsGridInside(pgrid_inside, pgrid, 0));
    }

}


/** function to check if greatest distance from a station to a 3D grid
                is less than the horizontal extent of a 2D grid and
                if depth range of the 3D grid is less than that of the 2D grid
                and if station is within dist_horiz_min -> dist_horiz_max of grid location
                xcent, ycent */

/* return 1 if yes, return 0 otherwise */

int IsDistStaGridOK(GridDesc* pgrid_3D, SourceDesc* station,
        double dist_horiz_min, double dist_horiz_max, double xcent, double ycent) {

    double distance;


    /* check distances from grid location xcent, ycent */
    if (dist_horiz_min < dist_horiz_max) {
        // check if   dist_horiz_min < DIST < dist_horiz_max
        if (dist_horiz_min > SMALL_DOUBLE) {
            if ((distance = GetEpiDist(station, xcent, ycent)) < dist_horiz_min)
                return (-2);
        }
        if (dist_horiz_max > SMALL_DOUBLE) {
            if ((distance = GetEpiDist(station, xcent, ycent)) > dist_horiz_max)
                return (-2);
        }
    } else {
        // check if   DIST < dist_horiz_max && DIST > dist_horiz_min
        if (dist_horiz_min > SMALL_DOUBLE && dist_horiz_max > SMALL_DOUBLE) {
            distance = GetEpiDist(station, xcent, ycent);
            if (distance < dist_horiz_min && distance > dist_horiz_max)
                return (-2);
        }
    }

    return (1);

}


/** function to check if greatest distance from a station to a 3D grid
                is less than the horizontal extent of a 2D grid and
                if depth range of the 3D grid is less than that of the 2D grid
                and if station is within dist_horiz_min -> dist_horiz_max of grid location
                xcent, ycent */

/* return 1 if yes, return 0 otherwise */

int IsGrid2DBigEnough(GridDesc* pgrid_3D, GridDesc* pgrid_2D,
        SourceDesc* station,
        double dist_horiz_min, double dist_horiz_max, double xcent, double ycent) {

    double extent_horiz, distance;

    double xmin_in, xmax_in, ymin_in, ymax_in, zmin_in, zmax_in;
    double ymin, ymax, zmin, zmax;



    /* check distances from grid location xcent, ycent */
    if (dist_horiz_min < dist_horiz_max) {
        // check if   dist_horiz_min < DIST < dist_horiz_max
        if (dist_horiz_min > SMALL_DOUBLE) {
            if ((distance = GetEpiDist(station, xcent, ycent)) < dist_horiz_min)
                return (-2);
        }
        if (dist_horiz_max > SMALL_DOUBLE) {
            if ((distance = GetEpiDist(station, xcent, ycent)) > dist_horiz_max)
                return (-2);
        }
    } else {
        // check if   DIST < dist_horiz_max && DIST > dist_horiz_min
        if (dist_horiz_min > SMALL_DOUBLE && dist_horiz_max > SMALL_DOUBLE) {
            distance = GetEpiDist(station, xcent, ycent);
            if (distance < dist_horiz_min && distance > dist_horiz_max)
                return (-2);
        }
    }

    // 20160926 AJL - moved this block here from before previous block
    // GLOBAL
    if (GeometryMode == MODE_GLOBAL) {
        return (1); //INGV do not worry about grid size
    }


    /* get extent of 2D grid */
    ymin = pgrid_2D->origy;
    ymax = pgrid_2D->origy + (double) (pgrid_2D->numy - 1) * pgrid_2D->dy;
    extent_horiz = ymax - ymin;
    zmin = pgrid_2D->origz;
    zmax = pgrid_2D->origz + (double) (pgrid_2D->numz - 1) * pgrid_2D->dz;

    /* get bounds of 3D grid */
    xmin_in = pgrid_3D->origx;
    xmax_in = pgrid_3D->origx +
            (double) (pgrid_3D->numx - 1) * pgrid_3D->dx;
    ymin_in = pgrid_3D->origy;
    ymax_in = pgrid_3D->origy +
            (double) (pgrid_3D->numy - 1) * pgrid_3D->dy;
    zmin_in = pgrid_3D->origz;
    zmax_in = pgrid_3D->origz +
            (double) (pgrid_3D->numz - 1) * pgrid_3D->dz;

    // adjust for GLOBAL (grid units = deg, GetEpiDist units = km)
    if (GeometryMode == MODE_GLOBAL)
        extent_horiz /= KM2DEG;

    /* check each horizontal corner */
    if ((distance = GetEpiDist(station, xmin_in, ymin_in)) > extent_horiz)
        return (-1);
    if ((distance = GetEpiDist(station, xmin_in, ymax_in)) > extent_horiz)
        return (-1);
    if ((distance = GetEpiDist(station, xmax_in, ymax_in)) > extent_horiz)
        return (-1);
    if ((distance = GetEpiDist(station, xmax_in, ymin_in)) > extent_horiz)
        return (-1);
    /* check depth range */
    if (zmin_in < zmin || zmax_in > zmax)

        return (-3);


    return (1);
}

/** function to calculate epicentral (horizontal) distance from a source to an x-y position */

double GetEpiDist(SourceDesc* psrce, double xval, double yval) {
    double xtmp, ytmp;

    if (GeometryMode == MODE_GLOBAL) {
        return (GCDistance(yval, xval, psrce->y, psrce->x));
        //return(EllipsoidDistance(yval, xval, psrce->y, psrce->x));
    } else {
        xtmp = xval - psrce->x;
        ytmp = yval - psrce->y;

        return (sqrt(xtmp * xtmp + ytmp * ytmp));
    }
}

/** function to calculate epicentral (horizontal) distance from a station to an x-y position */

double GetEpiDistSta(StationDesc* psta, double xval, double yval) {
    double xtmp, ytmp;

    if (GeometryMode == MODE_GLOBAL) {
        return (GCDistance(yval, xval, psta->y, psta->x));
        //return(EllipsoidDistance(yval, xval, psta->y, psta->x));
    } else {
        xtmp = xval - psta->x;
        ytmp = yval - psta->y;

        return (sqrt(xtmp * xtmp + ytmp * ytmp));
    }
}

/** function to calculate epicentral (horizontal) azimuth from a x-y position to a source */

double GetEpiAzim(SourceDesc* psrce, double xval, double yval) {
    double xtmp, ytmp, azim;

    if (GeometryMode == MODE_GLOBAL) {
        return (GCAzimuth(yval, xval, psrce->y, psrce->x));
    } else {
        xtmp = psrce->x - xval;
        ytmp = psrce->y - yval;
        azim = atan2(xtmp, ytmp) / cRPD;
        if (azim < 0.0)
            azim += 360.0;

        return (azim);
    }
}

/** function to calculate epicentral (horizontal) azimuth from a
        x-y position to a station */

double GetEpiAzimSta(StationDesc* psta, double xval, double yval) {
    double xtmp, ytmp, azim;

    if (GeometryMode == MODE_GLOBAL) {
        return (GCAzimuth(yval, xval, psta->y, psta->x));
    } else {
        xtmp = psta->x - xval;
        ytmp = psta->y - yval;
        azim = atan2(xtmp, ytmp) / cRPD;
        if (azim < 0.0)
            azim += 360.0;

        return (azim);
    }
}

/** function to calculate distance between 2 XYZ points */

double Dist2D(double x1, double x2, double y1, double y2) {

    if (GeometryMode == MODE_GLOBAL) {
        return (GCDistance(y1, x1, y2, x2));
        //return(EllipsoidDistance(yval, xval, psta->y, psta->x));
    } else {
        double dx, dy;
        dx = x1 - x2;
        dy = y1 - y2;
        return (sqrt(dx * dx + dy * dy));
    }
}

/** function to calculate distance between 2 XYZ points */

double Dist3D(double x1, double x2, double y1, double y2,
        double z1, double z2) {

    double epi_dist = Dist2D(x1, x2, y1, y2);
    double dz = z1 - z2;

    return (sqrt(epi_dist * epi_dist + dz * dz));
}

/** function to calculate average inter-station distance */

double calcAveInterStationDistance(SourceDesc *stations, int numStations) {

    int ndist, n, m = 0;
    double x, y;
    double x2, y2;
    double dist_sum, dist_ave;


    // calculate average distance
    dist_sum = 0.0;
    ndist = 0;
    for (n = 0; n < numStations; n++) {
        // 20101021 AJL = Bug fix
        // check if station has ignored reading for this event
        if ((stations + n)->ignored)
            continue;
        x = (stations + n)->x;
        y = (stations + n)->y;
        if (!stationLocationIsKnown(x, y)) // station location not known
            continue;
        for (m = 0; m < n; m++) {
            if ((stations + m)->ignored)
                continue;
            x2 = (stations + m)->x;
            y2 = (stations + m)->y;
            if (!stationLocationIsKnown(x2, y2)) // station location not known
                continue;
            dist_sum += GetEpiDist((stations + m), x, y);
            ndist++;
        }
        //printf("Sta: %s  x %f y %f  dist_sum %f\n", (stations + n)->label, x, y, dist_sum);
    }
    if (ndist <= 0)
        return (-1.0);
    dist_ave = dist_sum / (double) ndist;

    return (dist_ave);

}

/** function to check if station location is known */

int stationLocationIsKnown(double x, double y) {

    // station location not known
    if (x == 0.0 && y == 0.0)
        return (0);

    // station location not known
    if (x < -LARGE_DOUBLE / 2.0 || y < -LARGE_DOUBLE / 2.0)

        return (0);

    // station location known
    return (1);

}

/** function to generate a grid from an OctTree structure */

int ConvertOctTree2Grid(Tree3D* tree, double dx, double dy, double dz, char *grid_type, GridDesc * pgrid_out) {

    int ix, iy, iz;
    Vect3D coords;
    OctNode* node;
    GRID_FLOAT_TYPE value;

    pgrid_out->numx = 1 + (int) ((double) tree->numx * tree->ds.x / dx);
    pgrid_out->numy = 1 + (int) ((double) tree->numy * tree->ds.y / dy);
    pgrid_out->numz = 1 + (int) ((double) tree->numz * tree->ds.z / dz);
    pgrid_out->origx = tree->orig.x;
    pgrid_out->origy = tree->orig.y;
    pgrid_out->origz = tree->orig.z;
    pgrid_out->dx = dx;
    pgrid_out->dy = dy;
    pgrid_out->dz = dz;
    if (grid_type != NULL) {
        strcpy(pgrid_out->chr_type, grid_type);
        convert_grid_type(pgrid_out, 1);
    } else {
        pgrid_out->type = tree->data_code;
        convert_grid_type(pgrid_out, 0);
    }

    /* allocate grid buffer */
    if (message_flag >= 4)
        display_grid_param(pgrid_out);
    pgrid_out->buffer = AllocateGrid(pgrid_out);
    if (pgrid_out->buffer == NULL) {
        nll_puterr("ERROR: allocating memory for 3D PDF grid buffer.");
        exit(EXIT_ERROR_MEMORY);
    }

    /* create grid array access pointers */

    pgrid_out->array = CreateGridArray(pgrid_out);
    if (pgrid_out->array == NULL) {
        nll_puterr("ERROR: creating array for accessing 3D PDF grid buffer.");
        exit(EXIT_ERROR_MEMORY);
    }


    /*
    // use values for center of target cell
    coords.x = pgrid_out->origx + dx / 2.0;	// center of cell
    for (ix = 0; ix < pgrid_out->numx; ix++) {
            if (ix == pgrid_out->numx - 1)
                    coords.x -= dx / 2.0 + dx / 1000.0;
            coords.y = pgrid_out->origy + dy / 2.0;
            for (iy = 0; iy < pgrid_out->numy; iy++) {
                    if (iy == pgrid_out->numy - 1)
                            coords.y -= dy / 2.0 + dy / 1000.0;
                    coords.z = pgrid_out->origz + dz / 2.0;
                    for (iz = 0; iz < pgrid_out->numz; iz++) {
                            if (iz == pgrid_out->numz - 1)
                                    coords.z -= dz / 2.0 + dz / 1000.0;
                            node = getLeafNodeContaining(tree, coords);
                            value = (GRID_FLOAT_TYPE) node->value;
                            pgrid_out->array[ix][iy][iz] = value;

                            coords.z += dz;
                    }
                    coords.y += dy;
            }
            coords.x += dx;
    }
     */
    // use values for corner of target cell
    coords.x = pgrid_out->origx; // corner of cell
    for (ix = 0; ix < pgrid_out->numx; ix++) {
        if (ix == pgrid_out->numx - 1)
            coords.x -= dx / 1000.0;
        coords.y = pgrid_out->origy;
        for (iy = 0; iy < pgrid_out->numy; iy++) {
            if (iy == pgrid_out->numy - 1)
                coords.y -= dy / 1000.0;
            coords.z = pgrid_out->origz;
            for (iz = 0; iz < pgrid_out->numz; iz++) {

                if (iz == pgrid_out->numz - 1)
                    coords.z -= dz / 1000.0;
                node = getLeafNodeContaining(tree, coords);
                value = (GRID_FLOAT_TYPE) node->value;
                ((GRID_FLOAT_TYPE ***) pgrid_out->array)[ix][iy][iz] = value;

                coords.z += dz;
            }
            coords.y += dy;
        }
        coords.x += dx;
    }

    return (0);

}

/** function to swap bytes in grid buffer ***/

int swapBytes(GRID_FLOAT_TYPE *buffer, long bufsize) {

    GRID_FLOAT_TYPE *bufpos;
    char ctmp;
    unsigned short itmp;

    union {
        char cval[4];
        unsigned short ival[2];
        float fval;
    }
    bufval;


    bufpos = buffer;
    while (bufpos < buffer + bufsize) {

        bufval.fval = *bufpos;
        ctmp = bufval.cval[0];
        bufval.cval[0] = bufval.cval[1];
        bufval.cval[1] = ctmp;
        ctmp = bufval.cval[2];
        bufval.cval[2] = bufval.cval[3];
        bufval.cval[3] = ctmp;
        itmp = bufval.ival[0];
        bufval.ival[0] = bufval.ival[1];
        bufval.ival[1] = itmp;
        *bufpos = bufval.fval;
        bufpos++;
    }

    return (0);
}

/** function to read entire grid buffer from disk ***/

int ReadGrid3dBuf(GridDesc* pgrid, FILE * fpio) {

    long readsize;



    // 20161021 AJL  readsize = pgrid->numx * pgrid->numy * pgrid->numz * sizeof (GRID_FLOAT_TYPE);
    readsize = pgrid->buffer_size;

    /* read from grid file to buffer */

    int ireturn;
    if ((ireturn = fread((char *) pgrid->buffer, readsize, 1, fpio)) != 1) {
        printf("DEBUG: pgrid->buffer %ld, readsize %ld, fpio %ld, ireturn %d\n", (long) pgrid->buffer, readsize, (long) fpio, ireturn);
        nll_puterr2("ERROR: reading grid file", pgrid->title);
        return (-1);
    }

    if (pgrid->iSwapBytes)
        swapBytes(pgrid->buffer, readsize / sizeof (float));

    return (0);
}

/** function to read y-z sheet of grid buffer from disk ***/

int ReadGrid3dBufSheet(GRID_FLOAT_TYPE* sheetbuf, GridDesc* pgrid_disk,
        FILE* fpio, int ix) {

    long offset;
    long readsize;


    /* check indexes in range */

    if (ix < 0 || ix >= pgrid_disk->numx) {
        sprintf(MsgStr,
                "WARNING: grid file x-sheet index %d out of range (%d,%d)",
                ix, 0, pgrid_disk->numx - 1);
        return (-1);
    }


    /* calculate offset in bytes */

    offset = sizeof (GRID_FLOAT_TYPE) * (ix * (pgrid_disk->numy * pgrid_disk->numz));
    fseek(fpio, offset, SEEK_SET);


    /* calculate number of bytes to read */

    readsize = pgrid_disk->numy * pgrid_disk->numz * sizeof (GRID_FLOAT_TYPE);


    /* read sheet from grid file to buffer */

    if (fread((char *) sheetbuf, readsize, 1, fpio) != 1) {
        nll_puterr("ERROR: reading x-sheet grid file.");
        return (-1);
    }

    if (pgrid_disk->iSwapBytes)
        swapBytes(sheetbuf, readsize / sizeof (GRID_FLOAT_TYPE));

    return (0);
}

/** function to read grid header file ***/

int ReadGrid3dHdr(GridDesc* pgrid, SourceDesc* psrce, char* filename, char* file_type) {

    FILE *fpio;
    char fname[FILENAME_MAX];



    /* read header file */

    sprintf(fname, "%s.%s.hdr", filename, file_type);
    if ((fpio = fopen(fname, "r")) == NULL) {
        if (message_flag >= 1)
            nll_puterr2("ERROR: opening grid header file: %s", fname);
        return (-1);
    }
    NumFilesOpen++;


    if (ReadGrid3dHdr_grid_description(fpio, pgrid, fname) < 0) {
        fclose(fpio);
        NumFilesOpen--;
        return (-1);
    }

    if (strncmp(file_type, "time", 4) == 0
            || strncmp(file_type, "angle", 4) == 0)
        fscanf(fpio, "%s %lf %lf %lf\n",
            psrce->label, &(psrce->x), &(psrce->y), &(psrce->z));


    //  check for optional header lines
    char line[MAXLINE_LONG];
    char tag[MAXLINE_LONG];

    // check for map projection
    strcpy(pgrid->mapProjStr, "");
    rewind(fpio);
    while (fgets(line, MAXLINE_LONG, fpio) != NULL) {
        int istat = sscanf(line, "%s", tag);
        if (istat == 1 && strcmp(tag, "TRANSFORM") == 0) {
            strcpy(pgrid->mapProjStr, line);
        }
    }

    // check if cascading grid
    pgrid->flagGridCascading = IS_NOT_CASCADING;
    int num_z_merge_depths;
    rewind(fpio);
    while (fgets(line, MAXLINE_LONG, fpio) != NULL) {
        int istat = sscanf(line, "%s %d", tag, &num_z_merge_depths);
        if (istat == 2 && strcmp(tag, "CASCADING_GRID") == 0) {
            setCascadingGrid(pgrid);
            pgrid->gridDesc_Cascading.num_z_merge_depths = num_z_merge_depths;
            if (pgrid->gridDesc_Cascading.num_z_merge_depths > MAX_NUM_Z_MERGE_DEPTHS) {
                pgrid->gridDesc_Cascading.num_z_merge_depths = MAX_NUM_Z_MERGE_DEPTHS;
                sprintf(MsgStr, "ERROR: too many cascading grid Z merge depths, only using first %d depths.",
                        pgrid->gridDesc_Cascading.num_z_merge_depths);
                nll_puterr(MsgStr);
            }
            // 20170207 AJL - z_merge_depths moved to fixed array to ease memory management
            //pgrid->gridDesc_Cascading.z_merge_depths = (double*) malloc((size_t) num_z_merge_depths * sizeof (double));
            char doubling_depths[1024];
            sscanf(line, "%*s %*d %s", doubling_depths);
            char *str_pos = strtok(doubling_depths, ",");
            int n = 0;
            while (str_pos != NULL) {
                pgrid->gridDesc_Cascading.z_merge_depths[n] = atof(str_pos);
                //printf("DEBUG: CASCADING_GRID doubling depth added: %s %f\n", str_pos, pgrid->gridDesc_Cascading.z_merge_depths[n]);
                n++;
                str_pos = strtok(NULL, ",");
            }
        }
    }

    fclose(fpio);
    NumFilesOpen--;

    return (0);
}

/** function to read grid header file grid description line ***/

int ReadGrid3dHdr_grid_description(FILE *fpio, GridDesc *pgrid, char *fname) {

    char line[MAXLINE_LONG];
    if (fgets(line, MAXLINE_LONG, fpio) == NULL) {
        nll_puterr2("ERROR: reading grid header file: ", fname);
        return (-1);
    }
    strcpy(pgrid->float_type, "FLOAT");
    sscanf(line, "%d %d %d  %lf %lf %lf  %lf %lf %lf %s %s",
            &(pgrid->numx), &(pgrid->numy), &(pgrid->numz),
            &(pgrid->origx), &(pgrid->origy), &(pgrid->origz),
            &(pgrid->dx), &(pgrid->dy), &(pgrid->dz), pgrid->chr_type, pgrid->float_type);

#ifdef GRID_FLOAT_TYPE_DOUBLE
    if (strcmp(pgrid->float_type, "DOUBLE") != 0) {
        nll_puterr("ERROR: Global grid float type is set to DOUBLE, but grid file type is not DOUBLE. (see compiler flag GRID_FLOAT_TYPE_DOUBLE)");
        return (-1);
    }
#else
    if (strcmp(pgrid->float_type, "FLOAT") != 0) {
        nll_puterr("ERROR: Global grid float type is set to FLOAT, but grid file type is not FLOAT. (see compiler flag GRID_FLOAT_TYPE_DOUBLE)");

        return (-1);
    }
#endif

    return (0);

}

/** function to open grid file and read header ***/

int OpenGrid3dFile(char *fname, FILE **fp_grid, FILE **fp_hdr,
        GridDesc* pgrid, char* file_type, SourceDesc* psrce, int iSwapBytes) {

    char fn_grid[FILENAME_MAX], fn_hdr[FILENAME_MAX];

    /* open grid file and header file */

    sprintf(fn_grid, "%s.buf", fname);
    if (message_flag >= 3) {
        sprintf(MsgStr, "Opening Grid File: %s", fn_grid);
        nll_putmsg(3, MsgStr);
    }
    if ((*fp_grid = fopen(fn_grid, "r")) == NULL) {
        if (message_flag >= 3) {
            sprintf(MsgStr, "WARNING: cannot open grid buffer file: %s", fn_grid);
            nll_putmsg(3, MsgStr);
        }
        //return(-1);	// sometimes only header is wanted
    } else {
        NumGridBufFilesOpen++;
        NumFilesOpen++;
    }
    sprintf(fn_hdr, "%s.hdr", fname);
    if ((*fp_hdr = fopen(fn_hdr, "r")) == NULL) {
        if (message_flag >= 3) {
            sprintf(MsgStr,
                    "WARNING: cannot open grid header file: %s", fn_hdr);
            nll_putmsg(3, MsgStr);
        }
        if (*fp_grid != NULL) {
            fclose(*fp_grid);
            NumGridBufFilesOpen--;
            NumFilesOpen--;
        }
        return (-1);
    }
    NumGridHdrFilesOpen++;
    NumFilesOpen++;

    // initialize key fields
    pgrid->array = NULL;
    pgrid->buffer = NULL;


    /* read header file */

    pgrid->iSwapBytes = iSwapBytes;
    if (ReadGrid3dHdr_grid_description(*fp_hdr, pgrid, fn_hdr) < 0) {
        fclose(*fp_hdr);
        NumGridBufFilesOpen--;
        NumFilesOpen--;
        if (*fp_grid != NULL) {
            fclose(*fp_grid);
            NumGridBufFilesOpen--;
            NumFilesOpen--;
        }
        return (-1);
    }

    // make sure that dx for 2D grids is non-zero
    if (pgrid->numx == 1)
        pgrid->dx = 1.0;


    convert_grid_type(pgrid, 1);
    if (message_flag >= 4)
        display_grid_param(pgrid);

    if (psrce != NULL && (strncmp(file_type, "time", 4) == 0
            || strncmp(file_type, "angle", 4) == 0))
        fscanf(*fp_hdr, "%s %lf %lf %lf\n",
            psrce->label, &(psrce->x), &(psrce->y), &(psrce->z));

    // save filename as grid identifier
    strcpy(pgrid->title, fname);

    //  check for optional header lines
    char line[MAXLINE_LONG];
    char tag[MAXLINE_LONG];

    // check for map projection
    strcpy(pgrid->mapProjStr, "");
    rewind(*fp_hdr);
    while (fgets(line, MAXLINE_LONG, *fp_hdr) != NULL) {
        int istat = sscanf(line, "%s", tag);
        if (istat == 1 && strcmp(tag, "TRANSFORM") == 0) {
            strcpy(pgrid->mapProjStr, line);
        }
    }

    // check if cascading grid
    pgrid->flagGridCascading = IS_NOT_CASCADING;
    int num_z_merge_depths;
    rewind(*fp_hdr);
    while (fgets(line, MAXLINE_LONG, *fp_hdr) != NULL) {
        int istat = sscanf(line, "%s %d", tag, &num_z_merge_depths);
        if (istat == 2 && strcmp(tag, "CASCADING_GRID") == 0) {
            setCascadingGrid(pgrid);
            pgrid->gridDesc_Cascading.num_z_merge_depths = num_z_merge_depths;
            if (pgrid->gridDesc_Cascading.num_z_merge_depths > MAX_NUM_Z_MERGE_DEPTHS) {
                pgrid->gridDesc_Cascading.num_z_merge_depths = MAX_NUM_Z_MERGE_DEPTHS;
                sprintf(MsgStr, "ERROR: too many cascading grid Z merge depths, only using first %d depths.",
                        pgrid->gridDesc_Cascading.num_z_merge_depths);
                nll_puterr(MsgStr);
            }
            // 20170207 AJL - z_merge_depths moved to fixed array to ease memory management
            //pgrid->gridDesc_Cascading.z_merge_depths = (double*) malloc((size_t) num_z_merge_depths * sizeof (double));
            char doubling_depths[1024];
            sscanf(line, "%*s %*d %s", doubling_depths);
            char *str_pos = strtok(doubling_depths, ",");
            int n = 0;
            while (str_pos != NULL) {
                pgrid->gridDesc_Cascading.z_merge_depths[n] = atof(str_pos);
                //printf("DEBUG: CASCADING_GRID doubling depth added: %s %f\n", str_pos, pgrid->gridDesc_Cascading.z_merge_depths[n]);
                n++;
                str_pos = strtok(NULL, ",");
            }
        }
    }



    return (0);

}

/** function to close grid file and header ***/

void CloseGrid3dFile(GridDesc* pgrid, FILE **fp_grid, FILE **fp_hdr) {

    // 20170207 AJL - z_merge_depths moved to fixed array to ease memory management
    /*if (pgrid != NULL) {
        if (isCascadingGrid(pgrid)) {
            if (pgrid->gridDesc_Cascading.z_merge_depths != NULL) {
                free(pgrid->gridDesc_Cascading.z_merge_depths);
                pgrid->gridDesc_Cascading.z_merge_depths = NULL;
                NumAllocations--;
            }
            pgrid->gridDesc_Cascading.z_merge_depths = NULL;
        }
    }*/

    if (*fp_grid != NULL) {
        fclose(*fp_grid);
        *fp_grid = NULL;
        NumGridBufFilesOpen--;
        NumFilesOpen--;
    }

    if (*fp_hdr != NULL) {
        fclose(*fp_hdr);
        *fp_hdr = NULL;
        NumGridHdrFilesOpen--;
        NumFilesOpen--;
    }

}

/** function to read values from 2D or 3D grid file
 *
 *  The required grids must be present as *.buf and *.hyp disk files.
 */

GRID_FLOAT_TYPE * ReadGridFile
(
        // calling parameters
        GRID_FLOAT_TYPE* values, // GRID_FLOAT_TYPE array to return grid values, must be of size nvalues or larger
        char *fname, // grid file path and root name
        char* file_type, // grid file type (e.g. "mod", "time" or "angle")
        double* xloc, double* yloc, double* zloc, // arrays of x, y, z coordinates in grid units of target point, arrays must be of size nvalues or larger
        int nvalues, // number of values to read
        int iSwapBytes, // flag to indicate if hi and low bytes of input velocity grid file should be swapped (1=swap, 0=no swap)
        SourceDesc* psrceIn // station source to use instead of srce in grid file (use only for station=DEFAULT 2D grids, use NULL otherwise)
        ) {
    int istat, i;
    FILE *fp_grid, *fp_hdr;
    GridDesc gdesc;
    SourceDesc srce;
    SourceDesc* psrce_use;
    double yval_grid;


    // initialize values to invalid
    for (i = 0; i < nvalues; i++)
        values[i] = -VERY_LARGE_FLOAT;


    // open grid file
    if ((istat = OpenGrid3dFile(fname, &fp_grid, &fp_hdr, &gdesc, file_type, &srce, iSwapBytes)) < 0) {
        if (message_flag >= 3) {
            sprintf(MsgStr, "WARNING: cannot open grid file: %s", fname);
            nll_putmsg(3, MsgStr);
        }
        return (values);
    }

    if (gdesc.type == GRID_TIME_2D || gdesc.type == GRID_ANGLE_2D) {
        // 2D grid (1D model)
        psrce_use = &srce;
        if (psrceIn != NULL) {
            psrce_use = psrceIn;
        }
        //        printf("DEBUG: SOURCE: Name: %s  Loc:  latlon: %d  X(east) %lg  Y(north) %lg  Z(pos DOWN) %lg  lat %lg lon %lg\n",
        //                psrce_use->label, psrce_use->is_coord_latlon,
        //                psrce_use->x, psrce_use->y, psrce_use->z,
        //                psrce_use->dlat, psrce_use->dlong
        //                );
        for (i = 0; i < nvalues; i++) {
            yval_grid = GetEpiDist(psrce_use, xloc[i], yloc[i]);
            if (GeometryMode == MODE_GLOBAL)
                yval_grid *= KM2DEG;
            values[i] = ReadAbsInterpGrid2d(fp_grid, &gdesc, yval_grid, zloc[i]);
            //            printf("DEBUG: yval_grid=%f xloc[i]=%f yloc[i]=%f zloc[i]=%f values[i]=%f\n", yval_grid, xloc[i], yloc[i], zloc[i], values[i]);
        }
    } else {
        // 3D grid
        for (i = 0; i < nvalues; i++) {
            // get GRID_FLOAT_TYPE value on grid
            values[i] = ReadAbsInterpGrid3d(fp_grid, &gdesc, xloc[i], yloc[i], zloc[i], 0);
        }
    }


    // close grid file
    CloseGrid3dFile(&gdesc, &fp_grid, &fp_hdr);

    return (values);

}

/** function to read cascading grid data from disk or array at index location
 *
 * 20161019 AJL - added
 *
 *  ix_casc, iy_casc, iz_casc are indices in cascading grid index units
 */

GRID_FLOAT_TYPE ReadCascadingGrid3dValue(FILE *fpgrid, int ix_casc, int iy_casc, int iz_casc, GridDesc * pgrid) {

    long offset;
    GRID_FLOAT_TYPE fvalue;

    // get fvalue
    if (fpgrid != NULL) {
        // calculate offset in bytes from array offset
        offset = ((GRID_FLOAT_TYPE **) pgrid->array[ix_casc][iy_casc] + (long) iz_casc) - (GRID_FLOAT_TYPE **) pgrid->array[0][0];
        offset *= sizeof (GRID_FLOAT_TYPE);
        //offset *= 2;
        //if (ix_casc < 4 && iy_casc < 4 && iz_casc < 4)
        //    printf("ixyz %d,%d,%d  ixyz_casc %d,%d,%d  offset %ld/%ld = (GRID_FLOAT_TYPE **) pgrid->array[ix_casc][iy_casc] %ld + iz_casc %d - (GRID_FLOAT_TYPE **) pgrid->array[0][0] %ld\n",
        //       ix, iy, iz, ix_casc, iy_casc, iz_casc, offset, offset * sizeof (GRID_FLOAT_TYPE), (GRID_FLOAT_TYPE **) pgrid->array[ix_casc][iy_casc], iz_casc, (GRID_FLOAT_TYPE **) pgrid->array[0][0]);
        fseek(fpgrid, offset, SEEK_SET);
        // read fvalue
        if (fread(&fvalue, sizeof (GRID_FLOAT_TYPE), 1, fpgrid) != 1) {
            sprintf(MsgStr, "ERROR: reading cascading grid value at: x%d y%d z%d (offset %ld buf_size %ld diff %ld) in: %s",
                    ix_casc, iy_casc, iz_casc, offset, pgrid->buffer_size, offset - pgrid->buffer_size, pgrid->title);
            nll_puterr(MsgStr);
            return (-VERY_LARGE_FLOAT);
        }
        if (pgrid->iSwapBytes)
            swapBytes(&fvalue, 1);
    } else {
        fvalue = ((GRID_FLOAT_TYPE ***) pgrid->array)[ix_casc][iy_casc][iz_casc];
    }

    return (fvalue);
}

/** function to read cascading grid data from disk or array at interpolated point corresponding to regular grid index location
 *
 * 20161021 AJL - added
 *
 *  ix_dbl, iy_dbl, iz_dbl are indices in virtual, regular grid equivalent to this cascading grid
 */

GRID_FLOAT_TYPE ReadGrid3dValue_Cascading_Interp(FILE *fpgrid, double ix_dbl, double iy_dbl, double iz_dbl, GridDesc * pgrid, int clean_casc_allocs) {

    int ix = (int) ix_dbl;
    int iy = (int) iy_dbl;
    int iz = (int) iz_dbl;

    // check indexes in range
    // ix, iy, iz are indices in virtual, regular grid equivalent to this cascading grid
    if (ix < 0 || ix >= pgrid->numx || iy < 0 || iy >= pgrid->numy || iz < 0 || iz >= pgrid->numz) {
        //nll_puterr("WARNING: grid file index out of range.");
        return (-VERY_LARGE_FLOAT);
    }

    // need grid array to find offset
    int created_grid_array_cascading = 0;
    int allocated_grid_cascading = 0;
    if (pgrid->array == NULL) {
        if (pgrid->buffer == NULL) {
            // prepare grid memory without allocating buffer, required to initialize cascading grid array
            AllocateGrid_Cascading(pgrid, 0);
            allocated_grid_cascading = 1;
        }
        pgrid->array = CreateGridArray_Cascading(pgrid);
        created_grid_array_cascading = 1;
    }

    int xyz_scale = pgrid->gridDesc_Cascading.xyz_scale[iz];

    // set upper x,y,z indices in cascading grid
    int iz0_casc = pgrid->gridDesc_Cascading.zindex[iz];
    int ix0_casc = ix / xyz_scale;
    int iy0_casc = iy / xyz_scale;


    // interpolate within cascading grid cell

    // check if at change in grid cell size, gets complicated...
    int xy_scale_use = xyz_scale;
    int ix0_casc_up, iy0_casc_up, ix0_casc_dn, iy0_casc_dn;
    int ix1_casc_up, iy1_casc_up, ix1_casc_dn, iy1_casc_dn;
    int iz1_casc;
    int irescale = 0;
    if (iz < pgrid->numz - 2) {
        // find first regular z index in next cascading grid level
        int iz_test = iz + 1;
        // 20170207 AJL while (iz_test < pgrid->numz && pgrid->gridDesc_Cascading.zindex[iz_test] == iz0_casc) {
        while (iz_test < pgrid->numz - 1 && pgrid->gridDesc_Cascading.zindex[iz_test] == iz0_casc) {
            iz_test++;
        }
        // check if scale changes in next cascading grid level
        irescale = pgrid->gridDesc_Cascading.xyz_scale[iz_test] > xyz_scale;
    }
    if (irescale) { // change of grid size
        //printf("DEBUG: rescale! iz %d scale %d, iz+1 %d, scale %d\n",iz, xyz_scale, iz + 1, pgrid->gridDesc_Cascading.xyz_scale[iz + 1]);
        xy_scale_use = 2 * xyz_scale;
        // set upper and lower x,y indices in cascading grid
        // x
        // make sure upper values are at even index, align with cell below
        ix0_casc_up = 2 * (ix0_casc / 2);
        // upper index increases by 2
        ix1_casc_up = ix0_casc_up + 2;
        int ixmax = (pgrid->numx - 1) / xyz_scale + ((pgrid->numx - 1) % xyz_scale == 0 ? 0 : 1); // includes fractional final casc grid node if xyz_scale != 1
        if (ix1_casc_up > ixmax) {
            ix1_casc_up = ixmax;
        }
        // lower index in larger cell size increases by 1
        ix0_casc_dn = ix0_casc_up / 2;
        ix1_casc_dn = ix0_casc_dn + 1;
        ixmax = (pgrid->numx - 1) / xy_scale_use + ((pgrid->numx - 1) % xy_scale_use == 0 ? 0 : 1); // includes fractional final casc grid node if xyz_scale != 1
        if (ix1_casc_dn > ixmax) {
            ix1_casc_dn = ixmax;
        }
        // y
        // make sure upper values are at even index, align with cell below
        iy0_casc_up = 2 * (iy0_casc / 2);
        // upper index increases by 2
        iy1_casc_up = iy0_casc_up + 2;
        int iymax = (pgrid->numy - 1) / xyz_scale + ((pgrid->numy - 1) % xyz_scale == 0 ? 0 : 1); // includes fractional final casc grid node if xyz_scale != 1
        if (iy1_casc_up > iymax) {
            iy1_casc_up = iymax;
        }
        // lower index in larger cell size increases by 1
        iy0_casc_dn = iy0_casc_up / 2;
        iy1_casc_dn = iy0_casc_dn + 1;
        iymax = (pgrid->numy - 1) / xy_scale_use + ((pgrid->numy - 1) % xy_scale_use == 0 ? 0 : 1); // includes fractional final casc grid node if xyz_scale != 1
        if (iy1_casc_dn > iymax) {
            iy1_casc_dn = iymax;
        }
    } else {
        // set upper and lower x,y indices in cascading grid
        // x
        ix0_casc_up = ix0_casc_dn = ix0_casc;
        ix1_casc_up = ix0_casc + 1;
        if (ix1_casc_up >= pgrid->numx) {
            ix1_casc_up = pgrid->numx - 1;
        }
        ix1_casc_dn = ix1_casc_up;
        // y
        iy0_casc_up = iy0_casc_dn = iy0_casc;
        iy1_casc_up = iy0_casc + 1;
        if (iy1_casc_up >= pgrid->numy) {
            iy1_casc_up = pgrid->numy - 1;
        }
        iy1_casc_dn = iy1_casc_up;
    }
    iz1_casc = iz0_casc + 1;
    if (iz1_casc > pgrid->gridDesc_Cascading.zindex[pgrid->numz - 1]) {
        iz1_casc = pgrid->gridDesc_Cascading.zindex[pgrid->numz - 1];
    }


    DOUBLE xdiff, ydiff, zdiff;
    // x
    int lastx = ((pgrid->numx - 1) / xy_scale_use) * xy_scale_use; // reg index of end of last full casc grid cell
    if (ix > lastx) { // casc grid cell is truncated in x relative to cell of size xy_scale_use reg grid cells
        xdiff = (ix_dbl - (DOUBLE) lastx) / (DOUBLE) (pgrid->numx - 1 - lastx);
        if (ix0_casc_up != ix0_casc_dn && iy == 133) {
            printf("xy_scale_use %d, xdiff %f = (DOUBLE) (ix %d - lastx %d) / (DOUBLE) (pgrid->numx %d - 1 - lastx %d)\n", xy_scale_use, xdiff, ix, lastx, pgrid->numx, lastx);
            printf("ix0_casc_up/dn %d/%d, ix1_casc_up/dn %d/%d, iz0_casc %d, iz1_casc %d\n", ix0_casc_up, ix0_casc_dn, ix1_casc_up, ix1_casc_dn, iz0_casc, iz1_casc);
            printf("iy0_casc_up %d, iy0_casc_dn %d, iy1_casc_up %d, iy1_casc_dn %d, iz0_casc %d, iz1_casc %d\n", iy0_casc_up, iy0_casc_dn, iy1_casc_up, iy1_casc_dn, iz0_casc, iz1_casc);
        }
    } else {
        //xdiff = (DOUBLE) (ix % xy_scale_use) / (DOUBLE) xy_scale_use;
        xdiff = fmod(ix_dbl, (DOUBLE) xy_scale_use) / (DOUBLE) xy_scale_use;
    }
    // y
    int lasty = ((pgrid->numy - 1) / xy_scale_use) * xy_scale_use; // reg index of end of last full casc grid cell
    if (iy > lasty) { // casc grid cell is truncated in y relative to cell of size xy_scale_use reg grid cells
        ydiff = (iy_dbl - (DOUBLE) lasty) / (DOUBLE) (pgrid->numy - 1 - lasty);
        //printf("xy_scale_use %d, ydiff %f = (DOUBLE) (iy %d - lasty %d) / (DOUBLE) (pgrid->numy %d - 1 - lasty %d)\n", xy_scale_use, ydiff, iy, lasty, pgrid->numy, lasty);
    } else {
        //ydiff = (DOUBLE) (iy % xy_scale_use) / (DOUBLE) xy_scale_use;
        ydiff = fmod(iy_dbl, (DOUBLE) xy_scale_use) / (DOUBLE) xy_scale_use;
    }
    //xdiff = (DOUBLE) (ix % xy_scale_use) / (DOUBLE) xy_scale_use;
    //ydiff = (DOUBLE) (iy % xy_scale_use) / (DOUBLE) xy_scale_use;
    // find regular grid z steps to reach upper limit
    int iz_test = iz;
    while (iz_test > 0 && pgrid->gridDesc_Cascading.zindex[iz_test - 1] == iz0_casc) {
        iz_test--;
    }
    // 20161026  zdiff = (DOUBLE) (iz - iz_test) / (DOUBLE) xy_scale_use;
    zdiff = (iz_dbl - (DOUBLE) iz_test) / (DOUBLE) xyz_scale;

    GRID_FLOAT_TYPE value;

    if (xdiff < 0.0 || xdiff > 1.0) {
        value = -VERY_LARGE_FLOAT;
        goto cleanup;
    }
    if (ydiff < 0.0 || ydiff > 1.0) {
        value = -VERY_LARGE_FLOAT;
        goto cleanup;
    }
    if (zdiff < 0.0 || zdiff > 1.0) {
        value = -VERY_LARGE_FLOAT;
        goto cleanup;
    }

    // location at grid node
    /*if (xdiff + ydiff + zdiff < SMALL_FLOAT) {
        value = ReadCascadingGrid3dValue(fpgrid, ix0_casc, iy0_casc, iz0_casc, pgrid);
        return (value);
    }*/

    DOUBLE vval000, vval001, vval010, vval011, vval100, vval101, vval110, vval111;
    vval000 = ReadCascadingGrid3dValue(fpgrid, ix0_casc_up, iy0_casc_up, iz0_casc, pgrid);
    vval001 = ReadCascadingGrid3dValue(fpgrid, ix0_casc_dn, iy0_casc_dn, iz1_casc, pgrid);
    vval010 = ReadCascadingGrid3dValue(fpgrid, ix0_casc_up, iy1_casc_up, iz0_casc, pgrid);
    vval011 = ReadCascadingGrid3dValue(fpgrid, ix0_casc_dn, iy1_casc_dn, iz1_casc, pgrid);
    vval100 = ReadCascadingGrid3dValue(fpgrid, ix1_casc_up, iy0_casc_up, iz0_casc, pgrid);
    vval101 = ReadCascadingGrid3dValue(fpgrid, ix1_casc_dn, iy0_casc_dn, iz1_casc, pgrid);
    vval110 = ReadCascadingGrid3dValue(fpgrid, ix1_casc_up, iy1_casc_up, iz0_casc, pgrid);
    vval111 = ReadCascadingGrid3dValue(fpgrid, ix1_casc_dn, iy1_casc_dn, iz1_casc, pgrid);

    // check for invalid / mask nodes
    if (vval000 < 0.0 || vval010 < 0.0 || vval100 < 0.0 || vval110 < 0.0
            || vval001 < 0.0 || vval011 < 0.0 || vval101 < 0.0 || vval111 < 0.0) {
        value = -VERY_LARGE_FLOAT;
        goto cleanup;
    }
    value = InterpCubeLagrange(xdiff, ydiff, zdiff,
            vval000, vval001, vval010, vval011,
            vval100, vval101, vval110, vval111);

cleanup:

    // clean up
    if (clean_casc_allocs) {
        if (allocated_grid_cascading) {
            FreeGrid_Cascading(pgrid);
        }
        if (created_grid_array_cascading) {
            DestroyGridArray(pgrid);
        }
    }

    return (value);
}

/** function to read cascading grid data from disk or array at upper xyz indices corresponding to regular grid index location
 *
 * 20161019 AJL - added
 *
 *  ix, iy, iz are indices in virtual, regular grid equivalent to this cascading grid
 */

GRID_FLOAT_TYPE ReadGrid3dValue_Cascading(FILE *fpgrid, int ix, int iy, int iz, GridDesc * pgrid) {

    // check indexes in range
    // ix, iy, iz are indices in virtual, regular grid equivalent to this cascading grid
    if (ix < 0 || ix >= pgrid->numx || iy < 0 || iy >= pgrid->numy || iz < 0 || iz >= pgrid->numz) {
        //nll_puterr("WARNING: grid file index out of range.");
        return (-VERY_LARGE_FLOAT);
    }

    // need grid array to find offset
    if (pgrid->array == NULL) {
        if (pgrid->buffer == NULL) {
            // prepare grid memory without allocating buffer, required to initialize cascading grid array
            AllocateGrid_Cascading(pgrid, 0);
        }
        pgrid->array = CreateGridArray_Cascading(pgrid);
    }

    // set x,y,z indices in cascading grid
    int iz_casc = pgrid->gridDesc_Cascading.zindex[iz];
    int ix_casc = ix / pgrid->gridDesc_Cascading.xyz_scale[iz];
    int iy_casc = iy / pgrid->gridDesc_Cascading.xyz_scale[iz];


    // get fvalue
    GRID_FLOAT_TYPE fvalue = ReadCascadingGrid3dValue(fpgrid, ix_casc, iy_casc, iz_casc, pgrid);

    return (fvalue);
}

/** function to read grid data from disk or array at index location ***/

GRID_FLOAT_TYPE ReadGrid3dValue(FILE *fpgrid, int ix, int iy, int iz, GridDesc * pgrid, int clean_casc_allocs) {

    if (isCascadingGrid(pgrid)) {
        //return (ReadGrid3dValue_Cascading(fpgrid, ix, iy, iz, pgrid));
        return (ReadGrid3dValue_Cascading_Interp(fpgrid, (double) ix, (double) iy, (double) iz, pgrid, clean_casc_allocs));
    }


    int numyz;
    long offset;
    GRID_FLOAT_TYPE fvalue;

    /* check indexes in range */

    if (ix < 0 || ix >= pgrid->numx || iy < 0 || iy >= pgrid->numy
            || iz < 0 || iz >= pgrid->numz) {
        //nll_puterr("WARNING: grid file index out of range.");
        return (-VERY_LARGE_FLOAT);
    }

    /* get fvalue */

    if (fpgrid != NULL) {
        /* calculate offset in bytes */
        numyz = pgrid->numy * pgrid->numz;
        offset = sizeof (GRID_FLOAT_TYPE) * (ix * numyz + iy * pgrid->numz + iz);
        fseek(fpgrid, offset, SEEK_SET);
        /* read fvalue */
        if (fread(&fvalue, sizeof (GRID_FLOAT_TYPE), 1, fpgrid) != 1) {
            sprintf(MsgStr,
                    "ERROR: reading grid value: %s: ix%d iy=%d iz=%d", pgrid->title, ix, iy, iz);
            nll_puterr(MsgStr);
            return (-VERY_LARGE_FLOAT);
        }
        if (pgrid->iSwapBytes)
            swapBytes(&fvalue, 1);
    } else {

        fvalue = ((GRID_FLOAT_TYPE ***) pgrid->array)[ix][iy][iz];
    }


    return (fvalue);
}

/** function to read grid data from disk or array at absolute location ***/

GRID_FLOAT_TYPE ReadAbsGrid3dValue(FILE *fpgrid, GridDesc* pgrid,
        double xloc, double yloc, double zloc, int ifloor) {
    int ix, iy, iz;
    GRID_FLOAT_TYPE fvalue;
    double shift;


    if (ifloor)
        shift = 0.0;
    else
        shift = 0.5;


    /* calculate nearest grid location */

    ix = (int) (shift + (xloc - pgrid->origx) / pgrid->dx);
    iy = (int) (shift + (yloc - pgrid->origy) / pgrid->dy);
    iz = (int) (shift + (zloc - pgrid->origz) / pgrid->dz);


    /* read grid data */

    fvalue = ReadGrid3dValue(fpgrid, ix, iy, iz, pgrid, 0);

    return (fvalue);
}



/** function to find value inside a cube using Lagrange interpolation***/

/* 	0.0 <= vvalKLM <= 1.0 */

/*	returns interp value at (xdiff, ydiff, zdiff) */

DOUBLE InterpCubeLagrange(DOUBLE xdiff, DOUBLE ydiff, DOUBLE zdiff,
        DOUBLE vval000, DOUBLE vval001, DOUBLE vval010, DOUBLE vval011,
        DOUBLE vval100, DOUBLE vval101, DOUBLE vval110, DOUBLE vval111) {

    DOUBLE value;
    DOUBLE oneMinusXdiff;
    DOUBLE oneMinusYdiff;
    DOUBLE oneMinusZdiff;

    oneMinusXdiff = 1.0 - xdiff;
    oneMinusYdiff = 1.0 - ydiff;
    oneMinusZdiff = 1.0 - zdiff;

    /*value = vval000 * (oneMinusXdiff) * oneMinusYdiff * (oneMinusZdiff)
            + vval001 * (oneMinusXdiff) * (oneMinusYdiff) * zdiff
            + vval010 * (oneMinusXdiff) * ydiff * (oneMinusZdiff)
            + vval011 * (oneMinusXdiff) * ydiff * zdiff
            + vval100 * xdiff * (oneMinusYdiff) * (oneMinusZdiff)
            + vval101 * xdiff * (oneMinusYdiff) * zdiff
            + vval110 * xdiff * ydiff * (oneMinusZdiff)
            + vval111 * xdiff * ydiff * zdiff;*/

    value = oneMinusXdiff * (
            oneMinusYdiff * (vval000 * oneMinusZdiff + vval001 * zdiff)
            + ydiff * (vval010 * oneMinusZdiff + vval011 * zdiff)
            )
            + xdiff * (
            oneMinusYdiff * (vval100 * oneMinusZdiff + vval101 * zdiff)
            + ydiff * (vval110 * oneMinusZdiff + vval111 * zdiff)
            );

    return (value);

}

/** function to find angles inside a cube */

float InterpCubeAngles(DOUBLE xdiff, DOUBLE ydiff, DOUBLE zdiff,
        DOUBLE vval000, DOUBLE vval001, DOUBLE vval010, DOUBLE vval011,
        DOUBLE vval100, DOUBLE vval101, DOUBLE vval110, DOUBLE vval111) {

    int nx, ny, nz;
    double azim[2][2][2], dip[2][2][2], azim_interp, dip_interp;
    int iqual[2][2][2], iqual_low;
    float value;
    TakeOffAngles angles;
    double azim_ref, azim_test;


    /* decode angles on vertices of cube */

    SetAnglesFloat(&angles, vval000);
    iqual[0][0][0] = GetTakeOffAngles(&angles,
            &(azim[0][0][0]), &(dip[0][0][0]), &(iqual[0][0][0]));
    SetAnglesFloat(&angles, vval001);
    iqual[0][0][1] = GetTakeOffAngles(&angles,
            &(azim[0][0][1]), &(dip[0][0][1]), &(iqual[0][0][1]));
    SetAnglesFloat(&angles, vval010);
    iqual[0][1][0] = GetTakeOffAngles(&angles,
            &(azim[0][1][0]), &(dip[0][1][0]), &(iqual[0][1][0]));
    SetAnglesFloat(&angles, vval011);
    iqual[0][1][1] = GetTakeOffAngles(&angles,
            &(azim[0][1][1]), &(dip[0][1][1]), &(iqual[0][1][1]));
    SetAnglesFloat(&angles, vval100);
    iqual[1][0][0] = GetTakeOffAngles(&angles,
            &(azim[1][0][0]), &(dip[1][0][0]), &(iqual[1][0][0]));
    SetAnglesFloat(&angles, vval101);
    iqual[1][0][1] = GetTakeOffAngles(&angles,
            &(azim[1][0][1]), &(dip[1][0][1]), &(iqual[1][0][1]));
    SetAnglesFloat(&angles, vval110);
    iqual[1][1][0] = GetTakeOffAngles(&angles,
            &(azim[1][1][0]), &(dip[1][1][0]), &(iqual[1][1][0]));
    SetAnglesFloat(&angles, vval111);
    iqual[1][1][1] = GetTakeOffAngles(&angles,
            &(azim[1][1][1]), &(dip[1][1][1]), &(iqual[1][1][1]));


    // check for lowest quality angles
    // 20130823 AJL - bug fix:
    //    correct azimuths to avoid discontinuity at 0/360 deg

    iqual_low = 999;
    azim_ref = azim[0][0][0];
    for (nx = 0; nx < 2; nx++) {
        for (ny = 0; ny < 2; ny++) {
            for (nz = 0; nz < 2; nz++) {
                //printf("%d %d %d az %f dip %f q %d\n", nx, ny, nz, azim[nx][ny][nz], dip[nx][ny][nz], iqual[nx][ny][nz]);
                if (iqual[nx][ny][nz] < iqual_low)
                    iqual_low = iqual[nx][ny][nz];
                azim_test = azim[nx][ny][nz] - azim_ref;
                if (azim_test < -90.0) {
                    azim[nx][ny][nz] += 360.0;
                } else if (azim_test > 90.0) {
                    azim[nx][ny][nz] -= 360.0;
                }
            }
        }
    }


    /* determine angles to return */

    if (iqual_low < ANGLE_QUALITY_CUTOFF) {
        /* if lowest quality is too low, use nearest node */
        value = vval000;
    } else {
        /* otherwise interpolate */
        azim_interp = InterpCubeLagrange(xdiff, ydiff, zdiff,
                azim[0][0][0], azim[0][0][1], azim[0][1][0],
                azim[0][1][1], azim[1][0][0], azim[1][0][1],
                azim[1][1][0], azim[1][1][1]);
        // 20130823 AJL - bug fix:
        if (azim_interp < 0.0) {
            azim_interp += 360.0;
        } else if (azim_interp > 360.0) {

            azim_interp -= 360.0;
        }
        dip_interp = InterpCubeLagrange(xdiff, ydiff, zdiff,
                dip[0][0][0], dip[0][0][1], dip[0][1][0],
                dip[0][1][1], dip[1][0][0], dip[1][0][1],
                dip[1][1][0], dip[1][1][1]);
        angles = SetTakeOffAngles(azim_interp, dip_interp, iqual_low);
        value = angles.fval;
    }


    return (value);

}




/** function to read grid data from disk or buffer at absolute location
                with interpolation ***/

/* read from file if fpgrid != NULL, otherwise read from grid buffer */

GRID_FLOAT_TYPE ReadAbsInterpGrid3d(FILE *fpgrid, GridDesc* pgrid, double xloc, double yloc, double zloc, int clean_casc_allocs) {

    DOUBLE xoff, yoff, zoff;
    xoff = (xloc - pgrid->origx) / pgrid->dx;
    yoff = (yloc - pgrid->origy) / pgrid->dy;
    zoff = (zloc - pgrid->origz) / pgrid->dz;

    // check if cascading grid
    if (isCascadingGrid(pgrid)) {
        return (ReadGrid3dValue_Cascading_Interp(fpgrid, xoff, yoff, zoff, pgrid, clean_casc_allocs));
    }

    int ix0, ix1, iy0, iy1, iz0, iz1;
    GRID_FLOAT_TYPE value;
    DOUBLE vval000, vval001, vval010, vval011, vval100, vval101, vval110, vval111;
    DOUBLE xdiff, ydiff, zdiff;
    int numx, numy, numz, numyz;
    GRID_FLOAT_TYPE *buffer;

    buffer = (GRID_FLOAT_TYPE *) pgrid->buffer;
    numx = pgrid->numx;
    numy = pgrid->numy;
    numz = pgrid->numz;
    numyz = numy * numz;


    /* calculate grid locations on edge of solid containing point */

    ix0 = (int) (xoff - VERY_SMALL_DOUBLE);
    iy0 = (int) (yoff - VERY_SMALL_DOUBLE);
    iz0 = (int) (zoff - VERY_SMALL_DOUBLE);

    ix1 = (ix0 < numx - 1) ? ix0 + 1 : ix0;
    iy1 = (iy0 < numy - 1) ? iy0 + 1 : iy0;
    iz1 = (iz0 < numz - 1) ? iz0 + 1 : iz0;

    /*	if (ix1 < 0 || ix1 >= numx) fprintf(stderr, "GRID INDEX ERROR: ix1 %d (%d-%d) xloc %lf ix0 %d ix1 %d\n", ix1, 0, numx, xloc, ix0, ix1);
            if (iy1 < 0 || iy1 >= numy) fprintf(stderr, "GRID INDEX ERROR: iy1 %d (%d-%d) yloc %lf iy0 %d iy1 %d\n", iy1, 0, numy, yloc, iy0, iy1);
            if (iz1 < 0 || iz1 >= numz) fprintf(stderr, "GRID INDEX ERROR: iz1 %d (%d-%d) zloc %lf iz0 %d iz1 %d\n", iz1, 0, numz, zloc, iz0, iz1);
            if (ix0 < 0 || ix0 >= numx) fprintf(stderr, "GRID INDEX ERROR: ix0 %d (%d-%d) xloc %lf ix0 %d ix1 %d\n", ix0, 0, numx, xloc, ix0, ix1);
            if (iy0 < 0 || iy0 >= numy) fprintf(stderr, "GRID INDEX ERROR: iy0 %d (%d-%d) yloc %lf iy0 %d iy1 %d\n", iy0, 0, numy, yloc, iy0, iy1);
            if (iz0 < 0 || iz0 >= numz) fprintf(stderr, "GRID INDEX ERROR: iz0 %d (%d-%d) zloc %lf iz0 %d iz1 %d\n", iz0, 0, numz, zloc, iz0, iz1);
     */

    /*	if (ix1 < 0 || ix1 >= numx) ix1 = ix0;
            if (iy1 < 0 || iy1 >= numy) iy1 = iy0;
            if (iz1 < 0 || iz1 >= numz) iz1 = iz0;

            if (ix0 < 0 || ix0 >= numx) ix0 = ix1;
            if (iy0 < 0 || iy0 >= numy) iy0 = iy1;
            if (iz0 < 0 || iz0 >= numz) iz0 = iz1;
     */

    xdiff = xoff - (DOUBLE) ix0;
    ydiff = yoff - (DOUBLE) iy0;
    zdiff = zoff - (DOUBLE) iz0;

    if (xdiff < 0.0 || xdiff > 1.0)
        return (-VERY_LARGE_FLOAT);
    if (ydiff < 0.0 || ydiff > 1.0)
        return (-VERY_LARGE_FLOAT);
    if (zdiff < 0.0 || zdiff > 1.0)
        return (-VERY_LARGE_FLOAT);

    /* location at grid node */

    if (xdiff + ydiff + zdiff < SMALL_FLOAT) {
        if (fpgrid != NULL)
            value = ReadGrid3dValue(fpgrid, ix0, iy0, iz0, pgrid, 0);
        else
            value = *(buffer + ix0 * numyz + iy0 * numz + iz0);
        return (value);
    }


    /* read vertex values from grid file or array */

    if (fpgrid != NULL) {
        vval000 = ReadGrid3dValue(fpgrid, ix0, iy0, iz0, pgrid, 0);
        vval001 = ReadGrid3dValue(fpgrid, ix0, iy0, iz1, pgrid, 0);
        vval010 = ReadGrid3dValue(fpgrid, ix0, iy1, iz0, pgrid, 0);
        vval011 = ReadGrid3dValue(fpgrid, ix0, iy1, iz1, pgrid, 0);
        vval100 = ReadGrid3dValue(fpgrid, ix1, iy0, iz0, pgrid, 0);
        vval101 = ReadGrid3dValue(fpgrid, ix1, iy0, iz1, pgrid, 0);
        vval110 = ReadGrid3dValue(fpgrid, ix1, iy1, iz0, pgrid, 0);
        vval111 = ReadGrid3dValue(fpgrid, ix1, iy1, iz1, pgrid, 0);
    } else {
        /*
        vval000 = pgrid->array[ix0][iy0][iz0];
        vval001 = pgrid->array[ix0][iy0][iz1];
        vval010 = pgrid->array[ix0][iy1][iz0];
        vval011 = pgrid->array[ix0][iy1][iz1];
        vval100 = pgrid->array[ix1][iy0][iz0];
        vval101 = pgrid->array[ix1][iy0][iz1];
        vval110 = pgrid->array[ix1][iy1][iz0];
        vval111 = pgrid->array[ix1][iy1][iz1];
         */
        vval000 = *(buffer + ix0 * numyz + iy0 * numz + iz0);
        vval001 = *(buffer + ix0 * numyz + iy0 * numz + iz1);
        vval010 = *(buffer + ix0 * numyz + iy1 * numz + iz0);
        vval011 = *(buffer + ix0 * numyz + iy1 * numz + iz1);
        vval100 = *(buffer + ix1 * numyz + iy0 * numz + iz0);
        vval101 = *(buffer + ix1 * numyz + iy0 * numz + iz1);
        vval110 = *(buffer + ix1 * numyz + iy1 * numz + iz0);
        vval111 = *(buffer + ix1 * numyz + iy1 * numz + iz1);
    }

    /* interpolate values */

    if (pgrid->type == GRID_ANGLE || pgrid->type == GRID_ANGLE_2D) {
        value = InterpCubeAngles(xdiff, ydiff, zdiff,
                vval000, vval001, vval010, vval011,
                vval100, vval101, vval110, vval111);
    } else {
        // INGV
        // check for invalid / mask nodes
        if (vval000 < 0.0 || vval010 < 0.0 || vval100 < 0.0 || vval110 < 0.0
                || vval001 < 0.0 || vval011 < 0.0 || vval101 < 0.0 || vval111 < 0.0)

            return (-VERY_LARGE_FLOAT);
        value = InterpCubeLagrange(xdiff, ydiff, zdiff,
                vval000, vval001, vval010, vval011,
                vval100, vval101, vval110, vval111);
    }


    return (value);
}


/** function to read grid data from disk or buffer at absolute location with interpolation ***/

/* 2D version - ix assumed = 0 */

/* read from file if fpgrid != NULL, otherwise read from grid buffer */

DOUBLE ReadAbsInterpGrid2d(FILE *fpgrid, GridDesc* pgrid, double yloc, double zloc) {

    int ix0, ix1, iy0, iy1, iz0, iz1;
    DOUBLE value;
    DOUBLE vval00, vval01, vval10, vval11;

    DOUBLE yoff, zoff;
    DOUBLE ydiff, zdiff;

    /* calculate grid locations on edge of solid containing point */

    ix0 = 0;
    yoff = (yloc - pgrid->origy) / pgrid->dy;
    iy0 = (int) (yoff - VERY_SMALL_DOUBLE);
    zoff = (zloc - pgrid->origz) / pgrid->dz;
    iz0 = (int) (zoff - VERY_SMALL_DOUBLE);

    ix1 = 0;
    iy1 = (iy0 < pgrid->numy - 1) ? iy0 + 1 : iy0;
    iz1 = (iz0 < pgrid->numz - 1) ? iz0 + 1 : iz0;


    ydiff = yoff - (DOUBLE) iy0;
    zdiff = zoff - (DOUBLE) iz0;

    // DEBUG
    /*
    printf("%d %d %d  %lf %lf %lf  %lf %lf %lf %s\n",
            pgrid->numx, pgrid->numy, pgrid->numz,
            pgrid->origx, pgrid->origy, pgrid->origz,
            pgrid->dx, pgrid->dy, pgrid->dz, pgrid->chr_type);
    if (iy0 < 0 || iy1 >= pgrid->numy )
            printf("ERROR: ReadAbsInterpGrid2d: (iy0 < 0 || iy1 >= pgrid->numy )\n");
    if (iz0 < 0 || iz1 >= pgrid->numz )
            printf("ERROR: ReadAbsInterpGrid2d: (iz0 < 0 || iz1 >= pgrid->numz )\n");
    if (ydiff < 0.0 || ydiff > 1.0)
            printf("ERROR: ReadAbsInterpGrid2d: (ydiff < 0.0 || ydiff > 1.0 )\n");
    if (zdiff < 0.0 || zdiff > 1.0)
            printf("ERROR: ReadAbsInterpGrid2d: (zdiff < 0.0 || zdiff > 1.0 )\n");
     */
    // End - DEBUG

    //INGV
    if (iy0 < 0 || iy1 >= pgrid->numy)
        return (-VERY_LARGE_DOUBLE);
    if (iz0 < 0 || iz1 >= pgrid->numz)
        return (-VERY_LARGE_DOUBLE);

    if (ydiff < 0.0 || ydiff > 1.0)
        return (-VERY_LARGE_DOUBLE);
    if (zdiff < 0.0 || zdiff > 1.0)
        return (-VERY_LARGE_DOUBLE);

    /* location at grid node */

    if (ydiff + zdiff < SMALL_FLOAT) {
        if (fpgrid != NULL)
            value = ReadGrid3dValue(fpgrid, ix0, iy0, iz0, pgrid, 0);
        else
            value = ((GRID_FLOAT_TYPE ***) pgrid->array)[ix0][iy0][iz0];
        return (value);
    }


    /* read vertex values from grid file or array */

    if (fpgrid != NULL) {
        vval00 = ReadGrid3dValue(fpgrid, ix0, iy0, iz0, pgrid, 0);
        vval01 = ReadGrid3dValue(fpgrid, ix0, iy0, iz1, pgrid, 0);
        vval10 = ReadGrid3dValue(fpgrid, ix0, iy1, iz0, pgrid, 0);
        vval11 = ReadGrid3dValue(fpgrid, ix0, iy1, iz1, pgrid, 0);
    } else {
        vval00 = ((GRID_FLOAT_TYPE ***) pgrid->array)[ix0][iy0][iz0];
        vval01 = ((GRID_FLOAT_TYPE ***) pgrid->array)[ix0][iy0][iz1];
        vval10 = ((GRID_FLOAT_TYPE ***) pgrid->array)[ix0][iy1][iz0];
        vval11 = ((GRID_FLOAT_TYPE ***) pgrid->array)[ix0][iy1][iz1];
    }

    // DEBUG
    /*
    if (vval00 < 0.0 || vval01 < 0.0 || vval10 < 0.0 || vval11 < 0.0)
            printf("ERROR: ReadAbsInterpGrid2d: (vval00 < 0.0 || vval01 < 0.0 || vval10 < 0.0 || vval11 < 0.0)\n");
     */
    // End - DEBUG

    // INGV
    // check for invalid / mask nodes
    if (vval00 < 0.0 || vval01 < 0.0 || vval10 < 0.0 || vval11 < 0.0)
        return (-VERY_LARGE_DOUBLE);

    /* interpolate values */

    value = InterpSquareLagrange(ydiff, zdiff,
            vval00, vval01, vval10, vval11);

    return (value);
}


/** function to find value inside a square using Lagrange interpolation */

/* 	0.0 <= vvalKLM <= 1.0 */

/*	returns interp value at (xdiff, zdiff) */

DOUBLE InterpSquareLagrange(DOUBLE xdiff, DOUBLE zdiff,
        DOUBLE vval00, DOUBLE vval01, DOUBLE vval10, DOUBLE vval11) {

    DOUBLE value;

    value = vval00 * (1.0 - xdiff) * (1.0 - zdiff)
            + vval01 * (1.0 - xdiff) * zdiff
            + vval10 * xdiff * (1.0 - zdiff)
            + vval11 * xdiff * zdiff;

    return (value);

}

/** function to write hypocenter/arrivals to output */

int WriteLocation(FILE *fpio, HypoDesc* phypo, ArrivalDesc* parrivals,
        int narrivals, char* filename,
        int iWriteArrivals, int iWriteEndLoc, int iWriteMinimal,
        GridDesc* pgrid, int n_proj) {

    return (_WriteLocation(fpio, phypo, parrivals,
            narrivals, filename,
            iWriteArrivals, iWriteEndLoc, iWriteMinimal,
            pgrid, n_proj, IO_ARRIVAL_ALL));
}

/** function to write arrivals to output */

int WritePhases(FILE *fpio, HypoDesc* phypo, ArrivalDesc* parrivals,
        int narrivals, char* filename,
        int iWriteArrivals, int iWriteEndLoc, int iWriteMinimal,
        GridDesc* pgrid, int n_proj, int io_arrival_mode) {

    return (_WriteLocation(fpio, phypo, parrivals,
            narrivals, filename,
            iWriteArrivals, iWriteEndLoc, iWriteMinimal,
            pgrid, n_proj, io_arrival_mode));
}

/** function to write hypocenter/arrivals to output */

int _WriteLocation(FILE *fpio, HypoDesc* phypo, ArrivalDesc* parrivals,
        int narrivals, char * filename,
        int iWriteArrivals, int iWriteEndLoc, int iWriteMinimal,
        GridDesc* pgrid, int n_proj, int io_arrival_mode) {

    int ifile = 0, narr;
    ArrivalDesc* parr;


    /* write hypocenter to file */

    if (fpio == NULL) {
        if ((fpio = fopen(filename, "w")) == NULL) {
            nll_puterr2("ERROR: opening hypocenter output file", filename);
            return (-1);
        }
        NumFilesOpen++;
        ifile = 1;
    }

    /* write hypocenter parameters */

    if (iWriteMinimal) {
        fprintf(fpio, "NLLOC \"%s\" \"%s\" \" \"\n",
                phypo->fileroot, phypo->locStat);
    } else {
        fprintf(fpio, "NLLOC \"%s\" \"%s\" \"%s\"\n",
                phypo->fileroot, phypo->locStat, phypo->locStatComm);
        fprintf(fpio, "PUBLIC_ID %s\n", phypo->public_id); // 20190823 AJL - added
        fprintf(fpio, "SIGNATURE \"%s\"\n", phypo->signature);
        fprintf(fpio, "COMMENT \"%s\"\n", phypo->comment);
        if (pgrid != NULL)
            fprintf(fpio, "GRID  %d %d %d  %lg %lg %lg  %lg %lg %lg %s\n",
                pgrid->numx, pgrid->numy, pgrid->numz,
                pgrid->origx, pgrid->origy, pgrid->origz,
                pgrid->dx, pgrid->dy, pgrid->dz, pgrid->chr_type);
        else
            fprintf(fpio, "GRID  %d %d %d  %lg %lg %lg  %lg %lg %lg %s\n",
                -1, -1, -1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, "NULLGRID");
        fprintf(fpio, "SEARCH %s\n", phypo->searchInfo);
        fprintf(fpio,
                "HYPOCENTER  x %lg y %lg z %lg  OT %lg  ix %d iy %d iz %d %s\n",
                phypo->x, phypo->y, phypo->z, (double) phypo->sec,
                phypo->ix, phypo->iy, phypo->iz,
                phypo->type // 20170811 AJL - added to allow saving of expectation hypocenter results instead of maximum likelihood
                );
    }
    fprintf(fpio,
            "GEOGRAPHIC  OT %4.4d %2.2d %2.2d  %2.2d %2.2d %lf  Lat %lf Long %lf Depth %lf\n",
            // 20151224  "GEOGRAPHIC  OT %4.4d %2.2d %2.2d  %2.2d %2.2d %f  Lat %f Long %f Depth %lg\n",
            //"GEOGRAPHIC  OT %4.4d %2.2d %2.2d  %2.2d %2.2d %lg  Lat %lg Long %lg Depth %lg\n",
            phypo->year, phypo->month, phypo->day,
            phypo->hour, phypo->min, (double) phypo->sec,
            phypo->dlat, phypo->dlong, phypo->depth);
    fprintf(fpio,
            "QUALITY  Pmax %Lg MFmin %lg MFmax %lg RMS %lg Nphs %d Gap %lg Dist %lg Mamp %5.2lf %d Mdur %5.2lf %d\n",
            phypo->probmax, phypo->misfit, phypo->grid_misfit_max,
            phypo->rms, phypo->nreadings, phypo->gap,
            GeometryMode == MODE_GLOBAL ? phypo->dist * KM2DEG : phypo->dist,
            phypo->amp_mag, phypo->num_amp_mag,
            phypo->dur_mag, phypo->num_dur_mag
            );

    fprintf(fpio,
            "VPVSRATIO  VpVsRatio %lg  Npair %d  Diff %lg\n",
            phypo->VpVs, phypo->nVpVs, phypo->tsp_min_max_diff
            );

    if (!iWriteMinimal) {

        fprintf(fpio,
                // 20190903 AJL  "STATISTICS  ExpectX %lg Y %lg Z %lg",
                "STATISTICS  ExpectX %f Y %f Z %f", // increase precision, important for spherical (lat, lon) output
                phypo->expect.x, phypo->expect.y, phypo->expect.z);
        fprintf(fpio, "  CovXX %lg XY %lg XZ %lg",
                phypo->cov.xx, phypo->cov.xy, phypo->cov.xz);
        fprintf(fpio, " YY %lg YZ %lg",
                phypo->cov.yy, phypo->cov.yz);
        fprintf(fpio, " ZZ %lg",
                phypo->cov.zz);
        fprintf(fpio, " EllAz1  %lg Dip1  %lg Len1  %lg",
                phypo->ellipsoid.az1, phypo->ellipsoid.dip1,
                phypo->ellipsoid.len1);
        fprintf(fpio, " Az2  %lg Dip2  %lg Len2  %lg",
                phypo->ellipsoid.az2, phypo->ellipsoid.dip2,
                phypo->ellipsoid.len2);
        fprintf(fpio, " Len3  %le\n", phypo->ellipsoid.len3);

        fprintf(fpio,
                "STAT_GEOG  ExpectLat %lf Long %lf Depth %lf\n",
                // 20151224  "STAT_GEOG  ExpectLat %f Long %f Depth %f\n",
                phypo->expect_dlat, phypo->expect_dlong, phypo->expect.z);

        // 20170811 AJL - added to allow saving of expectation hypocenter results instead of maximum likelihood
        if (strcmp(phypo->type, HYPO_TYPE_EXPECTATION) == 0) {
            fprintf(fpio,
                    "MAXIMUM_LIKELIHOOD  MaxLikeLat %lf Long %lf Depth %lf OT %lf\n",
                    phypo->max_like_dlat, phypo->max_like_dlong, phypo->max_like.z, (double) phypo->max_like_sec);
        }

        fprintf(fpio, "%s\n", MapProjStr[n_proj]);

        // 20100519 AJL Added new hypo line
        /*
        int associatedPhaseCount; //   - Number of associated phases, regardless of their use for origin computation.
        // [->nreadings] int usedPhaseCount; // QML - Number of defining phases, i. e., phase observations that were actually used for computing
        // the origin. Note that there may be more than one defining phase per station.
        int associatedStationCount; // QML - Number of stations at which the event was observed.
        int usedStationCount; // QML - Number of stations from which data was used for origin computation.
        int depthPhaseCount; // QML - Number of depth phases (typically pP, sometimes sP) used in depth computation.
        // [->rms] double standardError; // QML - RMS of the travel time residuals of the arrivals used for the origin computation. Unit: s
        // [->gap] double azimuthalGap; // QML - Largest azimuthal gap in station distribution as seen from epicenter. Unit: deg
        // [->gap_secondary] double secondaryAzimuthalGap; // QML - Secondary azimuthal gap in station distribution, i. e., the largest azimuthal gap a station closes. Unit: deg
        char groundTruthLevel[8]; // QML - String describing ground-truth level, e. g. GT0, GT5, etc.
        double minimumDistance; // QML - Epicentral distance of station closest to the epicenter. Unit: km
        double maximumDistance; // QML - Epicentral distance of station farthest from the epicenter. Unit: km
        double medianDistance; // QML - Median epicentral distance of used stations. Unit: km
         */
        fprintf(fpio,
                "QML_OriginQuality  assocPhCt %d  usedPhCt %d  assocStaCt %d  usedStaCt %d  depthPhCt %d  stdErr %lg  azGap %lg  secAzGap %lg",
                phypo->associatedPhaseCount, phypo->nreadings, phypo->associatedStationCount, phypo->usedStationCount, phypo->depthPhaseCount,
                phypo->rms, phypo->gap, phypo->gap_secondary
                );
        fprintf(fpio,
                "  gtLevel %s  minDist %lg maxDist %lg medDist %lg\n",
                phypo->groundTruthLevel, phypo->minimumDistance, phypo->maximumDistance, phypo->medianDistance
                );

        // 20100519 AJL Added new hypo line
        /*
        // QML fields added for compatibility with QuakeML OriginUncertainty attributes (AJL 201005)
        // preferredDescription enum 0..1 – Element (OriginUncertaintyDescription)
        double horizontalUncertainty; // QML - Circular confidence region, given by single value of horizontal uncertainty. Unit: km
        double minHorizontalUncertainty; // QML - Semi-major axis of confidence ellipse. Unit: km
        double maxHorizontalUncertainty; // QML - Semi-minor axis of confidence ellipse. Unit: km
        double azimuthMaxHorizontalUncertainty; // QML - Azimuth of major axis of confidence ellipse. Unit: km
         */
        double azMaxHorUnc = phypo->ellipse.az1 + 90.0;
        if (azMaxHorUnc >= 360.0)
            azMaxHorUnc -= 360.0;
        if (azMaxHorUnc >= 180.0)
            azMaxHorUnc -= 180.0;
        fprintf(fpio,
                "QML_OriginUncertainty  horUnc %lg  minHorUnc %lg  maxHorUnc %lg  azMaxHorUnc %lg\n",
                -1.0, // 20100617 AJL - horizontalUncertainty: not clear what this is, set undefined
                phypo->ellipse.len1, phypo->ellipse.len2, azMaxHorUnc
                );


        // 20150602 AJL Added new hypo line
        /*
        // QML fields added for compatibility with QuakeML ConfidenceEllipsoid attributes (AJL 201005)
        semiMajorAxisLength Largest uncertainty, corresponding to the semi-major axis of the confidence ellipsoid. Unit: m
        semiMinorAxisLength Smallest uncertainty, corresponding to the semi-minor axis of the confidence ellipsoid. Unit: m
        semiIntermediateAxisLength Uncertainty in direction orthogonal to major and minor axes of the confidence ellipsoid. Unit: m
        majorAxisPlunge Plunge angle of major axis of confidence ellipsoid. Corresponds to Tait-Bryan angle phi. Unit: deg
        majorAxisAzimuth Azimuth angle of major axis of confidence ellipsoid. Corresponds to Tait-Bryan angle psi. Unit: deg
        majorAxisRotation This angle describes a rotation about the confidence ellipsoid’s major axis which is required
           to define the direction of the ellipsoid’s minor axis. Corresponds to Tait-Bryan angle θ. Unit: deg
         */
        double semiMajorAxisLength;
        double semiMinorAxisLength;
        double semiIntermediateAxisLength;
        double majorAxisAzimuth;
        double majorAxisPlunge;
        double majorAxisRotation;
        if (nllEllipsiod2QMLConfidenceEllipsoid(
                &(phypo->ellipsoid),
                &semiMajorAxisLength, &semiMinorAxisLength, &semiIntermediateAxisLength,
                &majorAxisAzimuth, &majorAxisPlunge, &majorAxisRotation) < 0) {
            semiMajorAxisLength = -1.0;
            semiMinorAxisLength = -1.0;
            semiIntermediateAxisLength = -1.0;
            majorAxisAzimuth = -1.0;
            majorAxisPlunge = -1.0;
            majorAxisRotation = -1.0;
        }
        fprintf(fpio,
                "QML_ConfidenceEllipsoid  semiMajorAxisLength %lg  semiMinorAxisLength %lg  semiIntermediateAxisLength %lg  majorAxisPlunge %lg  majorAxisAzimuth %lg  majorAxisRotation %lg\n",
                semiMajorAxisLength, semiMinorAxisLength, semiIntermediateAxisLength, majorAxisPlunge, majorAxisAzimuth, majorAxisRotation
                );


        // 20100617 AJL Added new hypo line
        // SED-ETH fields added for compatibility with legacy SED location quality indicators (AJL 201006)
        // SED_Origin errx 33.0775  erry 7.5298  errz 27.5298  diffMaxLikeExpect 4.9584 quality B
        if (phypo->qualitySED != '\0') {
            fprintf(fpio,
                    "SED_Origin  errx %lg  erry %lg  errz %lg  diffMaxLikeExpect %lg  quality %c\n",
                    sqrt(phypo->cov.xx), sqrt(phypo->cov.yy), sqrt(phypo->cov.zz), phypo->diffMaxLikeExpect, phypo->qualitySED
                    );
        }


        /* write mechanism */
        fprintf(fpio,
                "FOCALMECH  Hyp  %f %f %f",
                phypo->dlat, phypo->dlong, phypo->depth);
        fprintf(fpio,
                " Mech  %lg %lg %lg",
                phypo->focMech.dipDir, phypo->focMech.dipAng,
                phypo->focMech.rake);
        fprintf(fpio, " mf  %lg nObs %d",
                phypo->focMech.misfit, phypo->focMech.nObs);
        fprintf(fpio, "\n");


        /* write differential loc parameters */
        if (nll_mode == MODE_DIFFERENTIAL
                || phypo->event_id >= 0) { // 20110620 AJL - preserve event id if available
            fprintf(fpio,
                    "DIFFERENTIAL  Nhyp %ld", phypo->event_id);
            fprintf(fpio, "\n");
        }


    }


    /* write arrival parameters */
    /* !!! Remember to modify the sscanf in function GetHypLoc() when modifying this fprintf */

    if (iWriteArrivals) {

        fprintf(fpio,
                "PHASE ID Ins Cmp On Pha  FM Date     HrMn   Sec     Err  ErrMag    Coda      Amp       Per");

        if (PhaseFormat == FORMAT_PHASE_2)
            fprintf(fpio,
                "       PriorWt");

        if (io_arrival_mode == IO_ARRIVAL_ALL) {
            fprintf(fpio,
                    "  >   TTpred    Res       Weight    StaLoc(X  Y         Z)        SDist    SAzim  RAz  RDip RQual    Tcorr ");
            if (PhaseFormat == FORMAT_PHASE_2)
                fprintf(fpio,
                    "      TTerr");

        }
        fprintf(fpio, "\n");

        for (narr = 0; narr < narrivals; narr++) {
            parr = parrivals + narr;
            if (nll_mode == MODE_DIFFERENTIAL) {
                //printf("phypo->event_id %ld  parr->dd_event_id_1 %ld  parr->dd_event_id_2 %ld\n", phypo->event_id, parr->dd_event_id_1, parr->dd_event_id_2);
                if (parr->flag_ignore || !(phypo->event_id == parr->dd_event_id_1
                        || phypo->event_id == parr->dd_event_id_2))
                    continue;
            }
            WriteArrival(fpio, parr, io_arrival_mode);
        }
        fprintf(fpio, "END_PHASE\n");
    }


    /* write end line and blank line */
    if (iWriteEndLoc)
        fprintf(fpio, "END_NLLOC\n\n");

    if (ifile) {

        fclose(fpio);
        NumFilesOpen--;
    }

    return (0);
}

/** function to read hypocenter/arrival parameters */

int GetHypLoc(FILE *fpio, const char* filein, HypoDesc* phypo,
        ArrivalDesc* parrivals, int *pnarrivals, int iReadArrivals,
        GridDesc* pgrid, int n_proj) {

    int istat, ifile = 0;
    int lineLength;
    char fn_in[FILENAME_MAX];
    char line[MAXLINE_LONG], *pstr, *pstr2 = NULL;
    double hypo_sec, templat, templong;
    ArrivalDesc* parr;


    /* open hypocenter file */

    if (fpio == NULL) {
        if ((pstr = strstr(filein, ".hyp")) == NULL
                || ((pstr - filein) < strlen(filein) - 4))
            sprintf(fn_in, "%s.hyp", filein);
        else
            sprintf(fn_in, "%s", filein);
        if ((fpio = fopen(fn_in, "r")) == NULL) {
            nll_puterr2("ERROR: opening hypocenter file", fn_in);
            return (-1);
        }
        NumFilesOpen++;
        ifile = 1;
    }


    /* read hypocenter parameters */

    /* search for first line of hypocenter description */
    do {
        if (fgets(line, MAXLINE_LONG, fpio) == NULL)
            goto eof_exit;
    } while (strncmp("NLLOC", line, 5));
    phypo->fileroot[0] = '\0';
    if ((pstr = strchr(line, '"')) != NULL) {
        if ((pstr2 = strchr(pstr + 1, '"')) != NULL) {
            strncpy(phypo->fileroot, pstr + 1, pstr2 - pstr - 1);
            phypo->fileroot[pstr2 - pstr - 1] = '\0';
        }
    }
    phypo->locStat[0] = '\0';
    if ((pstr = strchr(pstr2 + 1, '"')) != NULL) {
        if ((pstr2 = strchr(pstr + 1, '"')) != NULL) {
            strncpy(phypo->locStat, pstr + 1, pstr2 - pstr - 1);
            phypo->locStat[pstr2 - pstr - 1] = '\0';
        }
    }
    phypo->locStatComm[0] = '\0';
    if ((pstr = strchr(pstr2 + 1, '"')) != NULL) {
        if ((pstr2 = strchr(pstr + 1, '"')) != NULL) {
            strncpy(phypo->locStatComm, pstr + 1,
                    pstr2 - pstr - 1);
            phypo->locStatComm[pstr2 - pstr - 1] = '\0';
        }
    }


    /* read hypocenter description lines until END_NLLOC reached */

    while (1) {

        /* read next line */
        if (fgets(line, MAXLINE_LONG, fpio) == NULL)
            goto eof_exit;

        /* end of NLLoc hypocenter descrition */
        if (strncmp(line, "END_NLLOC", 9) == 0)
            break;

        lineLength = TrimString(line);

        /* blank line */
        if (lineLength == 0)
            continue;

        /* identify and read line */

        if (strncmp(line, "PUBLIC_ID", 9) == 0) {
            // 20190823 AJL - added
            if (sscanf(line, "%*s %s", phypo->public_id) == EOF)
                goto eof_exit;
        } else if (strncmp(line, "SIGNATURE", 9) == 0) {
            /* SIGNATURE */
            strcpy(phypo->signature, strchr(line, '"') + 1);
            *(strchr(phypo->signature, '"')) = '\0';
        } else if (strncmp(line, "COMMENT", 7) == 0) {
            /* COMMENT */
            strcpy(phypo->comment, strchr(line, '"') + 1);
            *(strchr(phypo->comment, '"')) = '\0';
        } else if (strncmp(line, "GRID", 4) == 0) {
            /* GRID */
            if (pgrid != NULL) {
                if (sscanf(line,
                        "%*s  %d %d %d  %lf %lf %lf  %lf %lf %lf %s",
                        &pgrid->numx, &pgrid->numy, &pgrid->numz,
                        &pgrid->origx, &pgrid->origy, &pgrid->origz,
                        &pgrid->dx, &pgrid->dy, &pgrid->dz,
                        pgrid->chr_type)
                        == EOF)
                    goto eof_exit;
            }
        } else if (strncmp(line, "SEARCH", 6) == 0) {
            /* SEARCH */
            line[MAXLINE_LONG - 1] = '\0';
            strcpy(phypo->searchInfo, strchr(line, ' ') + 1);
        } else if (strncmp(line, "HYPOCENTER", 10) == 0) {
            /* HYPOCENTER */
            if (sscanf(line,
                    "%*s %*s %lf %*s %lf %*s %lf %*s %lf %*s %d %*s %d %*s %d %s",
                    &(phypo->x), &(phypo->y), &(phypo->z),
                    &hypo_sec,
                    &(phypo->ix), &(phypo->iy), &(phypo->iz),
                    phypo->type // 20170811 AJL - added to allow saving of expectation hypocenter results instead of maximum likelihood
                    )
                    == EOF)
                goto eof_exit;
            phypo->sec = (long double) hypo_sec;
        } else if (strncmp(line, "GEOGRAPHIC", 10) == 0) {
            /* GEOGRAPHIC */
            if (sscanf(line,
                    "%*s %*s %d %d %d   %d %d %lf %*s %lf %*s %lf %*s %lf",
                    &(phypo->year), &(phypo->month), &(phypo->day),
                    &(phypo->hour), &(phypo->min), &hypo_sec,
                    &(phypo->dlat), &(phypo->dlong),
                    &(phypo->depth)) == EOF)
                goto eof_exit;
            phypo->sec = (long double) hypo_sec;
        } else if (strncmp(line, "QUALITY", 7) == 0) {
            /* QUALITY */
            if (sscanf(line,
                    "%*s %*s %Lf %*s %lf %*s %lf %*s %lf %*s %d %*s %lf %*s %lf %*s %lf %d %*s %lf %d",
                    &phypo->probmax, &phypo->misfit,
                    &phypo->grid_misfit_max,
                    &phypo->rms, &phypo->nreadings, &phypo->gap,
                    &phypo->dist,
                    &phypo->amp_mag, &phypo->num_amp_mag,
                    &phypo->dur_mag, &phypo->num_dur_mag
                    ) == EOF)
                goto eof_exit;
            if (GeometryMode == MODE_GLOBAL)
                phypo->dist /= KM2DEG; // always store in memory in km
        } else if (strncmp(line, "QUALITY2", 8) == 0) { // 20100519 AJL Added new hypo line
            /* QUALITY2 */
            if (sscanf(line,
                    "%*s %*s %lf",
                    &phypo->gap_secondary
                    ) == EOF)
                goto eof_exit;
        } else if (strncmp(line, "VPVSRATIO", 9) == 0) {
            /* VPVSRATIO */
            if (sscanf(line,
                    "%*s %*s %lf %*s %d %*s %lf",
                    &phypo->VpVs, &phypo->nVpVs, &phypo->tsp_min_max_diff
                    ) == EOF)
                goto eof_exit;
        } else if (strncmp(line, "STATISTICS", 10) == 0) {
            /* STATISTICS */
            if (sscanf(line,
                    "%*s %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf",
                    &phypo->expect.x, &phypo->expect.y,
                    &phypo->expect.z,
                    &phypo->cov.xx, &phypo->cov.xy, &phypo->cov.xz,
                    &phypo->cov.yy, &phypo->cov.yz, &phypo->cov.zz,
                    &phypo->ellipsoid.az1, &phypo->ellipsoid.dip1,
                    &phypo->ellipsoid.len1,
                    &phypo->ellipsoid.az2, &phypo->ellipsoid.dip2,
                    &phypo->ellipsoid.len2,
                    &phypo->ellipsoid.len3) == EOF)
                goto eof_exit;
            phypo->cov.yx = phypo->cov.xy;
            phypo->cov.zx = phypo->cov.xz;
            phypo->cov.zy = phypo->cov.yz;
        } else if (strncmp(line, "STAT_GEOG", 9) == 0) {
            /* STATISTICS */
            if (sscanf(line,
                    "%*s %*s %lf %*s %lf %*s %lf",
                    &phypo->expect_dlat, &phypo->expect_dlong, &phypo->expect.z) == EOF)
                goto eof_exit;
        } else if (strncmp(line, "MAXIMUM_LIKELIHOOD", 18) == 0) {
            /* MAXIMUM_LIKELIHOOD */
            // 20170811 AJL - added to allow saving of expectation hypocenter results instead of maximum likelihood
            if (sscanf(line,
                    "%*s %*s %lf %*s %lf %*s %lf %*s %lf",
                    &phypo->max_like_dlat, &phypo->max_like_dlong, &phypo->max_like.z, &phypo->max_like_sec) == EOF)
                goto eof_exit;
        } else if (strncmp(line, "ELLIPSOID", 9) == 0) {
            /* ELLIPSOID */
            if (sscanf(line,
                    "%*s %*s %lf %lf %lf  %*s %lf %lf %lf  %*s %lf %lf %lf  %*s %lf",
                    &templat, &templong,
                    &phypo->expect.z,
                    &phypo->ellipsoid.az1, &phypo->ellipsoid.dip1,
                    &phypo->ellipsoid.len1,
                    &phypo->ellipsoid.az2, &phypo->ellipsoid.dip2,
                    &phypo->ellipsoid.len2,
                    &phypo->ellipsoid.len3
                    ) == EOF)
                goto eof_exit;
            latlon2rect(0, templat, templong,
                    &phypo->expect.x, &phypo->expect.y);
        } else if (strncmp(line, "FOCALMECH", 9) == 0) {
            /* FOCALMECH */
            /* Hyp dlat dlong depth Mech dipDir dipAng rake mf misfit nObs nObs */
            if (sscanf(line,
                    "%*s %*s %lf %lf %lf %*s %lf %lf %lf %*s %lf %*s %d",
                    &phypo->focMech.dlat, &phypo->focMech.dlong,
                    &phypo->focMech.depth,
                    &phypo->focMech.dipDir, &phypo->focMech.dipAng,
                    &phypo->focMech.rake,
                    &phypo->focMech.misfit, &phypo->focMech.nObs
                    ) == EOF)
                goto eof_exit;
        } else if (strncmp(line, "TRANSFORM", 9) == 0) {
            /* TRANSFORM */
            if (n_proj >= 0)
                line[MAXLINE_LONG - 1] = '\0';
            strcpy(MapProjStr[n_proj], line);
        } else if (strncmp(line, "DIFFERENTIAL", 12) == 0) {

            /* DIFFERENTIAL */
            /* Nhyp n */
            if (sscanf(line, "%*s %*s %ld", &phypo->event_id) == EOF)
                goto eof_exit;


        } else if (strncmp(line, "QML_OriginQuality", 12) == 0) { // 20181017 AJL - added

            if (sscanf(line,
                    "QML_OriginQuality  assocPhCt %d  usedPhCt %*d  assocStaCt %d  usedStaCt %d  depthPhCt %d"
                    "  stdErr %*lg  azGap %*lg  secAzGap %lg  gtLevel %s  minDist %lg maxDist %lg medDist %lg",
                    &phypo->associatedPhaseCount, &phypo->associatedStationCount, &phypo->usedStationCount, &phypo->depthPhaseCount,
                    &phypo->gap_secondary, phypo->groundTruthLevel, &phypo->minimumDistance, &phypo->maximumDistance, &phypo->medianDistance) == EOF)
                goto eof_exit;

        } else if (strncmp(line, "QML_OriginUncertainty", 12) == 0) { // 20181017 AJL - added

            double azMaxHorUnc;
            if (sscanf(line,
                    "QML_OriginUncertainty  horUnc %*lg  minHorUnc %lg  maxHorUnc %lg  azMaxHorUnc %lg",
                    // 20100617 AJL - horizontalUncertainty: not clear what this is, ignore
                    &phypo->ellipse.len1, &phypo->ellipse.len2, &azMaxHorUnc) == EOF)
                goto eof_exit;
            phypo->ellipse.az1 = azMaxHorUnc - 90.0;
            if (phypo->ellipse.az1 <= 0.0)
                phypo->ellipse.az1 += 360.0;

        } else if (strncmp(line, "PHASE", 5) == 0) {

            /* read arrival/phase parameters */
            /* !!! Remember to modify this sscanf when modifying
                    fprintf format  in function WriteLocation() */

            if (pnarrivals != NULL)
                *pnarrivals = 0;

            /* if requested to read arrivals */
            if (iReadArrivals) {

                /* read arrivals */
                while ((pstr = fgets(line, MAXLINE_LONG, fpio))
                        != NULL
                        && strncmp("END_PHASE", line, 9)) {
                    if (*pnarrivals >= MAX_NUM_ARRIVALS) {
                        sprintf(MsgStr, "WARNING: maximum number of arrivals (%d) exceeded.", MAX_NUM_ARRIVALS);
                        nll_puterr(MsgStr);
                        (*pnarrivals)--;
                        break;
                    }
                    parr = parrivals + *pnarrivals;
                    istat = ReadArrival(line, parr, IO_ARRIVAL_ALL);
                    (*pnarrivals)++;
                };
            }/* not requested to read arrivals */
            else {
                while ((pstr = fgets(line, MAXLINE_LONG, fpio))
                        != NULL
                        && strncmp("END_PHASE", line, 9))
                    ;
            }

            if (pstr == NULL)
                goto eof_exit;
        } else if (strncmp(line, "SCATTER", 7) == 0) {

            /* skip scatter points */
            while ((pstr = fgets(line, MAXLINE_LONG, fpio))
                    != NULL
                    && strncmp("END_SCATTER", line, 11))
                ;
            if (pstr == NULL)
                goto eof_exit;
        } else if (strncmp(line, "QML_", 4) == 0) {
            /* QML_ */
            ;
        } else {

            nll_putmsg(1,
                    "WARNING: unrecognized line in NLLOC hypocenter description:");
            sprintf(MsgStr, "   <%s>", line);
            nll_putmsg(1, MsgStr);
        }

    } /* while(1) */



    /* normal return */
    if (ifile) {
        fclose(fpio);
        NumFilesOpen--;
    }
    return (0);


    /* return at end of file */
eof_exit:
    if (ifile) {

        fclose(fpio);
        NumFilesOpen--;
    }
    return (EOF);

}



/** function to read arrival */

/* returns:	1    if only observation part of phase read
                2    if observation and calculated parts of phase read
                EOF  if EOF or error occurs before any values read
                -1   otherwise
 *
 * NOTE: if this function is changed, then Early-est->trace_processing/timedomain_processing/timedomain_processing_data.c->get_next_pick_nll() should be modified.
 */

int ReadArrival(char* line, ArrivalDesc* parr, int iReadType) {

    int istat, istat2;
    long int idate, ihrmin;
    char *line_calc;
    static char label[10 * ARRIVAL_LABEL_LEN];

    // new values NLL PHASE_2 format
    // 20060629 AJL - Added
    double apriori_weight;
    // 20070326 AJL - Added
    double tt_error;


    /* dummy values for unsupported fields */
    strcpy(parr->network, ARRIVAL_NULL_STR);

    /* read observation part of phase line */

    istat = sscanf(line, "%s %s %s %s %s %s %ld %ld %lf %s %lf %lf %lf %lf",
            label,
            parr->inst,
            parr->comp,
            parr->onset,
            parr->phase,
            parr->first_mot,
            /*&parr->quality, */
            &idate, &ihrmin,
            &(parr->sec),
            parr->error_type,
            // AJL 2000628 bug fix - replaced following line
            //&parr->error,
            &(parr->error),
            &(parr->coda_dur),
            &(parr->amplitude),
            &(parr->period)
            );

    // check for QUAL error type and convert to GAU error using LOCQUAL2ERR
    // 20160727 AJL - added
    if (strcmp(parr->error_type, "QUAL") == 0) {
        parr->quality = (int) lround(parr->error);
        Qual2Err(parr);
    }

    // DEBUG
    /*printf("%s %s %s %s %s %s %ld %ld %lf %s %lf %lf %lf %lf %lf\n",
    parr->label,
    parr->inst,
    parr->comp,
    parr->onset,
    parr->phase,
    parr->first_mot,
    //&parr->quality,
    idate, ihrmin,
    (parr->sec),
    parr->error_type, parr->error,
    (parr->coda_dur),
    (parr->amplitude),
    (parr->period),
    parr->apriori_weight
                  );
     */

    // new values NLL PHASE_2 format
    // 20060629 AJL - Added

    // test input of new values

    istat2 = sscanf(line, "%*s %*s %*s %*s %*s %*s %*d %*d %*f %*s %*f %*f %*f %*f %lf",
            &apriori_weight
            );
    if (istat2 == 1) {
        parr->apriori_weight = apriori_weight;
    } else {
        parr->apriori_weight = 1.0;
    }

    // DEBUG
    /*printf(" %lf\n",
    parr->apriori_weight
    );
     */

    // AJL 20050324 $$$ XXX
    // work around for bug in SeisGram2K Cumul_XX uncertainty
    // parr->error /= 2.0;


    strncpy(parr->label, label, ARRIVAL_LABEL_LEN - 1);


    if (istat == EOF)
        return (istat);
    if (istat != 14)
        return (-1);

    /* decode data and time integers */
    parr->year = idate / 10000;
    idate = idate % 10000;
    parr->month = idate / 100;
    parr->day = idate % 100;
    parr->hour = ihrmin / 100;
    parr->min = ihrmin % 100;

    /* set null quality value */
    parr->quality = -1;


    if (iReadType != IO_ARRIVAL_ALL)
        return (1);

    /* read calculated part of phase line */

    if ((line_calc = strchr(line, '>')) == NULL)
        return (1);

    istat = sscanf(line_calc + 1,
            "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf",
            &parr->pred_travel_time, &parr->residual,
            &parr->weight,
            &parr->station.x, &parr->station.y,
            &parr->station.z,
            &parr->dist, &parr->azim,
            &parr->ray_azim, &parr->ray_dip, &parr->ray_qual,
            &parr->delay
            );
    parr->station.is_coord_xyz = 1; // 20200206 AJL - Bug fix.
    parr->station.is_coord_latlon = 0; // 20200206 AJL - Bug fix.
    if (istat == EOF)
        return (istat);
    //if (istat != 11)  // delay added
    if (istat < 11)
        return (-1);

    // new values NLL PHASE_2 format
    // 20070326 AJL - Added

    // test input of new values

    istat2 = sscanf(line, "%*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*d %*f %lf",
            &tt_error
            );
    if (istat2 == 1) {
        parr->tt_error = tt_error;
    } else {
        parr->tt_error = -1.0;
    }


    if (GeometryMode == MODE_GLOBAL)
        parr->dist /= KM2DEG; // always store in memory in km

    /* convert azimuths to grid coords direction */
    parr->azim = latlon2rectAngle(0, parr->azim);
    parr->ray_azim = latlon2rectAngle(0, parr->ray_azim);

    return (2);

}

/** function to write arrival */

int WriteArrival(FILE* fpio, ArrivalDesc* parr, int iWriteType) {

    int istat;
    long int idate, ihrmin;
    double sta_azim, ray_azim;


    /* code data and time integers */
    idate = parr->year * 10000 + parr->month * 100
            + parr->day;
    ihrmin = parr->hour * 100 + parr->min;


    if (PhaseFormat == FORMAT_PHASE_1) {
        // write observation part of phase line

        istat = fprintf(fpio,
                "%-6.6s %-4.4s %-4.4s %-1.1s %-6.6s %-1.1s %8.8ld %4.4ld %9.4lf %-3.3s %9.2le %9.2le %9.2le %9.2le",
                parr->label,
                parr->inst,
                parr->comp,
                parr->onset,
                parr->phase,
                parr->first_mot,
                /*parr->quality, */
                idate, ihrmin,
                parr->sec,
                parr->error_type, parr->error,
                parr->coda_dur,
                parr->amplitude,
                parr->period
                );
        if (istat < 0)
            return (-1);

    } else if (PhaseFormat == FORMAT_PHASE_2) {
        // write observation part of FORMAT_PHASE_2 phase line
        istat = fprintf(fpio,
                // 20110107 AJL "%-6s %-4s %-4s %-1s %-6s %-1s %8.8ld %4.4ld %9.4lf %-3s %9.2le %9.2le %9.2le %9.2le %9.4lf",
                "%-12s %-4s %-4s %-1s %-6s %-1s %8.8ld %4.4ld %9.4lf %-3s %9.2le %9.2le %9.2le %9.2le %9.4lf",
                parr->label,
                parr->inst,
                parr->comp,
                parr->onset,
                parr->phase,
                parr->first_mot,
                /*parr->quality, */
                idate, ihrmin,
                parr->sec,
                parr->error_type, parr->error,
                parr->coda_dur,
                parr->amplitude,
                parr->period,
                parr->apriori_weight
                );
        if (istat < 0)
            return (-1);
    }



    /* write calculated part of phase line */

    if (iWriteType == IO_ARRIVAL_ALL) {

        /* convert ray azimuth to geographic direction */

        sta_azim = rect2latlonAngle(0, parr->azim);
        ray_azim = rect2latlonAngle(0, parr->ray_azim);

        if (PhaseFormat == FORMAT_PHASE_1) {

            // write calculated part of phase line

            istat = fprintf(fpio,
                    " > %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %6.2lf %5.1lf %5.1lf %2d  %9.4lf",
                    parr->pred_travel_time, parr->residual,
                    parr->weight,
                    parr->station.x, parr->station.y,
                    parr->station.z,
                    GeometryMode == MODE_GLOBAL ? parr->dist * KM2DEG : parr->dist,
                    sta_azim,
                    ray_azim, parr->ray_dip, parr->ray_qual,
                    parr->delay
                    );
            if (istat < 0)
                return (-1);

        } else if (PhaseFormat == FORMAT_PHASE_2) {

            // write calculated part of FORMAT_PHASE_2 phase line
            istat = fprintf(fpio,
                    " > %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %6.2lf %5.1lf %5.1lf %2d  %9.4lf %9.4lf",
                    parr->pred_travel_time, parr->residual,
                    parr->weight,
                    parr->station.x, parr->station.y,
                    parr->station.z,
                    GeometryMode == MODE_GLOBAL ? parr->dist * KM2DEG : parr->dist,
                    sta_azim,
                    ray_azim, parr->ray_dip, parr->ray_qual,
                    parr->delay,
                    parr->tt_error
                    );
            if (istat < 0)
                return (-1);
        }

    }

    istat = fprintf(fpio, "\n");
    if (istat < 0)

        return (-1);

    return (0);

}

/** function to clone (deep copy) a HypoDesc */

HypoDesc * cloneHypoDesc(HypoDesc * phypo_orig) {

    HypoDesc *phypo = NULL;

    // clone HypoDesc

    phypo = (HypoDesc*) malloc(sizeof (HypoDesc));
    if (phypo == NULL) { // memory allocation problem
        nll_puterr("ERROR: cloneHypoDesc(): allocating memory for hypocenter.\n");
    } else {

        *phypo = *phypo_orig;
    }

    /*
    // clone member structures
    //Vect3D expect;
    phypo->expect = (Vect3D*) malloc(sizeof(Vect3D));
     *(phypo->expect) = *(phypo_orig->expect);
    //Mtrx3D cov
    phypo->cov = (Vect3D*) malloc(sizeof(Mtrx3D));
     *(phypo->cov) = *(phypo_orig->cov);
    //Ellipsoid3D ellipsoid;
    phypo->ellipsoid = (Ellipsoid3D*) malloc(sizeof(Ellipsoid3D));
     *(phypo->ellipsoid) = *(phypo_orig->ellipsoid);
    //FocalMech focMech;
    phypo->focMech = (Vect3D*) malloc(sizeof(FocalMech));
     *(phypo->focMech) = *(phypo_orig->focMech);
     */

    return (phypo);

}

/** function to clone (deep copy) an array of ArrivalDesc's */

ArrivalDesc * cloneArrivalDescArray(ArrivalDesc* parrivals_orig, int narrivals) {

    int i;
    ArrivalDesc* parrivals;

    parrivals = (ArrivalDesc*) calloc(narrivals, sizeof (ArrivalDesc));
    if (parrivals == NULL) { // memory allocation problem
        nll_puterr("ERROR: cloneArrivalDescArray(): allocating memory for arrivals.\n");
    } else {

        for (i = 0; i < narrivals; i++)
            parrivals[i] = parrivals_orig[i];
    }

    return (parrivals);

}

/** function to clone (deep copy) a GridDesc */

GridDesc * cloneGridDesc(GridDesc * pgrid_orig) {

    GridDesc* pgrid;

    pgrid = (GridDesc*) malloc(sizeof (GridDesc));
    if (pgrid == NULL) { // memory allocation problem
        nll_puterr("ERROR: cloneGridDesc(): allocating memory for grid.\n");
    } else {

        *pgrid = *pgrid_orig;
    }

    return (pgrid);

}

/** function to write arrival in hypo71/hypoellipse format */

int WriteArrivalHypo(FILE* fpio, ArrivalDesc* arrival, int iwriteEOL) {

    int istat = 0;
    int pha_qual;


    /* write phase */

    pha_qual = (arrival->quality >= 0 && arrival->quality <= 4) ?
            arrival->quality : Err2Qual(arrival);
    if (pha_qual < 0)
        pha_qual = 0; /* !! not necessarily a good choice */


    if (iwriteEOL)
        istat = fprintf(fpio, "\n");

    /* P phase */
    if (strcmp(arrival->phase, "P") == 0) {

        /* write P arrival output */
        istat = fprintf(fpio, "%4.4s", arrival->label);
        istat = fprintf(fpio, "%1s", arrival->onset);
        istat = fprintf(fpio, "%1s", arrival->phase);
        istat = fprintf(fpio, "%1s", arrival->first_mot);
        istat = fprintf(fpio, "%1.1d", pha_qual);
        istat = fprintf(fpio, " %2.2d", arrival->year % 100);
        istat = fprintf(fpio, "%2.2d", arrival->month);
        istat = fprintf(fpio, "%2.2d", arrival->day);
        istat = fprintf(fpio, "%2.2d", arrival->hour);
        istat = fprintf(fpio, "%2.2d", arrival->min);
        istat = fprintf(fpio, "%5.2f", arrival->sec);

        /* S phase */
    } else if (strcmp(arrival->phase, "S") == 0) {

        /* write S phase output (NOTE: assumes written directly after corresponding P) */
        istat = fprintf(fpio, "       %5.2f", arrival->sec);
        istat = fprintf(fpio, " %1s ", arrival->phase);
        istat = fprintf(fpio, "%1.1d", pha_qual);

    }

    if (istat < 0)

        return (-1);

    return (0);

}

/** function to generate GMT JVAL from map transformation parameters */

double getGMTJVAL(int n_proj, char* jval_string, double xlen, double vxmax, double vxmin,
        double ylen, double vymax, double vymin) {

    double gmt_scale;

    jval_string[0] = '\0';

    if (map_itype[n_proj] == MAP_TRANS_GLOBAL) {

        /* -Jmscale or -JMwidth (Mercator [C])
           Give scale along equator (1:xxxx or inch/degree).
         */

        gmt_scale = xlen / (vxmax - vxmin);

        sprintf(jval_string, "-Jm%lf", gmt_scale);

        return (gmt_scale);

    } else if (map_itype[n_proj] == MAP_TRANS_SIMPLE || map_itype[n_proj] == MAP_TRANS_SDC || map_itype[n_proj] == MAP_TRANS_NONE) {

        /* -Jmscale or -JMwidth (Mercator [C])
           Give scale along equator (1:xxxx or inch/degree).
         */

        gmt_scale = xlen / (vxmax - vxmin);

        sprintf(jval_string, "-Jm%lf", gmt_scale);

        return (gmt_scale);

    } else if (map_itype[n_proj] == MAP_TRANS_LAMBERT) {

        /* -JLlon0/lat0/lat1/lat2/width (Lambert [C])
           Give origin, 2 standard parallels, and scale along
           these (1:xxxx or inch/degree)
         */

        gmt_scale = (ylen / (vymax - vymin)); // * (xmaxrect0 - xminrect0) / (xmaxrect - xminrect);
        sprintf(jval_string, "-JL%lf/%lf/%lf/%lf/%lf",
                map_orig_long[n_proj], map_orig_lat[n_proj],
                map_lambert_1st_std_paral[n_proj], map_lambert_2nd_std_paral[n_proj],
                xlen);

        return (gmt_scale);

    } else if (map_itype[n_proj] == MAP_TRANS_TM) {

        /* -Jtlon0/[lat0/]scale or -JTlon0/[lat0/]width (Transverse Mercator [C])
           Give the central meridian lon0, central parallel lat0 (optional), and scale (1:xxxx or UNIT/degree).
         */

        gmt_scale = (ylen / (vymax - vymin)); // * (xmaxrect0 - xminrect0) / (xmaxrect - xminrect);
        sprintf(jval_string, "-JT%lf/%lf/%lf",
                map_orig_long[n_proj], map_orig_lat[n_proj], xlen);

        return (gmt_scale);

    } else if (map_itype[n_proj] == MAP_TRANS_AZ_EQUID) {

        /* -Je|E<lon0>/<lat0>[/<horizon>]/<scale>|<width> (Azimuthal Equidistant)
             <lon0>/<lat0> is the center of the projection.
             <horizon> is max distance from center of the projection (<= 180, default 180).
             <scale> can also be given as <radius>/<lat>, where <radius> is the distance
             in cm to the oblique parallel <lat>.
         */

        gmt_scale = (ylen / (vymax - vymin)); // * (xmaxrect0 - xminrect0) / (xmaxrect - xminrect);
        sprintf(jval_string, "-JE%lf/%lf/180/%lf",
                map_orig_long[n_proj], map_orig_lat[n_proj], xlen);

        return (gmt_scale);

    }

    return (-1.0);

}

/** function to convert between coord systems */

int convertCoordsRect(int proj_index_from, int proj_index_to, double x, double y, double *pxnew, double *pynew) {

    double dlat, dlong;

    if (proj_index_from < 0 || proj_index_to < 0)
        return (-1);

    if (proj_index_from == proj_index_to) {
        *pxnew = x;
        *pynew = y;
        return (0);
    }

    rect2latlon(proj_index_from, x, y, &dlat, &dlong);
    latlon2rect(proj_index_to, dlat, dlong, pxnew, pynew);

    return (0);

}



/** function to convert lat/long to rectangular km coord */

/* rotation about rect coord origin */

int latlon2rect(int n_proj, double dlat, double dlong, double* pxrect, double* pyrect) {

    double xtemp, ytemp;

    double xlt1;


    if (map_itype[n_proj] == MAP_TRANS_GLOBAL) {
        *pxrect = dlong;
        *pyrect = dlat;
        return (0);

    } else if (map_itype[n_proj] == MAP_TRANS_NONE) {
        *pxrect = dlong;
        *pyrect = dlat;
        return (0);

    } else if (map_itype[n_proj] == MAP_TRANS_SIMPLE) {
        xtemp = dlong - map_orig_long[n_proj];
        if (xtemp > 180.0)
            xtemp -= 360.0;
        if (xtemp < -180.0)
            xtemp += 360.0;
        xtemp = xtemp * c111 * cos(cRPD * dlat);
        ytemp = (dlat - map_orig_lat[n_proj]) * c111;
        *pxrect = xtemp * map_cosang[n_proj] - ytemp * map_sinang[n_proj];
        *pyrect = ytemp * map_cosang[n_proj] + xtemp * map_sinang[n_proj];
        return (0);

    } else if (map_itype[n_proj] == MAP_TRANS_SDC) {

        /*
        c  now convert lat and lon differences to km
        x=xln-pln
        c     x=pln-xln
        y=plt-xlt
        xlt1=atan(drlt*tan(drad*(plt+xlt)/120.))
        x=x*xlnkm*cos(xlt1)
        y=y*xltkm
         */

        xtemp = dlong - map_orig_long[n_proj];
        if (xtemp > 180.0)
            xtemp -= 360.0;
        if (xtemp < -180.0)
            xtemp += 360.0;
        ytemp = dlat - map_orig_lat[n_proj];

        xlt1 = atan(MAP_TRANS_SDC_DRLT * tan(DE2RA * (dlat + map_orig_lat[n_proj]) / 2.0));
        xtemp = xtemp * map_sdc_xlnkm[n_proj] * cos(xlt1);
        ytemp = ytemp * map_sdc_xltkm[n_proj];

        *pxrect = xtemp * map_cosang[n_proj] - ytemp * map_sinang[n_proj];
        *pyrect = ytemp * map_cosang[n_proj] + xtemp * map_sinang[n_proj];
        return (0);

    } else if (map_itype[n_proj] == MAP_TRANS_LAMBERT) {

        lamb(n_proj, dlong, dlat, &xtemp, &ytemp);
        xtemp /= 1000.0; /* m -> km */
        ytemp /= 1000.0; /* m -> km */
        *pxrect = xtemp * map_cosang[n_proj] - ytemp * map_sinang[n_proj];
        *pyrect = ytemp * map_cosang[n_proj] + xtemp * map_sinang[n_proj];

        return (0);

    } else if (map_itype[n_proj] == MAP_TRANS_TM) {

        tm(n_proj, dlong, dlat, &xtemp, &ytemp);
        xtemp /= 1000.0; /* m -> km */
        ytemp /= 1000.0; /* m -> km */
        *pxrect = xtemp * map_cosang[n_proj] - ytemp * map_sinang[n_proj];
        *pyrect = ytemp * map_cosang[n_proj] + xtemp * map_sinang[n_proj];

        return (0);

    } else if (map_itype[n_proj] == MAP_TRANS_AZ_EQUID) {

        azeqdist(n_proj, dlong, dlat, &xtemp, &ytemp);
        xtemp /= 1000.0; /* m -> km */
        ytemp /= 1000.0; /* m -> km */
        *pxrect = xtemp * map_cosang[n_proj] - ytemp * map_sinang[n_proj];
        *pyrect = ytemp * map_cosang[n_proj] + xtemp * map_sinang[n_proj];

        return (0);

    }

    return (-1);

}

/** function to convert rectangular km coord to lat/long */

int rect2latlon(int n_proj, double xrect, double yrect, double* pdlat, double* pdlong) {

    double xtemp, ytemp;

    double xlt1;

    if (map_itype[n_proj] == MAP_TRANS_GLOBAL) {
        *pdlat = yrect;
        *pdlong = xrect;
        if (*pdlong < -180.0)
            *pdlong += 360.0;
        else if (*pdlong > 180.0)
            *pdlong -= 360.0;

        return (0);

    } else if (map_itype[n_proj] == MAP_TRANS_NONE) {
        *pdlat = yrect;
        *pdlong = xrect;
        /* 20170321 AJL - bug fix
        if (*pdlong < -180.0)
         *pdlong += 360.0;
        else if (*pdlong > 180.0)
         *pdlong -= 360.0;*/

        return (0);

    } else if (map_itype[n_proj] == MAP_TRANS_SIMPLE) {
        xtemp = xrect * map_cosang[n_proj] + yrect * map_sinang[n_proj];
        ytemp = yrect * map_cosang[n_proj] - xrect * map_sinang[n_proj];
        *pdlat = map_orig_lat[n_proj] + ytemp / c111;
        *pdlong = map_orig_long[n_proj] + xtemp / (c111 * cos(cRPD * *pdlat));
        // 20121005 AJL - prevent longitude outside of -180 -> 180 deg range
        if (*pdlong < -180.0)
            *pdlong += 360.0;
        else if (*pdlong > 180.0)
            *pdlong -= 360.0;

        return (0);

    } else if (map_itype[n_proj] == MAP_TRANS_SDC) {
        xtemp = xrect * map_cosang[n_proj] + yrect * map_sinang[n_proj];
        ytemp = yrect * map_cosang[n_proj] - xrect * map_sinang[n_proj];

        /*
        c
        fy=fy/xltkm
        plt=xlt+fy
        c
        xlt1=atan(rlt*tan(rad*(plt+xlt)/120.))
        fx=fx/(xlnkm*cos(xlt1))
        pln=xln+fx
         */

        ytemp = ytemp / map_sdc_xltkm[n_proj];
        *pdlat = map_orig_lat[n_proj] + ytemp;
        xlt1 = atan(MAP_TRANS_SDC_DRLT * tan(DE2RA * (*pdlat + map_orig_lat[n_proj]) / 2.0));
        xtemp = xtemp / (map_sdc_xlnkm[n_proj] * cos(xlt1));
        *pdlong = map_orig_long[n_proj] + xtemp;
        // 20121005 AJL - prevent longitude outside of -180 -> 180 deg range
        if (*pdlong < -180.0)
            *pdlong += 360.0;
        else if (*pdlong > 180.0)
            *pdlong -= 360.0;

        return (0);

    } else if (map_itype[n_proj] == MAP_TRANS_LAMBERT) {

        xtemp = xrect * map_cosang[n_proj] + yrect * map_sinang[n_proj];
        ytemp = yrect * map_cosang[n_proj] - xrect * map_sinang[n_proj];
        ilamb(n_proj, pdlong, pdlat, xtemp * 1000.0, ytemp * 1000.0);
        // 20121005 AJL - prevent longitude outside of -180 -> 180 deg range
        if (*pdlong < -180.0)
            *pdlong += 360.0;
        else if (*pdlong > 180.0)
            *pdlong -= 360.0;

        return (0);

    } else if (map_itype[n_proj] == MAP_TRANS_TM) {

        xtemp = xrect * map_cosang[n_proj] + yrect * map_sinang[n_proj];
        ytemp = yrect * map_cosang[n_proj] - xrect * map_sinang[n_proj];
        itm(n_proj, pdlong, pdlat, xtemp * 1000.0, ytemp * 1000.0);
        // 20121005 AJL - prevent longitude outside of -180 -> 180 deg range
        if (*pdlong < -180.0)
            *pdlong += 360.0;
        else if (*pdlong > 180.0)
            *pdlong -= 360.0;

        return (0);

    } else if (map_itype[n_proj] == MAP_TRANS_AZ_EQUID) {

        xtemp = xrect * map_cosang[n_proj] + yrect * map_sinang[n_proj];
        ytemp = yrect * map_cosang[n_proj] - xrect * map_sinang[n_proj];
        iazeqdist(n_proj, pdlong, pdlat, xtemp * 1000.0, ytemp * 1000.0);
        if (*pdlong < -180.0)
            *pdlong += 360.0;
        else if (*pdlong > 180.0)
            *pdlong -= 360.0;

        return (0);

    }

    return (-1);
}

/** function to convert rectangular km angle to lat/long angle */

double rect2latlonAngle(int n_proj, double rectAngle) {
    double angle;

    if (map_itype[n_proj] == MAP_TRANS_SIMPLE ||
            map_itype[n_proj] == MAP_TRANS_SDC ||
            map_itype[n_proj] == MAP_TRANS_LAMBERT ||
            map_itype[n_proj] == MAP_TRANS_TM ||
            map_itype[n_proj] == MAP_TRANS_AZ_EQUID) {
        angle = rectAngle - map_rot[n_proj];
        if (angle < 0.0)
            angle += 360.0;
        else if (angle > 360.0)
            angle -= 360.0;
        return (angle);
    } else
        return (rectAngle);
}

/** function to convert lat/long km angle to rectangular angle */

double latlon2rectAngle(int n_proj, double latlonAngle) {
    double angle;

    if (map_itype[n_proj] == MAP_TRANS_SIMPLE ||
            map_itype[n_proj] == MAP_TRANS_SDC ||
            map_itype[n_proj] == MAP_TRANS_LAMBERT ||
            map_itype[n_proj] == MAP_TRANS_TM ||
            map_itype[n_proj] == MAP_TRANS_AZ_EQUID) {
        angle = latlonAngle + map_rot[n_proj];
        if (angle < 0.0)
            angle += 360.0;
        else if (angle > 360.0)
            angle -= 360.0;
        return (angle);
    } else
        return (latlonAngle);
}

/** function to convert source location parameters between coord systems */

int ConvertSourceLoc(int n_proj, SourceDesc *source, int numSources, int toXY, int toLatLon) {

    int istat = 0, nsource;
    SourceDesc *srce_in;


    for (nsource = 0; nsource < numSources; nsource++) {

        srce_in = source + nsource;

        istat = ConvertASourceLocation(n_proj, srce_in, toXY, toLatLon);

    }

    return (istat);

}

/** function to convert source location parameters between coord systems */

int ConvertASourceLocation(int n_proj, SourceDesc *srce_in, int toXY, int toLatLon) {

    int istat = 0;

    if (toXY && srce_in->is_coord_latlon && !srce_in->is_coord_xyz) {
        istat = latlon2rect(n_proj, srce_in->dlat, srce_in->dlong,
                &(srce_in->x), &(srce_in->y));
        srce_in->is_coord_xyz = 1; // 20200206 AJL - Bug fix.
        srce_in->z = srce_in->depth;
    }
    if (toLatLon && srce_in->is_coord_xyz && !srce_in->is_coord_latlon) {

        istat = rect2latlon(n_proj, srce_in->x, srce_in->y,
                &(srce_in->dlat), &(srce_in->dlong));
        srce_in->is_coord_latlon = 1; // 20200206 AJL - Bug fix.
        srce_in->depth = srce_in->z;
    }

    return (istat);

}

/** function to set model coordinates mode to rect or latlon */

int SetModelCoordsMode(int num_surfaces) {
    double xloc, yloc;

    /* if surfaces read, assume lat/long */
    if (num_surfaces > 0) {
        ModelCoordsMode = COORDS_LATLON;
        /* check geographic transformation */
        if (rect2latlon(0, 0.0, 0.0, &yloc, &xloc) < 0) {
            nll_puterr(
                    "FATAL ERROR: geographic transformation required with SURFACE options,\n\tbut transformation (TRANS) not initialized.");
            exit(-1);
        }
    } else
        ModelCoordsMode = COORDS_RECT;

    return (0);

}




/*** date functions */

/** function to convert character month to integer month */

int Month2Int(char* cmonth) {
    int i;

    for (i = 0; i < strlen(cmonth); i++)
        cmonth[i] = toupper(cmonth[i]);

    if (strcmp(cmonth, "JAN") == 0)
        return (1);
    if (strcmp(cmonth, "FEB") == 0)
        return (2);
    if (strcmp(cmonth, "MAR") == 0)
        return (3);
    if (strcmp(cmonth, "APR") == 0)
        return (4);
    if (strcmp(cmonth, "MAY") == 0)
        return (5);
    if (strcmp(cmonth, "JUN") == 0)
        return (6);
    if (strcmp(cmonth, "JUL") == 0)
        return (7);
    if (strcmp(cmonth, "AUG") == 0)
        return (8);
    if (strcmp(cmonth, "SEP") == 0)
        return (9);
    if (strcmp(cmonth, "OCT") == 0)
        return (10);
    if (strcmp(cmonth, "NOV") == 0)
        return (11);
    if (strcmp(cmonth, "DEC") == 0)
        return (12);

    nll_puterr2("ERROR: unrecognized charcter month", cmonth);

    return (0);

}


/*** date functions */

static char daytab[2][13] = {
    {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
    {0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}
};

/** function to set day of year */

int DayOfYear(int year, int month, int day) {
    int i, leap;

    leap = (year % 4 == 0 && year % 100 != 0) || year % 400 == 0;
    for (i = 1; i < month; i++)
        day += daytab[leap][i];

    return day;

}

/** function to set month / day from day of year */

void MonthDay(int year, int yearday, int* pmonth, int* pday) {
    int i, leap;

    leap = (year % 4 == 0 && year % 100 != 0) || year % 400 == 0;

    for (i = 1; yearday > daytab[leap][i]; i++)
        yearday -= daytab[leap][i];
    *pmonth = i;
    *pday = yearday;

}

/** function to construct current date/time string */

char* CurrTimeStr(void) {

    static char timestr[MAXLINE];
    time_t curr_time;

    curr_time = time(NULL);

    strftime(timestr, (size_t) MAXLINE, "%d%b%Y %Hh%Mm%S", localtime(&curr_time));

    return (timestr);

}

/** function to check for and expand wild card characters in filenames
        and to return a list of equivalent files */

int ExpandWildCards(char* fileFilter, char fileList[][FILENAME_MAX], int maxNumFiles) {

    int nfiles = 0;
    char *pchr;


    /* check for no '*' or '?' character */

    if ((pchr = strchr(fileFilter, '*')) == NULL && (pchr = strchr(fileFilter, '?')) == NULL) {
        strcpy(fileList[0], fileFilter);
        nfiles = 1;
        return (nfiles);
    }


    // get directory and filename
    char directory[FILENAME_MAX];
    if ((pchr = strrchr(fileFilter, '/')) != NULL) {
        strncpy(directory, fileFilter, pchr - fileFilter);
        directory[pchr - fileFilter] = '\0';
        strcpy(ExpandWildCards_pattern, pchr + 1);
    } else {
        strcpy(directory, ".");
        strcpy(ExpandWildCards_pattern, fileFilter);
    }
    /*DEBUG*///printf("directory: <%s>\n", directory);
    /*DEBUG*///printf("filename: <%s>\n", ExpandWildCards_pattern);


    /* expand wildcard file names into list of files */

    nfiles = 0;
    struct dirent **namelist;
    int n;
    n = scandir(directory, &namelist, fnmatch_wrapper, alphasort);
    //n = scandir(directory, &namelist, NULL, alphasort);
    /*DEBUG*///printf("directory: <%s>  n %d\n", directory, n);
    if (n <= 0) {
        nll_puterr2("ERROR: expanding wildcard filenames in: ", fileFilter);
        return (-1);
    } else if (n == 0) {
        nll_puterr2("ERROR: empty directory: expanding wildcard filenames in: ", fileFilter);
        return (-1);
    } else if (n > maxNumFiles) { // 20111011 AJL - added this block to catch excess number of wildcard files
        sprintf(MsgStr,
                "ERROR: too many files: expanding wildcard filenames in: %s, max number of files = %d",
                fileFilter, maxNumFiles);
        nll_puterr(MsgStr);
        return (-1);
    } else {
        while (--n >= 0) {

            sprintf(fileList[n], "%s/%s", directory, namelist[n]->d_name);
            /*DEBUG*///printf("%s -> %s\n", namelist[n]->d_name,fileList[n]);
            nfiles++;
            free(namelist[n]);
            namelist[n] = NULL;
        }
        free(namelist);
        namelist = NULL;
    }

    return (nfiles);


}

/** function to wrap fnmatch */

int fnmatch_wrapper(const struct dirent * entry) {
    int match = 0;
    int flags = 0;
    match = fnmatch(ExpandWildCards_pattern, entry->d_name, flags);

    /*DEBUG*///printf("ExpandWildCards_pattern: <%s>  entry->d_name: <%s>  match %d\n", ExpandWildCards_pattern, entry->d_name, match);

    return (!match);

}

/** function to check for and expand wild card characters in filenames
        and to return a list of equivalent files */

int ExpandWildCards_OLD(char* fileName, char fileList[][FILENAME_MAX], int maxNumFiles) {
    int istat;
    int nfiles = 0;
    char system_str[MAXLINE];
    char list_file[FILENAME_MAX] = "filelist.tmp";
    char *pchr;
    FILE* fpio;


    /* check for no '*' or '?' character */

    if ((pchr = strchr(fileName, '*')) == NULL
            && (pchr = strchr(fileName, '?')) == NULL) {
        strcpy(fileList[0], fileName);
        nfiles = 1;
        return (nfiles);
    }


    /* expand wildcard file names into list of files */

    sprintf(system_str, "ls %s > %s", fileName, list_file);
    system(system_str);

    if ((fpio = fopen(list_file, "r")) == NULL) {
        nll_puterr2("ERROR: opening fileList temporary file: ", list_file);
        return (-1);
    }
    NumFilesOpen++;

    nfiles = 0;
    while (nfiles < maxNumFiles &&
            (istat = fscanf(fpio, "%s", fileList[nfiles])) != EOF
            && istat == 1) {
        nfiles++;
    }


    fclose(fpio);
    NumFilesOpen--;

    return (nfiles);


}

/** function to sort arrivals by obs_time field */

int SortArrivalsTime(ArrivalDesc* arrival, int num_arrivals) {

    qsort((void *) arrival, (size_t) num_arrivals, sizeof (ArrivalDesc),
            (int (*)(const void *, const void *)) CmpArrivalsTime);

    return (0);

}

/** function to compare arrivals by obs_time field */

int CmpArrivalsTime(const ArrivalDesc *keyval, const ArrivalDesc * datum) {
    if (keyval->obs_time < datum->obs_time)
        return (-1);
    if (keyval->obs_time > datum->obs_time)

        return (1);

    return (0);
}

/** function to sort arrivals by flag_ignore field */

int SortArrivalsIgnore(ArrivalDesc* arrival, int num_arrivals) {

    qsort((void *) arrival, (size_t) num_arrivals, sizeof (ArrivalDesc),
            (int (*)(const void *, const void *)) CmpArrivalsIgnore);

    return (0);

}

/** function to compare arrivals by flag_ignore field */

int CmpArrivalsIgnore(const ArrivalDesc *keyval, const ArrivalDesc * datum) {
    if (keyval->flag_ignore < datum->flag_ignore)
        return (-1);
    if (keyval->flag_ignore > datum->flag_ignore)

        return (1);

    return (0);
}

/** function to sort arrivals by dist field */

int SortArrivalsDist(ArrivalDesc* arrival, int num_arrivals) {

    qsort((void *) arrival, (size_t) num_arrivals, sizeof (ArrivalDesc),
            (int (*)(const void *, const void *)) CmpArrivalsDist);

    return (0);

}

/** function to compare arrivals by dist field */

int CmpArrivalsDist(const ArrivalDesc *keyval, const ArrivalDesc * datum) {
    if (keyval->dist < datum->dist)
        return (-1);
    if (keyval->dist > datum->dist)
        return (1);

    // if same distance, assume same sta sort by time   // 20161004 AJL - added
    if (keyval->obs_time < datum->obs_time)
        return (-1);
    if (keyval->obs_time > datum->obs_time)

        return (1);

    return (0);
}

/** function to sort array of doubles */

int SortDoubles(double* array, int num_elements) {

    qsort((void *) array, (size_t) num_elements, sizeof (double),
            (int (*)(const void *, const void *)) CmpDoubles);

    return (0);

}

/** function to compare doubles */

int CmpDoubles(const double *keyval, const double *datum) {
    if (*keyval < *datum)
        return (-1);
    if (*keyval > *datum)

        return (1);

    return (0);
}

/** function to set angle values in take-off angles union */

TakeOffAngles SetTakeOffAngles(double azim, double dip, int iqual) {
    TakeOffAngles angles;

    // 20201221 AJL - Bug fix: cannot represent a negative number as unsigned short - duh!
    if (azim < 0.0) {
        angles.ival[1] = ANGLES_DIP_REVERSE;
    } else {
        angles.ival[1] = (unsigned short) (0.5 + 10.0 * azim);
    }
    angles.ival[0] = (unsigned short) iqual
            + (unsigned short) ANGLES_OFFSET
            * (unsigned short) (0.5 + 10.0 * dip);

    return (angles);
}

/** function to set float values in take-off angles union */

void SetAnglesFloat(TakeOffAngles* pangles, float fvalue) {

    pangles->fval = fvalue;
}

/** function to get values in take-off angles union */

int GetTakeOffAngles(TakeOffAngles *pangles, double *pazim, double *pdip, int *piqual) {
    *pazim = ((double) pangles->ival[1]) / 10.0;
    *pdip = ((double) (pangles->ival[0] / (int) ANGLES_OFFSET)) / 10.0;
    *piqual = (int) pangles->ival[0] % (int) ANGLES_OFFSET;

    //printf("i0 %d  i1 %d\n", pangles->ival[0], pangles->ival[1]);

    return (*piqual);
}

/** function to read take-off angles from file */

int ReadTakeOffAnglesFile(char *fname, double xloc, double yloc, double zloc,
        double *pazim, double *pdip, int *piqual, double sta_azim, int iSwapBytes) {
    int istat;
    FILE *fp_grid, *fp_hdr;
    float fvalue;
    GridDesc gdesc;
    TakeOffAngles angles;

    //printf("DEBUG: ReadTakeOffAnglesFile: fname %s\n", fname);

    /* open angle grid file */
    if ((istat = OpenGrid3dFile(fname, &fp_grid, &fp_hdr, &gdesc, "angle", NULL, iSwapBytes)) < 0) {
        if (message_flag >= 3) {
            sprintf(MsgStr, "WARNING: cannot open angle grid file, ignoring angles: %s", fname);
            nll_putmsg(3, MsgStr);
            //printf("DEBUG: ReadTakeOffAnglesFile: WARNING: cannot open angle grid file, ignoring angles: %s", fname);
        }
        angles = SetTakeOffAngles(0.0, 0.0, 0);
        GetTakeOffAngles(&angles, pazim, pdip, piqual);
        return (-1);
    }

    /* get angles float value on grid */
    fvalue = ReadAbsInterpGrid3d(fp_grid, &gdesc, xloc, yloc, zloc, 0);

    /* get angles */
    SetAnglesFloat(&angles, fvalue);
    GetTakeOffAngles(&angles, pazim, pdip, piqual);
    //printf("DEBUG: ReadTakeOffAnglesFile: x y z %f %f %f  raz %f  rdip %f rq  %d\n", xloc, yloc, zloc, *pazim, *pdip, *piqual);

    /* determine azimuth (2D grids) */
    if (gdesc.type == GRID_ANGLE_2D) {
        //printf("DEBUG: *pazim %f <= ANGLES_DIP_MAX %f\n", *pazim, ANGLES_DIP_MAX);
        // 20201221 AJL - Bug fix: cannot represent a negative number as unsigned short - duh!
        //if (*pazim > 0.0)
        if (*pazim <= ANGLES_DIP_MAX)
            *pazim = sta_azim;
        else { // reverse azimuth
            *pazim = sta_azim - 180.0;
            if (*pazim < 0.0)
                *pazim += 360.0;
        }
    }

    /* close angle grid file */
    CloseGrid3dFile(&gdesc, &fp_grid, &fp_hdr);

    return (0);


}

/** function to generate normally distributed deviate with zero mean and unit variance */

double normal_dist_deviate() {

    static int iset = 0;
    static float gset;
    double fac, r, v1, v2;

    if (iset == 0) {
        do {
            v1 = get_rand_double(-1.0, 1.0);
            v2 = get_rand_double(-1.0, 1.0);
            r = v1 * v1 + v2*v2;
        } while (r >= 1.0);
        fac = sqrt(-2.0 * log(r) / r);
        gset = v1 * fac;
        iset = 1;
        return v2 * fac;
    } else {
        iset = 0;

        return gset;
    }
}


/** function to test normal_dist_deviate function */

#define NUM_BIN 21
#define WIDTH 3.0

void test_normal_dist_deviate() {
    long nmax = 210000;
    long n, m;
    long ibin[NUM_BIN];
    double test;
    double dbin, binmax[NUM_BIN];
    int sum_one_std;

    dbin = 2.0 / (float) NUM_BIN;
    for (n = 0; n < NUM_BIN; n++) {
        ibin[n] = 0;
        binmax[n] = WIDTH * ((double) (n + 1) * dbin - 1.0);
    }


    for (n = 0; n < nmax; n++) {
        test = normal_dist_deviate();
        m = 0;
        while (test > binmax[m] && m < NUM_BIN - 1)
            m++;
        ibin[m]++;
    }

    sum_one_std = 0;
    fprintf(stdout,
            "\nnormal_dist_deviate function test (samples= %ld)\n", nmax);
    fprintf(stdout, "  Bin -Inf,%lf  N=%ld\n", binmax[0], ibin[0]);
    for (n = 1; n < NUM_BIN - 1; n++) {
        fprintf(stdout,
                "  Bin %lf,%lf  N=%ld\n", binmax[n - 1], binmax[n], ibin[n]);

        if (binmax[n - 1] >= -1.0 && binmax[n] <= 1.0)
            sum_one_std += ibin[n];
    }
    fprintf(stdout,
            "  Bin %lf,Inf  N=%ld\n", binmax[NUM_BIN - 2],
            ibin[NUM_BIN - 1]);
    fprintf(stdout,
            "Percent in range (-1,1) %lf\n", (double) sum_one_std / (double) nmax);

}

/** function to generate and display "traditional" statistics from grid file */

int GenTraditionStats(GridDesc *pgrid, Vect3D *pexpect, Mtrx3D *pcov,
        FILE * fpgrid) {
    int istat;


    /* allocate grid buffer */

    pgrid->buffer = AllocateGrid(pgrid);
    if (pgrid->buffer == NULL) {
        nll_puterr(
                "ERROR: allocating memory for 3D PDF grid buffer.");
        exit(EXIT_ERROR_MEMORY);
    }

    /* create grid array access pointers */

    pgrid->array = CreateGridArray(pgrid);
    if (pgrid->array == NULL) {
        nll_puterr(
                "ERROR: creating array for accessing 3D PDF grid buffer.");
        exit(EXIT_ERROR_MEMORY);
    }


    /* read PDF grid */

    if ((istat = ReadGrid3dBuf(pgrid, fpgrid)) < 0) {
        nll_puterr("ERROR: reading PDF grid from disk.");
        exit(EXIT_ERROR_IO);
    }


    /* calculate expectation */

    *pexpect = CalcExpectation(pgrid, NULL);
    if (message_flag >= 3) {
        sprintf(MsgStr, "EXPECTATION { x %lf  y %lf  z %lf }",
                pexpect->x, pexpect->y, pexpect->z);
        nll_putmsg(3, MsgStr);
    }

    /* calculate covariance matrix */

    *pcov = CalcCovariance(pgrid, pexpect, NULL);
    if (message_flag >= 3) {
        sprintf(MsgStr, "COVARIANCE: {");
        nll_putmsg(3, MsgStr);
        sprintf(MsgStr, "   xx: %lf  xy: %lf  xz: %lf",
                pcov->xx, pcov->xy, pcov->xz);
        nll_putmsg(3, MsgStr);
        sprintf(MsgStr, "   yx: %lf  yy: %lf  yz: %lf",
                pcov->yx, pcov->yy, pcov->yz);
        nll_putmsg(3, MsgStr);
        sprintf(MsgStr, "   zx: %lf  zy: %lf  zz: %lf",
                pcov->zx, pcov->zy, pcov->zz);
        nll_putmsg(3, MsgStr);
        sprintf(MsgStr, "}");
        nll_putmsg(3, MsgStr);
    }

    /* clean up */

    FreeGrid(pgrid);
    DestroyGridArray(pgrid);

    return (0);

}

/** function to calculate the expectation (mean) of a PDF grid */

Vect3D CalcExpectation(GridDesc* pgrid, FILE * fpgrid) {

    int ix, iy, iz;

    GRID_FLOAT_TYPE val;
    double volume;
    Vect3D expect = {0.0, 0.0, 0.0};


    /* cannot calculate for misfit grid */
    if (pgrid->type == GRID_MISFIT) {
        expect.x = expect.y = expect.z = -LARGE_DOUBLE;
        return (expect);
    }


    for (ix = 0; ix < pgrid->numx; ix++) {
        for (iy = 0; iy < pgrid->numy; iy++) {
            for (iz = 0; iz < pgrid->numz; iz++) {

                if (fpgrid != NULL)
                    val = ReadGrid3dValue(fpgrid,
                        ix, iy, iz, pgrid, 0);
                else
                    val = ((GRID_FLOAT_TYPE ***) pgrid->array)[ix][iy][iz];

                expect.x += (double) val * (double) ix;
                expect.y += (double) val * (double) iy;
                expect.z += (double) val * (double) iz;

            }
        }
    }

    volume = pgrid->dx * pgrid->dy * pgrid->dz;
    expect.x = pgrid->origx + expect.x * pgrid->dx * volume;
    expect.y = pgrid->origy + expect.y * pgrid->dy * volume;
    expect.z = pgrid->origz + expect.z * pgrid->dz * volume;

    return (expect);
}

/** function to calculate the covariance a PDF grid */

Mtrx3D CalcCovariance_OLD(GridDesc* pgrid, Vect3D* pexpect, FILE * fpgrid) {

    int ix, iy, iz;

    double val;
    double x, y, z, xx, xy, xz, yy, yz, zz;
    double volume;

    Mtrx3D cov = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};



    /* cannot calculate for misfit grid */
    if (pgrid->type == GRID_MISFIT) {
        cov.xx = cov.xy = cov.xz =
                cov.yx = cov.yy = cov.yz =
                cov.zx = cov.zy = cov.zz
                = -LARGE_DOUBLE;
        return (cov);
    }



    /* calculate covariance following eq. (6-12), T & V, 1982 */


    for (ix = 0; ix < pgrid->numx; ix++) {
        x = pgrid->origx + (double) ix * pgrid->dx;
        xx = x * x;

        for (iy = 0; iy < pgrid->numy; iy++) {
            y = pgrid->origy + (double) iy * pgrid->dy;
            yy = y * y;
            xy = x * y;

            for (iz = 0; iz < pgrid->numz; iz++) {
                z = pgrid->origz + (double) iz * pgrid->dz;
                xz = x * z;
                yz = y * z;
                zz = z * z;

                if (fpgrid != NULL)
                    val = ReadGrid3dValue(fpgrid, ix, iy, iz, pgrid, 0);
                else
                    val = ((GRID_FLOAT_TYPE ***) pgrid->array)[ix][iy][iz];

                if (val < 0.0) {
                    printf("ERROR: CalcCovariance: Grid value < 0: ixyz= %d %d %d  value= %g\n", ix, iy, iz, val);
                    continue;
                }

                cov.xx += (double) val * xx;
                cov.xy += (double) val * xy;
                cov.xz += (double) val * xz;

                cov.yy += (double) val * yy;
                cov.yz += (double) val * yz;

                cov.zz += (double) val * zz;

            }
        }
    }

    volume = pgrid->dx * pgrid->dy * pgrid->dz;

    //printf("DEBUG: cov.yy = cov.yy(%g) * volume(%g) (= %g) - pexpect->y(%g) * pexpect->y (= %g)\n", cov.yy, volume, cov.yy * volume, pexpect->y, pexpect->y * pexpect->y);
    cov.xx = cov.xx * volume - pexpect->x * pexpect->x;
    cov.xy = cov.xy * volume - pexpect->x * pexpect->y;
    cov.xz = cov.xz * volume - pexpect->x * pexpect->z;

    cov.yx = cov.xy;
    cov.yy = cov.yy * volume - pexpect->y * pexpect->y;
    cov.yz = cov.yz * volume - pexpect->y * pexpect->z;

    cov.zx = cov.xz;
    cov.zy = cov.yz;
    cov.zz = cov.zz * volume - pexpect->z * pexpect->z;

    //printf("DEBUG: CalcCovariance: volume= %g  cov.yy= %g\n", volume, cov.yy);

    return (cov);
}

/** function to calculate the covariance a PDF grid */

Mtrx3D CalcCovariance(GridDesc* pgrid, Vect3D* pexpect, FILE * fpgrid) {

    int ix, iy, iz;

    double val;
    double x, y, z, xx, xy, xz, yy, yz, zz;
    double volume;

    Mtrx3D cov = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};



    /* cannot calculate for misfit grid */
    if (pgrid->type == GRID_MISFIT) {
        cov.xx = cov.xy = cov.xz =
                cov.yx = cov.yy = cov.yz =
                cov.zx = cov.zy = cov.zz
                = -LARGE_DOUBLE;
        return (cov);
    }



    /* calculate covariance following eq. (6-12), T & V, 1982 */


    for (ix = 0; ix < pgrid->numx; ix++) {
        x = pgrid->origx + (double) ix * pgrid->dx - pexpect->x;
        xx = x * x;

        for (iy = 0; iy < pgrid->numy; iy++) {
            y = pgrid->origy + (double) iy * pgrid->dy - pexpect->y;
            yy = y * y;
            xy = x * y;

            for (iz = 0; iz < pgrid->numz; iz++) {
                z = pgrid->origz + (double) iz * pgrid->dz - pexpect->z;
                xz = x * z;
                yz = y * z;
                zz = z * z;

                if (fpgrid != NULL)
                    val = ReadGrid3dValue(fpgrid, ix, iy, iz, pgrid, 0);
                else
                    val = ((GRID_FLOAT_TYPE ***) pgrid->array)[ix][iy][iz];

                if (val < 0.0) {
                    printf("ERROR: CalcCovariance: Grid value < 0: ixyz= %d %d %d  value= %g\n", ix, iy, iz, val);
                    continue;
                }

                cov.xx += (double) val * xx;
                cov.xy += (double) val * xy;
                cov.xz += (double) val * xz;

                cov.yy += (double) val * yy;
                cov.yz += (double) val * yz;

                cov.zz += (double) val * zz;

            }
        }
    }

    volume = pgrid->dx * pgrid->dy * pgrid->dz;

    //printf("DEBUG: cov.yy = cov.yy(%g) * volume(%g) (= %g) - pexpect->y(%g) * pexpect->y (= %g)\n", cov.yy, volume, cov.yy * volume, pexpect->y, pexpect->y * pexpect->y);
    cov.xx = cov.xx * volume;
    cov.xy = cov.xy * volume;
    cov.xz = cov.xz * volume;

    cov.yx = cov.xy;
    cov.yy = cov.yy * volume;
    cov.yz = cov.yz * volume;

    cov.zx = cov.xz;
    cov.zy = cov.yz;
    cov.zz = cov.zz * volume;

    //printf("DEBUG: CalcCovariance: volume= %g  cov.yy= %g\n", volume, cov.yy);

    return (cov);
}

/** function to calculate the covariance of a set of samples */

Mtrx3D CalcCovarianceSamples(float* fdata, int nSamples, Vect3D * pexpect) {

    //printf("GeometryMode %s\n", GeometryMode == MODE_GLOBAL ? "MODE_GLOBAL" : "MODE_RECT");
    if (GeometryMode == MODE_GLOBAL)
        return (CalcCovarianceSamplesGlobal(fdata, nSamples, pexpect));

    else
        return (CalcCovarianceSamplesRect(fdata, nSamples, pexpect));

}

/** function to read expectation and covariance from a hyp file */

int ReadHypStatistics(FILE **pfpio, char* fnroot_in,
        Vect3D* pmax_like, Vect3D* pexpect,
        Mtrx3D* pcov, Ellipsoid3D *pellipsoid,
        ArrivalDesc* parrivals, int *pnarrivals) {

    char fn_in[FILENAME_MAX];
    static HypoDesc hypo;


    /* open hypocenter file if necessary */

    if (*pfpio == NULL) {
        sprintf(fn_in, "%s.hyp", fnroot_in);
        if ((*pfpio = fopen(fn_in, "r")) == NULL) {
            nll_puterr("ERROR: opening hypocenter file.");
            return (EOF);
        }
        NumFilesOpen++;
    }


    /* read next hypocenter */

    if (GetHypLoc(*pfpio, fnroot_in, &hypo, parrivals, pnarrivals, 1, NULL, 0) != EOF) {
        pmax_like->x = hypo.x;
        pmax_like->y = hypo.y;
        pmax_like->z = hypo.z;
        *pexpect = hypo.expect;
        *pcov = hypo.cov;
        *pellipsoid = hypo.ellipsoid;
        return (0);
    }


    /* end of file */

    fclose(*pfpio);
    NumFilesOpen--;

    return (EOF);

}

/** function to read mechanism from a hyp file */

int ReadFocalMech(FILE **pfpio, char* fnroot_in,
        FocalMech* pfocalMech,
        ArrivalDesc* parrivals, int *pnarrivals) {

    char fn_in[FILENAME_MAX];
    static HypoDesc hypo;


    /* open hypocenter file if necessary */

    if (*pfpio == NULL) {
        sprintf(fn_in, "%s.hyp", fnroot_in);
        if ((*pfpio = fopen(fn_in, "r")) == NULL) {
            nll_puterr("ERROR: opening hypocenter file.");
            return (EOF);
        }
        NumFilesOpen++;
    }


    /* read next hypocenter */

    if (GetHypLoc(*pfpio, fnroot_in, &hypo, parrivals, pnarrivals, 1, NULL, 0) != EOF) {
        *pfocalMech = hypo.focMech;
        return (0);
    }


    /* end of file */

    fclose(*pfpio);
    NumFilesOpen--;

    return (EOF);

}

/** function to read arrivals with first-motions readings from a hyp file
 *
 * returns ArrivalDesc array containing only arrivals with first motion readings
 * readings are converted to "+" or "-"
 *
 */

int ReadFirstMotionArrivals(FILE **pfpio, char* fnroot_in, ArrivalDesc* parrivals, int *pnarrivals) {

    char fn_in[FILENAME_MAX];
    static HypoDesc hypo;

    // open hypocenter file if necessary

    if (*pfpio == NULL) {
        sprintf(fn_in, "%s.hyp", fnroot_in);
        if ((*pfpio = fopen(fn_in, "r")) == NULL) {
            nll_puterr("ERROR: opening hypocenter file.");
            return (EOF);
        }
        NumFilesOpen++;
    }

    // read next hypocenter

    if (GetHypLoc(*pfpio, fnroot_in, &hypo, parrivals, pnarrivals, 1, NULL, 0) != EOF) {

        // compact arrival array to contain only those with first-motion readings
        int naccept = 0;
        ArrivalDesc* parr = parrivals;
        for (int narr = 0; narr < *pnarrivals; narr++) {
            if (strstr("CcUu+", parr->first_mot)) { // follows fmamp conventions in fmamp/read_input.c
                strcpy(parr->first_mot, "+");
                parrivals[naccept] = *parr;
                naccept++;
            } else if (strstr("DdRr-", parr->first_mot)) { // follows fmamp conventions in fmamp/read_input.c
                strcpy(parr->first_mot, "-");
                parrivals[naccept] = *parr;
                naccept++;
            }
            parr++;
        }
        *pnarrivals = naccept;

        return (0);
    }


    // end of file

    fclose(*pfpio);
    NumFilesOpen--;

    return (EOF);

}

/** function to trim leading and trailing blanks from character strings */

int TrimString(char* line) {
    char *line0, *line1, *pendchr;

    /* find end of string */
    pendchr = strchr(line, '\0');
    if (pendchr == NULL)
        return (-1);

    /* remove leading spaces */
    while (isspace(line[0])) {
        line0 = line + 1;
        line1 = line;
        do
            *(line1++) = *line0; while (*(line0++) != '\0');
    }

    /* re-find end of string */
    pendchr = strchr(line, '\0');
    if (pendchr == NULL)
        return (-1);

    /* remove trailing spaces */
    // AJL 20071018 while (--pendchr > line && isspace(*pendchr))
    while (--pendchr > line && (isspace(*pendchr) || *pendchr == '\t' || *pendchr == '\n' || *pendchr == '\r'))
        *pendchr = '\0';

    return (pendchr - line);

}

/** function to check for blank line */

int LineIsBlank(char *line) {
    char* pstr;

    pstr = line;
    while (*pstr) {
        if (!isspace(*pstr++))

            return (0);
    }

    return (1);
}

/** function to read FORTRAN character format fields */

int ReadFortranString(char* line, int istart, int ilen, char* string_in) {
    char chrtmp[MAXLINE], *chrpos;


    /* check for null char before end of read */
    chrpos = strchr(line, '\0');
    if (chrpos - line < istart + ilen - 1)
        return (-1);

    /* copy field to temp string */
    strncpy(chrtmp, line + (istart - 1), ilen);
    chrtmp[ilen] = '\0';

    /* check for blank field */
    /*	for (n = 0; n < ilen; n++) {
                    if (chrtmp[n] != ' ')
                            break;
            }
            if (n == ilen) {
                    for (n = 0; n < ilen; n++)
                            string_in[n] = ' ';
                    return(1);
            }
     */

    strncpy(string_in, chrtmp, ilen);
    string_in[ilen] = '\0';

    return (1);

}

/** function to read FORTRAN integer format fields */

int ReadFortranInt(char* line, int istart, int ilen, int* pintval) {
    char chrtmp[MAXLINE], *chrpos;
    int istat, n;


    /* check for null char before end of read */
    chrpos = strchr(line, '\0');
    if (chrpos - line < istart + ilen - 1)
        return (-1);

    /* copy field to temp string */
    strncpy(chrtmp, line + (istart - 1), ilen);
    chrtmp[ilen] = '\0';

    /* check for blank field */
    for (n = 0; n < ilen; n++) {
        if (chrtmp[n] != ' ')
            break;
    }
    if (n == ilen) {
        *pintval = 0;
        return (1);
    }


    istat = sscanf(chrtmp, "%d", pintval);

    return (istat);

}

/** function to read FORTRAN real format fields */

int ReadFortranReal(char* line, int istart, int ilen, double* pdblval) {
    char chrtmp[MAXLINE], *chrpos;
    int istat, n;


    /* check for null char before end of read */
    chrpos = strchr(line, '\0');
    if (chrpos - line < istart + ilen - 1)
        return (-1);

    /* copy field to temp string */
    strncpy(chrtmp, line + (istart - 1), ilen);
    chrtmp[ilen] = '\0';

    /* check for blank field */
    for (n = 0; n < ilen; n++) {
        if (chrtmp[n] != ' ')
            break;
    }
    if (n == ilen) {
        *pdblval = 0.0;
        return (1);
    }


    istat = sscanf(chrtmp, "%lf", pdblval);

    return (istat);

}

/** function to convert phase error to hypo71 style quality */

int Err2Qual(ArrivalDesc * arrival) {
    int errLevel;

    /* find error level with error >= arrival error */
    for (errLevel = 0; errLevel < NumQuality2ErrorLevels; errLevel++) {
        if (arrival->error <= Quality2Error[errLevel]) {
            arrival->quality = errLevel;

            return (errLevel);
        }
    }

    /* not found or no quality/error values available */
    return (-1);
}

/** function to convert quality to error */

void Qual2Err(ArrivalDesc * arrival) {

    /* set error fields */
    strcpy(arrival->error_type, "GAU");
    if (arrival->quality >= 0 &&
            arrival->quality < NumQuality2ErrorLevels) {
        arrival->error = Quality2Error[arrival->quality];
    } else {
        arrival->error = Quality2Error[NumQuality2ErrorLevels - 1];
        nll_puterr("WARNING: invalid arrival quality.");
    }

}

/** function to read Quality2Error values ***/

int GetQuality2Err(char* line1) {

    int istat, ierr, nlev;
    double qual2err;
    char frmt[MAXLINE] = "%lf";
    char frmttmp[MAXLINE];


    while ((istat = sscanf(line1, frmt, &qual2err)) == 1) {
        Quality2Error[NumQuality2ErrorLevels++] = qual2err;
        sprintf(frmttmp, "%%*f %s", frmt);
        strcpy(frmt, frmttmp);
    }


    if (message_flag >= 2) {
        sprintf(MsgStr, "NLLoc LOCQUAL2ERR:");
        nll_putmsg(2, MsgStr);
    }
    ierr = 0;
    for (nlev = 0; nlev < NumQuality2ErrorLevels; nlev++) {
        if (message_flag >= 2) {
            sprintf(MsgStr, " %d ->  %lf", nlev, Quality2Error[nlev]);
            nll_putmsg(2, MsgStr);
        }
        if ((istat = checkRangeDouble("QUAL2ERR", "Quality2Error",
                Quality2Error[nlev], 1, 0.0, 0, 0.0)) != 0)
            ierr = -1;
    }

    if (ierr < 0)

        return (-1);

    return (0);
}

/** function to check if phase_in matches phase_check
true if phase_in == phase_check
true if phase_in is in phase ID list for phase_check */

int IsPhaseID(char *phase_in, char *phase_check) {
    int npha_id;
    char test_str[PHASE_LABEL_LEN + 4];

    //printf("phase_in %s, phase_check %s, ", phase_in, phase_check);
    /* check for blank phase_in */
    if (strstr("              ", phase_in) != NULL) {
        //printf(" NO - BLANK\n");
        return (0);
    }

    // check if phase_in == phase_check (AJL 20041201 added)
    if (strcmp(phase_in, phase_check) == 0) {
        //printf(" YES (%s %s)\n", phase_in, phase_check);
        return (1);
    }

    // check if phase_in is in phase ID list for phase_check
    removeSpace(phase_in);
    sprintf(test_str, " %s ", phase_in);
    for (npha_id = 0; npha_id < NumPhaseID; npha_id++) {
        //printf(" compare (<%s> <%s>)\n", PhaseID[npha_id].id_string, test_str);
        if (strcmp(PhaseID[npha_id].phase, phase_check) == 0) {
            if (strstr(PhaseID[npha_id].id_string, test_str) != NULL) {
                //printf(" YES (%s %s)\n", PhaseID[npha_id].id_string,  phase_in);

                return (1);
            }
        }
    }

    //printf(" NO\n");
    return (0);
}

/** function to evaluate phase ID from phase ID lists */
// 20161004 AJL - moved here from NLLocLib.c

int EvalPhaseID(char *phase_out, char *phase_in) {
    int npha_id;

    for (npha_id = 0; npha_id < NumPhaseID; npha_id++) {
        if (IsPhaseID(phase_in, PhaseID[npha_id].phase)) {
            if (phase_out != NULL)
                strcpy(phase_out, PhaseID[npha_id].phase);
            return (1);
        }
    }

    if (phase_out != NULL)
        strcpy(phase_out, phase_in);

    return (0); // returned unchanged
}

/** function remove whitespace from string */

void removeSpace(char *str) {

    int n, m;

    n = 0;
    while (str[n] != '\0' && n < 1000000)
        if (isspace(str[n])) {
            m = n;
            while (str[m] != '\0' && m < 1000000) {
                str[m] = str[m + 1];
                m++;
            }
        } else {

            n++;
        }

}

/** function to read fpfit summary record to HypoDesc structure */

int ReadFpfitSum(FILE *fp_in, HypoDesc * phypo) {

    int istat;
    char *cstat;
    double mag, dtemp;

    double deg, dmin;
    char strNS[2], strMagType[2];

    static char line[MAXLINE_LONG];


    /* read next line */
    cstat = fgets(line, MAXLINE_LONG, fp_in);
    if (cstat == NULL)
        return (EOF);

    /* read hypocenter parameters */

    istat = 0;
    istat += ReadFortranInt(line, 1, 2, &phypo->year);
    if (phypo->year < 20)
        phypo->year += 2000;
    if (phypo->year < 100)
        phypo->year += 1900;
    istat += ReadFortranInt(line, 3, 2, &phypo->month);
    istat += ReadFortranInt(line, 5, 2, &phypo->day);
    istat += ReadFortranInt(line, 8, 2, &phypo->hour);
    istat += ReadFortranInt(line, 10, 2, &phypo->min);
    istat += ReadFortranReal(line, 12, 6, &phypo->sec);

    istat += ReadFortranReal(line, 18, 3, &deg);
    istat += ReadFortranString(line, 21, 1, strNS);
    istat += ReadFortranReal(line, 22, 5, &dmin);
    phypo->dlat = deg + dmin / 60.0;
    if (strncmp(strNS, "S", 1) == 0)
        phypo->dlat = -phypo->dlat;
    istat += ReadFortranReal(line, 27, 4, &deg);
    istat += ReadFortranString(line, 31, 1, strNS);
    istat += ReadFortranReal(line, 32, 5, &dmin);
    phypo->dlong = deg + dmin / 60.0;
    if (strncmp(strNS, "W", 1) == 0)
        phypo->dlong = -phypo->dlong;

    istat += ReadFortranReal(line, 37, 7, &phypo->depth);

    istat += ReadFortranReal(line, 46, 5, &mag);

    istat += ReadFortranInt(line, 51, 3, &phypo->nreadings);
    istat += ReadFortranReal(line, 54, 4, &dtemp);
    phypo->gap = (int) 0.5 + dtemp;
    istat += ReadFortranReal(line, 58, 5, &phypo->dist);
    istat += ReadFortranReal(line, 63, 5, &phypo->rms);

    /* hypoinverse horiz and vertical error are converted to ellipsoid,
                    this is not correct statistically  */
    istat += ReadFortranReal(line, 68, 5, &(phypo->ellipsoid.len1));
    phypo->ellipsoid.az1 = 0.0;
    phypo->ellipsoid.dip1 = 0.0;
    phypo->ellipsoid.len2 = phypo->ellipsoid.len1;
    phypo->ellipsoid.az2 = 90.0;
    phypo->ellipsoid.dip2 = 0.0;
    istat += ReadFortranReal(line, 73, 5, &(phypo->ellipsoid.len3));

    istat += ReadFortranString(line, 80, 1, strMagType);

    /* focal mechanism parameters */

    istat += ReadFortranReal(line, 82, 3, &phypo->focMech.dipDir);
    istat += ReadFortranReal(line, 86, 2, &phypo->focMech.dipAng);
    istat += ReadFortranReal(line, 88, 4, &phypo->focMech.rake);
    istat += ReadFortranReal(line, 94, 4, &phypo->focMech.misfit);
    istat += ReadFortranInt(line, 99, 3, &phypo->focMech.nObs);
    istat += ReadFortranReal(line, 103, 5, &phypo->focMech.misfit90);
    istat += ReadFortranReal(line, 109, 4, &phypo->focMech.staDist);
    istat += ReadFortranReal(line, 114, 4, &phypo->focMech.ratioMH);
    istat += ReadFortranReal(line, 120, 2, &phypo->focMech.conf90strike);
    istat += ReadFortranReal(line, 123, 2, &phypo->focMech.conf90dip);
    istat += ReadFortranReal(line, 126, 2, &phypo->focMech.conf90rake);
    istat += ReadFortranString(line, 128, 1, phypo->focMech.convFlag);
    istat += ReadFortranString(line, 129, 1, phypo->focMech.multSolFlag);

    return (istat);



}













//========================================================
// DD
// see hypoDD doc (Open File Report 01-113)

/** function to read hypoDD Initial hypocenter input format to HypoDesc structure */

/* example
19850124   2195871   37.8832  -122.2415      9.800  1.4    0.15    0.51   0.02      38542
19911126  14274555   37.8738  -122.2432      9.950  1.4    0.22    0.53   0.09     238298
 */

/*
One event per line: DATE, TIME, LAT, LON, DEP, MAG, EH, EV, RMS, ID
Parameter Description:
        DATE Concatenated origin date (YYYYMMDD; year, month and day).
        TIME Concatenated origin time (HHMMSSSS; hour, minute and seconds).
        LAT, LON, DEP Latitude, longitude, depth (km).
        MAG Magnitude (set to 0.0 if not available).
        EH, EV Horizontal and vertical error (km) of original location (set to 0.0 if not available).
        RMS Residual rms (s) (set to 0.0 is not available).
        ID Event identification [I9].
 */


int ReadHypoDDInitHypo(FILE *fp_in, HypoDesc *phypo, int n_proj) {

    int istat;
    char *cstat;
    double err_horiz, err_vert;

    static char line[MAXLINE_LONG];


    /* read next line */
    cstat = fgets(line, MAXLINE_LONG, fp_in);
    if (cstat == NULL)
        return (EOF);
    // DEBUG
    printf("%s\n", line);

    /* read hypocenter parameters */

    istat = 0;
    istat += ReadFortranInt(line, 1, 4, &phypo->year);
    istat += ReadFortranInt(line, 5, 2, &phypo->month);
    istat += ReadFortranInt(line, 7, 2, &phypo->day);
    istat += ReadFortranInt(line, 11, 2, &phypo->hour);
    istat += ReadFortranInt(line, 13, 2, &phypo->min);
    if (line[16] == '.') {
        // non-standard HypoDD: decimal point in seconds
        // YMD HMS str1 str2 lat lon depth mag errH errZ rms event_id
        // 20200414 101739.37 XXX XXX 31.705758 -104.035262 6.996 2.7 -1 -1 -1 1001
        istat += ReadFortranReal(line, 15, 5, &phypo->sec);
    } else {
        istat += ReadFortranReal(line, 15, 4, &phypo->sec);
        phypo->sec /= 100.0;
    }

    istat += sscanf(line, "%*s %*s %lf %lf %lf %lf %lf %lf %lf %ld",
            &phypo->dlat, &phypo->dlong, &phypo->depth, &phypo->amp_mag,
            &err_horiz, &err_vert, &phypo->rms, &phypo->event_id);

    // DEBUG
    printf("%d %d %d %d %d %f %s %s %f %f %f %f %f %f %f %ld\n", phypo->year, phypo->month, phypo->day, phypo->hour, phypo->min, phypo->sec, "X", "X", phypo->dlat, phypo->dlong, phypo->depth, phypo->amp_mag, err_horiz, err_vert, phypo->rms, phypo->event_id);

    //phypo->depth -= 8.0;

    /* hypoDD horiz and vertical error are converted to ellipsoid,
                    this is not correct statistically  */
    phypo->ellipsoid.az1 = 0.0;
    phypo->ellipsoid.dip1 = 0.0;
    phypo->ellipsoid.len1 = phypo->ellipsoid.len2 = err_horiz;
    phypo->ellipsoid.az2 = 90.0;
    phypo->ellipsoid.dip2 = 0.0;
    phypo->ellipsoid.len3 = err_vert;


    /* set additional hypocenter fields */
    latlon2rect(n_proj, phypo->dlat, phypo->dlong, &(phypo->x), &(phypo->y));
    phypo->z = phypo->depth;
    phypo->dotime = 0.0;

    return (istat == 14 ? 1 : -1);


}


/** function to write hypoDD Initial hypocenter input format to out */

/* example
19850124   2195871   37.8832  -122.2415      9.800  1.4    0.15    0.51   0.02      38542  x y z
19911126  14274555   37.8738  -122.2432      9.950  1.4    0.22    0.53   0.09     238298
 */


int WriteHypoDDInitHypo(FILE *fp_out, HypoDesc * phypo) {

    fprintf(fp_out,
            "%4.4d%2.2d%2.2d  %2.2d%2.2d%4.4d %9.4f %10.4f %10.3f %4.1f %7.2f %7.2f %6.2f %10ld  %10.4f %10.4f %10.4f\n",
            phypo->year, phypo->month, phypo->day,
            phypo->hour, phypo->min, (int) (phypo->sec * 100.0 + 0.5),
            phypo->dlat, phypo->dlong, phypo->depth, phypo->amp_mag,
            phypo->ellipsoid.len1, phypo->ellipsoid.len3, phypo->rms, phypo->event_id,
            phypo->x, phypo->y, phypo->z);

    return (0);


}



/** function to write hypoDD Cross correlation differential time input (e.g. file dt.cc) format to out */

/* example:
#    28136    46442     -0.174000
#   402093   402094     -0.009000
NCCMO        0.000000 -1.994625  ??
NCCMO    -0.038245201    0.63    P
NCCRP    -0.062200069    0.78    P
NCJBG    -0.020649910    0.60    P
NCJPS    -0.011000156    0.50    P
 */

int WriteHypoDDXCorrDiff(FILE *fp_out, int num_arrivals, ArrivalDesc *arrival, HypoDesc * hypos) {

    int narr;
    long dd_event_id_1, dd_event_id_2;
    double dd_otime_corr;
    ArrivalDesc *parr;


    dd_event_id_1 = dd_event_id_2 = -1;

    for (narr = 0; narr < num_arrivals; narr++) {

        parr = arrival + narr;

        if (parr->flag_ignore)
            continue;

        if (parr->dd_event_id_1 != dd_event_id_1 || parr->dd_event_id_2 != dd_event_id_2) {

            dd_event_id_1 = parr->dd_event_id_1;
            dd_event_id_2 = parr->dd_event_id_2;
            // set new otime corr OTC based on current perturbatios to otimes
            // original OTC was incorporated into dd_dtime when reading dt file.
            dd_otime_corr = hypos[parr->dd_event_index_1].dotime
                    - hypos[parr->dd_event_index_2].dotime;
            fprintf(fp_out, "# %8ld %8ld %13.6lf\n",
                    dd_event_id_1, dd_event_id_2, dd_otime_corr);
        }

        fprintf(fp_out, "%-8s %12lf %7lf %4s\n",
                parr->label, parr->dd_dtime, parr->weight, parr->phase);

    }


    return (0);

}

/** function to write arrival */

int WriteDiffArrival(FILE* fpio, HypoDesc* hypos, ArrivalDesc* parr, int iWriteType) {

    int istat;
    double sta_azim, ray_azim;
    double dd_otime_corr;


    /* write observation part of phase line */

    // set new otime corr OTC based on current perturbatios to otimes
    // original OTC was incorporated into dd_dtime when reading dt file.
    dd_otime_corr = hypos[parr->dd_event_index_1].dotime
            - hypos[parr->dd_event_index_2].dotime;
    istat = fprintf(fpio,
            "%-6.6s %-4.4s %-4.4s %-6.6s %8ld %8ld %9.5lf %9.5lf",
            parr->label,
            parr->inst,
            parr->comp,
            parr->phase,
            parr->dd_event_id_1, parr->dd_event_id_2,
            parr->dd_dtime - dd_otime_corr,
            parr->weight
            );
    if (istat < 0)
        return (-1);


    if (iWriteType == IO_ARRIVAL_ALL) {

        /* convert ray aziumuth to geographic direction */

        sta_azim = rect2latlonAngle(0, parr->azim);
        ray_azim = rect2latlonAngle(0, parr->ray_azim);

        /* write calculated part of phase line */

        istat = fprintf(fpio,
                " > %9.5lf %9.5lf %9.4lf %9.4lf %9.4lf %9.4lf %6.2lf %5.1lf %5.1lf %2d",
                parr->pred_travel_time, parr->residual,
                parr->station.x, parr->station.y,
                parr->station.z,
                GeometryMode == MODE_GLOBAL ? parr->dist * KM2DEG : parr->dist,
                sta_azim,
                ray_azim, parr->ray_dip, parr->ray_qual
                );

        if (istat < 0)
            return (-1);

    }

    istat = fprintf(fpio, "\n");
    if (istat < 0)

        return (-1);

    return (0);

}

/** function to generate take-off angles from travel time grid
                                using a numerical gradient algorithm */

int CalcAnglesGradient(GridDesc* ptgrid, GridDesc* pagrid, int angle_mode, int grid_mode) {

    int ix, iy, iz, edge_flagx = 0, edge_flagy = 0, iflag2D = 0;
    double origx, origy, origz;
    double dx, dy, dz, dvol;
    double xlow = 0.0, xhigh = 0.0;
    double azim, dip;
    int iqual;

    TakeOffAngles angles = AnglesNULL;


    /* write message */
    sprintf(MsgStr, "Generating take-off angle grid...");
    nll_putmsg(1, MsgStr);


    if (grid_mode == GRID_TIME_2D) {
        iflag2D = 1;
        xlow = xhigh = 0.0;
    }

    /* estimate take-off angles from numerical gradients */

    origx = pagrid->origx;
    origy = pagrid->origy;
    origz = pagrid->origz;
    dx = pagrid->dx;
    dy = pagrid->dy;
    dz = pagrid->dz;
    dvol = dx * dy * dz;

    for (ix = 0; ix < pagrid->numx; ix++) {
        /* 2D grids, store angles in ix = 0 sheet */
        if (ix == 1 && iflag2D)
            edge_flagx = 1;
        if ((ix == 0 || ix == pagrid->numx - 1)
                && grid_mode == GRID_TIME)
            edge_flagx = 1;
        for (iy = 0; iy < pagrid->numy; iy++) {
            if (iy == 0 || iy == pagrid->numy - 1)
                edge_flagy = 1;
            for (iz = 0; iz < pagrid->numz; iz++) {

                /* no calculation for edges of grid */
                if (edge_flagx || edge_flagy
                        || iz == 0 || iz == pagrid->numz - 1) {
                    ((GRID_FLOAT_TYPE ***) pagrid->array)[ix][iy][iz] = AnglesNULL.fval;
                    continue;
                }

                if (!iflag2D) {
                    xlow = ((GRID_FLOAT_TYPE ***) ptgrid->array)[ix - 1][iy][iz];
                    xhigh = ((GRID_FLOAT_TYPE ***) ptgrid->array)[ix + 1][iy][iz];
                }
                angles = GetGradientAngles(
                        ((GRID_FLOAT_TYPE ***) ptgrid->array)[ix][iy][iz],
                        xlow,
                        xhigh,
                        ((GRID_FLOAT_TYPE ***) ptgrid->array)[ix][iy - 1][iz],
                        ((GRID_FLOAT_TYPE ***) ptgrid->array)[ix][iy + 1][iz],
                        /* intentional reversal of z
                                signs to get pos = up */
                        ((GRID_FLOAT_TYPE ***) ptgrid->array)[ix][iy][iz + 1],
                        ((GRID_FLOAT_TYPE ***) ptgrid->array)[ix][iy][iz - 1],
                        dx, dy, dz, iflag2D,
                        &azim, &dip, &iqual);
                if (angle_mode == ANGLE_MODE_YES) {
                    ((GRID_FLOAT_TYPE ***) pagrid->array)[ix][iy][iz] = angles.fval;
                } else if (angle_mode == ANGLE_MODE_INCLINATION) {
                    ((

                            GRID_FLOAT_TYPE ***) pagrid->array)[ix][iy][iz] = dip;
                }

            }
            edge_flagy = 0;
        }
        edge_flagx = 0;
    }


    return (0);

}

/** function to generate take-off angles from time grid node values */

TakeOffAngles GetGradientAngles(double vcent, double xlow, double xhigh,
        double ylow, double yhigh, double zlow, double zhigh,
        double dx, double dy, double dz, int iflag2D,
        double *pazim, double *pdip, int *piqual) {

    double grad_low, grad_high, gradx, grady, gradz, azim, dip;
    int iqualx, iqualy, iqualz, iqual, iflip;
    TakeOffAngles angles = AnglesNULL;



    /* calculate gradient of travel time and quality in Z direction */
    grad_low = (vcent - zlow) / dz;
    grad_high = (zhigh - vcent) / dz;
    iqualz = CalcAnglesQuality(grad_low, grad_high);
    gradz = (grad_low + grad_high) / 2.0;
    gradz = -gradz; /* reverse sign to get take-off angle */

    /* calculate gradient of travel time and quality in Y direction */
    grad_low = (vcent - ylow) / dy;
    grad_high = (yhigh - vcent) / dy;
    iqualy = CalcAnglesQuality(grad_low, grad_high);
    grady = (grad_low + grad_high) / 2.0;
    grady = -grady; /* reverse sign to get take-off angle */

    /* thats all for 2D grids */
    if (iflag2D) {
        /* calculate dip angle (range of 0 (down) to 180 (up)) */
        //dip = atan2(grady, -gradz) / cRPD;
        dip = atan2(fabs(grady), -gradz) / cRPD;
        //dip = -dip; // TEST!
        iflip = 0;
        // 20201221 AJL - Bug fix: ray depart direction depend only on grady
        if (grady > 0.0) { // ray departs source away from station
            iflip = 1;
        }
        // 20201221 AJL - Bug fix: with atan2(<positive>, real) dip is restricted to 0-180 (atan2 restricted to 0-PI)
        /*if (dip > 180.0) {
            dip = dip - 180.0;
            iflip = 1;
        } else if (dip < 0.0) {
            dip = -dip;
            iflip = 1;
        }*/
        /* calculate azimuth polarity (1 or -1) relative to pos Y dir */
        azim = iflip ? -1.0 : 1.0;
        /* find combined quality - weighted average of component qual */
        iqual = (fabs(grady) * (double) iqualy
                + fabs(gradz) * (double) iqualz)
                / (fabs(grady) + fabs(gradz));
        /* set angles */
        angles = SetTakeOffAngles(azim, dip, iqual);
        *pazim = azim;
        *pdip = dip;
        *piqual = iqual;
        return (angles);
    }

    /* calculate gradient of travel time and quality in X direction */
    grad_low = (vcent - xlow) / dx;
    grad_high = (xhigh - vcent) / dx;
    iqualx = CalcAnglesQuality(grad_low, grad_high);
    gradx = (grad_low + grad_high) / 2.0;
    gradx = -gradx; /* reverse sign to get take-off angle */

    /* find combined quality - weighted average of component qual */
    iqual = (fabs(gradx) * (double) iqualx
            + fabs(grady) * (double) iqualy
            + fabs(gradz) * (double) iqualz)
            / (fabs(gradx) + fabs(grady) + fabs(gradz));

    /* calculate dip angle (range of 0 (down) to 180 (up)) */
    dip = atan2(sqrt(gradx * gradx + grady * grady), -gradz) / cRPD;
    /* calculate azimuth angle (0 to 360) */
    azim = atan2(gradx, grady) / cRPD;
    if (azim < 0.0)
        azim += 360.0;
    angles = SetTakeOffAngles(azim, dip, iqual);

    // return double angles values
    *pazim = azim;
    *pdip = dip;
    *piqual = iqual;

    return (angles);

}



/** function to estimate quality of take-off angle determination */

/* quality is:	0 if sign of A = grad_low and B = grad_high differ
                0->10 as (2AB / (AA + BB)) -> 1;
 */

int CalcAnglesQuality(double grad_low, double grad_high) {

    double ratio;

    /* if both gradients are zero, return highest quality */
    if (fabs(grad_low) + fabs(grad_high) < SMALL_DOUBLE)
        return (10);

    /* calculate quality */
    ratio = 2.0 * grad_low * grad_high /
            (grad_low * grad_low + grad_high * grad_high);
    return (ratio > 0.0 ? (int) (10.01 * ratio) : 0);

}




// END - DD
//========================================================



// 20190509 AJL - following function moved here from NLLocLib.c

/** function to generate sample (scatter) of location PDF
 *
 * IMPORTANT: assumes integral of dvol * prob_den over grid = 1.0
 */

int GenEventScatterGrid(GridDesc* ptgrid, HypoDesc* phypo, ScatterParams* pscat, char* filename) {

    FILE *fpio;
    char fname[4 * FILENAME_MAX];

    int ix, iy, iz;
    double origx, origy, origz;
    double dx, dy, dz, dvol;
    double xval, yval, zval;
    float fdata[4];
    int tot_npoints = 0;
    double xnpoints;
    double probmax, prob_den;



    /* return if no scatter samples requested */
    if (pscat->npts < 1)
        return (0);

    sprintf(fname, "%s.scat", filename);

    /* write message */
    if (message_flag >= 3) {
        nll_putmsg(3, "");
        nll_putmsg2(3, "Generating event scatter file:", fname);
    }

    /* open scatter file */

    // 20190509 AJL  sprintf(fname, "%s.loc.scat", filename);
    if ((fpio = fopen(fname, "w")) == NULL) {
        nll_puterr("ERROR: opening scatter output file.");
        return (-1);
    } else {
        NumFilesOpen++;
    }
    /* skip header record (used later to store number of samples taken) */
    fseek(fpio, 4 * sizeof (float), SEEK_SET);

    /* generate N=Scatter->npts events with prob P=prob_den/probmax at  */
    /*	uniformly-randomly chosen grid locations */

    origx = ptgrid->origx;
    origy = ptgrid->origy;
    origz = ptgrid->origz;
    dx = ptgrid->dx;
    dy = ptgrid->dy;
    dz = ptgrid->dz;
    dvol = dx * dy * dz;
    probmax = phypo->probmax;

    for (ix = 0; ix < ptgrid->numx; ix++) {
        for (iy = 0; iy < ptgrid->numy; iy++) {
            for (iz = 0; iz < ptgrid->numz; iz++) {
                //printf("DEBUG: message_flag %d, ix, iy, iz %d %d %d  %d\n", message_flag, ix, iy, iz, tot_npoints);

                prob_den = ((GRID_FLOAT_TYPE ***) ptgrid->array)[ix][iy][iz];

                xnpoints = (double) pscat->npts * dvol * prob_den;

                xval = origx + (double) ix * dx;
                yval = origy + (double) iy * dy;
                zval = origz + (double) iz * dz;

                while (xnpoints > 0.0) {

                    if (xnpoints > 1.0 || xnpoints - (double) ((int) xnpoints) > get_rand_double(0.0, 1.0)) {
                        fdata[0] = (float) (xval + get_rand_double(-dx / 2.0, dx / 2.0));
                        fdata[1] = (float) (yval + get_rand_double(-dy / 2.0, dy / 2.0));
                        fdata[2] = (float) (zval + get_rand_double(-dz / 2.0, dz / 2.0));
                        fdata[3] = (float) prob_den;
                        fwrite(fdata, sizeof (float), 4, fpio);
                        tot_npoints++;
                    }

                    xnpoints -= 1.0;

                }

            }
        }
    }


    /* write header information */
    fseek(fpio, 0, SEEK_SET);
    fwrite(&tot_npoints, sizeof (int), 1, fpio);
    fdata[0] = (float) probmax;
    fwrite(fdata, sizeof (float), 1, fpio);

    fclose(fpio);
    NumFilesOpen--;

    /* write message */
    if (message_flag >= 3) {
        sprintf(MsgStr, "  %d points generated.", tot_npoints);
        nll_putmsg(3, MsgStr);
        sprintf(MsgStr, "  (%d points requested, dvol= %lf, probmax=%lf)",
                pscat->npts, dvol, probmax);
        nll_putmsg(3, MsgStr);
    }

    return (0);

}

/** function to calculate the integral over volume of a PDF grid */

double IntegrateGrid(GridDesc* pgrid, int flag_normalize) {

    double integral = 0.0;

    int ix, iy, iz;
    double dvolume = pgrid->dx * pgrid->dy * pgrid->dz;

    // integrate
    for (ix = 0; ix < pgrid->numx; ix++) {
        for (iy = 0; iy < pgrid->numy; iy++) {
            for (iz = 0; iz < pgrid->numz; iz++) {
                integral += dvolume * ((GRID_FLOAT_TYPE ***) pgrid->array)[ix][iy][iz];
            }
        }
    }
    //printf("DEBUG:  IntegrateGrid: integral=%f\n", integral);

    // normalize
    if (flag_normalize && fabs(integral) > FLT_MIN) {
        for (ix = 0; ix < pgrid->numx; ix++) {
            for (iy = 0; iy < pgrid->numy; iy++) {
                for (iz = 0; iz < pgrid->numz; iz++) {
                    ((GRID_FLOAT_TYPE ***) pgrid->array)[ix][iy][iz] /= integral;
                }
            }
        }
    }

    return (integral);
}



// 20200122 AJL - following function moved here from NLLocLib.c

/** function to add arrival station to station list */

int addToStationList(SourceDesc *stations, int numStations, ArrivalDesc *arrival, int nArrivals, int iuse_phaseid_in_label) {

    int i, n, nAdded = 0;


    // set default ignored flag, weight and nullname
    for (n = 0; n < MAX_NUM_ARRIVALS; n++) {
        (stations + n)->ignored = 1;
        (stations + n)->station_weight = 1.0;
        //strcpy((stations + n)->label, "");
    }

    char arrival_label[ARRIVAL_LABEL_LEN];

    for (i = 0; i < nArrivals; i++) {

        // 20200122 AJL - added
        if (iuse_phaseid_in_label) {
            sprintf(arrival_label, "%s#%s", (arrival + i)->label, (arrival + i)->phase);
        } else {
            strcpy(arrival_label, (arrival + i)->label);
        }
        // find station in list
        n = 0;
        for (; n < numStations; n++) {
            if (strcmp((stations + n)->label, arrival_label) == 0) {
                break; // already in list
            }
        }

        // not found, add station to list
        if (n == numStations) {
            if (numStations >= MAX_NUM_ARRIVALS) {
                sprintf(MsgStr, "ERROR: addToStationList: numStations (%d) >= MAX_NUM_ARRIVALS (%d): cannot add station %s ",
                        numStations, MAX_NUM_ARRIVALS, arrival_label);
                nll_puterr(MsgStr);
                //return (0); // 20101209 AJL - bug fix
                continue;
            }
            *(stations + n) = (arrival + i)->station;
            strcpy((stations + n)->label, arrival_label);
            nAdded++;
            numStations++;
            if (message_flag >= 4) {
                sprintf(MsgStr, "Added to station list: %s (%lf,%lf,%lf)",
                        (stations + n)->label, (stations + n)->x, (stations + n)->y,
                        (stations + n)->z);
                nll_putmsg(4, MsgStr);
            }
        }

        if (!((arrival + i)->flag_ignore))
            (stations + n)->ignored = 0;

    }

    return (numStations);


}

/** function to write station list to file */

int WriteStationList(FILE* fpio, SourceDesc *stations, int numStations) {

    int n;

    for (n = 0; n < numStations; n++) {

        fprintf(fpio, "%s %lf %lf %lf\n",
                (stations + n)->label, (stations + n)->x, (stations + n)->y, (stations + n)->z);
    }

    return (0);

}

/** function to read phase identification values ***/

int GetPhaseID(char* line1) {
    int istat, ilen;
    char *substr, *cpos;


    if (NumPhaseID >= MAX_NUM_PHASE_ID) {
        nll_puterr(
                "LOCPHASEID: WARNING: maximum number of PhaseIDs reached, ignoring phase ID.");
        return (-1);
    }

    if ((istat = sscanf(line1, "%s", PhaseID[NumPhaseID].phase)) != 1)
        return (-1);

    substr = strstr(line1, PhaseID[NumPhaseID].phase);

    /* save phase id values with spaces at each end */
    if ((cpos = strchr(substr, '\n')) != NULL)
        *cpos = '\0';
    sprintf(PhaseID[NumPhaseID].id_string, " %s ", substr + 1);

    if ((ilen = strlen(PhaseID[NumPhaseID].id_string)) == 0)
        return (-1);

    sprintf(MsgStr, "LOCPHASEID:");
    nll_putmsg(3, MsgStr);
    sprintf(MsgStr, "  Phase: %s  PhaseID: <%s>",
            PhaseID[NumPhaseID].phase,
            PhaseID[NumPhaseID].id_string);
    nll_putmsg(3, MsgStr);

    NumPhaseID++;

    return (0);
}

