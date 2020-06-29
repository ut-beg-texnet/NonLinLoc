/*
 * Copyright (C) 1999-2008 Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
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


/*   NLLoc_func_test.c

        Program to demonstrate running NLLoc through a function call.


 */

/*-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
161 Allee du Micocoulier, 06370 Mouans-Sartoux, France
tel: +33(0)493752502  e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/


/*
        history:	(see also http://alomax.net/nlloc -> Updates)

        ver 01    17DEC2007  AJL  Original version
        ver 02    29NOV2010  AJL  Added multiple calls of NLLoc funciton if multiple observation files specified

        see NLLoc1.c and NLLocLib.c


.........1.........2.........3.........4.........5.........6.........7.........8

 */





#define PNAME  "NLLoc_func_test"

#include "GridLib.h"
#include "ran1/ran1.h"
#include "velmod.h"
#include "GridMemLib.h"
#include "calc_crust_corr.h"
#include "phaseloclist.h"
#include "otime_limit.h"
#include "NLLocLib.h"




/** program to demonstrate running NLLoc through a function call
 *
 *  This demonstration program uses a standard NLL control file specified by the program argument <control file>.
 *  The control file is parsed into a set of character strings in memory which are passed to the
 *  function invocation of NLLoc.
 *
 *  Some data needed for location may be on disk files, e.g
 *      station locations (if specified by an INCLUDE statement in the NLL control file),
 *      observations files (if the <obs file> parameter is "-" or not present, then the observations files specified
 *         in the LOCFILES statement in the NLL control file will be read,
 *      travel-time grids as specified in LOCFILES statement in the NLL control file are always
 *         read from disk files.
 *
 *  Most location result information is returned from the function invocation of NLLoc
 *  though pointers in memory, this information is written to disk files by this demonstration program.
 *  Depending on the LOCHYPOUT parameter in the control file, the location results may also be written to disk
 *  within the function invocation of NLLoc. e.g. for no disk output use: LOCHYPOUT NONE
 */


#define NARGS_MIN 2
#define ARG_DESC "<control file> [<obs file>]"

int main(int argc, char *argv[]) {

    int istat = 0;
    char line[4 * MAXLINE];


    char pid_main[255]; // string process id (for CUSTOM_ETH)
    char **param_line_array = NULL;
    char **obs_line_array = NULL;
    LocNode *loc_list_head = NULL; // root node of location list
    int return_locations, return_oct_tree_grid, return_scatter_sample;





    // set program name
    strcpy(prog_name, PNAME);

    // check command line for correct usage
    if (argc < NARGS_MIN) {
        disp_usage(prog_name, ARG_DESC);
        return (EXIT_ERROR_USAGE);
    }

    // set file names
    // control file name
    char fn_control_main[MAXLINE];
    strcpy(fn_control_main, argv[1]);


    /** ===========================================================================
     *  Loop over all observation files specified in command line
     */
    int id_filename = 0;
    int narg;
    for (narg = 1; narg < argc; narg++) {

        if (narg == 1 && argc > 2) // observation file(s) specified on command line
            continue;

        // obs file name
        char fn_obs[MAXLINE];
        if (narg == 1) // no observation files specified on command line
            strcpy(fn_obs, "-");
        else
            strcpy(fn_obs, argv[narg]);
        sprintf(MsgStr, "========> Running observation file: %d: %s", narg - 1, fn_obs);
        nll_putmsg(1, MsgStr);

        int n_param_lines = 0;
        int n_obs_lines = 0;

        /** ===========================================================================
         *  Convert nll control file specified on command line to array of string in memory.
         *  A production program might construct these control strings entirely in memory.
         */

        FILE* fp_control;

        if ((fp_control = fopen(fn_control_main, "r")) == NULL) {
            nll_puterr("FATAL ERROR: opening control file.");
            return (EXIT_ERROR_FILEIO);
        } else {
            NumFilesOpen++;
        }

        param_line_array = (char **) calloc(1000, sizeof (char *));
        while (fp_control != NULL && fgets(line, 4 * MAXLINE, fp_control) != NULL) {
            param_line_array[n_param_lines] = (char *) malloc(4 * MAXLINE);
            strcpy(param_line_array[n_param_lines], line);
            n_param_lines++;
        }
        fclose(fp_control);
        NumFilesOpen--;



        /** ===========================================================================
         *  Convert obs file specified on command line to array of string in memory.
         *  A production program might construct these control strings entirely in memory.
         */

        if (strcmp(fn_obs, "-") != 0) {

            FILE* fp_obs;

            if ((fp_obs = fopen(fn_obs, "r")) == NULL) {
                nll_puterr("FATAL ERROR: opening observations file.");
                return (EXIT_ERROR_FILEIO);
            } else {
                NumFilesOpen++;
            }

            obs_line_array = (char **) calloc(1000, sizeof (char *));
            while (fp_obs != NULL && fgets(line, 4 * MAXLINE, fp_obs) != NULL) {
                obs_line_array[n_obs_lines] = (char *) malloc(4 * MAXLINE);
                strcpy(obs_line_array[n_obs_lines], line);
                n_obs_lines++;
            }
            fclose(fp_obs);
            NumFilesOpen--;
        }



        /** ===========================================================================
         *  Call function invocation of NLLoc.
         *
         *  NOTE: the parameter loc_list_head in NLLoc() returns the location results:
         *  LocNode **loc_list_head - pointer to pointer to head of list of LocNodes containing Location's for
         *  located events (see phaseloclist.h)
         *  *loc_list_head must be initialized to NULL on first call to NLLoc()
         */

        return_locations = 1;
        return_oct_tree_grid = 1;
        return_scatter_sample = 1;
        istat = NLLoc(pid_main, NULL, (char **) param_line_array, n_param_lines, (char **) obs_line_array, n_obs_lines, return_locations, return_oct_tree_grid, return_scatter_sample, &loc_list_head);



        /** ===========================================================================
         *  Write location results to disk for each returned event.
         *  A production program might scan and process these location results entirely in memory.
         */

        LocNode* locNode = NULL;
        char frootname[FILENAME_MAX];
        char fname[FILENAME_MAX];

        // loop over returned location results
        int id = 0;
        while ((locNode = getLocationFromLocList(loc_list_head, id)) != NULL) {

            sprintf(frootname, "out/%3.3d", id_filename);
            sprintf(fname, "%s.loc.hyp", frootname);

            // NOTE: the angles in locNode->plocation->ellipsoid, ellipse and
            //    phase/arrival angles in locNode->plocation->parrivals are in NLL internal coordinates
            //    and must be transformed to geographic directions if there is a non-zero rotAngle in TRANS.
            //    For example:
            //          ellipsoid_az1 = rect2latlonAngle(0, phypo->ellipsoid.az1);
            //          ellipsoid_az2 = rect2latlonAngle(0, phypo->ellipsoid.az2);
            //    These two transformations ARE NOT applied in WriteLocation() used below!
            //          sta_azim = rect2latlonAngle(0, parr->azim);
            //          ray_azim = rect2latlonAngle(0, parr->ray_azim);
            //    These two transformations ARE applied in WriteLocation() used below.

            // write NLLoc Hypocenter-Phase file to disk
            if ((istat = WriteLocation(NULL, locNode->plocation->phypo, locNode->plocation->parrivals,
                    locNode->plocation->narrivals, fname, 1, 1, 0, locNode->plocation->pgrid, 0)) < 0) {
                nll_puterr2("ERROR: writing location to event file: %s", fname);
            }

            // write NLLoc location Grid Header file to disk
            if ((istat = WriteGrid3dHdr(locNode->plocation->pgrid, NULL, frootname, "loc")) < 0) {
                nll_puterr2("ERROR: writing grid header to disk: %s", frootname);
            }

            // write NLLoc location Oct tree structure of locaiton likelihood values to disk
            if (return_oct_tree_grid) {
                sprintf(fname, "%s.loc.octree", frootname);
                FILE *fpio;
                if ((fpio = fopen(fname, "w")) != NULL) {
                    istat = writeTree3D(fpio, locNode->plocation->poctTree);
                    fclose(fpio);
                    sprintf(MsgStr, "Oct tree structure written to file : %d nodes", istat);
                    nll_putmsg(1, MsgStr);
                }
            }

            // write NLLoc binary Scatter file to disk
            if (return_scatter_sample) {
                sprintf(fname, "%s.loc.scat", frootname);
                FILE *fpio;
                if ((fpio = fopen(fname, "w")) != NULL) {
                    // write scatter file header informaion
                    fseek(fpio, 0, SEEK_SET);
                    fwrite(&(locNode->plocation->phypo->nScatterSaved), sizeof (int), 1, fpio);
                    float ftemp = (float) locNode->plocation->phypo->probmax;
                    fwrite(&ftemp, sizeof (float), 1, fpio);
                    // skip header record
                    fseek(fpio, 4 * sizeof (float), SEEK_SET);
                    // write scatter samples
                    fwrite(locNode->plocation->pscatterSample, 4 * sizeof (float), locNode->plocation->phypo->nScatterSaved, fpio);
                    fclose(fpio);
                }
            }

            id++;
            id_filename++;
        }

        // clean up
        freeLocList(loc_list_head, 1);
        loc_list_head = NULL;
        int i;
        for (i = 0; i < n_param_lines; i++)
            free(param_line_array[i]);
        free(param_line_array);
        param_line_array = NULL;
        n_param_lines = 0;
        for (i = 0; i < n_obs_lines; i++)
            free(obs_line_array[i]);
        free(obs_line_array);
        obs_line_array = NULL;
        n_obs_lines = 0;


        /** ===========================================================================
         *  END of oop over all observation files specified in command line
         */
    }

    return (istat);

}



