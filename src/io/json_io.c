/*
 * Copyright (C) 1999-2021 Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
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


/* json_io.c

   NonLinLoc json read / write functions.

 */


/*-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
Mouans-Sartoux, France
e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/


/*
        history:

        ver 01    12Oct2021  AJL  Original version

 */




#ifdef EXTERN_MODE
#define EXTERN_TXT extern
#else
#define EXTERN_TXT
#endif


#include "jRead.h"
#include "jWrite.h"

#include "GridLib.h"
#include "ran1/ran1.h"
#include "velmod.h"
#include "GridMemLib.h"
#include "calc_crust_corr.h"
#include "phaseloclist.h"
#include "otime_limit.h"
#include "NLLocLib.h"

#include "json_io.h"

#define FILE_OK 0
#define FILE_NOT_EXIST 1
#define FILE_TO_LARGE 2
#define FILE_READ_ERROR 3

static size_t json_data_buffer_size;
static char *json_data_buffer = NULL;
static struct jReadElement theArray;
static const char *pArray;
static int curent_index;

static char strtmp[1024];

/** function to read contents of file to string
 *
 * from: https://stackoverflow.com/questions/174531/how-to-read-the-content-of-a-file-to-a-string-in-c
 *
 */
char *c_read_file(FILE *fp_in, int *err, size_t *f_size) {

    char *buffer;
    size_t length;
    size_t read_length;

    if (fp_in) {
        fseek(fp_in, 0, SEEK_END);
        length = ftell(fp_in);
        fseek(fp_in, 0, SEEK_SET);

        // 1 GiB; best not to load a whole large file in one string
        if (length > 1073741824) {
            *err = FILE_TO_LARGE;

            return NULL;
        }

        buffer = (char *) malloc(length + 1);

        if (length) {
            read_length = fread(buffer, 1, length, fp_in);

            if (length != read_length) {
                free(buffer);
                *err = FILE_READ_ERROR;

                return NULL;
            }
        }

        fclose(fp_in);

        *err = FILE_OK;
        buffer[length] = '\0';
        *f_size = length;

    } else {
        *err = FILE_NOT_EXIST;

        return NULL;
    }

    return buffer;
}

/** function to read next arrival in INGV_JSON
 *
 *
 */
int jReadNextArrival_INGV(HypoDesc* phypo, FILE* fp_obs, ArrivalDesc *arrival, int nfirst) {

    // INGV json format 20211012
    // assumes a single event in each input file
    /*
    {
        "data": {
            "hyp2000_conf": [
            ...
            ],
            "model": [
            ...
            ],
            "output": "json",
            "phases": [
                {
                    "net": "IV",
                    "sta": "MTRA",
                    "cha": "EHZ",
                    "loc": "--",
                    "arrival_time": "2021-10-11 07:15:29.570",
                    "isc_code": "P",
                    "firstmotion": null,
                    "emersio": null,
                    "weight": 1,
                    "amplitude": null,
                    "ampType": "1"
                },
            ...
            ]
        }
    }
     */

    // read json file to char buffer on first invocation for an event, free buffer at end of event

    if (nfirst) {

        // read json to string
        int err;
        json_data_buffer = c_read_file(fp_obs, &err, &json_data_buffer_size);

        //printf("DEBUG: json_data_buffer:\n%s\n", json_data_buffer);
        if (json_data_buffer == NULL) {
            return (OBS_FILE_END_OF_INPUT);
        }

        if (err) {
            // process error
            if (err == FILE_NOT_EXIST) {
                nll_puterr("ERROR: in jReadNextArrival_INGV: FILE_NOT_EXIST");
            } else if (err == FILE_TO_LARGE) {
                nll_puterr("ERROR: in jReadNextArrival_INGV: FILE_TO_LARGE");
            } else if (err == FILE_READ_ERROR) {
                nll_puterr("ERROR: in jReadNextArrival_INGV: FILE_READ_ERROR");
            }
        }

        // read array of phases
        jRead(json_data_buffer, "{'data'{'phases'", &theArray);
        //printf("DEBUG: theArray:\n%s\n", (char *) theArray.pValue);
        //printf("DEBUG: theArray.dataType:%d JREAD_ARRAY %d theArray.elements %d\n", theArray.dataType, JREAD_ARRAY, theArray.elements);
        if (theArray.dataType == JREAD_ARRAY) {
            curent_index = 0;
            pArray = (char *) theArray.pValue;
        }
    }

    //printf("DEBUG: curent_index: %d theArray.elements %d nfirst %d\n", curent_index, theArray.elements, nfirst);
    //printf("DEBUG: curent_index: %d nfirst %d\n", curent_index, nfirst);
    if (curent_index < theArray.elements) {
        struct jReadElement arrayElement;
        pArray = jReadArrayStep(pArray, &arrayElement);
        struct jReadElement phaseElement;
        // net
        jRead(arrayElement.pValue, "{'net'", &phaseElement);
        if (phaseElement.dataType == JREAD_STRING) {
            strncpy(arrival->label, phaseElement.pValue, phaseElement.bytelen);
            arrival->label[phaseElement.bytelen] = '\0';
            strncat(arrival->label, "_", 1);
        }
        // sta
        jRead(arrayElement.pValue, "{'sta'", &phaseElement);
        if (phaseElement.dataType == JREAD_STRING) {
            strncat(arrival->label, phaseElement.pValue, phaseElement.bytelen);
        }
        // location / instrument
        jRead(arrayElement.pValue, "{'loc'", &phaseElement);
        if (phaseElement.dataType == JREAD_STRING) {
            strncpy(arrival->inst, phaseElement.pValue, phaseElement.bytelen);
        }
        // channel
        jRead(arrayElement.pValue, "{'cha'", &phaseElement);
        if (phaseElement.dataType == JREAD_STRING) {
            strncpy(arrival->comp, phaseElement.pValue, phaseElement.bytelen);
        }
        // onset
        jRead(arrayElement.pValue, "{'emersio'", &phaseElement);
        if (phaseElement.dataType == JREAD_STRING) {
            strncpy(arrival->onset, phaseElement.pValue, phaseElement.bytelen);
        }
        // phase
        jRead(arrayElement.pValue, "{'isc_code'", &phaseElement);
        if (phaseElement.dataType == JREAD_STRING) {
            strncpy(arrival->phase, phaseElement.pValue, phaseElement.bytelen);
        }
        // first motion
        jRead(arrayElement.pValue, "{'firstmotion'", &phaseElement);
        if (phaseElement.dataType == JREAD_STRING) {
            strncpy(arrival->first_mot, phaseElement.pValue, phaseElement.bytelen);
        }
        // time
        jRead(arrayElement.pValue, "{'arrival_time'", &phaseElement);
        //                 "arrival_time": "2021-10-11 07:15:33.510",
        //printf("DEBUG: arrival_time: %s\n", (char *) phaseElement.pValue);
        int istat = sscanf(phaseElement.pValue, "%4d-%2d-%2d %2d:%2d:%lf",
                &arrival->year, &arrival->month, &arrival->day, &arrival->hour, &arrival->min, &arrival->sec);
        //printf("DEBUG: arrival_time istat: %d\n", istat);
        if (istat < 6) { // no arrival time, skip phase
            curent_index++;
            return (OBS_FILE_INVALID_PHASE);
        }
        //
        // quality
        jRead(arrayElement.pValue, "{'weight'", &phaseElement);
        //printf("DEBUG: phaseElement.dataType %d phaseElement.bytelen %d \n", phaseElement.dataType, phaseElement.bytelen);
        if (phaseElement.dataType == JREAD_NUMBER) {
            strncpy(strtmp, phaseElement.pValue, phaseElement.bytelen);
            //printf("DEBUG: weight %s phaseElement.bytelen %d\n", strtmp, phaseElement.bytelen);
            arrival->quality = atoi(strtmp);
            Qual2Err(arrival);
        }
        //
        // amplitude
        jRead(arrayElement.pValue, "{'amp'", &phaseElement);
        //printf("DEBUG: phaseElement.dataType %d phaseElement.bytelen %d \n", phaseElement.dataType, phaseElement.bytelen);
        if (phaseElement.dataType == JREAD_NUMBER) {
            strncpy(strtmp, phaseElement.pValue, phaseElement.bytelen);
            //printf("DEBUG: amplitude %s phaseElement.bytelen %d\n", strtmp, phaseElement.bytelen);
            arrival->amplitude = atof(strtmp);
        }
        //
        WriteArrival(stdout, arrival, IO_ARRIVAL_OBS);

        curent_index++;
    }


    // reached end of event
    if (curent_index >= theArray.elements) {
        if (json_data_buffer != NULL) {
            free(json_data_buffer);
        }
        json_data_buffer = NULL;
        return (OBS_FILE_END_OF_INPUT);
    }

    return (1);

}

/** function to read nll-control JSON format
 *
 *
 */
char **read_nll_control_json(FILE* fp_control, int *pn_param_lines) {

    char **param_line_array; // array of NLLoc control file lines (set to NULL if fn_control_main not NULL)
    int n_param_lines = 0; // number of elements (parameter lines) in array param_line_array (use 0 if fn_control_main not NULL)

    // read json to string
    int err;
    json_data_buffer = c_read_file(fp_control, &err, &json_data_buffer_size);

    // read array of nll control statements
    jRead(json_data_buffer, "{'nll-control'{'statements'", &theArray);
    if (theArray.elements < 1) {
        // try with extended nesting
        jRead(json_data_buffer, "{'data'{'nll-control'{'statements'", &theArray);
    }
    //printf("DEBUG: theArray:\n%s\n", (char *) theArray.pValue);
    printf("DEBUG: theArray.dataType:%d JREAD_ARRAY %d theArray.elements %d\n", theArray.dataType, JREAD_ARRAY, theArray.elements);
    if (theArray.elements < 1) {
        nll_puterr("ERROR: read_nll_control_json(): no JSON \"nll-control\" \"statements\" elements found.");
        return (NULL);
    }
    if (theArray.dataType == JREAD_ARRAY) {
        pArray = (char *) theArray.pValue;
    } else {
        nll_puterr("ERROR: read_nll_control_json(): JSON \"nll-control\" \"statements\" array not found.");
        return (NULL);
    }

    // allocate statement array
    n_param_lines = theArray.elements;
    if ((param_line_array = (char **) calloc(n_param_lines, sizeof (char*))) == NULL) {
        nll_puterr("ERROR: read_nll_control_json(): allocating memory for param_line_array.");
        return (NULL);
    }
    //allocate space for each string
    for (int i = 0; i < n_param_lines; i++) {
        param_line_array[i] = (char *) malloc(MAXLINE_LONG * sizeof (char));
        //printf("DEBUG: i:%d %s\n", i, param_line_array[curent_index]);
    }

    // load each statement to param_line_array
    curent_index = 0;
    while (curent_index < n_param_lines) {
        struct jReadElement arrayElement;
        pArray = jReadArrayStep(pArray, &arrayElement);
        //printf("DEBUG: curent_index:%d %s\n", curent_index, param_line_array[curent_index]);
        if (arrayElement.dataType == JREAD_STRING) {
            strncpy(param_line_array[curent_index], arrayElement.pValue, arrayElement.bytelen);
        } else {
            nll_puterr("ERROR: read_nll_control_json(): JSON \"nll-control\" \"statements\" element is not a string.");
            for (int i = 0; i < n_param_lines; i++) {
                free(param_line_array[i]);
            }
            free(param_line_array);
            return (NULL);
        }
        curent_index++;
    }

    /* DEBUG
    for (int n = 0; n < n_param_lines; n++) {
        printf("DEBUG: statements %d : %s\n", n, param_line_array[n]);
    }
    // END DEBUG */

    *pn_param_lines = n_param_lines;
    return (param_line_array);
}