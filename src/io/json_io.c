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




#define JW_GLOBAL_CONTROL_STRUCT
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
#define FILE_DOES_NOT_EXIST 1
#define FILE_TO_LARGE 2
#define FILE_READ_ERROR 3

static size_t json_data_buffer_size;
static char *json_data_buffer = NULL;
static struct jReadElement theArray;
static const char *pArray;
static int curent_index;

static char strtmp[1024];

/** function to allocate data buffer
 *
 */
int json_alloc_data_buffer(size_t buf_size) {

    // check that json_data_buffer is NULL or was previously freed, THIS SHOULD NOT HAPPEN!
    if (json_data_buffer != NULL) {
        nll_puterr("ERROR: json_alloc_data_buffer(): json_data_buffer not NULL, was previously freed, THIS SHOULD NOT HAPPEN!");
        return (-2);
    }

    // 1 GiB; best not to load a whole large file in one string
    if (buf_size > 1073741824) {
        nll_puterr("ERROR: json_io: json_alloc_data_buffer(): FILE_TO_LARGE");
        return (-1);
    }

    json_data_buffer = (char *) malloc(buf_size + 1);
    //printf("DEBUG: json_alloc_data_buffer: size: %ld\n", buf_size + 1);

    return (0);

}

/** function to free data buffer
 *
 */
void json_free_data_buffer() {

    if (json_data_buffer != NULL) {
        free(json_data_buffer);
    }
    json_data_buffer = NULL;

}

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
            nll_puterr("ERROR: json_io: c_read_file(): reading JSON file: FILE_TO_LARGE");
            return NULL;
        }

        buffer = (char *) malloc(length + 1);

        if (length) {
            read_length = fread(buffer, 1, length, fp_in);

            if (length != read_length) {
                free(buffer);
                *err = FILE_READ_ERROR;
                nll_puterr("ERROR: json_io: c_read_file(): reading JSON file: FILE_READ_ERROR");
                return NULL;
            }
        }

        fclose(fp_in);

        *err = FILE_OK;
        buffer[length] = '\0';
        *f_size = length;

    } else {
        *err = FILE_DOES_NOT_EXIST;
        nll_puterr("ERROR: json_io: c_read_file(): reading JSON file: FILE_DOES_NOT_EXIST");
        return NULL;
    }

    return buffer;
}

/** function to read next arrival in INGV_JSON
 *
 *
 */
int json_read_next_arrival_INGV(HypoDesc* phypo, FILE* fp_obs, ArrivalDesc *arrival, int nfirst) {

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

        // check that json_data_buffer is NULL or was previously freed, THIS SHOULD NOT HAPPEN!
        if (json_data_buffer != NULL) {
            nll_puterr("ERROR: json_read_nll_control(): json_data_buffer not NULL, was previously freed, THIS SHOULD NOT HAPPEN!");
            return (OBS_FILE_INTERNAL_ERROR);
        }

        // read json to string
        int err;
        json_data_buffer = c_read_file(fp_obs, &err, &json_data_buffer_size);
        //printf("DEBUG: json_data_buffer:\n%s\n", json_data_buffer);
        if (json_data_buffer == NULL) {
            nll_puterr("ERROR: json_read_nll_control(): reading JSON file");
            return (OBS_FILE_INTERNAL_ERROR);
        }

        if (err) {
            // process error
            if (err == FILE_DOES_NOT_EXIST) {
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
            json_free_data_buffer();
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
        json_free_data_buffer();
        return (OBS_FILE_END_OF_INPUT);
    }

    return (1);

}

/** function to read nll-control JSON format
 *
 *
 */
char **json_read_nll_control(FILE* fp_control, int *pn_param_lines) {

    char **param_line_array; // array of NLLoc control file lines (set to NULL if fn_control_main not NULL)
    int n_param_lines = 0; // number of elements (parameter lines) in array param_line_array (use 0 if fn_control_main not NULL)

    // check that json_data_buffer is NULL or was previously freed, THIS SHOULD NOT HAPPEN!
    if (json_data_buffer != NULL) {
        nll_puterr("ERROR: json_read_nll_control(): json_data_buffer not NULL, was previously freed, THIS SHOULD NOT HAPPEN!");
        goto cleanup_error;
    }

    // read json to string
    int err;
    json_data_buffer = c_read_file(fp_control, &err, &json_data_buffer_size);
    if (json_data_buffer == NULL) {
        nll_puterr("ERROR: json_read_nll_control(): reading JSON file");
        return (NULL);
    }

    // read array of nll control statements
    jRead(json_data_buffer, "{'nll-control'{'statements'", &theArray);
    if (theArray.elements < 1) {
        // try with extended nesting
        jRead(json_data_buffer, "{'data'{'nll-control'{'statements'", &theArray);
    }
    //printf("DEBUG: theArray:\n%s\n", (char *) theArray.pValue);
    //printf("DEBUG: theArray.dataType:%d JREAD_ARRAY %d theArray.elements %d\n", theArray.dataType, JREAD_ARRAY, theArray.elements);
    if (theArray.elements < 1) {
        nll_puterr("ERROR: json_read_nll_control(): no JSON \"nll-control\" \"statements\" elements found.");
        goto cleanup_error;
    }
    if (theArray.dataType == JREAD_ARRAY) {
        pArray = (char *) theArray.pValue;
    } else {
        nll_puterr("ERROR: json_read_nll_control(): JSON \"nll-control\" \"statements\" array not found.");
        goto cleanup_error;
    }

    // allocate statement array
    n_param_lines = theArray.elements;
    if ((param_line_array = (char **) calloc(n_param_lines, sizeof (char*))) == NULL) {
        nll_puterr("ERROR: json_read_nll_control(): allocating memory for param_line_array.");
        goto cleanup_error;
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
            nll_puterr("ERROR: json_read_nll_control(): JSON \"nll-control\" \"statements\" element is not a string.");
            for (int i = 0; i < n_param_lines; i++) {
                free(param_line_array[i]);
            }
            free(param_line_array);
            goto cleanup_error;
        }
        curent_index++;
    }

    /* DEBUG
    for (int n = 0; n < n_param_lines; n++) {
        printf("DEBUG: statements %d : %s\n", n, param_line_array[n]);
    }
    // END DEBUG */

    *pn_param_lines = n_param_lines;
    json_free_data_buffer();
    return (param_line_array);

cleanup_error:
    json_free_data_buffer();
    return (NULL);

}



// helper functions for json_write_NLL_location

/** helper function for json_write_NLL_location
 *
 * special processing for PHASE block
 *
 */
char *json_write_NLL_location_PHASE(char *nextLine) {

    /*
     PHASE ID Ins Cmp On Pha  FM Date     HrMn   Sec     Err  ErrMag    Coda      Amp       Per       PriorWt  >   TTpred    Res       Weight    StaLoc(X  Y         Z)        SDist    SAzim  RAz  RDip RQual    Tcorr       TTerr
    IV_NRCA      --   HHZ  ? P      ? 20211011 0715   28.7300 GAU  1.00e+04 -1.00e+00 -1.00e+00 -1.00e+00    1.0000 >    3.4586   -0.2814    0.0000   50.2232   92.7749   -0.9270   11.8189 324.97 325.0 130.2  9     0.0000    0.0000
    ...
    END_PHASE
     */

    static char *phase_hdr[] = {"ID", "Ins", "Cmp", "On", "Pha", "FM", "Date", "HrMn", "Sec", "Err", "ErrMag", "Coda", "Amp", "Per", "PriorWt", ">", "TTpred", "Res", "Weight", "StaLoc(X", "Y", "Z)", "SDist", "SAzim", "RAz", "RDip", "RQual", "Tcorr", "TTerr"};

    char datetime[1024] = "\0";

    // start phase array
    jwObj_array("PHASE");

    // skip PHASE header line
    if (nextLine) *nextLine = '\n'; // then restore newline-char, just to be tidy
    char *currline = nextLine ? (nextLine + 1) : NULL;

    while (currline) {
        nextLine = strchr(currline, '\n');
        if (nextLine) *nextLine = '\0'; // temporarily terminate the current line
        if (strlen(currline) > 0 && currline[0] != '\n') { // not a blank line
            // parse line into tokens
            char *token = strtok(currline, " ");
            if (token != NULL && strcmp(token, "END_PHASE") == 0) {
                if (nextLine) *nextLine = '\n'; // then restore newline-char, just to be tidy
                currline = nextLine ? (nextLine + 1) : NULL;
                break;
            }
            int hdr_index = 0;
            long itest;
            double dtest;
            while (token != NULL) {
                if (hdr_index == 0) {
                    jwArr_object();
                }
                // check for line types that need special processing
                if (strcmp(phase_hdr[hdr_index], "Date") == 0) {
                    // date-time
                    strncpy(datetime, token, 4); // year
                    datetime[4] = '\0';
                    strncat(strcat(datetime, "-"), token + 4, 2); // month
                    strncat(strcat(datetime, "-"), token + 2, 2); // day
                    token = strtok(NULL, " ");
                    hdr_index++;
                    if (token != NULL) {
                        strncat(strcat(datetime, "T"), token, 2); // hour
                        strncat(strcat(datetime, ":"), token + 2, 2); // min
                        token = strtok(NULL, " ");
                        hdr_index++;
                        if (token != NULL) {
                            strncat(strcat(datetime, ":"), token, 16); // sec
                            jwObj_string("time", datetime);
                            token = strtok(NULL, " ");
                            hdr_index++;
                        }
                    }
                } else if (strcmp(phase_hdr[hdr_index], ">") == 0) {
                    // skip obs > calc separator
                    token = strtok(NULL, " ");
                    hdr_index++;
                } else if (strncmp(phase_hdr[hdr_index], "StaLoc", 6) == 0) {
                    // station location
                    jwObj_object("StaLoc");
                    token = strtok(NULL, " ");
                    hdr_index++;
                    if (token != NULL) {
                        jwObj_raw("X", token);
                        token = strtok(NULL, " ");
                        hdr_index++;
                        if (token != NULL) {
                            jwObj_raw("Y", token);
                            token = strtok(NULL, " ");
                            hdr_index++;
                            if (token != NULL) {
                                jwObj_raw("Z", token);
                                jwEnd();
                                token = strtok(NULL, " ");
                                hdr_index++;
                            }
                        }
                    }
                } else {
                    // standard processing
                    // write without quotes in integer or double value
                    if (sscanf(token, "%ld", &itest) == 1 || sscanf(token, "%lf", &dtest) == 1) {
                        jwObj_raw(phase_hdr[hdr_index], token);
                    } else {
                        jwObj_string(phase_hdr[hdr_index], token);
                    }
                    token = strtok(NULL, " ");
                    hdr_index++;
                }
            }
            if (hdr_index > 0) {
                jwEnd();
            }
        }

        if (nextLine) *nextLine = '\n'; // then restore newline-char, just to be tidy
        currline = nextLine ? (nextLine + 1) : NULL;
    }

    // end array
    jwEnd();

    return (nextLine);

}

/** helper function for json_write_NLL_location
 *
 * special processing for NLLOC line
 *
 */
void json_write_NLL_location_NLLOC(char *currline) {

    // NLLOC "out/test_20211012/loc/Italy.20211011.071528.grid0" "LOCATED" "Location completed."

    // JSON object already created for this line type
    char *token = strtok(NULL, "\"");
    if (token != NULL) {
        jwObj_string("fileroot", token);
        token = strtok(NULL, "\"");
        if (token != NULL) {
            token = strtok(NULL, "\""); // extra "
            if (token != NULL) {
                jwObj_string("locStat", token);
                token = strtok(NULL, "\""); // extra "
                if (token != NULL) {
                    token = strtok(NULL, "\"");
                    if (token != NULL) {
                        jwObj_string("locStatComm", token);
                    }
                }
            }
        }
    }
}

/** helper function for json_write_NLL_location
 *
 * special processing for PUBLIC_ID line
 *
 */
void json_write_NLL_location_PUBLIC_ID(char *currline) {

    // PUBLIC_ID None

    // JSON object already created for this line type
    char *token = strtok(NULL, "");
    if (token != NULL) {
        jwObj_string("public_id", token);
    }
}

/** helper function for json_write_NLL_location
 *
 * special processing for SIGNATURE line
 *
 */
void json_write_NLL_location_SIGNATURE(char *currline) {

    // SIGNATURE "B2_ItalyLoc Project   obs:Phases_Example.json   NLLoc:v7.00.13(17Jan2022)   run:07Feb2022 15h20m43"

    // JSON object already created for this line type
    char *token = strtok(NULL, "\"");
    if (token != NULL) {
        jwObj_string("signature", token);
    }
}

/** helper function for json_write_NLL_location
 *
 * special processing for COMMENT line
 *
 */
void json_write_NLL_location_COMMENT(char *currline) {

    // COMMENT "JSON test Italy"

    // JSON object already created for this line type
    char *token = strtok(NULL, "\"");
    if (token != NULL) {
        jwObj_string("comment", token);
    }
}

/** helper function for json_write_NLL_location
 *
 * special processing for GRID line
 *
 */
void json_write_NLL_location_GRID(char *currline) {

    // GRID  1251 1201 86  -625 -600 -5  1 1 1 PROB_DENSITY

    // JSON object already created for this line type
    char *token = strtok(NULL, " ");
    if (token != NULL) {
        jwObj_raw("numx", token);
        token = strtok(NULL, " ");
        if (token != NULL) {
            jwObj_raw("numy", token);
            token = strtok(NULL, " ");
            if (token != NULL) {
                jwObj_raw("numz", token);
                token = strtok(NULL, " ");
                if (token != NULL) {
                    jwObj_raw("origx", token);
                    token = strtok(NULL, " ");
                    if (token != NULL) {
                        jwObj_raw("origy", token);
                        token = strtok(NULL, " ");
                        if (token != NULL) {
                            jwObj_raw("origz", token);
                            token = strtok(NULL, " ");
                            if (token != NULL) {
                                jwObj_raw("dx", token);
                                token = strtok(NULL, " ");
                                if (token != NULL) {
                                    jwObj_raw("dy", token);
                                    token = strtok(NULL, " ");
                                    if (token != NULL) {
                                        jwObj_raw("dz", token);
                                        token = strtok(NULL, " ");
                                        if (token != NULL) {
                                            jwObj_string("chr_type", token);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

/** helper function for json_write_NLL_location
 *
 * special processing for GEOGRAPHIC line
 *
 */
void json_write_NLL_location_GEOGRAPHIC(char *currline) {

    // GEOGRAPHIC  OT 2021 10 11  07 15 25.552831  Lat 42.745959 Long 13.196256 Depth 12.464844

    char datetime[1024] = "\0";

    // JSON object already created for this line type
    char *token = strtok(NULL, " ");
    if (token != NULL) {
        // "OT" ignore
        token = strtok(NULL, " ");
        if (token != NULL) {
            strncpy(datetime, token, 4); // year
            token = strtok(NULL, " ");
            if (token != NULL) {
                strncat(strcat(datetime, "-"), token, 2); // month
                token = strtok(NULL, " ");
                if (token != NULL) {
                    strncat(strcat(datetime, "-"), token, 2); // day
                    token = strtok(NULL, " ");
                    if (token != NULL) {
                        strncat(strcat(datetime, "T"), token, 2); // hour
                        token = strtok(NULL, " ");
                        if (token != NULL) {
                            strncat(strcat(datetime, ":"), token, 2); // min
                            token = strtok(NULL, " ");
                            if (token != NULL) {
                                strncat(strcat(datetime, ":"), token, 16); // sec
                                jwObj_string("origin_time", datetime);
                                token = strtok(NULL, " ");
                                char *token2;
                                long itest;
                                double dtest;
                                while (token != NULL) {
                                    token2 = strtok(NULL, " ");
                                    if (token2 != NULL) {
                                        // write without quotes in integer or double value
                                        if (sscanf(token2, "%ld", &itest) == 1 || sscanf(token2, "%lf", &dtest) == 1) {
                                            jwObj_raw(token, token2);
                                        } else {
                                            jwObj_string(token, token2);
                                        }
                                    }
                                    token = strtok(NULL, " ");
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

/** helper function for json_write_NLL_location
 *
 * special processing for QUALITY line
 *
 */
void json_write_NLL_location_QUALITY(char *currline) {

    // QUALITY  Pmax 1.7719e+18 MFmin 8.43777 MFmax 8.95647 RMS 0.168386 Nphs 10 Gap 63.9267 Dist 16.4412 Mamp -9.90 0 Mdur -9.90 0

    // JSON object already created for this line type
    char *token = strtok(NULL, " ");
    char *token2;
    long itest;
    double dtest;
    int nobject = 0;
    while (token != NULL && nobject < 7) {
        token2 = strtok(NULL, " ");
        if (token2 != NULL) {
            // write without quotes in integer or double value
            if (sscanf(token2, "%ld", &itest) == 1 || sscanf(token2, "%lf", &dtest) == 1) {
                jwObj_raw(token, token2);
            } else {
                jwObj_string(token, token2);
            }
        }
        nobject++;
        token = strtok(NULL, " ");
    }

    //  Mamp -9.90 0 Mdur -9.90 0
    char cnum_rdgs_name[1024] = "\0";
    char *clastname;
    if (token != NULL) {
        // Mamp
        clastname = token;
        token = strtok(NULL, " ");
        if (token != NULL) {
            jwObj_raw(clastname, token);
            token = strtok(NULL, " ");
            if (token != NULL) {
                strncpy(cnum_rdgs_name, clastname, 32);
                strncat(cnum_rdgs_name, "_Nrdgs", 32);
                jwObj_raw(cnum_rdgs_name, token);
            }
        }
    }
    token = strtok(NULL, " ");
    if (token != NULL) {
        // Mdur
        clastname = token;
        token = strtok(NULL, " ");
        if (token != NULL) {
            jwObj_raw(clastname, token);
            token = strtok(NULL, " ");
            if (token != NULL) {

                strncpy(cnum_rdgs_name, clastname, 32);
                strncat(cnum_rdgs_name, "_Nrdgs", 32);
                jwObj_raw(cnum_rdgs_name, token);
            }
        }
    }
}

/** helper function for json_write_NLL_location
 *
 * special processing for FOCALMECH line
 *
 */
void json_write_NLL_location_FOCALMECH(char *currline) {

    // FOCALMECH  Hyp  42.745959 13.196256 12.464844 Mech  0 0 0 mf  0 nObs -1

    // JSON object already created for this line type
    char *token = strtok(NULL, " ");
    if (token != NULL) {
        // Hyp
        jwObj_object(token);
        token = strtok(NULL, " ");
        if (token != NULL) {
            jwObj_raw("Lat", token);
            token = strtok(NULL, " ");
            if (token != NULL) {
                jwObj_raw("Lon", token);
                token = strtok(NULL, " ");
                if (token != NULL) {
                    jwObj_raw("Depth", token);
                    jwEnd();
                    token = strtok(NULL, " ");
                    if (token != NULL) {
                        // Mech
                        jwObj_object(token);
                        token = strtok(NULL, " ");
                        if (token != NULL) {
                            jwObj_raw("dipDir", token);
                            token = strtok(NULL, " ");
                            if (token != NULL) {
                                jwObj_raw("dipAng", token);
                                token = strtok(NULL, " ");
                                if (token != NULL) {
                                    jwObj_raw("rake", token);
                                    jwEnd();
                                    token = strtok(NULL, " ");
                                    char *token2;
                                    long itest;
                                    double dtest;
                                    while (token != NULL) {
                                        token2 = strtok(NULL, " ");
                                        if (token2 != NULL) {
                                            // write without quotes in integer or double value
                                            if (sscanf(token2, "%ld", &itest) == 1 || sscanf(token2, "%lf", &dtest) == 1) {
                                                jwObj_raw(token, token2);
                                            } else {
                                                jwObj_string(token, token2);
                                            }
                                        }
                                        token = strtok(NULL, " ");
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

/** function to read nll-control JSON format
 *
 *
 */
int json_write_NLL_location(char *file_str_loc_in, size_t stream_size, FILE * fp_json_out) {

    int err;

    // check that json_data_buffer is NULL or was previously freed, THIS SHOULD NOT HAPPEN!
    if (json_data_buffer != NULL) {
        nll_puterr("ERROR: json_read_nll_control(): json_data_buffer not NULL, was previously freed, THIS SHOULD NOT HAPPEN!");
        return (OBS_FILE_INTERNAL_ERROR);
    }

    // allocate json_data_buffer
    size_t buf_size;
    if (stream_size > 0) {
        buf_size = 10 * stream_size;
    } else {
        buf_size = 100000000; // very large value 100MB
    }
    if (json_alloc_data_buffer(buf_size) != 0) {
        nll_puterr("ERROR: json_read_nll_control(): allocating json_data_buffer");
        return (-1);
    }

    /*
     {
    "nll-hypocenter-phase": {
        "lines": [
            "CONTROL 1 54321",
     */

    json_data_buffer_size = buf_size;
    jwOpen(json_data_buffer, json_data_buffer_size, JW_OBJECT, JW_PRETTY); // start root object

    jwObj_object("nll-hypocenter-phase");
    {

        jwObj_object("data");
        {
            char *currline = file_str_loc_in;
            //printf("DEBUG: TP -3 <%s>\n", currline);
            while (currline) {
                char *nextLine = strchr(currline, '\n');
                if (nextLine) *nextLine = '\0'; // temporarily terminate the current line
                //printf("DEBUG: currline <%s>\n", currline);
                if (strlen(currline) > 0 && currline[0] != '\n') { // not a blank line
                    // parse line into tokens
                    char *token = strtok(currline, " ");
                    if (token != NULL && strcmp(token, "END_NLLOC") == 0) {
                        if (nextLine) *nextLine = '\n'; // then restore newline-char, just to be tidy
                        currline = nextLine ? (nextLine + 1) : NULL;
                        break;
                    }
                    if (token != NULL) { // line tag
                        // special processing for PHASE block
                        if (strcmp(token, "PHASE") == 0) {
                            nextLine = json_write_NLL_location_PHASE(nextLine);
                            currline = nextLine ? (nextLine + 1) : NULL;
                        } else {
                            jwObj_object(token);
                            // check for line types that need special processing
                            if (strcmp(token, "NLLOC") == 0) {
                                json_write_NLL_location_NLLOC(currline);
                            } else if (strcmp(token, "PUBLIC_ID") == 0) {
                                json_write_NLL_location_PUBLIC_ID(currline);
                            } else if (strcmp(token, "SIGNATURE") == 0) {
                                json_write_NLL_location_SIGNATURE(currline);
                            } else if (strcmp(token, "COMMENT") == 0) {
                                json_write_NLL_location_COMMENT(currline);
                            } else if (strcmp(token, "GRID") == 0) {
                                json_write_NLL_location_GRID(currline);
                            } else if (strcmp(token, "GEOGRAPHIC") == 0) {
                                json_write_NLL_location_GEOGRAPHIC(currline);
                            } else if (strcmp(token, "QUALITY") == 0) {
                                json_write_NLL_location_QUALITY(currline);
                            } else if (strcmp(token, "FOCALMECH") == 0) {
                                json_write_NLL_location_FOCALMECH(currline);
                            } else {
                                char *line_name = token;
                                // other lines
                                token = strtok(NULL, " ");
                                if (token != NULL && (strcmp(line_name, "SEARCH") == 0 || strcmp(line_name, "TRANSFORM") == 0)) {
                                    // first token is unnamed type
                                    jwObj_string("type", token);
                                    token = strtok(NULL, " ");
                                }
                                char *token2;
                                long itest;
                                double dtest;
                                while (token != NULL) {
                                    token2 = strtok(NULL, " ");
                                    if (token2 != NULL) {
                                        // write without quotes in integer or double value
                                        if (sscanf(token2, "%ld", &itest) == 1 || sscanf(token2, "%lf", &dtest) == 1) {
                                            jwObj_raw(token, token2);
                                        } else {
                                            jwObj_string(token, token2);
                                        }
                                    }
                                    token = strtok(NULL, " ");
                                }
                            }
                            jwEnd();
                        }
                    }
                }
                if (nextLine) *nextLine = '\n'; // then restore newline-char, just to be tidy
                currline = nextLine ? (nextLine + 1) : NULL;
            }

            jwEnd(); // end of data object, back to root object
        }

        jwEnd();
    }

    err = jwClose(); // close and get error code
    if (err != JWRITE_OK)
        nll_puterr2("ERROR: json_write_NLL_location(): closing JSON buffer", jwErrorToString(err));

    // write to output
    fprintf(fp_json_out, "%s", json_data_buffer);

    json_free_data_buffer();

    return (0);
}