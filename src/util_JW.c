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


/* util.c

   AJL utility functions.

*/


/*
	by Anthony Lomax
	Geosciences Azur, Valbonne, France
*/


/*	history:

	21 SEP 1998	AJL	Extracted from GridLib.c
*/



#ifdef EXTERN_MODE
#define	EXTERN_TXT extern
#else
#define EXTERN_TXT
#endif

#include "util.h"
FILE* NLL_Message_Stream = NULL;
FILE* NLL_Error_Stream = NULL;

/*** function to copy file by Jan Wiszniowski 2022-01-31*/
void copy_file(char* in_name, char* out_name)
{
#ifdef _GNU_SOURCE
    static char* RF = "r"
    static char* WF = "w"
#elif defined (WIN32)
    static char* RF = "rb";
    static char* WF = "wb";
#else
#error SYSTEM???
#endif

    FILE* in_file = NULL;
    FILE* out_file = NULL;
    char buf[1024];
    int bytes_in;
    in_file = fopen(in_name, RF);
    if (in_file == NULL) {
        return;
    }
    out_file = fopen(out_name, RF);
    if (out_file == NULL) {
        fclose(in_file);
        return;
    }
    while ((bytes_in = fread(buf, 1, 1024, in_file)) > 0)
    {
        fwrite(buf, 1, bytes_in, out_file);
    }
    fclose(in_file);
    fclose(out_file);
}

/*** function to display correct command line usage */

void disp_usage(const char * progname, const char * options)
{
	if (NLL_Error_Stream) fprintf(NLL_Error_Stream, "Usage: %s %s\n", progname, options);
}



/*** function to display error message */

void nll_puterr(const char *pm)
{
	if (NLL_Error_Stream) {
		fprintf(NLL_Error_Stream, "[0]%s: %s\n", prog_name, pm);
		fflush(NLL_Error_Stream);
	}
}



/*** function to display error message */

void nll_puterr2(const char *pmessage1, const char *pmessage2)
{
	if (NLL_Error_Stream) {
		fprintf(NLL_Error_Stream, "[0]%s: %s: %s\n", prog_name, pmessage1, pmessage2);
		fflush(NLL_Error_Stream);
	}
}



/*** function to display message */

void nll_putmsg(int imsg_level, const char *pm)
{
	if (NLL_Message_Stream)	fprintf(NLL_Message_Stream, "[%d]%s\n", imsg_level+1, pm);
}



/*** function to display message */

void nll_putmsg2(int imsg_level, const char *pmessage1, const char *pmessage2)
{
	if (NLL_Message_Stream) fprintf(NLL_Message_Stream, "[%d]%s: %s\n", imsg_level+1, pmessage1, pmessage2);
}



/*** function to display program name, version, date */

void DispProgInfo()
{
	sprintf(MsgStr, "%s (%s v%s %s) %s",
		prog_name, package_name, prog_ver, prog_date, prog_copyright);
	nll_putmsg(1, MsgStr);
}



/*** function to check that int val is in range */

int checkRangeInt(const char * name, const char * param, int val,
		int checkMin, int min, int checkMax, int max)
{
	int stat = 0;

	if (checkMin && val < min) {
		sprintf(MsgStr,
			"ERROR: %s param %s: value: %d is less than min value: %d",
			name, param, val, min);
		nll_puterr(MsgStr);
		stat = -1;
	}

	if (checkMax && val > max) {
		sprintf(MsgStr,
			"ERROR: %s param %s: value: %d is greater than max value: %d",
			name, param, val, max);
		nll_puterr(MsgStr);
		stat = 1;
	}

	return(stat);

}



/*** function to check that double val is in range */

int checkRangeDouble(const char * name, const char * param, double val,
		int checkMin, double min, int checkMax, double max)
{
	int stat = 0;

	if (checkMin && val < min - VERY_SMALL_DOUBLE) {
		sprintf(MsgStr,
			"ERROR: %s param %s: value: %lf is less than min value: %lf",
			name, param, val, min);
		nll_puterr(MsgStr);
		stat = -1;
	}

	if (checkMax && val > max + VERY_SMALL_DOUBLE) {
		sprintf(MsgStr,
			"ERROR: %s param %s: value: %lf is greater than max value: %lf",
			name, param, val, max);
		nll_puterr(MsgStr);
		stat = 1;
	}

	return(stat);

}
