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



#include <string.h>

#include "util.h"

char package_name[MAXLINE];
char prog_name[MAXLINE];
char prog_ver[MAXLINE];
char prog_date[MAXLINE];
char prog_copyright[MAXLINE];
int message_flag;
char MsgStr[100 * MAXLINE];



/*** function to copy file by Jan Wiszniowski 2022-01-31*/
void copy_file(char* in_name, char* out_name)
{
#if defined(WIN32) || defined(WIN64)
    static char* RF = "rb";
    static char* WF = "wb";
#else
    static char* RF = "r";
    static char* WF = "w";
#endif


    FILE* in_file = NULL;
    FILE* out_file = NULL;
    char buf[1024];
    int bytes_in;
    in_file = fopen(in_name, RF);
    if (in_file == NULL) {
        return;
    }
    out_file = fopen(out_name, WF);
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
	fprintf(stderr, "Usage: %s %s\n", progname, options);
}



/*** function to display error message */

void nll_puterr(const char *pm)
{
	fprintf(stderr, "%s: %s\n", prog_name, pm);
        fflush(stderr);
}



/*** function to display error message */

void nll_puterr2(const char *pmessage1, const char *pmessage2)
{
	fprintf(stderr, "%s: %s: %s\n", prog_name, pmessage1, pmessage2);
        fflush(stderr);
}



/*** function to display message */

void nll_putmsg(int imsg_level, const char *pm)
{
	if (message_flag >= imsg_level)
		fprintf(stdout, "%s\n", pm);
}



/*** function to display message */

void nll_putmsg2(int imsg_level, const char *pmessage1, const char *pmessage2)
{
	if (message_flag >= imsg_level)
		fprintf(stdout, "%s: %s\n", pmessage1, pmessage2);
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



/*** function to copy a word read with scanf "%s" into dest, cut with a warning if too long */
// 20261127 Bug fix (Mustafa Comoglu & Claude): scanf "%s" has no length limit, so a too-long word
//   overflowed dest. Callers now read the word into a buffer as long as the input line, then call this.
// An empty word means "%s" did not read it (it never stores an empty word); dest is then left unchanged, as before.
// what names the input in the warning; what == NULL cuts without a warning (for lines the caller will reject anyway).
// Returns 1 if the word was cut, else 0.

int CopyWordCut(char *dest, size_t dest_size, const char *word, const char *what)
{
	if (word[0] == '\0')
		return(0);

	size_t len = strlen(word);
	if (len < dest_size) {
		memcpy(dest, word, len + 1);
		return(0);
	}

	memcpy(dest, word, dest_size - 1);
	dest[dest_size - 1] = '\0';
	if (what != NULL) {
		fprintf(stderr, "%s: WARNING: %s: \"%s\" is longer than %zu characters, CUT to \"%s\"\n",
			prog_name, what, word, dest_size - 1, dest);
		fflush(stderr);
	}

	return(1);
}
