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
