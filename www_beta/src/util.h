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


/* util.h

   AJL utility functions.

*/


/*
	by Anthony Lomax
	Geosciences Azur, Valbonne, France
*/



#include <stdlib.h>
#include <stdio.h>


#ifndef MAXLINE
#define MAXLINE 101
#endif

#ifndef VERY_SMALL_DOUBLE
#define VERY_SMALL_DOUBLE 1.0e-30
#endif

EXTERN_TXT char package_name[MAXLINE];
EXTERN_TXT char prog_name[MAXLINE];
EXTERN_TXT char prog_ver[MAXLINE];
EXTERN_TXT char prog_date[MAXLINE];
EXTERN_TXT char prog_copyright[MAXLINE];
EXTERN_TXT int message_flag;
EXTERN_TXT char MsgStr[10 * MAXLINE];


/*** function to display command usage */
void disp_usage(const char * , const char *);

/*** function to display error message */
void nll_puterr(const char *);

/*** function to display error message */
void nll_puterr2(const char *, const char *);

/*** function to display message */
void nll_putmsg(int , const char *);

/*** function to display message */
void nll_putmsg2(int imsg_level, const char *pmessage1, const char *pmessage2);

/*** function to display program information */
void DispProgInfo();

/*** function to check that int val is in range */
int checkRangeInt(const char * name, const char * param, int val,
	int checkMin, int min, int checkMax, int max);

/*** function to check that double val is in range */
int checkRangeDouble(const char * name, const char * param, double val,
	int checkMin, double min, int checkMax, double max);

/* misc structures */

/* 3D vector or point */
// AEH/AJL 20080709 - already in geometry.h
/*
struct Vect3D {
	double x;
	double y;
	double z;
};
*/




