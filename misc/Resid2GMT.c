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


/*   Resid2GMT.c

	Program to generate GMT psxy data from NLLoc .stat file

	output in GMT readable format

*/

//------------------------------------------------------------/
// Anthony Lomax           | email: lomax@faille.unice.fr     /
// UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax /
// 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        /
// 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        /
//------------------------------------------------------------/


/*
	history:

	ver 01    18Dec1998  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

*/



#define PNAME  "Resid2GMT"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#define MAXLINE 100


// globals


// functions

int Convert2GMT(char* fn_outroot);


Program to generate GMT psxy data from NLLoc .stat file


main(int argc, char *argv[])
{

	int istat;
	int c;
	char fn_outroot[MAXLINE];


	/* read each set of residuals */

	while (1)
	{
		/* get residual group name */
		while ((istat = scanf("%s", fn_outroot)) != 1 && istat != EOF)
			;
printf("FILENAME: %s\n", fn_outroot);

		if (istat == EOF)
			break;

		while((c = getchar()) != '*' && c != EOF)
			;
		if (Convert2GMT(fn_outroot) < 0)
			break;
	}


	exit(0);

}



/** function to convert residuals list to gmt psxy (residual vs. Nsta) and pstext file  */

int Convert2GMT(char* fn_outroot)
{
	int istat, nRdg;
	int c = 0;
	char code[MAXLINE];
	char fname[MAXLINE];
	FILE *fp_xy_out, *fp_text_out;

	char staName[20];
	char phase[20];
	double resid;



	sprintf(fname, "%s.resid.xy", fn_outroot);
	if ((fp_xy_out = fopen(fname, "w")) == NULL) {
		fprintf(stderr, 
			"ERROR: Cannot open xy output file <%s>.\n",
				fname);
		return(-1);
    	}
	sprintf(fname, "%s.resid.text", fn_outroot);
	if ((fp_text_out = fopen(fname, "w")) == NULL) {
		fprintf(stderr, 
			"ERROR: Cannot open text output file <%s>.\n",
				fname);
		return(-1);
    	}


	/* read residuals */

	nRdg = 0;
	while (scanf("%s", code) == 1 && strncmp(code, "*", 1))
	{
		scanf("%s %s %*s %lf", 
			staName, phase, &(resid));
		fprintf(fp_xy_out, "%d %lf\n", nRdg + 1, resid);
		fprintf(fp_text_out, "%d 1.0 10 0 4 2 %s\n", 
				nRdg + 1, staName);
		fprintf(fp_text_out, "%d -1.0 10 0 4 2 %s\n", 
				nRdg + 1, phase);

		while((c = getchar()) != '\n' && c !=EOF)
			;

		nRdg++;
	}


	fclose(fp_xy_out);
	fclose(fp_text_out);

	return(0);

}









//------------------------------------------------------------/
// Anthony Lomax           | email: lomax@faille.unice.fr     /
// UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax /
// 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        /
// 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        /
//------------------------------------------------------------/

