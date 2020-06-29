/*
 * Copyright (C) 2002 Anthony Lomax <anthony@alomax.net>
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


/*   grad2layer.c

	Program to convert gradient model (NonLinLoc) to layered model (doubDiff)

*/

/*------------------------------------------------------------*/
/* Anthony Lomax           | email: anthony@alomax.net        */
/* Scientific Software     | web: www.alomax.net              */
/* Mouans-Sartoux, FRANCE  |                                  */
/*------------------------------------------------------------*/


/*
	history:

	ver 01    15FEB2002  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

*/



#define PNAME  "grad2layer"

/* defines */

/* globals  */


/*** Program to convert gradient model to layered model */

main(int argc, char *argv[])
{

	int i, index;
	double depth;
	

// ---------- describe model here ----------

	double ratio = 1.73;
	
	int nSegments = 6;
	double grad_vel[][3] = {
		{ 0.0, 3.45, 0.7100},
		{ 1.0, 4.16, 0.4000},
		{ 3.0, 4.96, 0.1933},
		{ 6.0, 5.54, 0.0538},
		{14.0, 5.97, 0.0455},
		{25.0, 7.98, 0.00  }
	};
	
	double layer_thickness = 0.5;
	int num_layers = 21;
	
// ---------- end - describe model here ----------

	printf("* NLAY RATIO\n");
	printf("%d %lf\n", num_layers, ratio);
	printf("* TOP\n");
	for (i = 0; i < num_layers; i++)
		printf("%.1lf ", layer_thickness * (double) i);
	printf("\n* VEL\n");
	index = 0;
	for (i = 0; i < num_layers; i++) {
		depth = layer_thickness * (double) i + 0.000001;
		while (index < nSegments - 1 && depth > grad_vel[index + 1][0])
			index++;
		printf("%.2lf ", grad_vel[index][1] + grad_vel[index][2] * (depth - grad_vel[index ][0]));
	}
	printf("\n*\n");

}



/*------------------------------------------------------------/ */
/* Anthony Lomax           | email: lomax@faille.unice.fr     / */
/* UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax / */
/* 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        / */
/* 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        / */
/*------------------------------------------------------------/ */

