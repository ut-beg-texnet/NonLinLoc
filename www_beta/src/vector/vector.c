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


/*  Vector mathematics functions
	Anthony Lomax

	Nov 1992
*/


#include <math.h>
#include "vector.h"

/*** function to calculate dot product of 3-D vectors (A o B) */

double dot_product_3d(a1, a2, a3, b1, b2, b3)		    
double a1, a2, a3, b1, b2, b3;
{
	return(a1 * b1 + a2 * b2 + a3 * b3);
}


/*** function to calculate cross product of 3-D vectors (A X B) */

int cross_product_3d(a1, a2, a3, b1, b2, b3, p1, p2, p3)		    
double a1, a2, a3, b1, b2, b3;
double *p1, *p2, *p3;
{
	*p1 = a2 * b3 - a3 * b2;
	*p2 = a3 * b1 - a1 * b3;
	*p3 = a1 * b2 - a2 * b1;

	return(0);
}


/*** function to calculate magnitude of 3-D vectors ||A|| */

double magnitude_3d(a1, a2, a3)		    
double a1, a2, a3;
{
	return(sqrt(a1 * a1 + a2 * a2 + a3 * a3));
}
