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


#include <stdio.h>
#include <stdlib.h>


/*------------------------------------------------------------/ */
/* function declarations */
/*------------------------------------------------------------/ */

int json_read_next_arrival_INGV(HypoDesc* phypo, FILE* fp_obs, ArrivalDesc *arrival, int nfirst);
char **json_read_nll_control(FILE* fp_control, int *pn_param_lines);
int json_write_NLL_location(char *file_str_loc_in, size_t stream_size, FILE *fp_json_out);

/*------------------------------------------------------------/ */


