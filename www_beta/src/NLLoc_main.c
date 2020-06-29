/*
 * Copyright (C) 1999-2007 Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
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


/*   NLLoc.c

	Program to do global search earthquake location in 3-D models

*/

/*-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
161 Allee du Micocoulier, 06370 Mouans-Sartoux, France
tel: +33(0)493752502  e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/


/*
	history:	(see also http://alomax.net/nlloc -> Updates)

	ver 01    17DEC2007  AJL  Created NLLoc_main and function call NLLoc()

	see NLLoc1.c and NLLocLib.c


.........1.........2.........3.........4.........5.........6.........7.........8

*/





#ifdef CUSTOM_ETH
#define PNAME  "NLLoc(ETH)"
#else
#define PNAME  "NLLoc"
#endif

#include "GridLib.h"
#include "ran1/ran1.h"
#include "velmod.h"
#include "GridMemLib.h"
#include "calc_crust_corr.h"
#include "phaseloclist.h"
#include "otime_limit.h"
#include "NLLocLib.h"

#ifdef CUSTOM_ETH
#include "custom_eth/eth_functions.h"
#endif


// function declarations


/** program to perform global search event locations */

#ifdef CUSTOM_ETH
#define NARGS_MIN 3
#define ARG_DESC "<control file> <snap_pid> <snap_param_file>"
#else
#define NARGS_MIN 2
#define ARG_DESC "<control file>"
#endif

int main(int argc, char *argv[])
{

	int istat;

	char fn_control_main[MAXLINE];	// control file name
	char pid_main[255];	// string process id


	// set program name
	strcpy(prog_name, PNAME);

	// check command line for correct usage
	if (argc < NARGS_MIN) {
		disp_usage(prog_name, ARG_DESC);
		return(EXIT_ERROR_USAGE);
	}

	// set control file
	strcpy(fn_control_main, argv[1]);

#ifdef CUSTOM_ETH
	// SH 02/27/2004 added snap_pid
	if (argc > 2)
		strcpy(pid_main, argv[2]);
	else
		strcpy(pid_main, "000");
	/* SH 02AUG2004 not needed any more
	// AJL 20040527 added snap param file
	if (argc > 3)
		strcpy(snap_param_file, argv[3]);
	else
		strcpy(snap_param_file, "snap_param.txt");
	*/
#endif

	// run NLLoc
	istat = NLLoc(pid_main, fn_control_main, NULL, -1, NULL, -1, 0, 0, 0, NULL);

	return(istat);

}




