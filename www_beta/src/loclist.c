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


/*   loclist.c

	loclist management and search functions

*/

/*-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
161 Allee du Micocoulier, 06370 Mouans-Sartoux, France
tel: +33(0)493752502  e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/


/*
	history:

	ver 01    27JUN2006  AJL  Original version

*/



#define EXTERN_MODE 1

#include "GridLib.h"
#include "phaseloclist.h"
//#include "ran1.h"


/** function to create a new LocNode from a NLL HypoDesc structure - sorted by id */

LocNode *addLocationToLocList(LocNode **phead, Location *plocation, int id)
{

	int i;
	double first_phase_time;
	double first_phase_time_min;
	LocNode* addr;
	LocNode* addr_test;


	addr = (LocNode*) malloc(sizeof(LocNode));
	if (addr == NULL) {	/* memory allocation problem */
		printf( "loclist: ERROR: allocating memory for LocNode.\n");
		return(addr);
	}

	// initalize members
	addr->plocation = plocation;
	addr->id = id;
	first_phase_time_min = VERY_LARGE_DOUBLE;
	for (i = 0; i < plocation->narrivals; i++) {
		if ((first_phase_time = getPhaseTimeValue(plocation->parrivals + i)) < first_phase_time_min)
			first_phase_time_min = first_phase_time;
	}
	addr->first_phase_time = first_phase_time_min;
	// NOT NEEDED? addr->phsNodeList = phsNodeList;

	if (*phead == NULL) {		/* no elements in list */
		*phead = addr->next = addr->prev = addr;
	} else {
		/* assume locations are added approximately in order, search from end of list backwards */
		addr_test = (*phead)->prev;
		while (addr->id  < addr_test->id ) {
			addr_test = addr_test->prev;
			if (addr_test == (*phead)->prev)
				break;
		}
		addr->prev = addr_test;
		addr->next = addr_test->next;
		addr->prev->next = addr;
		addr->next->prev = addr;
		/* check if arrival belongs at head of list */
		if (addr->next == *phead && (addr->id  < (*phead)->id ))
			*phead = addr;
	}


	return(addr);
}



/** function to remove a LocNode from a LocList */

LocNode *removeLocationFromLocList(LocNode *head, LocNode* addr, int freeHypoDesc)
{

	if (freeHypoDesc)
		freeLocation(addr->plocation);

	if (addr == head) {
		if (addr->next == head)	// removed last remaining node
			head = NULL;
		else 	// removed head, make next head
			head = addr->next;
	}

	addr->prev->next = addr->next;
	addr->next->prev = addr->prev;

	free(addr);

	return(head);

}



/** function to get a LocNode with specified id from a LocList */

LocNode *getLocationFromLocList(LocNode *head, int id)
{
	LocNode* addr;

	if (head == NULL)
		return(NULL);

	addr = head;
	do {
		if (addr->id == id)
			return(addr);
	} while ((addr = addr->next) != head);

	return(NULL);

}



/** function to free an entire LocList */

int freeLocList(LocNode *head, int freeHypoDesc)
{

	if (head == NULL)
		return(0);

	int n = 0;

	while ((head = removeLocationFromLocList(head, head, freeHypoDesc)) != NULL)
		n++;

	return(n + 1);

}



/** function to create a new Location */

Location *newLocation(HypoDesc *phypo, ArrivalDesc* parrivals, int narrivals, GridDesc *pgrid, Tree3D* poctTree, float *pscatterSample)
{

	Location *plocation;


	plocation = (Location*) malloc(sizeof(Location));
	if (plocation == NULL) {	/* memory allocation problem */
		printf( "loclist: ERROR: allocating memory for location.\n");
		return(NULL);
	}
	plocation->phypo = phypo;
	plocation->parrivals = parrivals;
	plocation->narrivals = narrivals;
	plocation->pgrid = pgrid;
	plocation->poctTree = poctTree;
	plocation->pscatterSample = pscatterSample;

	return(plocation);

}



/** function to free a Location */

void freeLocation(Location *plocation)
{

	if (plocation == NULL)
		return;

	free(plocation->phypo);
	free(plocation->parrivals);
	free(plocation->pgrid);
	freeTree3D(plocation->poctTree, 1);
	free(plocation->pscatterSample);
	free(plocation);

}


/** function to create a new LocNode from a NLL-HypocenterPhase file */
/*
LocNode *addNLLLocationFileToLocList(LocNode **phead, char *filein, int id, PhsNode *phsNodeList)
{

	int i;
	HypoDesc *phypo = NULL;
	ArrivalDesc* parrivals;
	int narrivals;
	GridDesc *pgrid;
	Location *plocation;


	phypo = (HypoDesc*) malloc(sizeof(HypoDesc));
	pgrid = (GridDesc*) malloc(sizeof(GridDesc));

	if ( GetHypLoc(NULL, filein, phypo, arrivals_tmp, &narrivals, 1, pgrid, -1) != 0 ) {
		free(phypo);
		free(pgrid);
		return(NULL);
	}

	parrivals = (ArrivalDesc*) calloc(narrivals, sizeof(ArrivalDesc));
	if (parrivals == NULL) {	// memory allocation problem
		printf( "loclist: ERROR: allocating memory for arrivals.\n");
	} else {
		for (i = 0; i < narrivals; i++)
			parrivals[i] = arrivals_tmp[i];
	}

	plocation = newLocation(phypo, parrivals, narrivals, pgrid);
	if (plocation == NULL) {	// memory allocation problem
		printf( "loclist: ERROR: constructin new location.\n");
		return(NULL);
	}

	return(addLocationToLocList(phead, plocation, id));

}
*/


/** function to display location list */

int writeLocList(LocNode *head, FILE *out)
{
	LocNode *addr;

	if ((addr = head) == NULL) {
		printf("LocList:  EMPTY.\n");
		return(0);
	}

	fprintf(out, "LocList:\n");
	do {
		writeLocNode(addr, out);
	} while ((addr = addr->next) != head);
	fprintf(out, "\n");

	return(0);

}


/** function to display location */

void writeLocNode(LocNode *addr, FILE *out)
{

	int iWriteArrivals = 0;
	int iWriteEndLoc = 1;
	int iWriteMinimal = 1;

	fprintf(out, "N=%d ", addr->id); /*DEBUG*/
	WriteLocation(out, addr->plocation->phypo, addr->plocation->parrivals,
		      addr->plocation->narrivals, NULL,
		      iWriteArrivals, iWriteEndLoc, iWriteMinimal,
		      addr->plocation->pgrid, -1);

}



/** function to convert hypocenter date/time to double time value */

static struct tm time_1970 = {0, 0, 1, 1, 0, 70};
static time_t time_1970_seconds = LONG_MIN;
static time_t TIME_T_INVALID = LONG_MIN;

double getLocTimeValue(HypoDesc *phypo)
{

	struct tm hypo_time;
	time_t time_seconds;
	double time_value;

	/* initalize time_1970_seconds */
	if (time_1970_seconds == TIME_T_INVALID)
		time_1970_seconds = mktime(&time_1970);


	hypo_time.tm_year = (*phypo).year - 1900;
	hypo_time.tm_mon = (*phypo).month - 1;
	hypo_time.tm_mday = (*phypo).day;
	hypo_time.tm_hour = (*phypo).hour;
	hypo_time.tm_min = (*phypo).min;
	hypo_time.tm_sec = 0;

	time_seconds = mktime(&hypo_time);

	if (time_seconds == -1)
		return(INVALID_DOUBLE);

	time_value = difftime(time_seconds, time_1970_seconds) + (*phypo).sec;


	fprintf(stdout, "Hypo: time_1970_seconds=%f ", (double) time_1970_seconds); /*DEBUG*/
	fprintf(stdout, "time_seconds=%f ", (double) time_seconds); /*DEBUG*/
	fprintf(stdout, "(*phypo).sec=%f ", (*phypo).sec); /*DEBUG*/
	fprintf(stdout, "time_value=%f ", time_value); /*DEBUG*/
	fprintf(stdout, "\n"); /*DEBUG*/


	return(time_value);

}



/** function to find all locations with first phase in a specified time window */

static LocNode *locNodesTmp[MAX_NUM_LOCATIONS];

LocNode **findLocsWithFirstPhaseInTimeWindow(LocNode *head, double tmin, double tmax)
{

	int nfound, i;
	LocNode *addr;
	LocNode **locNodes;


	if ((addr = head) == NULL) {
		return(NULL);
	}


	nfound = 0;
	do {
		if (addr->first_phase_time >= tmin && addr->first_phase_time <= tmax)
			locNodesTmp[nfound++] = addr;
	} while ((addr = addr->next) != head);

	if (nfound <= 0)
		return(NULL);

	locNodes = (LocNode**) calloc(nfound + 1, sizeof(LocNode*));
	if (locNodes == NULL) {	/* memory allocation problem */
		printf( "loclist: ERROR: allocating memory for locNodes.\n");
		return(NULL);
	}
	for (i = 0; i < nfound; i++)
		locNodes[i] = locNodesTmp[i];

	locNodes[nfound] = NULL;

	return(locNodes);

}




/** function to find a phase in the arrivals of a Location */

ArrivalDesc *findArrivalInLocation(Location *plocation, ArrivalDesc *arrivalKey)
{

	int n;
	double phaseTime;


	phaseTime = getPhaseTimeValue(arrivalKey);

	for (n = 0; n < plocation->narrivals; n++)
	{
		if (compareArrivals(plocation->parrivals + n, arrivalKey, 0) &&
				  compareTimeValue(getPhaseTimeValue(plocation->parrivals + n), phaseTime) == 0)
			return(plocation->parrivals + n);
	}

	return(NULL);


}






/** function to test if an arrivals occurs earlier in time than one or more arrivals used for location (apriori_weight > 0) */

int arrivalBeforeLastLocateArrival(ArrivalDesc *arrivalKey, Location *plocation)
{

	int n;
	double phaseTime;


	phaseTime = getPhaseTimeValue(arrivalKey);

	for (n = 0; n < plocation->narrivals; n++)
	{
		if ((plocation->parrivals[n].apriori_weight > 0.0) &&
				   (compareTimeValue(getPhaseTimeValue(plocation->parrivals + n), phaseTime) > 0))
			return(1);
	}

	return(0);


}




