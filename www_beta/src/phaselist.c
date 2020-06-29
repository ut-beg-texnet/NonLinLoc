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


/*   phaselist.c

	phaselist management and search functions

*/

/*-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
161 Allee du Micocoulier, 06370 Mouans-Sartoux, France
tel: +33(0)493752502  e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/


/*
	history:

	ver 01    23JUN2006  AJL  Original version

*/



#define EXTERN_MODE 1

#include "GridLib.h"
#include "phaseloclist.h"
//#include "ran1.h"


/** function to create a copy of a PhsNode */

PhsNode *copyPhsNode(PhsNode *original)
{

	int n;
	PhsNode* addr;

	addr = (PhsNode*) malloc(sizeof(PhsNode));
	if (addr == NULL) {	/* memory allocation problem */
		printf( "phaselist: ERROR: allocating memory for PhsNode.\n");
		return(NULL);
	}

	addr->prev = addr->next = NULL;
	addr->parrival = (ArrivalDesc*) malloc(sizeof(ArrivalDesc));
	*(addr->parrival) = *(original->parrival);
	addr->id = original->id;
	addr->phase_time = original->phase_time;

	addr->passoc_locations_size = original->passoc_locations_size;
	addr->passoc_locations = (int*) calloc(addr->passoc_locations_size, sizeof(int));
	if (addr->passoc_locations  == NULL) {	/* memory allocation problem */
		printf( "phaselist: ERROR: allocating memory for PhsNode->passoc_locations .\n");
		return(NULL);
	}
	n = -1;
	do {
		n++;
		*(addr->passoc_locations + n) = *(original->passoc_locations + n);
	} while (*(original->passoc_locations + n) != -1);

	return(addr);
}



/** function to create a new PhsNode from a NLL-Phase string */

PhsNode *addNLLPhaseStringToPhaseList(PhsNode **phead, char *phase_string, int id, int addDuplicates)
{

	int istat;
	ArrivalDesc *parrival = NULL;
	PhsNode *phsNode;


	parrival = (ArrivalDesc*) malloc(sizeof(ArrivalDesc));

	istat = ReadArrival(phase_string, parrival, IO_ARRIVAL_OBS);
	if (istat != 1) {
		free(parrival);
		return(NULL);
	}

	phsNode = addArrivalToPhaseList(phead, parrival, id, addDuplicates);
	if (phsNode == DUPLICATE_PAHSE_FOUND) {
		free(parrival);
	}

	return(phsNode);

}


/** function to get a PhsNode with specified id from a PhaseList */

PhsNode *getPhsNodeFromPhaseList(PhsNode *head, int id)
{
	PhsNode* addr;

	if (head == NULL)
		return(NULL);

	addr = head;
	do {
		if (addr->id == id)
			return(addr);
	} while ((addr = addr->next) != head);

	return(NULL);

}




/** function to create a new PhsNode from a NLL-Phase ArrivalDesc structure and insert in a PhsNode list*/

PhsNode *addArrivalToPhaseList(PhsNode **phead, ArrivalDesc* parrival, int id, int addDuplicates)
{

	PhsNode* addr;

	if (!addDuplicates) {
		if(findPhase(*phead, parrival) != NULL)
			return(DUPLICATE_PAHSE_FOUND);
	}

	addr = (PhsNode*) malloc(sizeof(PhsNode));
	if (addr == NULL) {	/* memory allocation problem */
		printf( "phaselist: ERROR: allocating memory for PhsNode.\n");
		return(addr);
	}

	addr->parrival = parrival;
	addr->id = id;
	addr->phase_time = getPhaseTimeValue(parrival);
	addr->passoc_locations_size = INIT_NUM_ASSOC_LOC_PER_PHASE;
	addr->passoc_locations = (int*) calloc(addr->passoc_locations_size, sizeof(int));
	if (addr->passoc_locations == NULL) {	/* memory allocation problem */
		printf( "phaselist: ERROR: allocating memory for PhsNode->passoc_locations .\n");
		return(NULL);
	}
	addr->passoc_locations[0] = -1;

	*phead = addPhsNodeToPhaseList(*phead, addr);

	return(addr);

}



/** function to add a PhsNode to a PhsNode list */

PhsNode *addPhsNodeToPhaseList(PhsNode *head, PhsNode* addr)
{

	PhsNode* addr_test;


	if (head == NULL) {		/* no elements in list */
		head = addr->next = addr->prev = addr;
	} else {
		/* assume arrival are added approximately in order, search from end of list backwards */
		addr_test = head->prev;
		while (addr->phase_time < addr_test->phase_time) {
			addr_test = addr_test->prev;
			if (addr_test == head->prev)
				break;
		}
		addr->prev = addr_test;
		addr->next = addr_test->next;
		addr->prev->next = addr;
		addr->next->prev = addr;
		/* check if arrival belongs at head of list */
		if (addr->next == head && (addr->phase_time <= head->phase_time))
			head = addr;
	}

	return(head);
}



/** function to remove a PhsNode from a PhaseList */

PhsNode *removeNodeFromPhaseList(PhsNode *head, PhsNode* addr, int freeArrivalDesc)
{

	if (freeArrivalDesc)
		free(addr->parrival);

	addr->prev->next = addr->next;
	addr->next->prev = addr->prev;

	if (addr == head) {
		if (addr->next == head)	// removed last node
			head = NULL;
		else 	// removed head
			head = addr->next;
	}

	free(addr);

	return(head);

}



/** function to free an entire PhaseList */

int freePhaseList(PhsNode *head, int freeArrivalDesc)
{

	int n = 0;

	while (head != NULL) {
		head = removeNodeFromPhaseList(head, head, freeArrivalDesc);
		n++;
	}

	return(n + 1);

}




/** function to display phase list */

int writePhaseList(PhsNode *head, FILE *out)
{
	int n;
	PhsNode *addr;

	if ((addr = head) == NULL) {
		printf("PhaseList:  EMPTY.\n");
		return(0);
	}

	fprintf(out, "PhaseList:\n");
	do {
		fprintf(out, "N=%d ", addr->id); /*DEBUG*/
		fprintf(out, "t=%f assoc=", addr->phase_time); /*DEBUG*/
		n = 0;
		do {
			fprintf(out, "%d,", *(addr->passoc_locations + n)); /*DEBUG*/
		} while (*(addr->passoc_locations + n++) != -1);
		fprintf(out, "  "); /*DEBUG*/
		WriteArrival(out, addr->parrival, IO_ARRIVAL_OBS);
	} while ((addr = addr->next) != head);
	fprintf(out, "\n");

	return(0);

}


/** function to convert arrival date/time to double time value */

static struct tm time_1970 = {0, 0, 1, 1, 0, 70, 0, 0, 0};
static time_t time_1970_seconds = LONG_MIN;
static time_t TIME_T_INVALID = LONG_MIN;

double getPhaseTimeValue(ArrivalDesc *parrival)
{

	struct tm arrival_time;
	time_t time_seconds;
	double time_value;

	/* initalize time_1970_seconds */
	if (time_1970_seconds == TIME_T_INVALID)
		time_1970_seconds = mktime(&time_1970);


	arrival_time.tm_year = parrival->year - 1900;
	arrival_time.tm_mon = parrival->month - 1;
	arrival_time.tm_mday = parrival->day;
	arrival_time.tm_hour = parrival->hour;
	arrival_time.tm_min = parrival->min;
	arrival_time.tm_sec = 0;
	arrival_time.tm_isdst = 0;

	time_seconds = mktime(&arrival_time);

	if (time_seconds == -1)
		return(INVALID_DOUBLE);

	time_value = difftime(time_seconds, time_1970_seconds) + parrival->sec;


	/*
	fprintf(stdout, "Phase: time_1970_seconds=%f ", (double) time_1970_seconds);
	fprintf(stdout, "time_seconds=%f ", (double) time_seconds);
	fprintf(stdout, "parrival->sec=%f ", parrival->sec);
	fprintf(stdout, "time_value=%f ", time_value);
	fprintf(stdout, "\n");
	*/

	return(time_value);

}


/** function to compare time values */

#define TIME_TOLERANCE 0.0001

int compareTimeValue(double t1, double t2)
{

	double diff;

	diff = t1 - t2;

	if (diff > TIME_TOLERANCE)
		return(1);
	if (diff < -TIME_TOLERANCE)
		return(-1);

	return(0);

}




/** function to add or remove a location from associated location list in a PhsNode */

int addRemoveLocationInAssocLocationsList(PhsNode *addr, int locID, int addLocID)
{

	int index, isize, ifound, newsize;
	int *plist, *pnewlist;


	plist = addr->passoc_locations;
	isize = addr->passoc_locations_size;

	// search list for locID
	ifound = 0;
	index = 0;
	for ( ; index < isize; index++) {
		if (plist[index] < 0)
			break;
		if (plist[index] == locID) {
			ifound = 1;
			break;
		}
	}

	if (ifound && addLocID)
		return(0);

	// reached end of list (passoc_locations[isize - 1] must be -1)
	if (index == isize) {
		printf( "phaselist: ERROR: PhsNode->passoc_locations list not terminated by -1. (1) : sizeof(plist) / sizeof(int) %d\n", isize);
		WriteArrival(stdout, addr->parrival, IO_ARRIVAL_OBS);
		return(-1);
	} else if (index == isize - 1) {
		if (addLocID) {
			// adding, need to enlarge list
			newsize = isize + INIT_NUM_ASSOC_LOC_PER_PHASE;
			pnewlist  = (int*) realloc(plist, newsize * sizeof(int));
			if (pnewlist == NULL) {	/* memory allocation problem */
				printf( "phaselist: ERROR: re-allocating memory for PhsNode->passoc_locations .\n");
				return(-1);
			}
			plist = addr->passoc_locations = pnewlist;
			// update size of array
			isize = addr->passoc_locations_size = newsize;
		} else {
			// removing, not found
			return(0);
		}
	}

	// need to add or remove
	if (addLocID) {
		// adding, not found
		plist[index] = locID;
		plist[index + 1] = -1;
	} else if (ifound) {
		// removing, found - compact list
		index--;
		do {
			index++;
			if (index == isize) {
				printf( "phaselist: ERROR: PhsNode->passoc_locations list not terminated by -1. (2): ");
				WriteArrival(stdout, addr->parrival, IO_ARRIVAL_OBS);
				return(-1);
			}
			plist[index] = plist[index + 1];
		} while (plist[index] >= 0);
	}

	return(1);

}




/** function to remove assocated flag fro a location from all phases in a phaselist */

int removeLocationAssociation(PhsNode *head, int locID, double tmin_nomimnal, double tmax_nomimnal)
{

	PhsNode *addr;

	double tmin, tmax;

	if ((addr = head) == NULL) {
		return(0);
	}

	// enlarge window to guarantee check of all possible associated phases
	tmin = tmin_nomimnal - (tmax_nomimnal - tmin_nomimnal) / 2.0;
	tmax = tmax_nomimnal + (tmax_nomimnal - tmin_nomimnal) / 2.0;

	// remove this location from all arrivals
	addr = head;
	do {
		if (addr->phase_time >= tmin) {
			addRemoveLocationInAssocLocationsList(addr, locID, 0);
		}
	} while ((addr = addr->next) != head  && addr->phase_time <= tmax); // !! assume sorted on phase_time

	return(0);

}


/** function to update associated location list for phases in a phase list */

int updateAssociatedLocation(PhsNode *head, Location *plocation, int locID, double tmin, double tmax)
{

	int i, nUpdated;
	double phaseTime;
	ArrivalDesc *parrival;
	PhsNode *addr;


	if ((addr = head) == NULL) {
		return(0);
	}

	// remove this location from all arrivals
	removeLocationAssociation(head, locID, tmin, tmax);


	// add this location to assocated arrivals
	nUpdated = 0;
	// loop over arrivals in associated location
	for (i = 0; i < plocation->narrivals; i++) {
		parrival = plocation->parrivals + i;
		phaseTime = getPhaseTimeValue(parrival);
		// loop over phases in phaselist within tmin->tmax
		addr = head;
		do {
			if (addr->phase_time >= tmin) {
				if (compareArrivals(addr->parrival, parrival, 0)
					&& compareTimeValue(getPhaseTimeValue(addr->parrival), phaseTime) == 0)
				{
					if (parrival->apriori_weight > VERY_SMALL_DOUBLE) {
						addRemoveLocationInAssocLocationsList(addr, locID, 1);
						nUpdated++;
					}
				}
			}
		} while ((addr = addr->next) != head  && addr->phase_time <= tmax); // !! assume sorted on phase_time
	}

	return(nUpdated);


}




/** function to find all phases with time in a specified time window
	if associatedFlag == -1   find unassociated
 	                      0   find all
  	                      N   find associated with location id=N (not yet implemented)
 */


PhsNode *findPhaseInTimeWindow(PhsNode *head, double tmin, double tmax, int associatedFlag)
{

	int nfound;
	PhsNode *addr;
	PhsNode *phsNodeList = NULL;


	if ((addr = head) == NULL) {
		return(NULL);
	}

	nfound = 0;
	do {
		if (addr->phase_time >= tmin && addr->phase_time <= tmax) {
			if (!associatedFlag || (associatedFlag == -1 && addr->passoc_locations[0] < 0)) {
				phsNodeList =
					addPhsNodeToPhaseList(phsNodeList, copyPhsNode(addr));
				nfound++;
			}
		}
	} while ((addr = addr->next) != head);

	if (nfound <= 0)
		return(NULL);

	return(phsNodeList);

}





/** function to find a phase */

PhsNode *findPhase(PhsNode *head, ArrivalDesc *arrivalKey)
{

	PhsNode *addr;
	double phaseTime;

	if ((addr = head) == NULL) {
		return(NULL);
	}

	phaseTime = getPhaseTimeValue(arrivalKey);

	do {
		if (compareArrivals(addr->parrival, arrivalKey, 0) &&
				  compareTimeValue(getPhaseTimeValue(addr->parrival), phaseTime) == 0)
			return(addr);
	} while ((addr = addr->next) != head);

	return(NULL);


}


/** function to find a phase */

#define CTA_VERBOSE 0

int compareArrivals(ArrivalDesc *parrival, ArrivalDesc *arrivalKey, int compareTimes)
{

	if (CTA_VERBOSE > 1) fprintf(stderr, "compareArrivals: 0\n");
	if (strcmp_to_null(parrival->label, arrivalKey->label) != 0)
		return(0);
	if (CTA_VERBOSE) fprintf(stderr, "compareArrivals: 1 %s %s\n", parrival->network, arrivalKey->network);
	if (strcmp_to_null(parrival->network, arrivalKey->network) != 0)
		return(0);
	if (CTA_VERBOSE) fprintf(stderr, "compareArrivals: 2 %s %s\n", parrival->inst, arrivalKey->inst);
	if (strcmp_to_null(parrival->inst, arrivalKey->inst) != 0)
		return(0);
	if (CTA_VERBOSE) fprintf(stderr, "compareArrivals: 3 %s %s\n",parrival->comp, arrivalKey->comp);
	if (strcmp_to_null(parrival->comp, arrivalKey->comp) != 0)
		return(0);
	if (CTA_VERBOSE) fprintf(stderr, "compareArrivals: 4 %s %s\n", parrival->phase, arrivalKey->phase);
	if (strcmp_to_null(parrival->phase, arrivalKey->phase) != 0)
		return(0);
	if (CTA_VERBOSE) fprintf(stderr, "compareArrivals: 5 %s %s\n", parrival->onset, arrivalKey->onset);
	if (strcmp_to_null(parrival->onset, arrivalKey->onset) != 0)
		return(0);
	if (CTA_VERBOSE) fprintf(stderr, "compareArrivals: 6 %s %s\n", parrival->first_mot, arrivalKey->first_mot);
	if (strcmp_to_null(parrival->first_mot, arrivalKey->first_mot) != 0)
		return(0);
	if (CTA_VERBOSE) fprintf(stderr, "compareArrivals: 7 %f %f\n", getPhaseTimeValue(parrival), getPhaseTimeValue(arrivalKey));
	if (compareTimes && compareTimeValue(getPhaseTimeValue(parrival), getPhaseTimeValue(arrivalKey)) != 0)
		return(0);
	if (CTA_VERBOSE) fprintf(stderr, "compareArrivals: 8 SUCCESS!\n");

	return(1);

}




/** function to compare two strings up to first '\0' char  */

int strcmp_to_null(char *s1, char *s2)
{
	int len, len1, len2;

	len1 = strlen(s1);
	len2 = strlen(s2);

	len = len1 < len2 ? len1 : len2;

	return(strncmp(s1, s2, len));

}









