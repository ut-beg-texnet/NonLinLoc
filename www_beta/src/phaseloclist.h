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


/*  phaseloclist.h

	include file for phaselist and loclist

*/



/*-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
161 Allee du Micocoulier, 06370 Mouans-Sartoux, France
tel: +33(0)493752502  e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>

#ifdef EXTERN_MODE
#define	EXTERN_TXT extern
#else
#define EXTERN_TXT
#endif



/* misc defines */

#ifndef SMALL_DOUBLE
#define SMALL_DOUBLE 1.0e-20
#endif
#ifndef LARGE_DOUBLE
#define LARGE_DOUBLE 1.0e20
#endif
#ifndef VERY_SMALL_DOUBLE
#define VERY_SMALL_DOUBLE 1.0e-30
#endif
#ifndef VERY_LARGE_DOUBLE
#define VERY_LARGE_DOUBLE 1.0e30
#endif


#define INVALID_DOUBLE -VERY_LARGE_DOUBLE

#define INIT_NUM_ASSOC_LOC_PER_PHASE 8

#define MAX_NUM_LOCATIONS 1000
#define MAX_NUM_PHASES_PER_LOC 4096

#define DUPLICATE_PAHSE_FOUND ((PhsNode *) -9876)



/*------------------------------------------------------------/ */
/* structures */
/*------------------------------------------------------------/ */


/* phaselist */

/* phaselist node */

typedef struct phasenode
{
	struct phasenode *prev;	/* pointer to previous vertex */
	struct phasenode *next;	/* pointer to next vertex */
	int id;			/* vertex identification */
	double phase_time;	/* phase time in seconds - sort value */
	ArrivalDesc* parrival;	/* phase arrival data (observation part only is initialized) */
	int *passoc_locations;	/* id's of associated locations, -1 if none */
	int passoc_locations_size;	/* size of passoc_locations array */
} PhsNode;



/* loclist */

/* location structure */

typedef struct location
{
	HypoDesc* phypo;	/* location data */
	ArrivalDesc* parrivals;	/* phase arrival data */
		/* WARNING: for nll_assoc this data is temporary, definitive associated arrivals specified by PhsNode.passoc_locations location id list!!! */
	int narrivals;		/* num of phase arrival data */
	GridDesc* pgrid;	/* location grid */
	Tree3D* poctTree;	/* location oct-tree structure */
	float* pscatterSample;	/* location scatter sample data */

} Location;


/* loclist node structure*/

typedef struct locnode
{
	struct locnode *prev;	/* pointer to previous vertex */
	struct locnode *next;	/* pointer to next vertex */
	int id;			/* vertex identification - sort value */
	Location *plocation;	/* location data */
	double first_phase_time;	/* time value for earliest phase */
	// NOT NEEDED? PhsNode *phsNodeList;	/* list of phase nodes used for location */
} LocNode;


/* */
/*------------------------------------------------------------/ */



/*------------------------------------------------------------/ */
/* globals  */
/*------------------------------------------------------------/ */

EXTERN_TXT ArrivalDesc arrivals_tmp[MAX_NUM_PHASES_PER_LOC];

/* */
/*------------------------------------------------------------/ */




/*------------------------------------------------------------/ */
/* function declarations */
/*------------------------------------------------------------/ */

/* phaselist */

PhsNode *addArrivalToPhaseList(PhsNode **phead, ArrivalDesc* parrival, int id, int addDuplicates);
PhsNode *addPhsNodeToPhaseList(PhsNode *phead, PhsNode* addr);
PhsNode *removeArrivalFromPhaseList(PhsNode *head, PhsNode* addr, int freeArrivalDesc);
int freePhaseList(PhsNode *head, int freeArrivalDesc);
PhsNode *addNLLPhaseStringToPhaseList(PhsNode **phead, char *phase_string, int id, int addDuplicates);
int writePhaseList(PhsNode *head, FILE *out);
double getPhaseTimeValue(ArrivalDesc *parrival);
int compareTimeValue(double t1, double t2);
int addRemoveLocationInAssocLocationsList(PhsNode *addr, int locID, int addLocID);
int updateAssociatedLocation(PhsNode *head, Location *plocation, int locID, double tmin, double tmax);
PhsNode *findPhaseInTimeWindow(PhsNode *head, double tmin, double tmax, int associatedFlag);
PhsNode *findPhase(PhsNode *head, ArrivalDesc *arrivalKey);
int compareArrivals(ArrivalDesc *parrival, ArrivalDesc *arrivalKey, int compareTimes);
int strcmp_to_null(char *s1, char *s2);
// AJL 20070323-
PhsNode *getPhsNodeFromPhaseList(PhsNode *head, int id);
int removeLocationAssociation(PhsNode *head, int locID, double tmin_nomimnal, double tmax_nomimnal);


/* loclist */

Location *newLocation(HypoDesc *phypo, ArrivalDesc* parrivals, int narrivals, GridDesc *pgrid, Tree3D* poctTree, float *pscatterSample);
LocNode *addLocationToLocList(LocNode **phead, Location *plocation, int id);
LocNode *removeLocationFromLocList(LocNode *head, LocNode* addr, int freeHypoDesc);
int freeLocList(LocNode *head, int freeHypoDesc);
void freeLocation(Location *plocation);
//LocNode *addNLLLocationFileToLocList(LocNode **phead, char *filein, int id, PhsNode **phsNodeList);
int writeLocList(LocNode *head, FILE *out);
void writeLocNode(LocNode *addr, FILE *out);
double getLocTimeValue(HypoDesc *phypo);
LocNode **findLocsWithFirstPhaseInTimeWindow(LocNode *LocListHead, double tmin, double tmax);
ArrivalDesc *findArrivalInLocation(Location *plocation, ArrivalDesc *arrivalKey);
int arrivalBeforeLastLocateArrival(ArrivalDesc *arrivalKey, Location *plocation);
// AJL 20070323-
LocNode *getLocationFromLocList(LocNode *head, int id);


/* */
/*------------------------------------------------------------/ */


