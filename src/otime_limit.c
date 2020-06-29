#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>

#include "otime_limit.h"





/** otime limit class */

OtimeLimit* new_OtimeLimit(int data_id, double time, double otime, int polarity, double dist_range, double time_range) {

	OtimeLimit* otime_limit = calloc(1, sizeof(OtimeLimit));

	otime_limit->data_id = data_id;
	otime_limit->time = time;
	otime_limit->otime = otime;
	otime_limit->polarity = polarity;
	//otime_limit->pair = (OtimeLimit*) NULL;
	otime_limit->dist_range = dist_range;
	otime_limit->time_range = time_range;

	return(otime_limit);

}


/** add an origin time limit to an OtimeLimit list
 * list will be sorted by increasing otime->time
 */

#define SIZE_INCREMENT_OTIME_LIMIT_LIST 128

void addOtimeLimitToList(OtimeLimit* otimeLimit, OtimeLimit*** potime_limit_list, int* pnum_otime_limit) {


	OtimeLimit** newOtimeLimitList = NULL;

	if (*potime_limit_list == NULL) {		// list not yet created
		*potime_limit_list = calloc(SIZE_INCREMENT_OTIME_LIMIT_LIST, sizeof(OtimeLimit*));
		*pnum_otime_limit = 0;
	}
	else if (*pnum_otime_limit != 0 && (*pnum_otime_limit % SIZE_INCREMENT_OTIME_LIMIT_LIST) == 0) {	// list will be too small
		newOtimeLimitList = calloc(*pnum_otime_limit + SIZE_INCREMENT_OTIME_LIMIT_LIST, sizeof(OtimeLimit*));
		int n;
		for (n = 0; n < *pnum_otime_limit; n++)
			newOtimeLimitList[n] = (*potime_limit_list)[n];
		free(*potime_limit_list);
		*potime_limit_list = newOtimeLimitList;
	}

	// find first limit later than new limit
	int ninsert;
	for (ninsert = 0; ninsert < *pnum_otime_limit; ninsert++)
		if ((*potime_limit_list)[ninsert]->time > otimeLimit->time)
			break;
	// shift later data
	if (ninsert < *pnum_otime_limit) {
		int m;
		for (m = *pnum_otime_limit - 1; m >= ninsert; m--)
			(*potime_limit_list)[m + 1] = (*potime_limit_list)[m];
	}
	// insert new OtimeLimit
	(*potime_limit_list)[ninsert] = otimeLimit;
	(*pnum_otime_limit)++;

}


/** clean up list memory */

void free_OtimeLimitList(OtimeLimit*** potime_limit_list, int* pnum_otime_limit)
{
	if (*potime_limit_list == NULL)
		return;

	int n;
	for (n = 0; n < *pnum_otime_limit; n++)
		free_OtimeLimit(*(*potime_limit_list + n));

	free(*potime_limit_list);
	*potime_limit_list = NULL;

	*pnum_otime_limit = 0;

}




/** clean up data memory */

void free_OtimeLimit(OtimeLimit* otime_limit)
{
	if (otime_limit == NULL)
		return;

	free(otime_limit);
	otime_limit = NULL;

}



