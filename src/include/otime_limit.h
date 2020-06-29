/* 
 * File:   otime_limit.h
 * Author: anthony
 *
 * Created on 14 September 2009, 15:52
 */

#ifndef _OTIME_LIMIT_H
#define	_OTIME_LIMIT_H

#ifdef	__cplusplus
extern "C" {
#endif


typedef struct Otime_Limit
{
	int data_id;
	double time;		// time of limit
	double otime;		// estimated origin time
	int polarity;		// = +1 if start time of otime limit, -1 if end time
        //struct Otime_Limit* pair; // paired limit
	double dist_range;	// distance range (effective cell size) estimated from (otime_end - otime_start) / slowness
	double time_range;	// time range (otime_end - otime_start)

}
OtimeLimit;

OtimeLimit* new_OtimeLimit(int data_id, double time, double otime, int polarity, double dist_range, double time_range) ;
void addOtimeLimitToList(OtimeLimit* otimeLimit, OtimeLimit*** potime_limit_list, int* pnum_otime_limit);
void free_OtimeLimitList(OtimeLimit*** otime_limit_list, int* pnum_otime_limit);
void free_OtimeLimit(OtimeLimit* otime_limit);



#ifdef	__cplusplus
}
#endif

#endif	/* _OTIME_LIMIT_H */

