#!/bin/bash

# 20241022 - Anthony Lomax, ALomax Scientific

# script to run GridCascadingDecimate

# For a description of GridCascadingDecimate, see section "3.2. Hypocenter location method and parameterization strategy" in:
# Latorre, D., Di Stefano, R., Castello, B., Michele, M., & Chiaraluce, L. (2023). An updated view of the Italian seismicity from probabilistic location in 3D velocity models: The 1981–2018 Italian catalog of absolute earthquake locations (CLASS). Tectonophysics, 846, 229664. https://doi.org/10.1016/j.tecto.2022.229664
#
# In order to use the Italian 3D Tomographic Model for FD traveltime computation with the P &L code, we have resampled both P and S tomographic grids into finer grids of cubic cells with mesh spacing of 2 km and constant slowness inside each cell. The topography is not included in our parameterization, but the velocity model grid is extended up to the cells that include the seismic stations. By following the NLL procedure, P and S traveltime grids need to be computed for each seismic station and stored before solving the inversion problem. The main issue concerning earthquake locations at the Italian scale based on 3D, FD traveltime computation is the huge dimension of the traveltime grids involved in the inversion procedure, since the Italian seismicity is distributed in a target volume of about 1200 × 1200 × 800 km3 and recorded at more than 400 stations.Therefore, using cubic cells having 2 km × 2 km × 2 km size, we should manage hundreds of very large grids (up to 2,88 × 108 nodes for each station) making the hypocenter location of the entire catalog excessively time-consuming. To tackle this problem, we first calculate the traveltimes with the P&L code in the original finer grids of cubic cells with mesh spacing of 2 km to maintain precision, and then we decimatethe computed traveltime grids for increasing depths by preserving the finer, 2 km-size grids in the upper crust (from the Earth’s surface to 20 km depth) where we need a more detailed representation of the 3D structure, and doubling incrementally the cell size at specified depths, in the lower crust (4 km-size grid from 20 to 40 km depth),at lithospheric depths (8 km-size grid from 40 to 100 km depth),and in the uppermantle (16 km-size grid from 100 to 600 km depth). The conversion of the traveltime grids has been performed with the cascading grid module included in the NLL software (http://alomax.free.fr/nlloc). Using our parametrization, we estimated that typical errors due to the conversion from full regular grids to cascading grids are in a range of ±0.1 s, lower than average picking errors estimated on P and S waves input data. On the contrary, resampled grids allow us to reduce by ~90% the size of each station-phase traveltime grid, making the earthquake location of the entire catalog easier to manage.




# USER VARIABLES TO EDIT =================================================================

MODEL_GRID_ROOT=out/model/BIG_3D_MODEL_NAME
TIME_GRIDS_PATH=out/time
GRID_NAME=BIG_3D_MODEL_NAME

STATION_SRCE_FILE=my_stations_GTSRCE.in

TRANS="TRANS  TRANS_MERC WGS-84  0.0 15.0  0.0"

# GridCascadingDecimate Arguments: <GridCascadingDecimate> 
# Usage: GridCascadingDecimate <doubling_depths> <input grid(s)> <output grid path> [<flag_use_mean_value_in_casc_cell>]
#    doubling_depths - comma separated, increasing (approx) depths to double cell size
DEPTHS=20,40,100  # doubling_depths

PHASE=P
#PHASE=S

DELETE_ORIGINAL_TIME_GRIDS=NO
#DELETE_ORIGINAL_TIME_GRIDS=YES  # delete large, original time grids

# END - USER VARIABLES TO EDIT =================================================================




echo "Phase is ${PHASE}"
mkdir ${TIME_GRIDS_PATH}
mkdir ${TIME_GRIDS_PATH}_casc

{
	# GTSRCE ATN LATLON 38.15950 15.46470 0.0 1.130
	while read GTSRCE_KEY STATION LATLON_KEY LAT LON DEPTH ELEV ; do
	
		echo "===================================================================="
		echo "Processing: ${GTSRCE_KEY} ${STATION} ${LATLON_KEY} ${LAT} ${LON} ${DEPTH} ${ELEV}"

		# check if casc grid already generated
		if [ -f ${TIME_GRIDS_PATH}_casc/casc.${GRID_NAME}.${PHASE}.${STATION}.time.buf ] ; then
			echo "Cascading time grid exists, not regenerated: ${TIME_GRIDS_PATH}_casc/casc.${GRID_NAME}.${PHASE}.${STATION}.time.buf"
		else

			# check if time grid already generated
			if [ -f ${TIME_GRIDS_PATH}/${GRID_NAME}.${PHASE}.${STATION}.time.buf ] ; then
				echo "Original time grid exists, not regenerated: ${TIME_GRIDS_PATH}/${GRID_NAME}.${PHASE}.${STATION}.time.buf"
			else
				cat > nll_temp.in << END
CONTROL 1 54321
${TRANS}
GTMODE GRID3D ANGLES_NO
GT_PLFD  1.0e-3  0
${GTSRCE_KEY} ${STATION} ${LATLON_KEY} ${LAT} ${LON} ${DEPTH} ${ELEV}
GTFILES  ${MODEL_GRID_ROOT}  ${TIME_GRIDS_PATH}/${GRID_NAME}   ${PHASE}
END
				COMMAND="nice Grid2Time nll_temp.in"
				echo ${COMMAND}
				time ${COMMAND}
			fi
		
			COMMAND="nice GridCascadingDecimate ${DEPTHS} ${TIME_GRIDS_PATH}/${GRID_NAME}.${PHASE}.${STATION}.time.buf ${TIME_GRIDS_PATH}_casc/"
			echo ${COMMAND}
			time ${COMMAND}
		
			if [ ${DELETE_ORIGINAL_TIME_GRIDS} == YES ]; then
				# delete large, original time grids
				COMMAND="rm ${TIME_GRIDS_PATH}/${GRID_NAME}.${PHASE}.${STATION}.time.*"
				echo ${COMMAND}
				{COMMAND}
			fi
			
		fi
	
		echo ""

	done

} < ${STATION_SRCE_FILE}


