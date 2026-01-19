#!/bin/bash

# 20231024 - Anthony Lomax, ALomax Scientific

# script to run NLLoc in parallel

# log output is written to run_nll_relocations.log


# USER VARIABLES TO EDIT =================================================================


# NLLoc control file settings ------------------------------------------------------------

# identifier for a particular run
RUN_NAME=20231212A

# existing NLLoc control file for locations
# IMPORTANT: comment out any INCLUDE statements not needed for all configurations!
CONTROL_FILE=my_nll_control_file.in

# output file root name
PROJECT_NAME=MyProjectName_2023

# existing travel-time grid files path/root
TIME_ROOT=path_to_time_grids/root_name
SWAP_BYTES=0	# Grid2Time travel-time files
#SWAP_BYTES=1	# TauPToolkit travel-time files

# observations settings
OBS_FILE_TYPE=NLLOC_OBS
LOC_OBS="obs/*.hyp"

# set different configurations for NLLoc
if [ YES == YES ]; then
	# initial locations
	mv run_nll_relocations.log run_ssst_relocations_OLD.log   # clear log file
	LOCCOM="My Locations 2023"
	OUT_ROOT="out/${RUN_NAME}/loc"
	INCLUDE=""
else
	# simple static station corrections
	LOCCOM="My Locations 2023 corr1"
	OUT_ROOT="out/${RUN_NAME}/loc_corr1"
	INCLUDE="INCLUDE out/${RUN_NAME}/loc/${PROJECT_NAME}.sum.grid0.loc.stat_totcorr"
fi

# specify the number NLLoc to be run in parallel (e.g. up to the number of physical or virtual CPU cores available)
NUM_CORES=7


# less important NLLoc control file settings ------------------------------------------------------------



# misc less important settings ------------------------------------------------------------

# visualisation command
SV_CMD="java net.alomax.seismicity.Seismicity"


# END - USER VARIABLES TO EDIT =================================================================




# NLLoc control file for parallel runs
SKELETON_CONF=tmp/NLL.cluster_SKELETON_SSST.conf
# comment out LOCFILES in original loc control file and output to SKELETON_CONF
sed '/LOCFILES/ s/^#*/#/' ${CONTROL_FILE} > ${SKELETON_CONF}

# ----------------------------------------------------------------
# run NLLoc in parallel

mkdir tmp
mkdir -p ${OUT_ROOT}

cat << END >> run_nll_relocations.log
Running NLLoc:
${SV_CMD} ${OUT_ROOT}/${PROJECT_NAME}*sum*hyp &
END

CONTROL_FILE_TMP=tmp/${PROJECT_NAME}_nll.in
cp ${SKELETON_CONF} ${CONTROL_FILE_TMP}

# get a list of all obs files
OBS_LIST=$(echo ${LOC_OBS})
for ENTRY in ${OBS_LIST}; do
	echo "${ENTRY}"
done  > tmp/obs.txt
# count the total number of obs files
COUNT=$(wc -l < tmp/obs.txt) 
rm -r tmp/obs_*
rm -r tmp/obsfiles_*
# split obs file list into NUM_CORES sub-lists
echo "split -l $((1 + ${COUNT} / ${NUM_CORES})) tmp/obs.txt tmp/obs_"
split -l $((1 + ${COUNT} / ${NUM_CORES})) tmp/obs.txt tmp/obs_
# produces temp/obs_aa, temp/obs_ab, etc
# run NLLoc for each sub-list
declare -i INDEX=0
for SPLIT_OBS_FILE in tmp/obs_* ; do
	echo "Running: ${INDEX} ${SPLIT_OBS_FILE}"
	cp ${CONTROL_FILE_TMP} ${CONTROL_FILE_TMP}_${INDEX}
	# copy obs files in sub-list to temp obs file directory
	mkdir tmp/obsfiles_${INDEX}
	{
		while read OFILE; do
			cp -p ${OFILE} tmp/obsfiles_${INDEX}
		done
	} < ${SPLIT_OBS_FILE}
	cat << END >> ${CONTROL_FILE_TMP}_${INDEX}
LOCCOM ${PROJECT_NAME} ${RUN_NAME} ${LOCCOM}
LOCFILES tmp/obsfiles_${INDEX}/* ${OBS_FILE_TYPE}  ${TIME_ROOT}  ${OUT_ROOT}_${INDEX}/${PROJECT_NAME}  ${SWAP_BYTES}
${INCLUDE}
END
	mkdir -p ${OUT_ROOT}_${INDEX}
	NLLoc ${CONTROL_FILE_TMP}_${INDEX} &
	PIDS[${INDEX}]=$!
	INDEX=INDEX+1
done
# wait for all PIDS
for PID in ${PIDS[*]}; do
	wait $PID
	status=$?
	echo "Finished: PID=${PID} status=${status} ================================="
done

# assemble unique NLL output sum files
cp -a ${OUT_ROOT}_*/. ${OUT_ROOT}/
cat ${OUT_ROOT}_*/${PROJECT_NAME}.sum.grid0.loc.hyp > ${OUT_ROOT}/${PROJECT_NAME}.sum.grid0.loc.hyp
echo "" > ${OUT_ROOT}/WARNING.concatenated_output_of_multiple_NLLoc_runs.WARNING
cat ${OUT_ROOT}_*/${PROJECT_NAME}.sum.grid0.loc.stations > ${OUT_ROOT}/${PROJECT_NAME}.sum.grid0.loc.stations
cat ${OUT_ROOT}_*/${PROJECT_NAME}.sum.grid0.loc.stat > ${OUT_ROOT}/${PROJECT_NAME}.sum.grid0.loc.stat
cat ${OUT_ROOT}_*/${PROJECT_NAME}.sum.grid0.loc.stat_totcorr > ${OUT_ROOT}/${PROJECT_NAME}.sum.grid0.loc.stat_totcorr

# cleanup
rm -r ${OUT_ROOT}_?
rm -r ${OUT_ROOT}_??

cat << END >> nll_parallel.list
${OUT_ROOT}/${PROJECT_NAME}.sum.grid0.loc.hyp
END

echo "Done!"
echo "${OUT_ROOT}/${PROJECT_NAME}.sum.grid0.loc.hyp"

# END run NLLoc in parallel
# ----------------------------------------------------------------
