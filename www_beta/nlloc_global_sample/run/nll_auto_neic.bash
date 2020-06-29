# set volatile variables
NUM_OBS_MIN=5
NUM_OBS_MAX=100
VISUALIZE=YES
ITERATION_METHOD=EDT
RUN_ROOT=/home/anthony/nlloc_tmp/neic/loc

# set quasi-static variables
# general ===========================
declare -i ITERATION_MAX=5
# NLLoc ===========================
NLL_CONTROL_FILE=nll.in
#OBS_FILE=./obs/20040224.0227_MOROCCO.neic
OBS_FILE=/windows/C/download/seisdata/neic/obs/20040317_crete.html
OBS_FILE_FORMAT=NEIC
TT_GRID_FILES=/home/anthony/nlloc_tmp/taup/ak135/ak135
#TT_GRID_FILES=/home/anthony/nlloc_tmp/taup/global/iasp91
#RUN_ROOT=${RUN_ROOT}_IPO
#
# global
LOC_SEARCH="LOCSEARCH  OCT 36 18 10 0.1 50000 5000"
LOC_GRID="LOCGRID  361 181 1001  -180.0 -90.0 0.0  1.0 1.0 1.0   PROB_DENSITY  SAVE"
#
# PhsAssoc ===========================
declare -i MAX_RESIDUAL0=10
declare -i MAX_RESIDUAL1=5

I_SWAP_BYTES=1
I_LOOK_FOR_STATION_GRIDS=0



# functions =============================================================================

# function to write skeleton NLL control file
write_nll_skeleton ()
{
cat > $1 << END
# =============================================================================
#  Automatically generated NonLinLoc programs control file
#  $(date)
# =============================================================================
END
cat >> $1 << END
# see documentation http://alomax.free.fr/nlloc/soft2.37/control.html#_NLLoc_
CONTROL 1 54321
TRANS GLOBAL
INCLUDE run/sta_list_neic.in
LOCSIG Anthony Lomax
LOCCOM NEIC eq locations (IASPEI91 model)
LOCHYPOUT SAVE_NLLOC_ALL
${LOC_SEARCH}
${LOC_GRID}
LOCGAU 1.0 5.0
LOCPHASEID  P   P p Pn Pdiff PKP PKiKP PKIKP
LOCPHASEID  S   S s Sn Sdiff SKS SKiKS SKIKS
LOCQUAL2ERR 1.0 2.0 4.0 8.0 99999.9
LOCPHSTAT 9999.0 -1 9999.0 1.0 1.0
LOCANGLES ANGLES_NO 5
LOCMAG MD_FMAG -1.465 2.22 0.0 0.0 0.0
LOCALIAS PAD_ X_PAD   0 0 0   9999 99 99
END
}

# function to run SeismicityViewer
run_seismicity ()
{
echo
echo java net.alomax.seismicity.Seismicity  $1/$2.*.*.grid0.loc.hyp
java -mx100m net.alomax.seismicity.Seismicity  $1/$2.*.*.grid0.loc.hyp > /dev/null &
}

# END - functions =============================================================================



# main script =============================================================================

# clean
rm -r ${RUN_ROOT}
mkdir ${RUN_ROOT}
cp $0 ${RUN_ROOT}

if [ ${ITERATION_METHOD} = EDT ]; then
	# run initial RMS locations
	RUN_ID_OUT=RMS_0
	# write NLL control file skeleton
	write_nll_skeleton ${NLL_CONTROL_FILE}
	# write run specific NLL control file statements
cat >> ${NLL_CONTROL_FILE} << END
LOCFILES ${OBS_FILE} ${OBS_FILE_FORMAT}  ${TT_GRID_FILES}  ${RUN_ROOT}/${RUN_ID_OUT} 1
END
cat >> ${NLL_CONTROL_FILE} << END
LOCMETH GAU_ANALYTIC 1.0e6 ${NUM_OBS_MIN} ${NUM_OBS_MAX} -1 -1.80 6
END
	# run NLL
	nice NLLoc ${NLL_CONTROL_FILE}
	# run SeismicityViewer
	if [ ${VISUALIZE} = YES ]; then
		run_seismicity ${RUN_ROOT} ${RUN_ID_OUT}
	fi
fi


# run initial locations
if [ ${ITERATION_METHOD} = EDT ]; then
	RUN_ID_OUT=EDT_0
else
	RUN_ID_OUT=RMS_0
fi
# write NLL control file skeleton
write_nll_skeleton ${NLL_CONTROL_FILE}
# write run specific NLL control file statements
cat >> ${NLL_CONTROL_FILE} << END
LOCFILES ${OBS_FILE} ${OBS_FILE_FORMAT}  ${TT_GRID_FILES}  ${RUN_ROOT}/${RUN_ID_OUT} 1
END
if [ ${ITERATION_METHOD} = EDT ]; then
cat >> ${NLL_CONTROL_FILE} << END
LOCMETH EDT_TEST 1.0e6 ${NUM_OBS_MIN} ${NUM_OBS_MAX} -1 -1.80 6
END
else
cat >> ${NLL_CONTROL_FILE} << END
LOCMETH GAU_ANALYTIC 1.0e6 ${NUM_OBS_MIN} ${NUM_OBS_MAX} -1 -1.80 6
END
fi
# run NLL
nice NLLoc ${NLL_CONTROL_FILE}
# run SeismicityViewer
if [ ${VISUALIZE} = YES ]; then
	run_seismicity ${RUN_ROOT} ${RUN_ID_OUT}
fi


# iterate PhsAssoc then re-locate

CHECK_AGAIN=YES
NUM_OBS_MIN_ITER=4
I_WRITE_ALL=1

declare -i I_COUNT=1
declare -i MAX_RESIDUAL=${MAX_RESIDUAL0}

while [ $((${I_COUNT} <= ${ITERATION_MAX})) = 1 ] && [ ${CHECK_AGAIN} = YES ]; do

	echo
	echo "--- Iteration: " ${I_COUNT}/${ITERATION_MAX}, from ${RUN_ID_IN} to ${RUN_ID_OUT}

if [ ${ITERATION_METHOD} = EDT ]; then
	RUN_ID_IN=${RUN_ID_OUT}
	RUN_ID_OUT=EDT_${I_COUNT}
else
	RUN_ID_IN=${RUN_ID_OUT}
	RUN_ID_OUT=RMS_${I_COUNT}
fi

	# associate phases
	if [ $((${I_COUNT} > 1)) = 1 ]; then
		rm ${RUN_ROOT}/${RUN_ID_IN}.*.*.grid0.loc.phs_assoc
	fi
	MAX_RESIDUAL=$((${MAX_RESIDUAL1}+(${MAX_RESIDUAL0}-${MAX_RESIDUAL1})*(${ITERATION_MAX}-${I_COUNT})/(${ITERATION_MAX}-1)))
	echo MAX_RESIDUAL=${MAX_RESIDUAL}
	if PhsAssoc ${MAX_RESIDUAL} ${I_SWAP_BYTES} ${I_LOOK_FOR_STATION_GRIDS} ${I_WRITE_ALL} \
		${TT_GRID_FILES} ${RUN_ROOT}/${RUN_ID_IN}.*.*.grid0.loc.hyp
	then
		CHECK_AGAIN=NO
	else
		CHECK_AGAIN=YES
#		if [ $((${I_COUNT} > 1)) = 1 ]; then
#			rm ${RUN_ROOT}/${RUN_ID_IN}.*.*.grid0.loc.hyp
#		fi
		# write NLL control file skeleton
		write_nll_skeleton ${NLL_CONTROL_FILE}
		# write run specific NLL control file statements
cat >> ${NLL_CONTROL_FILE} << END
LOCFILES ${RUN_ROOT}/${RUN_ID_IN}.*.*.grid0.loc.phs_assoc NLLOC_OBS  ${TT_GRID_FILES}  \
	${RUN_ROOT}/${RUN_ID_OUT} 1
END
if [ ${ITERATION_METHOD} = EDT ]; then
cat >> ${NLL_CONTROL_FILE} << END
LOCMETH EDT_TEST 1.0e6 ${NUM_OBS_MIN} ${NUM_OBS_MAX} -1 -1.80 6
END
else
cat >> ${NLL_CONTROL_FILE} << END
LOCMETH GAU_ANALYTIC 1.0e6 ${NUM_OBS_MIN} ${NUM_OBS_MAX} -1 -1.80 6
END
fi
		# run NLL
		nice NLLoc ${NLL_CONTROL_FILE}
	fi

	I_COUNT=I_COUNT+1

	I_WRITE_ALL=0

done

# run SeismicityViewer
if [ ${VISUALIZE} = YES ]; then
	run_seismicity ${RUN_ROOT} ${RUN_ID_OUT}
fi


# run final RMS locations
if [ ${ITERATION_METHOD} = EDT ]; then
	RUN_ID_IN=${RUN_ID_OUT}
	RUN_ID_OUT=RMS_N
fi
# write NLL control file skeleton
write_nll_skeleton ${NLL_CONTROL_FILE}
# write run specific NLL control file statements
cat >> ${NLL_CONTROL_FILE} << END
LOCFILES ${RUN_ROOT}/${RUN_ID_IN}.*.*.grid0.loc.hyp NLLOC_OBS  ${TT_GRID_FILES}  \
	${RUN_ROOT}/${RUN_ID_OUT} 1
END
cat >> ${NLL_CONTROL_FILE} << END
LOCMETH GAU_ANALYTIC 1.0e6 ${NUM_OBS_MIN} ${NUM_OBS_MAX} -1 -1.80 6
END
# run NLL
nice NLLoc ${NLL_CONTROL_FILE}
# run SeismicityViewer
if [ ${VISUALIZE} = YES ]; then
	run_seismicity ${RUN_ROOT} ${RUN_ID_OUT}
fi



# END - main script =============================================================================
