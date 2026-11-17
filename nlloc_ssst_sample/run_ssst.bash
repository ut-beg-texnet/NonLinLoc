# 20210414 - Anthony Lomax, ALomax Scientific

# script to run each iteration of SSST with multiple instances of Loc2SSST run in parallel

CONTROL_FILE=$1
 
 
declare -i INDEX=0

# NOTE: TYPICALY, ONLY EDIT HERE:
# specify N station sets, where N is the number of Loc2SSST instances to be run in parallel (e.g. the number of pyhsical CPU cores available)
# 	separate each station set by a blank space " "
# 	separate stations in each set with the equals sign "="
# see stations.xlxs to get station codes
#
if [ NO = YES ]; then
	# to randomize order:
	# orig.txt contains one NET_STA codes per line (e.g. column of stations.ods)
	cat orig.txt | sort -R > tmp/rand.txt
	#
	# station list into NUM_CORES station sets:
	COUNT=$(wc -l < tmp/rand.txt) 
	NUM_CORES=9
	rm tmp/rand_sta_??; split -l $((1 + ${COUNT} / ${NUM_CORES})) tmp/rand.txt tmp/rand_sta_
	rm tmp/rand_sta_ssst.txt; for SFILE in tmp/rand_sta_??; do cat ${SFILE} >> tmp/rand_sta_ssst.txt; echo " " >> tmp/rand_sta_ssst.txt; done;
# replace \n with = in tmp/rand_sta_ssst.txt, remove extra, empty station set at end
fi
#
for STA_SET in \
NC_PPT=NC_PSC=NC_PWK=NC_PMM=NC_PAG= =NC_PSM=NC_PADB=BP_JCNB=NC_PHP=BP_MMNB= =NC_PCA=NC_PMPB=NC_PHF=NC_PL11=NC_PFR= =NC_PSR=BP_SCYB=NC_PMR=NC_PVC=NC_PCHC= =NC_PJC=BP_LCCB=BP_CCRB=NC_PST=NC_PHA= =NC_PPO=BP_SMNB=NC_PKE=NC_PTR=NC_PGH= =NC_PDR=NC_PSA=NC_PMC=NC_PAR=NC_PHOB= =BP_VCAB=NC_PAN=NC_PBS=
do
	
	echo "Running: ${INDEX} ${STA_SET}"
	
	cp ${CONTROL_FILE} tmp/ssst_${INDEX}.in
	cat << END >> tmp/ssst_${INDEX}.in
LSSTATIONS ${STA_SET}
END
	Loc2ssst tmp/ssst_${INDEX}.in &
    PIDS[${INDEX}]=$!
    	
	INDEX=INDEX+1
	
done

# wait for all PIDS
for PID in ${PIDS[*]}; do
    wait $PID
	status=$?
	echo "Finished: PID=${PID} status=${status} ================================="
done
