#  --- GMT defaults
gmtset  PAPER_MEDIA A4   PAGE_ORIENTATION landscape   MEASURE_UNIT cm    ANOT_FONT_SIZE 10p    HEADER_FONT_SIZE 18p    LABEL_FONT_SIZE 12p    DEGREE_FORMAT 5


#########
PHASE=P
#########

echo PHASE: ${PHASE}

MODEL=ak135
echo MODEL: ${MODEL}
OUTPATH=./${MODEL}
echo OUTPATH: ${OUTPATH}

DATAFILE0=${OUTPATH}/${MODEL}.${PHASE}.DEFAULT.time
PSFILE=${DATAFILE0}.ps
DATAFILE=${DATAFILE0}.buf

CPTFILE=${DATAFILE0}.cpt
makecpt -T0/2500/50 > ${CPTFILE}

STEP_SIZE=0.1/2.5
RVAL=-R0/180/5678/6378
xyz2grd  ${DATAFILE} -G${DATAFILE}.grd -I${STEP_SIZE} -ZfLTw ${RVAL}

JVAL=-Jpa0.002/90
grdimage  ${DATAFILE}.grd -C${CPTFILE} ${RVAL} ${JVAL}  -K -X2 -Y5 > ${PSFILE}
grdcontour  ${DATAFILE}.grd -C${CPTFILE} ${RVAL} ${JVAL}   -A- -O -K >> ${PSFILE}
psbasemap  ${RVAL} ${JVAL}  -Ba10:Distance:/a200:Depth::.${DATAFILE0}:WESN -O -K >> ${PSFILE}

#  --- scale
psscale -D12.75/-2/25/1.0h -C${CPTFILE} -B100:"Travel Time":/:sec: -X0 -Y0  -O >> ${PSFILE}
