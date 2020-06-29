echo "============"
echo "Script to create spherical, layered-model, travel-time grids for NonLinLoc - Global mode"
echo "see http://alomax.net/nlloc"
echo
echo "IMPORTANT:  Requires:"
echo "   1. Java - the command \"java\" must be on your path"
echo "   2. TauP Toolkit classes must be on your java classpath- see: http://www.seis.sc.edu/software/TauP"
echo "   3. The additional class edu.sc.seis.TauP.TauP_Table_NLL must be on your java classpath"
echo "     - this class is available at http://alomax.net/nlloc/java"

echo ============ whole earth

MODEL=ak135
echo MODEL: ${MODEL}
OUTPATH=./${MODEL}
echo OUTPATH: ${OUTPATH}

mkdir ${OUTPATH}

java edu.sc.seis.TauP.TauP_Create -tvel ${MODEL}.tvel

# NLL_GRID=<num cells in distance>,<distance min in deg>,<distance max in deg>,<num cells in depth>,<depth min in km>,<depth max in km>
#  Set the parameters with regards to the area (cylindrical slice of the spherical earth) that you want to study.  The default is for the whole earth.  The grid is 2D, so it is not critical if you make it larger than needed.
NLL_GRID=1801,0.0,180.0,281,0.0,700.0
echo $NLL_GRID

echo First arriving P and S

java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph P,p,Pn,Pdiff,PKP,PKiKP,PKIKP -o ${OUTPATH}/${MODEL}
java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph S,s,Sn,Sdiff,SKS,SKiKS,SKIKS -o ${OUTPATH}/${MODEL}


echo Pg and Sg

#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph Pg -o ${OUTPATH}/${MODEL}
#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph Sg -o ${OUTPATH}/${MODEL}


echo pP and sP

#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph pP -o ${OUTPATH}/${MODEL}
#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph sP -o ${OUTPATH}/${MODEL}


echo PcP and ScP and PcS and ScS, etc

#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph PcP -o ${OUTPATH}/${MODEL}
#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph pPcP -o ${OUTPATH}/${MODEL}
#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph sPcP -o ${OUTPATH}/${MODEL}
#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph ScP -o ${OUTPATH}/${MODEL}
#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph pScP -o ${OUTPATH}/${MODEL}
#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph sScP -o ${OUTPATH}/${MODEL}
#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph PcS -o ${OUTPATH}/${MODEL}
#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph pPcS  -o ${OUTPATH}/${MODEL}
#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph sPcS  -o ${OUTPATH}/${MODEL}
#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph ScS -o ${OUTPATH}/${MODEL}
#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph pScS -o ${OUTPATH}/${MODEL}
#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph sScS -o ${OUTPATH}/${MODEL}


echo PSP, SKS, SKP and PKS, etc

#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph PKP -o ${OUTPATH}/${MODEL}
#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph SKS -o ${OUTPATH}/${MODEL}
#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph pSKP -o ${OUTPATH}/${MODEL}
#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph sSKP -o ${OUTPATH}/${MODEL}
#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph PKS -o ${OUTPATH}/${MODEL}
#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph pPKS -o ${OUTPATH}/${MODEL}
#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph sPKS -o ${OUTPATH}/${MODEL}


echo PP and PS and SP and SS

#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph PP -o ${OUTPATH}/${MODEL}
#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph pPP -o ${OUTPATH}/${MODEL}
#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph sPP -o ${OUTPATH}/${MODEL}
#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph PS -o ${OUTPATH}/${MODEL}
#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph pPS -o ${OUTPATH}/${MODEL}
#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph sPS -o ${OUTPATH}/${MODEL}
#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph SP -o ${OUTPATH}/${MODEL}
#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph pSP -o ${OUTPATH}/${MODEL}
#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph sSP -o ${OUTPATH}/${MODEL}
#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph SS -o ${OUTPATH}/${MODEL}
#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph pSS -o ${OUTPATH}/${MODEL}
#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph sSS -o ${OUTPATH}/${MODEL}
#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph SSS -o ${OUTPATH}/${MODEL}
#java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph SSSS -o ${OUTPATH}/${MODEL}

echo 410, 660, etc

### BAD? #java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph p410S -o ${OUTPATH}/${MODEL}
### BAD? #java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph P410S -o ${OUTPATH}/${MODEL}
### BAD? #java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph Pv410S -o ${OUTPATH}/${MODEL}
### BAD? #java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph s410P -o ${OUTPATH}/${MODEL}
### BAD? #java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph S410P -o ${OUTPATH}/${MODEL}
### BAD? #java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph Sv410P -o ${OUTPATH}/${MODEL}
### BAD? #java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph p660S -o ${OUTPATH}/${MODEL}
### BAD? #java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph P660S -o ${OUTPATH}/${MODEL}
### BAD? #java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph Pv660S -o ${OUTPATH}/${MODEL}
### BAD? #java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph s660P -o ${OUTPATH}/${MODEL}
### BAD? #java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph S660P -o ${OUTPATH}/${MODEL}
### BAD? #java edu.sc.seis.TauP.TauP_Table_NLL -nll ${NLL_GRID}  -mod ${MODEL} -ph Sv660P -o ${OUTPATH}/${MODEL}


