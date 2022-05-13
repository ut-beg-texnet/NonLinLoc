#!/bin/bash

TEST_IT_LOCAL=YES
TEST_IT_GLOBAL=NO
TEST_IT_FUNC=NO

# use GMT4
export PATH=/opt/gmt/bin:$PATH

#NLLOC_VER=7.00.03
NLLOC_VER=beta
#NLLOC_VER=MacOSX_TEST
#NLLOC_VER=SC3_TEST
#NLLOC_VER=NLDiffLoc_TEST

OS_IDENT=x86_64-apple-darwin14

SEISMICIY_JAR=/Users/anthony/www/seismicity/beta/SeismicityViewer50.jar



PWD=$(pwd)

NLL_DIR=${PWD}
SOURCE_DIR=${NLL_DIR}/src/
WEB_DIR=distributions/www_${NLLOC_VER}

echo "PATH=${PATH}"
MYBIN=../bin
PATH=${PWD}/${WEB_DIR}/bin:${PATH}
echo "PATH=${PATH}"

#PS_VIEWER=ghostview

#  create/clean directories./b

rm -rf ${WEB_DIR}
mkdir ${WEB_DIR}
mkdir ${WEB_DIR}/src


# copy files

echo "copying *.c source ..."
cp -p ${SOURCE_DIR}/NLLoc_main.c ${SOURCE_DIR}/NLLoc_func_test.c ${SOURCE_DIR}/ttime_func_test.c \
	${SOURCE_DIR}/mag_func_test.c ${SOURCE_DIR}/NLLoc1.c ${SOURCE_DIR}/NLLocLib.c \
	${SOURCE_DIR}/PhsAssoc.c ${SOURCE_DIR}/Vel2Grid1.c ${SOURCE_DIR}/Grid2Time1.c ${SOURCE_DIR}/Time2Angles1.c \
	${SOURCE_DIR}/Grid2GMT.c ${SOURCE_DIR}/LocSum.c ${SOURCE_DIR}/scat2latlon.c ${SOURCE_DIR}/Time2EQ1.c ${SOURCE_DIR}/GridLib.c \
	${SOURCE_DIR}/GridMemLib.c ${SOURCE_DIR}/GridGraphLib.c ${SOURCE_DIR}/util.c \
	${SOURCE_DIR}/velmod.c ${SOURCE_DIR}/map_project.c \
	${SOURCE_DIR}/NLDiffLoc.c ${SOURCE_DIR}/Loc2ddct.c \
	${SOURCE_DIR}/geo.c ${SOURCE_DIR}/Time_3d_NLL.c ${SOURCE_DIR}/hypoe2hyp.c \
	${SOURCE_DIR}/fpfit2hyp.c ${SOURCE_DIR}/calc_crust_corr.c ${SOURCE_DIR}/oct2grid.c \
	${SOURCE_DIR}/phaselist.c ${SOURCE_DIR}/loclist.c ${SOURCE_DIR}/otime_limit.c \
	${SOURCE_DIR}/Vel2Grid3D.c ${SOURCE_DIR}/interface2fmm.c ${SOURCE_DIR}/fmm2grid.c ${SOURCE_DIR}/GridCascadingDecimate.c \
	${WEB_DIR}/src

echo "copying *.h source ..."
cp -p ${SOURCE_DIR}/NLLocLib.h ${SOURCE_DIR}/GridLib.h ${SOURCE_DIR}/GridMemLib.h \
	${SOURCE_DIR}/GridGraphLib.h \
	${SOURCE_DIR}/map_project.h ${SOURCE_DIR}/geo.h \
	${SOURCE_DIR}/util.h ${SOURCE_DIR}/velmod.h \
	${SOURCE_DIR}/calc_crust_corr.h ${SOURCE_DIR}/crust_corr_model.h \
	${SOURCE_DIR}/crust_type.h ${SOURCE_DIR}/crust_type_key.h ${SOURCE_DIR}/phaseloclist.h \
	${SOURCE_DIR}/otime_limit.h  ${WEB_DIR}/src

# alomax libraries .c and .h
for ALOMAX_LIB in octtree ran1 vector alomax_matrix matrix_statistics ; do
mkdir ${WEB_DIR}/src/${ALOMAX_LIB}
cp -p ${SOURCE_DIR}/${ALOMAX_LIB}/*.h  ${WEB_DIR}/src/${ALOMAX_LIB}
cp -p ${SOURCE_DIR}/${ALOMAX_LIB}/*.c  ${WEB_DIR}/src/${ALOMAX_LIB}
done
# alomax libraries .h
for ALOMAX_LIB in geometry ; do
mkdir ${WEB_DIR}/src/${ALOMAX_LIB}
cp -p ${SOURCE_DIR}/${ALOMAX_LIB}/*.h  ${WEB_DIR}/src/${ALOMAX_LIB}
done

echo "copying Makefile ..."
cp -p ${SOURCE_DIR}/Makefile ${WEB_DIR}/src

echo "copying README ..."
cp -p ${NLL_DIR}/README.md ${WEB_DIR}/src/README.txt

echo "copying CHANGE_NOTES ..."
cp -p ${NLL_DIR}/CHANGE_NOTES.txt ${WEB_DIR}/src


# sample location

echo "copying ${NLL_DIR}/nlloc_sample ..."
cp -pr ${NLL_DIR}/nlloc_sample ${WEB_DIR}


# sample global location

echo "copying ${NLL_DIR}/nlloc_global_sample ..."
rm -rf ${NLL_DIR}/nlloc_global_sample/loc/*
rm -rf ${NLL_DIR}/nlloc_global_sample/taup/ak135/*
cp -pr ${NLL_DIR}/nlloc_global_sample ${WEB_DIR}
mkdir ${WEB_DIR}/nlloc_global_sample
mkdir ${WEB_DIR}/nlloc_global_sample/java
mkdir ${WEB_DIR}/nlloc_global_sample/java/edu
mkdir ${WEB_DIR}/nlloc_global_sample/java/edu/sc
mkdir ${WEB_DIR}/nlloc_global_sample/java/edu/sc/seis
mkdir ${WEB_DIR}/nlloc_global_sample/java/edu/sc/seis/TauP
cp -pr /Users/anthony/java/edu/sc/seis/TauP/*.java ${WEB_DIR}/nlloc_global_sample/java/edu/sc/seis/TauP
#${JAVA_PATH_DISTRIBUTION}/javac ${WEB_DIR}/nlloc_global_sample/java/edu/sc/seis/TauP/TauP_Table_NLL.java
cp -pr /Users/anthony/java/edu/sc/seis/TauP/StdModels ${WEB_DIR}/nlloc_global_sample/java/edu/sc/seis/TauP
${JAVA_PATH_DISTRIBUTION}/javac ${WEB_DIR}/nlloc_global_sample/java/edu/sc/seis/TauP/*.java
cp -p ${WEB_DIR}/nlloc_global_sample/java/edu/sc/seis/TauP/TauP_Table_NLL.* /Users/anthony/www/nlloc/java
${JAVA_PATH_DISTRIBUTION}/jar -cf ${WEB_DIR}/nlloc_global_sample/TauP_NLL.jar -C ${WEB_DIR}/nlloc_global_sample/java/ .
cp -p ${WEB_DIR}/nlloc_global_sample/TauP_NLL.jar /Users/anthony/www/nlloc/java
cp -p ${SEISMICIY_JAR} ${WEB_DIR}/nlloc_global_sample/SeismicityViewer.jar


# sample nll_func

echo "copying ${NLL_DIR}/nll_func ..."
cp -pr ${NLL_DIR}/nll_func ${WEB_DIR}
echo "copying ${NLL_DIR}/ttime_func ..."
cp -pr ${NLL_DIR}/ttime_func ${WEB_DIR}
echo "copying ${NLL_DIR}/mag_func ..."
cp -pr ${NLL_DIR}/mag_func ${WEB_DIR}



cd ${WEB_DIR}


# compile

pwd
mkdir bin
cd src
pwd

echo "compiling ..."

export MYBIN=../bin
make clean
make -R distrib
if [ ${TEST_IT_FUNC} = YES ]; then
make -R NLLoc_func_test_
make -R ttime_func_test_
make -R mag_func_test_
fi
rm *.o
rm */*.o
cd ..
pwd


################################

if [ ${TEST_IT_LOCAL} = YES ]; then

	# apply and test sample location

	cd nlloc_sample
	pwd
	./run_nlloc_sample.sh

	# clean
	rm gmt/*

	# save original output files
	mkdir original_output
	rm -rf original_output/*
	mkdir original_output/loc
	mv loc/* original_output/loc
	mkdir original_output/model
	mv model/* original_output/model
	mkdir original_output/time
	mv time/* original_output/time
	mkdir original_output/obs
	mv obs/* original_output/obs
	mkdir original_output/obs_synth
	mv obs_synth/* original_output/obs_synth
	cd ..
	pwd

fi





################################

if [ ${TEST_IT_GLOBAL} = YES ]; then

	# test sample global location

	cd nlloc_global_sample
	pwd
	./run_global_sample.sh

	# clean
#rm -rf loc/*
	rm -rf taup/ak135/*

	cd ..
	pwd

fi





################################

if [ ${TEST_IT_FUNC} = YES ]; then

	export PATH=../bin:${PATH}

	# nll_func
	cd nll_func
	pwd
	./run_func.sh
	# clean
	rm -rf out/*
	cd ..
	pwd

# ttime_func
	cd ttime_func
	pwd
	./run_func.sh
	cd ..
	pwd

# mag_func
	cd mag_func
	pwd
	./run_func.sh
	cd ..
	pwd

fi



################################

# make compressed tar files

echo "making compressed tar files ..."
pwd

export TAR_DIR=./tar
mkdir ${TAR_DIR}
rm -rf ${TAR_DIR}/*

tar czf ${TAR_DIR}/NLL${NLLOC_VER}_src.tgz src
# zip version
#tar cf ${TAR_DIR}/NLL${NLLOC_VER}_src.tar src
#zip ${TAR_DIR}/NLL${NLLOC_VER}_src.tar.zip ${TAR_DIR}/NLL${NLLOC_VER}_src.tar


tar cvzf ${TAR_DIR}/NLL${NLLOC_VER}_bin_${OS_IDENT}.tgz bin

tar czf ${TAR_DIR}/NLL${NLLOC_VER}_samples.tgz nlloc*sample/*

# custom source files
ETH_DIR=custom_eth
cd ../..
tar cvzf ${WEB_DIR}/${TAR_DIR}/NLL${NLLOC_VER}_${ETH_DIR}_src.tgz ${ETH_DIR}
cd ${WEB_DIR}

# nll_func, ttime_func, mag_func
tar czf ${TAR_DIR}/NLL${NLLOC_VER}_func.tgz nll_func ttime_func mag_func


mkdir /Users/anthony/www/nlloc/soft${NLLOC_VER}
mkdir /Users/anthony/www/nlloc/soft${NLLOC_VER}/tar
cp ${TAR_DIR}/* /Users/anthony/www/nlloc/soft${NLLOC_VER}/tar


echo "copying README ..."
cp -p ${NLL_DIR}/README.txt /Users/anthony/www/nlloc/soft${NLLOC_VER}/tar
echo "copying CHANGE_NOTES ..."
cp -p ${NLL_DIR}/CHANGE_NOTES.txt /Users/anthony/www/nlloc/soft${NLLOC_VER}/tar
echo "copying NLDiffLoc users_guide ..."
cp -p /Users/anthony/soft/NLLoc/doc/NLDiffLoc*_users_guide.pdf /Users/anthony/www/nlloc/soft${NLLOC_VER}/tar

cd ..
pwd

echo "-"
echo "Done"

# Twitter message
echo "Twitter message"
echo "NonLinLoc: New RELEASE version available: NLL${NLLOC_VER}: http://alomax.net//nlloc"
echo "NonLinLoc: New BETA version available: NLL${NLLOC_VER}: http://alomax.net//nlloc/softbeta/tar (minor changes)"

echo "IMPORTANT: For distribution, remember to run ~/java_distribution/jdist_TauP_Table_NLL_.bash"

