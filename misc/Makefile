# Makefile for NonLinLoc software package
#
# Invocation:
#     Solaris: make all
#     Linux:   make -R all


#  IMPORTANT: change here to set the directory for binary executable files
ifdef MYBIN
BINDIR=${MYBIN}
else
# with the following, binary executables will be placed in bin subdirectory of your home directory
#BINDIR=~/bin/
# if in doubt, use the following - binary executables will be placed in the current directory
BINDIR=.
endif

CC=gcc

#
GRID_FLOAT=
# double precision grid files (not yet supported fully for NonLinLoc)
#GRID_FLOAT=-D GRID_FLOAT_TYPE_DOUBLE

SRC_PATH=../src/
INC=-I../src/

# Linux, Mac OS X ?
CCFLAGS_BASIC =  -Wall  ${GRID_FLOAT}  $(INC)
# Mac OS X ?
#CCFLAGS_BASIC =  -Wall -D__APPLE__
#
# optimized
export CCFLAGS = -O3 $(CCFLAGS_BASIC)
#
# profile
#export CCFLAGS= -pg $(CCFLAGS_BASIC)
#
# debug - gdb, valgrind, ...
export CCFLAGS = $(CCFLAGS_BASIC) -g
# gdb:
# gdb --args exe_name
# gdb-apple --args exe_name
# lldb -- exe_name
# valgrind --leak-check=yes --dsymutil=yes exe_name <args>
# valgrind --leak-check=full --show-reachable=yes exe_name <args>



all : HypDiff_ alberto2hyp_ ingv_list2hyp_ hypoinv2hyp_ hypoe2hyp_ fpfit2hyp_ diffTime_ MergeFpfit_ MergeFmamp_ GridXY2XYZ_ isc_isf2hyp_ hyp2bull_ GridProc_

#GRID_LIB_OBJS=${SRC_PATH}/GridLib.o ${SRC_PATH}/util.o ${SRC_PATH}/geo.o ${SRC_PATH}/octtree/octtree.o ${SRC_PATH}/alomax_matrix/alomax_matrix.o ${SRC_PATH}/alomax_matrix/alomax_matrix_svd.o ${SRC_PATH}/matrix_statistics/matrix_statistics.o ${SRC_PATH}/vector/vector.o ${SRC_PATH}/ran1/ran1.o ${SRC_PATH}/map_project.o
GRID_LIB_OBJS=../src/CMakeFiles/GRID_LIB_OBJS.dir/GridLib.c.o ../src/CMakeFiles/GRID_LIB_OBJS.dir/util.c.o ../src/CMakeFiles/GRID_LIB_OBJS.dir/geo.c.o ../src/CMakeFiles/GRID_LIB_OBJS.dir/octtree/octtree.c.o ../src/CMakeFiles/GRID_LIB_OBJS.dir/alomax_matrix/alomax_matrix.c.o ../src/CMakeFiles/GRID_LIB_OBJS.dir/alomax_matrix/eigv.c.o ../src/CMakeFiles/GRID_LIB_OBJS.dir/alomax_matrix/alomax_matrix_svd.c.o ../src/CMakeFiles/GRID_LIB_OBJS.dir/matrix_statistics/matrix_statistics.c.o ../src/CMakeFiles/GRID_LIB_OBJS.dir/vector/vector.c.o ../src/CMakeFiles/GRID_LIB_OBJS.dir/ran1/ran1.c.o ../src/CMakeFiles/GRID_LIB_OBJS.dir/map_project.c.o

NLLOC_LIB_OBJS=../src/CMakeFiles/calc_crust_corr.dir/calc_crust_corr.c.o ../src/CMakeFiles/velmod.dir/velmod.c.o ../src/CMakeFiles/NLLOC_LIB_OBJS.dir/GridMemLib.c.o ../src/CMakeFiles/NLLOC_LIB_OBJS.dir/phaselist.c.o ../src/CMakeFiles/NLLOC_LIB_OBJS.dir/loclist.c.o ../src/CMakeFiles/NLLOC_LIB_OBJS.dir/otime_limit.c.o


PVER=1

# --------------------------------------------------------------------------
# HypDiff
#
OBJS4 = HypDiff.o $(GRID_LIB_OBJS)
HypDiff_ : $(BINDIR)/HypDiff
$(BINDIR)/HypDiff : $(OBJS4)
	${CC} $(OBJS4) $(CCFLAGS) -o $(BINDIR)/HypDiff -lm
HypDiff.o : HypDiff.c ${SRC_PATH}/include/GridLib.h
	${CC} $(CCFLAGS) -c HypDiff.c
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# alberto2hyp
#
OBJS5 = alberto2hyp.o $(GRID_LIB_OBJS)
alberto2hyp_ : $(BINDIR)/alberto2hyp
$(BINDIR)/alberto2hyp : $(OBJS5)
	${CC} $(OBJS5) $(CCFLAGS) -o $(BINDIR)/alberto2hyp -lm
alberto2hyp.o : alberto2hyp.c ${SRC_PATH}/include/GridLib.h
	${CC} $(CCFLAGS) -c alberto2hyp.c
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# ingv_list2hyp
#
OBJS5 = ingv_list2hyp.o $(GRID_LIB_OBJS)
ingv_list2hyp_ : $(BINDIR)/ingv_list2hyp
$(BINDIR)/ingv_list2hyp : $(OBJS5)
	${CC} $(OBJS5) $(CCFLAGS) -o $(BINDIR)/ingv_list2hyp -lm
ingv_list2hyp.o : ingv_list2hyp.c ${SRC_PATH}/include/GridLib.h
	${CC} $(CCFLAGS) -c ingv_list2hyp.c
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# hypoinv2hyp
#
OBJS6 = hypoinv2hyp.o $(GRID_LIB_OBJS)
hypoinv2hyp_ : $(BINDIR)/hypoinv2hyp
$(BINDIR)/hypoinv2hyp : $(OBJS6)
	${CC} $(OBJS6) $(CCFLAGS) -o $(BINDIR)/hypoinv2hyp -lm
hypoinv2hyp.o : hypoinv2hyp.c ${SRC_PATH}/include/GridLib.h
	${CC} $(CCFLAGS) -c hypoinv2hyp.c
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# hypoe2hyp
#
OBJS7=hypoe2hyp.o $(GRID_LIB_OBJS)
hypoe2hyp_ : ${BINDIR}/hypoe2hyp
${BINDIR}/hypoe2hyp : ${OBJS7}
	${CC} ${OBJS7} ${CCFLAGS} -o ${BINDIR}/hypoe2hyp -lm
hypoe2hyp.o : hypoe2hyp.c ${SRC_PATH}/include/GridLib.h
	${CC} ${CCFLAGS} -c hypoe2hyp.c
# --------------------------------------------------------------------------



# --------------------------------------------------------------------------
# fpfit2hyp
#
OBJS8=fpfit2hyp.o $(GRID_LIB_OBJS)
fpfit2hyp_ : ${BINDIR}/fpfit2hyp
${BINDIR}/fpfit2hyp : ${OBJS8}
	${CC} ${OBJS8} ${CCFLAGS} -o ${BINDIR}/fpfit2hyp -lm
fpfit2hyp.o : fpfit2hyp.c ${SRC_PATH}/include/GridLib.h
	${CC} ${CCFLAGS} -c fpfit2hyp.c
# --------------------------------------------------------------------------



# --------------------------------------------------------------------------
# MergeFpfit
#
OBJS8B=MergeFpfit.o $(GRID_LIB_OBJS)
MergeFpfit_ : ${BINDIR}/MergeFpfit
${BINDIR}/MergeFpfit : ${OBJS8B}
	${CC} ${OBJS8B} ${CCFLAGS} -o ${BINDIR}/MergeFpfit -lm
MergeFpfit.o : MergeFpfit.c ${SRC_PATH}/include/GridLib.h
	${CC} ${CCFLAGS} -c MergeFpfit.c
# --------------------------------------------------------------------------



# --------------------------------------------------------------------------
# MergeFmamp
#
OBJS8C=MergeFmamp.o $(GRID_LIB_OBJS)
MergeFmamp_ : ${BINDIR}/MergeFmamp
${BINDIR}/MergeFmamp : ${OBJS8C}
	${CC} ${OBJS8C} ${CCFLAGS} -o ${BINDIR}/MergeFmamp -lm
MergeFmamp.o : MergeFmamp.c ${SRC_PATH}/include/GridLib.h
	${CC} ${CCFLAGS} -c MergeFmamp.c
# --------------------------------------------------------------------------



# --------------------------------------------------------------------------
# diffTime
#
OBJS9 = diffTime.o $(GRID_LIB_OBJS)
diffTime_ : $(BINDIR)/diffTime
$(BINDIR)/diffTime : $(OBJS9)
	${CC} $(OBJS9) $(CCFLAGS) -o $(BINDIR)/diffTime -lm
diffTime.o : diffTime.c ${SRC_PATH}/include/GridLib.h
	${CC} $(CCFLAGS) -c diffTime.c
# --------------------------------------------------------------------------


# --------------------------------------------------------------------------
# GridXY2XYZ
#
OBJS11 = GridXY2XYZ.o $(GRID_LIB_OBJS)
GridXY2XYZ_ : $(BINDIR)/GridXY2XYZ
$(BINDIR)/GridXY2XYZ : $(OBJS11)
	${CC} $(OBJS11) $(CCFLAGS) -o $(BINDIR)/GridXY2XYZ -lm
GridXY2XYZ.o : GridXY2XYZ.c ${SRC_PATH}/include/GridLib.h
	${CC} $(CCFLAGS) -c GridXY2XYZ.c
# --------------------------------------------------------------------------


# --------------------------------------------------------------------------
# Loc2coul
#
OBJS12 = Loc2coul.o $(GRID_LIB_OBJS) ../okada/okada.o
Loc2coul_ : $(BINDIR)/Loc2coul
$(BINDIR)/Loc2coul : $(OBJS12)
	${CC} $(OBJS12) $(CCFLAGS) -o $(BINDIR)/Loc2coul -lm
Loc2coul.o : Loc2coul.c ${SRC_PATH}/include/GridLib.h
	${CC} $(CCFLAGS) -c Loc2coul.c
# --------------------------------------------------------------------------


# --------------------------------------------------------------------------
# FindEquivPhases
#
OBJS13 = FindEquivPhases.o $(GRID_LIB_OBJS)
FindEquivPhases_ : $(BINDIR)/FindEquivPhases
$(BINDIR)/FindEquivPhases : $(OBJS13)
	${CC} $(OBJS13) $(CCFLAGS) -o $(BINDIR)/FindEquivPhases -lm
FindEquivPhases.o : FindEquivPhases.c ${SRC_PATH}/include/GridLib.h
	${CC} $(CCFLAGS) -c FindEquivPhases.c
# --------------------------------------------------------------------------


# --------------------------------------------------------------------------
# Corr2MapGMT
#
OBJS14 = Corr2MapGMT.o $(GRID_LIB_OBJS)
Corr2MapGMT_ : $(BINDIR)/Corr2MapGMT
$(BINDIR)/Corr2MapGMT : $(OBJS14)
	${CC} $(OBJS14) $(CCFLAGS) -o $(BINDIR)/Corr2MapGMT -lm
Corr2MapGMT.o : Corr2MapGMT.c ${SRC_PATH}/include/GridLib.h
	${CC} $(CCFLAGS) -c Corr2MapGMT.c
# --------------------------------------------------------------------------


# --------------------------------------------------------------------------
# DDcc2MapGMT
#
OBJS14B = DDcc2MapGMT.o
DDcc2MapGMT_ : $(BINDIR)/DDcc2MapGMT
$(BINDIR)/DDcc2MapGMT : $(OBJS14B)
	${CC} $(OBJS14B) $(CCFLAGS) -o $(BINDIR)/DDcc2MapGMT -lm
DDcc2MapGMT.o : DDcc2MapGMT.c
	${CC} $(CCFLAGS) -c DDcc2MapGMT.c
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# isc_isf2hyp
#
OBJS15 = isc_isf2hyp.o $(GRID_LIB_OBJS)
isc_isf2hyp_ : $(BINDIR)/isc_isf2hyp
$(BINDIR)/isc_isf2hyp : $(OBJS15)
	${CC} $(OBJS15) $(CCFLAGS) -o $(BINDIR)/isc_isf2hyp -lm
isc_isf2hyp.o : isc_isf2hyp.c ${SRC_PATH}/include/GridLib.h
	${CC} $(CCFLAGS) -c isc_isf2hyp.c
# --------------------------------------------------------------------------


# --------------------------------------------------------------------------
# isc_isf2hyp
#
OBJS16 = nll2qml.o $(GRID_LIB_OBJS)
nll2qml_ : $(BINDIR)/nll2qml
$(BINDIR)/nll2qml : $(OBJS16)
	${CC} $(OBJS16) $(CCFLAGS) -o $(BINDIR)/nll2qml -lm `xml2-config --cflags` `xml2-config --libs`
nll2qml.o : nll2qml.c ${SRC_PATH}/include/GridLib.h
	${CC} $(CCFLAGS) -c nll2qml.c
# --------------------------------------------------------------------------



# --------------------------------------------------------------------------
# hyp2bull
#
OBJS17 = hyp2bull.o $(GRID_LIB_OBJS)
hyp2bull_ : $(BINDIR)/hyp2bull
$(BINDIR)/hyp2bull : $(OBJS17)
	${CC} $(OBJS17) $(CCFLAGS) -o $(BINDIR)/hyp2bull -lm
hyp2bull.o : hyp2bull.c ${SRC_PATH}/include/GridLib.h
	${CC} $(CCFLAGS) -c hyp2bull.c
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# GridProc
#
OBJS18 = GridProc.o $(GRID_LIB_OBJS)
GridProc_ : $(BINDIR)/GridProc
$(BINDIR)/GridProc : $(OBJS18)
	${CC} $(OBJS18) $(CCFLAGS) -o $(BINDIR)/GridProc -lm
GridProc.o : GridProc.c ${SRC_PATH}/include/GridLib.h
	${CC} $(CCFLAGS) -c GridProc.c
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# geo_test
#
OBJS19 = geo_test.o geo.o
geo_test_ : $(BINDIR)/geo_test
$(BINDIR)/geo_test : $(OBJS19)
	${CC} $(OBJS19) $(CCFLAGS) -o $(BINDIR)/geo_test -lm
geo_test.o : geo_test.c geo.h
	${CC} $(CCFLAGS) -c geo_test.c
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# compare_geog_trans
#
OBJS20 = compare_geog_trans.o $(GRID_LIB_OBJS)
compare_geog_trans_ : $(BINDIR)/compare_geog_trans
$(BINDIR)/compare_geog_trans : $(OBJS20)
	${CC} $(OBJS20) $(CCFLAGS) -o $(BINDIR)/compare_geog_trans -lm
compare_geog_trans.o : compare_geog_trans.c ${SRC_PATH}/include/GridLib.h
	${CC} $(CCFLAGS) -c compare_geog_trans.c
# --------------------------------------------------------------------------


# --------------------------------------------------------------------------
# read_grid2d
#
OBJS21 = read_grid2d.o $(GRID_LIB_OBJS)
read_grid2d_ : $(BINDIR)/read_grid2d
$(BINDIR)/read_grid2d : $(OBJS21)
	${CC} $(OBJS21) $(CCFLAGS) -o $(BINDIR)/read_grid2d -lm
read_grid2d.o : read_grid2d.c ${SRC_PATH}/include/GridLib.h
	${CC} $(CCFLAGS) -c read_grid2d.c
# --------------------------------------------------------------------------


# --------------------------------------------------------------------------
# bazSP2latlon
#
OBJS21 = bazSP2latlon.o $(GRID_LIB_OBJS)
bazSP2latlon_ : $(BINDIR)/bazSP2latlon
$(BINDIR)/bazSP2latlon : $(OBJS21)
	${CC} $(OBJS21) $(CCFLAGS) -o $(BINDIR)/bazSP2latlon -lm
bazSP2latlon.o : bazSP2latlon.c ${SRC_PATH}/include/GridLib.h
	${CC} $(CCFLAGS) -c bazSP2latlon.c
# --------------------------------------------------------------------------


# --------------------------------------------------------------------------
# hyp2pretty_phases_AlaskaHistorical
#
OBJS22 = hyp2pretty_phases_AlaskaHistorical.o $(GRID_LIB_OBJS)
hyp2pretty_phases_AlaskaHistorical_ : $(BINDIR)/hyp2pretty_phases_AlaskaHistorical
$(BINDIR)/hyp2pretty_phases_AlaskaHistorical : $(OBJS22)
	${CC} $(OBJS22) $(CCFLAGS) -o $(BINDIR)/hyp2pretty_phases_AlaskaHistorical -lm
hyp2pretty_phases_AlaskaHistorical.o : hyp2pretty_phases_AlaskaHistorical.c ${SRC_PATH}/include/GridLib.h
	${CC} $(CCFLAGS) -c hyp2pretty_phases_AlaskaHistorical.c
# --------------------------------------------------------------------------


# --------------------------------------------------------------------------
# hyp_gt_2_corr
#
OBJS23 = hyp_gt_2_corr.o ${SRC_PATH}/NLLoc${PVER}.o ${SRC_PATH}/NLLocLib.o $(NLLOC_LIB_OBJS) $(GRID_LIB_OBJS)
hyp_gt_2_corr_ : $(BINDIR)/hyp_gt_2_corr
$(BINDIR)/hyp_gt_2_corr : $(OBJS23)
	${CC} $(OBJS23) $(CCFLAGS) -o $(BINDIR)/hyp_gt_2_corr -lm
hyp_gt_2_corr.o : hyp_gt_2_corr.c ${SRC_PATH}/include/GridLib.h
	${CC} $(CCFLAGS) -c hyp_gt_2_corr.c
# --------------------------------------------------------------------------


# --------------------------------------------------------------------------
# General compile rule
#

%.o : %.c
	$(CC) -c $(CCFLAGS)   $(OPTIONS) $< -o $@

#
# --------------------------------------------------------------------------



# --------------------------------------------------------------------------
# Misc
#

clean :
	rm -f *.o

#
# --------------------------------------------------------------------------



