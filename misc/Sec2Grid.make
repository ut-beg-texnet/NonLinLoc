BINDIR=${MYBIN}
BINDIR=/u/durance/home/lomax/bin
INCLUDE_DIR=./
CCFLAGS=-O3
#CCFLAGS=-g


all : Sec2Grid
distrib : Sec2Grid


# --------------------------------------------------------------------------
# Sec2Grid
#
OBJS2=GridLib.o util.o nrutil.o nrmatrix.o map_project.o ran1.o do_spline.o velmod.o 
#OBJS2=GridLib.o util.o nrutil.o nrmatrix.o map_project.o ran1.o spline_funct.o velmod.o libcommon.a
#OBJS2=GridLib.o util.o nrutil.o nrmatrix.o map_project.o ran1.o velmod.o libcommon.a
PVER=
Sec2Grid : ${BINDIR}/Sec2Grid
${BINDIR}/Sec2Grid : Sec2Grid${PVER}.o ${OBJS2}
	gcc Sec2Grid${PVER}.o  ${OBJS2} ${CCFLAGS} -o ${BINDIR}/Sec2Grid -lm
Sec2Grid${PVER}.o : Sec2Grid${PVER}.c GridLib.h
	gcc ${CCFLAGS} -c Sec2Grid${PVER}.c
# --------------------------------------------------------------------------


# --------------------------------------------------------------------------
# Librarires
#
GridLib.o : GridLib.c GridLib.h
	gcc -c ${CCFLAGS} GridLib.c

util.o : ${INCLUDE_DIR}util.c ${INCLUDE_DIR}util.h
	gcc -c ${CCFLAGS} ${INCLUDE_DIR}util.c

nrutil.o : ${INCLUDE_DIR}nrutil.c ${INCLUDE_DIR}nrutil.h
	gcc -c ${CCFLAGS} ${INCLUDE_DIR}nrutil.c

nrmatrix.o : ${INCLUDE_DIR}nrmatrix.c
	gcc -c ${CCFLAGS} ${INCLUDE_DIR}nrmatrix.c

ran1.o : ${INCLUDE_DIR}ran1.c ${INCLUDE_DIR}ran1.h
	gcc -c ${CCFLAGS} ${INCLUDE_DIR}ran1.c

vector.o : ${INCLUDE_DIR}vector.c ${INCLUDE_DIR}vector.h
	gcc -c ${CCFLAGS} ${INCLUDE_DIR}vector.c
velmod.o : ${INCLUDE_DIR}velmod.c
	gcc -c ${CCFLAGS} ${INCLUDE_DIR}velmod.c

SPLINE_DIR=./
do_spline.o : ${SPLINE_DIR}do_spline.c
	gcc -c ${CCFLAGS} ${SPLINE_DIR}do_spline.c

GMT_INCLUDE=./
map_project.o : ${GMT_INCLUDE}map_project.c
	gcc -c ${CCFLAGS} ${GMT_INCLUDE}map_project.c
#
# --------------------------------------------------------------------------
