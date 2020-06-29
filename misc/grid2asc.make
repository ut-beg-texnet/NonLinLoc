# Makefile for NonLinLoc software package
#

BINDIR=${MYBIN}
#BINDIR=/home/anthony/bin
INCLUDE_DIR=./
#CCFLAGS=-O4
CCFLAGS=-g

GRID_LIB_OBJS=GridLib.o util.o geo.o octtree/octtree.o alomax_matrix/alomax_matrix.o alomax_matrix/alomax_matrix_svd.o matrix_statistics/matrix_statistics.o vector/vector.o ran1/ran1.o map_project.o

all : grid2asc_

# --------------------------------------------------------------------------
# grid2asc
#
OBJS1=${GRID_LIB_OBJS}
grid2asc_ : ${BINDIR}/grid2asc
${BINDIR}/grid2asc : grid2asc.o ${OBJS1}
	gcc grid2asc.o ${OBJS1} ${CCFLAGS} -o ${BINDIR}/grid2asc -lm
grid2asc.o : grid2asc.c GridLib.h
	gcc -c ${CCFLAGS} grid2asc.c
# --------------------------------------------------------------------------



# --------------------------------------------------------------------------
# Librarires
#
GridLib.o : GridLib.c GridLib.h
	gcc -c ${CCFLAGS} GridLib.c

util.o : ${INCLUDE_DIR}util.c ${INCLUDE_DIR}util.h
	gcc -c ${CCFLAGS} ${INCLUDE_DIR}util.c

octtree.o : ${INCLUDE_DIR}octtree.c ${INCLUDE_DIR}octtree.h
	gcc -c ${CCFLAGS} ${INCLUDE_DIR}octtree.c

nrutil.o : ${INCLUDE_DIR}nrutil.c ${INCLUDE_DIR}nrutil.h
	gcc -c ${CCFLAGS} ${INCLUDE_DIR}nrutil.c

nrmatrix.o : ${INCLUDE_DIR}nrmatrix.c
	gcc -c ${CCFLAGS} ${INCLUDE_DIR}nrmatrix.c

ran1.o : ${INCLUDE_DIR}ran1.c ${INCLUDE_DIR}ran1.h
	gcc -c ${CCFLAGS} ${INCLUDE_DIR}ran1.c

vector.o : ${INCLUDE_DIR}vector.c ${INCLUDE_DIR}vector.h
	gcc -c ${CCFLAGS} ${INCLUDE_DIR}vector.c

GMT_INCLUDE=./
map_project.o : ${GMT_INCLUDE}map_project.c
	gcc -c ${CCFLAGS} ${GMT_INCLUDE}map_project.c
#
# --------------------------------------------------------------------------
