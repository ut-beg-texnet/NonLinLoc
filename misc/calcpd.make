# Makefile for NonLinLoc software package
#
# Invocation:
#     Solaris: make all
#     Linux:   make -R all

BINDIR=${MYBIN}
BINDIR=/home/anthony/bin
INCLUDE_DIR=./
CCFLAGS=-O4
#CCFLAGS=-g


# --------------------------------------------------------------------------
# calcpd
#
OBJS1=GridLib.o util.o nrutil.o nrmatrix.o map_project.o
calcpd : ${BINDIR}/calcpd
${BINDIR}/calcpd : calcpd.o ${OBJS1}
	gcc calcpd.o ${OBJS1} ${CCFLAGS} -o ${BINDIR}/calcpd -lm
calcpd.o : calcpd.c GridLib.h ${INCLUDE_DIR}nrutil.h
	gcc -c ${CCFLAGS} calcpd.c
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
