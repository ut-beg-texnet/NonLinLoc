BINDIR=${MYBIN}
INCLUDE_DIR=./
CCFLAGS=-O3
#CCFLAGS=-O3 -pg
#CCFLAGS=-g


# --------------------------------------------------------------------------
# HypAssoc
#
OBJS9=HypAssoc.o GridLib.o geo.o util.o nrutil.o nrmatrix.o map_project.o
HypAssoc : ${BINDIR}/HypAssoc
${BINDIR}/HypAssoc : ${OBJS9}
	gcc ${OBJS9} ${CCFLAGS} -o ${BINDIR}/HypAssoc -lm
HypAssoc.o : HypAssoc.c GridLib.h
	gcc ${CCFLAGS} -c HypAssoc.c
# --------------------------------------------------------------------------
