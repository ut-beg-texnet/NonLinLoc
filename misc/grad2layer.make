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
# grad2layer
#

grad2layer : ${BINDIR}/grad2layer
${BINDIR}/grad2layer : grad2layer.o
	gcc grad2layer.o ${CCFLAGS} -o ${BINDIR}/grad2layer -lm
grad2layer.o : grad2layer.c
	gcc -c ${CCFLAGS} grad2layer.c
# --------------------------------------------------------------------------

