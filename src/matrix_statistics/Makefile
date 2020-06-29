
# Build environment can be configured the following
# environment variables:
#   CC : Specify the C compiler to use
#   CFLAGS : Specify compiler options to use

#  IMPORTANT: change here to set the directory for binary executable files
ifdef MYBIN
BINDIR=${MYBIN}
else
# with the following, binary executables will be placed in bin subdirectory of your home directory
#BINDIR=~/bin/
# if in doubt, use the following - binary executables will be placed in the current directory
BINDIR=.
endif


# Options specific for GCC
ifndef CC
CC = gcc
endif
#
ifndef CCFLAGS
#
CCFLAGS_BASIC =  -Wall -fPIC -I..
#
# optimized
CCFLAGS = -O3 $(CCFLAGS_BASIC)
#
# profile
#CCFLAGS=-O3 -pg $(CCFLAGS_BASIC)
#
# debug - gdb, valgrind, ...
#CCFLAGS = $(CCFLAGS_BASIC) -g
# valgrind --leak-check=yes  exe_name <args>
# valgrind --leak-check=full --show-reachable=yes exe_name exe_name <args>
endif


OBJS =  matrix_statistics.o

all: $(OBJS)

clean:
	rm -f *.o
	rm -f *.a

test: ${BINDIR}/test_ellipsoid

# Implicit rule for building object files
%.o: %.c
	$(CC) -c $< $(CCFLAGS)

# --------------------------------------------------------------------------
# test_ellipsoid
#
OBJ_TEST1= test_ellipsoid.o $(OBJS) ../vector/vector.o ../alomax_matrix/alomax_matrix.o ../alomax_matrix/alomax_matrix_svd.o
${BINDIR}/test_ellipsoid : ${OBJ_TEST1}
	$(CC) ${OBJ_TEST1} ${CCFLAGS} -o ${BINDIR}/test_ellipsoid -lm
# --------------------------------------------------------------------------
