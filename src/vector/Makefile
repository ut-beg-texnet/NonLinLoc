# Makefile for picker test
#
# Invocation:
#     Solaris: make -f Make_picker_func_test
#     Linux:   gmake -Rf Make_picker_func_test

#  IMPORTANT: change here to set the directory for binary executable files
ifdef MYBIN
BINDIR=${MYBIN}
else
# with the following, binary executables will be placed in bin subdirectory of your home directory
#BINDIR=~/bin/
# if in doubt, use the following - binary executables will be placed in the current directory
BINDIR=.
endif
INCLUDE_DIR=./

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



OBJS=vector.o

all : ${OBJS}

%.o : %.c
	$(CC) -c $(CCFLAGS) $< -o $@


clean :
	rm -f *.o
	rm -f *.a

