#BINDIR=${MYBIN}
BINDIR=/u/durance/home/lomax/bin
INCLUDE_DIR=./
CCFLAGS=-O4
#CCFLAGS=-g

asc2grid: asc2grid.o GridLib.o
	gcc asc2grid.o GridLib.o  util.o  nrutil.o nrmatrix.o ran1.o map_project.o -o $(BINDIR)/asc2grid -lm

asc2grid.o : asc2grid.c GridLib.h
	gcc $(CCFLAGS) -c asc2grid.c

GridLib.o : GridLib.c GridLib.h
	gcc -c $(CCFLAGS) GridLib.c

