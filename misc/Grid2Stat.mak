BINDIR = $(MYBIN)
AJL_INCLUDE = ../include/
PVER = 1
CCFLAGS = -O4

Grid2Stat: Grid2Stat.o $(GRID_DIR)GridLib.o $(GRID_DIR)GridGraphLib.o
	gcc $(CCFLAGS) Grid2Stat.o $(GRID_DIR)GridLib.o $(GRID_DIR)GridGraphLib.o -o $(BINDIR)/Grid2Stat -lm

Grid2Stat.o : Grid2Stat.c GridLib.h GridGraphLib.h
	gcc $(CCFLAGS) -c Grid2Stat.c

GridLib.o : GridLib.c GridLib.h
	gcc -c $(CCFLAGS) GridLib.c -o GridLib.o

GridGraphLib.o : GridGraphLib.c GridGraphLib.h
	gcc -c $(CCFLAGS) GridGraphLib.c -o GridGraphLib.o

