BINDIR = $(MYBIN)
AJL_INCLUDE = ../include/
PVER = 1
CCFLAGS = -O4

GridProc: GridProc.o GridLib.o
	gcc GridProc.o GridLib.o -o $(BINDIR)/GridProc -lm

GridProc.o : GridProc.c GridLib.h GridGraphLib.h
	gcc $(CCFLAGS) -c GridProc.c

GridLib.o : GridLib.c GridLib.h
	gcc -c $(CCFLAGS) GridLib.c -o GridLib.o

