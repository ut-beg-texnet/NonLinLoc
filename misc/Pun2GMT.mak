BINDIR = $(MYBIN)
AJL_INCLUDE = ../include/
PVER = 1
CCFLAGS = -O4

Pun2GMT: Pun2GMT.o GridLib.o
	gcc Pun2GMT.o GridLib.o -o $(BINDIR)/Pun2GMT -lm

Pun2GMT.o : Pun2GMT.c GridLib.h GridGraphLib.h
	gcc $(CCFLAGS) -c Pun2GMT.c

GridLib.o : GridLib.c GridLib.h
	gcc -c $(CCFLAGS) GridLib.c -o GridLib.o

