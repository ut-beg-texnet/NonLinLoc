BINDIR = $(MYBIN)
INCLUDE_DIR = ./
CCFLAGS = -O4
#CCFLAGS = -ggdb


all : Loc2Mod
distrib : Loc2Mod


# --------------------------------------------------------------------------
# Loc2Mod
#
OBJS2 = GridLib.o util.o nrutil.o map_project.o nrmatrix.o ran1.o
PVER = 
Loc2Mod : $(BINDIR)/Loc2Mod
$(BINDIR)/Loc2Mod : Loc2Mod$(PVER).o $(OBJS2)
	gcc Loc2Mod$(PVER).o  $(OBJS2) $(CCFLAGS) -o $(BINDIR)/Loc2Mod -lm
Loc2Mod$(PVER).o : Loc2Mod$(PVER).c GridLib.h
	gcc $(CCFLAGS) -c Loc2Mod$(PVER).c
# --------------------------------------------------------------------------


# --------------------------------------------------------------------------
# Librarires
#
GridLib.o : GridLib.c GridLib.h
	gcc -c $(CCFLAGS) GridLib.c

util.o : $(INCLUDE_DIR)util.c $(INCLUDE_DIR)util.h
	gcc -c $(CCFLAGS) $(INCLUDE_DIR)util.c

nrutil.o : $(INCLUDE_DIR)nrutil.c $(INCLUDE_DIR)nrutil.h
	gcc -c $(CCFLAGS) $(INCLUDE_DIR)nrutil.c

nrmatrix.o : $(INCLUDE_DIR)nrmatrix.c
	gcc -c $(CCFLAGS) $(INCLUDE_DIR)nrmatrix.c

ran1.o : $(INCLUDE_DIR)ran1.c $(INCLUDE_DIR)ran1.h
	gcc -c $(CCFLAGS) $(INCLUDE_DIR)ran1.c

vector.o : $(INCLUDE_DIR)vector.c $(INCLUDE_DIR)vector.h
	gcc -c $(CCFLAGS) $(INCLUDE_DIR)vector.c
velmod.o : $(INCLUDE_DIR)velmod.c
	gcc -c $(CCFLAGS) $(INCLUDE_DIR)velmod.c

SPLINE_DIR = ./
do_spline.o : $(SPLINE_DIR)do_spline.c
	gcc -c $(CCFLAGS) $(SPLINE_DIR)do_spline.c

GMT_INCLUDE = ./
map_project.o : $(GMT_INCLUDE)map_project.c
	gcc -c $(CCFLAGS) $(GMT_INCLUDE)map_project.c
#
# --------------------------------------------------------------------------
