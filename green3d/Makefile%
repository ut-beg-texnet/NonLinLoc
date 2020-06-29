BINDIR = $(MYBIN)
F77FLAGS = -native -Bstatic -O4 -e

SUB  = addcel.o addrai.o calderive.o ajuste.o addseg.o modelrese.o \
       cellule.o critere.o derive.o raisat.o calcray.o modelray.o \
       deter.o expanh3.o green3.o limit.o derivebis.o inv3x3.o centra.o \
       indice.o infronti.o intcel.o inisource.o desvit.o cleanfront.o \
       intetra.o lenteur2.o modelvit.o celsat.o update.o cleanmap.o \
       newpoint.o propage.o runge4.o sortie.o stock.o segsat.o refpos.o \
       subcel.o testseg.o valtetra.o volume.o extirp.o main.o tabu.o \
       green3a.o

$(BINDIR)/green3 : $(SUB)
	f77 $(F77FLAGS) $(SUB) -o $(BINDIR)/green3

.SUFFIXES: .f .o
.f.o:
	f77 $(F77FLAGS) -c $*.f
