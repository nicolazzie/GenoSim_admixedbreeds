COMP=gfortran
OPTS=-O3 #  -g -fbounds-check -Wall
OPTMOD=-c
VPATH=progs

all: NUMKIN PROGPAR RNG SUBSGEN SUBSPEC SPECIES CLEAN

#DIRS: 
#	-mkdir -p ./GTYPES ./RESULTS

NUMKIN: Numeric_Kinds.f90
	$(COMP) $(OPTS) $(OPTMOD) $(VPATH)/Numeric_Kinds.f90

PROGPAR: Prog_Param.f90
	$(COMP) $(OPTS) $(OPTMOD) $(VPATH)/Prog_Param.f90

RNG: RNG.f90
	$(COMP) $(OPTS) $(OPTMOD) $(VPATH)/RNG.f90 

SUBSGEN: Subs_General.f90
	$(COMP) $(OPTS) $(OPTMOD) $(VPATH)/Subs_General.f90

SUBSPEC: Subs_Special.f90	
	$(COMP) $(OPTS) $(OPTMOD) $(VPATH)/Subs_Special.f90

SPECIES: 
	$(COMP) $(OPTS) $(VPATH)/Species_Data.f90 -o GenoSim.obj \
	Numeric_Kinds.o Prog_Param.o RNG.o Subs_General.o Subs_Special.o 
CLEAN:
	rm -f *o *.mod

