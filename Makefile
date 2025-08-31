SYSTEM     = x86-64_sles10_4.1
LIBFORMAT  = static_pic
DEPS = CG.h gnrl.h lap.h Param.h system.h
#------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where CPLEX and CONCERT are installed.
#
#------------------------------------------------------------

CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio125/cplex
CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio125/concert
# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------

CCC = g++ -Ofast 

# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------

CCOPT = -m64 -O -fPIC -fstrict-aliasing -fexceptions -DIL_STD -pthread

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
LOCALDIR = .

CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -pthread


all:
	make all_cpp

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -I${LOCALDIR} 


CPP_EX = main

all_cpp: $(CPP_EX)

# ------------------------------------------------------------

clean :
	/bin/rm -rf $(CPP_EX)
	/bin/rm -rf *.mps *.ord *.sos *.lp *.sav *.net *.msg *.log *.clp *.o

# ------------------------------------------------------------
#
# The examples
#
main: main.o JV/lap.o
	$(CCC) $(CCFLAGS) main.o -o main $(CCLNFLAGS)


main.o: main.cc
	$(CCC) -c $(CCFLAGS) main.cc -o main.o


JV/lap.o: JV/lap.cpp JV/system.o 
	$(CCC) -c $(CCFLAGS) JV/lap.cpp -o JV/lap.o


JV/system.o: JV/system.cpp
	$(CCC) -c $(CCFLAGS) JV/system.cpp -o JV/system.o


# Local Variables:
# mode: makefile
# End:
