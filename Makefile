SYSTEM     ?= x86-64_sles10_4.1
LIBFORMAT  = static_pic
CPLEXVERSION ?= 2211

#------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where CPLEX and CONCERT are installed.
#
#------------------------------------------------------------

CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio$(CPLEXVERSION)/cplex
CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio$(CPLEXVERSION)/concert

# ---------------------------------------------------------------------
# Compiler selection
# ---------------------------------------------------------------------

CCC = g++

# ---------------------------------------------------------------------
# Compiler options
# ---------------------------------------------------------------------

CXXFLAGS = -std=c++17 -O2 -Wall -m64 -fPIC -fstrict-aliasing \
           -fexceptions -fopenmp -mavx2 -DIL_STD -pthread

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex \
            -L$(CONCERTLIBDIR) -lconcert -lm -fopenmp -pthread

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

INCLUDES = -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -Isrc -I.

# ---------------------------------------------------------------------
# Source files and objects
# ---------------------------------------------------------------------

SRCDIR = src
OBJDIR = obj

SRCS = $(SRCDIR)/main.cc \
       $(SRCDIR)/shortest_path.cc \
       $(SRCDIR)/reduction.cc \
       $(SRCDIR)/column_gen.cc

OBJS = $(patsubst $(SRCDIR)/%.cc,$(OBJDIR)/%.o,$(SRCS))

TARGET = mdvsp

# ---------------------------------------------------------------------
# Rules
# ---------------------------------------------------------------------

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CCC) $(CXXFLAGS) $(OBJS) -o $@ $(CCLNFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cc | $(OBJDIR)
	$(CCC) -c $(CXXFLAGS) $(INCLUDES) $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR)

clean:
	rm -rf $(TARGET) $(OBJDIR)
	rm -f *.mps *.ord *.sos *.lp *.sav *.net *.msg *.log *.clp

.PHONY: all clean
