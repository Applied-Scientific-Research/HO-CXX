APLLES_DIR ?= ../ApllesSolver
AMGCL_DIR ?= ../AMGCL

F90	= gfortran
F90FLAGS= -g -ggdb -gdwarf-2 -O2 -fopenmp -ffixed-line-length-none #-ffpe-summary='none'
CXX	= g++
CXXFLAGS= -g -ggdb -gdwarf-2 -O2 -std=c++11 -fopenmp

LDFLAGS =
LIBS    =

CXXSRC	= $(wildcard ../cxx/*.cpp)
OBJS	= $(CXXSRC:.cpp=.o)
DEPS	= $(CXXSRC:.cpp=.d)

#include $(APLLES_DIR)/user/Makefile.aplles
include $(APLLES_DIR)/lib/Makefile.aplles

F90FLAGS += $(APLLES_F90_MODULE)
CXXFLAGS += -I $(APLLES_DIR)/user
CXXFLAGS += -I ../cxx
CXXFLAGS += -I $(AMGCL_DIR)/include
CXXFLAGS += -I $(APLLES_DIR)/include
CXXFLAGS += -Wfatal-errors

all: seq omp

-include $(DEPS)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $< -c -o $@
	$(CXX) $(CXXFLAGS) -MM -MT $@ $< > $*.d

%.o : %.F90
	$(F90) $(F90FLAGS) $< -c -o $@

%.o : %.f90
	$(F90) $(F90FLAGS) $< -c -o $@

#VTM.o: $(APLLES_MOD)

seq: VTM.o $(OBJS) $(APLLES_OBJS)
	$(F90) $(F90FLAGS) -o VTM $< $(OBJS) $(APLLES_OBJS) $(LDFLAGS) $(APLLES_LIBS) $(APLLES_LINK_FLAGS)

omp: VTM.o $(OBJS) $(APLLES_OBJS)
	$(F90) $(F90FLAGS) -o VTM_mp $< $(OBJS) $(APLLES_OBJS) $(LDFLAGS) $(APLLES_MT_LIBS) $(APLLES_LINK_FLAGS)

clean:
	/bin/rm -fv *.o $(OBJS) *.mod $(DEPS)
