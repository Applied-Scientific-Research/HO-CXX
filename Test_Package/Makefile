APLLES_DIR ?= ../ApllesSolver
AMGCL_DIR ?= ../AMGCL

OPT	= -O3 -mtune=native
OMP	= -fopenmp
F90	= gfortran
F90FLAGS= -g -ggdb -gdwarf-2 $(OPT) -ffixed-line-length-none #-ffpe-summary='none'
CXX	= g++
CXXFLAGS= -g -ggdb -gdwarf-2 $(OPT) -std=c++17

CXXFLAGS += -I /usr/include/eigen3

LDFLAGS = $(OMP)
LIBS    =

CXXSRC	= ../cxx/laplacian.cpp
OBJS	= $(CXXSRC:.cpp=.o)
DEPS	= $(CXXSRC:.cpp=.d)

OBJS_MP	= $(CXXSRC:.cpp=_mp.o)
DEPS_MP	= $(CXXSRC:.cpp=_mp.d)

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

%_mp.o : %.cpp
	$(CXX) $(OMP) $(CXXFLAGS) $< -c -o $@
	$(CXX) $(OMP) $(CXXFLAGS) -MM -MT $@ $< > $*_mp.d

%_mp.o : %.F90
	$(F90) $(OMP) $(F90FLAGS) $< -c -o $@

%_mp.o : %.f90
	$(F90) $(OMP) $(F90FLAGS) $< -c -o $@

seq: VTM.o $(OBJS) $(APLLES_OBJS)
	$(F90) $(F90FLAGS) -o VTM $< $(OBJS) $(APLLES_OBJS) $(LDFLAGS) $(APLLES_LIBS) $(APLLES_LINK_FLAGS)

omp: VTM_mp.o $(OBJS_MP) $(APLLES_OBJS)
	$(F90) $(OMP) $(F90FLAGS) -o VTM_mp $< $(OBJS_MP) $(APLLES_OBJS) $(LDFLAGS) $(APLLES_MT_LIBS) $(APLLES_LINK_FLAGS)

clean:
	/bin/rm -fv *.o $(OBJS) $(OBJS_MP) *.mod $(DEPS)
