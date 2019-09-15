APLLES_DIR ?= ../ApllesSolver

F90	= gfortran
F90FLAGS= -g -ggdb -gdwarf-2 -O2 -fopenmp -ffixed-line-length-none #-ffpe-summary='none'
CXX	= g++
CXXFLAGS= -g -ggdb -gdwarf-2 -O2 -std=c++11

LDFLAGS =
LIBS    =
OBJS	= ../cxx/laplacian.o

#include $(APLLES_DIR)/user/Makefile.aplles
include $(APLLES_DIR)/lib/Makefile.aplles

F90FLAGS += $(APLLES_F90_MODULE)
CXXFLAGS += -I $(APLLES_DIR)/user

all: seq

%.o : %.cxx
	$(CXX) $(CXXFLAGS) $< -c -o $@

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
	/bin/rm -fv *.o $(OBJS) *.mod
