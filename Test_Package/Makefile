APLLES_DIR ?= ../../ApllesSolver
AMGCL_INC ?= ../../amgcl/
EIGEN3_INC ?= /usr/include/eigen3/
HYPRE_DIR ?= ../../hypre/

OPT	= -O3
OMP	= -fopenmp
DBG	= -g #-ggdb -gdwarf-2
F90	?= gfortran
F90FLAGS= $(DBG) $(OPT)
MODFLAGS=
CXX	?= g++
CXXFLAGS= $(DBG) $(OPT) -std=c++17
LDFLAGS	= $(OMP)
LIBS	=

ifeq ($(F90),gfortran)
  F90FLAGS += -mtune=native -ffixed-line-length-none #-ffpe-summary='none'
endif
ifeq ($(CXX),g++)
  CXXFLAGS += -mtune=native
endif

ifeq ($(F90),ifort)
  F90FLAGS += -xHost
endif
ifeq ($(CXX),icpc)
  CXXFLAGS += -xHost -Wno-unknown-pragmas
endif

CXXSRC	= ../cxx/laplacian.cpp
OBJS	= $(CXXSRC:.cpp=.o)
DEPS	= $(CXXSRC:.cpp=.d)

OBJS_MP	= $(CXXSRC:.cpp=_mp.o)
DEPS_MP	= $(CXXSRC:.cpp=_mp.d)

include $(APLLES_DIR)/lib/Makefile.aplles

#APLLES_LIBS = -laplles_debug

F90FLAGS += $(APLLES_F90_MODULE)
CXXFLAGS += -I $(APLLES_DIR)/user
CXXFLAGS += -I ../cxx
#CXXFLAGS += -I $(AMGCL_INC) -DENABLE_AMGCL
#CXXFLAGS += -I $(APLLES_DIR)/include
#CXXFLAGS += -I $(EIGEN3_INC) -DENABLE_EIGEN
CXXFLAGS += -I $(APLLES_DIR)/include -DENABLE_APLLES
CXXFLAGS += -Wfatal-errors

CXXFLAGS += -I $(HYPRE_DIR)/include -DENABLE_HYPRE
LDFLAGS  += -lHYPRE -L $(HYPRE_DIR)/lib

ifeq ($(F90),gfortran)
  F90FLAGS += -I. -J.
endif
ifeq ($(F90),ifort)
  F90FLAGS += -I. -module .
endif

all: seq omp

MODS = aplles_interface_fortran.o
MODS_MP += aplles_interface_fortran_mp.o

aplles_interface_fortran.o: $(APLLES_DIR)/user/aplles_interface_fortran.F90
	$(F90) $(F90FLAGS) $(MODFLAGS) $< -c -o $@

aplles_interface_fortran_mp.o: $(APLLES_DIR)/user/aplles_interface_fortran.F90
	$(F90) $(F90FLAGS) $(MODFLAGS) $< -c -o $@

-include $(DEPS) $(DEPS_MP)

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

seq: $(MODS) VTM.o $(OBJS) $(APLLES_OBJS)
	$(F90) $(F90FLAGS) -o VTM $^ $(LDFLAGS) $(LIBS) $(APLLES_LIBS) $(APLLES_LINK_FLAGS)

omp: $(MODS_MP) VTM_mp.o $(OBJS_MP) $(APLLES_OBJS)
	$(F90) $(OMP) $(F90FLAGS) -o VTM_mp $^ $(LDFLAGS) $(LIBS) $(APLLES_MT_LIBS) $(APLLES_LINK_FLAGS)

clean:
	/bin/rm -fv *.o $(OBJS) $(OBJS_MP) *.mod $(DEPS)
