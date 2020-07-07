CXX=g++
CXXOPTS=-O2 -fopenmp
INC=-I./include
LIBS=-L./libs -laplles_mp -lxml2 -llapack

all : 2DVortTrans.bin

2DVortTrans.bin : aplles_interface_fortran.F90 aplles_interface.cxx VTM.f90
	$(CXX) $(CXXOPTS) $(INC) -lgfortran $^ $(LIBS) -o $@

clean :
	rm -f 2DVortTrans.bin
