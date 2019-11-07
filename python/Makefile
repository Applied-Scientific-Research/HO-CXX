CXX = g++
CXXFLAGS = -g -fpic -O3 -mtune=native -std=c++11
SHLIB = -shared
LDFLAGS =

MKL_INCS = -DWITH_MKL -I${MKLROOT}/include
MKL_LIBS = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed \
             -lmkl_intel_lp64 \
             -lmkl_gnu_thread \
             -lmkl_core \
             -lmkl_def \
             -lmkl_avx2 \
             -lgomp -lpthread -lm -ldl -Wl,-rpath=${MKLROOT}/lib/intel64

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $<

all: libsp_solver_impl.so libcblock_sgs_impl.so libmkl_sparse_impl.so

libsp_solver_impl.so: sp_solver_impl.o
	$(CXX) $(CXXFLAGS) $(SHLIB) -o $@ $^ $(LDFLAGS)

libmkl_sparse_impl.so: mkl_sparse_impl.cpp
	$(CXX) $(CXXFLAGS) $(SHLIB) $(MKL_INCS) -o $@ $^ $(LDFLAGS) $(MKL_LIBS)

libcblock_sgs_impl.so: cblock_sgs.o
	$(CXX) $(CXXFLAGS) $(SHLIB) -o $@ $^ $(LDFLAGS)

clean:
	/bin/rm -fv *.o *.so
