CXX = g++
CXXFLAGS = -g -fpic -O3 -mtune=native -std=c++11
SHLIB = -shared

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $<

all: libsp_dgetrs.so libcblock_sgs_impl.so

libsp_dgetrs.so: sp_dgetrs.o
	$(CXX) $(CXXFLAGS) $(SHLIB) -o $@ $^

libcblock_sgs_impl.so: cblock_sgs.o
	$(CXX) $(CXXFLAGS) $(SHLIB) -o $@ $^

clean:
	/bin/rm -fv *.o *.so
