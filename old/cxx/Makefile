CXX	= g++
CXXFLAGS= -g -ggdb -gdwarf-2 -O2 -std=c++17 -fopenmp
CXXFLAGS += -Wfatal-errors

CXXSRC	= $(wildcard test_*.cpp)
EXES	= $(CXXSRC:.cpp=.exe)
OBJS	= $(CXXSRC:.cpp=.o)
DEPS	= $(CXXSRC:.cpp=.d)

all: $(EXES)

-include $(DEPS)

%.o : %.cpp

%.exe : %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@
	$(CXX) $(CXXFLAGS) -MM -MT $@ $< > $*.d

clean:
	/bin/rm -fv $(DEPS) $(EXES)
