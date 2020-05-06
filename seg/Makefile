DEBUG=0
CXX=g++
OMP=-fopenmp -DOMP
CXXFLAGS=-g -O3 -Wall -DDEBUG=$(DEBUG) -std=c++11
LDFLAGS=-lm

CXXFILES = main.cpp
HFILES = mst.hpp par.hpp

all: seg-omp

seg-omp: $(CXXFILES) $(HFILES)
	$(CXX) $(CXXFLAGS) $(OMP) -o seg-omp $(CXXFILES) $(LDFLAGS)

clean:
	rm -f seg-omp