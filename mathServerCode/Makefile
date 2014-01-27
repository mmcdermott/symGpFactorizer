CXX = /shared/local/sw/gcc/gcc/bin/g++
WARNINGFLAGS = -pedantic -Wall
LINKERFLAGS = -Wl,-rpath,/shared/local/sw/gcc/gcc/lib64
CXXFLAGS = -gdwarf-2 -gstrict-dwarf -std=c++11 -g $(WARNINGFLAGS) $(LINKERFLAGS)

TARGETS = cobFactorizer

all: $(TARGETS)

cobFactorizer: symGp.o tableau.o Vector.o Matrix.o main.o
	$(CXX) $^ $(CXXFLAGS) -o $@ 

profiler: symGp.o tableau.o Vector.o Matrix.o main.o
	$(CXX) $^ $(CXXFLAGS) -pg -o cobFactorizerProfilable

clean: 
	rm -f *.o cobFactorizer
