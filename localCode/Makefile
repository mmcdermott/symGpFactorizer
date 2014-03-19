CXX = g++
WARNINGFLAGS = #-pedantic -Wall
LINKERFLAGS = #
CXXFLAGS = -std=c++11 -O2 $(WARNINGFLAGS) $(LINKERFLAGS)

TARGETS = cobFactorizer

all: $(TARGETS)

cobFactorizer: symGp.o tableau.o Vector.o Matrix.o main.o
	$(CXX) $^ $(CXXFLAGS) -o $@ 

profiler: symGp.o tableau.o Vector.o Matrix.o main.o
	$(CXX) $^ $(CXXFLAGS) -pg -o cobFactorizerProfilable

clean: 
	rm -f *.o cobFactorizer
