CXX = g++
WARNINGFLAGS = #-pedantic -Wall
CXXFLAGS = -std=c++11 -g -I ~/boost_1_54_0 -l boost_system -l boost_filesystem $(WARNINGFLAGS)

TARGETS = cobFactorizer

all: $(TARGETS)

cobFactorizer: symGp.o tableau.o Vector.o Matrix.o main.o
	$(CXX) $^ $(CXXFLAGS) -o $@ 

profiler: symGp.o tableau.o Vector.o Matrix.o main.o
	$(CXX) $^ $(CXXFLAGS) -pg -o cobFactorizerProfilable

clean: 
	rm -f *.o cobFactorizer
