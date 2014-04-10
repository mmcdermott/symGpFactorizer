CXX = g++
WARNINGFLAGS = #-pedantic -Wall
LINKERFLAGS = #
OPTIMIZATION = #-O2
CXXFLAGS = -std=c++11 $(OPTIMIZATION) $(WARNINGFLAGS) $(LINKERFLAGS)

TARGETS = cobFactorizer

all: $(TARGETS)

cobFactorizer: symGp.o tableau.o Vector.o Matrix.o main.o
	$(CXX) $^ $(CXXFLAGS) -o $@ 

profiler: symGp.o tableau.o Vector.o Matrix.o main.o
	$(CXX) $^ $(CXXFLAGS) -pg -o cobFactorizerProfilable

rationals: Rational.o rationalTest.o
	$(CXX) $^ $(CXXFLAGS) -o rationalTest

sqrteds: Rational.o Sqrt.o sqrtTest.o
	$(CXX) $^ $(CXXFLAGS) -o sqrtTest

clean: 
	rm -f *.o cobFactorizer
