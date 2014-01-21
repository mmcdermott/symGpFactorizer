CXX = g++
WARNINGFLAGS = -pedantic -Wall
CXXFLAGS = -std=c++11 -g  $(WARNINGFLAGS)

TARGETS = matrixTest

all: $(TARGETS)

matrixTest: myMatrix.o main.o
	$(CXX) $^ $(CXXFLAGS) -o $@ 

clean: 
	rm -f *.o matrixTest
