ICXX = icpc
CXXFLAGS = -std=c++17 -march=native -fopenmp
BOOSTPATH = ../lib/boost_1_74_0/ # customize to local boost installation
INCLUDE = -isystem $(BOOSTPATH) -I ./include/
LIB = -L $(BOOSTPATH)
SRC = ./src/
TEST =

CPPFILES = $(SRC)main.cpp \
					 $(SRC)CTimeDependentHjbSolver.cpp\
					 $(SRC)CTimeGrid.cpp\
					 $(SRC)CFoodDensity.cpp\
					 $(SRC)CTimeDependentTracer.cpp\
					 $(SRC)CEikonalSolver.cpp\
					 $(SRC)CTerrain.cpp

main: $(CPPFILES)
	$(CXX) $(CXXFLAGS) -O3 -o $@ $^ $(INCLUDE) $(LIB)

debug: $(CPPFILES)
	$(CXX) $(CXXFLAGS) -O1 -g -o $@ $^ $(INCLUDE) $(LIB)

clean:
	rm main

realclean: clean
	rm output/*
