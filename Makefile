# Makefile for ismail_code.cpp

CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++11
ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)

SRC = calibration.cpp
OBJ = $(SRC:.cpp=.o)
EXE = calibration

all: $(EXE)

$(EXE): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(ROOTLIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(ROOTFLAGS)

clean:
	rm -f $(OBJ) $(EXE)