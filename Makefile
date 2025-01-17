# Makefile for calibration.cpp

CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++11
ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)

SRCS = calibration.cpp calibration2.cpp plot_overlaid_histograms.cpp feature_space.cpp landau_sim.cpp
OBJS = $(SRC:.cpp=.o)
EXES = calibration calibration2 plot_overlaid_histograms feature_space landau_sim

all: $(EXES)

calibration: calibration.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(ROOTLIBS)

calibration2: calibration2.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(ROOTLIBS)

plot_overlaid_histograms: plot_overlaid_histograms.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(ROOTLIBS)

feature_space: feature_space.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(ROOTLIBS)

landau_sim: landau_sim.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(ROOTLIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(ROOTFLAGS)

clean:
	rm -f $(OBJS) $(EXES)
