CXX = g++
CXXFLAGS = -std=c++11 -Wall

all: matrix_mult schedule

matrix_mult: matrix_mult_1D.cpp
	$(CXX) $(CXXFLAGS) -o matrix_mult matrix_mult_1D.cpp

schedule: schedule.cpp
	$(CXX) $(CXXFLAGS) -o schedule schedule.cpp

clean:
	rm -f matrix_mult schedule