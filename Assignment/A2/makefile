# Variables
CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++17
SRC_BFS = BFSTopo.cpp alg_graphs.cpp
SRC_EDGE = EdgeDetect.cpp alg_graphs.cpp
OBJ_BFS = $(SRC_BFS:.cpp=.o)
OBJ_EDGE = $(SRC_EDGE:.cpp=.o)
TARGET_BFS = BFSTopo
TARGET_EDGE = EdgeDetect

# Default target
all: $(TARGET_BFS) $(TARGET_EDGE)

# Link object files to create the BFSTopo executable
$(TARGET_BFS): $(OBJ_BFS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Link object files to create the EdgeDetect executable
$(TARGET_EDGE): $(OBJ_EDGE)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile source files into object files
%.o: %.cpp alg_graphs.h
	$(CXX) $(CXXFLAGS) -c $<

# Clean up compiled files
clean:
	rm -f $(OBJ_BFS) $(OBJ_EDGE) $(TARGET_BFS) $(TARGET_EDGE)

# Phony targets
.PHONY: all clean