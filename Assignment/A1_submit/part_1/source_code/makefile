# Makefile for compiling foursum with alg_stopwatch

# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -std=c++17 -Wall -Wextra -O2

# Target executable name
TARGET = foursum

# Source files
SRCS = four_sum.cpp alg_stopwatch.cpp

# Header files
HEADERS = alg_stopwatch.h

# Object files (replace .cpp with .o in SRCS)
OBJS = $(SRCS:.cpp=.o)

# Default target
all: $(TARGET)

# Link the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

# Compile individual object files
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up generated files
clean:
	rm -f $(OBJS) $(TARGET)

# Phony targets
.PHONY: all clean
