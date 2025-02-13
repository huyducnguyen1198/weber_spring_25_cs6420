# Compiler
CXX = g++

# Compiler Flags
CXXFLAGS = -std=c++17 -Wall -Wextra -O2

# Target Executable
TARGET = best_case_quick

# Source File
SRC = best_case_quick.cpp

# Build Target
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

# Clean Rule
clean:
	rm -f $(TARGET)

# PHONY Targets
.PHONY: clean