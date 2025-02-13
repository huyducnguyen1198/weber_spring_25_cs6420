# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -Wall -Wextra -O2 -std=c++17

# Target executables
TARGET = rabin_karp
KMP_TARGET = kmp

# Source files for each target
SRCS = rabin_karp.cpp
KMP_SRCS = KMP.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)
KMP_OBJS = $(KMP_SRCS:.cpp=.o)

# Default rule: builds both executables
all: $(TARGET) $(KMP_TARGET)

# Rule to build only Rabin-Karp
rk: $(TARGET)  # Now correctly referencing rabin_karp

# Rule to build only KMP
kmp: $(KMP_TARGET)  # Now correctly referencing kmp

# Linking the Rabin-Karp executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Linking the KMP executable
$(KMP_TARGET): $(KMP_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compiling C++ source files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJS) $(KMP_OBJS) $(TARGET) $(KMP_TARGET)

# Phony targets
.PHONY: all clean kmp rk