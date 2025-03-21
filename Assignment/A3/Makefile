# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -Wall -Wextra -O2 -std=c++17

# Target executables
TARGET = rabin_karp
KMP_TARGET = kmp  # Renamed from 'kmp' to avoid conflicts

# Source files for each target
SRCS = rabin_karp.cpp
KMP_SRCS = kmp.cpp

# Object files
OBJS = rabin_karp.o
KMP_OBJS = kmp.o

# Default rule: builds both executables
all: $(TARGET) $(KMP_TARGET)

# Rule to build only Rabin-Karp
rk: $(TARGET)

# Rule to build only KMP
kmp: $(KMP_TARGET)

# Linking the Rabin-Karp executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Linking the KMP executable (renamed)
$(KMP_TARGET): $(KMP_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Explicit compilation rules
rabin_karp.o: rabin_karp.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

kmp.o: kmp.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJS) $(KMP_OBJS) $(TARGET) $(KMP_TARGET)

# Phony targets
.PHONY: all clean kmp rk