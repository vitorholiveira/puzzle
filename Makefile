# Compiler
CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++23 -O3 -march=native

# Executable name
TARGET = main

# Source files
SRCS = main.cpp puzzle.cpp
OBJS = $(SRCS:.cpp=.o)

# Default rule
all: $(TARGET)

# Link object files to create executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile .cpp into .o with dependency generation
%.o: %.cpp puzzle.h
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

# Include dependencies
-include $(OBJS:.o=.d)

# Clean build files
clean:
	rm -f $(OBJS) $(TARGET) $(OBJS:.o=.d)

# Run program
run: $(TARGET)
	./$(TARGET)
