# Compiler
CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++23

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

# Compile .cpp into .o
%.o: %.cpp puzzle.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean build files
clean:
	rm -f $(OBJS) $(TARGET)

# Run program
run: $(TARGET)
	./$(TARGET)
