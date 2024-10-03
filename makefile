# Makefile

# Compiler and flags
CC = g++
CFLAGS = -std=c++17 -O3 -fPIC -MMD -MP
LDFLAGS = -shared

# Target library
TARGET = libmc.so

# Source files
SOURCES = MC_front.cpp MC.cpp BOX.cpp Objects.cpp

# Object files
OBJECTS = $(SOURCES:.cpp=.o)

# Default target
all: $(TARGET)

# Rule to create the shared library
$(TARGET): $(OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $(OBJECTS)

# Rule to compile .cpp files into .o files
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

# Clean up generated files
clean:
	rm -f $(OBJECTS) $(SOURCES:.cpp=.d) $(TARGET)

# Include dependency files
-include $(SOURCES:.cpp=.d)
