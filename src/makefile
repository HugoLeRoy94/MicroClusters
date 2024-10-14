# Makefile

# Compiler and flags
CC = g++
CFLAGS = -std=c++17 -O3 -MMD -MP
SHARED_CFLAGS = $(CFLAGS) -fPIC
LDFLAGS = -shared

# Target library
TARGET = libmc.so

# Source files
SOURCES = MC_front.cpp MC.cpp BOX.cpp Objects.cpp Move.cpp

# Object files for the shared library
SHARED_OBJECTS = $(SOURCES:.cpp=.so.o)

# Object files for the executable
OBJECTS = $(SOURCES:.cpp=.o)

# Default target
all: $(TARGET)

# Rule to create the shared library
$(TARGET): $(SHARED_OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $(SHARED_OBJECTS)

# Rule to compile .cpp files into .so.o files for the shared library
%.so.o: %.cpp
	$(CC) $(SHARED_CFLAGS) -c $< -o $@

# Rule to compile .cpp files into .o files for the executable
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

# Rule to compile main.cpp into main.o
main.o: main.cpp
	$(CC) $(CFLAGS) -c $< -o $@

# Executable target
exec: simulation

simulation: main.o $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^

# Clean up generated files
clean:
	rm -f $(OBJECTS) $(SHARED_OBJECTS) main.o main.d $(SOURCES:.cpp=.d) $(TARGET) simulation

# Include dependency files
-include $(SOURCES:.cpp=.d) main.d
