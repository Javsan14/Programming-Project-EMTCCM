# MAKEFILE FOR PROGRAMMING_PROJECT.F90
# Compiler and flags
FC = gfortran
FLAGS = -g -O2 -Wall
LIBS = -llapack

# Files
SOURCES = programming_project.f90
OBJECTS = $(SOURCES:.f90=.o)
EXECUTABLE = program

# Default rule
all: $(EXECUTABLE)

# Link the program
$(EXECUTABLE): $(OBJECTS)
	$(FC) $(FLAGS) -o $@ $(OBJECTS) $(LIBS)

# Compile source files
%.o: %.f90
	$(FC) $(FLAGS) -c $< -o $@

# Clean up
clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

# Rebuild
rebuild: clean all
