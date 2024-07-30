# Define the compiler
FC = gfortran

# Define the compiler flags
FFLAGS = -Jrelease -O3 -flto -Irelease

# Define the name of the executable
EXE = main

# Source directory to access the source files
SRCDIR = src

# Directory to put all the files and executables
CONFIG = release

# Define the source files with path
BARE_SRCS = parameters.f90 mod_global.f90 marsaglia.f90 commands.f90 mod_io.f90 boltzmann.f90 moldynamics.f90 coulomb_explosion.f90 main.f90

SRCS = $(addprefix $(SRCDIR)/,$(BARE_SRCS))

# Define the object files with path
BARE_OBJS = $(SRCS:$(SRCDIR)/%.f90=%.o)
OBJS = $(addprefix $(CONFIG)/,$(BARE_OBJS))

# Define the module files with path
BARE_MODS = $(SRCS:$(SRCDIR)/%.f90=%.o)
MODS = $(addprefix $(CONFIG)/,$(BARE_OBJS))

# The default rule: build the executable
all: $(CONFIG)/$(EXE)

# Rule to compile the program
$(CONFIG)/parameters.o: $(SRCDIR)/parameters.f90
	@mkdir -p $(CONFIG)
	$(FC) $(FFLAGS) -c $(SRCDIR)/parameters.f90 -o $(CONFIG)/parameters.o

# Rule to compile the program
$(CONFIG)/mod_global.o: $(SRCDIR)/mod_global.f90
	$(FC) $(FFLAGS) -c $(SRCDIR)/mod_global.f90 -o $(CONFIG)/mod_global.o

# Rule to compile the program
$(CONFIG)/marsaglia.o: $(SRCDIR)/marsaglia.f90
	$(FC) $(FFLAGS) -c $(SRCDIR)/marsaglia.f90 -o $(CONFIG)/marsaglia.o

# Rule to compile the program
$(CONFIG)/commands.o: $(SRCDIR)/commands.f90
	$(FC) $(FFLAGS) -c $(SRCDIR)/commands.f90 -o $(CONFIG)/commands.o

# Rule to compile the program
$(CONFIG)/mod_io.o: $(SRCDIR)/mod_io.f90
	$(FC) $(FFLAGS) -c $(SRCDIR)/mod_io.f90 -o $(CONFIG)/mod_io.o

# Rule to compile the program
$(CONFIG)/boltzmann.o: $(SRCDIR)/boltzmann.f90
	$(FC) $(FFLAGS) -c $(SRCDIR)/boltzmann.f90 -o $(CONFIG)/boltzmann.o

# Rule to compile the program
$(CONFIG)/moldynamics.o: $(SRCDIR)/moldynamics.f90
	$(FC) $(FFLAGS) -c $(SRCDIR)/moldynamics.f90 -o $(CONFIG)/moldynamics.o

# Rule to compile the program
$(CONFIG)/coulomb_explosion.o: $(SRCDIR)/coulomb_explosion.f90
	$(FC) $(FFLAGS) -c $(SRCDIR)/coulomb_explosion.f90 -o $(CONFIG)/coulomb_explosion.o

# Rule to compile the program
$(CONFIG)/main.o: $(SRCDIR)/main.f90
	$(FC) $(FFLAGS) -c $(SRCDIR)/main.f90 -o $(CONFIG)/main.o


# Rule to link object files to create the executable
$(CONFIG)/$(EXE): $(OBJS)
	$(FC) $(LINKFLAGS) -o $(CONFIG)/$(EXE) $(OBJS) 


# Rule to clean up the directory
clean:
	rm -f release/*

cleaner:
	rm -rf release

cleanest:
	rm -rf release
	rm -rf output

clean_output:
	rm -rf output