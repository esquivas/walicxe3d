# ==============================================================================
# Walicxe 3D Makefile
# ==============================================================================

# Name of the compiled binary
PROGRAM= Walicxe3D

# Non-MPI Compiler
# Supported options: gfortran, ifort
# For MPI, mpfi90 is used
#COMPILER= ifort
COMPILER= gfortran


#  Search Path for Makefile, object files are searched in the order below
VPATH = . : ../source

# Additional user compiler flags
## ifort-compatible flags
#USER_FLAGS= -O2 -warn nounused -nogen-interfaces
#USER_FLAGS= -traceback -warn all -check all,noarg_temp_created -nogen-interfaces
## gfortran-compatible flags
USER_FLAGS= -Wall -pedantic -fbounds-check -g -fbacktrace
#USER_FLAGS= -g -fbacktrace
## Generic flags
USER_FLAGS= -O3

# ============================================= #
# Compilation Time Parameters (Y=on, N=off)     #
# Warning: the following is all case sensitive! #
# ============================================= #

# Use MPI parallelization?
MPI= Y

# Use double precision for all real variables?
DOUBLEP= Y

# Enable passive magnetic (passive *or* active)
BFIELD  = Y

# ==============================================================================

# Additional USER modules to compile
MODULES_USER= brio-wu.o

# Independently compilable modules in the MAIN CODE
MODULES_MAIN= \
constants.o \
parameters.o \
globals.o \
tictoc.o \
clean_quit.o \
hydro_core.o \
utils.o \
hilbert.o \
amr.o

# Dependent source files and modules in the MAIN CODE
OBJECTS_MAIN= \
loadbalance.o \
uniformISM.o  \
snr.o \
winds.o \
orbits.o \
user.o \
cooling.o \
init.o \
output.o \
boundary.o \
lax.o \
flux_cd.o \
hll.o \
hllc.o \
hlle.o \
hlld.o \
sources.o \
godunov.o \
hydro_solver.o \
report.o  \
main.o

# List of modules and objects to compile the Column Density facility
OBJECTS_COLDENS= \
constants.o \
utils.o \
coldens.o

# List of modules and objects to compile the Data Extractor facility
OBJECTS_EXTRACT= \
utils.o \
extract.o

# ==============================================================================
# THERE SHOULD BE NO NEED TO MODIFY BELOW THIS LINE
# ==============================================================================

# Build compiler flags
CFLAGS = $(USER_FLAGS) -cpp

# MPI
ifeq ($(MPI),Y)
CFLAGS += -DMPIP
endif

# Double precision (compiler dependent)
ifeq ($(DOUBLEP),Y)
CFLAGS += -DDOUBLEP
ifeq ($(COMPILER),ifort)
CFLAGS += -r8 -Difort
endif
ifeq ($(COMPILER),gfortran)
CFLAGS += -fdefault-real-8 -Dgfortran
endif
endif

# MPI
ifeq ($(BFIELD),Y)
CFLAGS += -DBFIELD
endif

# Set mpif90 as compiler if MPI, otherwise use specified
ifeq ($(MPI),Y)
COMPILER = mpif90
endif

# Join object lists
OBJECTS_ALL = ${MODULES_MAIN} ${MODULES_USER} ${OBJECTS_MAIN}

# ==============================================================================

$(PROGRAM) : prebuild ${OBJECTS_ALL}
	@echo Linking object files ...
	@echo Linking flags used: $(LINKFLAGS)
	@$(COMPILER) $(CFLAGS) $(OBJECTS_ALL) -o $@
	@echo Cleaning up ...
	@rm -f *.o *.mod
	@echo "Done! (`date`)"

coldens : prebuild $(OBJECTS_COLDENS)
	@echo Linking object files ...
	@$(COMPILER) $(CFLAGS) $(OBJECTS_COLDENS) -o $@
	@echo Cleaning up ...
	@rm -f *.o *.mod
	@echo "Done! (`date`)"

extract : prebuild $(OBJECTS_EXTRACT)
	@echo Linking object files ...
	@$(COMPILER) $(CFLAGS) $(OBJECTS_EXTRACT) -o $@
	@echo Cleaning up ...
	@rm -f *.o *.mod
	@echo "Done! (`date`)"

prebuild :
	@echo "Walicxe3D build started `date`"
	@echo Compiling flags used: $(CFLAGS)
%.o : %.f90
	@echo Compiling $^ ...
	@$(COMPILER) $(CFLAGS) -c $^ -o $@

clean :
	rm -f $(PROGRAM) *.o *.mod

cleanall :
	rm -f *.o *.mod
	rm -f $(PROGRAM)
	rm -f $(PROGRAM).*
	rm -f coldens
	rm -f extract
	rm -f DATA/*.bin
	rm -f DATA/*.vtk
	rm -f DATA/*.dat
	rm -f DATA/*.log
	rm -f DATA/*.visit
