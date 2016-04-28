######################################################
# Makefile for Fortran programs
#  Author: Elmo Tempel
#  Last chaged: (28.04.2016)
#####################################################

# increase stack size: ulimit -s 65532
# ulimit -s shows the current stack size

# Fortran source files and program name:
PROG = bisous_program
# object files
OBJ = program.o constants.o parameters.o cylinders.o cubesort.o quicksort.o \
      galaxies.o utilities.o data_term.o interaction_term.o config.o \
	  operations.o mcmc.o simulation.o preparation.o read_cylinders.o \
	  post_processing.o fil_catalogue.o # fits_op.o
# Variable declaration (Fortran Compilator):
FC = ifort # for intel fortran compiler
#FC = gfortran  # for gfortran compiler
OPT = -O3 #-ffree-line-length-none # -ffree-line-length-none should be commented out for gfortran
DBG = -g -O0 -error-limit 10 -warn all
OMP = #-openmp
LIBS =  #-lcfitsio -lrecipes_f90
# choose flags for compiler
FLAGS = $(OPT) $(OMP) $(LIBS)

# Variable declaration (system directories):
DIR_INSTALL = ./
DIR_SRC = ./src

all: $(PROG)

$(OBJ): %.o: $(DIR_SRC)/%.f90
	$(FC) -c $< $(FLAGS)

$(PROG): $(OBJ)
	$(FC) $(OBJ) -o $(DIR_INSTALL)/$(PROG) $(FLAGS)

.PHONY: clean
clean:
	rm -f *.o *.mod

# define the dependencies
program.o : constants.o parameters.o cylinders.o cubesort.o galaxies.o data_term.o \
            interaction_term.o operations.o mcmc.o simulation.o preparation.o \
			read_cylinders.o post_processing.o fil_catalogue.o
cylinders.o : constants.o parameters.o
galaxies.o : constants.o parameters.o cubesort.o utilities.o
data_term.o : constants.o parameters.o cylinders.o galaxies.o utilities.o
interaction_term.o : constants.o parameters.o cylinders.o galaxies.o
operations.o : constants.o cylinders.o galaxies.o data_term.o interaction_term.o \
               parameters.o cubesort.o
mcmc.o : constants.o parameters.o galaxies.o cylinders.o data_term.o interaction_term.o \
         operations.o
simulation.o : constants.o parameters.o galaxies.o cylinders.o data_term.o interaction_term.o \
         operations.o mcmc.o utilities.o preparation.o
preparation.o : constants.o parameters.o galaxies.o cylinders.o data_term.o interaction_term.o \
         operations.o utilities.o
parameters.o : constants.o config.o
config.o : constants.o
utilities.o : constants.o
cubesort.o : constants.o quicksort.o
quicksort.o : constants.o
constants.o :
#fits_op.o : constants.o
read_cylinders.o : constants.o cylinders.o utilities.o operations.o
post_processing.o : constants.o parameters.o cylinders.o galaxies.o cubesort.o read_cylinders.o \
                    utilities.o fil_catalogue.o # fits_op.o
fil_catalogue.o : constants.o parameters.o cylinders.o galaxies.o cubesort.o









