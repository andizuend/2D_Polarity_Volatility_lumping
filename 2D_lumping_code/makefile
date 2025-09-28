## makefile for 2D_Pol_Vol_lumping framework with 	##
## source files located in different subfolders 	##

# compiler settings:
FC = gfortran
FCFLAGS = -O3 -Wall -g -fcheck=all -ffree-line-length-none -Wno-uninitialized -Wno-unused

# linker settings:
LDFLAGS = 

# executable program file name:
program = 2D_Pol_Vol.out

# source files and their relative paths are listed 
# in src_list.txt in order accounting for dependencies:
src_list := ./src_list.txt
srcs := $(file < $(src_list))

# object files:
objs = $(patsubst %.f90, %.o, $(srcs))
this_folder = ../2D_lumping_code


all: $(program)

$(program): $(objs)
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $< -o $(this_folder)/$@

.phony: clean
clean:
	rm -f *.mod *.smod $(this_folder)/*.o $(this_folder)/Fortran_code/*.o $(this_folder)/Fortran_code/Fortran_OS_commands/*.o $(this_folder)/Fortran_code/FortranCode_Lumping_Schemes/*.o $(program)