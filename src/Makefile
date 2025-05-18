# Makefile

FC = gfortran
#FFLAGS = -O3 -Wall -Wextra -std=f2008 -pedantic
FFLAGS = -O0 -g -fbacktrace -ffpe-trap=invalid,overflow,zero -fcheck=all -Wall -Wextra -std=f2008 -pedantic

OBJ = \
	  diffusion.o \
	  fileio.o \
	  geometry.o \
	  input.o \
	  kind.o \
	  linalg.o \
	  main.o \
	  output.o \
	  power.o \
	  transport.o \
	  xs.o

.PHONY : default all clean

default : all

all : siren.x

siren.x : $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	$(RM) $(OBJ)
	$(RM) siren.x

diffusion.o : kind.o linalg.o xs.o
fileio.o :
geometry.o : kind.o
input.o : kind.o fileio.o
kind.o :
linalg.o : kind.o
main.o : xs.o input.o geometry.o output.o power.o transport.o
output.o : fileio.o kind.o
power.o : kind.o xs.o
transport.o : kind.o linalg.o xs.o
xs.o : kind.o fileio.o
