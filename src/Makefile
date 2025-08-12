# Makefile

FC = gfortran
FFLAGS = -O3 -Wall -Wextra -std=f2008 -pedantic
#FFLAGS = -O0 -g -fbacktrace -ffpe-trap=invalid,overflow,zero -fcheck=all -Wall -Wextra -std=f2008 -pedantic

OBJ = \
	  analytic.o \
	  constant.o \
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

analytic.o: kind.o
constant.o : kind.o
diffusion.o : kind.o linalg.o output.o xs.o
fileio.o :
geometry.o : kind.o
input.o : fileio.o kind.o output.o
kind.o :
linalg.o : kind.o
main.o : xs.o input.o geometry.o output.o power.o transport.o
output.o : fileio.o kind.o
power.o : kind.o xs.o
transport.o : kind.o linalg.o xs.o
xs.o : fileio.o kind.o output.o
