# Makefile

FC = gfortran
FFLAGS = -O3 -Wall -Wextra -std=f2008 -pedantic
#FFLAGS = -O0 -g -fbacktrace -ffpe-trap=invalid,overflow,zero -fcheck=all -Wall -Wextra -std=f2008 -pedantic

OBJ = \
	  analytic.o \
	  constant.o \
	  diffusion.o \
	  diffusion_block.o \
	  exception_handler.o \
	  fileio.o \
	  geometry.o \
	  input.o \
	  kind.o \
	  linalg.o \
	  main.o \
	  output.o \
	  power.o \
	  timer.o \
	  transport.o \
	  transport_block.o \
	  xs.o

MOD = $(OBJ:o=mod)

LIBS = -llapack -lblas

.PHONY : default all clean

default : all

all : siren.x

siren.x : $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	$(RM) $(OBJ)
	$(RM) $(MOD)
	$(RM) siren.x

analytic.o: constant.o exception_handler.o fileio.o kind.o linalg.o output.o xs.o
constant.o : kind.o
diffusion.o : exception_handler.o kind.o linalg.o output.o timer.o xs.o
diffusion_block.o : exception_handler.o kind.o linalg.o output.o timer.o xs.o
exception_handler.o : kind.o output.o
fileio.o :
geometry.o : output.o kind.o
input.o : exception_handler.o fileio.o kind.o output.o
kind.o :
linalg.o : kind.o
main.o : kind.o xs.o input.o geometry.o diffusion.o diffusion_block.o transport.o transport_block.o output.o power.o analytic.o exception_handler.o timer.o
output.o : fileio.o kind.o
power.o : kind.o xs.o
timer.o : kind.o output.o
transport.o : exception_handler.o kind.o linalg.o output.o timer.o xs.o
transport_block.o : exception_handler.o kind.o linalg.o output.o timer.o xs.o
xs.o : exception_handler.o fileio.o kind.o output.o
