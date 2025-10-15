module kind
use iso_fortran_env, only : real64, real128, int32, int64
implicit none (external)

private

public :: rk, ik

integer, parameter :: rk = real64
integer, parameter :: ik = int32

endmodule kind
