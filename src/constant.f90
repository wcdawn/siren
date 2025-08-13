module constant
use kind
implicit none

private

public :: pi

real(rk), parameter :: pi = 4.0_rk * atan(1.0_rk)

endmodule constant
