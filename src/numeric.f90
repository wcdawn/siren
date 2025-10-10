module numeric
use kind
implicit none

private
public :: deriv

contains

  ! return first derivative using a second-order estimate on non-uniform grid
  ! derivative is returned at coordinate x2
  ! this function also works in the case of a non-uniform grid
  !
  ! f1    f2    f3
  ! |-----|-----|
  ! x1    x2    x3
  ! h_left|h_right
  !
  ! See:
  ! - numpy.gradient https://numpy.org/doc/stable/reference/generated/numpy.gradient.html
  ! - Forgnberg, Bengt. "Generation of Finite Difference Formulas on Aribtrarily Spaced Grids."
  !   Mathematics of Computation. Volume 51, Number 184. October 1988.
  real(rk) pure function deriv(x1, x2, x3, f1, f2, f3)
    real(rk), intent(in) :: x1, x2, x3 ! x-coordinate
    real(rk), intent(in) :: f1, f2, f3 ! function value

    real(rk) :: h_left, h_right
    real(rk) :: a, b, c

    h_left = x2 - x1
    h_right = x3 - x2

    a = -h_right**2  / (h_left * h_right * (h_left + h_right))
    c = h_left**2  / (h_left * h_right * (h_left + h_right))
    b = -(a+c)

    deriv = a * f1 + b * f2 + c * f3
  endfunction deriv

endmodule numeric
