module linalg
use kind
implicit none

private

public :: trid, norm

contains

  real(rk) function norm(ell, x)
    integer(ik), intent(in) :: ell
    real(rk), intent(in) :: x(:)
    integer(ik) :: i
    real(rk) :: xsum
    select case(ell)
      case (-1) ! infinity norm
        norm = maxval(abs(x))
      case (1)
        norm = sum(abs(x))
      case (2)
        norm = sqrt(sum(abs(x)*abs(x)))
      case default
        xsum = 0_rk
        do i = 1,size(x)
          xsum = xsum + abs(x(i))**ell
        enddo ! i = 1,size(x)
        xsum = xsum**(1_rk/ell)
        norm = xsum
    endselect
  endfunction norm

  subroutine trid(n, sub, dia, sup, b, x)
    integer(ik), intent(in) :: n
    real(rk), intent(inout) :: sub(:), dia(:), sup(:)
    real(rk), intent(inout) :: b(:)
    real(rk), intent(inout) :: x(:)

    integer :: i
    real(rk) :: w

    do i = 2,n
      w = sub(i-1)/dia(i-1)
      dia(i) = dia(i) - w*sup(i-1)
      b(i) = b(i) - w*b(i-1)
    enddo

    x(n) = b(n)/dia(n)
    do i = n-1,1,-1
      x(i) = (b(i) - sup(i)*x(i+1))/dia(i)
    enddo

    return
  endsubroutine trid

endmodule linalg
