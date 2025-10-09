module linalg
use kind
implicit none

private

public :: trid, norm, trid_block, inv

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
  endsubroutine trid

  subroutine inv(n, a, ainv)
    integer(ik), intent(in) :: n
    real(rk), intent(in) :: a(:,:) ! (n,n)
    real(rk), intent(out) :: ainv(:,:) ! (n,n)

    integer :: info

    ! scratch space used by LAPACK
    integer, allocatable :: ipiv(:)
    real(rk), allocatable :: work(:)

    allocate(ipiv(n))
    allocate(work(n))

    ainv = a
    call dgetrf(n, n, ainv, n, ipiv, info)
    if (info /= 0) then
      stop 'failure from dgetrf'
    endif
    call dgetri(n, ainv, n, ipiv, work, n, info)
    if (info /= 0) then
      stop 'failure from dgetri'
    endif

    deallocate(ipiv, work)
  endsubroutine inv

  ! block tri-diagonal system of n blocks, each nprime*nprime
  subroutine trid_block(n, nprime, sub, dia, sup, b, x)
    integer(ik), intent(in) :: n, nprime
    real(rk), intent(inout) :: sub(:,:,:), dia(:,:,:), sup(:,:,:) ! (nprime, nprime, n)
    real(rk), intent(inout) :: b(:,:) ! (nprime,n)
    real(rk), intent(inout) :: x(:,:) ! (nprime,n)

    integer(ik) :: i

    real(rk), allocatable :: w(:,:) ! (nprime,nprime)
    real(rk), allocatable :: y(:) ! (nprime)

    ! TODO there is a trade to be made here...
    ! I think that it would kill my memory if I had to actually store all of these
    ! inverses, but what do I know.
    ! I guess fundamentally, this is probably where it's time to move to a
    ! "matrix-free" method of evaluation, because at this point I have already stored
    ! the nprime*nprime*n*3 numbers...
    real(rk), allocatable :: a(:,:)

    allocate(w(nprime,nprime))
    allocate(y(nprime))
    allocate(a(nprime,nprime))

    do i = 2,n
      call inv(nprime, dia(:,:,i-1), a)
      w = matmul(sub(:,:,i-1), a)
      dia(:,:,i) = dia(:,:,i) - matmul(w, sup(:,:,i-1))
      b(:,i) = b(:,i) - matmul(w, b(:,i-1))
    enddo

    call inv(nprime, dia(:,:,n), a)
    x(:,n) = matmul(a, b(:,n))
    do i = n-1,1,-1
      call inv(nprime, dia(:,:,i), a)
      y = matmul(sup(:,:,i), x(:,i+1))
      y = b(:,i) - y
      x(:,i) = matmul(a, y)
    enddo

    deallocate(w, a, y)
  endsubroutine trid_block

endmodule linalg
