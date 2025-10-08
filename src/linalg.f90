module linalg
use kind
implicit none

private

public :: trid, norm, trid_block

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

  ! block tri-diagonal system of n blocks, each nprime*nprime
  subroutine trid_block(n, nprime, sub, dia, sup, b, x)
    integer(ik), intent(in) :: n, nprime
    real(rk), intent(inout) :: sub(:,:,:), dia(:,:,:), sup(:,:,:) ! (nprime, nprime, n)
    real(rk), intent(inout) :: b(:,:) ! (nprime,n)
    real(rk), intent(inout) :: x(:,:) ! (nprime,n)

    integer(ik) :: i
    integer :: info

    real(rk), allocatable :: w(:,:) ! (nprime,nprime)
    real(rk), allocatable :: y(:) ! (nprime)

    ! scratch space used by LAPACK
    integer, allocatable :: ipiv(:)
    real(rk), allocatable :: work(:)

    ! TODO there is a trade to be made here...
    ! I think that it would kill my memory if I had to actually store all of these
    ! inverses, but what do I know.
    ! I guess fundamentally, this is probably where it's time to move to a
    ! "matrix-free" method of evaluation, because at this point I have already stored
    ! the nprime*nprime*n*3 numbers...
    real(rk), allocatable :: inv(:,:)

    allocate(w(nprime,nprime))
    allocate(inv(nprime,nprime))

    allocate(ipiv(nprime))
    allocate(work(nprime))

    allocate(inv(nprime,nprime))

    do i = 2,n

      inv = dia(:,:,i)
      call dgetrf(nprime, nprime, inv, ipiv, info)
      if (info /= 0) then
        stop 'failure from dgetrf'
      endif
      call dgetri(nprime, inv, nprime, ipiv, work, n, info)
      if (info /= 0) then
        stop 'failure from dgetri'
      endif

      w = matmul(sub(:,:,i-1), inv)
      dia(:,:,i) = dia(:,:,i) - matmul(w, sup(:,:,i-1))
      b(:,i) = b(:,i) - matmul(w, b(:,i-1))
    enddo

    inv = dia(:,:,n)
    call dgetrf(nprime, nprime, inv, ipiv, info)
    if (info /= 0) then
      stop 'failure from dgetrf'
    endif
    call dgetri(nprime, inv, nprime, ipiv, work, n, info)
    if (info /= 0) then
      stop 'failure from dgetri'
    endif

    x(:,n) = matmul(inv, b(:,n))
    do i = n-1,1,-1

      inv = dia(:,:,i)
      call dgetrf(nprime, nprime, inv, ipiv, info)
      if (info /= 0) then
        stop 'failure from dgetrf'
      endif
      call dgetri(nprime, inv, nprime, ipiv, work, n, info)
      if (info /= 0) then
        stop 'failure from dgetri'
      endif

      y = matmul(sup(:,:,i), x(:,i+1))
      y = b(:,i) - y
      x(:,i) = matmul(inv, y)
    enddo

    deallocate(w, inv)
    deallocate(ipiv, work)
    deallocate(inv)
  endsubroutine trid_block

endmodule linalg
