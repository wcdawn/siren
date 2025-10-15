module linalg
use kind, only : rk, ik
implicit none (external)

private

public :: trid, norm, trid_block, inv, solve, geneig

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

    external :: dgetrf
    external :: dgetri

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

  subroutine solve(n, a, b, x)
    integer(ik), intent(in) :: n
    real(rk), intent(in) :: a(:,:) ! (n,n)
    real(rk), intent(in) :: b(:) ! (n)
    real(rk), intent(out) :: x(:) ! (n)

    integer :: info
    integer(ik), allocatable :: ipiv(:)

    real(rk), allocatable :: acpy(:,:)
    real(rk), allocatable :: bcpy(:)

    external :: dgesv

    allocate(acpy(n,n))
    allocate(bcpy(n))
    allocate(ipiv(n))

    acpy = a
    bcpy = b

    call dgesv(n, 1, acpy, n, ipiv, bcpy, n, info)
    if (info /= 0) then
      stop 'failure from dgesv'
    endif
    x = bcpy

    deallocate(ipiv)
    deallocate(acpy)
    deallocate(bcpy)
  endsubroutine solve

  ! block tri-diagonal system of n blocks, each nprime*nprime
  ! credit to textbook "Computational Mathematics: Models, Methods, and Analysis with MATLAB and MPI"
  ! by R. E. White
  ! This StackOverflow post was also helpful
  ! https://scicomp.stackexchange.com/questions/6701/how-to-solve-block-tridiagonal-matrix-using-thomas-algorithm
  subroutine trid_block(n, nprime, sub, dia, sup, b, x)
    integer(ik), intent(in) :: n, nprime
    real(rk), intent(in) :: sub(:,:,:), dia(:,:,:), sup(:,:,:) ! (nprime, nprime, n)
    real(rk), intent(in) :: b(:,:) ! (nprime,n)
    real(rk), intent(out) :: x(:,:) ! (nprime,n)

    integer(ik) :: i

    real(rk), allocatable :: gmat(:,:,:) ! (nprime, nprime, n)
    real(rk), allocatable :: yvec(:,:) ! (nprime, n)

    real(rk), allocatable :: invmat(:,:), auxmat(:,:) ! (nprime, nprime)

    allocate(gmat(nprime,nprime,n))
    allocate(yvec(nprime,n))
    allocate(invmat(nprime,nprime), auxmat(nprime,nprime))

    ! compute values of gmat
    auxmat = dia(:,:,1)
    call inv(nprime, auxmat, invmat)
    gmat(:,:,1) = matmul(invmat, sup(:,:,1))
    yvec(:,1) = matmul(invmat, b(:,1))
    do i = 2,n-1
      auxmat = dia(:,:,i) - matmul(sub(:,:,i-1), gmat(:,:,i-1))
      call inv(nprime, auxmat, invmat)
      gmat(:,:,i) = matmul(invmat, sup(:,:,i))
      yvec(:,i) = matmul(invmat, b(:,i) - matmul(sub(:,:,i-1), yvec(:,i-1)))
    enddo
    ! unroll for i = n
    auxmat = dia(:,:,n) - matmul(sub(:,:,n-1), gmat(:,:,n-1))
    call inv(nprime, auxmat, invmat)
    yvec(:,n) = matmul(invmat, b(:,n) - matmul(sub(:,:,n-1), yvec(:,n-1)))

    ! backward substitution
    x(:,n) = yvec(:,n)
    do i = n-1,1,-1
      x(:,i) = yvec(:,i) - matmul(gmat(:,:,i), x(:,i+1))
    enddo

    deallocate(invmat, auxmat)
    deallocate(yvec)
    deallocate(gmat)
  endsubroutine trid_block

  subroutine geneig(n, a, b, eigval, eigvec)
    integer(ik), intent(in) :: n
    real(rk), intent(in) :: a(:,:), b(:,:)
    complex(rk), intent(out) :: eigval(:)
    real(rk), intent(out) :: eigvec(:,:)

    real(rk), allocatable :: acpy(:,:), bcpy(:,:)
    real(rk), allocatable :: alphar(:), alphai(:)
    real(rk), allocatable :: beta(:)
    real(rk), allocatable :: vr(:,:), vl(:,:)

    integer(ik) :: lwork
    integer(ik), parameter :: lwork_factor = 10
    real(rk), allocatable :: work(:)
    integer :: info

    integer(ik) :: i

    external :: dggev

    allocate(acpy(n,n), bcpy(n,n))
    allocate(alphar(n), alphai(n))
    allocate(beta(n))
    allocate(vr(n,n), vl(n,n))
    acpy = a
    bcpy = b

    lwork = lwork_factor * n
    allocate(work(lwork))

    call dggev('N', 'V', n, acpy, n, bcpy, n, alphar, alphai, beta, &
      vl, n, vr, n, work, lwork, info)
    if (info /= 0) then
      stop 'failure in dgeev'
    endif

    do i = 1,n
      if (abs(beta(i)) > epsilon(1.0_rk)) then
        eigval(i) = cmplx(alphar(i), alphai(i), rk)/beta(i)
      endif
    enddo
    eigvec = vr

    deallocate(work)
    deallocate(vr, vl)
    deallocate(beta)
    deallocate(alphar, alphai)
    deallocate(acpy, bcpy)
  endsubroutine geneig

endmodule linalg
