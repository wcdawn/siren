module linalg
use kind, only : rk, ik
implicit none

private

public :: trid, norm, trid_block, inv, solve, geneig, &
  trid_sor, trid_conjugate_gradient

contains

  real(rk) function norm(x, ell)
    real(rk), intent(in) :: x(:)
    integer(ik), intent(in) :: ell
    integer(ik) :: i
    real(rk) :: xsum
    select case(ell)
      case (-1) ! infinity norm
        xsum = 0.0_rk
        !$omp parallel do default(none) private(i) shared(x) reduction(max:xsum)
        do i = 1,size(x)
          xsum = max(xsum, abs(x(i)))
        enddo ! i = 1,size(x)
        !$omp end parallel do
        norm = xsum
      case (1)
        xsum = 0.0_rk
        !$omp parallel do default(none) private(i) shared(x) reduction(+:xsum)
        do i = 1,size(x)
          xsum = xsum + abs(x(i))
        enddo
        !$omp end parallel do
        norm = xsum
      case (2)
        xsum = 0.0_rk
        !$omp parallel do default(none) private(i) shared(x) reduction(+:xsum)
        do i = 1,size(x)
          xsum = xsum + abs(x(i))**2
        enddo ! i = 1,size(x)
        !$omp end parallel do
        norm = sqrt(xsum)
      case default
        xsum = 0.0_rk
        do i = 1,size(x)
          xsum = xsum + abs(x(i))**ell
        enddo ! i = 1,size(x)
        norm = xsum**(1_rk/ell)
    endselect
  endfunction norm

  real(rk) function dot(n, x, y)
    integer(ik), intent(in) :: n
    real(rk), intent(in) :: x(:), y(:) ! (n)
    integer(ik) :: i
    real(rk) :: xsum
    xsum = 0.0_rk
    !$omp parallel do default(none) private(i) shared(n, x, y) reduction(+:xsum)
    do i = 1,n
      xsum = xsum + x(i)*y(i)
    enddo
    !$omp end parallel do
    dot = xsum
  endfunction dot

  subroutine copy(n, x, y)
    integer(ik), intent(in) :: n
    real(rk), intent(in) :: x(:) ! (n)
    real(rk), intent(out) :: y(:) ! (n)
    integer(ik) :: i
    !$omp parallel do default(none) private(i) shared(n, x, y)
    do i = 1,n
      y(i) = x(i)
    enddo ! i = 1,n
    !$omp endparallel do
  endsubroutine copy

  subroutine set(n, xx, x)
    integer, intent(in) :: n
    real(rk), intent(in) :: xx
    real(rk), intent(out) :: x(:) ! (n)
    integer(ik) :: i
    !$omp parallel do default(none) private(i) shared(n, x, xx)
    do i = 1,n
      x(i) = xx
    enddo
    !$omp end parallel do
  endsubroutine set

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
    use exception_handler, only : exception_fatal
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
      call exception_fatal('Failure from DGETRF.')
    endif
    call dgetri(n, ainv, n, ipiv, work, n, info)
    if (info /= 0) then
      call exception_fatal('Failure from DGETRI.')
    endif

    deallocate(ipiv, work)
  endsubroutine inv

  subroutine solve(n, a, b, x)
    use exception_handler, only : exception_fatal
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
      call exception_fatal('Failure from DGESV.')
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
    use exception_handler, only : exception_fatal
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
      call exception_fatal('Failure in DGGEV.')
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

  subroutine axpy(n, alpha, x, y, z)
    integer(ik), intent(in) :: n
    real(rk), intent(in) :: alpha
    real(rk), intent(in) :: x(:), y(:) ! (n)
    real(rk), intent(inout) :: z(:) ! (n)
    integer(ik) :: i
    !$omp parallel do default(none) private(i) shared(alpha, x, y, z)
    do i = 1,n
      z(i) = alpha*x(i) + y(i)
    enddo ! i = 1,n
    !$omp end parallel do
  endsubroutine axpy

  subroutine trid_matvec(n, sub, dia, sup, x, z)
    integer(ik), intent(in) :: n
    real(rk), intent(in) :: sub(:) ! (n-1)
    real(rk), intent(in) :: dia(:) ! (n)
    real(rk), intent(in) :: sup(:) ! (n-1)
    real(rk), intent(in) :: x(:) ! (n)
    real(rk), intent(out) :: z(:) ! (n)
    integer(ik) :: i
    if (n == 1) then
      z(1) = dia(1)*x(1)
    else
      z(1) = dia(1)*x(1) + sup(1)*x(2)
      !$omp parallel do default(none) private(i) shared(sub, dia, sup, x, z)
      do i = 2,n-1
        z(i) = sub(i-1)*x(i-1) + dia(i)*x(i) + sup(i)*x(i+1)
      enddo ! i = 1,n
      !$omp end parallel do
      z(n) = sub(n-1)*x(n-1) + dia(n)*x(n)
    endif
  endsubroutine trid_matvec

  subroutine trid_conjugate_gradient(n, sub, dia, sup, b, maxit, atol, rtol, x, verbose)
    use exception_handler, only : exception_fatal
    integer(ik), intent(in) :: n
    real(rk), intent(in) :: sub(:) ! (n-1)
    real(rk), intent(in) :: dia(:) ! (n)
    real(rk), intent(in) :: sup(:) ! (n-1)
    real(rk), intent(in) :: b(:) ! (n)
    integer(ik), intent(in) :: maxit
    real(rk), intent(in) :: atol, rtol
    real(rk), intent(out) :: x(:) ! (n)
    logical, intent(in), optional :: verbose

    integer(ik) :: iter
    real(rk) :: alpha, beta
    real(rk) :: initial_norm, nom

    real(rk), allocatable :: r(:)
    real(rk), allocatable :: w(:), p(:)

    character(1024) :: line

    allocate(r(n))

    call trid_matvec(n, sub, dia, sup, x, r) ! r = A*x
    call axpy(n, -1.0_rk, r, b, r)

    nom = norm(r, 2)
    initial_norm = nom
    if (nom < atol) then
      ! if the initial residual is sufficiently small, then return x0 as the result
      deallocate(r)
      if (present(verbose)) then
        if (verbose) then
          write(*,*) 'quitting early'
        endif
      endif
      return
    endif

    allocate(w(n), p(n))

    call copy(n, r, p)
    
    do iter = 1,maxit
      call trid_matvec(n, sub, dia, sup, p, w)
      alpha = norm(r, 2)**2 / dot(n, p, w)
      call axpy(n, alpha, p, x, x)
      call axpy(n, -alpha, w, r, r)
      nom = norm(r,2)
      if (present(verbose)) then
        if (verbose) then
          write(*,'(a,i0,a,es13.6,a,es13.6)') 'it=', iter, ' norm=', nom, ' improvement=', nom/initial_norm
        endif
      endif
      if ((nom < rtol*initial_norm) .or. (nom < atol)) then
        exit
      endif
      call axpy(n, beta, p, r, p)
    enddo ! iter = 1,max_iter

    if (iter > maxit) then
      write(line, '(a,i0,a,a,es8.1,a,es8.1,a,es8.1,a,es8.1)') &
        'Failed to converge CG after ', maxit, ' iterations.', &
        ' rtol_target=', rtol, ' rtol_achieved=', nom/initial_norm, &
        ' atol_target=', atol, ' atol_achieved', nom
      call exception_fatal(line)
    endif

    deallocate(r)
    deallocate(w, p)
  endsubroutine trid_conjugate_gradient

  subroutine trid_sor(n, sub, dia, sup, b, maxit, rtol, omega, x)
    use exception_handler, only : exception_fatal
    integer(ik), intent(in) :: n
    real(rk), intent(in) :: sub(:) ! (n-1)
    real(rk), intent(in) :: dia(:) ! (n)
    real(rk), intent(in) :: sup(:) ! (n-1)
    real(rk), intent(in) :: b(:) ! (n)
    integer(ik), intent(in) :: maxit
    real(rk), intent(in) :: rtol
    real(rk), intent(in) :: omega
    real(rk), intent(out) :: x(:) ! (n)

    integer(ik) :: i, iter
    real(rk) :: xdif, xmax, xold, xnew

    character(1024) :: line

    xdif = 0.0_rk
    xmax = 0.0_rk

    do iter = 1,maxit
      xdif = 0.0_rk
      xmax = 0.0_rk
      
      xold = x(1)
      xnew = sup(1)*x(2)
      xnew = (b(1) - xnew) / dia(1)
      x(1) = xold + omega * (xnew - xold)
      xdif = max(xdif, abs(x(1) - xold))
      xmax = max(xmax, abs(x(1)))

      ! NOTE: this cannot be parallelized or else it approaches a Jacobi iteration
      ! as the key to the Gauss-Siedel iteration (and thereby SOR) is using the
      ! updated vector from the previous step
      do i = 2,n-1
        xold = x(i)
        xnew = sub(i-1)*x(i-1) + sup(i)*x(i+1)
        xnew = (b(i) - xnew) / dia(i)
        x(i) = xold + omega * (xnew - xold)
        xdif = max(xdif, abs(x(i) - xold))
        xmax = max(xmax, abs(x(i)))
      enddo ! i = 2,n-1

      xold = x(n)
      xnew = sub(n-1)*x(n-1)
      xnew = (b(n) - xnew) / dia(n)
      x(n) = xold + omega * (xnew - xold)
      xdif = max(xdif, abs(x(n) - xold))
      xmax = max(xmax, abs(x(n)))

      if (xdif/xmax < rtol) then
        exit
      endif
    enddo ! iter = 1,maxit

    if (iter > maxit) then
      write(line, '(a,i0,a,a,es8.1,a,es8.1,a,f5.2)') &
        'Failed to converge SOR after ', maxit, ' iterations.', &
        ' target=', rtol, ' achieved=', xdif/xmax, ' omega=', omega
      call exception_fatal(line)
    endif
  endsubroutine trid_sor

endmodule linalg
