module analytic
use kind, only : rk, ik
implicit none

private
public :: analytic_error

! twogroup
real(rk), parameter :: Lx_twogroup = 1e2_rk

! tworeg
real(rk), parameter :: Bf_tworeg = 0.01716859736590981012_rk ! [1/cm] LF=80cm
real(rk), parameter :: Lf_tworeg = 80.0_rk
real(rk), parameter :: Lr_tworeg = 100.0_rk

! analytic p3
real(rk), parameter :: Lx_p3 = 1e2_rk

contains

  subroutine analytic_error(analytic_name, fname, nx, ngroup, pnorder, xslib, dx, phi, keff)
    use xs, only : XSLibrary
    use linalg, only : norm
    use output, only : output_write
    use exception_handler, only : exception_fatal
    character(*), intent(in) :: analytic_name
    character(*), intent(in) :: fname
    integer(ik), intent(in) :: nx, ngroup, pnorder
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(in) :: dx(:)
    real(rk), intent(in) :: phi(:,:,:) ! (nx,ngroup,pnorder+1)
    real(rk), intent(in) :: keff

    ! assume that the coordinate system starts at xleft==0.0
    ! we don't really know in general since all we have are deltas
    ! but we need an exact coordinate for the error analysis
    real(rk), parameter :: x0 = 0.0_rk

    real(rk) :: xleft
    integer(ik) :: i, g, n

    real(rk), allocatable :: x(:)
    real(rk), allocatable :: phi_analysis(:,:,:), phi_exact(:,:,:), phi_diff(:,:,:)
    real(rk) :: phi_zero
    real(rk) :: keff_exact, keff_diff
    real(rk) :: ratio_exact, ratio_siren, ratio_diff
    real(rk), allocatable :: linferr(:,:) ! (ngroup,pnorder+1)

    character(1024) :: line

    allocate(x(nx))
    xleft = x0
    do i = 1,nx
      x(i) = xleft + 0.5d0 * dx(i)
      xleft = xleft + dx(i)
    enddo

    ! work on a copy so that we can change normalization
    allocate(phi_analysis(nx,ngroup,pnorder+1))
    phi_analysis = phi
    ! extrapolate to estimate phi(0) for the first group, scalar flux
    ! this is used to normalize such that phi(0) = 1.0 for first group, scalar flux
    ! this is usually equivalent to phi0 == 1 from the analytic solutions
    phi_zero = -(phi(2,1,1)-phi(1,1,1)) * dx(1) / (dx(1) + dx(2)) + phi(1,1,1)
    phi_analysis = phi_analysis / phi_zero

    ! get phi_exact
    allocate(phi_exact(nx,ngroup,pnorder+1))
    select case (analytic_name)
      case ('onegroup')
        call analytic_fun_onegroup(x, phi_exact)
      case ('twogroup')
        call analytic_fun_twogroup(x, xslib, phi_exact)
      case ('tworeg')
        call analytic_fun_tworeg(x, xslib, phi_exact)
      case ('analytic_p3')
        call analytic_fun_p3(x, xslib, phi_exact)
      case ('analytic_pn')
        call analytic_fun_pn(x, xslib, phi_exact)
      case default
        call exception_fatal('unknown analytic_name: ' // trim(adjustl(analytic_name)))
    endselect

    allocate(phi_diff(nx,ngroup,pnorder+1))
    allocate(linferr(ngroup,pnorder+1))
    phi_diff = phi_exact - phi_analysis
    ! NOTE: using the two-norm would be reasonable.
    ! However, one would need to be careful since the mesh is not uniform.
    do n = 1,pnorder+1
      do g = 1,ngroup
        linferr(g,n) = norm(-1, phi_diff(:,g,n))
      enddo ! g = 1,ngroup
    enddo ! n = 1,pnorder+1

    ! get keff_exact
    select case (analytic_name)
      case ('onegroup')
        keff_exact = analytic_keff_onegroup(xslib)
      case ('twogroup')
        keff_exact = analytic_keff_twogroup(xslib)
      case ('tworeg')
        keff_exact = analytic_keff_tworeg(xslib)
      case ('analytic_p3')
        keff_exact = analytic_keff_p3(xslib)
      case ('analytic_pn')
        call analytic_fun_pn(x, xslib, phi_exact, keff_exact)
      case default
        call exception_fatal('second -- unknown analytic_name: ' // trim(adjustl(analytic_name)))
    endselect
    keff_diff = keff_exact - keff

    call output_write("=== ANALYSIS OF ERROR ===")
    write(line, '(a,f22.20)') 'keff_exact = ', keff_exact
    call output_write(line)
    write(line, '(a,f15.6,a)') 'keff_diff = ', keff_diff*1e5_rk, ' [pcm]'
    call output_write(line)

    if (analytic_name == 'twogroup') then
      ratio_exact = analytic_ratio_twogroup(xslib)
      ratio_siren = phi(1,2,1) / phi(1,1,1)
      ratio_diff = ratio_exact - ratio_siren
      write(line, '(a,f23.20)') 'ratio_exact = ', ratio_exact
      call output_write(line)
      write(line, '(a,f23.20)') 'ratio_siren = ', ratio_siren
      call output_write(line)
      write(line, '(a,es13.6)') 'ratio_diff = ', ratio_diff
      call output_write(line)
    elseif (analytic_name == 'analytic_p3') then
      ratio_exact = analytic_ratio_p3(xslib)
      ratio_siren = phi(1,1,3) / phi(1,1,1)
      ratio_diff = ratio_exact - ratio_siren
      write(line, '(a,f23.20)') 'ratio_exact = ', ratio_exact
      call output_write(line)
      write(line, '(a,f23.20)') 'ratio_siren = ', ratio_siren
      call output_write(line)
      write(line, '(a,es13.6)') 'ratio_diff = ', ratio_diff
      call output_write(line)
    endif

    do n = 1,pnorder+1
      do g = 1,ngroup
        write(line, '(a,i0,a,i0,a,es13.6)') 'linferr_n', n-1, '_g', g, ' = ', linferr(g,n)
        call output_write(line)
      enddo ! g = 1,ngroup
    enddo ! n = 1,pnorder+1

    call output_write('writing ' // trim(adjustl(fname)))
    call analytic_output_csv(fname, nx, ngroup, pnorder, x, phi_analysis, phi_exact, phi_diff)

    call output_write('')

    deallocate(linferr)
    deallocate(phi_analysis, phi_exact, phi_diff)
    deallocate(x)
  endsubroutine analytic_error

  subroutine analytic_output_csv(fname, nx, ngroup, pnorder, x, phi_normalized, phi_exact, phi_diff)
    use fileio, only : fileio_open_write
    character(*), intent(in) :: fname
    integer(ik), intent(in) :: nx, ngroup, pnorder
    real(rk), intent(in) :: x(:)
    real(rk), intent(in) :: phi_normalized(:,:,:)
    real(rk), intent(in) :: phi_exact(:,:,:)
    real(rk), intent(in) :: phi_diff(:,:,:)

    integer(ik) :: i, g, n

    integer, parameter :: iounit = 17

    call fileio_open_write(fname, iounit)

    write(iounit, '(a)', advance='no') 'x [cm]'
    do n = 1,pnorder+1
      do g = 1,ngroup
        write(iounit, '(a,a,i0,a,i0)', advance='no') ',', 'phisiren_n', n-1, '_g', g
      enddo ! g = 1,ngroup
    enddo ! n = 1,pnorder+1
    do n = 1,pnorder+1
      do g = 1,ngroup
        write(iounit, '(a,a,i0,a,i0)', advance='no') ',', 'phiexact_n', n-1, '_g', g
      enddo ! g = 1,ngroup
    enddo ! n = 1,pnorder+1
    do n = 1,pnorder+1
      do g = 1,ngroup
        write(iounit, '(a,a,i0,a,i0)', advance='no') ',', 'phidiff_n', n-1, '_g', g
      enddo ! g = 1,ngroup
    enddo ! n = 1,pnorder+1
    write(iounit,*)

    do i = 1,nx
      write(iounit, '(es23.16)', advance='no') x(i)

      do n = 1,pnorder+1
        do g = 1,ngroup
          write(iounit, '(a,es23.16)', advance='no') ',', phi_normalized(i,g,n)
        enddo ! g = 1,ngroup
      enddo ! n = 1,pnorder+1

      do n = 1,pnorder+1
        do g = 1,ngroup
          write(iounit, '(a,es23.16)', advance='no') ',', phi_exact(i,g,n)
        enddo ! g = 1,ngroup
      enddo ! n = 1,pnorder+1

      do n = 1,pnorder+1
        do g = 1,ngroup
          write(iounit, '(a,es23.16)', advance='no') ',', phi_diff(i,g,n)
        enddo ! g = 1,ngroup
      enddo ! n = 1,pnorder+1

      write(iounit,*)

    enddo ! i = 1,nx

    close(iounit)
  endsubroutine analytic_output_csv

  subroutine analytic_fun_twogroup(x, xslib, exact)
    use xs, only : XSLibrary
    use constant, only : pi
    real(rk), intent(in) :: x(:)
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(out) :: exact(:,:,:)

    real(rk), parameter :: phi0 = 1.0_rk

    exact(:,1,1) = phi0 * cos(pi * x/ Lx_twogroup)
    exact(:,2,1) = analytic_ratio_twogroup(xslib) * exact(:,1,1)
  endsubroutine analytic_fun_twogroup

  real(rk) pure function analytic_ratio_twogroup(xslib)
    use xs, only : XSLibrary
    use constant, only : pi
    type(XSLibrary), intent(in) :: xslib
    real(rk) :: rem2
    rem2 = xslib%mat(1)%sigma_t(2) - xslib%mat(1)%scatter(2,2,1)
    analytic_ratio_twogroup = &
      xslib%mat(1)%scatter(1,2,1) / (xslib%mat(1)%diffusion(2) * (pi/Lx_twogroup)**2 + rem2)
  endfunction analytic_ratio_twogroup

  real(rk) pure function analytic_keff_twogroup(xslib)
    use xs, only : XSLibrary
    use constant, only : pi
    type(XSLibrary), intent(in) :: xslib
    real(rk) :: rem1
    rem1 = xslib%mat(1)%sigma_t(1) - xslib%mat(1)%scatter(1,1,1)
    analytic_keff_twogroup = &
      (xslib%mat(1)%nusf(1) + xslib%mat(1)%nusf(2)*analytic_ratio_twogroup(xslib)) &
      / (xslib%mat(1)%diffusion(1) * (pi/Lx_twogroup)**2 + rem1)
  endfunction analytic_keff_twogroup

  subroutine analytic_fun_tworeg(x, xslib, exact)
    use xs, only : XSLibrary
    real(rk), intent(in) :: x(:)
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(out) :: exact(:,:,:)

    integer(ik) :: i

    real(rk), parameter :: phi0 = 1.0_rk
    real(rk) :: kr

    kr = &
      sqrt((xslib%mat(2)%sigma_t(1) - xslib%mat(2)%scatter(1,1,1)) &
      / xslib%mat(2)%diffusion(1))

    do i = 1,size(x)
      if (x(i) <= Lf_tworeg) then
        exact(i,1,1) = phi0 * cos(Bf_tworeg * x(i))
      else
        exact(i,1,1) = phi0 * cos(Bf_tworeg * Lf_tworeg) &
          * (cosh(kr*(x(i)-Lf_tworeg)) - sinh(kr*(x(i)-Lf_tworeg))/tanh(kr*(Lr_tworeg-Lf_tworeg)))
      endif
    enddo
  endsubroutine analytic_fun_tworeg

  real(rk) pure function analytic_keff_tworeg(xslib)
    use xs, only : XSLibrary
    type(XSLibrary), intent(in) :: xslib
    analytic_keff_tworeg = xslib%mat(1)%nusf(1) &
      / (xslib%mat(1)%diffusion(1) * Bf_tworeg**2 &
      + (xslib%mat(1)%sigma_t(1) - xslib%mat(1)%scatter(1,1,1)))
  endfunction analytic_keff_tworeg

  subroutine analytic_fun_p3(x, xslib, exact)
    use xs, only : XSLibrary
    use constant, only : pi
    real(rk), intent(in) :: x(:)
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(out) :: exact(:,:,:)

    real(rk) :: ratio

    real(rk), parameter :: phi0 = 1.0_rk

    ratio = analytic_ratio_p3(xslib)

    exact(:,1,1) = phi0 * cos(pi * x / Lx_p3) ! p0
    exact(:,1,2) = &
      (1.0_rk/3.0_rk) * phi0 * pi/Lx_p3 * (1.0_rk + 2.0_rk*ratio) * sin(pi * x / Lx_p3) &
      / xslib%mat(1)%sigma_t(1)
    exact(:,1,3) = ratio * exact(:,1,1) ! p2
    exact(:,1,4) = &
      (3.0_rk/7.0_rk) * phi0 * pi/Lx_p3 * ratio * sin(pi * x / Lx_p3) &
      / xslib%mat(1)%sigma_t(1)
  endsubroutine analytic_fun_p3

  real(rk) pure function analytic_ratio_p3(xslib)
    use xs, only : XSLibrary
    use constant, only : pi
    type(XSLibrary), intent(in) :: xslib
    real(rk), parameter :: bsq = (pi/Lx_p3)**2
    analytic_ratio_p3 = (-2.0_rk/15.0_rk/xslib%mat(1)%sigma_t(1)*bsq) &
      / (11.0_rk/21.0_rk/xslib%mat(1)%sigma_t(1)*bsq + xslib%mat(1)%sigma_t(1))
  endfunction analytic_ratio_p3

  real(rk) pure function analytic_keff_p3(xslib)
    use xs, only : XSLibrary
    use constant, only : pi
    type(XSLibrary), intent(in) :: xslib
    real(rk), parameter :: bsq = (pi/Lx_p3)**2
    real(rk) :: rem
    rem = xslib%mat(1)%sigma_t(1) - xslib%mat(1)%scatter(1,1,1)
    analytic_keff_p3 = xslib%mat(1)%nusf(1) &
      / (1.0_rk/3.0_rk/xslib%mat(1)%sigma_t(1)*bsq &
      + 2.0_rk/3.0_rk/xslib%mat(1)%sigma_t(1)*bsq * analytic_ratio_p3(xslib) &
      + rem)
  endfunction analytic_keff_p3

  subroutine analytic_fun_onegroup(x, phi_exact)
    use constant, only : pi
    real(rk), intent(in) :: x(:)
    real(rk), intent(out) :: phi_exact(:,:,:)
    real(rk), parameter :: phi0 = 1.0_rk
    phi_exact(:,1,1) = phi0 * cos(pi * x / Lx_twogroup)
  endsubroutine analytic_fun_onegroup

  real(rk) function analytic_keff_onegroup(xslib)
    use xs, only : XSLibrary
    use constant, only : pi
    type(XSLibrary), intent(in) :: xslib
    real(rk), parameter :: bsq = (pi / Lx_twogroup)**2
    real(rk) :: rem
    rem = xslib%mat(1)%sigma_t(1) - xslib%mat(1)%scatter(1,1,1)
    analytic_keff_onegroup = xslib%mat(1)%nusf(1) &
      / (xslib%mat(1)%diffusion(1) * bsq + rem)
  endfunction analytic_keff_onegroup

  subroutine analytic_pn_fission_matrix(xsmat, pnorder, fmat)
    use xs, only : XSMaterial
    type(XSMaterial), intent(in) :: xsmat
    integer(ik), intent(in) :: pnorder
    real(rk), intent(out) :: fmat(:,:)

    integer(ik) :: g, gprime

    fmat = 0.0_rk
    do g = 1,xsmat%ngroup
      do gprime = 1,xsmat%ngroup
        fmat(gprime, g) = xsmat%chi(gprime) * xsmat%nusf(g)
      enddo ! gprime = 1,xsmat%ngroup
    enddo ! g = 1,xsmat%ngroup
  endsubroutine analytic_pn_fission_matrix

  subroutine analytic_pn_transport_matrix(xsmat, bsq, pnorder, amat)
    use xs, only : XSMaterial
    use linalg, only : inv
    type(XSMaterial), intent(in) :: xsmat
    real(rk), intent(in) :: bsq
    integer(ik), intent(in) :: pnorder
    real(rk), intent(out) :: amat(:,:)

    real(rk), allocatable :: mat(:,:), invmat(:,:)
    real(rk), allocatable :: removal(:,:,:), inv_removal(:,:,:)

    real(rk), allocatable :: Dmat(:,:), Pmat(:,:), Nmat(:,:), Tmat(:,:)
    real(rk), allocatable :: Dhat(:,:), Phat(:,:), Nhat(:,:), That(:,:)

    integer(ik) :: g, n, idxn
    integer(ik) :: neq, rank
    integer(ik) :: lo, hi
    real(rk) :: xidxn

    allocate(removal(xsmat%ngroup,xsmat%ngroup,pnorder+1))
    allocate(inv_removal(xsmat%ngroup,xsmat%ngroup,pnorder+1))
    allocate(mat(xsmat%ngroup, xsmat%ngroup), invmat(xsmat%ngroup, xsmat%ngroup))
    do n = 1,pnorder+1
      idxn = n-1

      mat = 0.0_rk
      do g = 1,xsmat%ngroup
        mat(g,g) = xsmat%sigma_t(g)
      enddo ! g = 1,xsmat%ngroup

      if (idxn+1 <= xsmat%nmoment) then
        ! we have a scattering moment
        mat = mat - transpose(xsmat%scatter(:,:,idxn+1))
        call inv(xsmat%ngroup, mat, invmat)
        removal(:,:,idxn+1) = mat
        inv_removal(:,:,idxn+1) = invmat
      else
        invmat = mat
        do g = 1,xsmat%ngroup
          invmat(g,g) = 1.0_rk / invmat(g,g)
        enddo ! g = 1,xsmat%ngroup
        removal(:,:,idxn+1) = mat
        inv_removal(:,:,idxn+1) = invmat
      endif
    enddo ! n = 1,pnorder+1
    deallocate(mat, invmat)

    allocate(Dmat(xsmat%ngroup,xsmat%ngroup))
    allocate(Pmat(xsmat%ngroup,xsmat%ngroup))
    allocate(Nmat(xsmat%ngroup,xsmat%ngroup))
    allocate(Tmat(xsmat%ngroup,xsmat%ngroup))

    neq = int((pnorder+1)/2)
    rank = neq * xsmat%ngroup
    allocate(Dhat(rank, rank))
    allocate(Phat(rank, rank))
    allocate(Nhat(rank, rank))
    allocate(That(rank, rank))

    do n = 1,neq
      idxn = 2*(n-1)
      xidxn = real(idxn, rk)

      lo = xsmat%ngroup*(n-1)+1
      hi = xsmat%ngroup*n

      Tmat = removal(:,:,idxn+1)
      That(lo:hi,lo:hi) = Tmat

      Dmat = (xidxn+1.0_rk)*(xidxn+1.0_rk)/((2.0_rk*xidxn+1.0_rk)*(2.0_rk*xidxn+3.0_rk))*inv_removal(:,:,idxn+1+1)
      if ((idxn - 1) > 0) then
        Dmat = Dmat + xidxn*xidxn/((2.0_rk*xidxn+1.0_rk)*(2.0_rk*xidxn-1.0_rk))*inv_removal(:,:,idxn+1-1)
      endif
      Dhat(lo:hi,lo:hi) = Dmat

      if ((idxn - 1) > 0) then
        Pmat = xidxn*(xidxn-1.0_rk)/((2.0_rk*xidxn+1.0_rk)*(2.0_rk*xidxn-1.0_rk))*inv_removal(:,:,idxn+1-1)
        Phat(lo:hi,xsmat%ngroup*(n-2)+1:xsmat%ngroup*(n-1)) = Pmat ! TODO transpose block?
      endif

      if ((idxn + 1) < pnorder-1) then
        Nmat = (xidxn+1.0_rk)*(xidxn+2.0_rk)/((2.0_rk*xidxn+1.0_rk)*(2.0_rk*xidxn+3.0_rk))*inv_removal(:,:,idxn+1+1)
        Nhat(lo:hi,xsmat%ngroup*n+1:xsmat%ngroup*(n+1)) = Nmat ! TODO transpose block?
      endif

    enddo ! n = 1,neq

    amat = (Dhat + Phat + Nhat) * bsq + That

    deallocate(removal, inv_removal)
    deallocate(Dmat, Pmat, Nmat, Tmat)
    deallocate(Dhat, Phat, Nhat, That)
  endsubroutine analytic_pn_transport_matrix

  subroutine analytic_fun_pn(x, xslib, exact, keff)
    use xs, only : XSLibrary
    use constant, only : pi
    use linalg, only : geneig, inv
    real(rk), intent(in) :: x(:)
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(out) :: exact(:,:,:)
    real(rk), intent(out), optional :: keff

    integer(ik) :: pnorder, neq, rank
    integer(ik) :: i, g, n, idxn
    real(rk) :: bsq, xidxn

    real(rk), allocatable :: a(:,:), f(:,:)
    complex(rk), allocatable :: eigval(:)
    real(rk), allocatable :: eigvec(:,:)

    integer(ik) :: idx
    real(rk) :: kcalc

    real(rk), allocatable :: amplitude(:,:)
    real(rk), allocatable :: vec(:)
    real(rk), allocatable :: removal(:,:), inv_removal(:,:)

    ! a bit iffy since this is taken with intent(out)
    pnorder = size(exact,3)

    neq = int((pnorder+1)/2)
    rank = neq * xslib%ngroup
    allocate(f(rank,rank))
    allocate(a(rank,rank))

    bsq = (pi / Lx_p3)**2 ! NOTE: hard wired

    call analytic_pn_fission_matrix(xslib%mat(1), pnorder, f)
    call analytic_pn_transport_matrix(xslib%mat(1), bsq, pnorder, a)

    allocate(eigval(rank))
    allocate(eigvec(rank,rank))
    call geneig(rank, a, f, eigval, eigvec)

    kcalc = huge(1.0_rk)
    idx = -1
    do i = 1,rank
      if ((abs(eigval(i)) < kcalc) .and. (abs(eigval(i)) > epsilon(1.0_rk)) .and. (abs(aimag(eigval(i))) < epsilon(1.0_rk))) then
        idx = i
        kcalc = abs(eigval(i))
        write(*,*) 1.0_rk/kcalc
      endif
    enddo ! i = 1,rank
    kcalc = 1.0_rk/kcalc

    if (present(keff)) then
      keff = kcalc
    endif

    allocate(amplitude(xslib%ngroup,pnorder+1))

    ! first pass: copy even moments
    do n = 1,neq
      idxn = 2*(n-1)
      amplitude(:,idxn+1) = eigvec(xslib%ngroup*(n-1)+1:xslib%ngroup*n,idx)
    enddo ! n = 1,neq

    ! second pass: compute odd moments
    allocate(vec(xslib%ngroup))
    allocate(removal(xslib%ngroup,xslib%ngroup))
    allocate(inv_removal(xslib%ngroup,xslib%ngroup))
    do n = 1,neq
      idxn = 2*(n-1)+1

      if (idxn+1 <= xslib%nmoment) then
        ! we have a scattering moment
        removal = 0.0_rk
        do g = 1,xslib%ngroup
          removal(g,g) = xslib%mat(1)%sigma_t(g)
        enddo ! g = 1,xslib%ngroup
        removal = removal - transpose(xslib%mat(1)%scatter(:,:,idxn+1))
        call inv(xslib%ngroup, removal, inv_removal)
      else
        removal = 0.0_rk
        inv_removal = 0.0_rk
        do g = 1,xslib%ngroup
          removal(g,g) = xslib%mat(1)%sigma_t(g)
          inv_removal(g,g) = 1.0_rk / xslib%mat(1)%sigma_t(g)
        enddo ! g = 1,xslib%ngroup
      endif

      xidxn = real(idxn, rk)
      vec = xidxn / (2.0_rk * xidxn + 1.0_rk) * amplitude(:,idxn+1-1)
      if ((idxn + 1) <= pnorder-1) then
        vec = vec + (xidxn + 1.0_rk) / (2.0_rk * xidxn + 1.0_rk) * amplitude(:,idxn+1+1)
      endif

      amplitude(:,idxn+1) = sqrt(bsq) * matmul(inv_removal, vec)
    enddo ! n = 1,neq
    deallocate(vec)
    deallocate(removal, inv_removal)

    amplitude = amplitude / amplitude(1,1)

    do n = 1,pnorder
      idxn = n-1
      if (mod(idxn, 2) == 0) then
        do g = 1,xslib%ngroup
          exact(:,g,idxn+1) = amplitude(g,idxn+1) * cos(sqrt(bsq) * x)
        enddo ! g = 1,xslib%ngroup
      else
        do g = 1,xslib%ngroup
          exact(:,g,idxn+1) = amplitude(g,idxn+1) * sin(sqrt(bsq) * x)
        enddo ! g = 1,xslib%ngroup
      endif
    enddo ! n = 1,pnorder

    deallocate(eigval, eigvec)
    deallocate(amplitude)
    deallocate(a, f)
  endsubroutine analytic_fun_pn

endmodule analytic
