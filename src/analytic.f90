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
      case ('twogroup')
        call analytic_fun_twogroup(x, xslib, phi_exact)
      case ('tworeg')
        call analytic_fun_tworeg(x, xslib, phi_exact)
      case ('analytic_p3')
        call analytic_fun_p3(x, xslib, phi_exact)
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
      case ('twogroup')
        keff_exact = analytic_keff_twogroup(xslib)
      case ('tworeg')
        keff_exact = analytic_keff_tworeg(xslib)
      case ('analytic_p3')
        keff_exact = analytic_keff_p3(xslib)
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

endmodule analytic
