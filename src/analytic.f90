module analytic
use kind
implicit none

private
public :: analytic_error

real(rk), parameter :: Lx_twogroup = 1e2_rk

contains

  subroutine analytic_error(analytic_name, nx, ngroup, pnorder, xslib, dx, phi, keff)
    use xs, only : XSLibrary
    use linalg, only : norm
    use output, only : output_write
    character(*), intent(in) :: analytic_name
    integer(ik), intent(in) :: nx, ngroup, pnorder
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(in) :: dx(:)
    real(rk), intent(in) :: phi(:,:,:) ! (nx,ngroup,neven)
    real(rk), intent(in) :: keff

    ! assume that the coordinate system starts at xleft==0.0
    ! we don't really know in general since all we have are deltas
    ! but we need an exact coordinate for the error analysis
    real(rk), parameter :: x0 = 0.0_rk

    real(rk) :: xleft
    integer(ik) :: neven
    integer(ik) :: i

    real(rk), allocatable :: x(:)
    real(rk), allocatable :: phi_analysis(:,:,:), phi_exact(:,:,:), phi_diff(:,:,:)
    real(rk) :: phi_zero
    real(rk) :: keff_exact, keff_diff
    real(rk) :: ratio_exact, ratio_siren, ratio_diff
    real(rk) :: linferr

    character(1024) :: line

    allocate(x(nx))
    xleft = x0
    do i = 1,nx
      x(i) = xleft + 0.5d0 * dx(i)
      xleft = xleft + dx(i)
    enddo

    ! work on a copy so that we can change normalization
    neven = max((pnorder + 1) / 2, 1)
    allocate(phi_analysis(nx,ngroup,neven))
    phi_analysis = phi
    ! extrapolate to estimate phi(0) for the first group, scalar flux
    ! this is used to normalize such that phi(0) = 1.0 for first group, scalar flux
    ! this is usually equivalent to phi0 == 1 from the analytic solutions
    phi_zero = -(phi(2,1,1)-phi(1,1,1)) * dx(1) / (dx(1) + dx(2)) + phi(1,1,1)
    phi_analysis = phi_analysis / phi_zero
    
    ! get phi_exact
    allocate(phi_exact(nx,ngroup,neven))
    select case (analytic_name)
      case ('twogroup')
        call analytic_fun_twogroup(x, xslib, phi_exact)
      case default
        ! TODO exception handling
        write(*,*) 'unknown analytic_name: ' // trim(adjustl(analytic_name))
        stop
    endselect

    allocate(phi_diff(nx,ngroup,neven))
    phi_diff = phi_exact - phi_analysis
    ! NOTE: using the two-norm would be reasonable.
    ! However, one would need to be careful since the mesh is not uniform.
    linferr = norm(-1, phi_diff(:,1,1))

    ! get keff_exact
    select case (analytic_name)
      case ('twogroup')
        keff_exact = analytic_keff_twogroup(xslib)
      case default
        write(*,*) 'second -- unknwon analytic_name :' // trim(adjustl(analytic_name))
        stop
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
    endif
    write(line, '(a,es13.6)') 'linferr = ', linferr
    call output_write(line)
    call output_write('')

    deallocate(phi_analysis, phi_exact, phi_diff)
    deallocate(x)
  endsubroutine analytic_error

  subroutine analytic_fun_twogroup(x, xslib, exact)
    use xs, only : XSLibrary
    use constant, only : pi
    real(rk), intent(in) :: x(:)
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(out) :: exact(:,:,:)

    real(rk) :: ratio

    real(rk), parameter :: phi0 = 1.0_rk

    exact(:,1,1) = phi0*cos(pi*x/Lx_twogroup)
    ratio = analytic_ratio_twogroup(xslib)
    exact(:,2,1) = ratio * exact(:,1,1)
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

endmodule analytic
