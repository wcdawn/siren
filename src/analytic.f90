module analytic
use kind
implicit none

private
public :: analytic_error

contains

  subroutine analytic_error(analytic_name, nx, ngroup, pnorder, xslib, dx, phi, keff)
    use xs, only : XSLibrary
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
    real(rk), allocatable :: phi_analysis(:,:,:), phi_exact(:,:,:)

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
    
    allocate(phi_exact(nx,ngroup,neven))
    select case (analytic_name)
      case ('twogroup')
        call analytic_fun_twogroup(x, xslib, phi_exact)
      case default
        ! TODO exception handling
        write(*,*) 'unknown analytic_name: ' // trim(adjustl(analytic_name))
        stop
    endselect

    deallocate(phi_analysis, phi_exact)
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
    real(rk), parameter :: Lx = 100.0_rk

    exact(:,1,1) = phi0*cos(pi*x/Lx)
    ratio = analytic_ratio_twogroup(xslib, Lx)
    exact(:,2,1) = ratio * exact(:,1,1)
  endsubroutine analytic_fun_twogroup

  real(rk) function analytic_ratio_twogroup(xslib, Lx)
    use xs, only : XSLibrary
    use constant, only : pi
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(in) :: Lx
    real(rk) :: rem2
    rem2 = xslib%mat(1)%sigma_t(2) - xslib%mat(1)%scatter(2,2,1)
    analytic_ratio_twogroup = &
      xslib%mat(1)%scatter(2,1,1) / (xslib%mat(1)%diffusion(2) * (pi/Lx)**2 + rem2)
    ! TODO
    write(*,'(es13.6)') analytic_ratio_twogroup
  endfunction analytic_ratio_twogroup

endmodule analytic
