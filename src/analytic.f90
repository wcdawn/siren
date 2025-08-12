module analytic
use kind
implicit none

private
public :: analytic_error

contains

  subroutine analytic_error(analytic_name, nx, ngroup, nmoment, dx, phi)
    character(*), intent(in) :: analytic_name
    integer(ik), intent(in) :: nx, ngroup, nmoment
    real(rk), intent(in) :: dx(:)
    real(rk), intent(in) :: phi(:,:,:) ! (nx,ngroup,nmoment)

    ! assume that the coordinate system starts at xleft==0.0
    ! we don't really know in general since all we have are deltas
    ! but we need an exact coordinate for the error analysis
    real(rk), parameter :: x0 = 0.0_rk

    real(rk) :: xleft
    integer(ik) :: i

    real(rk), allocatable :: x(:)

    allocate(x(nx))

    xleft = x0
    do i = 1,nx
      x(i) = xleft + 0.5d0 * dx(i)
      xleft = xleft + dx(i)
    enddo

    select case (analytic_name)
      case default
        write(*,*) 'unknown analytic_name: ' // trim(adjustl(analytic_name))
        stop
    endselect

    deallocate(x)
  endsubroutine analytic_error

endmodule analytic
