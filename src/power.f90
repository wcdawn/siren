module power
use kind, only : rk, ik
implicit none

private

public :: power_calculate, power_total

contains

  subroutine power_calculate(nx, mat_map, xslib, flux, power)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    integer(ik), intent(in) :: mat_map(:)
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(in) :: flux(:,:) ! (nx,ngroup)
    real(rk), intent(out) :: power(:) ! (nx)

    integer(ik) :: i
    integer(ik) :: mthis

    do i = 1,nx
      mthis = mat_map(i)
      if (xslib%mat(mthis)%is_fiss) then
        ! use sigma_f if available, otherwise use nusf (it must be available)
        ! using simga_f is the "fission reaction rate"
        ! using nusf is the "neutron production rate due to fission"
        if (allocated(xslib%mat(mthis)%sigma_f)) then
          power(i) = sum(xslib%mat(mthis)%sigma_f(:) * flux(i,:))
        else
          power(i) = sum(xslib%mat(mthis)%nusf(:) * flux(i,:))
        endif
      else
        power(i) = 0d0
      endif
    enddo
  endsubroutine power_calculate

  real(rk) pure function power_total(dx, power)
    real(rk), intent(in) :: dx (:)
    real(rk), intent(in) :: power(:)
    power_total = sum(dx*power)
  endfunction power_total


endmodule power
