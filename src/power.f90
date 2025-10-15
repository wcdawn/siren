module power
use kind, only : rk, ik
implicit none (external)

private

public :: power_calculate

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
        ! NOTE: using nusf for now (fission neutron production rate)
        ! could alternatively use sigma_f
        power(i) = sum(xslib%mat(mthis)%nusf(:) * flux(i,:))
      else
        power(i) = 0d0
      endif
    enddo
  endsubroutine power_calculate

endmodule power
