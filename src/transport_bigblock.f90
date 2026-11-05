module transport_bigblock
use kind, only : rk, ik
implicit none

private

public :: transport_bigblock_power_iteration

contains

  subroutine transport_bigblock_power_iteration( &
      nx, dx, mat_map, xs, boundary_left, boundary_right, albedo_coeff, &
      calc_type, k_tol, phi_tol, max_iter, pnorder, keff, sigma_tr, phi)
    use exception_handler, only : exception_fatal
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xs
    character(*), intent(in) :: boundary_left, boundary_right
    real(rk), intent(in) :: albedo_coeff ! A \in [-1,1]
    character(*), intent(in) :: calc_type
    real(rk), intent(in) :: k_tol, phi_tol
    integer(ik), intent(in) :: max_iter
    integer(ik), intent(in) :: pnorder
    real(rk), intent(out) :: keff
    real(rk), intent(out) :: sigma_tr(:,:,:) ! (nx, ngroup, nmoment)
    real(rk), intent(out) :: phi(:,:,:) ! (nx, ngroup, nmoment)


    if (boundary_left /= 'mirror') then
      call exception_fatal('In bigblock, boundary_left must be mirror.')
    elseif (boundary_right /= 'mirror') then
      call exception_fatal('In bigblock, boundary_right must be mirror.')
    elseif (calc_type /= 'eigenvalue') then
      call exception_fatal('In bigblock, calc_type must be eigenvalue.')
    endif

  endsubroutine transport_bigblock_power_iteration

endmodule transport_bigblock
