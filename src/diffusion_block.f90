module diffusion_block
use kind
implicit none

private
public :: diffusion_block_power_iteration

contains

  subroutine diffusion_block_power_iteration(nx, dx, mat_map, xslib, boundary_right, k_tol, phi_tol, max_iter, keff, flux)
    use xs, only : XSLibrary
    use linalg, only : trid_block
    use output, only : output_write
    use timer, only : timer_start, timer_stop
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    character(*), intent(in) :: boundary_right
    real(rk), intent(in) :: k_tol, phi_tol
    integer(ik), intent(in) :: max_iter
    real(rk), intent(out) :: keff
    real(rk), intent(out) :: flux(:,:) ! (nx,ngroup)
  endsubroutine diffusion_block_power_iteration

endmodule diffusion_block
