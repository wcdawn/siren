module geometry
use kind
implicit none

private

public :: uniform_refinement

contains

  subroutine uniform_refinement(nx, hx, mat_map)
    integer(ik), intent(inout) :: nx
    real(rk), intent(inout) :: hx
    integer(ik), allocatable, intent(inout) :: mat_map(:)

    integer(ik) :: i

    integer(ik) :: nx_old
    real(rk) :: hx_old
    integer(ik), allocatable :: mat_map_old(:)

    nx_old = nx
    hx_old = hx
    allocate(mat_map_old(nx_old))
    mat_map_old = mat_map

    nx = 2 * nx_old
    hx = 0.5d0 * hx_old

    deallocate(mat_map)
    allocate(mat_map(nx))

    do i = 1,nx_old
      mat_map(2*(i-1)+1) = mat_map_old(i)
      mat_map(2*(i-1)+2) = mat_map_old(i)
    enddo ! i = 1,nx_old

    deallocate(mat_map_old)
  endsubroutine uniform_refinement

endmodule geometry
