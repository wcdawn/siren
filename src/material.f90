module material
use kind, only : rk, ik
implicit none

private

public :: material_idx_from_name

contains

  subroutine material_idx_from_name(nx, xs, name, idx)
    use xs, only : XSLibrary
    integer, intent(in) :: nx
    type(XSLibrary), intent(in) :: xs
    character(*), intent(in) :: name(:)
    integer(ik), intent(out) :: idx(:)

    integer(ik) :: i, id

    do i = 1,nx
      read(name(i), *) id
      idx(i) = id
    enddo ! i = 1,nx
  endsubroutine material_idx_from_name

endmodule material
