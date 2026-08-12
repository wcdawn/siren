module material
use kind, only : rk, ik
implicit none

private

public :: material_idx_from_name

contains

  subroutine material_idx_from_name(nx, xs, name, idx)
    use xs, only : XSLibrary
    use exception_handler, only : exception_fatal
    integer, intent(in) :: nx
    type(XSLibrary), intent(in) :: xs
    character(*), intent(in) :: name(:)
    integer(ik), intent(out) :: idx(:)

    integer(ik) :: i, j

    idx = -1
    do i = 1,nx

      ! linear search
      do j = 1,xs%niso
        if (xs%mat(j)%name == name(i)) then
          idx(i) = j
          exit
        endif
      enddo ! j = 1,xs%niso

      if (idx(i) < 0) then
        call exception_fatal(&
          'Failed to identify material in mat_map: ' &
          // trim(adjustl(name(i))))
      endif

    enddo ! i = 1,nx
  endsubroutine material_idx_from_name

endmodule material
