module geometry
use kind, only : rk, ik
implicit none

private

public :: geometry_uniform_refinement, geometry_summary

contains

  subroutine geometry_uniform_refinement(nx, dx, mat_map)
    integer(ik), intent(inout) :: nx
    real(rk), allocatable, intent(inout) :: dx(:)
    integer(ik), allocatable, intent(inout) :: mat_map(:)

    integer(ik) :: i

    integer(ik) :: nx_old
    real(rk), allocatable :: dx_old(:)
    integer(ik), allocatable :: mat_map_old(:)

    nx_old = nx

    allocate(mat_map_old(nx_old))
    mat_map_old = mat_map

    allocate(dx_old(nx_old))
    dx_old = dx

    nx = 2 * nx_old

    deallocate(mat_map)
    allocate(mat_map(nx))

    deallocate(dx)
    allocate(dx(nx))

    do i = 1,nx_old
      mat_map(2*(i-1)+1) = mat_map_old(i)
      mat_map(2*(i-1)+2) = mat_map_old(i)
      dx(2*(i-1)+1) = 0.5d0*dx_old(i)
      dx(2*(i-1)+2) = 0.5d0*dx_old(i)
    enddo ! i = 1,nx_old

    deallocate(mat_map_old)
    deallocate(dx_old)
  endsubroutine geometry_uniform_refinement

  subroutine geometry_summary(nx, dx, mat_map)
    use output, only : output_write
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:)
    integer(ik), intent(in) :: mat_map(:)

    character(2**16) :: line, tmp
    integer(ik) :: i

    integer(ik), parameter :: max_nx_listing = 100

    call output_write('=== GEOMETRY ===')

    write(line, '(a,i0)') 'nx = ', nx
    call output_write(line)

    if (nx <= max_nx_listing) then
      write(line, '(a)') 'dx ='
      do i = 1,nx
        write(tmp, '(es9.2)') dx(i)
        line = trim(adjustl(line)) // ' ' // trim(adjustl(tmp))
      enddo ! i = 1,nx
      call output_write(line)
    else
      call output_write('dx list omitted for space')
    endif

    write(line, '(a,es13.6)') 'max_dx = ', maxval(dx)
    call output_write(line)
    write(line, '(a,es13.6)') 'min_dx = ', minval(dx)
    call output_write(line)

    if (nx <= max_nx_listing) then
      write(line, '(a)') 'mat_map ='
      do i = 1,nx
        write(tmp, '(i0)') mat_map(i)
        line = trim(adjustl(line)) // ' ' // trim(adjustl(tmp))
      enddo ! i = 1,nx
      call output_write(line)
    else
      call output_write('mat_map list omitted for space')
    endif

    call output_write('')
  endsubroutine geometry_summary

endmodule geometry
