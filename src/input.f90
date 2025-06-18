module input
use kind
implicit none

private

! required input
integer(ik) :: nx
real(rk), allocatable :: dx(:)
character(1024) :: xslib_fname
integer(ik), allocatable :: mat_map(:)
integer(ik) :: pnorder

! optional input (with defaults)
real(rk) :: k_tol = 1d-6, phi_tol = 1d-5
integer(ik) :: max_iter = 100
integer(ik) :: refine = 0

! optional
character(16) :: boundary_left = 'mirror', boundary_right = 'mirror'

public :: input_read, input_cleanup

public :: nx, dx, xslib_fname, mat_map
public :: k_tol, phi_tol, max_iter
public :: pnorder
public :: refine

contains

  subroutine input_read(fname)
    use fileio, only : fileio_open_read
    character(*), intent(in) :: fname

    integer, parameter :: iounit = 15
    character(*), parameter :: comment_char = '#'

    character(1024) :: line, card
    integer :: ios

    call fileio_open_read(trim(adjustl(fname)), iounit)
    ! TODO input echo

    do
      read(iounit, '(a)', iostat=ios) line
      if (ios /= 0) then
        exit
      endif

      if (line == '') then
        cycle
      endif

      line = trim(adjustl(line))
      if (line(1:1) == comment_char) then
        cycle
      endif

      read(line, *) card

      select case (card)
        case ('nx')
          read(line, *) card, nx
          allocate(mat_map(nx))
          allocate(dx(nx))
          mat_map = 0
        case ('dx')
          read(line, *) card, dx
        case ('xslib_fname')
          xslib_fname = ''
          read(line, *) card, xslib_fname
        case ('mat_map')
          read(line, *) card, mat_map
        case ('k_tol')
          read(line, *) card, k_tol
        case ('phi_tol')
          read(line, *) card, phi_tol
        case ('max_iter')
          read(line, *) card, max_iter
        case ('refine')
          read(line, *) card, refine
        case ('pnorder')
          read(line, *) card, pnorder
        case ('boundary_left')
          read(line, *) card, boundary_left
        case ('boundary_right')
          read(line, *) card, boundary_right
        case default
          write(*,*) 'card:', trim(adjustl(card))
          stop 'unknown input card'
      endselect
    enddo

    call input_summary()

    call input_check()

    close(iounit)
  endsubroutine input_read

  subroutine input_check()
    if (trim(adjustl(boundary_left)) /= 'mirror') then
      stop 'boundary_left must be set to mirror (for now)'
    elseif ((pnorder /= 0) .and. (trim(adjustl(boundary_right)) /= 'mirror')) then
      stop 'all boundaries must be set to mirror for PN (non-diffusion) calculation'
    endif
  endsubroutine input_check

  subroutine input_summary()
    use output, only : output_write
    character(1024) :: line

    call output_write('=== INPUT ===')
    call output_write('xslib_fname = "' // trim(adjustl(xslib_fname)) // '"')
    write(line, '(a,i0)') 'refine = ', refine
    call output_write(line)
    write(line, '(a,es13.6)') 'k_tol = ', k_tol
    call output_write(line)
    write(line, '(a,es13.6)') 'phi_tol = ', phi_tol
    call output_write(line)
    call output_write('boundary_left = *' // trim(adjustl(boundary_left)) // '*')
    call output_write('boundary_right = *' // trim(adjustl(boundary_right)) // '*')

    call output_write('')
  endsubroutine input_summary

  subroutine input_cleanup()
    if (allocated(mat_map)) then
      deallocate(mat_map)
    endif
    if (allocated(dx)) then
      deallocate(dx)
    endif
  endsubroutine input_cleanup

endmodule input
