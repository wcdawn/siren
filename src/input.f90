module input
use kind
implicit none

private

! required input
integer(ik) :: nx
real(rk) :: hx
character(1024) :: xslib_fname
integer(ik), allocatable :: mat_map(:)
integer(ik) :: pnorder

! optional input (with defaults)
real(rk) :: k_tol = 1d-6, phi_tol = 1d-5
integer(ik) :: max_iter = 100
integer(ik) :: refine = 0

public :: input_read, input_cleanup

public :: nx, hx, xslib_fname, mat_map
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

    integer(ik) :: i

    call fileio_open_read(trim(adjustl(fname)), iounit)

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
          mat_map = 0
        case ('hx')
          read(line, *) card, hx
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
        case default
          write(*,*) 'card:', trim(adjustl(card))
          stop 'unknown input card'
      endselect
    enddo

    write(*,*) '=== INPUT ==='
    write(*, '(a,i0)') 'nx = ', nx
    write(*, '(a,es13.6)') 'hx = ', hx
    write(*, '(a,a,a)') 'xslib_fname = "', trim(adjustl(xslib_fname)), '"'
    write(*, '(a)', advance='no') 'mat_map ='
    do i = 1,nx
      write(*, '(a,i0)', advance='no') ' ', mat_map(i)
    enddo
    write(*,*)
    write(*,'(a,i0)') 'refine = ', refine
    write(*,'(a,es13.6)') 'k_tol = ', k_tol
    write(*,'(a,es13.6)') 'phi_tol = ', phi_tol
    write(*,*)

    close(iounit)
  endsubroutine input_read

  subroutine input_cleanup()
    if (allocated(mat_map)) then
      deallocate(mat_map)
    endif
  endsubroutine input_cleanup

endmodule input
