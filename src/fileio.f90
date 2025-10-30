module fileio
implicit none

private

public :: fileio_open_read, fileio_open_write

contains

  subroutine fileio_open_read(fname, iounit)
    character(*), intent(in) :: fname
    integer, intent(in) :: iounit
    integer :: ios
    open(unit=iounit, file=trim(adjustl(fname)), action='read', status='old', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'failed to open file in fileio_open_read'
      write(*,*) 'fname=', trim(adjustl(fname)), ' iounit=', iounit
      stop
    endif
  endsubroutine fileio_open_read

  subroutine fileio_open_write(fname, iounit, replace_in)
    character(*), intent(in) :: fname
    integer, intent(in) :: iounit
    logical, intent(in), optional :: replace_in
    integer :: ios
    logical :: replace
    if (present(replace_in)) then
      replace = replace_in
    else
      replace = .true. ! default replace
    endif
    if (replace) then
      open(unit=iounit, file=fname, action='write', status='replace', iostat=ios)
    else
      open(unit=iounit, file=fname, action='write', status='new', iostat=ios)
    endif
    if (ios /= 0) then
      write(*,*) 'failed to open file in fileio_open_read'
      write(*,*) 'fname=', trim(adjustl(fname)), 'iounit=', iounit
      stop
    endif
  endsubroutine fileio_open_write

endmodule fileio
