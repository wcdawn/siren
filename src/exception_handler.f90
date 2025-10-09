module exception_handler
use kind

public :: exception_note, exception_warning, exception_fatal, exception_summary

private

integer(ik), parameter :: exception_lvl_note = 1
integer(ik), parameter :: exception_lvl_warning = 2
integer(ik), parameter :: exception_lvl_fatal = 3

type Exception_obj
  integer :: lvl
  character(1024) :: msg
endtype Exception_obj

type(Exception_obj), allocatable :: exception_stack(:)

contains

  subroutine exception_append(stack, x)
    type(Exception_obj), allocatable, intent(inout) :: stack(:)
    type(Exception_obj), intent(in) :: x
    integer(ik) :: length
    type(Exception_obj), allocatable :: stack_old(:)

    if (allocated(stack)) then
      length = size(stack)
      allocate(stack_old(length))
      stack_old = stack
      deallocate(stack)
    else
      length = 0
    endif

    allocate(stack(length+1))
    if (allocated(stack_old)) then
      stack(1:length) = stack_old
    endif
    stack(length+1) = x

    if (allocated(stack_old)) then
      deallocate(stack_old)
    endif
  endsubroutine exception_append

  subroutine exception_note(msg)
    use output, only : output_write
    character(*), intent(in) :: msg
    type(Exception_obj) :: obj
    obj%lvl = exception_lvl_note
    obj%msg = trim(adjustl(msg))
    call exception_append(exception_stack, obj)
    call output_write('NOTE :: ' // trim(adjustl(msg)))
  endsubroutine exception_note

  subroutine exception_warning(msg)
    use output, only : output_write
    character(*), intent(in) :: msg
    type(Exception_obj) :: obj
    obj%lvl = exception_lvl_warning
    obj%msg = trim(adjustl(msg))
    call exception_append(exception_stack, obj)
    call output_write('WARNING :: ' // trim(adjustl(msg)))
  endsubroutine exception_warning

  subroutine exception_fatal(msg)
    use output, only : output_write
    character(*), intent(in) :: msg
    type(Exception_obj) :: obj
    obj%lvl = exception_lvl_fatal
    obj%msg = trim(adjustl(msg))
    call exception_append(exception_stack, obj)
    call output_write('FATAL :: ' // trim(adjustl(msg)))
    call exception_summary()
    stop
  endsubroutine exception_fatal

  subroutine exception_summary()
    use output, only : output_write
    integer(ik) :: i
    integer(ik) :: count_note, count_warning, count_fatal
    character(1024) :: line

    call output_write('=== EXCEPTION SUMMARY ===')

    count_note = 0
    count_warning = 0
    count_fatal = 0

    if (allocated(exception_stack)) then
      ! right now, just printed in order that they occured
      do i = 1,size(exception_stack)
        select case (exception_stack(i)%lvl)
          case (exception_lvl_note)
            call output_write('NOTE :: ' // trim(adjustl(exception_stack(i)%msg)))
            count_note = count_note + 1
          case (exception_lvl_warning)
            call output_write('WARNING :: ' // trim(adjustl(exception_stack(i)%msg)))
            count_warning = count_warning + 1
          case (exception_lvl_fatal)
            call output_write('FATAL :: ' // trim(adjustl(exception_stack(i)%msg)))
            count_fatal = count_fatal + 1
        endselect
      enddo
    endif

    call output_write('Encountered:')
    write(line, '(2x,i0,a)') count_note, ' note(s)'
    call output_write(line)
    write(line, '(2x,i0,a)') count_warning, ' warning(s)'
    call output_write(line)
    write(line, '(2x,i0,a)') count_fatal, ' fatal error(s)'
    call output_write(line)
    call output_write('')

    if (allocated(exception_stack)) then
      deallocate(exception_stack)
    endif

  endsubroutine exception_summary

endmodule exception_handler
