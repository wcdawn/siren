module timer
use kind

public :: timer_init, timer_start, timer_stop, timer_summary

private

type Timer_obj
  character(32) :: name
  logical :: running
  real(rk) :: elapsed
  integer :: count
endtype Timer_obj

type(Timer_obj), allocatable :: timer_arr(:)
integer :: time_rate

contains

  subroutine timer_init()
    call system_clock(count_rate=time_rate)
    call timer_start('total')
  endsubroutine timer_init

  subroutine timer_start(name)
    character(*), intent(in) :: name
    integer :: i
    integer :: length
    type(Timer_obj), allocatable :: timer_arr_cpy(:)

    if (allocated(timer_arr)) then
      do i = 1,size(timer_arr)
        if (timer_arr(i)%name == name) then
          ! found the timer
          if (.not. timer_arr(i)%running) then
            ! if the timer is not running, start it
            ! if the timer was already running, do nothing
            call system_clock(count=timer_arr(i)%count)
            timer_arr(i)%running = .true.
            return
          endif
        endif
      enddo
    endif

    if (allocated(timer_arr)) then
      length = size(timer_arr)
      allocate(timer_arr_cpy(length))
      timer_arr_cpy = timer_arr
      deallocate(timer_arr)
    else
      length = 0
    endif

    allocate(timer_arr(length+1))
    if (allocated(timer_arr_cpy)) then
      timer_arr(1:length) = timer_arr_cpy
    endif

    timer_arr(length+1)%name = name
    call system_clock(count=timer_arr(length+1)%count)
    timer_arr(length+1)%running = .true.
    timer_arr(length+1)%elapsed = 0.0_rk

    if (allocated(timer_arr_cpy)) then
      deallocate(timer_arr_cpy)
    endif
  endsubroutine timer_start

  subroutine timer_stop_this(t)
    type(Timer_obj), intent(inout) :: t
    integer :: final_count
    if (t%running) then
      call system_clock(count=final_count)
      t%elapsed = t%elapsed + (final_count - t%count) / real(time_rate, rk)
      t%running = .false.
    endif
  endsubroutine timer_stop_this

  subroutine timer_stop(name)
    character(*), intent(in) :: name
    integer :: i
    do i = 1,size(timer_arr)
      if (timer_arr(i)%name == name) then
        if (timer_arr(i)%running) then
          call timer_stop_this(timer_arr(i))
        endif
        return
      endif
    enddo
  endsubroutine timer_stop

  subroutine timer_summary()
    use output, only : output_write
    integer :: i, max_len
    real(rk) :: total
    character(1024) :: line, format_str

    call output_write('=== TIMER SUMMARY ===')

    ! stop the total timer and store the total elapsed time to compute fractions/percentages
    max_len = 0
    total = 0.0_rk
    do i = 1,size(timer_arr)
      call timer_stop_this(timer_arr(i))
      if (timer_arr(i)%name == 'total') then
        total = timer_arr(i)%elapsed
      endif
      max_len = max(max_len, len(trim(adjustl(timer_arr(i)%name))))
    enddo

    format_str = ''
    write(format_str, '(a,i0,a)') '(a,', max_len+1-4, 'x'
    format_str = trim(adjustl(format_str)) // ',a)'
    write(line, format_str) 'name', 'dt [s]    [%]'
    call output_write(line)

    format_str = ''
    write(format_str, '(a,i0)') '(a', max_len+1
    if (total < 1e1_rk) then
      format_str = trim(adjustl(format_str)) // ', f6.4'
    elseif (total < 1e2_rk) then
      format_str = trim(adjustl(format_str)) // ', f6.3'
    elseif (total < 1e3_rk) then
      format_str = trim(adjustl(format_str)) // ', f6.2'
    else
      format_str = trim(adjustl(format_str)) // ', es10.3'
    endif
    format_str = trim(adjustl(format_str)) // ', 1x, f6.2)'

    do i = 1,size(timer_arr)
      write(line, format_str) timer_arr(i)%name, timer_arr(i)%elapsed, timer_arr(i)%elapsed/total*1e2_rk
      call output_write(line)
    enddo

    call output_write('')

    deallocate(timer_arr)
  endsubroutine timer_summary

endmodule timer
