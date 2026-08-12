module transient
use kind, only : rk, ik
implicit none

private

type DelayedNeutronData

  integer(ik) :: nd ! number of delayed neutron precursor families
  real(rk), allocatable :: beta(:) ! (nd) [-] Delayed neutron fractions
  real(rk), allocatable :: lamd(:) ! (nd) [1/s] Delayed neutron decay constants

  integer(ik) :: ng ! number of energy groups (for velocities)
  real(rk), allocatable :: vel(:) ! (ng) [cm/s] Neutron velocity

  real(rk) :: deltat ! [s] time-step size
  real(rk) :: tend ! [s] final time

  character(16) :: reference ! reference case (used to specify absorption xs)
endtype

public :: DelayedNeutronData
public :: transient_read, transient_summary, transient_cleanup

contains

  subroutine transient_read(fname, dnd)
    use fileio, only : fileio_open_read
    use exception_handler, only : exception_fatal
    character(*), intent(in) :: fname
    type(DelayedNeutronData), intent(out) :: dnd

    integer, parameter :: iounit = 15
    character(*), parameter :: comment_char = '#'

    character(1024) :: line, card
    integer :: ios

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
      backspace(iounit)

      select case (card)
        case ('nd')
          read(iounit, *) card, dnd%nd
          allocate(dnd%beta(dnd%nd))
          allocate(dnd%lamd(dnd%nd))
          dnd%beta = 0.0_rk
          dnd%lamd = 0.0_rk
        case ('beta')
          read(iounit, *) card, dnd%beta
        case ('lambda')
          read(iounit, *) card, dnd%lamd
        case ('ng')
          read(iounit, *) card, dnd%ng
          allocate(dnd%vel(dnd%ng))
          dnd%vel = 0.0_rk
        case ('vel')
          read(iounit, *) card, dnd%vel
        case ('deltat')
          read(iounit, *) card, dnd%deltat
        case ('tend')
          read(iounit, *) card, dnd%tend
        case ('reference')
          read(iounit, *) card, dnd%reference
        case default
          call exception_fatal('unknown transient input card: ' // &
            trim(adjustl(card)))
      endselect
    enddo

    close(iounit)
  endsubroutine transient_read

  subroutine transient_summary(dnd)
    use output, only : output_write
    type(DelayedNeutronData), intent(in) :: dnd

    character(2**16) :: line, tmp
    integer(ik) :: d, g

    call output_write('=== TRANSIENT ===')

    write(line, '(a,i0)') 'Number of delayed neutron precursor families: ', dnd%nd
    call output_write(line)

    line = ''
    do d = 1,dnd%nd
      write(tmp, '(es9.2)') dnd%beta(d)
      line = trim(adjustl(line)) // ' ' // trim(adjustl(tmp))
    enddo ! d = 1,dnd%nd
    call output_write('beta = ' // line)

    line = ''
    do d = 1,dnd%nd
      write(tmp, '(es9.2)') dnd%lamd(d)
      line = trim(adjustl(line)) // ' ' // trim(adjustl(tmp))
    enddo ! d = 1,dnd%nd
    call output_write('lambda [1/s] = ' // line)

    write(line, '(a,i0)') 'Number of energy groups in transient data: ', dnd%ng
    call output_write(line)
    
    line = ''
    do g = 1,dnd%ng
      write(tmp, '(es9.2)') dnd%vel(g)
      line = trim(adjustl(line)) // ' ' // trim(adjustl(tmp))
    enddo ! g = dnd%ng
    call output_write('velocity [cm/s] = ' // line)

    write(line, '(a,es9.2)') 'Time-step size [s] = ', dnd%deltat
    call output_write(line)
    write(line, '(a,es9.2)') 'Tend [s] = ', dnd%tend
    call output_write(line)
    
    call output_write('Transient reference case: ' &
      // trim(adjustl(dnd%reference)))

    call output_write('')
  endsubroutine transient_summary

  subroutine transient_cleanup(dnd)
    type(DelayedNeutronData), intent(out) :: dnd
    if (allocated(dnd%beta)) then
      deallocate(dnd%beta)
    endif
    if (allocated(dnd%lamd)) then
      deallocate(dnd%lamd)
    endif
    if (allocated(dnd%vel)) then
      deallocate(dnd%vel)
    endif
  endsubroutine transient_cleanup

endmodule transient
