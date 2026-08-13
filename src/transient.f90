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
public :: transient_read, transient_summary, transient_solve, transient_cleanup

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

  subroutine transient_init_precursors(nx, mat_map, xslib, dnd, flux, prec)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    type(DelayedNeutronData), intent(in) :: dnd
    real(rk), intent(in) :: flux(:,:) ! (nx,ngroup)
    real(rk), intent(out) :: prec(:,:) ! (nx,nd)
  endsubroutine transient_init_precursors

  subroutine transient_solve(nx, dx, mat_map, xslib, dnd, boundary_right, phi_tol, max_iter, keff, flux)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    type(DelayedNeutronData), intent(in) :: dnd
    character(*), intent(in) :: boundary_right
    real(rk), intent(in) :: phi_tol
    integer(ik), intent(in) :: max_iter
    real(rk), intent(in) :: keff
    real(rk), intent(inout) :: flux(:,:) ! (nx,ngroup)

    real(rk), allocatable :: sub(:,:), dia(:,:), sup(:,:) ! (nx,ng)
    real(rk), allocatable :: sub_copy(:), dia_copy(:), sup_copy(:) ! (nx)
    real(rk), allocatable :: fsource(:,:), upsource(:,:), downsource(:) ! (nx,ng), (nx,ng), (nx)
    real(rk), allocatable :: q(:) ! (nx)

    real(rk), allocatable :: prec(:,:) ! (nx,nd)

    real(rk), parameter :: omega = 1.9_rk ! over-relaxation

    allocate(sub(nx,xslib%ngroup))
    allocate(dia(nx,xslib%ngroup))
    allocate(sup(nx,xslib%ngroup))
    allocate(sub_copy(nx))
    allocate(dia_copy(nx))
    allocate(sup_copy(nx))
    sub = 0.0_rk
    dia = 0.0_rk
    sup = 0.0_rk
    sub_copy = 0.0_rk
    dia_copy = 0.0_rk
    sup_copy = 0.0_rk

    allocate(fsource(nx,xslib%ngroup))
    allocate(upsource(nx,xslib%ngroup))
    allocate(downsource(nx))
    allocate(q(nx))
    fsource = 0.0_rk
    upsource = 0.0_rk
    downsource = 0.0_rk
    q = 0.0_rk

    allocate(prec(nx,dnd%nd))
    prec = 0.0_rk
    call transient_init_precursors(nx, mat_map, xslib, dnd, flux, prec)

    deallocate(sub, dia, sup)
    deallocate(sub_copy, dia_copy, sup_copy)
    deallocate(fsource, upsource, downsource, q)
    deallocate(prec)
  endsubroutine transient_solve

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
