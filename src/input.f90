module input
use kind, only : rk, ik
implicit none (external)

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

! optional boundary conditions
! only allowed for very special problems
character(16) :: boundary_left = 'mirror', boundary_right = 'mirror'

! optional analytic reference solution
character(16) :: analytic_reference = 'none'

! optional solver option for pn solver
! right now, just chooses iteration order either classical "flip" or my version of "lupine"
character(16) :: pn_solver_opt = 'lupine'

! optional solver option for diffusion solver (to be extended to pn)
! choose whether to solve one-group at-a-time ("onegroup") or all at once ("block")
character(16) :: energy_solver_opt = 'onegroup'

! optional high-fidelity -> low-fidelity option to run diffusion after transport
logical :: high_low = .false.

! optional recalculate diffusion coefficient to 1/(3*sigma_t) for consistency with P1
logical :: force_consistent_diffusion = .false.

public :: input_read, input_cleanup

public :: nx, dx, xslib_fname, mat_map
public :: k_tol, phi_tol, max_iter
public :: pnorder
public :: refine
public :: boundary_left, boundary_right
public :: analytic_reference
public :: pn_solver_opt
public :: energy_solver_opt
public :: high_low, force_consistent_diffusion

contains

  subroutine input_read(fname)
    use fileio, only : fileio_open_read
    use exception_handler, only : exception_fatal
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
        case ('analytic_reference')
          read(line, *) card, analytic_reference
        case ('pn_solver_opt')
          read(line, *) card, pn_solver_opt
        case ('energy_solver_opt')
          read(line, *) card, energy_solver_opt
        case ('high_low')
          read(line, *) card, high_low
        case ('force_consistent_diffusion')
          read(line, *) card, force_consistent_diffusion
        case default
          call exception_fatal('unknown input card: ' // trim(adjustl(card)))
      endselect
    enddo

    call input_summary()

    call input_check()

    close(iounit)
  endsubroutine input_read

  subroutine input_check()
    use output, only : output_write
    use exception_handler, only : exception_fatal, exception_warning
    if (trim(adjustl(boundary_left)) /= 'mirror') then
      call exception_fatal('boundary_left must be set to mirror (for now)')
    elseif ((pnorder /= 0) .and. (trim(adjustl(boundary_right)) == 'zero')) then
      call exception_warning('Zero-flux boundary condition with PN transport ' // &
        'is probably not what you want. ' // &
        'This will set the scalar flux to identically zero at the boundary. ' // &
        'It is intended only for numerical benchmarking.')
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
    call output_write('analytic_reference = *' // trim(adjustl(analytic_reference)) // '*')

    if (pnorder /= 0) then
      call output_write('pn_solver_opt = *' // trim(adjustl(pn_solver_opt)) // '*')
    endif

    call output_write('energy_solver_opt = *' // trim(adjustl(energy_solver_opt)) // '*')

    write(line, '(a,l2)') 'high_low = ', high_low
    call output_write(line)
    write(line, '(a,l2)') 'force_consistent_diffusion = ', force_consistent_diffusion
    call output_write(line)

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
