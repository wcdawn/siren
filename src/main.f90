program siren
use kind
use xs, only : XSLibrary, xs_read_library, xs_cleanup
use input, only : input_read, input_cleanup, &
  xslib_fname, refine, nx, dx, mat_map, pnorder, boundary_right, &
  k_tol, phi_tol, max_iter, analytic_reference, pn_solver_opt, energy_solver_opt
use geometry, only : geometry_uniform_refinement, geometry_summary
use diffusion, only : diffusion_power_iteration
use diffusion_block, only : diffusion_block_power_iteration
use transport, only : sigma_tr, transport_power_iteration, transport_power_iteration_flip
use transport_block, only : transport_block_power_iteration
use output, only : output_open_file, output_close_file, output_write, &
  output_flux_csv, output_power_csv, output_phi_csv, output_transportxs_csv, output_matmap_csv
use power, only : power_calculate
use analytic, only : analytic_error
use exception_handler
use timer
implicit none

integer(ik) :: i
character(1024) :: input_fname
character(1024) :: fname_stub, &
  fname_flux, fname_phi, fname_power, fname_transportxs, fname_out, fname_analytic, &
  fname_matmap
character(1024) :: line
type(XSLibrary) :: xs

real(rk) :: keff
real(rk), allocatable :: phi(:,:,:) ! (nx, ngroup, pnorder)
real(rk), allocatable :: power(:)

if (command_argument_count() == 0) then
  stop 'missing input filename'
endif

call timer_init()

input_fname = ''
call get_command_argument(1, input_fname)

i = index(trim(adjustl(input_fname)), '.', back=.true.)
i = i - 1
fname_stub = input_fname(:i)
fname_out = trim(adjustl(fname_stub)) // '.out'
fname_flux = trim(adjustl(fname_stub)) // '_flux.csv'
fname_phi = trim(adjustl(fname_stub)) // '_phi.csv'
fname_power = trim(adjustl(fname_stub)) // '_power.csv'
fname_transportxs = trim(adjustl(fname_stub)) // '_transportxs.csv'
fname_analytic = trim(adjustl(fname_stub)) // '_analytic.csv'
fname_matmap = trim(adjustl(fname_stub)) // '_matmap.csv'

call output_open_file(trim(adjustl(fname_out)))

call output_write('begin SIREN')
call output_write('input file: ' // trim(adjustl(input_fname)))
call output_write('')

call timer_start('input_read')
call input_read(input_fname)
call output_write('(before refinement)')
call geometry_summary(nx, dx, mat_map)
call timer_stop('input_read')

call timer_start('xs_read')
call xs_read_library(xslib_fname, xs)
call timer_stop('xs_read')

if (refine > 0) then
  call timer_start('refinement')
  do i = 1,refine
    call geometry_uniform_refinement(nx, dx, mat_map)
  enddo
  call output_write('(after refinement)')
  call geometry_summary(nx, dx, mat_map)
  call timer_stop('refinement')
endif

call timer_start('solver')
allocate(phi(nx, xs%ngroup, pnorder+1))
if (pnorder == 0) then
  select case (energy_solver_opt)
    case ('block')
      call diffusion_block_power_iteration(nx, dx, mat_map, xs, boundary_right, k_tol, phi_tol, max_iter, keff, phi(:,:,1))
    case ('onegroup')
      call diffusion_power_iteration(nx, dx, mat_map, xs, boundary_right, k_tol, phi_tol, max_iter, keff, phi(:,:,1))
    case default
      call exception_fatal('unknown energy_solver_opt: ' // trim(adjustl(energy_solver_opt)))
  endselect
else
  select case (energy_solver_opt)
    case ('block')
      call transport_block_power_iteration(nx, dx, mat_map, xs, boundary_right, k_tol, phi_tol, max_iter, pnorder, keff, phi)
    case ('onegroup')
      select case (pn_solver_opt)
        case ('flip')
          call transport_power_iteration_flip(nx, dx, mat_map, xs, boundary_right, k_tol, phi_tol, max_iter, pnorder, keff, phi)
        case ('lupine')
          call transport_power_iteration(nx, dx, mat_map, xs, boundary_right, k_tol, phi_tol, max_iter, pnorder, keff, phi)
        case default
          call exception_fatal('unknown pn_solver_opt: ' // trim(adjustl(pn_solver_opt)))
      endselect
    case default
      call exception_fatal('unknown energy_solver_opt: ' // trim(adjustl(energy_solver_opt)))
  endselect
endif
call timer_stop('solver')

write(line, '(a,f22.20)') 'keff = ', keff
call output_write(line)
call output_write('')

if (analytic_reference /= 'none') then
  call analytic_error(analytic_reference, fname_analytic, nx, xs%ngroup, pnorder, xs, dx, phi, keff)
endif

call timer_start('output')
call output_write('writing ' // trim(adjustl(fname_flux)))
call output_flux_csv(trim(adjustl(fname_flux)), nx, xs%ngroup, dx, phi(:,:,1))

call output_write('writing ' // trim(adjustl(fname_phi)))
call output_phi_csv(trim(adjustl(fname_phi)), nx, xs%ngroup, pnorder, dx, phi)

call output_write('writing ' // trim(adjustl(fname_power)))
allocate(power(nx))
call power_calculate(nx, mat_map, xs, phi(:,:,1), power)
call output_power_csv(trim(adjustl(fname_power)), nx, dx, power)

if (allocated(sigma_tr)) then
  call output_write('writing ' // trim(adjustl(fname_transportxs)))
  call output_transportxs_csv(trim(adjustl(fname_transportxs)), nx, xs%ngroup, pnorder, dx, sigma_tr)
  deallocate(sigma_tr)
endif

call output_write('writing ' // trim(adjustl(fname_matmap)))
call output_matmap_csv(trim(adjustl(fname_matmap)), nx, dx, mat_map, xs%niso, xs%mat(:)%name)

call output_write('')
call timer_stop('output')

call timer_summary()
call exception_summary()

deallocate(power)
deallocate(phi)

call xs_cleanup(xs)
call input_cleanup()

call output_write('Normal Termination :)')
call output_write('end SIREN')

call output_close_file()

endprogram siren
