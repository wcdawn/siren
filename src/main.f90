program siren
use kind, only : rk, ik
use xs, only : XSLibrary, XSMaterial, xs_read_library, xs_cleanup
use input, only : input_read, input_cleanup, &
  xslib_fname, refine, nx, dx, mat_map, pnorder, boundary_right, &
  k_tol, phi_tol, max_iter, analytic_reference, pn_solver_opt, energy_solver_opt, &
  linear_solver_opt, krylov_max_iter, krylov_atol, krylov_rtol
use geometry, only : geometry_uniform_refinement, geometry_summary
use diffusion, only : diffusion_power_iteration
use diffusion_block, only : diffusion_block_power_iteration
use transport, only : transport_power_iteration, transport_power_iteration_flip
use transport_block, only : transport_block_power_iteration
use output, only : output_open_file, output_close_file, output_write, &
  output_flux_csv, output_power_csv, output_phi_csv, output_transportxs_csv, output_matmap_csv
use power, only : power_calculate
use analytic, only : analytic_error
use exception_handler, only : exception_fatal, exception_summary
use timer, only : timer_init, timer_start, timer_stop, timer_summary
implicit none

integer(ik) :: i
character(1024) :: input_fname
character(1024) :: fname_stub, &
  fname_flux, fname_phi, fname_power, fname_transportxs, fname_out, fname_analytic, &
  fname_matmap
character(1024) :: line
type(XSLibrary) :: xs

real(rk) :: keff
real(rk), allocatable :: sigma_tr(:,:,:) ! (nx, ngroup, pnorder+1)
real(rk), allocatable :: phi(:,:,:) ! (nx, ngroup, pnorder+1)
real(rk), allocatable :: power(:)

character(len=:), allocatable :: mat_name_list(:)

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

if ((pnorder /= 0) .and. (xs%nmoment /= 1) .and. (trim(adjustl(energy_solver_opt)) /= 'block')) then
  call exception_fatal(&
    'The onegroup PN solver is unstable for true anisotropic scattering. ' // &
    'You must use "energy_solver_opt block"')
endif

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
      call diffusion_block_power_iteration( &
        nx, dx, mat_map, xs, boundary_right, k_tol, phi_tol, max_iter, &
        keff, phi(:,:,1))
    case ('onegroup')
      call diffusion_power_iteration( &
        nx, dx, mat_map, xs, boundary_right, k_tol, phi_tol, max_iter, &
        linear_solver_opt, krylov_max_iter, krylov_atol, krylov_rtol, &
        keff, phi(:,:,1))
    case default
      call exception_fatal('unknown energy_solver_opt: ' // trim(adjustl(energy_solver_opt)))
  endselect
else
  allocate(sigma_tr(nx, xs%ngroup, pnorder+1))
  select case (energy_solver_opt)
    case ('block')
      call transport_block_power_iteration( &
        nx, dx, mat_map, xs, boundary_right, k_tol, phi_tol, max_iter, pnorder, &
        keff, sigma_tr, phi)
    case ('onegroup')
      select case (pn_solver_opt)
        case ('flip')
          call transport_power_iteration_flip( &
            nx, dx, mat_map, xs, boundary_right, k_tol, phi_tol, max_iter, pnorder, &
            keff, sigma_tr, phi)
        case ('lupine')
          call transport_power_iteration(&
            nx, dx, mat_map, xs, boundary_right, k_tol, phi_tol, max_iter, pnorder, &
            keff, sigma_tr, phi)
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
  call output_transportxs_csv( &
    trim(adjustl(fname_transportxs)), nx, xs%ngroup, pnorder, dx, sigma_tr)
  deallocate(sigma_tr)
endif

call output_write('writing ' // trim(adjustl(fname_matmap)))
allocate(character(len=len(xs%mat(1)%name)) :: mat_name_list(xs%niso))
do i = 1,xs%niso
  mat_name_list(i) = xs%mat(i)%name
enddo
call output_matmap_csv(trim(adjustl(fname_matmap)), nx, dx, mat_map, xs%niso, mat_name_list)
deallocate(mat_name_list)

call output_write('')
call timer_stop('output')

call timer_summary()
call exception_summary()

deallocate(power)
deallocate(phi)
if (allocated(sigma_tr)) deallocate(sigma_tr)

call xs_cleanup(xs)
call input_cleanup()

call output_write('Normal Termination :)')
call output_write('end SIREN')

call output_close_file()

endprogram siren
