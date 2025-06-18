program siren
use kind
use xs, only : XSLibrary, xs_read_library, xs_cleanup
use input, only : input_read, input_cleanup, &
  xslib_fname, refine, nx, hx, mat_map, pnorder, &
  k_tol, phi_tol, max_iter
use geometry, only : uniform_refinement
use diffusion, only : diffusion_power_iteration
use transport, only : sigma_tr, transport_power_iteration
use output, only : output_open_file, output_close_file, output_write, &
  output_flux_csv, output_power_csv, output_phi_csv, output_transportxs_csv
use power, only : power_calculate
implicit none

integer(ik) :: i, g
character(1024) :: input_fname
character(1024) :: fname_stub, fname_flux, fname_phi, fname_power, fname_transportxs, fname_out
character(1024) :: line
type(XSLibrary) :: xs

real(rk) :: keff
real(rk), allocatable :: phi(:,:,:) ! (nx, ngroup, pnorder)
real(rk), allocatable :: power(:)

if (command_argument_count() == 0) then
  stop 'missing input filename'
endif


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

call output_open_file(trim(adjustl(fname_out)))

call output_write('begin SIREN')
call output_write('input file: ' // trim(adjustl(input_fname)))

call input_read(input_fname)
call xs_read_library(xslib_fname, xs)

if (refine > 0) then
  do i = 1,refine
    call uniform_refinement(nx, hx, mat_map)
  enddo
  call output_write('=== REFINEMENT ===')
  write(line, '(a,i0)') 'levels = ', refine
  call output_write(line)
  write(line, '(a,i0)') 'refined nx = ', nx
  call output_write(line)
  write(line, '(a,es13.6)') 'refined hx = ', hx
  call output_write(line)
  call output_write('')
endif

allocate(phi(nx, xs%ngroup, pnorder+1))
if (pnorder == 0) then

  do i = 1,xs%niso
    do g = 1,xs%ngroup
      xs%mat(i)%diffusion(g) = 1d0/(3d0*xs%mat(i)%sigma_t(g))
    enddo
  enddo

  call diffusion_power_iteration(nx, hx, mat_map, xs, k_tol, phi_tol, max_iter, keff, phi(:,:,1))
else
  !call transport_outer_iteration(nx, hx, mat_map, xs, k_tol, phi_tol, max_iter, pnorder, keff, phi)

  phi = 1d0
  keff = 1d0
  call transport_power_iteration(nx, hx, mat_map, xs, k_tol, phi_tol, max_iter, pnorder, keff, phi)

endif

write(line, '(a,f22.20)') 'keff = ', keff
call output_write(line)
call output_write('')

call output_write('writing ' // trim(adjustl(fname_flux)))
call output_flux_csv(trim(adjustl(fname_flux)), nx, xs%ngroup, hx, phi(:,:,1))

call output_write('writing ' // trim(adjustl(fname_phi)))
call output_phi_csv(trim(adjustl(fname_phi)), nx, xs%ngroup, pnorder, hx, phi)

call output_write('writing ' // trim(adjustl(fname_power)))
allocate(power(nx))
call power_calculate(nx, mat_map, xs, phi(:,:,1), power)
call output_power_csv(trim(adjustl(fname_power)), nx, hx, power)

if (allocated(sigma_tr)) then
  call output_write('writing ' // trim(adjustl(fname_transportxs)))
  call output_transportxs_csv(trim(adjustl(fname_transportxs)), nx, xs%ngroup, pnorder, hx, sigma_tr)
  deallocate(sigma_tr)
endif

deallocate(power)
deallocate(phi)

call xs_cleanup(xs)
call input_cleanup()

call output_write('end SIREN')

call output_close_file()

endprogram siren
