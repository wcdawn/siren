program siren
use kind
use xs, only : XSLibrary, xs_read_library, xs_cleanup
use input, only : input_read, input_cleanup, &
  xslib_fname, refine, nx, hx, mat_map, &
  k_tol, phi_tol, max_iter
use geometry, only : uniform_refinement
use diffusion, only : diffusion_power_iteration
use output, only : output_flux_csv, output_power_csv
use power, only : power_calculate
implicit none

integer(ik) :: i
character(1024) :: input_fname
type(XSLibrary) :: xs

real(rk) :: keff
real(rk), allocatable :: flux(:,:)
real(rk), allocatable :: power(:)

write(*,*) 'begin SIREN'

if (command_argument_count() == 0) then
  stop 'missing input filename'
endif

input_fname = ''
call get_command_argument(1, input_fname)
write(*, '(a,a)') 'input file: ', trim(adjustl(input_fname))

call input_read(input_fname)
call xs_read_library(xslib_fname, xs)

if (refine > 0) then
  do i = 1,refine
    call uniform_refinement(nx, hx, mat_map)
  enddo
  write(*,*) '=== REFINEMENT ==='
  write(*,'(a,i0)') 'levels = ', refine
  write(*,'(a,i0)') 'refined nx = ', nx
  write(*,'(a,es13.6)') 'refined hx = ', hx
endif

allocate(flux(nx, xs%ngroup))
call diffusion_power_iteration(nx, hx, mat_map, xs, k_tol, phi_tol, max_iter, keff, flux)

write(*,'(a,f22.20)') 'keff = ', keff
write(*,*)

write(*,*) 'writing flux.csv'
call output_flux_csv('flux.csv', nx, xs%ngroup, hx, flux)

write(*,*) 'writing power.csv'
allocate(power(nx))
call power_calculate(nx, mat_map, xs, flux, power)
call output_power_csv('power.csv', nx, hx, power)

deallocate(power)
deallocate(flux)

call xs_cleanup(xs)
call input_cleanup()

write(*,*) 'end SIREN'

endprogram siren
