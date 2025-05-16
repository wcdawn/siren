module output
use kind
implicit none

private

public :: output_flux_csv

contains

  subroutine output_flux_csv(fname, nx, ng, hx, flux)
    use fileio, only : fileio_open_write
    character(*), intent(in) :: fname
    integer(ik), intent(in) :: nx, ng
    real(rk), intent(in) :: hx
    real(rk), intent(in) :: flux(:,:) ! (nx,ng)

    integer, parameter :: iounit = 17

    integer(ik) :: i, g

    call fileio_open_write(fname, iounit)

    write(iounit, '(a)', advance='no') 'x [cm]'
    do g = 1,ng
      write(iounit, '(a,a,i0)', advance='no') ',', 'flux_g', g
    enddo
    write(iounit,*)

    do i = 1,nx
      write(iounit, '(es23.16)', advance='no') hx*(i-0.5d0) ! cell center
      do g = 1,ng
        write(iounit, '(a,es23.16)', advance='no') ',', flux(i,g)
      enddo
      write(iounit,*)
    enddo

    close(iounit)
  endsubroutine output_flux_csv

endmodule output
