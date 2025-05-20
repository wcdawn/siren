module output
use kind
implicit none

private

public :: output_flux_csv, output_power_csv, output_phi_csv, output_transportxs_csv

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

  subroutine output_power_csv(fname, nx, hx, power)
    use fileio, only : fileio_open_write
    character(*), intent(in) :: fname
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: hx
    real(rk), intent(in) :: power(:) ! (nx)

    integer, parameter :: iounit = 17

    integer(ik) :: i

    call fileio_open_write(fname, iounit)

    write(iounit, '(a)') 'x [cm], power'
    do i = 1,nx
      write(iounit, '(es23.16,",",es23.16)') hx*(i-0.5d0), power(i)
    enddo

    close(iounit)
  endsubroutine output_power_csv

  subroutine output_phi_csv(fname, nx, ngroup, pnorder, hx, phi)
    use fileio, only : fileio_open_write
    character(*), intent(in) :: fname
    integer(ik), intent(in) :: nx
    integer(ik), intent(in) :: ngroup
    integer(ik), intent(in) :: pnorder
    real(rk), intent(in) :: hx
    real(rk), intent(in) :: phi(:,:,:) ! (nx,ngroup,pnorder+1)

    integer, parameter :: iounit = 17

    integer(ik) :: i, g, n

    call fileio_open_write(fname, iounit)

    write(iounit, '(a)', advance='no') 'x [cm]'
    do n = 1,pnorder+1
      do g = 1,ngroup
        write(iounit, '(a,a,i0,a,i0)', advance='no') ',', 'phi_n', n-1, '_g', g
      enddo
    enddo
    write(iounit,*)

    do i = 1,nx
      write(iounit, '(es23.16)', advance='no') hx*(i-0.5d0) ! cell center
      do n = 1,pnorder+1
        do g = 1,ngroup
          write(iounit, '(a,es23.16)', advance='no') ',', phi(i,g,n)
        enddo
      enddo
      write(iounit,*)
    enddo

    close(iounit)
  endsubroutine output_phi_csv

  subroutine output_transportxs_csv(fname, nx, ngroup, pnorder, hx, sigma_tr)
    use fileio, only : fileio_open_write
    character(*), intent(in) :: fname
    integer(ik), intent(in) :: nx
    integer(ik), intent(in) :: ngroup
    integer(ik), intent(in) :: pnorder
    real(rk), intent(in) :: hx
    real(rk), intent(in) :: sigma_tr(:,:,:) ! (nx,ngroup,pnorder+1)

    integer, parameter :: iounit = 17

    integer(ik) :: i, g, n

    call fileio_open_write(fname, iounit)

    write(iounit, '(a)', advance='no') 'x [cm]'
    do n = 1,pnorder+1
      do g = 1,ngroup
        write(iounit, '(a,a,i0,a,i0)', advance='no') ',', 'sigma_tr_n', n-1, '_g', g
      enddo
    enddo
    write(iounit,*)

    do i = 1,nx
      write(iounit, '(es23.16)', advance='no') hx*(i-0.5d0) ! cell center
      do n = 1,pnorder+1
        do g = 1,ngroup
          write(iounit, '(a,es23.16)', advance='no') ',', sigma_tr(i,g,n)
        enddo
      enddo
      write(iounit,*)
    enddo

    close(iounit)
  endsubroutine output_transportxs_csv

endmodule output
