module output
use kind
use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                          stdout=>output_unit, &
                                          stderr=>error_unit
implicit none

private

public :: output_open_file, output_close_file, output_write, &
  output_flux_csv, output_power_csv, output_phi_csv, output_transportxs_csv

integer, parameter, private :: output_file_unit = 99
integer, parameter, private :: output_list(2) = [ stdout, output_file_unit ]

contains

  subroutine output_open_file(fname)
    use fileio, only : fileio_open_write
    character(*), intent(in) :: fname
    call fileio_open_write(fname, output_file_unit)
  endsubroutine output_open_file

  subroutine output_close_file()
    close(output_file_unit)
  endsubroutine output_close_file

  subroutine output_write(str)
    character(*), intent(in) :: str
    integer :: i
    do i = 1,size(output_list)
      write(output_list(i), '(a)') trim(adjustl(str))
    enddo ! i = 1,size(output_list)
  endsubroutine output_write

  subroutine output_flux_csv(fname, nx, ng, dx, flux)
    use fileio, only : fileio_open_write
    character(*), intent(in) :: fname
    integer(ik), intent(in) :: nx, ng
    real(rk), intent(in) :: dx(:) ! (nx)
    real(rk), intent(in) :: flux(:,:) ! (nx,ng)

    integer, parameter :: iounit = 17

    integer(ik) :: i, g
    real(rk) :: xx

    call fileio_open_write(fname, iounit)

    write(iounit, '(a)', advance='no') 'x [cm]'
    do g = 1,ng
      write(iounit, '(a,a,i0)', advance='no') ',', 'flux_g', g
    enddo
    write(iounit,*)

    xx = 0d0
    do i = 1,nx
      write(iounit, '(es23.16)', advance='no') xx + 0.5d0*dx(i) ! cell center
      do g = 1,ng
        write(iounit, '(a,es23.16)', advance='no') ',', flux(i,g)
      enddo
      write(iounit,*)
      xx = xx + dx(i)
    enddo

    close(iounit)
  endsubroutine output_flux_csv

  subroutine output_power_csv(fname, nx, dx, power)
    use fileio, only : fileio_open_write
    character(*), intent(in) :: fname
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    real(rk), intent(in) :: power(:) ! (nx)

    integer, parameter :: iounit = 17

    integer(ik) :: i
    real(rk) :: xx

    call fileio_open_write(fname, iounit)

    write(iounit, '(a)') 'x [cm], power'
    xx = 0d0
    do i = 1,nx
      write(iounit, '(es23.16,",",es23.16)') xx + 0.5d0*dx(i), power(i)
      xx = xx + dx(i)
    enddo

    close(iounit)
  endsubroutine output_power_csv

  subroutine output_phi_csv(fname, nx, ngroup, pnorder, dx, phi)
    use fileio, only : fileio_open_write
    character(*), intent(in) :: fname
    integer(ik), intent(in) :: nx
    integer(ik), intent(in) :: ngroup
    integer(ik), intent(in) :: pnorder
    real(rk), intent(in) :: dx(:) ! (nx)
    real(rk), intent(in) :: phi(:,:,:) ! (nx,ngroup,pnorder+1)

    integer, parameter :: iounit = 17

    integer(ik) :: i, g, n
    real(rk) :: xx

    call fileio_open_write(fname, iounit)

    write(iounit, '(a)', advance='no') 'x [cm]'
    do n = 1,pnorder+1
      do g = 1,ngroup
        write(iounit, '(a,a,i0,a,i0)', advance='no') ',', 'phi_n', n-1, '_g', g
      enddo
    enddo
    write(iounit,*)

    xx = 0d0
    do i = 1,nx
      write(iounit, '(es23.16)', advance='no') xx + 0.5d0*dx(i) ! cell center
      do n = 1,pnorder+1
        do g = 1,ngroup
          write(iounit, '(a,es23.16)', advance='no') ',', phi(i,g,n)
        enddo
      enddo
      write(iounit,*)
      xx = xx + dx(i)
    enddo

    close(iounit)
  endsubroutine output_phi_csv

  subroutine output_transportxs_csv(fname, nx, ngroup, pnorder, dx, sigma_tr)
    use fileio, only : fileio_open_write
    character(*), intent(in) :: fname
    integer(ik), intent(in) :: nx
    integer(ik), intent(in) :: ngroup
    integer(ik), intent(in) :: pnorder
    real(rk), intent(in) :: dx(:) ! (nx)
    real(rk), intent(in) :: sigma_tr(:,:,:) ! (nx,ngroup,pnorder+1)

    integer, parameter :: iounit = 17

    integer(ik) :: i, g, n
    real(rk) :: xx

    call fileio_open_write(fname, iounit)

    write(iounit, '(a)', advance='no') 'x [cm]'
    do n = 1,pnorder+1
      do g = 1,ngroup
        write(iounit, '(a,a,i0,a,i0)', advance='no') ',', 'sigma_tr_n', n-1, '_g', g
      enddo
    enddo
    write(iounit,*)

    xx = 0d0
    do i = 1,nx
      write(iounit, '(es23.16)', advance='no') xx + 0.5d0*dx(i) ! cell center
      do n = 1,pnorder+1
        do g = 1,ngroup
          write(iounit, '(a,es23.16)', advance='no') ',', sigma_tr(i,g,n)
        enddo
      enddo
      write(iounit,*)
      xx = xx + dx(i)
    enddo

    close(iounit)
  endsubroutine output_transportxs_csv

endmodule output
