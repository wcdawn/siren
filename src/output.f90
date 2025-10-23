module output
use kind, only : rk, ik
use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                          stdout=>output_unit, &
                                          stderr=>error_unit
implicit none

private

public :: output_open_file, output_close_file, output_write, &
  output_flux_csv, output_power_csv, output_phi_csv, output_transportxs_csv, &
  output_matmap_csv, output_trid_matrix_csv, output_vec

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

  subroutine output_matmap_csv(fname, nx, dx, mat_map, niso, xsname)
    use fileio, only : fileio_open_write
    character(*), intent(in) :: fname
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    integer(ik), intent(in) :: niso
    character(*), intent(in) :: xsname(:) ! niso

    integer, parameter :: iounit = 17

    integer(ik) :: i
    real(rk) :: xlo, xhi

    call fileio_open_write(fname, iounit)

    write(iounit, '(a,i0)') 'niso ', niso
    do i = 1,niso
      write(iounit, '(a)') trim(adjustl(xsname(i)))
    enddo

    write(iounit, '(a)') 'xstart [cm] , xend [cm] , matid'
    xlo = 0d0
    do i = 1,nx
      xhi = xlo + dx(i)
      write(iounit, '(es23.16," , ",es23.16," , ",i0)') xlo, xhi, mat_map(i)
      xlo = xhi
    enddo

    close(iounit)
  endsubroutine output_matmap_csv

  subroutine output_trid_matrix_csv(fname, n, sub, dia, sup)
    use fileio, only : fileio_open_write
    character(*), intent(in) :: fname
    integer(ik), intent(in) :: n
    real(rk), intent(in) :: sub(:) ! (n-1)
    real(rk), intent(in) :: dia(:) ! (n)
    real(rk), intent(in) :: sup(:) ! (n-1)

    integer(ik) :: i

    integer, parameter :: iounit = 17

    call fileio_open_write(fname, iounit)

    write(iounit, '(a)') 'sub , dia , sup'
    write(iounit, '(es23.16, " , ", es23.16, " , ", es23.16)') &
      0.0_rk, dia(1), sup(1)
    do i = 2,n-2
      write(iounit, '(es23.16, " , ", es23.16, " , ", es23.16)') &
        sub(i-1), dia(i), sup(i)
    enddo ! i = 2,n-1
    write(iounit, '(es23.16, " , ", es23.16, " , ", es23.16)') &
      sub(n-1), dia(n), 0.0_rk

    close(iounit)
  endsubroutine output_trid_matrix_csv

  subroutine output_vec(fname, n, vec)
    use fileio, only : fileio_open_write
    character(*), intent(in) :: fname
    integer(ik), intent(in) :: n
    real(rk), intent(in) :: vec(:) ! (n)

    integer(ik) :: i

    integer, parameter :: iounit = 17

    call fileio_open_write(fname, iounit)

    write(iounit, '(a)') 'vec'
    do i = 1,n
      write(iounit, '(es23.16)') vec(i)
    enddo ! i = 1,n

    close(iounit)
  endsubroutine output_vec

endmodule output
