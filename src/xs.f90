module xs
use kind, only : rk, ik
implicit none
private

public :: XSMaterial, XSLibrary, xs_read_library, xs_cleanup

type XSMaterial
  integer(ik) :: ngroup
  integer(ik) :: nmoment
  character(16) :: name

  real(rk), allocatable :: diffusion(:)
  real(rk), allocatable :: transport(:)

  real(rk), allocatable :: sigma_a(:)
  real(rk), allocatable :: sigma_t(:)

  logical :: is_fiss

  real(rk), allocatable :: sigma_f(:)
  real(rk), allocatable :: nusf(:)
  real(rk), allocatable :: chi(:)

  real(rk), allocatable :: scatter(:,:,:) ! (ngroup, ngroup, nmoment)

endtype XSMaterial

type XSLibrary
  character(1024) :: filename
  integer(ik) :: ngroup
  integer(ik) :: niso
  integer(ik) :: nmoment
  type(XSMaterial), allocatable :: mat(:)
endtype XSLibrary

contains

  subroutine xs_read_library(fname, xslib)
    use fileio, only : fileio_open_read
    character(*), intent(in) :: fname
    type(XSLibrary), intent(out) :: xslib

    integer, parameter :: iounit = 13

    character(2**16) :: line, card
    integer :: ios

    integer(ik) :: i, pnt, mom

    pnt = 0

    call fileio_open_read(fname, iounit)
    xslib%filename = fname
    do

      read(iounit, '(a)', iostat=ios) line
      if (ios /= 0) then
        exit
      endif
      if (line == '') then
        cycle
      endif

      read(line, *) card

      card = adjustl(card)
      if (card(1:1) == '#') then
        ! comment
        cycle
      endif

      select case (card)
        case ('ngroup')
          read(line, *) card, xslib%ngroup
        case ('niso')
          read(line, *) card, xslib%niso
        case ('nmoment')
          read(line, *) card, xslib%nmoment
          xslib%nmoment = xslib%nmoment + 1 ! index begins at 1 in SIREN
        case ('name')
          if (pnt == 0) then
            ! allocate necessary space
            ! we must have encountered ngroup, nsio, and nmoment before now
            allocate(xslib%mat(xslib%niso))
            do i = 1,xslib%niso
              xslib%mat(i)%ngroup = xslib%ngroup
              xslib%mat(i)%nmoment = xslib%nmoment
            enddo
          endif
          pnt = pnt + 1
          read(line, *) card, xslib%mat(pnt)%name
          xslib%mat(pnt)%is_fiss = .false. ! initialize
        case ('diffusion')
          allocate(xslib%mat(pnt)%diffusion(xslib%ngroup))
          do i = 1,xslib%ngroup
            read(iounit, *) xslib%mat(pnt)%diffusion(i)
          enddo
          allocate(xslib%mat(pnt)%transport(xslib%ngroup))
          do i = 1,xslib%ngroup
            xslib%mat(pnt)%transport(i) = 1d0 / (3d0 * xslib%mat(pnt)%diffusion(i))
          enddo
        case ('sigma_a')
          allocate(xslib%mat(pnt)%sigma_a(xslib%ngroup))
          do i = 1,xslib%ngroup
            read(iounit, *) xslib%mat(pnt)%sigma_a(i)
          enddo
        case ('sigma_t')
          allocate(xslib%mat(pnt)%sigma_t(xslib%ngroup))
          do i = 1,xslib%ngroup
            read(iounit, *) xslib%mat(pnt)%sigma_t(i)
          enddo
        case ('sigma_f')
          xslib%mat(pnt)%is_fiss = .true.
          allocate(xslib%mat(pnt)%sigma_f(xslib%ngroup))
          do i = 1,xslib%ngroup
            read(iounit, *) xslib%mat(pnt)%sigma_f(i)
          enddo
        case ('nusf')
          xslib%mat(pnt)%is_fiss = .true.
          allocate(xslib%mat(pnt)%nusf(xslib%ngroup))
          do i = 1,xslib%ngroup
            read(iounit, *) xslib%mat(pnt)%nusf(i)
          enddo
        case ('chi')
          xslib%mat(pnt)%is_fiss = .true.
          allocate(xslib%mat(pnt)%chi(xslib%ngroup))
          do i = 1,xslib%ngroup
            read(iounit, *) xslib%mat(pnt)%chi(i)
          enddo
        case ('scatter')
          if (.not. allocated(xslib%mat(pnt)%scatter)) then
            allocate(xslib%mat(pnt)%scatter(xslib%ngroup, xslib%ngroup, xslib%nmoment))
          endif
          read(line, *) card, mom
          do i = 1,xslib%ngroup
            ! index begins at 1 in SIREN
            read(iounit, *) xslib%mat(pnt)%scatter(:,i,mom+1)
          enddo
        case default
          write(*,*) 'Unknown card:', trim(adjustl(card))
          stop 'error in xs_read_library'
      endselect
    enddo

    call xs_summary(xslib)

    close(iounit)
  endsubroutine xs_read_library

  subroutine xs_summary(xslib)
    use output, only : output_write
    type(XSLibrary), intent(in) :: xslib

    character(1024) :: line
    integer(ik) :: i

    call output_write("=== XSLIB ===")
    call output_write('filename: ' // trim(adjustl(xslib%filename)))
    write(line, '(a,i0)') 'niso = ', xslib%niso
    call output_write(line)
    write(line, '(a,i0)') 'ngroup = ', xslib%ngroup
    call output_write(line)
    write(line, '(a,i0)') 'nmoment = ', xslib%nmoment
    call output_write(line)
    do i = 1,xslib%niso
      if ((xslib%mat(i)%is_fiss) .and. &
        ((xslib%ngroup == 2) .or. (xslib%ngroup == 1))) then
        write(line, '(a,a,a,i0,a,l1,a,f8.6)') &
          'name = "', trim(adjustl(xslib%mat(i)%name)), '" , idx = ', i, &
          ' , fiss = ', xslib%mat(i)%is_fiss, ' , kinf = ', calc_kinf(xslib%mat(i))
      else
        write(line, '(a,a,a,i0,a,l1)') &
          'name = "', trim(adjustl(xslib%mat(i)%name)), '" , idx = ', i, ' , fiss = ', xslib%mat(i)%is_fiss
      endif
      call output_write(line)
    enddo
    call output_write('')

  endsubroutine xs_summary

  real(rk) pure function calc_kinf(xsmat)
    type(XSMaterial), intent(in) :: xsmat

    real(rk) :: ratio
    real(rk) :: rem1, rem2

    if (.not. xsmat%is_fiss) then
      calc_kinf = -1.0_rk
    else if (xsmat%ngroup == 1) then
      calc_kinf = xsmat%nusf(1) / (xsmat%sigma_t(1) - xsmat%scatter(1,1,1))
    else if (xsmat%ngroup == 2) then
      rem1 = xsmat%sigma_t(1) - xsmat%scatter(1,1,1)
      rem2 = xsmat%sigma_t(2) - xsmat%scatter(2,2,1)
      ratio = xsmat%scatter(1,2,1) / rem2
      calc_kinf = (xsmat%nusf(1) + xsmat%nusf(2)*ratio) / rem1
    else
      calc_kinf = -1.0_rk
    endif
  endfunction calc_kinf

  subroutine xs_cleanup(xslib)
    type(XSLibrary), intent(out) :: xslib
    if (allocated(xslib%mat)) then
      deallocate(xslib%mat)
    endif
  endsubroutine xs_cleanup

endmodule xs
