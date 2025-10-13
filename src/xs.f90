module xs
use kind
use exception_handler
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
            xslib%mat(pnt)%scatter = 0.0_rk
          endif
          read(line, *) card, mom
          do i = 1,xslib%ngroup
            ! index begins at 1 in SIREN
            read(iounit, *) xslib%mat(pnt)%scatter(:,i,mom+1)
          enddo
        case default
          call exception_fatal('Error in xs_read_library. Unknown card: ' // trim(adjustl(card)))
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
      if (xslib%mat(i)%is_fiss) then
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

  real(rk) function calc_kinf(xsmat)
    use linalg, only : geneig
    type(XSMaterial), intent(in) :: xsmat

    real(rk), allocatable :: a(:,:), f(:,:)
    complex(rk), allocatable :: eigval(:)
    real(rk), allocatable :: eigvec(:,:)
    real(rk), allocatable :: phi(:)
    integer(ik) :: g, gprime

    if (.not. xsmat%is_fiss) then
      calc_kinf = -1.0_rk
    else
      allocate(a(xsmat%ngroup,xsmat%ngroup), f(xsmat%ngroup,xsmat%ngroup))
      allocate(eigval(xsmat%ngroup), eigvec(xsmat%ngroup,xsmat%ngroup))

      do g = 1,xsmat%ngroup
        a(g,g) = xsmat%sigma_t(g)
      enddo ! g = 1,xsmat%ngroup
      a = a - transpose(xsmat%scatter(:,:,1))

      do g = 1,xsmat%ngroup
        do gprime = 1,xsmat%ngroup
          f(gprime,g) = xsmat%chi(gprime) * xsmat%nusf(g)
        enddo ! gprime = 1,xsmat%ngroup
      enddo ! g = 1,xsmat%ngroup

      call geneig(xsmat%ngroup, a, f, eigval, eigvec)

      calc_kinf = huge(1.0_rk)
      gprime = -1
      do g = 1,xsmat%ngroup
        if ((abs(eigval(g)) < calc_kinf) .and. (abs(eigval(g)) > epsilon(1.0_rk))) then
          gprime = g
          calc_kinf = 1.0_rk/abs(eigval(g))
        endif
      enddo

      ! NOTE: this isn't used now, but may be useful in the future
      allocate(phi(xsmat%ngroup))
      phi = eigvec(:,gprime)
      deallocate(phi)

      deallocate(a, f)
      deallocate(eigval, eigvec)
    endif
  endfunction calc_kinf

  subroutine xs_cleanup(xslib)
    type(XSLibrary), intent(out) :: xslib
    if (allocated(xslib%mat)) then
      deallocate(xslib%mat)
    endif
  endsubroutine xs_cleanup

endmodule xs
