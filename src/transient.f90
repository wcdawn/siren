module transient
use kind, only : rk, ik
implicit none

private

type DelayedNeutronData

  integer(ik) :: nd ! number of delayed neutron precursor families
  real(rk), allocatable :: beta(:) ! (nd) [-] Delayed neutron fractions
  real(rk), allocatable :: lamd(:) ! (nd) [1/s] Delayed neutron decay constants

  integer(ik) :: ng ! number of energy groups (for velocities)
  real(rk), allocatable :: vel(:) ! (ng) [cm/s] Neutron velocity

  real(rk) :: deltat ! [s] time-step size
  real(rk) :: tend ! [s] final time

  character(16) :: reference ! reference case (used to specify absorption xs)
endtype

public :: DelayedNeutronData
public :: transient_read, transient_summary, transient_solve, transient_cleanup

contains

  subroutine transient_read(fname, dnd)
    use fileio, only : fileio_open_read
    use exception_handler, only : exception_fatal
    character(*), intent(in) :: fname
    type(DelayedNeutronData), intent(out) :: dnd

    integer, parameter :: iounit = 15
    character(*), parameter :: comment_char = '#'

    character(1024) :: line, card
    integer :: ios

    call fileio_open_read(trim(adjustl(fname)), iounit)

    do
      read(iounit, '(a)', iostat=ios) line
      if (ios /= 0) then
        exit
      endif

      if (line == '') then
        cycle
      endif

      line = trim(adjustl(line))
      if (line(1:1) == comment_char) then
        cycle
      endif

      read(line, *) card
      backspace(iounit)

      select case (card)
        case ('nd')
          read(iounit, *) card, dnd%nd
          allocate(dnd%beta(dnd%nd))
          allocate(dnd%lamd(dnd%nd))
          dnd%beta = 0.0_rk
          dnd%lamd = 0.0_rk
        case ('beta')
          read(iounit, *) card, dnd%beta
        case ('lambda')
          read(iounit, *) card, dnd%lamd
        case ('ng')
          read(iounit, *) card, dnd%ng
          allocate(dnd%vel(dnd%ng))
          dnd%vel = 0.0_rk
        case ('vel')
          read(iounit, *) card, dnd%vel
        case ('deltat')
          read(iounit, *) card, dnd%deltat
        case ('tend')
          read(iounit, *) card, dnd%tend
        case ('reference')
          read(iounit, *) card, dnd%reference
        case default
          call exception_fatal('unknown transient input card: ' // &
            trim(adjustl(card)))
      endselect
    enddo

    close(iounit)
  endsubroutine transient_read

  subroutine transient_summary(dnd)
    use output, only : output_write
    type(DelayedNeutronData), intent(in) :: dnd

    character(2**16) :: line, tmp
    integer(ik) :: d, g

    call output_write('=== TRANSIENT ===')

    write(line, '(a,i0)') 'Number of delayed neutron precursor families: ', dnd%nd
    call output_write(line)

    line = ''
    do d = 1,dnd%nd
      write(tmp, '(es9.2)') dnd%beta(d)
      line = trim(adjustl(line)) // ' ' // trim(adjustl(tmp))
    enddo ! d = 1,dnd%nd
    call output_write('beta = ' // line)

    line = ''
    do d = 1,dnd%nd
      write(tmp, '(es9.2)') dnd%lamd(d)
      line = trim(adjustl(line)) // ' ' // trim(adjustl(tmp))
    enddo ! d = 1,dnd%nd
    call output_write('lambda [1/s] = ' // line)

    write(line, '(a,i0)') 'Number of energy groups in transient data: ', dnd%ng
    call output_write(line)
    
    line = ''
    do g = 1,dnd%ng
      write(tmp, '(es9.2)') dnd%vel(g)
      line = trim(adjustl(line)) // ' ' // trim(adjustl(tmp))
    enddo ! g = dnd%ng
    call output_write('velocity [cm/s] = ' // line)

    write(line, '(a,es9.2)') 'Time-step size [s] = ', dnd%deltat
    call output_write(line)
    write(line, '(a,es9.2)') 'Tend [s] = ', dnd%tend
    call output_write(line)
    
    call output_write('Transient reference case: ' &
      // trim(adjustl(dnd%reference)))

    call output_write('')
  endsubroutine transient_summary

  subroutine transient_init_precursors(nx, mat_map, xslib, dnd, kcrit, flux, prec)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    type(DelayedNeutronData), intent(in) :: dnd
    real(rk), intent(in) :: kcrit
    real(rk), intent(in) :: flux(:,:) ! (nx,ngroup)
    real(rk), intent(out) :: prec(:,:) ! (nx,nd)

    integer(ik) :: i, mthis
    real(rk) :: fsrc

    do i = 1,nx
      mthis = mat_map(i)
      if (xslib%mat(mthis)%is_fiss) then
        fsrc = sum(xslib%mat(mthis)%nusf(:) * flux(i,:))/kcrit
        prec(i,:) = dnd%beta/dnd%lamd * fsrc
      else
        prec(i,:) = 0.0_rk
      endif
    enddo ! i = 1,nx
  endsubroutine transient_init_precursors

  subroutine transient_update_precursors(nx, mat_map, xslib, dnd, kcrit, flux, &
    prec_old, prec)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    type(DelayedNeutronData), intent(in) :: dnd
    real(rk), intent(in) :: kcrit
    real(rk), intent(in) :: flux(:,:) ! (nx,ngroup)
    real(rk), intent(in) :: prec_old(:,:) ! (nx,nd)
    real(rk), intent(out) :: prec(:,:) ! (nx,nd)

    integer(ik) :: i, mthis
    real(rk) :: fsrc

    do i = 1,nx
      mthis = mat_map(i)
      if (xslib%mat(mthis)%is_fiss) then
        fsrc = sum(xslib%mat(mthis)%nusf(:) * flux(i,:))/kcrit
        prec(i,:) = dnd%deltat / (1.0_rk + dnd%deltat * dnd%lamd(:)) &
          * (dnd%beta(:)*fsrc + prec_old(i,:)/dnd%deltat)
      else
        prec(i,:) = 0.0_rk
      endif
    enddo ! i = 1,nx
  endsubroutine transient_update_precursors

  subroutine transient_solve(nx, dx, mat_map, xslib, dnd, &
      boundary_left, boundary_right, phi_tol, max_iter, keff, flux)
    use xs, only : XSLibrary
    use output, only : output_write
    use power, only : power_calculate, power_total
    use diffusion, only : diffusion_build_matrix, diffusion_build_fsource, &
      diffusion_build_upscatter, diffusion_build_downscatter
    use linalg, only : trid
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    type(DelayedNeutronData), intent(in) :: dnd
    character(*), intent(in) :: boundary_left, boundary_right
    real(rk), intent(in) :: phi_tol
    integer(ik), intent(in) :: max_iter
    real(rk), intent(in) :: keff
    real(rk), intent(inout) :: flux(:,:) ! (nx,ngroup)

    real(rk), allocatable :: sub(:,:), dia(:,:), sup(:,:) ! (nx,ngroup)
    real(rk), allocatable :: sub_copy(:), dia_copy(:), sup_copy(:) ! (nx)
    real(rk), allocatable :: fsource(:,:) ! (nx,ngroup)
    real(rk), allocatable :: upsource(:,:) ! (nx,ngroup)
    real(rk), allocatable :: downsource(:) ! (nx)
    real(rk), allocatable :: kinsource(:,:) ! (nx,ngroup)
    real(rk), allocatable :: q(:) ! (nx)

    real(rk), allocatable :: prec(:,:), prec_old(:,:) ! (nx,nd)

    real(rk), allocatable :: flux_old(:,:) ! (nx,ngroup)
    real(rk), allocatable :: power(:) ! (nx)

    real(rk), parameter :: omega = 1.9_rk ! over-relaxation
    real(rk), parameter :: power_init = 1.0_rk ! solving for relative power

    integer(ik) :: step, iter, g
    real(rk) :: tfinal
    real(rk) :: phi_conv

    character(1024) :: line

    allocate(sub(nx,xslib%ngroup))
    allocate(dia(nx,xslib%ngroup))
    allocate(sup(nx,xslib%ngroup))
    allocate(sub_copy(nx))
    allocate(dia_copy(nx))
    allocate(sup_copy(nx))
    sub = 0.0_rk
    dia = 0.0_rk
    sup = 0.0_rk
    sub_copy = 0.0_rk
    dia_copy = 0.0_rk
    sup_copy = 0.0_rk

    allocate(fsource(nx,xslib%ngroup))
    allocate(upsource(nx,xslib%ngroup))
    allocate(downsource(nx))
    allocate(kinsource(nx,xslib%ngroup))
    allocate(q(nx))
    fsource = 0.0_rk
    upsource = 0.0_rk
    downsource = 0.0_rk
    kinsource = 0.0_rk
    q = 0.0_rk

    allocate(flux_old(nx,xslib%ngroup))
    flux_old = 0.0_rk

    allocate(prec(nx,dnd%nd))
    allocate(prec_old(nx,dnd%nd))
    prec = 0.0_rk
    prec_old = 0.0_rk
    call transient_init_precursors(nx, mat_map, xslib, dnd, keff, flux, prec)

    ! compute intital power and normalize
    allocate(power(nx))
    power = 0.0_rk
    call power_calculate(nx, mat_map, xslib, flux, power)
    flux = flux * (power_init / power_total(dx, power))
    ! now, recompute with normalized flux
    call power_calculate(nx, mat_map, xslib, flux, power)

    call output_write('=== TRANSIENT CALCULATION ===')
    call output_write('elapt [s] , deltat [s] , rel. power')
    write(line, '(es13.6, " , ", es13.6, " , ", es13.6)') &
      0.0_rk, dnd%deltat, power_total(dx, power)
    call output_write(line)

    ! TIME LOOP
    step = 0
    do
      ! This works for uniform time steps.
      ! It is more accurate to use a multiplication rather than floating-point
      ! addition which will accumulate round-off.
      step = step + 1
      tfinal = step * dnd%deltat
      
      prec_old = prec

      ! TODO update cross sections

      ! iterative scheme is necesary for one-group at-a-time problem
      do iter = 1,max_iter
        flux_old = flux

        ! TODO build and solve
        ! TODO only the diagonal needs to be modified compared to the
        ! steady-state solution
        call diffusion_build_matrix(&
          nx, dx, mat_map, xslib, boundary_left, boundary_right, sub, dia, sup)

        call diffusion_build_fsource(nx, dx, mat_map, xslib, flux, fsource)
        call diffusion_build_upscatter(nx, dx, mat_map, xslib, flux, upsource)

        do g = 1,xslib%ngroup
          call diffusion_build_downscatter(nx, dx, mat_map, xslib, flux, g, &
            downsource)
          q = fsource(:,g)/keff + upsource(:,g) + kinsource(:,g) + downsource

          ! TODO there is some work to do here
          ! I think that the diagonal will be modified on every iteration, so we
          ! can just trample on it
          sub_copy = sub(:,g)
          dia_copy = dia(:,g)
          sup_copy = sup(:,g)
          call trid(nx, sub_copy, dia_copy, sup_copy, q, flux(:,g))

          ! over-relaxation
          flux(:,g) = flux_old(:,g) + omega * (flux(:,g) - flux_old(:,g))
        enddo ! g = 1,xslib%ngroup

        call transient_update_precursors(nx, mat_map, xslib, dnd, keff, flux, &
          prec_old, prec)

        phi_conv = maxval(abs(flux - flux_old))/maxval(flux)
        if (phi_conv < phi_tol) then
          exit
        endif
      enddo ! iter = 1,max_iter

      call power_calculate(nx, mat_map, xslib, flux, power)
      write(line, '(es13.6, " , ", es13.6, " , ", es13.6)') &
        tfinal, dnd%deltat, power_total(dx, power)
      call output_write(line)

      if (tfinal > dnd%tend) then
        call output_write('Transient Completed!')
        exit
      endif
    enddo ! TIME LOOP

    call output_write('')

    deallocate(sub, dia, sup)
    deallocate(sub_copy, dia_copy, sup_copy)
    deallocate(fsource, upsource, downsource, kinsource, q)
    deallocate(prec, prec_old)
    deallocate(power)
    deallocate(flux_old)
  endsubroutine transient_solve

  subroutine transient_cleanup(dnd)
    type(DelayedNeutronData), intent(out) :: dnd
    if (allocated(dnd%beta)) then
      deallocate(dnd%beta)
    endif
    if (allocated(dnd%lamd)) then
      deallocate(dnd%lamd)
    endif
    if (allocated(dnd%vel)) then
      deallocate(dnd%vel)
    endif
  endsubroutine transient_cleanup

endmodule transient
