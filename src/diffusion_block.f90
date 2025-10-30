module diffusion_block
use kind, only : rk, ik
implicit none

private
public :: diffusion_block_power_iteration

contains

  subroutine diffusion_block_build_matrix(nx, dx, mat_map, xslib, boundary_right, sub, dia, sup)
    use xs, only : XSLibrary
    use exception_handler, only : exception_fatal
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    character(*), intent(in) :: boundary_right
    real(rk), intent(out) :: sub(:,:,:) ! (ngroup,ngroup,nx-1)
    real(rk), intent(out) :: dia(:,:,:) ! (ngroup,ngroup,nx)
    real(rk), intent(out) :: sup(:,:,:) ! (ngroup,ngroup,nx-1)

    integer(ik) :: i, g
    integer(ik) :: mprev, mthis, mnext
    real(rk) :: dprev, dnext

    sub(1:xslib%ngroup,1:xslib%ngroup,1:nx-1) = 0d0
    dia(1:xslib%ngroup,1:xslib%ngroup,1:nx) = 0d0
    sup(1:xslib%ngroup,1:xslib%ngroup,1:nx-1) = 0d0

    ! BC at x=0, i=1 (mirror)
    mthis = mat_map(1)
    mnext = mat_map(2)
    do g = 1,xslib%ngroup
      dnext = 2 &
        * (xslib%mat(mthis)%diffusion(g) / dx(1) * xslib%mat(mnext)%diffusion(g) / dx(2)) &
        / (xslib%mat(mthis)%diffusion(g) / dx(1) + xslib%mat(mnext)%diffusion(g) / dx(2))
      dia(g,g,1) = dnext + xslib%mat(mthis)%sigma_t(g) * dx(1)
      sup(g,g,1) = -dnext
    enddo ! g = 1,xslib%ngroup

    do i = 2,nx-1

      mprev = mat_map(i-1)
      mthis = mat_map(i)
      mnext = mat_map(i+1)

      do g = 1,xslib%ngroup
        dprev = 2 &
          * (xslib%mat(mthis)%diffusion(g) / dx(i) * xslib%mat(mprev)%diffusion(g) / dx(i-1)) &
          / (xslib%mat(mthis)%diffusion(g) / dx(i) + xslib%mat(mprev)%diffusion(g) / dx(i-1))
        dnext = 2 &
          * (xslib%mat(mthis)%diffusion(g) / dx(i) * xslib%mat(mnext)%diffusion(g) / dx(i+1)) &
          / (xslib%mat(mthis)%diffusion(g) / dx(i) + xslib%mat(mnext)%diffusion(g) / dx(i+1))

        sub(g,g,i-1) = -dprev
        dia(g,g,i) = dprev + dnext + xslib%mat(mthis)%sigma_t(g) * dx(i)
        sup(g,g,i) = -dnext
      enddo ! g = 1,xslib%ngroup
    enddo ! i = 2,nx-1

    ! BC at x=L, i=N
    mprev = mat_map(nx-1)
    mthis = mat_map(nx)
    select case (boundary_right)
      case ('zero')
        do g = 1,xslib%ngroup
          dprev = 2 &
            * (xslib%mat(mthis)%diffusion(g) / dx(nx) * xslib%mat(mprev)%diffusion(g) / dx(nx-1)) &
            / (xslib%mat(mthis)%diffusion(g) / dx(nx) + xslib%mat(mprev)%diffusion(g) / dx(nx-1))
          sub(g,g,nx-1) = -dprev
          dia(g,g,nx) = dprev + xslib%mat(mthis)%sigma_t(g) * dx(nx) &
            + 2 * xslib%mat(mthis)%diffusion(g)/dx(nx)
        enddo ! g = 1,xslib%ngroup
      case ('mirror')
        do g = 1,xslib%ngroup
          dprev = 2 &
            * (xslib%mat(mthis)%diffusion(g) / dx(nx) * xslib%mat(mprev)%diffusion(g) / dx(nx-1)) &
            / (xslib%mat(mthis)%diffusion(g) / dx(nx) + xslib%mat(mprev)%diffusion(g) / dx(nx-1))
          sub(g,g,nx-1) = -dprev
          dia(g,g,nx) = dprev + xslib%mat(mthis)%sigma_t(g) * dx(nx)
        enddo ! g = 1,xslib%ngroup
      case default
        call exception_fatal('unknown boundary_right: ' // trim(adjustl(boundary_right)))
    endselect

    ! remove scattering
    do i = 1,nx
      mthis = mat_map(i)
      dia(:,:,i) = dia(:,:,i) - transpose(xslib%mat(mthis)%scatter(:,:,1)) * dx(i)
    enddo
  endsubroutine diffusion_block_build_matrix

  subroutine diffusion_block_build_fsource(nx, dx, mat_map, xslib, flux, fsource)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(in) :: flux(:,:) ! (ngroup,nx) NOTE this is "block" order
    real(rk), intent(out) :: fsource(:,:) ! (ngroup,nx)

    integer(ik) :: i, mthis

    do i = 1,nx
      mthis = mat_map(i)
      if (xslib%mat(mthis)%is_fiss) then
        fsource(:,i) = xslib%mat(mthis)%chi * sum(xslib%mat(mthis)%nusf*flux(:,i)) * dx(i)
      else
        fsource(:,i) = 0d0
      endif
    enddo ! i = 1,nx
  endsubroutine diffusion_block_build_fsource

  real(rk) function diffusion_block_fission_summation(nx, dx, mat_map, xslib, flux)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(in) :: flux(:,:) ! (ngroup,nx) NOTE this is in "block" order

    integer(ik) :: i, mthis
    real(rk) :: xsum

    xsum = 0d0
    do i = 1,nx
      mthis = mat_map(i)
      if (xslib%mat(mthis)%is_fiss) then
        xsum = xsum + sum(xslib%mat(mthis)%nusf * flux(:,i))*dx(i)
      endif
    enddo ! i = 1,nx
    diffusion_block_fission_summation = xsum
  endfunction diffusion_block_fission_summation

  subroutine diffusion_block_power_iteration( &
    nx, dx, mat_map, xslib, boundary_right, k_tol, phi_tol, max_iter, keff, flux)
    use xs, only : XSLibrary
    use linalg, only : trid_block
    use output, only : output_write
    use timer, only : timer_start, timer_stop
    use exception_handler, only : exception_warning
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    character(*), intent(in) :: boundary_right
    real(rk), intent(in) :: k_tol, phi_tol
    integer(ik), intent(in) :: max_iter
    real(rk), intent(out) :: keff
    real(rk), intent(out) :: flux(:,:) ! (nx,ngroup)

    integer(ik) :: iter
    integer(ik) :: i, g
    real(rk) :: fsum, fsum_old, k_old
    real(rk) :: delta_k, delta_phi

    ! (ngroup,ngroup,nx-1) , (nxgroup,ngroup,nx) , (ngroup,ngroup,nx-1)
    real(rk), allocatable :: sub(:,:,:), dia(:,:,:), sup(:,:,:)
    real(rk), allocatable :: sub_copy(:,:,:), dia_copy(:,:,:), sup_copy(:,:,:)

    real(rk), allocatable :: fsource(:,:) ! (ngroup,nx)

    ! NOTE: I'm storing copies here
    ! they're in transposed indices to make my block solves easier
    real(rk), allocatable :: flux_block(:,:), flux_old(:,:) ! (ngroup,nx)

    character(1024) :: line

    allocate(&
      sub(xslib%ngroup,xslib%ngroup,nx-1),&
      dia(xslib%ngroup,xslib%ngroup,nx),&
      sup(xslib%ngroup,xslib%ngroup,nx-1)&
    )
    allocate(sub_copy(xslib%ngroup,xslib%ngroup,nx-1), dia_copy(xslib%ngroup,xslib%ngroup,nx), &
      sup_copy(xslib%ngroup,xslib%ngroup,nx-1))

    call timer_start('diffusion_build_matrix')
    call diffusion_block_build_matrix(nx, dx, mat_map, xslib, boundary_right, sub, dia, sup)
    call timer_stop('diffusion_build_matrix')

    allocate(fsource(xslib%ngroup,nx))
    allocate(flux_block(xslib%ngroup,nx))
    allocate(flux_old(xslib%ngroup,nx))

    ! initialize
    flux_block = 1d0
    keff = 1d0
    fsum = 1d0

    call output_write("=== DIFFUSION BLOCK POWER ITERATION ===")

    do iter = 1,max_iter
      k_old = keff
      flux_old = flux_block
      fsum_old = fsum

      call timer_start('diffusion_source')
      call diffusion_block_build_fsource(nx, dx, mat_map, xslib, flux_block, fsource)
      call timer_stop('diffusion_source')

      fsource = fsource / keff

      call timer_start('diffusion_block_tridiagonal')
      sub_copy = sub
      dia_copy = dia
      sup_copy = sup
      call trid_block(nx, xslib%ngroup, sub_copy, dia_copy, sup_copy, fsource, flux_block)
      call timer_stop('diffusion_block_tridiagonal')

      fsum = diffusion_block_fission_summation(nx, dx, mat_map, xslib, flux_block)
      if (iter > 1) then
        keff = keff * fsum / fsum_old
      endif
      delta_k = abs(keff - k_old)
      delta_phi = maxval(abs(flux_block - flux_old))/maxval(flux_block)

      write(line, '(a,i4,a,es8.1,a,es8.1,a,f8.6)') &
        'it=', iter, ' dk=', delta_k, ' dphi=', delta_phi, ' keff=', keff
      call output_write(line)

      if ((delta_k < k_tol) .and. (delta_phi < phi_tol)) then
        call output_write('CONVERGENCE!')
        call output_write('')
        exit
      endif
    enddo ! iter = 1,max_iter

    if (iter > max_iter) then
      call exception_warning('failed to converge')
    endif

    ! copy before exit
    do i = 1,nx
      do g = 1,xslib%ngroup
        flux(i,g) = flux_block(g,i)
      enddo ! g = 1,xslib%ngroup
    enddo ! i = 1,nx

    deallocate(sub, dia, sup)
    deallocate(sub_copy, dia_copy, sup_copy)
    deallocate(fsource)
    deallocate(flux_block, flux_old)
  endsubroutine diffusion_block_power_iteration

endmodule diffusion_block
