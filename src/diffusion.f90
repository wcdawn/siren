module diffusion
use kind, only : rk, ik
implicit none

private

public :: diffusion_power_iteration

contains

  subroutine diffusion_build_matrix(nx, dx, mat_map, xslib, boundary_right, sub, dia, sup, diffusion_coeff)
    use xs, only : XSLibrary
    use exception_handler, only : exception_fatal
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    character(*), intent(in) :: boundary_right
    real(rk), intent(out) :: sub(:,:) ! (nx,ngroup)
    real(rk), intent(out) :: dia(:,:) ! (nx,ngroup)
    real(rk), intent(out) :: sup(:,:) ! (nx,ngroup)
    real(rk), intent(in), optional :: diffusion_coeff(:,:) ! (nx,ngroup)

    integer(ik) :: i, g
    integer(ik) :: mprev, mthis, mnext
    real(rk) :: dprev, dnext

    ! BC at x=0, i=1 (mirror)
    mthis = mat_map(1)
    mnext = mat_map(2)
    do g = 1,xslib%ngroup
      if (present(diffusion_coeff)) then
        dnext = 2 &
          * (diffusion_coeff(1,g) / dx(1) * diffusion_coeff(2,g) / dx(2)) &
          / (diffusion_coeff(1,g) / dx(1) + diffusion_coeff(2,g) / dx(2))
      else
        dnext = 2 &
          * (xslib%mat(mthis)%diffusion(g) / dx(1) * xslib%mat(mnext)%diffusion(g) / dx(2)) &
          / (xslib%mat(mthis)%diffusion(g) / dx(1) + xslib%mat(mnext)%diffusion(g) / dx(2))
      endif
      dia(1,g) = dnext + (xslib%mat(mthis)%sigma_t(g) - xslib%mat(mthis)%scatter(g,g,1)) * dx(1)
      sup(1,g) = -dnext
    enddo

    do g = 1,xslib%ngroup
      do i = 2,nx-1

        mprev = mat_map(i-1)
        mthis = mat_map(i)
        mnext = mat_map(i+1)

        if (present(diffusion_coeff)) then
          dprev = 2 &
            * (diffusion_coeff(i,g) / dx(i) * diffusion_coeff(i-1,g) / dx(i-1)) &
            / (diffusion_coeff(i,g) / dx(i) + diffusion_coeff(i-1,g) / dx(i-1))
          dnext = 2 &
            * (diffusion_coeff(i,g) / dx(i) * diffusion_coeff(i+1,g) / dx(i+1)) &
            / (diffusion_coeff(i,g) / dx(i) + diffusion_coeff(i+1,g) / dx(i+1))
        else
          dprev = 2 &
            * (xslib%mat(mthis)%diffusion(g) / dx(i) * xslib%mat(mprev)%diffusion(g) / dx(i-1)) &
            / (xslib%mat(mthis)%diffusion(g) / dx(i) + xslib%mat(mprev)%diffusion(g) / dx(i-1))
          dnext = 2 &
            * (xslib%mat(mthis)%diffusion(g) / dx(i) * xslib%mat(mnext)%diffusion(g) / dx(i+1)) &
            / (xslib%mat(mthis)%diffusion(g) / dx(i) + xslib%mat(mnext)%diffusion(g) / dx(i+1))
        endif

        sub(i-1,g) = -dprev
        dia(i,g) = dprev + dnext &
          + (xslib%mat(mthis)%sigma_t(g) - xslib%mat(mthis)%scatter(g,g,1)) * dx(i)
        sup(i,g) = -dnext

      enddo
    enddo

    ! BC at x=L, i=N
    mprev = mat_map(nx-1)
    mthis = mat_map(nx)
    select case(boundary_right)
      case ('zero')
        do g = 1,xslib%ngroup
          if (present(diffusion_coeff)) then
            dprev = 2 &
              * (diffusion_coeff(nx,g) / dx(nx) * diffusion_coeff(nx-1,g) / dx(nx-1)) &
              / (diffusion_coeff(nx,g) / dx(nx) + diffusion_coeff(nx-1,g) / dx(nx-1))
          else
            dprev = 2 &
              * (xslib%mat(mthis)%diffusion(g) / dx(nx) * xslib%mat(mprev)%diffusion(g) / dx(nx-1)) &
              / (xslib%mat(mthis)%diffusion(g) / dx(nx) + xslib%mat(mprev)%diffusion(g) / dx(nx-1))
          endif
          sub(nx-1,g) = -dprev
          dia(nx,g) = dprev &
            + (xslib%mat(mthis)%sigma_t(g) - xslib%mat(mthis)%scatter(g,g,1)) * dx(nx)
          if (present(diffusion_coeff)) then
            dia(nx,g) = dia(nx,g) + 2 * diffusion_coeff(nx,g) / dx(nx)
          else
            dia(nx,g) = dia(nx,g) + 2 * xslib%mat(mthis)%diffusion(g) / dx(nx)
          endif
        enddo
      case ('mirror')
        do g = 1,xslib%ngroup
          if (present(diffusion_coeff)) then
            dprev = 2 &
              * (diffusion_coeff(nx,g) / dx(nx) * diffusion_coeff(nx-1,g) / dx(nx-1)) &
              / (diffusion_coeff(nx,g) / dx(nx) + diffusion_coeff(nx-1,g) / dx(nx-1))
          else
            dprev = 2 &
              * (xslib%mat(mthis)%diffusion(g) / dx(nx) * xslib%mat(mprev)%diffusion(g) / dx(nx-1)) &
              / (xslib%mat(mthis)%diffusion(g) / dx(nx) + xslib%mat(mprev)%diffusion(g) / dx(nx-1))
          endif
          sub(nx-1,g) = -dprev
          dia(nx,g) = dprev &
            + (xslib%mat(mthis)%sigma_t(g) - xslib%mat(mthis)%scatter(g,g,1)) * dx(nx)
        enddo
      case default
        call exception_fatal(&
          'unknown boundary_right in diffusion_build_matrix' // trim(adjustl(boundary_right)))
    endselect

  endsubroutine diffusion_build_matrix

  subroutine diffusion_build_fsource(nx, dx, mat_map, xslib, flux, fsource)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(in) :: flux(:,:) ! (nx,ngroup)
    real(rk), intent(out) :: fsource(:,:) ! (nx,ngroup)

    integer(ik) :: i
    integer(ik) :: mthis

    do i = 1,nx
      mthis = mat_map(i)
      if (xslib%mat(mthis)%is_fiss) then
        fsource(i,:) = xslib%mat(mthis)%chi(:) * sum(xslib%mat(mthis)%nusf(:)*flux(i,:)) * dx(i)
      else
        fsource(i,:) = 0d0
      endif
    enddo ! i = 1,nx

  endsubroutine diffusion_build_fsource

  subroutine diffusion_build_upscatter(nx, dx, mat_map, xslib, flux, qup)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(in) :: flux(:,:) ! (nx,ngroup)
    real(rk), intent(out) :: qup(:,:) ! (nx,ngroup)

    integer(ik) :: i, g
    integer(ik) :: mthis

    do i = 1,nx
      mthis = mat_map(i)
      do g = 1,xslib%ngroup
        qup(i,g) = &
          sum(xslib%mat(mthis)%scatter(g+1:xslib%ngroup,g,1) &
          * flux(i,g+1:xslib%ngroup)) * dx(i)
      enddo
    enddo ! i = 1,nx

  endsubroutine diffusion_build_upscatter

  subroutine diffusion_build_downscatter(nx, dx, mat_map, xslib, flux, g, qdown)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(in) :: flux(:,:) ! (nx,ngroup)
    integer(ik), intent(in) :: g ! group index
    real(rk), intent(out) :: qdown(:) ! (nx,ngroup)

    integer(ik) :: i
    integer(ik) :: mthis

    do i = 1,nx
      mthis = mat_map(i)
      qdown(i) = sum(xslib%mat(mthis)%scatter(1:g-1,g,1) * flux(i,1:g-1)) * dx(i)
    enddo ! i = 1,nx

  endsubroutine diffusion_build_downscatter

  real(rk) function diffusion_fission_summation(nx, dx, mat_map, xslib, flux)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(in) :: flux(:,:) ! (nx,ngroup)

    integer(ik) :: i, mthis
    real(rk) :: xsum

    xsum = 0d0
    do i = 1,nx
      mthis = mat_map(i)
      if (xslib%mat(mthis)%is_fiss) then
        xsum = xsum + sum(xslib%mat(mthis)%nusf(:) * flux(i,:)) * dx(i)
      endif
    enddo ! i = 1,nx
    diffusion_fission_summation = xsum
  endfunction diffusion_fission_summation

  subroutine diffusion_power_iteration( &
    nx, dx, mat_map, xslib, boundary_right, k_tol, phi_tol, max_iter, keff, flux, transportxs)
    use xs, only : XSLibrary
    use linalg, only : trid
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
    real(rk), intent(in), optional :: transportxs(:,:) ! (nx,ngroup)

    real(rk), allocatable :: sub(:,:), dia(:,:), sup(:,:)
    real(rk), allocatable :: sub_copy(:), dia_copy(:), sup_copy(:)
    real(rk), allocatable :: fsource(:,:), upsource(:,:), downsource(:)
    real(rk), allocatable :: q(:)

    integer(ik) :: i, g
    integer(ik) :: iter
    real(rk) :: k_old, fsum, fsum_old
    real(rk), allocatable :: flux_old(:,:)
    real(rk) :: delta_k, delta_phi

    real(rk), allocatable :: diffusion_coeff(:,:) ! (nx,ngroup)

    character(1024) :: line

    logical, parameter :: matrix_dump = .false.

    allocate(sub(nx-1,xslib%ngroup), dia(nx,xslib%ngroup), sup(nx-1,xslib%ngroup))
    allocate(sub_copy(nx-1), dia_copy(nx), sup_copy(nx-1))

    if (present(transportxs)) then
      allocate(diffusion_coeff(nx,xslib%ngroup))
      call diffusion_populate_coeff(nx, mat_map, xslib, transportxs, diffusion_coeff)
    endif

    call timer_start('diffusion_build_matrix')
    if (allocated(diffusion_coeff)) then
      call diffusion_build_matrix(nx, dx, mat_map, xslib, boundary_right, sub, dia, sup, diffusion_coeff)
    else
      call diffusion_build_matrix(nx, dx, mat_map, xslib, boundary_right, sub, dia, sup)
      if (matrix_dump) then
        write(730, '(es13.6,",",es13.6,",",es13.6)') 0d0 , dia(1,1), sup(1,1)
        do i = 2,nx-1
          write(730, '(es13.6,",",es13.6,",",es13.6)') sub(i-1,1), dia(i,1), sup(i,1)
        enddo
        write(730, '(es13.6,",",es13.6,",",es13.6)') sub(nx-1,1), dia(nx,1), 0.0
      endif
    endif
    call timer_stop('diffusion_build_matrix')

    allocate(fsource(nx,xslib%ngroup))
    allocate(upsource(nx,xslib%ngroup))
    allocate(downsource(nx))
    allocate(q(nx))

    allocate(flux_old(nx,xslib%ngroup))

    ! initialization
    flux = 1d0
    keff = 1d0
    fsum = 1d0

    call output_write("=== DIFFUSION POWER ITERATION ===")

    do iter = 1,max_iter
      k_old = keff
      flux_old = flux
      fsum_old = fsum

      call timer_start('diffusion_source')
      call diffusion_build_fsource(nx, dx, mat_map, xslib, flux, fsource)
      call diffusion_build_upscatter(nx, dx, mat_map, xslib, flux, upsource)
      call timer_stop('diffusion_source')

      do g = 1,xslib%ngroup

        call timer_start('diffusion_source')
        call diffusion_build_downscatter(nx, dx, mat_map, xslib, flux, g, downsource)
        q = fsource(:,g)/keff + upsource(:,g) + downsource
        call timer_stop('diffusion_source')

        call timer_start('diffusion_tridiagonal')
        ! SOLVE
        ! need to store copies, trid uses them as scratch space
        sub_copy = sub(:,g)
        dia_copy = dia(:,g)
        sup_copy = sup(:,g)
        call trid(nx, sup_copy, dia_copy, sup_copy, q, flux(:,g))
        call timer_stop('diffusion_tridiagonal')

      enddo

      fsum = diffusion_fission_summation(nx, dx, mat_map, xslib, flux)
      if (iter > 1) keff = keff * fsum / fsum_old
      delta_k = abs(keff - k_old)
      delta_phi = maxval(abs(flux - flux_old))/maxval(flux)

      write(line, '(a,i4,a,es8.1,a,es8.1,a,f8.6)') &
        'it=', iter, ' dk=', delta_k, ' dphi=', delta_phi, ' keff=', keff
      call output_write(line)

      if ((delta_k < k_tol) .and. (delta_phi < phi_tol)) then
        call output_write('CONVERGENCE!')
        call output_write('')
        exit
      endif
    enddo

    if (iter > max_iter) then
      call exception_warning('failed to converge')
    endif

    deallocate(flux_old)
    deallocate(sub, dia, sup)
    deallocate(sub_copy, dia_copy, sup_copy)
    deallocate(fsource, upsource, downsource)
    deallocate(q)
    if (allocated(diffusion_coeff)) then
      deallocate(diffusion_coeff)
    endif
  endsubroutine diffusion_power_iteration

  subroutine diffusion_populate_coeff(nx, mat_map, xslib, transportxs, diffusion_coeff)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(in) :: transportxs(:,:) ! (nx,ngroup)
    real(rk), intent(out) :: diffusion_coeff(:,:) ! (nx,ngroup)

    integer(ik) :: i, g, mthis
    real(rk) :: xsmax, xsmin

    logical, parameter :: clamp = .true.

    ! copy
    diffusion_coeff(1:nx,1:xslib%ngroup) = transportxs(1:nx,1:xslib%ngroup)
    
    if (clamp) then
      do i = 1,nx
        mthis = mat_map(i)
        xsmax = maxval(xslib%mat(mthis)%sigma_t)
        xsmin = minval(xslib%mat(mthis)%sigma_t)
        do g = 1,xslib%ngroup
          diffusion_coeff(i,g) = min(diffusion_coeff(i,g), xsmax)
          diffusion_coeff(i,g) = max(diffusion_coeff(i,g), xsmin)
        enddo ! g = 1,xslib%ngroup
      enddo ! i = 1,nx
    endif

    diffusion_coeff = 1.0_rk / (3.0_rk * diffusion_coeff)
  endsubroutine diffusion_populate_coeff

endmodule diffusion
