module diffusion
use kind
implicit none

private

public :: diffusion_power_iteration

logical, parameter :: tworeg_benchmark = .false.

contains

  subroutine diffusion_build_matrix(nx, hx, mat_map, xslib, sub, dia, sup)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: hx
    integer(ik), intent(in) :: mat_map(:)
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(out) :: sub(:,:) ! (nx,ngroup)
    real(rk), intent(out) :: dia(:,:) ! (nx,ngroup)
    real(rk), intent(out) :: sup(:,:) ! (nx,ngroup)

    integer(ik) :: i, g
    integer(ik) :: mprev, mthis, mnext
    real(rk) :: dprev, dnext

    ! BC at x=0, i=1 (mirror)
    mthis = mat_map(1)
    mnext = mat_map(2)
    do g = 1,xslib%ngroup
      dnext = 2 * xslib%mat(mthis)%diffusion(g) * xslib%mat(mnext)%diffusion(g) / &
        (xslib%mat(mthis)%diffusion(g) + xslib%mat(mnext)%diffusion(g))
      dia(1,g) = dnext + (xslib%mat(mthis)%sigma_t(g) - xslib%mat(mthis)%scatter(g,g,1)) * hx**2
      sup(1,g) = -dnext
    enddo

    do g = 1,xslib%ngroup
      do i = 2,nx-1

        mprev = mat_map(i-1)
        mthis = mat_map(i)
        mnext = mat_map(i+1)

        dprev = 2 * xslib%mat(mthis)%diffusion(g) * xslib%mat(mprev)%diffusion(g) / &
          (xslib%mat(mthis)%diffusion(g) + xslib%mat(mprev)%diffusion(g))
        dnext = 2 * xslib%mat(mthis)%diffusion(g) * xslib%mat(mnext)%diffusion(g) / &
          (xslib%mat(mthis)%diffusion(g) + xslib%mat(mnext)%diffusion(g))

        sub(i-1,g) = -dprev
        dia(i,g) = dprev + dnext + (xslib%mat(mthis)%sigma_t(g) - xslib%mat(mthis)%scatter(g,g,1)) * hx**2
        sup(i,g) = -dnext

      enddo
    enddo

    ! BC at x=L, i=N
    if (tworeg_benchmark) then
      mprev = mat_map(nx-1)
      mthis = mat_map(nx)
      do g = 1,xslib%ngroup
        dprev = 2 * xslib%mat(mthis)%diffusion(g) * xslib%mat(mprev)%diffusion(g) / &
          (xslib%mat(mthis)%diffusion(g) + xslib%mat(mprev)%diffusion(g))
        sub(nx-1,g) = -dprev
        dia(nx,g) = dprev + (xslib%mat(mthis)%sigma_t(g) - xslib%mat(mthis)%scatter(g,g,1)) * hx**2 &
          + 2 * xslib%mat(mthis)%diffusion(g)
      enddo
    else
      mprev = mat_map(nx-1)
      mthis = mat_map(nx)
      do g = 1,xslib%ngroup
        dprev = 2 * xslib%mat(mthis)%diffusion(g) * xslib%mat(mprev)%diffusion(g) / &
          (xslib%mat(mthis)%diffusion(g) + xslib%mat(mprev)%diffusion(g))
        sub(nx-1,g) = -dprev
        dia(nx,g) = dprev + (xslib%mat(mthis)%sigma_t(g) - xslib%mat(mthis)%scatter(g,g,1)) * hx**2
      enddo
    endif

  endsubroutine diffusion_build_matrix

  subroutine diffusion_build_fsource(nx, hx, mat_map, xslib, flux, fsource)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: hx
    integer(ik), intent(in) :: mat_map(:)
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(in) :: flux(:,:) ! (nx,ngroup)
    real(rk), intent(out) :: fsource(:,:) ! (nx,ngroup)

    integer(ik) :: i
    integer(ik) :: mthis

    do i = 1,nx
      mthis = mat_map(i)
      if (xslib%mat(mthis)%is_fiss) then
        fsource(i,:) = xslib%mat(mthis)%chi(:) * sum(xslib%mat(mthis)%nusf(:)*flux(i,:)) * hx**2
      else
        fsource(i,:) = 0d0
      endif
    enddo ! i = 1,nx

  endsubroutine diffusion_build_fsource

  subroutine diffusion_build_upscatter(nx, hx, mat_map, xslib, flux, qup)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: hx
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(in) :: flux(:,:) ! (nx,ngroup)
    real(rk), intent(out) :: qup(:,:) ! (nx,ngroup)

    integer(ik) :: i, g
    integer(ik) :: mthis

    do i = 1,nx
      mthis = mat_map(i)
      do g = 1,xslib%ngroup
        qup(i,g) = sum(xslib%mat(mthis)%scatter(g+1:xslib%ngroup,g,1) * flux(i,g+1:xslib%ngroup))
      enddo
    enddo ! i = 1,nx
    qup = qup * hx**2

  endsubroutine diffusion_build_upscatter

  subroutine diffusion_build_downscatter(nx, hx, mat_map, xslib, flux, g, qdown)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: hx
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(in) :: flux(:,:) ! (nx,ngroup)
    integer(ik), intent(in) :: g ! group index
    real(rk), intent(out) :: qdown(:) ! (nx,ngroup)

    integer(ik) :: i
    integer(ik) :: mthis

    do i = 1,nx
      mthis = mat_map(i)
      qdown(i) = sum(xslib%mat(mthis)%scatter(1:g-1,g,1) * flux(i,1:g-1))
    enddo ! i = 1,nx
    qdown = qdown * hx**2

  endsubroutine diffusion_build_downscatter

  real(rk) function diffusion_fission_summation(nx, mat_map, xslib, flux)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(in) :: flux(:,:) ! (nx,ngroup)

    integer(ik) :: i, mthis
    real(rk) :: xsum

    xsum = 0d0
    do i = 1,nx
      mthis = mat_map(i)
      if (xslib%mat(mthis)%is_fiss) then
        xsum = xsum + sum(xslib%mat(mthis)%nusf(:) * flux(i,:))
      endif
    enddo ! i = 1,nx
    diffusion_fission_summation = xsum
    
  endfunction diffusion_fission_summation

  subroutine diffusion_power_iteration(nx, hx, mat_map, xslib, k_tol, phi_tol, max_iter, keff, flux)
    use xs, only : XSLibrary
    use linalg, only : trid
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: hx
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(in) :: k_tol, phi_tol
    integer(ik), intent(in) :: max_iter
    real(rk), intent(out) :: keff
    real(rk), intent(out) :: flux(:,:) ! (nx,ngroup)

    real(rk), allocatable :: sub(:,:), dia(:,:), sup(:,:)
    real(rk), allocatable :: sub_copy(:,:), dia_copy(:,:), sup_copy(:,:)
    real(rk), allocatable :: fsource(:,:), upsource(:,:), downsource(:)
    real(rk), allocatable :: q(:)

    integer(ik) :: g
    integer(ik) :: iter
    real(rk) :: k_old, fsum, fsum_old
    real(rk), allocatable :: flux_old(:,:)
    real(rk) :: delta_k, delta_phi

    allocate(sub(nx-1,xslib%ngroup), dia(nx,xslib%ngroup), sup(nx-1,xslib%ngroup))
    allocate(sub_copy(nx-1,xslib%ngroup), dia_copy(nx,xslib%ngroup), sup_copy(nx-1,xslib%ngroup))
    call diffusion_build_matrix(nx, hx, mat_map, xslib, sub, dia, sup)

    allocate(fsource(nx,xslib%ngroup))
    allocate(upsource(nx,xslib%ngroup))
    allocate(downsource(nx))
    allocate(q(nx))

    allocate(flux_old(nx,xslib%ngroup))

    ! initialization
    flux = 1d0
    keff = 1d0
    fsum = 1d0

    write(*,*) "=== DIFFUSION POWER ITERATION ==="

    do iter = 1,max_iter
      k_old = keff
      flux_old = flux
      fsum_old = fsum

      call diffusion_build_fsource(nx, hx, mat_map, xslib, flux, fsource)
      call diffusion_build_upscatter(nx, hx, mat_map, xslib, flux, upsource)

      do g = 1,xslib%ngroup
        call diffusion_build_downscatter(nx, hx, mat_map, xslib, flux, g, downsource)
        q = fsource(:,g)/keff + upsource(:,g) + downsource
        ! SOLVE
        ! need to store copies, trid uses them as scratch space
        sub_copy(:,g) = sub(:,g)
        dia_copy(:,g) = dia(:,g)
        sup_copy(:,g) = sup(:,g)
        call trid(nx, sup_copy(:,g), dia_copy(:,g), sup_copy(:,g), q, flux(:,g))
      enddo

      fsum = diffusion_fission_summation(nx, mat_map, xslib, flux)
      if (iter > 1) keff = keff * fsum / fsum_old
      delta_k = abs(keff - k_old)
      delta_phi = maxval(abs(flux - flux_old))/maxval(flux)

      write(*,'(a,i4,a,es8.1,a,es8.1,a,f8.6)') &
        'it=', iter, ' dx=', delta_k, ' dphi=', delta_phi, ' keff=', keff

      if ((delta_k < k_tol) .and. (delta_phi < phi_tol)) then
        write(*,*) 'CONVERGENCE!'
        write(*,*)
        exit
      endif
    enddo

    if (iter > max_iter) then
      write(*,*) 'WARNING: failed to converge'
    endif

    deallocate(flux_old)
    deallocate(sub, dia, sup)
    deallocate(sub_copy, dia_copy, sup_copy)
    deallocate(fsource, upsource, downsource)
    deallocate(q)

  endsubroutine diffusion_power_iteration

endmodule diffusion
