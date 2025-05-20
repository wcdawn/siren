module transport
use kind
implicit none

private

public :: transport_power_iteration, sigma_tr

integer(ik) :: neven
real(rk), allocatable :: sigma_tr(:,:,:) ! (nx, ngroup, nmoment)

! 0d0 < damping < 1d0  -- damped
! damping == 1d0 -- undamped
! 1d0 < damping < 2d0 -- accelerated

real(rk), parameter :: damping_transport = 0.5d0
real(rk), parameter :: damping_odd = 1d0
integer(ik), parameter :: dead_iter = 5

contains

  subroutine transport_build_matrix(nx, hx, mat_map, xslib, neven, sub, dia, sup)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: hx
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    integer, intent(in) :: neven
    real(rk), intent(out) :: sub(:,:,:), dia(:,:,:), sup(:,:,:) ! (nx, ngroup, neven)

    integer(ik) :: i, g, n
    integer(ik) :: mprev, mthis, mnext
    integer(ik) :: idxn
    real(rk) :: xn, xmul_next, xmul_prev
    real(rk) :: cprev, cthis, cnext
    real(rk) :: dprev, dnext

    ! BC at x=0, i=1
    mthis = mat_map(1)
    mnext = mat_map(2)
    do n = 1,neven
      idxn = 2*(n-1)
      xn = real(idxn, rk)
      xmul_next = (xn+1d0)**2 / ((2d0*xn+1d0)*(2d0*xn+3d0))
      xmul_prev = xn**2 / (4d0*xn**2 - 1d0)
      do g = 1,xslib%ngroup
        cthis = xmul_next / sigma_tr(1,g,idxn+1+1)
        cnext = xmul_next / sigma_tr(2,g,idxn+1+1)
        if (idxn > 1) then
          cthis = cthis + xmul_prev / sigma_tr(1,g,idxn+1-1)
          cnext = cnext + xmul_prev / sigma_tr(2,g,idxn+1-1)
        endif
        dnext = 2 * cthis*cnext / (cthis + cnext)
        dia(1,g,n) = + dnext + xslib%mat(mthis)%sigma_t(g) * hx**2
        if (idxn+1 < xslib%nmoment) then
          dia(1,g,n) = dia(1,g,n) - xslib%mat(mthis)%scatter(g,g,idxn+1) * hx**2
        endif
        sup(1,g,n) = -dnext
      enddo ! = g = 1,xslib%ngroup
    enddo ! n = 1,neven

    do n = 1,neven
      idxn = 2*(n-1)
      xn = real(idxn, rk)
      xmul_next = (xn+1d0)**2 / ((2d0*xn+1d0)*(2d0*xn+3d0))
      xmul_prev = xn**2 / (4d0*xn**2 - 1d0)
      do g = 1,xslib%ngroup
        do i = 2,nx-1

          mprev = mat_map(i-1)
          mthis = mat_map(i)
          mnext = mat_map(i+1)

          cprev = xmul_next / sigma_tr(i-1,g,idxn+1+1)
          cthis = xmul_next / sigma_tr(i  ,g,idxn+1+1)
          cnext = xmul_next / sigma_tr(i+1,g,idxn+1+1)
          if (idxn > 0) then
            cprev = cprev + xmul_prev / sigma_tr(i-1,g,idxn+1-1)
            cthis = cthis + xmul_prev / sigma_tr(i  ,g,idxn+1-1)
            cnext = cnext + xmul_prev / sigma_tr(i+1,g,idxn+1-1)
          endif

          dprev = 2 * cthis*cprev / (cthis + cprev)
          dnext = 2 * cthis*cnext / (cthis + cnext)

          sub(i-1,g,n) = -dprev
          dia(i,g,n) = dprev + dnext + xslib%mat(mthis)%sigma_t(g) * hx**2
          if (idxn+1 < xslib%nmoment) then
            dia(i,g,n) = dia(i,g,n) - xslib%mat(mthis)%scatter(g,g,idxn+1) * hx**2
          endif
          sup(i,g,n) = -dnext

        enddo ! i = 2,nx-1
      enddo ! g = 1,xslib%ngroup
    enddo ! n = 1,neven

    ! BC at x=L, i=N
    mprev = mat_map(nx-1)
    mthis = mat_map(nx)
    do n = 1,neven
      idxn = 2*(n-1)
      xn = real(idxn, rk)
      xmul_next = (xn+1d0)**2 / ((2d0*xn+1d0)*(2d0*xn+3d0))
      xmul_prev = xn**2 / (4d0*xn**2 - 1d0)
      do g = 1,xslib%ngroup
        cprev = xmul_next / sigma_tr(nx-1,g,idxn+1+1)
        cthis = xmul_next / sigma_tr(nx  ,g,idxn+1+1)
        if (idxn > 0) then
          cprev = cprev + xmul_prev / sigma_tr(nx-1,g,idxn+1-1)
          cthis = cthis + xmul_prev / sigma_tr(nx  ,g,idxn+1-1)
        endif
        dprev = 2 * cthis*cprev / (cthis + cprev)
        sub(nx-1,g,n) = -dprev
        dia(nx,g,n) = dprev + xslib%mat(mthis)%sigma_t(g) * hx**2
        if (idxn+1 < xslib%nmoment) then
          dia(nx,g,n) = dia(nx,g,n) - xslib%mat(mthis)%scatter(g,g,idxn+1) * hx**2
        endif
      enddo ! g = 1,xslib%ngroup
    enddo ! n = 1,neven

  endsubroutine transport_build_matrix

  subroutine transport_build_transportxs(nx, mat_map, xslib, pnorder, phi, sigma_tr)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    integer(ik), intent(in) :: mat_map(:)
    type(XSLibrary), intent(in) :: xslib
    integer(ik), intent(in) :: pnorder
    real(rk), intent(in) :: phi(:,:,:) ! (nx, ngroup, pnorder)
    real(rk), intent(out) :: sigma_tr(:,:,:) ! (nx, ngroup, pnorder)

    integer(ik) :: n, g, i
    integer(ik) :: mthis

    real(rk) :: xnew

    do n = 1,pnorder+1
      do g = 1,xslib%ngroup
        do i = 1,nx
          mthis = mat_map(i)

          if (n <= xslib%nmoment) then
            xnew = xslib%mat(mthis)%sigma_t(g) - sum(xslib%mat(mthis)%scatter(:,g,n)*phi(i,:,n))
          else
            xnew = xslib%mat(mthis)%sigma_t(g)
          endif

          sigma_tr(i,g,n) = sigma_tr(i,g,n) + damping_transport * (xnew - sigma_tr(i,g,n))
        enddo ! i = 1,nx
      enddo ! g = 1,ngroup
    enddo ! n = 1,pnorder+1
  endsubroutine transport_build_transportxs

  subroutine transport_odd_update(nx, hx, ng, pnorder, sigma_tr, phi)
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: hx
    integer(ik), intent(in) :: ng
    integer(ik), intent(in) :: pnorder
    real(rk), intent(in) :: sigma_tr(:,:,:) ! (nx, ngroup, pnorder)
    real(rk), intent(inout) :: phi(:,:,:) ! (nx, ngroup, pnorder)

    integer(ik) :: n, g, i
    integer(ik) :: idxn
    real(rk) :: xn, xmul_prev, xmul_next
    real(rk) :: xnew

    do n = 2,pnorder+1,2
      idxn = n-1
      xn = real(idxn, rk)
      xmul_prev = xn/(2d0*xn+1d0)
      xmul_next = (xn+1d0)/(2d0*n+1)
      do g = 1,ng
        do i = 2,nx-1
          if (n < pnorder+1) then
            xnew = -0.5d0/(sigma_tr(i,g,n)*hx) &
              * (xmul_prev * (phi(i+1,g,idxn+1-1) - 2d0*phi(i,g,idxn+1-1) + phi(i-1,g,idxn+1-1)) &
              + xmul_next * (phi(i+1,g,idxn+1+1) - 2d0*phi(i,g,idxn+1+1) + phi(i-1,g,idxn+1+1)))
          else
            xnew = -0.5d0/(sigma_tr(i,g,n)*hx) &
              * xmul_prev * (phi(i+1,g,idxn+1-1) - 2d0*phi(i,g,idxn+1-1) + phi(i-1,g,idxn+1-1))
          endif
          phi(i, g, idxn+1) = phi(i, g, idxn+1) + damping_odd * (xnew - phi(i, g, idxn+1))
        enddo ! i = 2,nx-1
        ! boundary conditions
        phi(1,g,idxn+1) = phi(2,g,idxn+1)/3d0
        phi(nx,g,idxn+1) = phi(nx-1,g,idxn+1)/3d0
      enddo ! g = 1,ng
    enddo ! n = 2,pnorder+1,2
  endsubroutine transport_odd_update

  subroutine transport_build_upscatter(nx, hx, mat_map, xslib, phi, idxn, qup)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: hx
    integer(ik), intent(in) :: mat_map(:)
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(in) :: phi(:,:,:) ! (nx, ngroup, pnorder)
    integer(ik), intent(in) :: idxn
    real(rk), intent(out) :: qup(:,:) ! (nx, ngroup) 

    integer(ik) :: i, g
    integer(ik) :: mthis

    if (idxn+1 > xslib%nmoment) then
      qup = 0d0
    else
      do i = 1,nx
        mthis = mat_map(i)
        do g = 1,xslib%ngroup
          qup(i,g) = sum(xslib%mat(mthis)%scatter(g+1:xslib%ngroup,g,idxn+1) * phi(i,g+1:xslib%ngroup,idxn+1))
        enddo
      enddo ! i = 1,nx
      qup = qup * hx**2
    endif
  endsubroutine transport_build_upscatter

  subroutine transport_build_downscatter(nx, hx, mat_map, xslib, phi, idxn, g, qdown)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: hx
    integer(ik), intent(in) :: mat_map(:)
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(in) :: phi(:,:,:) ! (nx, ngroup, pnorder)
    integer(ik), intent(in) :: idxn
    integer(ik), intent(in) :: g
    real(rk), intent(out) :: qdown(:) ! (nx) 

    integer(ik) :: i
    integer(ik) :: mthis

    if (idxn+1 > xslib%nmoment) then
      qdown = 0d0
    else
      do i = 1,nx
        mthis = mat_map(i)
        qdown(i) = sum(xslib%mat(mthis)%scatter(1:g-1,g,idxn+1) * phi(i,1:g-1,idxn+1))
      enddo ! i = 1,nx
      qdown = qdown * hx**2
    endif
  endsubroutine transport_build_downscatter

  ! NOTE: this is an exact copy of diffusion_build_fsource ...
  subroutine transport_build_fsource(nx, hx, mat_map, xslib, phi0, qfiss)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: hx
    integer(ik), intent(in) :: mat_map(:)
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(in) :: phi0(:,:) ! (nx, ngroup)
    real(rk), intent(out) :: qfiss(:,:) ! (nx, ngroup)

    integer(ik) :: i
    integer(ik) :: mthis

    do i = 1,nx
      mthis = mat_map(i)
      if (xslib%mat(mthis)%is_fiss) then
        qfiss(i,:) = xslib%mat(mthis)%chi(:) * sum(xslib%mat(mthis)%nusf(:) * phi0(i,:))
      else
        qfiss(i,:) = 0d0
      endif
    enddo ! i = 1,nx
    qfiss = qfiss * hx**2

  endsubroutine transport_build_fsource

  subroutine transport_build_next_source(nx, ngroup, sigma_tr, phi, qnext)
    integer(ik), intent(in) :: nx
    integer(ik), intent(in) :: ngroup
    real(rk), intent(in) :: sigma_tr(:,:,:) ! (nx, ngroup, pnorder)
    real(rk), intent(in) :: phi(:,:,:) ! (nx, ngroup, pnorder)
    real(rk), intent(out) :: qnext(:,:,:) ! (nx, ngroup, neven)

    integer(ik) :: i, g, n
    integer(ik) :: idxn
    real(rk) :: xn, xmul
    real(rk) :: kaprev, kathis, kanext
    real(rk) :: daprev, danext

    qnext(:,:,neven) = 0d0
    do n = 1,neven-1
      idxn = 2*(n-1)
      xn = real(idxn, rk)
      xmul = (xn+1d0)*(xn+2d0)/((2d0*xn+1d0)*(2d0*xn+3d0))
      do g = 1,ngroup
        qnext(1,g,n) = 0d0
        do i = 2,nx-1

          kaprev = xmul/sigma_tr(i-1,g,idxn+1+1)
          kathis = xmul/sigma_tr(i  ,g,idxn+1+1)
          kanext = xmul/sigma_tr(i+1,g,idxn+1+1)

          daprev = 2 * kathis*kaprev / (kathis + kaprev)
          danext = 2 * kathis*kanext / (kathis + kanext)

          qnext(i, g, n) = phi(i-1,g,idxn+1+2) * daprev - phi(i,g,idxn+1+2) * (daprev + danext) + phi(i+1,g,idxn+1+2) * danext

        enddo
        qnext(nx,g,n) = 0d0
      enddo ! g = 1,ngroup
    enddo ! n = 1,neven

  endsubroutine transport_build_next_source

  subroutine transport_build_prev_source(nx, ngroup, sigma_tr, phi, idxn, qprev)
    integer(ik), intent(in) :: nx
    integer(ik), intent(in) :: ngroup
    real(rk), intent(in) :: sigma_tr(:,:,:) ! (nx, ngroup, pnorder)
    real(rk), intent(in) :: phi(:,:,:) ! (nx, ngroup, pnorder)
    integer(ik), intent(in) :: idxn
    real(rk), intent(out) :: qprev(:,:) ! (nx, ngroup)

    integer(ik) :: g, i
    real(rk) :: xn, xmul
    real(rk) :: kbprev, kbthis, kbnext
    real(rk) :: dbprev, dbnext

    xn = real(idxn, rk)
    xmul = (xn**2-xn)/(4d0*xn**2 - 1d0)

    do g = 1,ngroup
      qprev(1,g) = 0d0
      do i = 2,nx-1

        kbprev = xmul/sigma_tr(i-1,g,idxn+1-1)
        kbthis = xmul/sigma_tr(i  ,g,idxn+1-1)
        kbnext = xmul/sigma_tr(i+1,g,idxn+1-1)

        dbprev = 2 * kbthis*kbprev / (kbthis + kbprev)
        dbnext = 2 * kbthis*kbnext / (kbthis + kbnext)

        qprev(i,g) = phi(i-1,g,idxn+1-2) * dbprev - phi(i,g,idxn+1-2) * (dbprev + dbnext) + phi(i+1,g,idxn+1-2) * dbnext

      enddo ! i = 2,nx-1
      qprev(nx,g) = 0d0
    enddo ! g = 1,ngroup

  endsubroutine transport_build_prev_source

  ! NOTE: this is an exact copy of diffusion_fission_sumation ...
  real(rk) function transport_fission_summation(nx, mat_map, xslib, phi0)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    integer(ik), intent(in) :: mat_map(:)
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(in) :: phi0(:,:) ! (nx, ngroup) -- scalar flux

    integer(ik) :: i, mthis
    real(rk) :: xsum

    xsum = 0d0
    do i = 1,nx
      mthis = mat_map(i)
      if (xslib%mat(mthis)%is_fiss) then
        xsum = xsum + sum(xslib%mat(mthis)%nusf(:) * phi0(i,:))
      endif
    enddo ! i = 1,nx
    transport_fission_summation = xsum

  endfunction transport_fission_summation

  subroutine transport_power_iteration(nx, hx, mat_map, xslib, k_tol, phi_tol, max_iter, pnorder, keff, phi)
    use xs, only : XSLibrary
    use linalg, only : trid
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: hx
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(in) :: k_tol, phi_tol 
    integer(ik), intent(in) :: max_iter
    integer(ik), intent(in) :: pnorder
    real(rk), intent(out) :: keff
    real(rk), intent(out) :: phi(:,:,:) ! (nx, ngroup, nmoment)

    ! matrix
    real(rk), allocatable :: sub(:,:,:), dia(:,:,:), sup(:,:,:) ! (nx, ngroup, neven)
    real(rk), allocatable :: sub_copy(:,:,:), dia_copy(:,:,:), sup_copy(:,:,:)
    ! neutron source
    real(rk), allocatable :: fsource(:,:) ! (nx, ngroup) -- all p0
    real(rk), allocatable :: upsource(:,:) ! (nx, ngroup) -- just this moment
    real(rk), allocatable :: downsource(:) ! (nx) -- just this moment & this group
    ! moment source
    real(rk), allocatable :: pn_next_source(:,:,:) ! (nx, ngroup, neven)
    real(rk), allocatable :: pn_prev_source(:,:) ! (nx, ngroup) -- just this moment
    ! combined source
    real(rk), allocatable :: q(:)

    integer(ik) :: iter
    real(rk) :: k_old, fsum, fsum_old
    real(rk), allocatable :: flux_old(:,:) ! (nx, ngroup) -- all p0
    real(rk) :: delta_k, delta_phi

    integer(ik) :: n, g, i
    integer(ik) :: mthis

    if (mod(pnorder,2) /= 1) then
      stop 'pnorder must be odd'
    endif

    neven = max((pnorder + 1) / 2, 1)

    allocate(sub(nx-1,xslib%ngroup,neven), dia(nx,xslib%ngroup,neven), sup(nx-1,xslib%ngroup,neven))
    allocate(sub_copy(nx-1,xslib%ngroup,neven), dia_copy(nx,xslib%ngroup,neven), sup_copy(nx-1,xslib%ngroup,neven))

    allocate(sigma_tr(nx,xslib%ngroup,pnorder+1))
    ! initialize to total xs
    ! reasonable? I'm not sure
    ! but necessary for damping
    do n = 1,pnorder+1
      do g = 1,xslib%ngroup
        do i = 1,nx
          mthis = mat_map(i)
          sigma_tr(i,g,n) = xslib%mat(mthis)%sigma_t(g)
        enddo ! i = 1,nx
      enddo ! g = 1,xslib%ngroup
    enddo ! n = 1,pnorder+1

    allocate(fsource(nx,xslib%ngroup))
    allocate(upsource(nx,xslib%ngroup))
    allocate(downsource(nx))
    allocate(pn_next_source(nx,xslib%ngroup,neven))
    allocate(pn_prev_source(nx,xslib%ngroup))
    allocate(q(nx))

    allocate(flux_old(nx,xslib%ngroup))

    keff = 1d0
    phi = 1d0
    fsum = 1d0

    write(*,*) "=== PN TRANSPORT POWER ITERATION ==="

    do iter = 1,max_iter
      k_old = keff
      flux_old = phi(:,:,1)
      fsum_old = fsum

      if (iter > dead_iter) then
        ! transport xs update
        call transport_build_transportxs(nx, mat_map, xslib, pnorder, phi, sigma_tr)
      endif

      ! even update
      ! matrix must be rebuilt every time since we update transport xs
      call transport_build_matrix(nx, hx, mat_map, xslib, neven, sub, dia, sup)

      call transport_build_next_source(nx, xslib%ngroup, sigma_tr, phi, pn_next_source)

      do n = 1,neven

        if (n == 1) then
          call transport_build_fsource(nx, hx, mat_map, xslib, phi(:,:,1), fsource)
        else ! (n > 1)
          call transport_build_prev_source(nx, xslib%ngroup, sigma_tr, phi, 2*(n-1), pn_prev_source)
        endif

        call transport_build_upscatter(nx, hx, mat_map, xslib, phi, 2*(n-1), upsource)

        do g = 1,xslib%ngroup

          call transport_build_downscatter(nx, hx, mat_map, xslib, phi, 2*(n-1), g, downsource)

          q = upsource(:,g) + downsource

          if (n == 1) then
            q = q + fsource(:,g)/keff
          else
            q = q + pn_prev_source(:,g)
          endif

          if (n < neven) then
            q = q + pn_next_source(:,g,n)
          endif

          ! SOLVE
          ! need to store copies, trid uses them as scratch space
          sub_copy(:,g,n) = sub(:,g,n)
          dia_copy(:,g,n) = dia(:,g,n)
          sup_copy(:,g,n) = sup(:,g,n)
          call trid(nx, sub_copy(:,g,n), dia_copy(:,g,n), sup_copy(:,g,n), q, phi(:,g,2*(n-1)+1))

        enddo ! g = 1,ngroup
      enddo ! n = 1,neven

      ! odd update
      call transport_odd_update(nx, hx, xslib%ngroup, pnorder, sigma_tr, phi)

      ! eigenvalue update
      fsum = transport_fission_summation(nx, mat_map, xslib, phi(:,:,1))
      if (iter > 1) keff = keff * fsum / fsum_old
      delta_k = abs(keff - k_old)
      ! only scalar flux, all groups
      delta_phi = maxval(abs(phi(:,:,1) - flux_old)) / maxval(phi(:,:,1))

      if ((keff < 0d0) .or. (keff > 2d0)) then
        write(*,*) 'keff', keff
        stop 'invalid keff'
      endif

      write(*,'(a,i4,a,es8.1,a,es8.1,a,f8.6)') &
        'it=', iter, ' dx=', delta_k, ' dphi=', delta_phi, ' keff=', keff

      if ((iter > dead_iter) .and. (delta_k < k_tol) .and. (delta_phi < phi_tol)) then
        write(*,*) 'CONVERGENCE!'
        write(*,*)
        exit
      endif

    enddo ! iter = 1,max_iter
    
    if (iter > max_iter) then
      write(*,*) 'WARNING: failed to converge'
    endif

    deallocate(sub, dia, sup)
    deallocate(sub_copy, dia_copy, sup_copy)
    deallocate(fsource, upsource, downsource, q)
    deallocate(pn_next_source, pn_prev_source)
    deallocate(flux_old)

  endsubroutine transport_power_iteration

endmodule transport
