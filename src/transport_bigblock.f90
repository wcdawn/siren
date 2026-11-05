module transport_bigblock
use kind, only : rk, ik
implicit none

private
public :: transport_bigblock_power_iteration

real(rk), allocatable :: invtransmat(:,:,:,:) ! (ngroup,ngroup,pnorder+1,xslib%niso)

contains

  subroutine transport_bigblock_build_matrix( &
      nx, dx, mat_map, xslib, &
      boundary_left, boundary_right, albedo_alpha, neven, sub, dia, sup)
    use xs, only : XSLibrary
    use exception_handler, only : exception_fatal
    use linalg, only : inv
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    character(*), intent(in) :: boundary_left, boundary_right
    real(rk), intent(in) :: albedo_alpha
    integer(ik), intent(in) :: neven
    real(rk), intent(out) :: sub(:,:,:) ! (dimlen, dimlen, nx-1)
    real(rk), intent(out) :: dia(:,:,:) ! (dimlen, dimlen, nx)
    real(rk), intent(out) :: sup(:,:,:) ! (dimlen, dimlen, nx-1)

    integer(ik) :: i, g, n
    integer(ik) :: mprev, mthis, mnext
    integer(ik) :: idxn, ilo, ihi, idx, jlo, jhi, klo, khi

    real(rk) :: xn, xmul_next, xmul_prev, xcup_next, xcup_prev

    real(rk), allocatable :: dhat_prev(:,:), dhat_this(:,:), dhat_next(:,:)
    real(rk), allocatable :: anext(:,:), aprev(:,:)

    real(rk), allocatable :: fhat_prev(:,:), fhat_this(:,:), fhat_next(:,:)
    real(rk), allocatable :: bnext(:,:), bprev(:,:)

    real(rk), allocatable :: ghat_prev(:,:), ghat_this(:,:), ghat_next(:,:)
    real(rk), allocatable :: cnext(:,:), cprev(:,:)

    real(rk), allocatable :: matinv(:,:)

    allocate( &
      dhat_prev(xslib%ngroup,xslib%ngroup), &
      dhat_this(xslib%ngroup,xslib%ngroup), &
      dhat_next(xslib%ngroup,xslib%ngroup) &
    )
    allocate( &
      fhat_prev(xslib%ngroup,xslib%ngroup), &
      fhat_this(xslib%ngroup,xslib%ngroup), &
      fhat_next(xslib%ngroup,xslib%ngroup) &
    )
    allocate( &
      ghat_prev(xslib%ngroup,xslib%ngroup), &
      ghat_this(xslib%ngroup,xslib%ngroup), &
      ghat_next(xslib%ngroup,xslib%ngroup) &
    )
    allocate( &
      anext(xslib%ngroup, xslib%ngroup), &
      aprev(xslib%ngroup, xslib%ngroup) &
    )
    allocate( &
      bnext(xslib%ngroup, xslib%ngroup), &
      bprev(xslib%ngroup, xslib%ngroup) &
    )
    allocate( &
      cnext(xslib%ngroup, xslib%ngroup), &
      cprev(xslib%ngroup, xslib%ngroup) &
    )
    allocate(matinv(xslib%ngroup,xslib%ngroup))

    ! BC at x=0, i=1
    mthis = mat_map(1)
    mnext = mat_map(2)
    do n = 1,neven
      idxn = 2*(n-1)
      xn = real(idxn, rk)
      xmul_next = (xn+1.0_rk)**2 / ((2.0_rk*xn+1.0_rk)*(2.0_rk*xn+3.0_rk))
      xmul_prev = xn**2 / (4.0_rk*xn**2 - 1.0_rk)
      xcup_next = (xn+1.0_rk)*(xn+2.0_rk)/((2.0_rk*xn+1.0_rk)*(2.0_rk*xn+3.0_rk))
      xcup_prev = (xn**2-xn)/(4.0_rk*xn**2 - 1.0_rk)
      ! XNEXT
      dhat_this = xmul_next * invtransmat(:,:,idxn+1+1,mthis)
      dhat_next = xmul_next * invtransmat(:,:,idxn+1+1,mnext)
      fhat_this = xcup_next * invtransmat(:,:,idxn+1+1,mthis)
      fhat_next = xcup_next * invtransmat(:,:,idxn+1+1,mnext)
      ! XPREV
      if (idxn > 0) then
        dhat_this = dhat_this + xmul_prev * invtransmat(:,:,idxn+1-1,mthis)
        dhat_next = dhat_next + xmul_prev * invtransmat(:,:,idxn+1-1,mnext)
        ghat_this = xcup_prev * invtransmat(:,:,idxn+1-1,mthis)
        ghat_next = xcup_prev * invtransmat(:,:,idxn+1-1,mnext)
      endif

      anext = dhat_next/dx(2) + dhat_this/dx(1)
      call inv(xslib%ngroup, anext, matinv)
      anext = 2.0_rk/dx(1)/dx(2) * matmul(matmul(dhat_this, matinv), dhat_next)

      bnext = fhat_next/dx(2) + fhat_this/dx(1)
      call inv(xslib%ngroup, bnext, matinv)
      bnext = 2.0_rk/dx(1)/dx(2) * matmul(matmul(fhat_this, matinv), fhat_next)

      if (idxn > 0) then
        cnext = ghat_next/dx(2) + ghat_this/dx(1)
        call inv(xslib%ngroup, cnext, matinv)
        cnext = 2.0_rk/dx(1)/dx(2) * matmul(matmul(ghat_this, matinv), ghat_next)
      endif

      ilo = 1 + xslib%ngroup * (n-1)
      ihi = xslib%ngroup * n
      sup(ilo:ihi,ilo:ihi,1) = -anext

      jlo = 1 + xslib%ngroup * (n+1-1)
      jhi = xslib%ngroup * (n+1)
      if (n < neven) then
        sup(ilo:ihi,jlo:jhi,1) = -bnext
      endif

      klo = 1 + xslib%ngroup * (n-1-1)
      khi = xslib%ngroup * (n-1)
      if (idxn > 0) then
        sup(ilo:ihi,klo:khi,1) = -cnext
      endif

      select case (boundary_left)
        case ('mirror')
          dia(ilo:ihi,ilo:ihi,1) = anext
          if (n < neven) then
            dia(ilo:ihi,jlo:jhi,1) = bnext
          endif
          if (idxn > 0) then
            dia(ilo:ihi,klo:khi,1) = cnext
          endif
        case default
          call exception_fatal(&
            'unknown boundary_left in transport_bigblock_buld_matrix: ' &
            // trim(adjustl(boundary_left)))
      endselect
    enddo ! n = 1,neven

    do i = 2,nx-1
      mprev = mat_map(i-1)
      mthis = mat_map(i)
      mnext = mat_map(i+1)

      do n = 1,neven
        idxn = 2*(n-1)
        xn = real(idxn, rk)
        xmul_next = (xn+1.0_rk)**2 / ((2.0_rk*xn+1.0_rk)*(2.0_rk*xn+3.0_rk))
        xmul_prev = xn**2 / (4.0_rk*xn**2 - 1.0_rk)
        xcup_next = (xn+1.0_rk)*(xn+2.0_rk)/((2.0_rk*xn+1.0_rk)*(2.0_rk*xn+3.0_rk))
        xcup_prev = (xn**2-xn)/(4.0_rk*xn**2 - 1.0_rk)

        ! XNEXT
        dhat_prev = xmul_next * invtransmat(:,:,idxn+1+1,mprev)
        dhat_this = xmul_next * invtransmat(:,:,idxn+1+1,mthis)
        dhat_next = xmul_next * invtransmat(:,:,idxn+1+1,mnext)
        fhat_prev = xcup_next * invtransmat(:,:,idxn+1+1,mprev)
        fhat_this = xcup_next * invtransmat(:,:,idxn+1+1,mthis)
        fhat_next = xcup_next * invtransmat(:,:,idxn+1+1,mnext)

        ! XPREV
        if (idxn > 0) then
          dhat_prev = dhat_prev + xmul_prev * invtransmat(:,:,idxn+1-1,mprev)
          dhat_this = dhat_this + xmul_prev * invtransmat(:,:,idxn+1-1,mthis)
          dhat_next = dhat_next + xmul_prev * invtransmat(:,:,idxn+1-1,mnext)
          ghat_prev = xcup_prev * invtransmat(:,:,idxn+1-1,mprev)
          ghat_this = xcup_prev * invtransmat(:,:,idxn+1-1,mthis)
          ghat_next = xcup_prev * invtransmat(:,:,idxn+1-1,mnext)
        endif

        aprev = dhat_prev/dx(i-1) + dhat_this/dx(i)
        call inv(xslib%ngroup, aprev, matinv)
        aprev = 2.0_rk/dx(i)/dx(i-1) * matmul(matmul(dhat_this, matinv), dhat_prev)
        anext = dhat_next/dx(i+1) + dhat_this/dx(i)
        call inv(xslib%ngroup, anext, matinv)
        anext = 2.0_rk/dx(i)/dx(i+1) * matmul(matmul(dhat_this, matinv), dhat_next)

        bprev = fhat_prev/dx(i-1) + fhat_this/dx(i)
        call inv(xslib%ngroup, bprev, matinv)
        bprev = 2.0_rk/dx(i)/dx(i-1) * matmul(matmul(fhat_this, matinv), fhat_prev)
        bnext = fhat_next/dx(i+1) + fhat_this/dx(i)
        call inv(xslib%ngroup, bnext, matinv)
        bnext = 2.0_rk/dx(i)/dx(i+1) * matmul(matmul(fhat_this, matinv), fhat_next)

        if (idxn > 0) then
          cprev = ghat_prev/dx(i-1) + ghat_this/dx(i)
          call inv(xslib%ngroup, cprev, matinv)
          cprev = 2.0_rk/dx(i)/dx(i-1) * matmul(matmul(ghat_this, matinv), ghat_prev)
          cnext = ghat_next/dx(i+1) + ghat_this/dx(i)
          call inv(xslib%ngroup, cnext, matinv)
          cnext = 2.0_rk/dx(i)/dx(i+1) * matmul(matmul(ghat_this, matinv), ghat_next)
        endif

        ilo = 1 + xslib%ngroup * (n-1)
        ihi = xslib%ngroup * n
        sub(ilo:ihi,ilo:ihi,i-1) = -aprev
        dia(ilo:ihi,ilo:ihi,i) = aprev + anext
        sup(ilo:ihi,ilo:ihi,i) = -anext

        jlo = 1 + xslib%ngroup * (n+1-1)
        jhi = xslib%ngroup * (n+1)
        if (n < neven) then
          sub(ilo:ihi,jlo:jhi,i-1) = -bprev
          dia(ilo:ihi,jlo:jhi,i) = bprev + bnext
          sup(ilo:ihi,jlo:jhi,i) = -bnext
        endif

        klo = 1 + xslib%ngroup * (n-1-1)
        khi = xslib%ngroup * (n-1)
        if (idxn > 0) then
          sub(ilo:ihi,klo:khi,i-1) = -cprev
          dia(ilo:ihi,klo:khi,i) = cprev + cnext
          sup(ilo:ihi,klo:khi,i) = -cnext
        endif

      enddo ! n = 1,neven
    enddo ! i = 2,nx-1

    ! BC at x=L, i=N
    mprev = mat_map(nx-1)
    mthis = mat_map(nx)
    do n = 1,neven
      idxn = 2*(n-1)
      xn = real(idxn, rk)
      xmul_next = (xn+1.0_rk)**2 / ((2.0_rk*xn+1.0_rk)*(2.0_rk*xn+3.0_rk))
      xmul_prev = xn**2 / (4.0_rk*xn**2 - 1.0_rk)
      xcup_next = (xn+1.0_rk)*(xn+2.0_rk)/((2.0_rk*xn+1.0_rk)*(2.0_rk*xn+3.0_rk))
      xcup_prev = (xn**2-xn)/(4.0_rk*xn**2 - 1.0_rk)
      ! XNEXT
      dhat_prev = xmul_next * invtransmat(:,:,idxn+1+1,mprev)
      dhat_this = xmul_next * invtransmat(:,:,idxn+1+1,mthis)
      fhat_this = xcup_next * invtransmat(:,:,idxn+1+1,mthis)
      fhat_next = xcup_next * invtransmat(:,:,idxn+1+1,mnext)
      ! XPREV
      if (idxn > 0) then
        dhat_prev = dhat_prev + xmul_prev * invtransmat(:,:,idxn+1-1,mprev)
        dhat_this = dhat_this + xmul_prev * invtransmat(:,:,idxn+1-1,mthis)
        ghat_this = xcup_prev * invtransmat(:,:,idxn+1-1,mthis)
        ghat_next = xcup_prev * invtransmat(:,:,idxn+1-1,mnext)
      endif

      aprev = dhat_prev/dx(nx-1) + dhat_this/dx(nx)
      call inv(xslib%ngroup, aprev, matinv)
      aprev = 2.0_rk/dx(nx)/dx(nx-1) * matmul(matmul(dhat_this, matinv), dhat_prev)

      bprev = fhat_prev/dx(nx-1) + fhat_this/dx(nx)
      call inv(xslib%ngroup, bprev, matinv)
      bprev = 2.0_rk/dx(nx)/dx(nx-1) * matmul(matmul(fhat_this, matinv), fhat_prev)

      if (idxn > 0) then
        cprev = ghat_prev/dx(nx-1) + ghat_this/dx(nx)
        call inv(xslib%ngroup, cprev, matinv)
        cprev = 2.0_rk/dx(nx)/dx(nx-1) * matmul(matmul(ghat_this, matinv), ghat_prev)
      endif

      ilo = 1 + xslib%ngroup * (n-1)
      ihi = xslib%ngroup * n
      sub(ilo:ihi,ilo:ihi,nx-1) = -aprev

      jlo = 1 + xslib%ngroup * (n+1-1)
      jhi = xslib%ngroup * (n+1)
      if (n < neven) then
        sub(ilo:ihi,jlo:jhi,nx-1) = -bprev
      endif

      klo = 1 + xslib%ngroup * (n-1-1)
      khi = xslib%ngroup * (n-1)
      if (idxn > 0) then
        sub(ilo:ihi,klo:khi,nx-1) = -cprev
      endif

      select case (boundary_right)
        case ('mirror')
          dia(ilo:ihi,ilo:ihi,nx) = aprev
          if (n < neven) then
            dia(ilo:ihi,jlo:jhi,nx) = bprev
          endif
          if (idxn > 0) then
            dia(ilo:ihi,klo:khi,nx) = cprev
          endif
        case ('zero')
          dia(ilo:ihi,ilo:ihi,nx) = aprev + 2.0_rk * dhat_this/dx(nx)
          if (n < neven) then
            dia(ilo:ihi,jlo:jhi,nx) = bprev + 2.0_rk * fhat_this/dx(nx)
          endif
          if (idxn > 0) then
            dia(ilo:ihi,klo:khi,nx) = cprev + 2.0_rk * ghat_this/dx(nx)
          endif
        case default
          call exception_fatal(&
            'unknown boundary_right in transport_bigblock_build_matrix: ' &
            // trim(adjustl(boundary_right)))
      endselect
    enddo ! n = 1,neven

    ! adjust diagonal block with scattering moments as appropriate
    do n = 1,neven
      idxn = 2*(n-1)
      do i = 1,nx
        mthis = mat_map(i)
        do g = 1,xslib%ngroup
          idx = g + xslib%ngroup * (n-1)
          dia(idx,idx,i) = dia(idx,idx,i) + xslib%mat(mthis)%sigma_t(g) * dx(i)
        enddo ! g = 1,xslib%ngroup
        if (idxn+1 <= xslib%nmoment) then
          ilo = 1 + xslib%ngroup * (n-1)
          ihi = xslib%ngroup * n
          dia(ilo:ihi,ilo:ihi,i) = dia(ilo:ihi,ilo:ihi,i) &
            - transpose(xslib%mat(mthis)%scatter(:,:,idxn+1)) * dx(i)
        endif
      enddo ! i = 1,nx
    enddo ! n = 1,neven

    deallocate(dhat_prev, dhat_this, dhat_next)
    deallocate(anext, aprev)
    deallocate(fhat_prev, fhat_this, fhat_next)
    deallocate(bnext, bprev)
    deallocate(ghat_prev, ghat_this, ghat_next)
    deallocate(cnext, cprev)
    deallocate(matinv)
  endsubroutine transport_bigblock_build_matrix

  subroutine transport_bigblock_build_fsource(nx, dx, mat_map, xslib, phi_bigblock, q)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(in) :: phi_bigblock(:,:) ! (dimlen,nx)
    real(rk), intent(inout) :: q(:,:) ! (dimlen,nx)

    integer(ik) :: i, mthis

    do i = 1,nx
      mthis = mat_map(i)
      if (xslib%mat(mthis)%is_fiss) then
        q(1:xslib%ngroup,i) = xslib%mat(mthis)%chi &
          * sum(xslib%mat(mthis)%nusf * phi_bigblock(1:xslib%ngroup,i)) * dx(i)
      endif
    enddo ! i = 1,nx
  endsubroutine transport_bigblock_build_fsource

  real(rk) pure function transport_bigblock_fission_summation(&
      nx, dx, mat_map, xslib, phi_bigblock)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(in) :: phi_bigblock(:,:) ! (dimlen, nx)

    integer(ik) :: i, mthis
    real(rk) :: xsum

    xsum = 0.0_rk
    do i = 1,nx
      mthis = mat_map(i)
      if (xslib%mat(mthis)%is_fiss) then
        xsum = xsum + &
          sum(xslib%mat(mthis)%nusf * phi_bigblock(1:xslib%ngroup,i)) * dx(i)
      endif
    enddo ! i = 1,nx
    transport_bigblock_fission_summation = xsum
  endfunction transport_bigblock_fission_summation

  subroutine transport_bigblock_power_iteration( &
      nx, dx, mat_map, xslib, boundary_left, boundary_right, albedo_coeff, &
      calc_type, k_tol, phi_tol, max_iter, pnorder, keff, sigma_tr, phi)
    use xs, only : XSLibrary
    use exception_handler, only : exception_fatal, exception_warning
    use output, only : output_write
    use transport_block, only : transport_invtransmat
    use timer, only : timer_start, timer_stop
    use albedo, only : albedo_calculate_alpha
    use linalg, only : trid_block
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    character(*), intent(in) :: boundary_left, boundary_right
    real(rk), intent(in) :: albedo_coeff ! A \in [-1,1]
    character(*), intent(in) :: calc_type
    real(rk), intent(in) :: k_tol, phi_tol
    integer(ik), intent(in) :: max_iter
    integer(ik), intent(in) :: pnorder
    real(rk), intent(out) :: keff
    real(rk), intent(out) :: sigma_tr(:,:,:) ! (nx, ngroup, nmoment)
    real(rk), intent(out) :: phi(:,:,:) ! (nx, ngroup, nmoment)

    integer(ik) :: neven, dimlen
    integer(ik) :: i, g, n, idx, idxn
    integer(ik) :: iter

    real(rk) :: k_old
    real(rk) :: fsum, fsum_old
    real(rk) :: delta_k, delta_phi
    real(rk) :: albedo_alpha ! \alpha \in [0,\infty)

    real(rk), allocatable :: sub(:,:,:) ! (dimlen, dimlen, nx-1)
    real(rk), allocatable :: dia(:,:,:) ! (dimlen, dimlen, nx)
    real(rk), allocatable :: sup(:,:,:) ! (dimlen, dimlen, nx-1)
    real(rk), allocatable :: q(:,:)   ! (dimlen, nx)
    real(rk), allocatable :: phi_bigblock(:,:), phi_bigblock_old(:,:) ! (dimlen, nx)

    character(1024) :: line

    if (mod(pnorder,2) /= 1) then
      call exception_fatal('pnorder must be odd')
    endif
    neven = (pnorder+1)/2

    dimlen = neven * xslib%ngroup
    allocate(&
      sub(dimlen,dimlen,nx-1), &
      dia(dimlen,dimlen,nx), &
      sup(dimlen,dimlen,nx-1) &
    )
    allocate(q(dimlen,nx))
    q = 0.0_rk

    allocate(phi_bigblock(dimlen,nx), phi_bigblock_old(dimlen,nx))

    ! The indexing in this method is:
    ! idx = g + NG * (n-1) + NG * NE * (i-1)
    ! g = 1,NG -- group index
    ! n = 1,NE -- even-moment index (idxn=2*(n-1))
    ! i = 1,NX -- cell index
    ! This is the same indexing as Eq. (56) in Hamilton & Evans (JCP).
    ! This indexing minimizes the bandwith for a problem that has many groups,
    ! relatively few (even) moments, and an arbitrary number of spatial cells.

    ! set scalar flux to "flat" 1.0
    ! set higher moments to zero
    phi_bigblock = 0.0_rk
    phi_bigblock(1:xslib%ngroup,:) = 1.0_rk

    keff = 1.0_rk
    fsum = 1.0_rk

    call output_write('=== PN TRANSPORT BIG-BLOCK POWER ITERATION ===')

    call timer_start('transport_invtransmat')
    allocate(invtransmat(xslib%ngroup,xslib%ngroup,pnorder+1,xslib%niso))
    do i = 1,xslib%niso
      do n = 1,pnorder+1
        call transport_invtransmat(xslib%mat(i), n-1, invtransmat(:,:,n,i))
      enddo ! n = 1,pnorder+1
    enddo ! i = 1,xslib%niso
    call timer_stop('transport_invtransmat')

    albedo_alpha = albedo_calculate_alpha(albedo_coeff)
    call output_write('building transport matrix')
    call timer_start('transport_build_matrix')
    call transport_bigblock_build_matrix(nx, dx, mat_map, xslib, &
      boundary_left, boundary_right, albedo_alpha, neven, sub, dia, sup)
    call timer_stop('transport_build_matrix')

    if (calc_type /= 'eigenvalue') then
      call exception_fatal('In bigblock, calc_type must be eigenvalue.')
    endif

    call output_write('beginning iterations')

    do iter = 1,max_iter
      k_old = keff
      phi_bigblock_old = phi_bigblock
      fsum_old = fsum

      call timer_start('transport_source')
      call transport_bigblock_build_fsource(nx, dx, mat_map, xslib, phi_bigblock, q)
      call timer_stop('transport_source')

      q = q / keff

      call timer_start('transport_block_tridiagonal')
      call trid_block(nx, dimlen, sub, dia, sup, q, phi_bigblock)
      call timer_stop('transport_block_tridiagonal')

      call timer_start('transport_convergence')
      fsum = transport_bigblock_fission_summation(nx, dx, mat_map, xslib, phi_bigblock)
      if (iter > 1) keff = keff * fsum / fsum_old
      delta_k = abs(keff - k_old)
      delta_phi = maxval(abs(phi_bigblock - phi_bigblock_old)) &
        / maxval(phi_bigblock_old)
      call timer_stop('transport_convergence')

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
    
    ! TODO calculate odd
    ! TODO calculate sigma_tr

    ! copy before exit to reorder
    do i = 1,nx
      do n = 1,neven
        idxn = 2*(n-1)
        do g = 1,xslib%ngroup
          idx = g + xslib%ngroup * (n-1)
          phi(i,g,idxn+1) = phi_bigblock(idx,i)
        enddo ! 1,xslib%ngroup
      enddo ! n = 1,neven
    enddo ! i = 1,nx

    deallocate(sub, dia, sup)
    deallocate(q)
    deallocate(invtransmat)
  endsubroutine transport_bigblock_power_iteration

endmodule transport_bigblock
