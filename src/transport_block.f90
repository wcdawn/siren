module transport_block
use kind
use exception_handler
implicit none

private
public :: transport_block_power_iteration

real(rk), allocatable :: invtransmat(:,:,:,:) ! (ngroup,ngroup,pnorder+1,xslib%niso)

contains

  subroutine transport_invtransmat(xsmat, idxn, mat)
    use xs, only : XSMaterial
    use linalg, only : inv
    type(XSMaterial), intent(in) :: xsmat
    integer(ik), intent(in) :: idxn
    real(rk), intent(out) :: mat(:,:) ! (ngroup,ngroup)

    integer(ik) :: g
    real(rk), allocatable :: matinv(:,:) ! (ngroup,ngroup)

    mat = 0d0

    do g = 1,xsmat%ngroup
      mat(g,g) = xsmat%sigma_t(g)
    enddo ! g = 1,xsmat%ngroup

    if (idxn+1 <= xsmat%nmoment) then
      ! we have a scattering moment
      mat = mat - transpose(xsmat%scatter(:,:,idxn+1))
      allocate(matinv(xsmat%ngroup,xsmat%ngroup))
      call inv(xsmat%ngroup, mat, matinv)
      mat = matinv
      deallocate(matinv)
    else
      do g = 1,xsmat%ngroup
        mat(g,g) = 1.0_rk / mat(g,g)
      enddo ! g = 1,xsmat%ngroup
    endif
  endsubroutine transport_invtransmat

  subroutine transport_block_build_matrix(nx, dx, mat_map, xslib, boundary_right, neven, sub, dia, sup)
    use xs, only : XSLibrary
    use linalg, only : inv
    use omp_lib, only : omp_get_num_threads, omp_get_thread_num
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    character(*), intent(in) :: boundary_right
    integer(ik), intent(in) :: neven
    real(rk), intent(out) :: sub(:,:,:,:) ! (ngroup,ngroup,nx-1,neven)
    real(rk), intent(out) :: dia(:,:,:,:) ! (ngroup,ngroup,nx,neven)
    real(rk), intent(out) :: sup(:,:,:,:) ! (ngroup,ngroup,nx-1,neven)

    integer(ik) :: i, g, n
    integer(ik) :: mprev, mthis, mnext
    integer(ik) :: idxn

    real(rk) :: xn, xmul_next, xmul_prev

    real(rk), allocatable :: dhat_prev(:,:,:), dhat_this(:,:,:), dhat_next(:,:,:)
    real(rk), allocatable :: anext(:,:,:), aprev(:,:,:)
    real(rk), allocatable :: matinv(:,:,:)

    integer :: nthread, myid

    !$omp parallel default(none) shared(nthread) private(myid)
    nthread = omp_get_num_threads()
    !$omp end parallel

    allocate(&
      dhat_prev(xslib%ngroup,xslib%ngroup,nthread), &
      dhat_this(xslib%ngroup,xslib%ngroup,nthread), &
      dhat_next(xslib%ngroup,xslib%ngroup,nthread)&
    )
    allocate(&
      anext(xslib%ngroup,xslib%ngroup,nthread), &
      aprev(xslib%ngroup,xslib%ngroup,nthread))
    allocate(matinv(xslib%ngroup,xslib%ngroup,nthread))

    myid = 1
    ! BC at x=0, i=1
    mthis = mat_map(1)
    mnext = mat_map(2)
    do n = 1,neven
      idxn = 2*(n-1)
      xn = real(idxn, rk)
      xmul_next = (xn+1.0_rk)**2 / ((2.0_rk*xn+1.0_rk)*(2.0_rk*xn+3.0_rk))
      xmul_prev = xn**2 / (4.0_rk*xn**2 - 1.0_rk)
      ! XNEXT
      dhat_this(:,:,myid) = xmul_next * invtransmat(:,:,idxn+1+1,mthis)
      dhat_next(:,:,myid) = xmul_next * invtransmat(:,:,idxn+1+1,mnext)
      ! XPREV
      if (idxn > 1) then
        dhat_this(:,:,myid) = dhat_this(:,:,myid) + xmul_prev * invtransmat(:,:,idxn+1-1,mthis)
        dhat_next(:,:,myid) = dhat_next(:,:,myid) + xmul_prev * invtransmat(:,:,idxn+1-1,mnext)
      endif
      anext(:,:,myid) = dhat_next(:,:,myid)/dx(2) + dhat_this(:,:,myid)/dx(1)
      call inv(xslib%ngroup, anext(:,:,myid), matinv(:,:,myid))
      anext(:,:,myid) = 2.0_rk/dx(1)/dx(2) &
        * matmul(matmul(dhat_this(:,:,myid), matinv(:,:,myid)), dhat_next(:,:,myid))
      dia(:,:,1,n) = anext(:,:,myid)
      sup(:,:,1,n) = -anext(:,:,myid)
    enddo ! n = 1,neven

    !$omp parallel do default(none) &
    !$omp shared(nx, neven) &
    !$omp private(myid, i, mprev, mthis, mnext) &
    !$omp shared (dx, mat_map, xslib) &
    !$omp private(n, idxn, xn, xmul_next, xmul_prev) &
    !$omp shared (dhat_prev, dhat_this, dhat_next, aprev, anext, matinv, invtransmat) &
    !$omp shared (sub, sup, dia)
    do i = 2,nx-1
      myid = omp_get_thread_num() + 1

      mprev = mat_map(i-1)
      mthis = mat_map(i)
      mnext = mat_map(i+1)

      do n = 1,neven
        idxn = 2*(n-1)

        xn = real(idxn, rk)
        xmul_next = (xn+1.0_rk)**2 / ((2.0_rk*xn+1.0_rk)*(2.0_rk*xn+3.0_rk))
        xmul_prev = xn**2 / (4.0_rk*xn**2 - 1.0_rk)

        ! XNEXT
        dhat_prev(:,:,myid) = xmul_next * invtransmat(:,:,idxn+1+1,mprev)
        dhat_this(:,:,myid) = xmul_next * invtransmat(:,:,idxn+1+1,mthis)
        dhat_next(:,:,myid) = xmul_next * invtransmat(:,:,idxn+1+1,mnext)

        ! XPREV
        if (idxn > 1) then
          dhat_prev(:,:,myid) = dhat_prev(:,:,myid) + xmul_prev * invtransmat(:,:,idxn+1-1,mprev)
          dhat_this(:,:,myid) = dhat_this(:,:,myid) + xmul_prev * invtransmat(:,:,idxn+1-1,mthis)
          dhat_next(:,:,myid) = dhat_next(:,:,myid) + xmul_prev * invtransmat(:,:,idxn+1-1,mnext)
        endif

        aprev(:,:,myid) = dhat_prev(:,:,myid)/dx(i-1) + dhat_this(:,:,myid)/dx(i)
        call inv(xslib%ngroup, aprev(:,:,myid), matinv(:,:,myid))
        aprev(:,:,myid) = 2.0_rk/dx(i)/dx(i-1) &
          * matmul(matmul(dhat_this(:,:,myid), matinv(:,:,myid)), dhat_prev(:,:,myid))

        anext(:,:,myid) = dhat_next(:,:,myid)/dx(i+1) + dhat_this(:,:,myid)/dx(i)
        call inv(xslib%ngroup, anext(:,:,myid), matinv(:,:,myid))
        anext(:,:,myid) = 2.0_rk/dx(i)/dx(i+1) &
          * matmul(matmul(dhat_this(:,:,myid), matinv(:,:,myid)), dhat_next(:,:,myid))

        sub(:,:,i-1,n) = -aprev(:,:,myid)
        dia(:,:,i,n) = aprev(:,:,myid) + anext(:,:,myid)
        sup(:,:,i,n) = -anext(:,:,myid)

      enddo ! n = 1,neven
    enddo ! i = 2,nx-1
    !$omp end parallel do

    ! BC at x=L, i=N
    myid = 1
    mprev = mat_map(nx-1)
    mthis = mat_map(nx)
    do n = 1,neven
      idxn = 2*(n-1)
      xn = real(idxn, rk)
      xmul_next = (xn+1.0_rk)**2 / ((2.0_rk*xn+1.0_rk)*(2.0_rk*xn+3.0_rk))
      xmul_prev = xn**2 / (4.0_rk*xn**2 - 1.0_rk)
      ! XNEXT
      dhat_prev(:,:,myid) = xmul_next * invtransmat(:,:,idxn+1+1,mprev)
      dhat_this(:,:,myid) = xmul_next * invtransmat(:,:,idxn+1+1,mthis)
      ! XPREV
      if (idxn > 1) then
        dhat_prev(:,:,myid) = dhat_prev(:,:,myid) + xmul_prev * invtransmat(:,:,idxn+1-1,mprev)
        dhat_this(:,:,myid) = dhat_this(:,:,myid) + xmul_prev * invtransmat(:,:,idxn+1-1,mthis)
      endif
      aprev(:,:,myid) = dhat_prev(:,:,myid)/dx(nx-1) + dhat_this(:,:,myid)/dx(nx)
      call inv(xslib%ngroup, aprev(:,:,myid), matinv(:,:,myid))
      aprev(:,:,myid) = 2.0_rk/dx(nx)/dx(nx-1) &
        * matmul(matmul(dhat_this(:,:,myid), matinv(:,:,myid)), dhat_prev(:,:,myid))
      select case (boundary_right)
        case ('mirror')
          sub(:,:,nx-1,n) = -aprev(:,:,myid)
          dia (:,:,nx,n) = aprev(:,:,myid)
        case ('zero')
          sub(:,:,nx-1,n) = -aprev(:,:,myid)
          dia(:,:,nx,n) = aprev(:,:,myid) + 2 * dhat_this(:,:,myid)/dx(nx)
        case default
          call exception_fatal('unknown boundary_right: ' // trim(adjustl(boundary_right)))
      endselect
    enddo ! n = 1,neven

    ! adjust diagonal block with scattering moments as appropriate
    do n = 1,neven
      idxn = 2*(n-1)
      do i = 1,nx
        mthis = mat_map(i)
        do g = 1,xslib%ngroup
          dia(g,g,i,n) = dia(g,g,i,n) + xslib%mat(mthis)%sigma_t(g) * dx(i)
        enddo
        if (idxn+1 <= xslib%nmoment) then
          dia(:,:,i,n) = dia(:,:,i,n) - transpose(xslib%mat(mthis)%scatter(:,:,idxn+1)) * dx(i)
        endif
      enddo ! i = 1,nx
    enddo ! n = 1,neven

    deallocate(matinv)
    deallocate(anext, aprev)
    deallocate(dhat_prev, dhat_this, dhat_next)
  endsubroutine transport_block_build_matrix

  subroutine transport_block_build_next_source(nx, dx, mat_map, xslib, boundary_right, neven, phi_block, pn_next_source)
    use xs, only : XSLibrary
    use linalg, only : inv
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    character(*), intent(in) :: boundary_right
    integer(ik), intent(in) :: neven
    real(rk), intent(in) :: phi_block(:,:,:) ! (ngroup,nx,pnorder+1)
    real(rk), intent(out) :: pn_next_source(:,:,:) ! (ngroup,nx,neven)

    integer(ik) :: i, n
    integer(ik) :: idxn
    integer(ik) :: mprev, mthis, mnext
    real(rk) :: xn, xmul

    real(rk), allocatable :: fhat_prev(:,:), fhat_this(:,:), fhat_next(:,:)
    real(rk), allocatable :: bnext(:,:), bprev(:,:)
    real(rk), allocatable :: matinv(:,:)

    allocate(&
      fhat_prev(xslib%ngroup,xslib%ngroup),&
      fhat_this(xslib%ngroup,xslib%ngroup),&
      fhat_next(xslib%ngroup,xslib%ngroup),&
    )
    allocate(bnext(xslib%ngroup,xslib%ngroup), bprev(xslib%ngroup,xslib%ngroup))
    allocate(matinv(xslib%ngroup,xslib%ngroup))

    ! NOTE: there is no "next" source after the last moment
    pn_next_source(1:xslib%ngroup,1:nx,neven) = 0d0

    ! BC at x=0, i=1
    mthis = mat_map(1)
    mnext = mat_map(2)
    do n = 1,neven-1
      idxn = 2*(n-1)
      xn = real(idxn, rk)
      xmul = (xn+1.0_rk)*(xn+2.0_rk)/((2.0_rk*xn+1.0_rk)*(2.0_rk*xn+3.0_rk))
      fhat_this = xmul * invtransmat(:,:,idxn+1+1,mthis)
      fhat_next = xmul * invtransmat(:,:,idxn+1+1,mnext)
      bnext = fhat_next/dx(2) + fhat_this/dx(1)
      call inv(xslib%ngroup, bnext, matinv)
      bnext = 2.0_rk/dx(1)/dx(2) * matmul(matmul(fhat_this, matinv), fhat_next)
      pn_next_source(:,1,n) = -matmul(bnext, phi_block(:,1,idxn+1+2)) &
        + matmul(bnext, phi_block(:,2,idxn+1+2))
    enddo ! n = 1,neven-1

    do i = 2,nx-1
      mprev = mat_map(i-1)
      mthis = mat_map(i)
      mnext = mat_map(i+1)
      do n = 1,neven-1
        idxn = 2*(n-1)
        xn = real(idxn, rk)
        xmul = (xn+1.0_rk)*(xn+2.0_rk)/((2.0_rk*xn+1.0_rk)*(2.0_rk*xn+3.0_rk))

        fhat_prev = xmul * invtransmat(:,:,idxn+1+1,mprev)
        fhat_this = xmul * invtransmat(:,:,idxn+1+1,mthis)
        fhat_next = xmul * invtransmat(:,:,idxn+1+1,mnext)

        bprev = fhat_prev/dx(i-1) + fhat_this/dx(i)
        call inv(xslib%ngroup, bprev, matinv)
        bprev = 2.0_rk/dx(i)/dx(i-1) * matmul(matmul(fhat_this, matinv), fhat_prev)

        bnext = fhat_next/dx(i+1) + fhat_this/dx(i)
        call inv(xslib%ngroup, bnext, matinv)
        bnext = 2.0_rk/dx(i)/dx(i+1) * matmul(matmul(fhat_this, matinv), fhat_next)

        pn_next_source(:,i,n) = matmul(bprev, phi_block(:,i-1,idxn+1+2)) &
          - matmul(bprev + bnext, phi_block(:,i,idxn+1+2)) &
          + matmul(bnext, phi_block(:,i+1,idxn+1+2))
      enddo
    enddo ! i = 2,nx-1

    ! BC at x=L, i=N
    mprev = mat_map(nx-1)
    mthis = mat_map(nx)
    do n = 1,neven-1
      idxn = 2*(n-1)
      xn = real(idxn, rk)
      xmul = (xn+1.0_rk)*(xn+2.0_rk)/((2.0_rk*xn+1.0_rk)*(2.0_rk*xn+3.0_rk))
      fhat_prev = xmul * invtransmat(:,:,idxn+1+1,mprev)
      fhat_this = xmul * invtransmat(:,:,idxn+1+1,mthis)
      bprev = fhat_prev/dx(nx-1) + fhat_this/dx(nx)
      call inv(xslib%ngroup, bprev, matinv)
      bprev = 2.0_rk/dx(nx-1)/dx(nx) * matmul(matmul(fhat_this, matinv), fhat_prev)
      select case (boundary_right)
        case ('mirror')
          pn_next_source(:,nx,n) = matmul(bprev, phi_block(:,nx-1,idxn+1+2)) &
            - matmul(bprev, phi_block(:,nx,idxn+1+2)) 
        case ('zero')
          pn_next_source(:,nx,n) = matmul(bprev, phi_block(:,nx-1,idxn+1+2)) &
            - matmul(bprev, phi_block(:,nx,idxn+1+2)) &
            - 2.0_rk/dx(nx) * matmul(fhat_this, phi_block(:,nx,idxn+1+2))
        case default
          call exception_fatal('unknown boundary_right in next_source: ' // trim(adjustl(boundary_right)))
      endselect
    enddo ! n = 1,neven

    deallocate(matinv)
    deallocate(bnext, bprev)
    deallocate(fhat_prev, fhat_this, fhat_next)
  endsubroutine transport_block_build_next_source

  subroutine transport_block_build_prev_source(nx, dx, mat_map, xslib, boundary_right, idxn, phi_block, pn_prev_source)
    use xs, only : XSLibrary
    use linalg, only : inv
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    character(*), intent(in) :: boundary_right
    integer(ik), intent(in) :: idxn
    real(rk), intent(in) :: phi_block(:,:,:) ! (ngroup,nx,pnorder+1)
    real(rk), intent(out) :: pn_prev_source(:,:) ! (ngroup,nx)

    integer(ik) :: i, n
    integer(ik) :: mprev, mthis, mnext

    real(rk) :: xn, xmul
    
    real(rk), allocatable :: ghat_prev(:,:), ghat_this(:,:), ghat_next(:,:)
    real(rk), allocatable :: cnext(:,:), cprev(:,:)
    real(rk), allocatable :: matinv(:,:)

    allocate(&
      ghat_prev(xslib%ngroup,xslib%ngroup),&
      ghat_this(xslib%ngroup,xslib%ngroup),&
      ghat_next(xslib%ngroup,xslib%ngroup),&
    )
    allocate(cnext(xslib%ngroup,xslib%ngroup), cprev(xslib%ngroup,xslib%ngroup))
    allocate(matinv(xslib%ngroup,xslib%ngroup))

    pn_prev_source(1:xslib%ngroup,1:nx) = 0d0
    if (idxn < 2) then
      return
    endif

    xn = real(idxn, rk)
    xmul = (xn**2-xn)/(4.0_rk*xn**2 - 1.0_rk)
    n = idxn/2 + 1

    ! BC at x=0, i=1
    mthis = mat_map(1)
    mnext = mat_map(2)
    ghat_this = xmul * invtransmat(:,:,idxn+1-1,mthis)
    ghat_next = xmul * invtransmat(:,:,idxn+1-1,mnext)
    cnext = ghat_next/dx(2) + ghat_this/dx(1)
    call inv(xslib%ngroup, cnext, matinv)
    cnext = 2.0_rk/dx(1)/dx(2) * matmul(matmul(ghat_this, matinv), ghat_next)
    pn_prev_source(:,1) = -matmul(cnext, phi_block(:,1,idxn+1-2)) &
      + matmul(cnext, phi_block(:,2,idxn+1-2))

    do i = 2,nx-1
      mprev = mat_map(i-1)
      mthis = mat_map(i)
      mnext = mat_map(i+1)

      ghat_prev = xmul * invtransmat(:,:,idxn+1-1,mprev)
      ghat_this = xmul * invtransmat(:,:,idxn+1-1,mthis)
      ghat_next = xmul * invtransmat(:,:,idxn+1-1,mnext)

      cprev = ghat_prev/dx(i-1) + ghat_this/dx(i)
      call inv(xslib%ngroup, cprev, matinv)
      cprev = 2.0_rk/dx(i)/dx(i-1) * matmul(matmul(ghat_this, matinv), ghat_prev)

      cnext = ghat_next/dx(i+1) + ghat_this/dx(i)
      call inv(xslib%ngroup, cnext, matinv)
      cnext = 2.0_rk/dx(i)/dx(i+1) * matmul(matmul(ghat_this, matinv), ghat_next)

      pn_prev_source(:,i) = matmul(cprev, phi_block(:,i-1,idxn+1-2)) &
        - matmul(cprev + cnext, phi_block(:,i,idxn+1-2)) &
        + matmul(cnext, phi_block(:,i+1,idxn+1-2))
    enddo ! i = 2,nx-1

    ! BC at x=L, i=N
    mprev = mat_map(nx-1)
    mthis = mat_map(nx)
    ghat_prev = xmul * invtransmat(:,:,idxn+1-1,mprev)
    ghat_this = xmul * invtransmat(:,:,idxn+1-1,mthis)
    cprev = ghat_prev/dx(nx-1) + ghat_this/dx(nx)
    call inv(xslib%ngroup, cprev, matinv)
    cprev = 2.0_rk/dx(nx-1)/dx(nx) * matmul(matmul(ghat_this, matinv), ghat_prev)
    select case (boundary_right)
      case ('mirror')
        pn_prev_source(:,nx) = matmul(cprev, phi_block(:,nx-1,idxn+1-2)) &
          - matmul(cprev, phi_block(:,nx,idxn+1-2))
      case ('zero')
        pn_prev_source(:,nx) = matmul(cprev, phi_block(:,nx-1,n-1)) &
          - matmul(cprev, phi_block(:,nx,idxn+1-2)) &
          - 2.0_rk/dx(nx) * matmul(ghat_this, phi_block(:,nx,idxn+1-2))
      case default
        call exception_fatal('unknown boundary_right in prev_source: ' // trim(adjustl(boundary_right)))
    endselect

    deallocate(matinv)
    deallocate(cnext, cprev)
    deallocate(ghat_prev, ghat_this, ghat_next)
  endsubroutine transport_block_build_prev_source

  subroutine transport_block_build_fsource(nx, dx, mat_map, xslib, phi_block, fsource)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(in) :: phi_block(:,:) ! (ngroup,nx) -- NOTE only scalar flux (zeorth moment)
    real(rk), intent(out) :: fsource(:,:) ! (ngroup,nx)

    integer(ik) :: i
    integer(ik) :: mthis

    do i = 1,nx
      mthis = mat_map(i)
      if (xslib%mat(mthis)%is_fiss) then
        fsource(:,i) = xslib%mat(mthis)%chi * sum(xslib%mat(mthis)%nusf * phi_block(:,i)) * dx(i)
      else
        fsource(:,i) = 0d0
      endif
    enddo ! i = 1,nx
  endsubroutine transport_block_build_fsource

  real(rk) function transport_block_fission_summation(nx, dx, mat_map, xslib, phi_block)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(in) :: phi_block(:,:) ! (ngroup,nx) -- NOTE only scalar flux (zeorth moment)

    integer(ik) :: i, mthis
    real(rk) :: xsum

    xsum = 0d0
    do i = 1,nx
      mthis = mat_map(i)
      if (xslib%mat(mthis)%is_fiss) then
        xsum = xsum + sum(xslib%mat(mthis)%nusf * phi_block(:,i)) * dx(i)
      endif
    enddo
    transport_block_fission_summation = xsum
  endfunction transport_block_fission_summation

  subroutine transport_block_power_iteration(&
    nx, dx, mat_map, xslib, boundary_right, k_tol, phi_tol, max_iter, pnorder, keff, sigma_tr,phi)
    use xs, only : XSLibrary
    use linalg, only : trid_block, count_inv
    use output, only : output_write
    use timer, only : timer_start, timer_stop
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    character(*), intent(in) :: boundary_right
    real(rk), intent(in) :: k_tol, phi_tol
    integer(ik), intent(in) :: max_iter
    integer(ik), intent(in) :: pnorder
    real(rk), intent(out) :: keff
    real(rk), intent(out) :: sigma_tr(:,:,:) ! (nx,ngroup,nmoment)
    real(rk), intent(out) :: phi(:,:,:) ! (nx,ngroup,nmoment)

    integer(ik) :: i, g, n
    integer(ik) :: neven, idxn
    integer(ik) :: iter

    real(rk) :: delta_k, delta_phi
    real(rk) :: fsum, fsum_old, k_old

    ! (ngroup,ngroup,nx-1,neven) , (ngroup,ngroup,nx,neven) , (ngroup,ngroup,nx-1,neven)
    real(rk), allocatable :: sub(:,:,:,:), dia(:,:,:,:), sup(:,:,:,:)
    ! (ngroup,ngroup,nx-1) , (ngroup,ngroup,nx) , (ngroup,ngroup,nx-1)
    real(rk), allocatable :: sub_copy(:,:,:), dia_copy(:,:,:), sup_copy(:,:,:)

    real(rk), allocatable :: fsource(:,:) ! (ngroup,nx)
    real(rk), allocatable :: pn_prev_source(:,:) ! (ngroup,nx) -- for this moment
    real(rk), allocatable :: pn_next_source(:,:,:) ! (ngroup,nx,neven) -- for all moments
    real(rk), allocatable :: q(:,:) ! (ngroup,nx)

    real(rk), allocatable :: phi_block(:,:,:), phi_old(:,:,:) ! (ngroup,nx,pnorder+1)

    character(1024) :: line

    if (mod(pnorder,2) /= 1) then
      call exception_fatal('pnorder must be odd')
    endif
    neven = (pnorder+1)/2

    allocate(sub(xslib%ngroup,xslib%ngroup,nx-1,neven),&
      dia(xslib%ngroup,xslib%ngroup,nx,neven),&
      sup(xslib%ngroup,xslib%ngroup,nx-1,neven))
    allocate(sub_copy(xslib%ngroup,xslib%ngroup,nx-1),&
      dia_copy(xslib%ngroup,xslib%ngroup,nx),&
      sup_copy(xslib%ngroup,xslib%ngroup,nx-1))

    allocate(fsource(xslib%ngroup,nx))
    allocate(pn_prev_source(xslib%ngroup,nx))
    allocate(pn_next_source(xslib%ngroup,nx,neven))
    allocate(q(xslib%ngroup,nx))

    allocate(phi_block(xslib%ngroup,nx,pnorder+1))
    allocate(phi_old(xslib%ngroup,nx,pnorder+1))

    phi_block(:,:,1) = 1d0
    phi_block(:,:,2:pnorder+1) = 0d0
    keff = 1d0
    fsum = 1d0

    call output_write('=== PN TRANSPORT BLOCK POWER ITERATION ===')

    call timer_start('transport_invtransmat')
    allocate(invtransmat(xslib%ngroup,xslib%ngroup,pnorder+1,xslib%niso))
    do i = 1,xslib%niso
      do n = 1,pnorder+1
        call transport_invtransmat(xslib%mat(i), n-1, invtransmat(:,:,n,i))
      enddo ! n = 1,pnorder+1
    enddo ! i = 1,xslib%niso
    call timer_stop('transport_invtransmat')

    call output_write('  building transport matrix')
    call timer_start('transport_build_matrix')
    call transport_block_build_matrix(nx, dx, mat_map, xslib, boundary_right, neven, sub, dia, sup)
    call timer_stop('transport_build_matrix')

    call output_write('  beginning iterations')
    do iter = 1,max_iter

      k_old = keff
      phi_old = phi_block
      fsum_old = fsum

      call timer_start('transport_pn_source')
      call transport_block_build_next_source(nx, dx, mat_map, xslib, boundary_right, neven, phi_block, pn_next_source)
      call timer_stop('transport_pn_source')

      do n = 1,neven

        idxn = 2*(n-1)

        q = pn_next_source(:,:,n)

        if (n == 1) then
          call timer_start('transport_fsource')
          call transport_block_build_fsource(nx, dx, mat_map, xslib, phi_block(:,:,1), fsource)
          q = q + fsource / keff
          call timer_stop('transport_fsource')
        else
          call timer_start('transport_pn_source')
          call transport_block_build_prev_source(nx, dx, mat_map, xslib, boundary_right, idxn, phi_block, pn_prev_source)
          q = q + pn_prev_source
          call timer_stop('transport_pn_source')
        endif

        call timer_start('transport_block_tridiagonal')
        sub_copy = sub(:,:,:,n)
        dia_copy = dia(:,:,:,n)
        sup_copy = sup(:,:,:,n)
        call trid_block(nx, xslib%ngroup, sub_copy, dia_copy, sup_copy, q, phi_block(:,:,idxn+1))
        call timer_stop('transport_block_tridiagonal')
      enddo ! n = 1,neven

      call timer_start('transport_convergence')
      fsum = transport_block_fission_summation(nx, dx, mat_map, xslib, phi_block(:,:,1))
      if (iter > 1) keff = keff * fsum / fsum_old
      delta_k = abs(keff - k_old)
      delta_phi = maxval(abs(phi_block - phi_old)) / maxval(phi_block)
      call timer_stop('transport_convergence')

      if ((keff < 0.0_rk) .or. (keff > 2.0_rk)) then
        write(*,*) 'invalid keff', keff
        stop
      endif

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

    write(line, '(a,i0)') 'number of dense matrix inversions: ', count_inv
    call output_write(line)

    call timer_start('calc_odd')
    call transport_block_calc_odd(nx, dx, mat_map, xslib, boundary_right, pnorder, phi_block)
    call timer_stop('calc_odd')
    call timer_start('calc_transportxs')
    call transport_block_calc_transportxs(nx, mat_map, xslib, pnorder, phi_block, sigma_tr)
    call timer_stop('calc_transportxs')
    ! copy before exit to reorder
    do i = 1,nx
      do g = 1,xslib%ngroup
        do n = 1,pnorder+1
          phi(i,g,n) = phi_block(g,i,n)
        enddo ! n = 1,neven
      enddo ! g = 1,xslib%ngroup
    enddo ! i = 1,nx

    deallocate(invtransmat)

    deallocate(fsource, pn_prev_source, pn_next_source, q)
    deallocate(phi_block, phi_old)
    deallocate(sub, dia, sup)
    deallocate(sub_copy, dia_copy, sup_copy)
  endsubroutine transport_block_power_iteration

  subroutine transport_block_calc_odd(nx, dx, mat_map, xslib, boundary_right, pnorder, phi_block)
    use xs, only : XSLibrary
    use numeric, only : deriv
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    character(*), intent(in) :: boundary_right
    integer(ik), intent(in) :: pnorder
    real(rk), intent(inout) :: phi_block(:,:,:) ! (ngroup,nx,pnorder+1)

    integer(ik) :: i, g, n
    integer(ik) :: idxn
    integer(ik) :: mthis

    real(rk) :: xn
    real(rk) :: xmul_next, xmul_prev
    real(rk) :: x1, x2, x3

    real(rk), allocatable :: dphi_prev(:), dphi_next(:) ! (ngroup)
    real(rk), allocatable :: trans(:,:) ! (ngroup,ngroup)

    allocate(dphi_prev(xslib%ngroup), dphi_next(xslib%ngroup))
    allocate(trans(xslib%ngroup,xslib%ngroup))

    do n = 2,pnorder+1,2
      idxn = n-1
      xn = real(idxn, rk)
      xmul_prev = xn/(2.0_rk*xn+1.0_rk)
      xmul_next = (xn+1.0_rk)/(2.0_rk*xn+1.0_rk)
      do i = 2,nx-1
        mthis = mat_map(i)
        x1 = -0.5_rk * (dx(i-1) + dx(i))
        x2 = 0.0_rk
        x3 = 0.5_rk * (dx(i+1) + dx(i))
        dphi_next = 0.0_rk
        if ((mat_map(i) == mat_map(i+1)) .and. (mat_map(i) == mat_map(i-1))) then
          ! central difference for interior
          do g = 1,xslib%ngroup
            dphi_prev(g) = deriv(x1, x2, x3, &
              phi_block(g,i-1,idxn+1-1), phi_block(g,i,idxn+1-1), phi_block(g,i+1,idxn+1-1))
            if (n < pnorder + 1) then
              dphi_next(g) = deriv(x1, x2, x3, &
                phi_block(g,i-1,idxn+1+1), phi_block(g,i,idxn+1+1), phi_block(g,i+1,idxn+1+1))
            endif
          enddo ! g = 1,xslib%ngroup
        elseif (mat_map(i) ==  mat_map(i+1)) then
          ! forward difference (first-order)
          dphi_prev = (phi_block(:,i+1,idxn+1-1) - phi_block(:,i,idxn+1-1))/(0.5_rk*(dx(i)+dx(i+1)))
          if (n < pnorder + 1) then
            dphi_next = (phi_block(:,i+1,idxn+1+1) - phi_block(:,i,idxn+1+1))/(0.5_rk*(dx(i)+dx(i+1)))
          endif
        elseif (mat_map(i) == mat_map(i-1)) then
          ! backward difference (first-order)
          dphi_prev = (phi_block(:,i,idxn+1-1) - phi_block(:,i-1,idxn+1-1))/(0.5_rk*(dx(i)+dx(i-1)))
          if (n < pnorder + 1) then
            dphi_next = (phi_block(:,i,idxn+1+1) - phi_block(:,i-1,idxn+1+1))/(0.5_rk*(dx(i)+dx(i-1)))
          endif
        else
          ! take a guess at central difference
          do g = 1,xslib%ngroup
            dphi_prev(g) = deriv(x1, x2, x3, &
              phi_block(g,i-1,idxn+1-1), phi_block(g,i,idxn+1-1), phi_block(g,i+1,idxn+1-1))
            if (n < pnorder + 1) then
              dphi_next(g) = deriv(x1, x2, x3, &
                phi_block(g,i-1,idxn+1+1), phi_block(g,i,idxn+1+1), phi_block(g,i+1,idxn+1+1))
            endif
          enddo ! g = 1,xslib%ngroup
        endif
        call transport_invtransmat(xslib%mat(mthis), idxn, trans)
        phi_block(:,i,idxn+1) = -matmul(trans, xmul_next * dphi_next + xmul_prev * dphi_prev)
      enddo ! i = 2,nx-1
      ! BC at x=0, i=1
      ! use the fact that odd moments must equal zero for mirror bc
      ! this stencil kind of extends to x3 because phi(2) was computed earlier
      phi_block(:,1,idxn+1) = phi_block(:,2,idxn+1) &
        * 0.5_rk * dx(1) / (dx(1) + 0.5_rk*dx(2))
      ! BC at x=L, i=N
      select case (boundary_right)
        case ('mirror')
          phi_block(:,nx,idxn+1) = phi_block(:,nx-1,idxn+1) &
            * 0.5_rk * dx(nx) / (dx(nx) + 0.5_rk*dx(nx-1))
        case ('zero')
          dphi_prev = -phi_block(:,nx-1,idxn+1-1)/(dx(nx) + 0.5_rk*(dx(nx-1)))
          if (n < pnorder + 1) then
            dphi_next = -phi_block(:,nx-1,idxn+1+1)/(dx(nx) + 0.5_rk*(dx(nx-1)))
          endif
          call transport_invtransmat(xslib%mat(mat_map(nx)), idxn, trans)
          phi_block(:,nx,idxn+1) = -matmul(trans, xmul_next * dphi_next + xmul_prev * dphi_prev)
        case default
          call exception_fatal('unknown boundary in odd_update: ' // trim(adjustl(boundary_right)))
      endselect
    enddo

    deallocate(trans)
    deallocate(dphi_prev, dphi_next)
  endsubroutine transport_block_calc_odd

  subroutine transport_block_calc_transportxs(nx, mat_map, xslib, pnorder, phi_block, sigma_tr)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    integer(ik), intent(in) :: mat_map(:) ! (dx)
    type(XSLibrary), intent(in) :: xslib
    integer(ik), intent(in) :: pnorder
    real(rk), intent(in) :: phi_block(:,:,:) ! (ngroup,nx,pnorder+1)

    ! NOTE the ordeirng is different here!
    ! This was done becasue I'm lazy...
    real(rk), intent(out) :: sigma_tr(:,:,:) ! (nx,ngroup,pnorder+1)

    integer(ik) :: i, g, n
    integer(ik) :: mthis
    
    do n = 1,pnorder+1
      if (n < xslib%nmoment) then
        do i = 1,nx
          mthis = mat_map(i)
          do g = 1,xslib%ngroup
            sigma_tr(i,g,n) = xslib%mat(mthis)%sigma_t(g) &
              - sum(xslib%mat(mthis)%scatter(:,g,n)*phi_block(:,i,n))/phi_block(g,i,n)
          enddo ! g = 1,xslib%ngroup
        enddo ! i = 1,nx
      else
        do i = 1,nx
          mthis = mat_map(i)
          sigma_tr(i,:,n) = xslib%mat(mthis)%sigma_t(:)
        enddo ! i = 1,nx
      endif
    enddo ! n = 1,pnorder+1
  endsubroutine transport_block_calc_transportxs

endmodule transport_block
