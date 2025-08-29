module transport
use kind
implicit none

private

public :: sigma_tr, transport_power_iteration

real(rk), allocatable :: sigma_tr(:,:,:) ! (nx, ngroup, nmoment)

contains

  subroutine transport_build_matrix(nx, dx, mat_map, xslib, boundary_right, neven, sub, dia, sup)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    character(*), intent(in) :: boundary_right
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
      xmul_next = (xn+1.0_rk)**2 / ((2.0_rk*xn+1.0_rk)*(2.0_rk*xn+3.0_rk))
      xmul_prev = xn**2 / (4.0_rk*xn**2 - 1.0_rk)
      do g = 1,xslib%ngroup
        cthis = xmul_next / sigma_tr(1,g,idxn+1+1)
        cnext = xmul_next / sigma_tr(2,g,idxn+1+1)
        if (idxn > 1) then
          cthis = cthis + xmul_prev / sigma_tr(1,g,idxn+1-1)
          cnext = cnext + xmul_prev / sigma_tr(2,g,idxn+1-1)
        endif
        dnext = 2 * cthis / dx(1) * cnext / dx(2) / (cthis / dx(1) + cnext / dx(2))
        dia(1,g,n) = + dnext + xslib%mat(mthis)%sigma_t(g) * dx(1)
        if (idxn+1 <= xslib%nmoment) then
          dia(1,g,n) = dia(1,g,n) - xslib%mat(mthis)%scatter(g,g,idxn+1) * dx(1)
        endif
        sup(1,g,n) = -dnext
      enddo ! = g = 1,xslib%ngroup
    enddo ! n = 1,neven

    do n = 1,neven
      idxn = 2*(n-1)
      xn = real(idxn, rk)
      xmul_next = (xn+1.0_rk)**2 / ((2.0_rk*xn+1.0_rk)*(2.0_rk*xn+3.0_rk))
      xmul_prev = xn**2 / (4.0_rk*xn**2 - 1.0_rk)
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

          dprev = 2 * cthis / dx(i) * cprev / dx(i-1) / (cthis / dx(i) + cprev / dx(i-1))
          dnext = 2 * cthis / dx(i) * cnext / dx(i+1) / (cthis / dx(i) + cnext / dx(i+1))

          sub(i-1,g,n) = -dprev
          dia(i,g,n) = dprev + dnext + xslib%mat(mthis)%sigma_t(g) * dx(i)
          if (idxn+1 <= xslib%nmoment) then
            dia(i,g,n) = dia(i,g,n) - xslib%mat(mthis)%scatter(g,g,idxn+1) * dx(i)
          endif
          sup(i,g,n) = -dnext

        enddo ! i = 2,nx-1
      enddo ! g = 1,xslib%ngroup
    enddo ! n = 1,neven

    ! BC at x=L, i=N
    mprev = mat_map(nx-1)
    mthis = mat_map(nx)
    select case (boundary_right)
      case ('mirror')
        do n = 1,neven
          idxn = 2*(n-1)
          xn = real(idxn, rk)
          xmul_next = (xn+1.0_rk)**2 / ((2.0_rk*xn+1.0_rk)*(2.0_rk*xn+3.0_rk))
          xmul_prev = xn**2 / (4.0_rk*xn**2 - 1.0_rk)
          do g = 1,xslib%ngroup
            cprev = xmul_next / sigma_tr(nx-1,g,idxn+1+1)
            cthis = xmul_next / sigma_tr(nx  ,g,idxn+1+1)
            if (idxn > 0) then
              cprev = cprev + xmul_prev / sigma_tr(nx-1,g,idxn+1-1)
              cthis = cthis + xmul_prev / sigma_tr(nx  ,g,idxn+1-1)
            endif
            dprev = 2 * cthis / dx(nx) * cprev / dx(nx-1) / (cthis / dx(nx) + cprev / dx(nx-1))
            sub(nx-1,g,n) = -dprev
            dia(nx,g,n) = dprev + xslib%mat(mthis)%sigma_t(g) * dx(nx)
            if (idxn+1 <= xslib%nmoment) then
              dia(nx,g,n) = dia(nx,g,n) - xslib%mat(mthis)%scatter(g,g,idxn+1) * dx(nx)
            endif
          enddo ! g = 1,xslib%ngroup
        enddo ! n = 1,neven
      case ('zero')
        do n = 1,neven
          idxn = 2*(n-1)
          xn = real(idxn, rk)
          xmul_next = (xn+1.0_rk)**2 / ((2.0_rk*xn+1.0_rk)*(2.0_rk*xn+3.0_rk))
          xmul_prev = xn**2 / (4.0_rk*xn**2 - 1.0_rk)
          do g = 1,xslib%ngroup
            cprev = xmul_next / sigma_tr(nx-1,g,idxn+1+1)
            cthis = xmul_next / sigma_tr(nx  ,g,idxn+1+1)
            if (idxn > 0) then
              cprev = cprev + xmul_prev / sigma_tr(nx-1,g,idxn+1-1)
              cthis = cthis + xmul_prev / sigma_tr(nx  ,g,idxn+1-1)
            endif
            dprev = 2 * cthis / dx(nx) * cprev / dx(nx-1) / (cthis / dx(nx) + cprev / dx(nx-1))
            sub(nx-1,g,n) = -dprev
            dia(nx,g,n) = dprev + xslib%mat(mthis)%sigma_t(g) * dx(nx) &
              + 2 * cthis / dx(nx)
            if (idxn+1 <= xslib%nmoment) then
              dia(nx,g,n) = dia(nx,g,n) - xslib%mat(mthis)%scatter(g,g,idxn+1) * dx(nx)
            endif
          enddo ! g = 1,xslib%ngroup
        enddo ! n = 1,neven
      case default
        write(*,*) 'unknown boundary_right in transport_build_matrix: ' // trim(adjustl(boundary_right))
        stop
    endselect


  endsubroutine transport_build_matrix

  subroutine transport_init_transportxs(nx, mat_map, xslib, pnorder, sigma_tr)
    use xs, only : XSLibrary
    integer, intent(in) :: nx
    integer, intent(in) :: mat_map(:)
    type(XSLibrary), intent(in) :: xslib
    integer, intent(in) :: pnorder
    real(rk), allocatable, intent(out) :: sigma_tr(:,:,:) ! (nx,ngroup,pnorder+1)

    integer(ik) :: i, g, n
    integer(ik) :: mthis

    allocate(sigma_tr(nx,xslib%ngroup,pnorder+1))

    do n = 1,pnorder+1
      if (n <= xslib%nmoment) then
        do g = 1,xslib%ngroup
          do i = 1,nx
            mthis = mat_map(i)
            sigma_tr(i,g,n) = xslib%mat(mthis)%sigma_t(g) - xslib%mat(mthis)%scatter(g,g,n)
          enddo ! i = 1,nx
        enddo ! g = 1,xslib%ngroup
      else
        do g = 1,xslib%ngroup
          do i = 1,nx
            mthis = mat_map(i)
            sigma_tr(i,g,n) = xslib%mat(mthis)%sigma_t(g)
          enddo ! i = 1,nx
        enddo ! g = 1,xslib%ngroup
      endif
    enddo ! n = 1,pnorder+1

  endsubroutine transport_init_transportxs

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

    do n = 1,pnorder+1
      if (n <= xslib%nmoment) then
        if (mod(n,2) == 0) then
          ! If n is even, idxn (the actual PN order) is odd... off-by-one
          do i = 1,nx
            mthis = mat_map(i)
            do g = 1,xslib%ngroup
              sigma_tr(i,g,n) = xslib%mat(mthis)%sigma_t(g) - sum(xslib%mat(mthis)%scatter(:,g,n)*phi(i,:,n))/phi(i,g,n)
            enddo ! g = 1,xslib%ngroup
          enddo ! i = 1,nx
        else
          ! I don't really care about the even transport xs.
          ! They are not used in the calculations anywhere
          do i = 1,nx
            mthis = mat_map(i)
            do g = 1,xslib%ngroup
              sigma_tr(i,g,n) = xslib%mat(mthis)%sigma_t(g) - sum(xslib%mat(mthis)%scatter(:,g,n)*phi(i,:,n))/phi(i,g,n)
            enddo ! xslib%ngroup
          enddo ! i = 1,nx
        endif
      else
        do i = 1,nx
          mthis = mat_map(i)
          do g = 1,xslib%ngroup
            sigma_tr(i,g,n) = xslib%mat(mthis)%sigma_t(g)
          enddo ! g = 1,xslib%ngroup
        enddo ! n = 1,nx
      endif
    enddo ! n = 1,pnorder+1
  endsubroutine transport_build_transportxs

  subroutine transport_odd_update(nx, dx, ng, boundary_right, pnorder, sigma_tr, phi)
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: ng
    character(*), intent(in) :: boundary_right
    integer(ik), intent(in) :: pnorder
    real(rk), intent(in) :: sigma_tr(:,:,:) ! (nx, ngroup, pnorder)
    real(rk), intent(inout) :: phi(:,:,:) ! (nx, ngroup, pnorder)

    integer(ik) :: n, g, i
    integer(ik) :: idxn
    real(rk) :: xn, xmul_prev, xmul_next
    real(rk) :: dphi_prev, dphi_next

    do n = 2,pnorder+1,2
      idxn = n-1
      xn = real(idxn, rk)
      xmul_prev = xn/(2.0_rk*xn+1.0_rk)
      xmul_next = (xn+1.0_rk)/(2.0_rk*xn+1.0_rk)
      if (n < pnorder+1) then
        do g = 1,ng
          do i = 2,nx-1
            ! central difference for interior
            dphi_prev = deriv(-0.5_rk*(dx(i-1)+dx(i)), 0.0_rk, +0.5_rk*(dx(i)+dx(i+1)), &
              phi(i-1,g,idxn+1-1), phi(i,g,idxn+1-1), phi(i+1,g,idxn+1-1))
            dphi_next = deriv(-0.5_rk*(dx(i-1)+dx(i)), 0.0_rk, +0.5_rk*(dx(i)+dx(i+1)), &
              phi(i-1,g,idxn+1+1), phi(i,g,idxn+1+1), phi(i+1,g,idxn+1+1))
            phi(i,g,idxn+1) = &
              - (xmul_prev * dphi_prev + xmul_next * dphi_next) &
              / sigma_tr(i,g,idxn+1)
          enddo ! i = 2,nx-1
          ! BC at x=0, i=1
          ! use the fact that odd moments must equal zero for mirror bc
          ! this stencil kind of extends to x3 because phi(2) was computed earlier
          phi(1,g,idxn+1) = phi(2,g,idxn+1) * 0.5_rk * dx(1) / (dx(1) + 0.5_rk*dx(2))
          ! BC at x=L, i=N
          select case (boundary_right)
            case ('mirror')
              phi(nx,g,idxn+1) = phi(nx-1,g,idxn+1) * 0.5_rk*dx(nx) / (dx(nx)+0.5_rk*dx(nx-1))
            case ('zero')
              dphi_prev = -phi(nx-1,g,idxn+1-1)/(dx(nx) + 0.5_rk*dx(nx-1))
              dphi_next = -phi(nx-1,g,idxn+1+1)/(dx(nx) + 0.5_rk*dx(nx-1))
              phi(nx,g,idxn+1) = &
                - (xmul_prev * dphi_prev + xmul_next * dphi_next) &
                / sigma_tr(nx,g,idxn+1)
            case default
              write(*,*) 'unknown boundary in odd_update: ' // trim(adjustl(boundary_right))
              stop
          endselect
        enddo ! g = 1,ng
      else
        do g = 1,ng
          do i = 2,nx-1
            ! central difference for interior
            dphi_prev = deriv(-0.5_rk*(dx(i-1)+dx(i)), 0.0_rk, +0.5_rk*(dx(i)+dx(i+1)), &
              phi(i-1,g,idxn+1-1), phi(i,g,idxn+1-1), phi(i+1,g,idxn+1-1))
            phi(i,g,idxn+1) = - xmul_prev * dphi_prev / sigma_tr(i,g,idxn+1)
          enddo ! i = 2,nx-1
          ! BC at x=0, i=1
          ! use the fact that odd moments must equal zero for mirror bc
          ! this stencil kind of extends to x3 because phi(2) was computed earlier
          phi(1,g,idxn+1) = phi(2,g,idxn+1) * 0.5_rk * dx(1) / (dx(1) + 0.5_rk*dx(2))
          ! BC at x=L, i=N
          select case (boundary_right)
            case ('mirror')
              phi(nx,g,idxn+1) = phi(nx-1,g,idxn+1) * 0.5_rk*dx(nx) / (dx(nx)+0.5_rk*dx(nx-1))
            case ('zero')
              dphi_prev = -phi(nx-1,g,idxn+1-1)/(dx(nx) + 0.5_rk*dx(nx-1))
              phi(nx,g,idxn+1) = - xmul_prev * dphi_prev / sigma_tr(nx,g,idxn+1)
            case default
              write(*,*) 'unknown boundary2 in odd_update: ' // trim(adjustl(boundary_right))
              stop
          endselect
        enddo ! g = 1,ng
      endif
    enddo ! n = 2,pnorder+1,2
  endsubroutine transport_odd_update

  subroutine transport_build_upscatter(nx, dx, mat_map, xslib, phi, idxn, qup)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(in) :: phi(:,:,:) ! (nx, ngroup, pnorder)
    integer(ik), intent(in) :: idxn
    real(rk), intent(out) :: qup(:,:) ! (nx, ngroup)

    integer(ik) :: i, g
    integer(ik) :: mthis

    if (idxn+1 > xslib%nmoment) then
      qup = 0.0_rk
    else
      do i = 1,nx
        mthis = mat_map(i)
        do g = 1,xslib%ngroup
          qup(i,g) = &
            sum(xslib%mat(mthis)%scatter(g+1:xslib%ngroup,g,idxn+1) * phi(i,g+1:xslib%ngroup,idxn+1)) * dx(i)
        enddo
      enddo ! i = 1,nx
    endif
  endsubroutine transport_build_upscatter

  subroutine transport_build_downscatter(nx, dx, mat_map, xslib, phi, idxn, g, qdown)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(in) :: phi(:,:,:) ! (nx, ngroup, pnorder)
    integer(ik), intent(in) :: idxn
    integer(ik), intent(in) :: g
    real(rk), intent(out) :: qdown(:) ! (nx)

    integer(ik) :: i
    integer(ik) :: mthis

    if (idxn+1 > xslib%nmoment) then
      qdown = 0.0_rk
    else
      do i = 1,nx
        mthis = mat_map(i)
        qdown(i) = &
          sum(xslib%mat(mthis)%scatter(1:g-1,g,idxn+1) * phi(i,1:g-1,idxn+1)) * dx(i)
      enddo ! i = 1,nx
    endif
  endsubroutine transport_build_downscatter

  ! NOTE: this is an exact copy of diffusion_build_fsource ...
  subroutine transport_build_fsource(nx, dx, mat_map, xslib, phi0, qfiss)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(in) :: phi0(:,:) ! (nx, ngroup)
    real(rk), intent(out) :: qfiss(:,:) ! (nx, ngroup)

    integer(ik) :: i
    integer(ik) :: mthis

    do i = 1,nx
      mthis = mat_map(i)
      if (xslib%mat(mthis)%is_fiss) then
        qfiss(i,:) = xslib%mat(mthis)%chi(:) * sum(xslib%mat(mthis)%nusf(:) * phi0(i,:)) * dx(i)
      else
        qfiss(i,:) = 0.0_rk
      endif
    enddo ! i = 1,nx

  endsubroutine transport_build_fsource

  subroutine transport_build_next_source(nx, dx, ngroup, boundary_right, neven, sigma_tr, phi, qnext)
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: ngroup
    character(*), intent(in) :: boundary_right
    integer(ik), intent(in) :: neven
    real(rk), intent(in) :: sigma_tr(:,:,:) ! (nx, ngroup, pnorder)
    real(rk), intent(in) :: phi(:,:,:) ! (nx, ngroup, pnorder)
    real(rk), intent(out) :: qnext(:,:,:) ! (nx, ngroup, neven)

    integer(ik) :: i, g, n
    integer(ik) :: idxn
    real(rk) :: xn, xmul
    real(rk) :: kaprev, kathis, kanext
    real(rk) :: daprev, danext

    qnext(:,:,neven) = 0.0_rk
    do n = 1,neven-1
      idxn = 2*(n-1)
      xn = real(idxn, rk)
      xmul = (xn+1.0_rk)*(xn+2.0_rk)/((2.0_rk*xn+1.0_rk)*(2.0_rk*xn+3.0_rk))
      do g = 1,ngroup

        ! BC at x=0, i=1
        kathis = xmul/sigma_tr(1,g,idxn+1+1)
        kanext = xmul/sigma_tr(2,g,idxn+1+1)
        danext = 2 * kathis / dx(1) * kanext / dx(2) / (kathis / dx(1) + kanext / dx(2))
        qnext(1,g,n) = -phi(1,g,idxn+1+2)*danext + danext*phi(2,g,idxn+1+2)

        do i = 2,nx-1

          kaprev = xmul/sigma_tr(i-1,g,idxn+1+1)
          kathis = xmul/sigma_tr(i  ,g,idxn+1+1)
          kanext = xmul/sigma_tr(i+1,g,idxn+1+1)

          daprev = 2 * kathis / dx(i) * kaprev / dx(i-1) / (kathis / dx(i) + kaprev / dx(i-1))
          danext = 2 * kathis / dx(i) * kanext / dx(i+1) / (kathis / dx(i) + kanext / dx(i+1))

          qnext(i, g, n) = phi(i-1,g,idxn+1+2) * daprev - phi(i,g,idxn+1+2) * (daprev + danext) + phi(i+1,g,idxn+1+2) * danext

        enddo

        ! BC at x=L, i=N
        select case (boundary_right)
          case ('mirror')
            kaprev = xmul/sigma_tr(nx-1,g,idxn+1+1)
            kathis = xmul/sigma_tr(nx  ,g,idxn+1+1)
            daprev = 2 * kathis / dx(nx) * kaprev / dx(nx-1) / (kathis / dx(nx) + kaprev / dx(nx-1))
            qnext(nx,g,n) = phi(nx-1,g,idxn+1+2)*daprev - phi(nx,g,idxn+1+2)*daprev
          case ('zero')
            kaprev = xmul/sigma_tr(nx-1,g,idxn+1+1)
            kathis = xmul/sigma_tr(nx  ,g,idxn+1+1)
            daprev = 2 * kathis / dx(nx) * kaprev / dx(nx-1) / (kathis / dx(nx) + kaprev / dx(nx-1))
            qnext(nx,g,n) = phi(nx-1,g,idxn+1+2)*daprev - phi(nx,g,idxn+1+2)*3*daprev
          case default
            write(*,*) 'unknown boundary in next_source: ' // trim(adjustl(boundary_right))
            stop
        endselect

      enddo ! g = 1,ngroup
    enddo ! n = 1,neven

  endsubroutine transport_build_next_source

  subroutine transport_build_prev_source(nx, dx, ngroup, boundary_right, sigma_tr, phi, idxn, qprev)
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: ngroup
    character(*), intent(in) :: boundary_right
    real(rk), intent(in) :: sigma_tr(:,:,:) ! (nx, ngroup, pnorder)
    real(rk), intent(in) :: phi(:,:,:) ! (nx, ngroup, pnorder)
    integer(ik), intent(in) :: idxn
    real(rk), intent(out) :: qprev(:,:) ! (nx, ngroup)

    integer(ik) :: g, i
    real(rk) :: xn, xmul
    real(rk) :: kbprev, kbthis, kbnext
    real(rk) :: dbprev, dbnext

    xn = real(idxn, rk)
    xmul = (xn**2-xn)/(4.0_rk*xn**2 - 1.0_rk)

    do g = 1,ngroup

      ! BC at x=0, i=1
      kbthis = xmul/sigma_tr(1,g,idxn+1-1)
      kbnext = xmul/sigma_tr(2,g,idxn+1-1)
      dbnext = 2 * kbthis / dx(1) * kbnext / dx(2) / (kbthis / dx(1) + kbnext / dx(2))
      qprev(1,g) = -phi(1,g,idxn+1-2)*dbnext + dbnext*phi(2,g,idxn+1-2)

      do i = 2,nx-1

        kbprev = xmul/sigma_tr(i-1,g,idxn+1-1)
        kbthis = xmul/sigma_tr(i  ,g,idxn+1-1)
        kbnext = xmul/sigma_tr(i+1,g,idxn+1-1)

        dbprev = 2 * kbthis / dx(i) * kbprev / dx(i-1) / (kbthis / dx(i) + kbprev / dx(i-1))
        dbnext = 2 * kbthis / dx(i) * kbnext / dx(i+1) / (kbthis / dx(i) + kbnext / dx(i+1))

        qprev(i,g) = phi(i-1,g,idxn+1-2) * dbprev - phi(i,g,idxn+1-2) * (dbprev + dbnext) + phi(i+1,g,idxn+1-2) * dbnext

      enddo ! i = 2,nx-1

      ! BC at x=L, i=N
      select case (boundary_right)
        case ('mirror')
          kbprev = xmul/sigma_tr(nx-1,g,idxn+1-1)
          kbthis = xmul/sigma_tr(nx  ,g,idxn+1-1)
          dbprev = 2 * kbthis / dx(nx) * kbprev / dx(nx-1) / (kbthis / dx(nx) + kbprev / dx(nx-1))
          qprev(nx,g) = phi(nx-1,g,idxn+1-2)*dbprev - phi(nx,g,idxn+1-2)*dbprev
        case ('zero')
          kbprev = xmul/sigma_tr(nx-1,g,idxn+1-1)
          kbthis = xmul/sigma_tr(nx  ,g,idxn+1-1)
          dbprev = 2 * kbthis / dx(nx) * kbprev / dx(nx-1) / (kbthis / dx(nx) + kbprev / dx(nx-1))
          qprev(nx,g) = phi(nx-1,g,idxn+1-2)*dbprev - phi(nx,g,idxn+1-2)*3*dbprev
        case default
          write(*,*) 'unknown boundary in prev_source: ' // trim(adjustl(boundary_right))
          stop
      endselect

    enddo ! g = 1,ngroup

  endsubroutine transport_build_prev_source

  ! NOTE: this is an exact copy of diffusion_fission_sumation ...
  real(rk) function transport_fission_summation(nx, dx, mat_map, xslib, phi0)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx) type(XSLibrary), intent(in) :: xslib
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(in) :: phi0(:,:) ! (nx, ngroup) -- scalar flux

    integer(ik) :: i, mthis
    real(rk) :: xsum

    xsum = 0.0_rk
    do i = 1,nx
      mthis = mat_map(i)
      if (xslib%mat(mthis)%is_fiss) then
        xsum = xsum + sum(xslib%mat(mthis)%nusf(:) * phi0(i,:)) * dx(i)
      endif
    enddo ! i = 1,nx
    transport_fission_summation = xsum

  endfunction transport_fission_summation

  subroutine transport_power_iteration(nx, dx, mat_map, xslib, boundary_right, k_tol, phi_tol, max_iter, pnorder, keff, phi)
    use xs, only : XSLibrary
    use linalg, only : trid
    use output, only : output_write
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    character(*), intent(in) :: boundary_right
    real(rk), intent(in) :: k_tol, phi_tol
    integer(ik), intent(in) :: max_iter
    integer(ik), intent(in) :: pnorder
    real(rk), intent(inout) :: keff
    real(rk), intent(inout) :: phi(:,:,:) ! (nx, ngroup, nmoment)

    ! matrix
    real(rk), allocatable :: sub(:,:,:), dia(:,:,:), sup(:,:,:) ! (nx, ngroup, neven)
    real(rk), allocatable :: sub_copy(:), dia_copy(:), sup_copy(:)
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
    real(rk), allocatable :: phi_old(:,:,:)
    real(rk) :: delta_k, delta_phi

    integer(ik) :: neven, idxn
    integer(ik) :: n, g

    character(1024) :: line

    if (mod(pnorder,2) /= 1) then
      stop 'pnorder must be odd'
    endif

    neven = max((pnorder + 1) / 2, 1)

    allocate(sub(nx-1,xslib%ngroup,neven), dia(nx,xslib%ngroup,neven), sup(nx-1,xslib%ngroup,neven))
    allocate(sub_copy(nx-1), dia_copy(nx), sup_copy(nx-1))

    allocate(fsource(nx,xslib%ngroup))
    allocate(upsource(nx,xslib%ngroup))
    allocate(downsource(nx))
    allocate(pn_next_source(nx,xslib%ngroup,neven))
    allocate(pn_prev_source(nx,xslib%ngroup))
    allocate(q(nx))

    allocate(phi_old(nx,xslib%ngroup,pnorder+1))

    keff = 1.0_rk
    phi = 1.0_rk
    fsum = 1.0_rk

    if (.not. allocated(sigma_tr)) then
      call transport_init_transportxs(nx, mat_map, xslib, pnorder, sigma_tr)
    endif

    call output_write('=== PN TRANSPORT POWER ITERATION ===')

    ! TODO I believe that the original FLIP algorithm had the nesting in the opposite order.
    ! I think I got away with it in LUPINE since it is pretty irrelevant for isotropic scattering.
    ! This could require somewhat intensive rewriting of routines.
    ! See Gelbard (1959) and Fletcher (1983).
    ! Particularly Fletcher, p. 36, bottom of column 1 and top of column 2.

    do iter = 1,max_iter
      k_old = keff
      phi_old = phi
      fsum_old = fsum

      ! update odd moments -> use odd moments for transport xs -> reconstruct transport matrices
      call transport_odd_update(nx, dx, xslib%ngroup, boundary_right, pnorder, sigma_tr, phi)
      call transport_build_transportxs(nx, mat_map, xslib, pnorder, phi, sigma_tr)
      call transport_build_matrix(nx, dx, mat_map, xslib, boundary_right, neven, sub, dia, sup)

      call transport_build_next_source(nx, dx, xslib%ngroup, boundary_right, neven, sigma_tr, phi, pn_next_source)

      do n = 1,neven

        idxn = 2*(n-1)

        if (n == 1) then
          call transport_build_fsource(nx, dx, mat_map, xslib, phi(:,:,1), fsource)
        else ! (n > 1)
          call transport_build_prev_source(nx, dx, xslib%ngroup, boundary_right, sigma_tr, phi, idxn, pn_prev_source)
        endif

        call transport_build_upscatter(nx, dx, mat_map, xslib, phi, idxn, upsource)

        do g = 1,xslib%ngroup

          call transport_build_downscatter(nx, dx, mat_map, xslib, phi, idxn, g, downsource)

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
          sub_copy = sub(:,g,n)
          dia_copy = dia(:,g,n)
          sup_copy = sup(:,g,n)
          call trid(nx, sub_copy, dia_copy, sup_copy, q, phi(:,g,idxn+1))

        enddo ! g = 1,ngroup
      enddo ! n = 1,neven

      ! eigenvalue update
      fsum = transport_fission_summation(nx, dx, mat_map, xslib, phi(:,:,1))
      if (iter > 1) keff = keff * fsum / fsum_old
      delta_k = abs(keff - k_old)
      ! NOTE: we look at convergence in all space, all groups, all moments!
      ! This is seriously overkill, but necessary to demonstrate that the odd moments are second-order convergent as well.
      ! Furthermore, we may expect that the higher-order moments are important for anisotropic scattering.
      delta_phi = maxval(abs(phi - phi_old)) / maxval(phi)

      if ((keff < 0.0_rk) .or. (keff > 2.0_rk)) then
        write(*,*) 'invalid keff', keff
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
      call output_write('WARNING: failed to converge')
    endif

    deallocate(sub, dia, sup)
    deallocate(sub_copy, dia_copy, sup_copy)
    deallocate(fsource, upsource, downsource, q)
    deallocate(pn_next_source, pn_prev_source)
    deallocate(phi_old)

  endsubroutine transport_power_iteration

  ! return first derivative using a second-order estimate on non-uniform grid
  ! derivative is returned at coordinate x2
  ! this function also works in the case of a non-uniform grid
  !
  ! f1    f2    f3
  ! |-----|-----|
  ! x1    x2    x3
  ! h_left|h_right
  !
  real(rk) pure function deriv(x1, x2, x3, f1, f2, f3)
    real(rk), intent(in) :: x1, x2, x3 ! x-coordinate
    real(rk), intent(in) :: f1, f2, f3 ! function value

    real(rk) :: h_left, h_right
    real(rk) :: a, b, c

    h_left = x2 - x1
    h_right = x3 - x2

    a = -h_right**2  / (h_left * h_right * (h_left + h_right))
    c = h_left**2  / (h_left * h_right * (h_left + h_right))
    b = -(a+c)

    deriv = a * f1 + b * f2 + c * f3
  endfunction deriv

endmodule transport
