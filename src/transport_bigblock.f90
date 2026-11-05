module transport_bigblock
use kind, only : rk, ik
implicit none

private
public :: transport_bigblock_power_iteration

real(rk), allocatable :: invtransmat(:,:,:,:) ! (ngroup,ngroup,pnorder+1,xslib%niso)

contains

  subroutine transport_bigblock_power_iteration( &
      nx, dx, mat_map, xslib, boundary_left, boundary_right, albedo_coeff, &
      calc_type, k_tol, phi_tol, max_iter, pnorder, keff, sigma_tr, phi)
    use xs, only : XSLibrary
    use exception_handler, only : exception_fatal
    use output, only : output_write
    use transport_block, only : transport_invtransmat
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
    integer(ik) :: i, g
    integer(ik) :: idx, idxn

    real(rk) :: fsum

    real(rk), allocatable :: sub(:,:) ! (dimlen,nx-1)
    real(rk), allocatable :: dia(:,:) ! (dimlen,nx)
    real(rk), allocatable :: sup(:,:) ! (dimlen,nx-1)
    real(rk), allocatable :: q(:,:) ! (dimlen,nx)
    real(rk), allocatable :: phi_bigblock(:,:), phi_bigblock_old(:,:) ! (dimlen,nx)

    if (mod(pnorder,2) /= 1) then
      call exception_fatal('pnorder must be odd')
    endif
    neven = (pnorder+1)/2

    dimlen = neven * xslib%ngroup
    allocate(&
      sub(dimlen,nx-1), &
      dia(dimlen,nx), &
      sup(dimlen,nx-1) &
    )
    allocate(q(dimlen,nx))

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
    do i =1,nx
      do g = 1,xslib%ngroup
        phi_bigblock(g,i) = 1.0_rk
      enddo ! g = 1,xslib%ngroup
    enddo ! i = 1,nx

    keff = 1.0_rk
    fsum = 1.0_rk

    call output_write('=== PN TRANSPORT BIG-BLOCK POWER ITERATION ===')

    if (boundary_left /= 'mirror') then
      call exception_fatal('In bigblock, boundary_left must be mirror.')
    elseif (boundary_right /= 'mirror') then
      call exception_fatal('In bigblock, boundary_right must be mirror.')
    elseif (calc_type /= 'eigenvalue') then
      call exception_fatal('In bigblock, calc_type must be eigenvalue.')
    endif

    deallocate(sub, dia, sup)
    deallocate(q)
  endsubroutine transport_bigblock_power_iteration

endmodule transport_bigblock
