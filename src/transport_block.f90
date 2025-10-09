module transport_block
use kind
use exception_handler
implicit none

private
public :: transport_block_power_iteration

contains

  subroutine transport_block_build_matrix(nx, dx, mat_map, xslib, boundary_right, neven, sub, dia, sup)
    use xs, only : XSLibrary
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    character(*), intent(in) :: boundary_right
    integer(ik), intent(in) :: neven
    real(rk), intent(out) :: sub(:,:,:,:) ! (ngroup,ngroup,nx-1,neven)
    real(rk), intent(out) :: dia(:,:,:,:) ! (ngroup,ngroup,nx,neven)
    real(rk), intent(out) :: sup(:,:,:,:) ! (ngroup,ngroup,nx-1,neven)
  endsubroutine transport_block_build_matrix

  subroutine transport_block_power_iteration(nx, dx, mat_map, xslib, boundary_right, k_tol, phi_tol, max_iter, pnorder, keff, phi)
    use xs, only : XSLibrary
    use linalg, only : trid_block
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
    real(rk), intent(out) :: phi(:,:,:) ! (nx,ngroup,nmoment)

    integer(ik) :: neven
    integer(ik) :: i, g, n

    real(rk) :: fsum, fsum_old, k_old

    ! (ngroup,ngroup,nx-1,neven) , (ngroup,ngroup,nx,neven) , (ngroup,ngroup,nx-1,neven)
    real(rk), allocatable :: sub(:,:,:,:), dia(:,:,:,:), sup(:,:,:,:)
    ! (ngroup,ngroup,nx-1) , (ngroup,ngroup,nx) , (ngroup,ngroup,nx-1)
    real(rk), allocatable :: sub_copy(:,:,:), dia_copy(:,:,:), sup_copy(:,:,:)

    real(rk), allocatable :: fsource(:,:) ! (ngroup,nx)
    real(rk), allocatable :: pn_prev_source(:,:) ! (ngroup,nx) -- for this moment
    real(rk), allocatable :: pn_next_source(:,:,:) ! (ngroup,nx,neven) -- for all moments
    real(rk), allocatable :: q(:,:) ! (ngroup,nx)

    real(rk), allocatable :: phi_block(:,:,:), phi_old(:,:,:) ! (ngroup,nx,nmoment)

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

    allocate(phi_block(xslib%ngroup,nx,neven))
    allocate(phi_old(xslib%ngroup,nx,neven))

    phi_block(:,:,1) = 1d0
    phi_block(:,:,2:neven) = 0d0

    keff = 1d0
    fsum = 1d0

    call output_write('=== PN TRANSPORT BLOCK POWER ITERATION ===')

    call timer_start('transport_build_matrix')
    call transport_block_build_matrix(nx, dx, mat_map, xslib, boundary_right, neven, sub, dia, sup)
    call timer_stop('transport_build_matrix')

    ! copy before exit
    do i = 1,nx
      do g = 1,xslib%ngroup
        do n = 1,neven
          phi(i,g,2*(n-1)+1) = phi_block(g,i,n)
        enddo ! n = 1,neven
      enddo ! g = 1,xslib%ngroup
    enddo ! i = 1,nx
    ! TODO need to compute odd moments and transport xs as an edit

    deallocate(fsource, pn_prev_source, pn_next_source, q)
    deallocate(phi_block, phi_old)
    deallocate(sub, dia, sup)
    deallocate(sub_copy, dia_copy, sup_copy)
  endsubroutine transport_block_power_iteration

endmodule transport_block
