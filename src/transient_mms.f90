module transient_mms
use kind, only : rk, ik
use constant, only : pi
use delayed_neutron_data, only : DelayedNeutronData
use xs, only : XSLibrary
implicit none

private

real(rk), parameter :: mms_L = 100.0_rk ! [cm] -- slab_length
real(rk), parameter :: mms_p0 = 1e-3_rk ! [W/cm^2] -- initial power
real(rk), parameter :: mms_coeff_a = -1184.0_rk/15.0_rk /mms_L**4
real(rk), parameter :: mms_coeff_b =  2368.0_rk/15.0_rk /mms_L**3
real(rk), parameter :: mms_coeff_c = -1486.0_rk/15.0_rk /mms_L**2
real(rk), parameter :: mms_coeff_d =   302.0_rk/15.0_rk /mms_L
real(rk), parameter :: mms_coeff_e = 0.0_rk
real(rk), parameter :: mms_coeff_fint = 163.0_rk/225.0_rk * mms_L
real(rk), parameter :: mms_nusf = 0.156_rk
real(rk), parameter :: mms_kappa_nu = 1.5e-11_rk
real(rk), parameter :: mms_kappasf = mms_kappa_nu * mms_nusf ! NOTE: hard-wired...
real(rk), parameter :: mms_phi0 = mms_p0/mms_kappasf * pi*0.5d0/mms_L
real(rk), parameter :: mms_alpha = 1e9_rk

public :: transient_mms_init, transient_mms_cleanup
public :: transient_mms_sigma_a

type(DelayedNeutronData) :: kin
type(XSLibrary) :: xslib
logical :: mms_negative_flag = .false.

contains

  subroutine transient_mms_init(xs, dnd)
    type(XSLibrary), intent(in) :: xs
    type(DelayedNeutronData), intent(in) :: dnd
    kin = dnd
    xslib = xs
    mms_negative_flag = .false.
  endsubroutine transient_mms_init

  subroutine transient_mms_cleanup()
    if (allocated(kin%beta)) then
      deallocate(kin%beta)
    endif
    if (allocated(kin%lamd)) then
      deallocate(kin%lamd)
    endif
    if (allocated(kin%vel)) then
      deallocate(kin%vel)
    endif
    if (allocated(xslib%mat)) then
      deallocate(xslib%mat)
    endif
  endsubroutine transient_mms_cleanup

  real(rk) pure function transient_mms_f(x)
    real(rk), intent(in) :: x
    transient_mms_f = mms_coeff_a*x**4 + mms_coeff_b*x**3 + mms_coeff_c*x**2 &
      + mms_coeff_d*x + mms_coeff_e
  endfunction transient_mms_f

  real(rk) pure function transient_mms_d2f_dx2(x)
    real(rk), intent(in) :: x
    transient_mms_d2f_dx2 = &
      12.0_rk*mms_coeff_a*x**2 + 6.0_rk*mms_coeff_b*x + 2.0_rk*mms_coeff_c
  endfunction transient_mms_d2f_dx2

  real(rk) pure function transient_mms_phi(x, t)
    real(rk), intent(in) :: x, t
    transient_mms_phi = mms_phi0 * sin(pi * x / mms_L) &
      + mms_alpha * t**2 * transient_mms_f(x)
  endfunction transient_mms_phi

  real(rk) pure function transient_mms_dphi_dt(x, t)
    real(rk), intent(in) :: x, t
    transient_mms_dphi_dt = 2.0_rk * mms_alpha * t * transient_mms_f(x)
  endfunction transient_mms_dphi_dt

  real(rk) pure function transient_mms_d2phi_dx2(x, t)
    real(rk), intent(in) :: x, t
    transient_mms_d2phi_dx2 = -mms_phi0 * (pi/mms_L)**2 * sin(pi*x/mms_L) &
      + mms_alpha * t**2 * transient_mms_d2f_dx2(x)
  endfunction transient_mms_d2phi_dx2

  real(rk) function transient_mms_kcrit()
    transient_mms_kcrit = xslib%mat(1)%nusf(1) &
      / (xslib%mat(1)%diffusion(1)*(pi/mms_L)**2 &
      + (xslib%mat(1)%sigma_t(1) - xslib%mat(1)%scatter(1,1,1)))
  endfunction transient_mms_kcrit

  real(rk) function transient_mms_c0(x)
    real(rk), intent(in) :: x
    transient_mms_c0 = kin%beta(1)/kin%lamd(1) &
      * xslib%mat(1)%nusf(1)/transient_mms_kcrit() * mms_phi0 * sin(pi*x/mms_L)
  endfunction transient_mms_c0

  real(rk) function transient_mms_c(x, t)
    real(rk), intent(in) :: x, t
    transient_mms_c = kin%beta(1) * xslib%mat(1)%nusf(1)/transient_mms_kcrit() * ( &
      mms_phi0 * sin(pi*x/mms_L) * (1.0_rk - exp(-kin%lamd(1)*t)) / kin%lamd(1) &
      + mms_alpha * transient_mms_f(x) &
      * (t**2 * kin%lamd(1)**2 - 2.0_rk * kin%lamd(1) * t + 2.0_rk - 2.0_rk * exp(-kin%lamd(1)*t)) &
      / kin%lamd(1)**3 &
      ) + exp(-kin%lamd(1)*t) * transient_mms_c0(x)
  endfunction transient_mms_c

  real(rk) function transient_mms_sigma_a(x, t)
    use exception_handler, only : exception_warning
    real(rk), intent(in) :: x, t
    transient_mms_sigma_a = ( &
      xslib%mat(1)%nusf(1)/transient_mms_kcrit() * (1.0_rk - kin%beta(1)) * transient_mms_phi(x, t) &
      + kin%lamd(1) * transient_mms_c(x, t) &
      - transient_mms_dphi_dt(x, t) / kin%vel(1) &
      + xslib%mat(1)%diffusion(1) * transient_mms_d2phi_dx2(x, t) &
      ) / transient_mms_phi(x, t)
    if ((transient_mms_sigma_a < 0d0) .and. (.not. mms_negative_flag)) then
      call exception_warning(&
        'A negative value of sigma_a has been encountered ' &
        // ' in the MMS calculation.')
      mms_negative_flag = .true.
    endif
  endfunction transient_mms_sigma_a

endmodule transient_mms
