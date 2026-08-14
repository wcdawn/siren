module albedo
use kind, only : rk
implicit none

private

public :: albedo_calculate_alpha, albedo_calculate_dstar

contains

  real(rk) pure function albedo_calculate_alpha(a)
    real(rk), intent(in) :: a
    albedo_calculate_alpha = 0.5_rk * (1.0_rk - a) / (1.0_rk + a)
    ! use some filtering to protect from numerical issues for vacuum bc
    albedo_calculate_alpha = min(albedo_calculate_alpha, 1e3_rk)
  endfunction albedo_calculate_alpha

  real(rk) pure function albedo_calculate_dstar(diffusion, alpha)
    real(rk), intent(in) :: diffusion ! diffusion coefficient [cm]
    real(rk), intent(in) :: alpha ! alpha albedo coefficient
    ! NOTE the input to this function is an alpha! Not a "big-A".
    albedo_calculate_dstar = diffusion / alpha
  endfunction albedo_calculate_dstar

endmodule albedo
