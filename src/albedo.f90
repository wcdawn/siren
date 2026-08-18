module albedo
use kind, only : rk
implicit none

private

public :: albedo_calculate_alpha

contains

  real(rk) pure function albedo_calculate_alpha(a)
    real(rk), intent(in) :: a
    real(rk), parameter :: albedo_max = 1e6_rk
    if (abs(1.0_rk + a) < 1e-8_rk) then
      albedo_calculate_alpha = albedo_max
    else
      albedo_calculate_alpha = 0.5_rk * (1.0_rk - a) / (1.0_rk + a)
    endif
    ! use some filtering to protect from numerical issues for vacuum bc
    albedo_calculate_alpha = min(albedo_calculate_alpha, albedo_max)
  endfunction albedo_calculate_alpha

endmodule albedo
