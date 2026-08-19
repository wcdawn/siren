module delayed_neutron_data
use kind, only : rk, ik
implicit none

private

type DelayedNeutronData

  integer(ik) :: nd ! number of delayed neutron precursor families
  real(rk), allocatable :: beta(:) ! (nd) [-] Delayed neutron fractions
  real(rk), allocatable :: lamd(:) ! (nd) [1/s] Delayed neutron decay constants

  integer(ik) :: ng ! number of energy groups (for velocities)
  real(rk), allocatable :: vel(:) ! (ng) [cm/s] Neutron velocity

  real(rk) :: deltat ! [s] time-step size
  real(rk) :: tend ! [s] final time

  character(16) :: reference ! reference case (used to specify absorption xs)

  integer(ik) :: ntime = 0
  real(rk), allocatable :: tedit(:)
endtype

public :: DelayedNeutronData

endmodule delayed_neutron_data
