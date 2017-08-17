module probdata_module

 ! Problem setup data
  use amrex_fort_module, only : rt => amrex_real
  real(rt)         :: rho1, rho2, pressure

  ! Problem number
  integer :: problem

  ! Uniform flow speed
  real(rt)         :: bulk_velocity

  real(rt)          :: g

end module probdata_module
