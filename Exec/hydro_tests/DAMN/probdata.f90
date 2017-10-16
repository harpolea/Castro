module probdata_module

  use network, only: nspec
 ! Problem setup data
  use amrex_fort_module, only : rt => amrex_real
  real(rt), save         :: h_in, h_out, damn_rad
  real(rt), save :: xn_zone(nspec)

  real(rt), save          :: p_ambient, dens_ambient, exp_energy, temp_ambient, e_ambient, g

  integer, save           :: swe_to_comp_level
  real(rt), save :: r_init
  integer, save ::  nsub

end module probdata_module
