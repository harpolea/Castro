module probdata_module

  use network, only: nspec
  use amrex_fort_module, only : rt => amrex_real

  real(rt), save :: p_ambient, dens_ambient, exp_energy, temp_ambient, e_ambient, g
  real(rt), save :: xn_zone(nspec)
  real(rt), save :: r_init
  integer, save ::  nsub
  integer, save :: swe_to_comp_level

end module probdata_module
