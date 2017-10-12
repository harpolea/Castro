module probdata_module

 ! Problem setup data
  use amrex_fort_module, only : rt => amrex_real
  real(rt)         :: h_in, h_out, damn_rad

  real(rt)          :: g, dens_ambient

  integer           :: swe_to_comp_level

end module probdata_module
