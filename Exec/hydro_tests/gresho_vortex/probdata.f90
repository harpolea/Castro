module probdata_module

      use network, only: nspec

      use amrex_fort_module, only : rt => amrex_real
      real(rt), save :: xn_zone(nspec)
      real(rt)        , save ::  t_r
      real(rt)        , save :: x_r, q_r, g, eos_K
      integer, save :: nsub
      integer, save           :: swe_to_comp_level

end module probdata_module
