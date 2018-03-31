module probdata_module

  use network, only: nspec
 ! Problem setup data
  use amrex_fort_module, only : rt => amrex_real
  real(rt), save         :: h_in, h_out, damn_rad
  real(rt), save :: xn_zone(nspec)

  real(rt), save          :: p_ambient, dens_ambient, exp_energy, temp_ambient, e_ambient, g, dens_incompressible

  integer, save           :: swe_to_comp_level
  integer, save ::  nsub

contains

  function rho_from_height(h, x) result(rho)
      ! Calculate density rho at depth x assuming adiabatic EoS given the fluid height h

      implicit none

      real(rt), intent(in) :: h, x
      real(rt) :: rho

      rho = dens_incompressible

  end function rho_from_height

  function p_from_height(h, x) result(p)
      ! Calculate pressure p at depth x assuming adiabatic EoS given the fluid height h

      implicit none

      real(rt), intent(in) :: h, x
      real(rt) :: p

      p = 0.5_rt * dens_incompressible * g * (h - x)**2

  end function p_from_height

  function height_from_p(p, x) result(h)
      ! Calculate the fluid height assuming adiabatic EoS given the pressure p at depth x

      implicit none

      real(rt), intent(in) :: p, x
      real(rt) :: h

      h = x + sqrt(2.0_rt * p / (dens_incompressible * g))

  end function height_from_p


  function rho_from_p(p) result(rho)
      ! Calculate the density rho assuming adiabatic EoS given pressure p

      implicit none

      real(rt), intent(in) :: p
      real(rt) :: rho

      rho = dens_incompressible

  end function rho_from_p

end module probdata_module
