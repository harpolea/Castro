module probdata_module

  use network, only: nspec
 ! Problem setup data
  use amrex_fort_module, only : rt => amrex_real
  real(rt), save         :: h_in, h_out
  real(rt), save :: xn_zone(nspec)

  real(rt), save          :: p_ambient, dens_ambient, exp_energy, temp_ambient, e_ambient, g, eos_K

  integer, save           :: swe_to_comp_level
  integer, save ::  nsub

contains

    function rho_from_height(h, x) result(rho)
        ! Calculate density rho at depth x assuming adiabatic EoS given the fluid height h
        use actual_eos_module, only: gamma_const

        implicit none

        real(rt), intent(in) :: h, x
        real(rt) :: rho

        rho = eos_K**(gamma_const/(gamma_const-1.0_rt)) * &
            ((gamma_const-1.0_rt)/gamma_const * g * &
            (h-x))**(1.0_rt / (gamma_const-1.0_rt))

    end function rho_from_height

    function p_from_height(h, x) result(p)
        ! Calculate pressure p at depth x assuming adiabatic EoS given the fluid height h
        use actual_eos_module, only: gamma_const

        implicit none

        real(rt), intent(in) :: h, x
        real(rt) :: p

        p = ((gamma_const - 1.0_rt) / gamma_const * g * eos_K * &
          (h - x))**(gamma_const / (gamma_const - 1.0_rt))

    end function p_from_height

    function height_from_p(p, x) result(h)
        ! Calculate the fluid height assuming adiabatic EoS given the pressure p at depth x
        use actual_eos_module, only: gamma_const

        implicit none

        real(rt), intent(in) :: p, x
        real(rt) :: h

        h = x + gamma_const / (gamma_const-1.0_rt) * &
            1.0_rt / (g * eos_K) * p**((gamma_const-1.0_rt)/gamma_const)

    end function height_from_p

    function p_from_rho(rho) result(p)
        ! Calculate the pressure p assuming adiabatic EoS given density rho
        use actual_eos_module, only: gamma_const

        implicit none

        real(rt) :: p
        real(rt), intent(in) :: rho

        p = (rho / eos_K)**gamma_const

    end function p_from_rho

    function rho_from_p(p) result(rho)
        ! Calculate the density rho assuming adiabatic EoS given pressure p
        use actual_eos_module, only: gamma_const

        implicit none

        real(rt), intent(in) :: p
        real(rt) :: rho

        rho = eos_K * p**(1.0_rt / gamma_const)

    end function rho_from_p

end module probdata_module
