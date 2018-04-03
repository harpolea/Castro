module probdata_module

      use network, only: nspec

      use amrex_fort_module, only : rt => amrex_real
      real(rt), save :: xn_zone(nspec)
      real(rt)        , save ::  t_r
      real(rt)        , save :: x_r, q_r, g, eos_K
      integer, save :: nsub
      integer, save           :: swe_to_comp_level

  contains

      function rho_from_height(h, x) result(rho)
          use actual_eos_module, only: gamma_const

          implicit none

          real(rt), intent(in) :: h, x
          real(rt) :: rho

          rho = eos_K**(gamma_const/(gamma_const-1.0_rt)) * &
              ((gamma_const-1.0_rt)/gamma_const * g * &
              (h-x))**(1.0_rt / (gamma_const-1.0_rt))

      end function rho_from_height

      function p_from_height(h, x) result(p)
          use actual_eos_module, only: gamma_const

          implicit none

          real(rt), intent(in) :: h, x
          real(rt) :: p

          p = ((gamma_const - 1.0_rt) / gamma_const * g * eos_K * &
            (h - x))**(gamma_const / (gamma_const - 1.0_rt))

      end function p_from_height

      function height_from_p(p, x) result(h)
          use actual_eos_module, only: gamma_const

          implicit none

          real(rt), intent(in) :: p, x
          real(rt) :: h

          h = x + gamma_const / (gamma_const-1.0_rt) * &
              1.0_rt / (g * eos_K) * p**((gamma_const-1.0_rt)/gamma_const)

      end function height_from_p

      function p_from_rho(rho) result(p)
          use actual_eos_module, only: gamma_const

          implicit none

          real(rt) :: p
          real(rt), intent(in) :: rho

          p = (rho / eos_K)**gamma_const

      end function p_from_rho

      function rho_from_p(p) result(rho)
          use actual_eos_module, only: gamma_const

          implicit none

          real(rt), intent(in) :: p
          real(rt) :: rho

          rho = eos_K * p**(1.0_rt / gamma_const)

      end function rho_from_p

end module probdata_module
