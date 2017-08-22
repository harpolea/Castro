module riemann_util_module

  use bl_types
  use bl_constants_module
  use amrex_fort_module, only : rt => amrex_real
  implicit none

contains

  pure function bc_test(idir, i, j, domlo, domhi) result (f)

    use prob_params_module, only : physbc_lo, physbc_hi, Symmetry, SlipWall, NoSlipWall

    integer, intent(in) :: idir, i, j, domlo(*), domhi(*)
    integer :: f

    ! Enforce that fluxes through a symmetry plane or wall are hard zero.
    f = 1

    if (idir == 1) then
       if (i == domlo(1) .and. &
            (physbc_lo(1) == Symmetry .or. &
             physbc_lo(1) == SlipWall .or. &
             physbc_lo(1) == NoSlipWall) ) then
          f = 0
       endif

       if (i == domhi(1)+1 .and. &
            (physbc_hi(1) == Symmetry .or. &
             physbc_hi(1) == SlipWall .or. &
             physbc_hi(1) == NoSlipWall) ) then
          f = 0
       endif
    end if

    if (idir == 2) then
       if (j == domlo(2) .and. &
            (physbc_lo(2) == Symmetry .or. &
             physbc_lo(2) == SlipWall .or. &
             physbc_lo(2) == NoSlipWall) ) then
          f = 0
       endif

       if (j == domhi(2)+1 .and. &
            (physbc_hi(2) == Symmetry .or. &
             physbc_hi(2) == SlipWall .or. &
             physbc_hi(2) == NoSlipWall) ) then
          f = 0
       end if
    endif

  end function bc_test


  subroutine grswe_cons_state(q, U, gamma_up)
    ! calculates the conserved state from the primitive variables
    use meth_params_module, only: NQ, QRHO, QU, QV, QW, &
         NVAR, URHO, UMX, UMY, UMZ
    use metric_module, only : calculate_scalar_W

    real(rt)        , intent(in)  :: q(NQ), gamma_up(9)
    real(rt)        , intent(out) :: U(NVAR)

    real(rt) :: W

    U(:) = 0.0d0

    call calculate_scalar_W(q(QU:QW), gamma_up, W)

    U(URHO) = q(QRHO) * W

    U(UMX:UMZ)  = q(QRHO) * W**2 * q(QU:QW)

  end subroutine grswe_cons_state


  subroutine grswe_compute_flux(idir, bnd_fac, q, U, F, alpha, beta, gamma_up)
    ! returns the GR flux in direction idir
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, &
         NQ, QRHO, QU, QV, QW

    integer, intent(in) :: idir, bnd_fac
    real(rt)        , intent(in) :: U(NVAR), q(NQ), alpha, beta(3), gamma_up(9)
    real(rt)        , intent(out) :: F(NVAR)

    integer :: ipassive, n
    real(rt)         :: u_flx

    if (idir==1) then
        u_flx = q(QU) * gamma_up(1) + q(QV) * gamma_up(2) + q(QW) * gamma_up(3)
    elseif (idir==2) then
        u_flx = q(QU) * gamma_up(4) + q(QV) * gamma_up(5) + q(QW) * gamma_up(6)
    else
        u_flx = q(QU) * gamma_up(7) + q(QV) * gamma_up(8) + q(QW) * gamma_up(9)
    endif

    u_flx = u_flx - beta(idir) / alpha

    if (bnd_fac == 0) then
       u_flx = ZERO
    endif

    F(:) = 0.0d0

    F(URHO) = U(URHO) * u_flx
    F(UMX) = U(UMX) * u_flx
    F(UMY) = U(UMY) * u_flx
    F(UMZ) = U(UMZ) * u_flx

    F(UMX-1+idir) = F(UMX-1+idir) + 0.5d0 * q(QRHO)**2

  end subroutine grswe_compute_flux

end module riemann_util_module
