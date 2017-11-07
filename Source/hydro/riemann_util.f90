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


  subroutine swe_cons_state(q, U)
    ! calculates the conserved state from the primitive variables
    use meth_params_module, only: NQ, QRHO, QU, QV, QW, QTEMP, &
         NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, QREINT, &
         npassive, upass_map, qpass_map, UFS
    use network, only : nspec

    real(rt)        , intent(in)  :: q(NQ)
    real(rt)        , intent(out) :: U(NVAR)

    !integer  :: ipassive, n, nq

    U(1:NVAR) = 0.0d0

    !U(:UMZ) = q(:UMZ)

    U(URHO) = q(QRHO)

    ! since we advect all 3 velocity components regardless of dimension, this
    U(UMX)  = q(QRHO) * q(QU)
    U(UMY)  = q(QRHO) * q(QV)
    U(UMZ)  = q(QRHO) * q(QW)

    ! don't care for swe but might help
    U(UEDEN) = q(QREINT) + 0.5d0*q(QRHO)*(q(QU)**2 + q(QV)**2 + q(QW)**2)
    U(UEINT) = q(QREINT)

    ! we don't care about T here, but initialize it to make NaN
    ! checking happy
    U(UTEMP) = q(QTEMP)

    U(UFS:UFS-1+nspec) = U(URHO) / nspec

  end subroutine swe_cons_state

  subroutine comp_cons_state(q, U)
    ! calculates the conserved state from the primitive variables
    use meth_params_module, only: NQ, QRHO, QU, QV, QW, QTEMP, &
         NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, QREINT, &
         npassive, upass_map, qpass_map, small_dens, small_temp, UFS
    use network, only : nspec

    real(rt)        , intent(in)  :: q(NQ)
    real(rt)        , intent(out) :: U(NVAR)

    !integer  :: ipassive, n, nq

    U(1:NVAR) = 0.0d0

    U(URHO) = q(QRHO)

    ! since we advect all 3 velocity components regardless of dimension, this
    ! will be general
    U(UMX)  = q(QRHO) * q(QU)
    U(UMY)  = q(QRHO) * q(QV)
    U(UMZ)  = q(QRHO) * q(QW)

    U(UEDEN) = q(QREINT) + 0.5d0*q(QRHO)*(q(QU)**2 + q(QV)**2 + q(QW)**2)
    U(UEINT) = q(QREINT)

    ! we don't care about T here, but initialize it to make NaN
    ! checking happy
    U(UTEMP) = q(QTEMP)

    U(UFS:UFS-1+nspec) = U(URHO) / nspec

    !write(*,*) "UMZ, QW, q(QW), U(UMZ), UEDEN, UEINT", UMZ, QW, q(QW), U(UMZ), UEDEN, UEINT

  end subroutine comp_cons_state


  subroutine swe_compute_flux(idir, bnd_fac, U, F)
    ! returns the GR flux in direction idir
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, &
         npassive, upass_map, UEDEN, UEINT, UTEMP
    use probdata_module, only: g

    integer, intent(in) :: idir, bnd_fac
    real(rt)        , intent(in) :: U(NVAR)
    real(rt)        , intent(out) :: F(NVAR)

    integer :: ipassive, n
    real(rt)         :: u_flx

    if (idir == 1) then
       u_flx = U(UMX) / U(URHO)
    elseif (idir == 2) then
       u_flx = U(UMY) / U(URHO)
    elseif (idir == 3) then
       u_flx = U(UMZ) / U(URHO)
    endif

    if (bnd_fac == 0) then
       u_flx = ZERO
    endif

    F(1:NVAR) = 0.0d0

    F(URHO) = U(URHO) * u_flx

    F(UMX) = U(UMX) * u_flx
    F(UMY) = U(UMY) * u_flx
    F(UMZ) = U(UMZ) * u_flx

    F(UMX-1+idir) = F(UMX-1+idir) + 0.5d0 * g * U(URHO)**2

    F(UEINT) = U(UEINT) * u_flx
    F(UEDEN) = (U(UEDEN) + 0.5d0 * g * U(URHO)**2) * u_flx

    F(UTEMP) = 0.0d0

  end subroutine swe_compute_flux

  subroutine comp_compute_flux(idir, bnd_fac, U, F, p)
    ! returns the GR flux in direction idir
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UTEMP, &
         UEDEN, UEINT, npassive, upass_map
    use probdata_module, only: g

    integer, intent(in) :: idir, bnd_fac
    real(rt)        , intent(in) :: U(NVAR)
    real(rt)        , intent(inout) :: p
    real(rt)        , intent(out) :: F(NVAR)

    integer :: ipassive, n
    real(rt)         :: u_flx

    if (idir == 1) then
       u_flx = U(UMX) / U(URHO)
    elseif (idir == 2) then
       u_flx = U(UMY) / U(URHO)
    elseif (idir == 3) then
       u_flx = U(UMZ) / U(URHO)
    endif

    if (bnd_fac == 0) then
       u_flx = ZERO
    endif
    F(1:NVAR) = 0.0d0

    F(URHO) = U(URHO) * u_flx

    F(UMX) = U(UMX) * u_flx
    F(UMY) = U(UMY) * u_flx
    F(UMZ) = U(UMZ) * u_flx

    F(UMX-1+idir) = F(UMX-1+idir) + p

    F(UEINT) = U(UEINT) * u_flx
    F(UEDEN) = (U(UEDEN) + p) * u_flx

    F(UTEMP) = 0.0d0

end subroutine comp_compute_flux

end module riemann_util_module
