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
    use meth_params_module, only: QVAR, QRHO, QU, QV, QW, &
         NVAR, URHO, UMX, UMY, UMZ, &
         npassive, upass_map, qpass_map

    real(rt)        , intent(in)  :: q(QVAR)
    real(rt)        , intent(out) :: U(NVAR)

    integer  :: ipassive, n, nq

    U(:) = 0.0d0

    U(:UMZ) = q(:UMZ)

    U(URHO) = q(QRHO)

    ! since we advect all 3 velocity components regardless of dimension, this
    ! will be general
    U(UMX)  = q(QRHO) * q(QU)
    U(UMY)  = q(QRHO) * q(QV)
    U(UMZ)  = q(QRHO) * q(QW)

    !do ipassive = 1, npassive
    !   n  = upass_map(ipassive)
    !   nq = qpass_map(ipassive)
    !   U(n) = q(QRHO)*q(nq)
    !enddo

end subroutine swe_cons_state


  subroutine swe_compute_flux(idir, bnd_fac, U, F)
    ! returns the GR flux in direction idir
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, &
         npassive, upass_map
    use probdata_module, only: g

    integer, intent(in) :: idir, bnd_fac
    real(rt)        , intent(in) :: U(NVAR)
    real(rt)        , intent(out) :: F(NVAR)

    integer :: ipassive, n
    real(rt)         :: u_flx

    if (idir == 1) then
       u_flx = U(UMX)
    elseif (idir == 2) then
       u_flx = U(UMY)
    elseif (idir == 3) then
       u_flx = U(UMZ)
    endif

    !if (bnd_fac == 0) then
    !   u_flx = ZERO
    !endif
    F(:) = 0.0d0

    F(URHO) = u_flx

    F(UMX) = U(UMX)*u_flx / U(URHO)
    F(UMY) = U(UMY)*u_flx / U(URHO)
    F(UMZ) = U(UMZ)*u_flx / U(URHO)

    F(UMX-1+idir) = F(UMX-1+idir) + 0.5d0 * g * U(URHO)**2

    !do ipassive = 1, npassive
    !   n = upass_map(ipassive)
    !   F(n) = U(n) * u_flx / U(URHO)
    !enddo

end subroutine swe_compute_flux

end module riemann_util_module
