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
         npassive, upass_map, qpass_map, UFS, UFA
    use network, only : nspec
    use metric_module, only : calculate_scalar_W

    real(rt)        , intent(in)  :: q(NQ)
    real(rt)        , intent(out) :: U(NVAR)

    real(rt) :: W

    !integer  :: ipassive, n, nq

    U(1:NVAR) = 0.0e0_rt

    call calculate_scalar_W(q(QU:QW), W)

    U(URHO) = q(QRHO) * W

    U(UMY:UMZ) = q(QRHO) * W**2 * q(QV:QW)
    ! U(UMX)  = 0.0e0_rt! q(QRHO) * q(QW)

    ! since we advect all 3 velocity components regardless of dimension, this
#if BL_SPACEDIM == 1
    U(UMY)  = 0.0e0_rt! q(QRHO) * q(QU)
#endif
#if BL_SPACEDIM <= 2
    U(UMZ)  = 0.0e0_rt!q(QRHO) * q(QV)
#endif

    ! don't care for swe but might help
    U(UEDEN) = q(QREINT) + 0.5e0_rt*q(QRHO)*(q(QU)**2 + q(QV)**2 + q(QW)**2)
    U(UEINT) = q(QREINT)

    ! we don't care about T here, but initialize it to make NaN
    ! checking happy
    U(UTEMP) = q(QTEMP)
    U(UFA) = 0.0e0_rt
    U(UFS:UFS-1+nspec) = U(URHO) / nspec

  end subroutine swe_cons_state

  subroutine comp_cons_state(q, U)
    ! calculates the conserved state from the primitive variables
    use meth_params_module, only: NQ, QRHO, QU, QV, QW, QTEMP, &
         NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, QREINT, &
         npassive, upass_map, qpass_map, UFS, UFA, QPRES, small_dens, small_temp
    use network, only : nspec
    use eos_module, only: eos, eos_init, initialized
    use eos_type_module, only: eos_input_re, eos_t, eos_input_rt, eos_input_rp
    use metric_module, only: calculate_scalar_W
    use probdata_module, only: dens_incompressible

    real(rt)        , intent(in)  :: q(NQ)
    real(rt)        , intent(out) :: U(NVAR)

    real(rt) :: W, rhoh, p, rho
    type(eos_t) :: eos_state

    if (.not. initialized) call eos_init(small_dens=small_dens, small_temp=small_temp)

    call calculate_scalar_W(q(QU:QW), W)
    !integer  :: ipassive, n, nq
    U(1:NVAR) = 0.0e0_rt

    ! INCOMPRESSIBLE
    rho = dens_incompressible

    U(URHO) = rho * W !q(QRHO)

    eos_state % rho  = rho
    eos_state % e    = q(QREINT) / rho
    eos_state % p = q(QPRES)
    eos_state % T = q(QTEMP)

    call eos(eos_input_re, eos_state)

    rhoh = eos_state % gam1 * q(QREINT) + rho
    p = eos_state % p

    ! since we advect all 3 velocity components regardless of dimension, this
    ! will be general
    U(UMX)  = rhoh * q(QU) * W**2
    U(UMY)  = rhoh * q(QV) * W**2
    U(UMZ)  = rhoh * q(QW) * W**2

#if BL_SPACEDIM == 1
    U(UMY) = 0.0e0_rt
#endif
#if BL_SPACEDIM <= 2
    U(UMZ) = 0.0e0_rt
#endif

    U(UEDEN) = rhoh * W**2 - p - U(URHO)
    U(UEINT) = q(QREINT)

    ! we don't care about T here, but initialize it to make NaN
    ! checking happy
    U(UTEMP) = q(QTEMP)
    U(UFA) = 0.0e0_rt
    U(UFS:UFS-1+nspec) = U(URHO) / nspec

  end subroutine comp_cons_state


  subroutine swe_compute_flux(idir, bnd_fac, q, U, F)
    ! returns the GR flux in direction idir
    use meth_params_module, only: NQ, NVAR, URHO, UMX, UMY, UMZ, &
         npassive, upass_map, UEDEN, UEINT, UTEMP, UFS, UFA, &
         QU, QV, QW, QRHO
    use probdata_module, only: g
    use network, only : nspec

    integer, intent(in) :: idir, bnd_fac
    real(rt)        , intent(in) :: q(NQ)
    real(rt)        , intent(in) :: U(NVAR)
    real(rt)        , intent(out) :: F(NVAR)

    integer :: ipassive, n
    real(rt)         :: u_flx

    if (idir == 1) then
       u_flx = q(QU)
    elseif (idir == 2) then
       u_flx = q(QV)
    elseif (idir == 3) then
       u_flx = q(QW)
    endif

    if (bnd_fac == 0) then
       u_flx = ZERO
    endif

    F(1:NVAR) = 0.0e0_rt

    F(URHO) = U(URHO) * u_flx

    F(UMX) = 0.0e0_rt!U(UMX) * u_flx
    F(UMY) = U(UMY) * u_flx
    F(UMZ) = U(UMZ) * u_flx

    F(UMX-1+idir) = F(UMX-1+idir) + 0.5e0_rt * g * q(QRHO)**2

#if BL_SPACEDIM == 1
    F(UMY) = 0.0e0_rt
#endif
#if BL_SPACEDIM <= 2
    F(UMZ) = 0.0e0_rt
#endif

    F(UEINT) = U(UEINT) * u_flx
    F(UEDEN) = (U(UEDEN) + 0.5e0_rt * g * q(QRHO)**2) * u_flx

    F(UTEMP) = 0.0e0_rt
    F(UFA) = U(UFA) * u_flx
    F(UFS:UFS-1+nspec) = U(UFS:UFS-1+nspec) * u_flx

  end subroutine swe_compute_flux

  subroutine comp_compute_flux(idir, bnd_fac, q, U, F, p)
    ! returns the GR flux in direction idir
    use meth_params_module, only: NQ, NVAR, URHO, UMX, UMY, UMZ, UTEMP, &
         UEDEN, UEINT, npassive, upass_map, UFS, UFA, QU, QV, QW
    use probdata_module, only: g
    use network, only : nspec

    integer, intent(in) :: idir, bnd_fac
    real(rt) , intent(in) :: q(NQ)
    real(rt)        , intent(in) :: U(NVAR)
    real(rt)        , intent(inout) :: p
    real(rt)        , intent(out) :: F(NVAR)

    integer :: ipassive, n
    real(rt)         :: u_flx

    if (idir == 1) then
       u_flx = q(QU)
    elseif (idir == 2) then
       u_flx = q(QV)
    elseif (idir == 3) then
       u_flx = q(QW)
    endif

    if (bnd_fac == 0) then
       u_flx = ZERO
    endif
    F(1:NVAR) = 0.0e0_rt

    ! NOTE: Made incompressible for now
    F(URHO) = 0.0e0_rt!U(URHO) * u_flx

    F(UMX) = U(UMX) * u_flx
    F(UMY) = U(UMY) * u_flx
    F(UMZ) = U(UMZ) * u_flx

    F(UMX-1+idir) = F(UMX-1+idir) + p

#if BL_SPACEDIM == 1
    F(UMY) = 0.0e0_rt
#endif
#if BL_SPACEDIM <= 2
    F(UMZ) = 0.0e0_rt
#endif

    F(UEINT) = U(UEINT) * u_flx
    F(UEDEN) = (U(UEDEN) + p) * u_flx

    F(UTEMP) = 0.0e0_rt
    F(UFA) = U(UFA) * u_flx
    F(UFS:UFS-1+nspec) = U(UFS:UFS-1+nspec) * u_flx

  end subroutine comp_compute_flux

  subroutine f_of_p(f, p, U)
    ! function used to recover the primitive variables. Calculates the pressure using the conserved variables and a guess of the pressure, then subtracts the pressure guess; the root of this function therefore gives the correct pressure.
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UEDEN
    use eos_module, only: eos
    use eos_type_module, only: eos_input_re, eos_t
    use probdata_module, only: dens_incompressible

    double precision, intent(in)  :: U(NVAR), p
    double precision, intent(out) :: f

    double precision :: ss, tpd, rho
    type (eos_t)     :: eos_state

    rho = dens_incompressible !INCOMPRESSIBLE

    tpd = U(UEDEN) + p + U(URHO)
    ss = sum(U(UMX:UMZ)**2) ! norm of S

    eos_state % rho  = rho! U(URHO) * sqrt(tpd**2 - ss) / tpd
    eos_state % e    = (sqrt(tpd**2 - ss) - &
        p * tpd / sqrt(tpd**2 - ss) - U(URHO)) / U(URHO)

    !write(*,*) "D, p, tau, ss", U(URHO), p, U(UEDEN), ss, U(UMX), U(UMY), U(UMZ)

    call eos(eos_input_re, eos_state)

    ! if (p < 1.0e-4_rt) then
    !     write(*,*) "p = ", eos_state % p, p, eos_state % rho, U(UEDEN), eos_state % e
    ! endif

    f = eos_state % p - p

  end subroutine f_of_p

  subroutine zbrent(p, x1, b, U)
    ! route finder using brent's method
    use meth_params_module, only: NVAR
    implicit none

    double precision, intent(out) :: p
    double precision, intent(in)  :: U(NVAR), x1
    double precision, intent(inout) :: b

    double precision, parameter :: TOL = 1.0e-6_rt
    integer, parameter :: ITMAX = 100

    double precision a, c, d, fa, fb, fc, fs, s
    logical mflag, con1, con2, con3, con4, con5
    integer i

    a = x1
    c = 0.0e0_rt
    d = 0.0e0_rt
    call f_of_p(fa, a, U)
    call f_of_p(fb, b, U)
    fc = 0.0e0_rt

    if (fa * fb >= 0.0e0_rt) then
        p = x1
        return
    end if

    if (abs(fa) < abs(fb)) then
        d = a
        a = b
        b = d

        d = fa
        fa = fb
        fb = d
    end if

    c = a
    fc = fa

    mflag = .true.

    do i = 1, ITMAX
        if (fa /= fc .and. fb /= fc) then
            s = a*fb*fc / ((fa-fb) * (fa-fc)) + b*fa*fc / ((fb-fa)*(fb-fc)) +&
                c*fa*fb / ((fc-fa)*(fc-fb))
        else
            s = b - fb * (b-a) / (fb-fa)
        end if

        con1 = .false.

        if (0.25e0_rt * (3.0e0_rt * a + b) < b) then
            if ( s < 0.25e0_rt * (3.0e0_rt * a + b) .or. s > b) then
                con1 = .true.
            end if
        else if (s < b .or. s > 0.25e0_rt  * (3.0e0_rt * a + b)) then
            con1 = .true.
        end if

        con2 = mflag .and. abs(s - b) >= 0.5e0_rt * abs(b-c)

        con3 = (.not. mflag) .and. abs(s-b) >= 0.5e0_rt * abs(c-d)

        con4 = mflag .and. abs(b-c) < TOL

        con5 = (.not. mflag) .and. abs(c-d) < TOL

        if (con1 .or. con2 .or. con3 .or. con4 .or. con5) then
            s = 0.5e0_rt * (a + b)
            mflag = .true.
        else
            mflag = .false.
        end if

        call f_of_p(fs, s, U)

        if (abs(fa) < abs(fb)) then
            d = a
            a = b
            b = d

            d = fa
            fa = fb
            fb = d
        end if

        d = c
        c = b
        fc = fb

        if (fa * fs < 0.0e0_rt) then
            b = s
            fb = fs
        else
            a = s
            fa = fs
        end if

        if (fb == 0.0e0_rt .or. fs == 0.0e0_rt .or. abs(b-a) < TOL) then
            p = b
            return
        end if

    end do

    p = x1

  end subroutine zbrent

end module riemann_util_module
