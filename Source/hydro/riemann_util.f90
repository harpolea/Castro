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


  subroutine gr_cons_state(q, U, gamma_up)
    ! calculates the conserved state from the primitive variables
    use meth_params_module, only: QVAR, QRHO, QU, QV, QW, QREINT, &
         NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, &
         npassive, upass_map, qpass_map, QTEMP, QPRES
    use eos_module, only: eos
    use eos_type_module, only: eos_input_re, eos_t, eos_input_rt
    use metric_module, only: calculate_scalar_W

    real(rt)        , intent(in)  :: q(QVAR), gamma_up(9)
    real(rt)        , intent(out) :: U(NVAR)

    integer  :: ipassive, n, nq
    real(rt) :: W, rhoh, p
    type (eos_t)     :: eos_state

    call calculate_scalar_W(q(QU:QW), gamma_up, W)

    U(:) = 0.0d0

    U(URHO) = q(QRHO) * W

    eos_state % rho  = q(QRHO)
    eos_state % e    = q(QREINT) / q(QRHO)
    eos_state % p = q(QPRES)
    eos_state % T = q(QTEMP)

    call eos(eos_input_re, eos_state)

    rhoh = eos_state % gam1 * q(QREINT) + q(QRHO)
    p = q(QPRES)!eos_state % p !q(QPRES)

    !write(*,*) "rho, erho, rhoh, erhoh, p, ep", q(QRHO), eos_state % rho, eos_state % h, rhoh, p, eos_state % p
    !write(*,*) "T, eT", q(QTEMP), eos_state % T

    ! since we advect all 3 velocity components regardless of dimension, this
    ! will be general
    U(UMX)  = rhoh * q(QU) * W**2
    U(UMY)  = rhoh * q(QV) * W**2
    U(UMZ)  = rhoh * q(QW) * W**2

    U(UEDEN) = rhoh * W**2 - p - U(URHO)
    U(UEINT) = q(QREINT)

    U(UTEMP) = eos_state % T!q(QTEMP)

    !do ipassive = 1, npassive
       !n  = upass_map(ipassive)
       !nq = qpass_map(ipassive)
       !U(n) = q(QRHO)*W*q(nq)
    !enddo

  end subroutine gr_cons_state


  subroutine gr_compute_flux(idir, bnd_fac, q, U, p, beta, alpha, F)
    ! returns the GR flux in direction idir
    use meth_params_module, only: QVAR, QRHO, QU, QV, QW, QREINT, &
        NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, &
         npassive, upass_map
    use prob_params_module, only : mom_flux_has_p

    integer, intent(in) :: idir, bnd_fac
    real(rt)        , intent(in)  :: q(QVAR)
    real(rt)        , intent(in) :: U(NVAR)
    real(rt)        , intent(in) :: p, beta(3), alpha
    real(rt)        , intent(out) :: F(NVAR)

    integer :: ipassive, n
    real(rt)         :: u_flx

    if (idir == 1) then
       u_flx = q(QU) - beta(1) / alpha
    elseif (idir == 2) then
       u_flx = q(QV) - beta(2) / alpha
    else
       u_flx = q(QW) - beta(3) / alpha
    endif

    if (bnd_fac == 0) then
       u_flx = ZERO
    endif

    F(:) = 0.0d0

    F(URHO) = U(URHO)*u_flx

    F(UMX) = U(UMX)*u_flx
    F(UMY) = U(UMY)*u_flx
    F(UMZ) = U(UMZ)*u_flx

    F(UMX-1+idir) = F(UMX-1+idir) + p

    F(UEINT) = U(UEINT)*u_flx
    F(UEDEN) = (U(UEDEN) + p)*u_flx + p*beta(idir)/alpha

    F(UTEMP) = 0.0d0

    !write(*,*) "FUEINT = ", F(UEINT), "FUDEN = ", F(UEDEN), "p = ", p
    !write(*,*) "F = ", F
    !stop

    !do ipassive = 1, npassive
       !n = upass_map(ipassive)
       !F(n) = U(n)*u_flx
    !enddo

end subroutine gr_compute_flux

subroutine f_of_p(f, p, U, gamma_up)
    ! function used to recover the primitive variables. Calculates the pressure using the conserved variables and a guess of the pressure, then subtracts the pressure guess; the root of this function therefore gives the correct pressure.
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UEDEN
    use eos_module, only: eos
    use eos_type_module, only: eos_input_re, eos_t
    use metric_module, only: calculate_norm
    implicit none

    double precision, intent(in)  :: U(NVAR), p, gamma_up(9)
    double precision, intent(out) :: f

    double precision :: ss, tpd
    type (eos_t)     :: eos_state

    tpd = U(UEDEN) + p + U(URHO)
    call calculate_norm(U(UMX:UMZ), gamma_up, ss)

    eos_state % rho  = U(URHO) * sqrt(tpd**2 - ss) / tpd
    eos_state % e    = (sqrt(tpd**2 - ss) - &
        p * tpd / sqrt(tpd**2 - ss) - U(URHO)) / U(URHO)

    !write(*,*) "D, p, tau, ss", U(URHO), p, U(UEDEN), ss, U(UMX), U(UMY), U(UMZ)

    call eos(eos_input_re, eos_state)

    !write(*,*) "p = ", eos_state % p, (eos_state % gam1 - 1.0d0) * eos_state % rho * eos_state % e

    f = eos_state % p - p

end subroutine f_of_p

subroutine zbrent(p, x1, b, U, gamma_up)
    ! route finder using brent's method
    use meth_params_module, only: NVAR
    implicit none

    double precision, intent(out) :: p
    double precision, intent(in)  :: U(NVAR), gamma_up(9), x1
    double precision, intent(inout) :: b

    double precision, parameter :: TOL = 1.0d-12
    integer, parameter :: ITMAX = 100

    double precision a, c, d, fa, fb, fc, fs, s
    logical mflag, con1, con2, con3, con4, con5
    integer i

    a = x1
    c = 0.0d0
    d = 0.0d0
    call f_of_p(fa, a, U, gamma_up)
    call f_of_p(fb, b, U, gamma_up)
    fc = 0.0d0

    if (fa * fb >= 0.0d0) then
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

        if (0.25d0 * (3.0d0 * a + b) < b) then
            if ( s < 0.25d0 * (3.0d0 * a + b) .or. s > b) then
                con1 = .true.
            end if
        else if (s < b .or. s > 0.25d0  * (3.0d0 * a + b)) then
            con1 = .true.
        end if

        con2 = mflag .and. abs(s - b) >= 0.5d0 * abs(b-c)

        con3 = (.not. mflag) .and. abs(s-b) >= 0.5d0 * abs(c-d)

        con4 = mflag .and. abs(b-c) < TOL

        con5 = (.not. mflag) .and. abs(c-d) < TOL

        if (con1 .or. con2 .or. con3 .or. con4 .or. con5) then
            s = 0.5d0 * (a + b)
            mflag = .true.
        else
            mflag = .false.
        end if

        call f_of_p(fs, s, U, gamma_up)

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

        if (fa * fs < 0.0d0) then
            b = s
            fb = fs
        else
            a = s
            fa = fs
        end if

        if (fb == 0.0d0 .or. fs == 0.0d0 .or. abs(b-a) < TOL) then
            p = b
            return
        end if

    end do

    p = x1

end subroutine zbrent

end module riemann_util_module
