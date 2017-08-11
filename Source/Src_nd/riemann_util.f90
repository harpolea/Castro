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


  pure subroutine gr_cons_state(q, U, gamma)

    use meth_params_module, only: QVAR, QRHO, QU, QV, QW, QREINT, &
         NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, &
         npassive, upass_map, qpass_map

    real(rt)        , intent(in)  :: q(QVAR), gamma
    real(rt)        , intent(out) :: U(NVAR)

    integer  :: ipassive, n, nq
    real(rt) :: W, rhoh, p

    W = 1.0d0 / sqrt(1.0d0 - sum(q(QU:QW)**2))

    U(URHO) = q(QRHO) * W
    rhoh = gamma * q(QREINT) / q(QRHO) + (1.0d0 - gamma) * q(QRHO) !W * (q(QRHO) + gamma * q(QREINT))
    p = (gamma - 1.0d0) * (q(QREINT) / q(QRHO) - q(QRHO))

    ! since we advect all 3 velocity components regardless of dimension, this
    ! will be general
    U(UMX)  = rhoh * q(QU) * W**2
    U(UMY)  = rhoh * q(QV) * W**2
    U(UMZ)  = rhoh * q(QW) * W**2

    U(UEDEN) = rhoh * W**2 - p - U(URHO) !q(QREINT) + HALF*q(QRHO)*(q(QU)**2 + q(QV)**2 + q(QW)**2)
    U(UEINT) = q(QREINT)

    ! we don't care about T here, but initialize it to make NaN
    ! checking happy
    U(UTEMP) = ZERO

    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)
       U(n) = q(QRHO)*W*q(nq)
    enddo

  end subroutine gr_cons_state


  pure subroutine gr_compute_flux(idir, bnd_fac, q, U, p, beta, alpha, F)

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
       u_flx = q(QU) - beta(1) / alpha !U(UMX)/U(URHO)
    elseif (idir == 2) then
       u_flx = q(QV) - beta(2) / alpha !U(UMY)/U(URHO)
    elseif (idir == 3) then
       u_flx = q(QW) - beta(3) / alpha !U(UMZ)/U(URHO)
    endif

    if (bnd_fac == 0) then
       u_flx = ZERO
    endif

    F(URHO) = U(URHO)*u_flx

    F(UMX) = U(UMX)*u_flx
    F(UMY) = U(UMY)*u_flx
    F(UMZ) = U(UMZ)*u_flx

    if (mom_flux_has_p(idir)%comp(UMX-1+idir)) then
       ! we do not include the pressure term in any non-Cartesian
       ! coordinate directions
       F(UMX-1+idir) = F(UMX-1+idir) + p
    endif

    F(UEINT) = U(UEINT)*u_flx
    F(UEDEN) = (U(UEDEN) + p)*u_flx + p*beta(idir)/alpha

    F(UTEMP) = ZERO

    do ipassive = 1, npassive
       n = upass_map(ipassive)
       F(n) = U(n)*u_flx
    enddo

end subroutine gr_compute_flux

subroutine f_of_p(f, p, U, gamma, gamma_up)
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UEDEN
    implicit none

    double precision, intent(in)  :: U(NVAR), p, gamma, gamma_up(9)
    double precision, intent(out) :: f

    double precision :: sq

    sq = sqrt((U(UEDEN) + p + U(URHO))**2 - U(UMX)**2*gamma_up(1) - &
        2.0d0 * U(UMX) * U(UMY) * gamma_up(2) - &
        2.0d0 * U(UMX) * U(UMZ) * gamma_up(3) - &
        U(UMY)**2 * gamma_up(5) - &
        2.0d0 * U(UMY) * U(UMZ) * gamma_up(6) - &
        U(UMZ)**2 * gamma_up(9))

    f = (gamma - 1.0d0) * sq / (U(UEDEN) + p + U(URHO)) * &
        (sq - p * (U(UEDEN) + p + U(URHO)) / sq - U(URHO)) - p

end subroutine f_of_p

subroutine zbrent(p, x1, b, U, gamma, gamma_up)
    ! route finder using brent's method
    use meth_params_module, only: NVAR
    implicit none

    double precision, intent(out) :: p
    double precision, intent(in)  :: U(NVAR), gamma, gamma_up(9), x1
    double precision, intent(inout) :: b

    double precision, parameter :: TOL = 1.0d-12
    integer, parameter :: ITMAX = 100

    double precision a, c, d, fa, fb, fc, fs, s
    logical mflag, con1, con2, con3, con4, con5
    integer i

    a = x1
    c = 0.0d0
    d = 0.0d0
    call f_of_p(fa, a, U, gamma, gamma_up)
    call f_of_p(fb, b, U, gamma, gamma_up)
    fc = 0.0d0

    if (fa * fb >= 0.0d0) then
        p = b
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

        call f_of_p(fs, s, U, gamma, gamma_up)

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
