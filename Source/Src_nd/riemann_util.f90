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


  pure subroutine wsqge(p,v,gam,gdot,gstar,pstar,wsq,csq,gmin,gmax)

    ! compute the lagrangian wave speeds -- this is the approximate
    ! version for the Colella & Glaz algorithm

    real(rt)        , intent(in) :: p,v,gam,gdot,pstar,csq,gmin,gmax
    real(rt)        , intent(out) :: wsq, gstar

    real(rt)        , parameter :: smlp1 = 1.e-10_rt
    real(rt)        , parameter :: small = 1.e-7_rt

    real(rt)         :: alpha, beta

    ! First predict a value of game across the shock

    ! CG Eq. 31
    gstar = (pstar-p)*gdot/(pstar+p) + gam
    gstar = max(gmin, min(gmax, gstar))

    ! Now use that predicted value of game with the R-H jump conditions
    ! to compute the wave speed.

    ! this is CG Eq. 34
    alpha = pstar - (gstar-ONE)*p/(gam-ONE)
    if (alpha == ZERO) alpha = smlp1*(pstar + p)

    beta = pstar + HALF*(gstar-ONE)*(pstar+p)

    wsq = (pstar-p)*beta/(v*alpha)

    if (abs(pstar - p) < smlp1*(pstar + p)) then
       wsq = csq
    endif
    wsq = max(wsq, (HALF*(gam-ONE)/gam)*csq)

    return
  end subroutine wsqge


  pure subroutine pstar_bisection(pstar_lo, pstar_hi, &
                                  ul, pl, taul, gamel, clsql, &
                                  ur, pr, taur, gamer, clsqr, &
                                  gdot, gmin, gmax, &
                                  pstar, gamstar, converged, pstar_hist_extra)

    ! we want to zero
    ! f(p*) = u*_l(p*) - u*_r(p*)
    ! we'll do bisection
    !
    ! this version is for the approximate Colella & Glaz
    ! version

    use meth_params_module, only : cg_maxiter, cg_tol

    real(rt)        , intent(inout) :: pstar_lo, pstar_hi
    real(rt)        , intent(in) :: ul, pl, taul, gamel, clsql
    real(rt)        , intent(in) :: ur, pr, taur, gamer, clsqr
    real(rt)        , intent(in) :: gdot, gmin, gmax
    real(rt)        , intent(out) :: pstar, gamstar
    logical, intent(out) :: converged
    real(rt)        , intent(out) :: pstar_hist_extra(:)

    real(rt)         :: pstar_c, ustar_l, ustar_r, f_lo, f_hi, f_c
    real(rt)         :: wl, wr, wlsq, wrsq

    integer :: iter


    ! lo bounds
    call wsqge(pl, taul, gamel, gdot,  &
               gamstar, pstar_lo, wlsq, clsql, gmin, gmax)

    call wsqge(pr, taur, gamer, gdot,  &
               gamstar, pstar_lo, wrsq, clsqr, gmin, gmax)

    wl = ONE / sqrt(wlsq)
    wr = ONE / sqrt(wrsq)

    ustar_l = ul - (pstar_lo - pstar)*wl
    ustar_r = ur + (pstar_lo - pstar)*wr

    f_lo = ustar_l - ustar_r

    ! hi bounds
    call wsqge(pl, taul, gamel, gdot,  &
               gamstar, pstar_hi, wlsq, clsql, gmin, gmax)

    call wsqge(pr, taur, gamer, gdot,  &
               gamstar, pstar_hi, wrsq, clsqr, gmin, gmax)

    wl = ONE / sqrt(wlsq)
    wr = ONE / sqrt(wrsq)

    ustar_l = ul - (pstar_hi - pstar)*wl
    ustar_r = ur + (pstar_hi - pstar)*wr

    f_hi = ustar_l - ustar_r

    ! bisection
    iter = 1
    do while (iter <= cg_maxiter)

       pstar_c = HALF * (pstar_lo + pstar_hi)

       pstar_hist_extra(iter) = pstar_c

       call wsqge(pl, taul, gamel, gdot,  &
                  gamstar, pstar_c, wlsq, clsql, gmin, gmax)

       call wsqge(pr, taur, gamer, gdot,  &
                  gamstar, pstar_c, wrsq, clsqr, gmin, gmax)

       wl = ONE / sqrt(wlsq)
       wr = ONE / sqrt(wrsq)

       ustar_l = ul - (pstar_c - pl)*wl
       ustar_r = ur - (pstar_c - pr)*wr

       f_c = ustar_l - ustar_r

       if ( HALF * abs(pstar_lo - pstar_hi) < cg_tol * pstar_c ) then
          converged = .true.
          exit
       endif

       if (f_lo * f_c < ZERO) then
          ! root is in the left half
          pstar_hi = pstar_c
          f_hi = f_c
       else
          pstar_lo = pstar_c
          f_lo = f_c
       endif
    enddo

    pstar = pstar_c

  end subroutine pstar_bisection


  subroutine HLL(ql, qr, cl, cr, idir, f)

    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, QPRES, QREINT, &
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, &
                                   npassive, upass_map, qpass_map
    use prob_params_module, only : mom_flux_has_p

    use amrex_fort_module, only : rt => amrex_real
    real(rt)        , intent(in) :: ql(QVAR), qr(QVAR), cl, cr
    real(rt)        , intent(inout) :: f(NVAR)
    integer, intent(in) :: idir

    integer :: ivel, ivelt, iveltt, imom, imomt, imomtt
    real(rt)         :: a1, a4, bd, bl, bm, bp, br
    real(rt)         :: cavg, uavg
    real(rt)         :: fl_tmp, fr_tmp
    real(rt)         :: rhod, rhoEl, rhoEr, rhol_sqrt, rhor_sqrt
    integer :: n, nq

    integer :: ipassive

    real(rt)        , parameter :: small = 1.e-10_rt

    select case (idir)
    case (1)
       ivel = QU
       ivelt = QV
       iveltt = QW

       imom = UMX
       imomt = UMY
       imomtt = UMZ

    case (2)
       ivel = QV
       ivelt = QU
       iveltt = QW

       imom = UMY
       imomt = UMX
       imomtt = UMZ

    case (3)
       ivel = QW
       ivelt = QU
       iveltt = QV

       imom = UMZ
       imomt = UMX
       imomtt = UMY

    end select

    rhol_sqrt = sqrt(ql(QRHO))
    rhor_sqrt = sqrt(qr(QRHO))

    rhod = ONE/(rhol_sqrt + rhor_sqrt)

    ! compute the average sound speed. This uses an approximation from
    ! E88, eq. 5.6, 5.7 that assumes gamma falls between 1
    ! and 5/3
    cavg = sqrt( (rhol_sqrt*cl**2 + rhor_sqrt*cr**2)*rhod + &
         HALF*rhol_sqrt*rhor_sqrt*rhod**2*(qr(ivel) - ql(ivel))**2 )

    ! Roe eigenvalues (E91, eq. 5.3b)
    uavg = (rhol_sqrt*ql(ivel) + rhor_sqrt*qr(ivel))*rhod

    a1 = uavg - cavg
    a4 = uavg + cavg

    ! signal speeds (E91, eq. 4.5)
    bl = min(a1, ql(ivel) - cl)
    br = max(a4, qr(ivel) + cr)

    bm = min(ZERO, bl)
    bp = max(ZERO, br)

    bd = bp - bm

    if (abs(bd) < small*max(abs(bm),abs(bp))) return

    bd = ONE/bd

    ! compute the fluxes according to E91, eq. 4.4b -- note that the
    ! min/max above picks the correct flux if we are not in the star
    ! region

    ! density flux
    fl_tmp = ql(QRHO)*ql(ivel)
    fr_tmp = qr(QRHO)*qr(ivel)

    f(URHO) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QRHO) - ql(QRHO))


    ! normal momentum flux.  Note for 1-d and 2-d non cartesian
    ! r-coordinate, we leave off the pressure term and handle that
    ! separately in the update, to accommodate different geometries
    fl_tmp = ql(QRHO)*ql(ivel)**2
    fr_tmp = qr(QRHO)*qr(ivel)**2
    if (mom_flux_has_p(idir)%comp(UMX-1+idir)) then
       fl_tmp = fl_tmp + ql(QPRES)
       fr_tmp = fr_tmp + qr(QPRES)
    endif

    f(imom) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QRHO)*qr(ivel) - ql(QRHO)*ql(ivel))


    ! transverse momentum flux
    fl_tmp = ql(QRHO)*ql(ivel)*ql(ivelt)
    fr_tmp = qr(QRHO)*qr(ivel)*qr(ivelt)

    f(imomt) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QRHO)*qr(ivelt) - ql(QRHO)*ql(ivelt))


    fl_tmp = ql(QRHO)*ql(ivel)*ql(iveltt)
    fr_tmp = qr(QRHO)*qr(ivel)*qr(iveltt)

    f(imomtt) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QRHO)*qr(iveltt) - ql(QRHO)*ql(iveltt))


    ! total energy flux
    rhoEl = ql(QREINT) + HALF*ql(QRHO)*(ql(ivel)**2 + ql(ivelt)**2 + ql(iveltt)**2)
    fl_tmp = ql(ivel)*(rhoEl + ql(QPRES))

    rhoEr = qr(QREINT) + HALF*qr(QRHO)*(qr(ivel)**2 + qr(ivelt)**2 + qr(iveltt)**2)
    fr_tmp = qr(ivel)*(rhoEr + qr(QPRES))

    f(UEDEN) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(rhoEr - rhoEl)


    ! eint flux
    fl_tmp = ql(QREINT)*ql(ivel)
    fr_tmp = qr(QREINT)*qr(ivel)

    f(UEINT) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QREINT) - ql(QREINT))


    ! passively-advected scalar fluxes
    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)

       fl_tmp = ql(QRHO)*ql(nq)*ql(ivel)
       fr_tmp = qr(QRHO)*qr(nq)*qr(ivel)

       f(n) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QRHO)*qr(nq) - ql(QRHO)*ql(nq))
    enddo

  end subroutine HLL


  pure subroutine cons_state(q, U)

    use meth_params_module, only: QVAR, QRHO, QU, QV, QW, QREINT, &
         NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, &
         npassive, upass_map, qpass_map

    real(rt)        , intent(in)  :: q(QVAR)
    real(rt)        , intent(out) :: U(NVAR)

    integer :: ipassive, n, nq

    U(URHO) = q(QRHO)

    ! since we advect all 3 velocity components regardless of dimension, this
    ! will be general
    U(UMX)  = q(QRHO)*q(QU)
    U(UMY)  = q(QRHO)*q(QV)
    U(UMZ)  = q(QRHO)*q(QW)

    U(UEDEN) = q(QREINT) + HALF*q(QRHO)*(q(QU)**2 + q(QV)**2 + q(QW)**2)
    U(UEINT) = q(QREINT)

    ! we don't care about T here, but initialize it to make NaN
    ! checking happy
    U(UTEMP) = ZERO

    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)
       U(n) = q(QRHO)*q(nq)
    enddo

  end subroutine cons_state


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

  pure subroutine compute_flux(idir, bnd_fac, U, p, F)

    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, &
         npassive, upass_map
    use prob_params_module, only : mom_flux_has_p

    integer, intent(in) :: idir, bnd_fac
    real(rt)        , intent(in) :: U(NVAR)
    real(rt)        , intent(in) :: p
    real(rt)        , intent(out) :: F(NVAR)

    integer :: ipassive, n
    real(rt)         :: u_flx

    if (idir == 1) then
       u_flx = U(UMX)/U(URHO)
    elseif (idir == 2) then
       u_flx = U(UMY)/U(URHO)
    elseif (idir == 3) then
       u_flx = U(UMZ)/U(URHO)
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
    F(UEDEN) = (U(UEDEN) + p)*u_flx

    F(UTEMP) = ZERO

    do ipassive = 1, npassive
       n = upass_map(ipassive)
       F(n) = U(n)*u_flx
    enddo

  end subroutine compute_flux

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
