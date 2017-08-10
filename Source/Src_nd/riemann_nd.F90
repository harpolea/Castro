module actual_riemann_module

  use amrex_fort_module, only : rt => amrex_real
  use bl_constants_module
  use meth_params_module, only : NQ, NVAR, NQAUX, &
                                 URHO, UMX, UMY, UMZ, &
                                 UEDEN, UEINT, UFS, UFX, &
                                 QRHO, QU, QV, QW, &
                                 QPRES, QGAME, QREINT, QFS, QFX, &
                                 QC, QGAMC, QCSML, &
                                 NGDNV, GDRHO, GDPRES, GDGAME, &
                                 npassive, upass_map, qpass_map, &
                                 small_dens, small_pres, small_temp, &
                                 use_eos_in_riemann
  use riemann_util_module


  implicit none

  private

  public :: hllc, gr_hll

  real(rt), parameter :: smallu = 1.e-12_rt
  real(rt), parameter :: small = 1.e-8_rt

contains


  subroutine HLLC(ql, qr, qpd_lo, qpd_hi, &
                  qaux, qa_lo, qa_hi, &
                  uflx, uflx_lo, uflx_hi, &
                  qint, q_lo, q_hi, &
                  idir, ilo, ihi, jlo, jhi, kc, kflux, k3d, &
                  domlo, domhi)


    ! this is an implementation of the HLLC solver described in Toro's
    ! book.  it uses the simplest estimate of the wave speeds, since
    ! those should work for a general EOS.  We also initially do the
    ! CGF Riemann construction to get pstar and ustar, since we'll
    ! need to know the pressure and velocity on the interface for the
    ! pdV term in the internal energy update.

    use mempool_module, only : bl_allocate, bl_deallocate
    use prob_params_module, only : physbc_lo, physbc_hi, &
                                   Symmetry, SlipWall, NoSlipWall

    implicit none

    integer, intent(in) :: qpd_lo(3), qpd_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: uflx_lo(3), uflx_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: idir, ilo, ihi, jlo, jhi
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) :: ql(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt), intent(in) :: qr(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)

    ! note: qaux comes in dimensioned as the fully box, so use k3d to
    ! index in z
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(inout) :: uflx(uflx_lo(1):uflx_hi(1),uflx_lo(2):uflx_hi(2),uflx_lo(3):uflx_hi(3),NVAR)
    real(rt), intent(inout) :: qint(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NGDNV)
    integer, intent(in) :: kc, kflux, k3d

    ! Note:
    !
    !  k3d: the k corresponding to the full 3d array -- it should be
    !       used for print statements or tests against domlo, domhi,
    !       etc
    !
    !  kc: the k corresponding to the 2-wide slab of k-planes, so in
    !      this routine it takes values only of 1 or 2
    !
    !  kflux: used for indexing the uflx array -- in the initial calls
    !         to cmpflx when uflx = {fx,fy,fxy,fyx,fz,fxz,fzx,fyz,fzy},
    !         kflux = kc, but in later calls, when uflx = {flux1,flux2,flux3},
    !         kflux = k3d

    integer :: i, j

    real(rt) :: rgdnv, regdnv
    real(rt) :: rl, ul, v1l, v2l, pl, rel
    real(rt) :: rr, ur, v1r, v2r, pr, rer
    real(rt) :: wl, wr, scr
    real(rt) :: rstar, cstar, estar, pstar, ustar
    real(rt) :: ro, uo, po, reo, co, gamco, entho
    real(rt) :: sgnm, spin, spout, ushock, frac
    real(rt) :: wsmall, csmall
    real(rt) :: cavg, gamcl, gamcr

    integer :: iu, iv1, iv2, sx, sy, sz
    logical :: special_bnd_lo, special_bnd_hi, special_bnd_lo_x, special_bnd_hi_x
    integer :: bnd_fac_x, bnd_fac_y, bnd_fac_z, bnd_fac
    real(rt) :: wwinv, roinv, co2inv

    real(rt) :: U_hllc_state(nvar), U_state(nvar), F_state(nvar)
    real(rt) :: S_l, S_r, S_c
    real(rt) :: gamma = 5.0d0 / 3.0d0, beta(3), alpha
    logical :: do_gr = .true.
    real(rt) :: sigma, l1, l2

    if (idir == 1) then
       iu = QU
       iv1 = QV
       iv2 = QW
       sx = 1
       sy = 0
       sz = 0
    else if (idir == 2) then
       iu = QV
       iv1 = QU
       iv2 = QW
       sx = 0
       sy = 1
       sz = 0
    else
       iu = QW
       iv1 = QU
       iv2 = QV
       sx = 0
       sy = 0
       sz = 1
    end if

    special_bnd_lo = (physbc_lo(idir) == Symmetry &
         .or.         physbc_lo(idir) == SlipWall &
         .or.         physbc_lo(idir) == NoSlipWall)
    special_bnd_hi = (physbc_hi(idir) == Symmetry &
         .or.         physbc_hi(idir) == SlipWall &
         .or.         physbc_hi(idir) == NoSlipWall)

    if (idir == 1) then
       special_bnd_lo_x = special_bnd_lo
       special_bnd_hi_x = special_bnd_hi
    else
       special_bnd_lo_x = .false.
       special_bnd_hi_x = .false.
    end if

    bnd_fac_z = 1
    if (idir == 3) then
       if ( k3d == domlo(3)   .and. special_bnd_lo .or. &
            k3d == domhi(3)+1 .and. special_bnd_hi ) then
          bnd_fac_z = 0
       end if
    end if

    do j = jlo, jhi

       bnd_fac_y = 1
       if (idir == 2) then
          if ( j == domlo(2)   .and. special_bnd_lo .or. &
               j == domhi(2)+1 .and. special_bnd_hi ) then
             bnd_fac_y = 0
          end if
       end if

       !dir$ ivdep
       do i = ilo, ihi

          rl = max(ql(i,j,kc,QRHO), small_dens)

          ! pick left velocities based on direction
          ul  = ql(i,j,kc,iu)
          v1l = ql(i,j,kc,iv1)
          v2l = ql(i,j,kc,iv2)

          pl  = max(ql(i,j,kc,QPRES ), small_pres)
          rel = ql(i,j,kc,QREINT)

          rr = max(qr(i,j,kc,QRHO), small_dens)

          ! pick right velocities based on direction
          ur  = qr(i,j,kc,iu)
          v1r = qr(i,j,kc,iv1)
          v2r = qr(i,j,kc,iv2)

          pr  = max(qr(i,j,kc,QPRES), small_pres)
          rer = qr(i,j,kc,QREINT)

          ! now we essentially do the CGF solver to get p and u on the
          ! interface, but we won't use these in any flux construction.
          csmall = max(qaux(i,j,k3d,QCSML), qaux(i-sx,j-sy,k3d-sz,QCSML) )
          cavg = HALF*(qaux(i,j,k3d,QC) + qaux(i-sx,j-sy,k3d-sz,QC))
          gamcl = qaux(i-sx,j-sy,k3d-sz,QGAMC)
          gamcr = qaux(i,j,k3d,QGAMC)

          wsmall = small_dens*csmall
          wl = max(wsmall, sqrt(abs(gamcl*pl*rl)))
          wr = max(wsmall, sqrt(abs(gamcr*pr*rr)))

          wwinv = ONE/(wl + wr)
          pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))*wwinv
          ustar = ((wl*ul + wr*ur) + (pl - pr))*wwinv

          pstar = max(pstar, small_pres)
          ! for symmetry preservation, if ustar is really small, then we
          ! set it to zero
          if (abs(ustar) < smallu*HALF*(abs(ul) + abs(ur))) then
             ustar = ZERO
          endif

          if (ustar > ZERO) then
             ro = rl
             uo = ul
             po = pl
             reo = rel
             gamco = gamcl

          else if (ustar < ZERO) then
             ro = rr
             uo = ur
             po = pr
             reo = rer
             gamco = gamcr
          else
             ro = HALF*(rl + rr)
             uo = HALF*(ul + ur)
             po = HALF*(pl + pr)
             reo = HALF*(rel + rer)
             gamco = HALF*(gamcl + gamcr)
          endif
          ro = max(small_dens, ro)

          roinv = ONE/ro
          co = sqrt(abs(gamco*po*roinv))
          co = max(csmall, co)
          co2inv = ONE/(co*co)

          rstar = ro + (pstar - po)*co2inv
          rstar = max(small_dens, rstar)

          entho = (reo + po)*co2inv/ro
          estar = reo + (pstar - po)*entho

          cstar = sqrt(abs(gamco*pstar/rstar))
          cstar = max(cstar, csmall)

          sgnm = sign(ONE, ustar)
          spout = co - sgnm*uo
          spin = cstar - sgnm*ustar
          ushock = HALF*(spin + spout)

          if (pstar-po > ZERO) then
             spin = ushock
             spout = ushock
          endif
          if (spout-spin == ZERO) then
             scr = small*cavg
          else
             scr = spout-spin
          endif
          frac = (ONE + (spout + spin)/scr)*HALF
          frac = max(ZERO, min(ONE, frac))

          rgdnv = frac*rstar + (ONE - frac)*ro
          regdnv = frac*estar + (ONE - frac)*reo

          qint(i,j,kc,iu) = frac*ustar + (ONE - frac)*uo
          qint(i,j,kc,GDPRES) = frac*pstar + (ONE - frac)*po
          qint(i,j,kc,GDGAME) = qint(i,j,kc,GDPRES)/regdnv + ONE


          ! now we do the HLLC construction

          ! Enforce that the fluxes through a symmetry plane or wall are hard zero.
          if ( special_bnd_lo_x .and. i== domlo(1) .or. &
               special_bnd_hi_x .and. i== domhi(1)+1 ) then
             bnd_fac_x = 0
          else
             bnd_fac_x = 1
          end if

          bnd_fac = bnd_fac_x*bnd_fac_y*bnd_fac_z

          sigma = cstar**2 / (gamco**2 * (1.0d0 - cstar**2))
          l1 = (ur - sqrt(sigma * (1.0d0 - ur**2 + sigma))) / (1.0d0 + sigma)
          l2 = (ul - sqrt(sigma * (1.0d0 - ul**2 + sigma))) / (1.0d0 + sigma)
          S_l = min(l1, l2)

          l1 = (ur + sqrt(sigma * (1.0d0 - ur**2 + sigma))) / (1.0d0 + sigma)
          l2 = (ul + sqrt(sigma * (1.0d0 - ul**2 + sigma))) / (1.0d0 + sigma)
          S_r = max(l1, l2)

          ! use the simplest estimates of the wave speeds
          !S_l = min(ul - sqrt(gamcl*pl/rl), ur - sqrt(gamcr*pr/rr))
          !S_r = max(ul + sqrt(gamcl*pl/rl), ur + sqrt(gamcr*pr/rr))

          ! estimate of the contact speed -- this is Toro Eq. 10.8
          S_c = (pr - pl + rl*ul*(S_l - ul) - rr*ur*(S_r - ur))/ &
               (rl*(S_l - ul) - rr*(S_r - ur))

          if (do_gr) then

              beta(:) = 0.0d0
              alpha = 1.0d0

              if (S_r <= ZERO) then
                 ! R region
                 call gr_cons_state(qr(i,j,kc,:), U_state, gamma)
                 call gr_compute_flux(idir, bnd_fac, qr(i,j,kc,:), U_state, pr, beta, alpha, F_state)

              else if (S_r > ZERO .and. S_c <= ZERO) then
                 ! R* region
                 call gr_cons_state(qr(i,j,kc,:), U_state, gamma)
                 call gr_compute_flux(idir, bnd_fac, qr(i,j,kc,:), U_state, pr, beta, alpha, F_state)

                 call gr_HLLC_state(idir, S_r, S_c, qr(i,j,kc,:), gamma, U_hllc_state)

                 ! correct the flux
                 F_state(:) = F_state(:) + S_r*(U_hllc_state(:) - U_state(:))

              else if (S_c > ZERO .and. S_l < ZERO) then
                 ! L* region
                 call gr_cons_state(ql(i,j,kc,:), U_state, gamma)
                 call gr_compute_flux(idir, bnd_fac, ql(i,j,kc,:), U_state, pl, beta, alpha, F_state)

                 call gr_HLLC_state(idir, S_l, S_c, ql(i,j,kc,:), gamma, U_hllc_state)

                 ! correct the flux
                 F_state(:) = F_state(:) + S_l*(U_hllc_state(:) - U_state(:))

              else
                 ! L region
                 call gr_cons_state(ql(i,j,kc,:), U_state, gamma)
                 call gr_compute_flux(idir, bnd_fac, ql(i,j,kc,:), U_state, pl, beta, alpha, F_state)

              endif

          else

              if (S_r <= ZERO) then
                 ! R region
                 call cons_state(qr(i,j,kc,:), U_state)
                 call compute_flux(idir, bnd_fac, U_state, pr, F_state)

              else if (S_r > ZERO .and. S_c <= ZERO) then
                 ! R* region
                 call cons_state(qr(i,j,kc,:), U_state)
                 call compute_flux(idir, bnd_fac, U_state, pr, F_state)

                 call HLLC_state(idir, S_r, S_c, qr(i,j,kc,:), U_hllc_state)

                 ! correct the flux
                 F_state(:) = F_state(:) + S_r*(U_hllc_state(:) - U_state(:))

              else if (S_c > ZERO .and. S_l < ZERO) then
                 ! L* region
                 call cons_state(ql(i,j,kc,:), U_state)
                 call compute_flux(idir, bnd_fac, U_state, pl, F_state)

                 call HLLC_state(idir, S_l, S_c, ql(i,j,kc,:), U_hllc_state)

                 ! correct the flux
                 F_state(:) = F_state(:) + S_l*(U_hllc_state(:) - U_state(:))

              else
                 ! L region
                 call cons_state(ql(i,j,kc,:), U_state)
                 call compute_flux(idir, bnd_fac, U_state, pl, F_state)

              endif

          endif

          uflx(i,j,kflux,:) = F_state(:)
       enddo
    enddo

  end subroutine HLLC


  subroutine gr_HLLC(ql, qr, qpd_lo, qpd_hi, &
                  qaux, qa_lo, qa_hi, &
                  uflx, uflx_lo, uflx_hi, &
                  qint, q_lo, q_hi, &
                  idir, ilo, ihi, jlo, jhi, kc, kflux, k3d, &
                  domlo, domhi)


    ! this is an implementation of the HLLC solver described in Toro's
    ! book.  it uses the simplest estimate of the wave speeds, since
    ! those should work for a general EOS.  We also initially do the
    ! CGF Riemann construction to get pstar and ustar, since we'll
    ! need to know the pressure and velocity on the interface for the
    ! pdV term in the internal energy update.

    use mempool_module, only : bl_allocate, bl_deallocate
    use prob_params_module, only : physbc_lo, physbc_hi, &
                                   Symmetry, SlipWall, NoSlipWall

    implicit none

    integer, intent(in) :: qpd_lo(3), qpd_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: uflx_lo(3), uflx_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: idir, ilo, ihi, jlo, jhi
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) :: ql(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt), intent(in) :: qr(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)

    ! note: qaux comes in dimensioned as the fully box, so use k3d to
    ! index in z
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(inout) :: uflx(uflx_lo(1):uflx_hi(1),uflx_lo(2):uflx_hi(2),uflx_lo(3):uflx_hi(3),NVAR)
    real(rt), intent(inout) :: qint(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NGDNV)
    integer, intent(in) :: kc, kflux, k3d

    ! Note:
    !
    !  k3d: the k corresponding to the full 3d array -- it should be
    !       used for print statements or tests against domlo, domhi,
    !       etc
    !
    !  kc: the k corresponding to the 2-wide slab of k-planes, so in
    !      this routine it takes values only of 1 or 2
    !
    !  kflux: used for indexing the uflx array -- in the initial calls
    !         to cmpflx when uflx = {fx,fy,fxy,fyx,fz,fxz,fzx,fyz,fzy},
    !         kflux = kc, but in later calls, when uflx = {flux1,flux2,flux3},
    !         kflux = k3d

    integer :: i, j

    real(rt) :: rgdnv, regdnv
    real(rt) :: rl, ul, v1l, v2l, pl, rel
    real(rt) :: rr, ur, v1r, v2r, pr, rer
    real(rt) :: wl, wr, scr
    real(rt) :: rstar, cstar, estar, pstar, ustar
    real(rt) :: ro, uo, po, reo, co, gamco, entho
    real(rt) :: sgnm, spin, spout, ushock, frac
    real(rt) :: wsmall, csmall
    real(rt) :: cavg, gamcl, gamcr

    integer :: iu, iv1, iv2, sx, sy, sz
    logical :: special_bnd_lo, special_bnd_hi, special_bnd_lo_x, special_bnd_hi_x
    integer :: bnd_fac_x, bnd_fac_y, bnd_fac_z, bnd_fac
    real(rt) :: wwinv, roinv, co2inv

    real(rt) :: U_hllc_state(nvar), U_state(nvar), F_state(nvar)
    real(rt) :: S_l, S_r, S_c
    real(rt) :: gamma = 5.0d0 / 3.0d0, beta(3), alpha
    real(rt) :: sigma, l1, l2

    if (idir == 1) then
       iu = QU
       iv1 = QV
       iv2 = QW
       sx = 1
       sy = 0
       sz = 0
    else if (idir == 2) then
       iu = QV
       iv1 = QU
       iv2 = QW
       sx = 0
       sy = 1
       sz = 0
    else
       iu = QW
       iv1 = QU
       iv2 = QV
       sx = 0
       sy = 0
       sz = 1
    end if

    special_bnd_lo = (physbc_lo(idir) == Symmetry &
         .or.         physbc_lo(idir) == SlipWall &
         .or.         physbc_lo(idir) == NoSlipWall)
    special_bnd_hi = (physbc_hi(idir) == Symmetry &
         .or.         physbc_hi(idir) == SlipWall &
         .or.         physbc_hi(idir) == NoSlipWall)

    if (idir == 1) then
       special_bnd_lo_x = special_bnd_lo
       special_bnd_hi_x = special_bnd_hi
    else
       special_bnd_lo_x = .false.
       special_bnd_hi_x = .false.
    end if

    bnd_fac_z = 1
    if (idir == 3) then
       if ( k3d == domlo(3)   .and. special_bnd_lo .or. &
            k3d == domhi(3)+1 .and. special_bnd_hi ) then
          bnd_fac_z = 0
       end if
    end if

    do j = jlo, jhi

       bnd_fac_y = 1
       if (idir == 2) then
          if ( j == domlo(2)   .and. special_bnd_lo .or. &
               j == domhi(2)+1 .and. special_bnd_hi ) then
             bnd_fac_y = 0
          end if
       end if

       !dir$ ivdep
       do i = ilo, ihi

          rl = max(ql(i,j,kc,QRHO), small_dens)

          ! pick left velocities based on direction
          ul  = ql(i,j,kc,iu)
          v1l = ql(i,j,kc,iv1)
          v2l = ql(i,j,kc,iv2)

          pl  = max(ql(i,j,kc,QPRES ), small_pres)
          rel = ql(i,j,kc,QREINT)

          rr = max(qr(i,j,kc,QRHO), small_dens)

          ! pick right velocities based on direction
          ur  = qr(i,j,kc,iu)
          v1r = qr(i,j,kc,iv1)
          v2r = qr(i,j,kc,iv2)

          pr  = max(qr(i,j,kc,QPRES), small_pres)
          rer = qr(i,j,kc,QREINT)

          ! now we essentially do the CGF solver to get p and u on the
          ! interface, but we won't use these in any flux construction.
          csmall = max(qaux(i,j,k3d,QCSML), qaux(i-sx,j-sy,k3d-sz,QCSML) )
          cavg = HALF*(qaux(i,j,k3d,QC) + qaux(i-sx,j-sy,k3d-sz,QC))
          gamcl = qaux(i-sx,j-sy,k3d-sz,QGAMC)
          gamcr = qaux(i,j,k3d,QGAMC)

          wsmall = small_dens*csmall
          wl = max(wsmall, sqrt(abs(gamcl*pl*rl)))
          wr = max(wsmall, sqrt(abs(gamcr*pr*rr)))

          wwinv = ONE/(wl + wr)
          pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))*wwinv
          ustar = ((wl*ul + wr*ur) + (pl - pr))*wwinv

          pstar = max(pstar, small_pres)
          ! for symmetry preservation, if ustar is really small, then we
          ! set it to zero
          if (abs(ustar) < smallu*HALF*(abs(ul) + abs(ur))) then
             ustar = ZERO
          endif

          if (ustar > ZERO) then
             ro = rl
             uo = ul
             po = pl
             reo = rel
             gamco = gamcl

          else if (ustar < ZERO) then
             ro = rr
             uo = ur
             po = pr
             reo = rer
             gamco = gamcr
          else
             ro = HALF*(rl + rr)
             uo = HALF*(ul + ur)
             po = HALF*(pl + pr)
             reo = HALF*(rel + rer)
             gamco = HALF*(gamcl + gamcr)
          endif
          ro = max(small_dens, ro)

          roinv = ONE/ro
          co = sqrt(abs(gamco*po*roinv))
          co = max(csmall, co)
          co2inv = ONE/(co*co)

          rstar = ro + (pstar - po)*co2inv
          rstar = max(small_dens, rstar)

          entho = (reo + po)*co2inv/ro
          estar = reo + (pstar - po)*entho

          cstar = sqrt(abs(gamco*pstar/rstar))
          cstar = max(cstar, csmall)

          sgnm = sign(ONE, ustar)
          spout = co - sgnm*uo
          spin = cstar - sgnm*ustar
          ushock = HALF*(spin + spout)

          if (pstar-po > ZERO) then
             spin = ushock
             spout = ushock
          endif
          if (spout-spin == ZERO) then
             scr = small*cavg
          else
             scr = spout-spin
          endif
          frac = (ONE + (spout + spin)/scr)*HALF
          frac = max(ZERO, min(ONE, frac))

          rgdnv = frac*rstar + (ONE - frac)*ro
          regdnv = frac*estar + (ONE - frac)*reo

          qint(i,j,kc,iu) = frac*ustar + (ONE - frac)*uo
          qint(i,j,kc,GDPRES) = frac*pstar + (ONE - frac)*po
          qint(i,j,kc,GDGAME) = qint(i,j,kc,GDPRES)/regdnv + ONE


          ! now we do the HLLC construction

          ! Enforce that the fluxes through a symmetry plane or wall are hard zero.
          if ( special_bnd_lo_x .and. i== domlo(1) .or. &
               special_bnd_hi_x .and. i== domhi(1)+1 ) then
             bnd_fac_x = 0
          else
             bnd_fac_x = 1
          end if

          bnd_fac = bnd_fac_x*bnd_fac_y*bnd_fac_z

          sigma = cstar**2 / (gamco**2 * (1.0d0 - cstar**2))
          l1 = (ur - sqrt(sigma * (1.0d0 - ur**2 + sigma))) / (1.0d0 + sigma)
          l2 = (ul - sqrt(sigma * (1.0d0 - ul**2 + sigma))) / (1.0d0 + sigma)
          S_l = min(l1, l2)

          l1 = (ur + sqrt(sigma * (1.0d0 - ur**2 + sigma))) / (1.0d0 + sigma)
          l2 = (ul + sqrt(sigma * (1.0d0 - ul**2 + sigma))) / (1.0d0 + sigma)
          S_r = max(l1, l2)

          ! use the simplest estimates of the wave speeds
          !S_l = min(ul - sqrt(gamcl*pl/rl), ur - sqrt(gamcr*pr/rr))
          !S_r = max(ul + sqrt(gamcl*pl/rl), ur + sqrt(gamcr*pr/rr))

          ! estimate of the contact speed -- this is Toro Eq. 10.8
          S_c = (pr - pl + rl*ul*(S_l - ul) - rr*ur*(S_r - ur))/ &
               (rl*(S_l - ul) - rr*(S_r - ur))

          beta(:) = 0.0d0
          alpha = 1.0d0

          if (S_r <= ZERO) then
             ! R region
             call gr_cons_state(qr(i,j,kc,:), U_state, gamma)
             call gr_compute_flux(idir, bnd_fac, qr(i,j,kc,:), U_state, pr, beta, alpha, F_state)

          else if (S_r > ZERO .and. S_c <= ZERO) then
             ! R* region
             call gr_cons_state(qr(i,j,kc,:), U_state, gamma)
             call gr_compute_flux(idir, bnd_fac, qr(i,j,kc,:), U_state, pr, beta, alpha, F_state)

             call gr_HLLC_state(idir, S_r, S_c, qr(i,j,kc,:), gamma, U_hllc_state)

             ! correct the flux
             F_state(:) = F_state(:) + S_r*(U_hllc_state(:) - U_state(:))

          else if (S_c > ZERO .and. S_l < ZERO) then
             ! L* region
             call gr_cons_state(ql(i,j,kc,:), U_state, gamma)
             call gr_compute_flux(idir, bnd_fac, ql(i,j,kc,:), U_state, pl, beta, alpha, F_state)

             call gr_HLLC_state(idir, S_l, S_c, ql(i,j,kc,:), gamma, U_hllc_state)

             ! correct the flux
             F_state(:) = F_state(:) + S_l*(U_hllc_state(:) - U_state(:))

          else
             ! L region
             call gr_cons_state(ql(i,j,kc,:), U_state, gamma)
             call gr_compute_flux(idir, bnd_fac, ql(i,j,kc,:), U_state, pl, beta, alpha, F_state)

          endif


          uflx(i,j,kflux,:) = F_state(:)
       enddo
    enddo

  end subroutine gr_HLLC

  subroutine gr_HLL(ql, qr, qpd_lo, qpd_hi, &
                  qaux, qa_lo, qa_hi, &
                  uflx, uflx_lo, uflx_hi, &
                  qint, q_lo, q_hi, &
                  idir, ilo, ihi, jlo, jhi, kc, kflux, k3d, &
                  domlo, domhi)


    ! this is an implementation of the HLLC solver described in Toro's
    ! book.  it uses the simplest estimate of the wave speeds, since
    ! those should work for a general EOS.  We also initially do the
    ! CGF Riemann construction to get pstar and ustar, since we'll
    ! need to know the pressure and velocity on the interface for the
    ! pdV term in the internal energy update.

    use mempool_module, only : bl_allocate, bl_deallocate
    use prob_params_module, only : physbc_lo, physbc_hi, &
                                   Symmetry, SlipWall, NoSlipWall

    implicit none

    integer, intent(in) :: qpd_lo(3), qpd_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: uflx_lo(3), uflx_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: idir, ilo, ihi, jlo, jhi
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) :: ql(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt), intent(in) :: qr(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)

    ! note: qaux comes in dimensioned as the fully box, so use k3d to
    ! index in z
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(inout) :: uflx(uflx_lo(1):uflx_hi(1),uflx_lo(2):uflx_hi(2),uflx_lo(3):uflx_hi(3),NVAR)
    real(rt), intent(inout) :: qint(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NGDNV)
    integer, intent(in) :: kc, kflux, k3d

    ! Note:
    !
    !  k3d: the k corresponding to the full 3d array -- it should be
    !       used for print statements or tests against domlo, domhi,
    !       etc
    !
    !  kc: the k corresponding to the 2-wide slab of k-planes, so in
    !      this routine it takes values only of 1 or 2
    !
    !  kflux: used for indexing the uflx array -- in the initial calls
    !         to cmpflx when uflx = {fx,fy,fxy,fyx,fz,fxz,fzx,fyz,fzy},
    !         kflux = kc, but in later calls, when uflx = {flux1,flux2,flux3},
    !         kflux = k3d

    integer :: i, j

    real(rt) :: rgdnv, regdnv
    real(rt) :: rl, ul, v1l, v2l, pl, rel, cl
    real(rt) :: rr, ur, v1r, v2r, pr, rer, cr
    real(rt) :: wl, wr, scr
    real(rt) :: rstar, cstar, estar, pstar, ustar
    real(rt) :: ro, uo, po, reo, co, gamco, entho
    real(rt) :: sgnm, spin, spout, ushock, frac
    real(rt) :: wsmall, csmall
    real(rt) :: cavg, gamcl, gamcr

    integer :: iu, iv1, iv2, sx, sy, sz
    logical :: special_bnd_lo, special_bnd_hi, special_bnd_lo_x, special_bnd_hi_x
    integer :: bnd_fac_x, bnd_fac_y, bnd_fac_z, bnd_fac
    real(rt) :: wwinv, roinv, co2inv

    real(rt) :: U_hll_state(nvar), U_state(nvar), F_state(nvar), Fr_state(nvar)
    real(rt) :: S_l, S_r, S_c
    real(rt) :: gamma = 5.0d0 / 3.0d0, beta(3), alpha
    real(rt) :: sigmal, sigmar, l1, l2

    if (idir == 1) then
       iu = QU
       iv1 = QV
       iv2 = QW
       sx = 1
       sy = 0
       sz = 0
    else if (idir == 2) then
       iu = QV
       iv1 = QU
       iv2 = QW
       sx = 0
       sy = 1
       sz = 0
    else
       iu = QW
       iv1 = QU
       iv2 = QV
       sx = 0
       sy = 0
       sz = 1
    end if

    special_bnd_lo = (physbc_lo(idir) == Symmetry &
         .or.         physbc_lo(idir) == SlipWall &
         .or.         physbc_lo(idir) == NoSlipWall)
    special_bnd_hi = (physbc_hi(idir) == Symmetry &
         .or.         physbc_hi(idir) == SlipWall &
         .or.         physbc_hi(idir) == NoSlipWall)

    if (idir == 1) then
       special_bnd_lo_x = special_bnd_lo
       special_bnd_hi_x = special_bnd_hi
    else
       special_bnd_lo_x = .false.
       special_bnd_hi_x = .false.
    end if

    bnd_fac_z = 1
    if (idir == 3) then
       if ( k3d == domlo(3)   .and. special_bnd_lo .or. &
            k3d == domhi(3)+1 .and. special_bnd_hi ) then
          bnd_fac_z = 0
       end if
    end if

    do j = jlo, jhi

       bnd_fac_y = 1
       if (idir == 2) then
          if ( j == domlo(2)   .and. special_bnd_lo .or. &
               j == domhi(2)+1 .and. special_bnd_hi ) then
             bnd_fac_y = 0
          end if
       end if

       !dir$ ivdep
       do i = ilo, ihi

          rl = max(ql(i,j,kc,QRHO), small_dens)

          ! pick left velocities based on direction
          ul  = ql(i,j,kc,iu)
          v1l = ql(i,j,kc,iv1)
          v2l = ql(i,j,kc,iv2)

          pl  = max(ql(i,j,kc,QPRES ), small_pres)
          rel = ql(i,j,kc,QREINT)

          rr = max(qr(i,j,kc,QRHO), small_dens)

          ! pick right velocities based on direction
          ur  = qr(i,j,kc,iu)
          v1r = qr(i,j,kc,iv1)
          v2r = qr(i,j,kc,iv2)

          pr  = max(qr(i,j,kc,QPRES), small_pres)
          rer = qr(i,j,kc,QREINT)

          ! now we essentially do the CGF solver to get p and u on the
          ! interface, but we won't use these in any flux construction.
          csmall = max(qaux(i,j,k3d,QCSML), qaux(i-sx,j-sy,k3d-sz,QCSML) )
          cavg = HALF*(qaux(i,j,k3d,QC) + qaux(i-sx,j-sy,k3d-sz,QC))
          gamcl = qaux(i-sx,j-sy,k3d-sz,QGAMC)
          gamcr = qaux(i,j,k3d,QGAMC)

          cl = sqrt(gamcl * pl / (rl + gamcl * pl / (gamcl - 1.0d0)))
          cr = sqrt(gamcr * pr / (rr + gamcr * pr / (gamcr - 1.0d0)))

          wsmall = small_dens*csmall
          wl = max(wsmall, sqrt(abs(gamcl*pl*rl)))
          wr = max(wsmall, sqrt(abs(gamcr*pr*rr)))

          wwinv = ONE/(wl + wr)
          pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))*wwinv
          ustar = ((wl*ul + wr*ur) + (pl - pr))*wwinv

          pstar = max(pstar, small_pres)
          ! for symmetry preservation, if ustar is really small, then we
          ! set it to zero
          if (abs(ustar) < smallu*HALF*(abs(ul) + abs(ur))) then
             ustar = ZERO
          endif

          if (ustar > ZERO) then
             ro = rl
             uo = ul
             po = pl
             reo = rel
             gamco = gamcl

          else if (ustar < ZERO) then
             ro = rr
             uo = ur
             po = pr
             reo = rer
             gamco = gamcr
          else
             ro = HALF*(rl + rr)
             uo = HALF*(ul + ur)
             po = HALF*(pl + pr)
             reo = HALF*(rel + rer)
             gamco = HALF*(gamcl + gamcr)
          endif
          ro = max(small_dens, ro)

          roinv = ONE/ro
          co = sqrt(abs(gamco*po / (ro + gamco * po / (gamco - 1.0d0))))
          co = max(csmall, co)
          co2inv = ONE/(co*co)

          rstar = ro + (pstar - po)*co2inv
          rstar = max(small_dens, rstar)

          entho = (reo + po)*co2inv/ro
          estar = reo + (pstar - po)*entho

          cstar = sqrt(abs(gamco*pstar/(rstar + gamco * pstar / (gamco - 1.0d0))))
          cstar = max(cstar, csmall)

          sgnm = sign(ONE, ustar)
          spout = co - sgnm*uo
          spin = cstar - sgnm*ustar
          ushock = HALF*(spin + spout)

          if (pstar-po > ZERO) then
             spin = ushock
             spout = ushock
          endif
          if (spout-spin == ZERO) then
             scr = small*cavg
          else
             scr = spout-spin
          endif
          frac = (ONE + (spout + spin)/scr)*HALF
          frac = max(ZERO, min(ONE, frac))

          rgdnv = frac*rstar + (ONE - frac)*ro
          regdnv = frac*estar + (ONE - frac)*reo

          qint(i,j,kc,iu) = frac*ustar + (ONE - frac)*uo
          qint(i,j,kc,GDPRES) = frac*pstar + (ONE - frac)*po
          qint(i,j,kc,GDGAME) = qint(i,j,kc,GDPRES)/regdnv + ONE


          ! now we do the HLLC construction

          ! Enforce that the fluxes through a symmetry plane or wall are hard zero.
          if ( special_bnd_lo_x .and. i== domlo(1) .or. &
               special_bnd_hi_x .and. i== domhi(1)+1 ) then
             bnd_fac_x = 0
          else
             bnd_fac_x = 1
          end if

          bnd_fac = bnd_fac_x*bnd_fac_y*bnd_fac_z

          sigmal = cl**2 / (gamcl**2 * (1.0d0 - cl**2))
          sigmar = cr**2 / (gamcr**2 * (1.0d0 - cr**2))
          l1 = (ur - sqrt(sigmar * (1.0d0 - ur**2 + sigmar))) / (1.0d0 + sigmar)
          l2 = (ul - sqrt(sigmal * (1.0d0 - ul**2 + sigmal))) / (1.0d0 + sigmal)
          S_l = min(l1, l2)

          l1 = (ur + sqrt(sigmar * (1.0d0 - ur**2 + sigmar))) / (1.0d0 + sigmar)
          l2 = (ul + sqrt(sigmal * (1.0d0 - ul**2 + sigmal))) / (1.0d0 + sigmal)
          S_r = max(l1, l2)

          ! use the simplest estimates of the wave speeds
          !S_l = min(ul - sqrt(gamcl*pl/rl), ur - sqrt(gamcr*pr/rr))
          !S_r = max(ul + sqrt(gamcl*pl/rl), ur + sqrt(gamcr*pr/rr))

          beta(:) = 0.0d0
          alpha = 1.0d0

          if (S_r <= ZERO) then
             ! R region
             call gr_cons_state(qr(i,j,kc,:), U_state, gamma)
             call gr_compute_flux(idir, bnd_fac, qr(i,j,kc,:), U_state, pr, beta, alpha, F_state)

         else if (S_r > ZERO .and. S_l <= ZERO) then
             ! * region
             call gr_cons_state(ql(i,j,kc,:), U_state, gamma)
             call gr_compute_flux(idir, bnd_fac, ql(i,j,kc,:), U_state, pl, beta, alpha, F_state)

             call gr_cons_state(qr(i,j,kc,:), U_hll_state, gamma)
             call gr_compute_flux(idir, bnd_fac, qr(i,j,kc,:), U_hll_state, pr, beta, alpha, Fr_state)

             ! correct the flux
             F_state(:) = (S_l * F_state(:) - S_r * Fr_state(:) + S_r * S_l * (U_hll_state(:) - U_state(:)))/ (S_r - S_l)

          else
             ! L region
             call gr_cons_state(ql(i,j,kc,:), U_state, gamma)
             call gr_compute_flux(idir, bnd_fac, ql(i,j,kc,:), U_state, pl, beta, alpha, F_state)

          endif

          uflx(i,j,kflux,:) = F_state(:)
       enddo
    enddo

  end subroutine gr_HLL


end module actual_riemann_module
