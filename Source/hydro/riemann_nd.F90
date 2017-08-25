module actual_riemann_module

  use amrex_fort_module, only : rt => amrex_real
  use bl_constants_module
  use meth_params_module, only : NQ, NVAR, NQAUX, &
                                 URHO, UMX, UMY, UMZ, &
                                 UEDEN, UFS, UFX, &
                                 QRHO, QU, QV, QW, &
                                 QPRES, QGAME, QREINT, QFS, QFX, &
                                 QC, QGAMC, QCSML, &
                                 npassive, upass_map, qpass_map, &
                                 small_dens, small_pres, small_temp, &
                                 use_eos_in_riemann
  use riemann_util_module

  implicit none

  private

  public :: gr_hll

  real(rt), parameter :: smallu = 1.e-12_rt
  real(rt), parameter :: small = 1.e-8_rt

contains

  subroutine gr_HLL(ql, qr, qpd_lo, qpd_hi, &
                  qaux, qa_lo, qa_hi, &
                  uflx, uflx_lo, uflx_hi, &
                  idir, lo, hi, & !ilo, ihi, jlo, jhi, &!kc, kflux, k3d, &
                  domlo, domhi, &
                  gamma_up, glo, ghi)


    ! this is an implementation of the HLL solver described in Mignone
    ! & Bodo 05. It uses the simplest estimate of the wave speeds, setting
    ! them both to +-1.
    use prob_params_module, only : physbc_lo, physbc_hi, &
                                   Symmetry, SlipWall, NoSlipWall
    use metric_module, only : calculate_alpha, calculate_beta

    implicit none

    integer, intent(in) :: qpd_lo(3), qpd_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: uflx_lo(3), uflx_hi(3)
    integer, intent(in) :: idir, lo(3), hi(3)
    integer, intent(in) :: domlo(3), domhi(3)
    integer, intent(in) :: glo(3), ghi(3)

    real(rt), intent(inout) :: ql(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt), intent(in) :: qr(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)

    ! note: qaux comes in dimensioned as the fully box, so use k3d to
    ! index in z
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(in) :: gamma_up(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), 9)

    real(rt), intent(inout) :: uflx(uflx_lo(1):uflx_hi(1),uflx_lo(2):uflx_hi(2),uflx_lo(3):uflx_hi(3),NVAR)
    !integer, intent(in) :: kc, kflux, k3d

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

    integer :: i, j,k

    real(rt) :: rgdnv, regdnv
    real(rt) :: rl, ul, v1l, v2l, pl, rel, cl
    real(rt) :: rr, ur, v1r, v2r, pr, rer, cr
    real(rt) :: wl, wr, scr
    real(rt) :: sgnm, spin, spout, ushock, frac
    real(rt) :: wsmall, csmall
    real(rt) :: cavg, gamcl, gamcr

    integer :: iu, iv1, iv2, sx, sy, sz
    logical :: special_bnd_lo, special_bnd_hi, special_bnd_lo_x, special_bnd_hi_x
    integer :: bnd_fac_x, bnd_fac_y, bnd_fac_z, bnd_fac

    real(rt) :: U_hll_state(NVAR), U_state(NVAR), F_state(NVAR), Fr_state(NVAR)
    real(rt) :: S_l, S_r, S_c
    real(rt) :: beta(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 3)
    real(rt) :: alpha(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    real(rt) :: sigmal, sigmar, l1, l2

    call calculate_alpha(alpha, lo, hi)
    call calculate_beta(beta, lo, hi)

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

    do k = lo(3), hi(3)

        bnd_fac_z = 1
        if (idir == 3) then
           if ( k == domlo(3)   .and. special_bnd_lo .or. &
                k == domhi(3)+1 .and. special_bnd_hi ) then
              bnd_fac_z = 0
           end if
        end if

        do j = lo(2), hi(2)

           bnd_fac_y = 1
           if (idir == 2) then
              if ( j == domlo(2)   .and. special_bnd_lo .or. &
                   j == domhi(2)+1 .and. special_bnd_hi ) then
                 bnd_fac_y = 0
              end if
           end if

           !dir$ ivdep
           do i = lo(1), hi(1)

              rl = max(ql(i,j,k,QRHO), small_dens)

              ! pick left velocities based on direction
              ul  = ql(i,j,k,iu)
              v1l = ql(i,j,k,iv1)
              v2l = ql(i,j,k,iv2)

              pl  = max(ql(i,j,k,QPRES ), small_pres)
              rel = ql(i,j,k,QREINT)

              rr = max(qr(i,j,k,QRHO), small_dens)

              ! pick right velocities based on direction
              ur  = qr(i,j,k,iu)
              v1r = qr(i,j,k,iv1)
              v2r = qr(i,j,k,iv2)

              pr  = max(qr(i,j,k,QPRES), small_pres)
              rer = qr(i,j,k,QREINT)

              ! find sound speeds
              csmall = max(qaux(i,j,k,QCSML), qaux(i-sx,j-sy,k-sz,QCSML) )
              cavg = HALF*(qaux(i,j,k,QC) + qaux(i-sx,j-sy,k-sz,QC))
              gamcl = qaux(i-sx,j-sy,k-sz,QGAMC)
              gamcr = qaux(i,j,k,QGAMC)

              cl = sqrt(gamcl * pl / (rl + gamcl * pl / (gamcl - 1.0d0)))
              cr = sqrt(gamcr * pr / (rr + gamcr * pr / (gamcr - 1.0d0)))

              wsmall = small_dens*csmall
              wl = max(wsmall, sqrt(abs(gamcl*pl*rl)))
              wr = max(wsmall, sqrt(abs(gamcr*pr*rr)))

              ! now we do the HLL construction

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

              !write(*,*) "Sl, Sr, pl, pr", S_l, S_r, pl, pr

              ! HACK FOR NOW
              !if (pl > 1.0d0) then
                  !pl = pr
                  !ql(i,j,k,:) = qr(i,j,k,:)
              !end if

              S_l = -1.0d0
              S_r = 1.0d0

              if (S_r <= ZERO) then
                 ! R region
                 call gr_cons_state(qr(i,j,k,:), U_state, gamma_up(i,j,k,:))
                 call gr_compute_flux(idir, bnd_fac, qr(i,j,k,:), U_state, pr, beta(i,j,k,:), alpha(i,j,k), F_state)

             else if (S_r > ZERO .and. S_l <= ZERO) then
                 ! * region
                 call gr_cons_state(ql(i,j,k,:), U_state, gamma_up(i,j,k,:))
                 call gr_compute_flux(idir, bnd_fac, ql(i,j,k,:), U_state, pl, beta(i,j,k,:), alpha(i,j,k), F_state)

                 call gr_cons_state(qr(i,j,k,:), U_hll_state, gamma_up(i,j,k,:))
                 call gr_compute_flux(idir, bnd_fac, qr(i,j,k,:), U_hll_state, pr, beta(i,j,k,:), alpha(i,j,k), Fr_state)

                 ! correct the flux
                 F_state(:) = (S_r * F_state(:) - S_l * Fr_state(:) + S_r * S_l * (U_hll_state(:) - U_state(:)))/ (S_r - S_l)

              else
                 ! L region
                 call gr_cons_state(ql(i,j,k,:), U_state, gamma_up(i,j,k,:))
                 call gr_compute_flux(idir, bnd_fac, ql(i,j,k,:), U_state, pl, beta(i,j,k,:), alpha(i,j,k), F_state)

              endif

              uflx(i,j,k,:) = F_state(:)
           enddo
        enddo
    enddo

    !stop

  end subroutine gr_HLL


end module actual_riemann_module
