module actual_riemann_module

  use amrex_fort_module, only : rt => amrex_real
  use bl_constants_module
  use meth_params_module, only : NQ, NVAR, NQAUX, &
                                 URHO, UMX, UMY, UMZ, &
                                 UFS, UFX, &
                                 QRHO, QU, QV, QW, &
                                 QFS, QFX, &
                                 npassive, upass_map, qpass_map, &
                                 small_dens, &
                                 use_eos_in_riemann
  use riemann_util_module

  implicit none

  private

  public :: swe_hll

  real(rt), parameter :: smallu = 1.e-12_rt
  real(rt), parameter :: small = 1.e-8_rt

contains

  subroutine swe_HLL(ql, qr, qpd_lo, qpd_hi, &
                  qaux, qa_lo, qa_hi, &
                  uflx, uflx_lo, uflx_hi, &
                  idir, ilo, ihi, jlo, jhi, kc, kflux, k3d, &
                  domlo, domhi)


    ! this is an implementation of the HLL solver described in Mignone
    ! & Bodo 05. It uses the simplest estimate of the wave speeds, setting
    ! them both to +-1.
    use prob_params_module, only : physbc_lo, physbc_hi, &
                                   Symmetry, SlipWall, NoSlipWall
    use metric_module

    implicit none

    integer, intent(in) :: qpd_lo(3), qpd_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: uflx_lo(3), uflx_hi(3)
    integer, intent(in) :: idir, ilo, ihi, jlo, jhi
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) :: ql(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt), intent(in) :: qr(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)

    ! note: qaux comes in dimensioned as the fully box, so use k3d to
    ! index in z
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(inout) :: uflx(uflx_lo(1):uflx_hi(1),uflx_lo(2):uflx_hi(2),uflx_lo(3):uflx_hi(3),NVAR)
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

    integer :: iu, iv1, iv2, sx, sy, sz, lo(3), hi(3)
    logical :: special_bnd_lo, special_bnd_hi, special_bnd_lo_x, special_bnd_hi_x
    integer :: bnd_fac_x, bnd_fac_y, bnd_fac_z, bnd_fac

    real(rt) :: U_hll_state(NVAR), U_state(NVAR), F_state(NVAR), Fr_state(NVAR)
    real(rt) :: S_l, S_r, S_c, Smax_l, Smax_r, Smax
    real(rt) :: gamma_upl(ilo:ihi, jlo:jhi, qpd_lo(3):qpd_hi(3), 9)
    real(rt) :: gamma_upr(ilo:ihi, jlo:jhi, qpd_lo(3):qpd_hi(3), 9)
    real(rt) :: alphal(ilo:ihi, jlo:jhi, qpd_lo(3):qpd_hi(3))
    real(rt) :: alphar(ilo:ihi, jlo:jhi, qpd_lo(3):qpd_hi(3))
    real(rt) :: betal(ilo:ihi, jlo:jhi, qpd_lo(3):qpd_hi(3), 3)
    real(rt) :: betar(ilo:ihi, jlo:jhi, qpd_lo(3):qpd_hi(3), 3)

    lo = [ilo, jlo, qpd_lo(3)]
    hi = [ihi, jhi, qpd_hi(3)]

    call calculate_alpha(lo, hi, alphal, lo, hi, ql(:,:,:,QRHO), qpd_lo, qpd_hi)
    call calculate_alpha(lo, hi, alphar, lo, hi, qr(:,:,:,QRHO), qpd_lo, qpd_hi)

    call calculate_beta(betal, lo, hi)
    call calculate_beta(betar, lo, hi)

    call calculate_gamma_up(lo, hi, gamma_upl, lo, hi, ql(:,:,:,QRHO), qpd_lo, qpd_hi)
    call calculate_gamma_up(lo, hi, gamma_upr, lo, hi, qr(:,:,:,QRHO), qpd_lo, qpd_hi)

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

    Smax_l = maxval(abs(ql(:,:,:,QU-1+idir))) + maxval(sqrt( ql(:,:,:,QRHO)))
    Smax_r = maxval(abs(qr(:,:,:,QU-1+idir))) + maxval(sqrt( qr(:,:,:,QRHO)))
    Smax = max(Smax_r, Smax_l)

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

          ! Enforce that the fluxes through a symmetry plane or wall are hard zero.
          if ( special_bnd_lo_x .and. i== domlo(1) .or. &
               special_bnd_hi_x .and. i== domhi(1)+1 ) then
             bnd_fac_x = 0
          else
             bnd_fac_x = 1
          end if

          bnd_fac = bnd_fac_x*bnd_fac_y*bnd_fac_z

          ! signal speeds
          S_l = -1.0d0!Smax
          S_r = 1.0d0!Smax

          if (S_r <= ZERO) then
             ! R region
             call grswe_cons_state(qr(i,j,kc,:), U_state, gamma_upr(i,j,kc,:))
             call grswe_compute_flux(idir, bnd_fac, qr(i,j,kc,:), U_state, F_state, alphar(i,j,kc), betar(i,j,kc,:), gamma_upr(i,j,kc,:))

         else if (S_r > ZERO .and. S_l <= ZERO) then
             ! * region
             call grswe_cons_state(ql(i,j,kc,:), U_state, gamma_upl(i,j,kc,:))
             call grswe_compute_flux(idir, bnd_fac, ql(i,j,kc,:), U_state, F_state, alphal(i,j,kc), betal(i,j,kc,:), gamma_upl(i,j,kc,:))

             call grswe_cons_state(qr(i,j,kc,:), U_hll_state, gamma_upr(i,j,kc,:))
             call grswe_compute_flux(idir, bnd_fac, qr(i,j,kc,:), U_hll_state, Fr_state, alphar(i,j,kc), betar(i,j,kc,:), gamma_upr(i,j,kc,:))

             ! correct the flux
             F_state(:) = (S_r * F_state(:) - S_l * Fr_state(:) + S_r * S_l * (U_hll_state(:) - U_state(:)))/ (S_r - S_l)

          else
             ! L region
             call grswe_cons_state(ql(i,j,kc,:), U_state, gamma_upl(i,j,kc,:))
             call grswe_compute_flux(idir, bnd_fac, ql(i,j,kc,:), U_state, F_state, alphal(i,j,kc), betal(i,j,kc,:), gamma_upl(i,j,kc,:))

          endif

          uflx(i,j,kflux,:) = F_state(:)
       enddo
    enddo

  end subroutine swe_HLL


end module actual_riemann_module
