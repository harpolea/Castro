module riemann_module

  use amrex_fort_module, only : rt => amrex_real
  use bl_constants_module
  use meth_params_module, only : NQ, NVAR, NQAUX, &
                                 URHO, UMX, UMY, UMZ, &
                                 UFS, UFX, &
                                 QRHO, QU, QV, QW, &
                                 QFS, QFX, QC, QCSML, QREINT, QPRES, QGAMC, &
                                 npassive, upass_map, qpass_map, &
                                 small_dens, small_pres, &
                                 use_eos_in_riemann
  use riemann_util_module

  implicit none

  private

  public :: swe_hll, comp_hll, swe_to_comp, comp_to_swe, cmpflx

  real(rt), parameter :: smallu = 1.e-12_rt
  real(rt), parameter :: small = 1.e-8_rt

contains


  subroutine cmpflx(level, qm, qp, qpd_lo, qpd_hi, &
                    flx, flx_lo, flx_hi, &
                    qaux, qa_lo, qa_hi, &
                    idir, lo, hi, domlo, domhi)

    use network, only: nspec, naux
    use probdata_module, only: swe_to_comp_level
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: level
    integer, intent(in) :: qpd_lo(3), qpd_hi(3)
    integer, intent(in) :: flx_lo(3), flx_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: lo(3), hi(3)

    integer, intent(in) :: idir
    integer, intent(in) :: domlo(3),domhi(3)

    real(rt), intent(inout) :: qm(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt), intent(inout) :: qp(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt), intent(inout) ::    flx(flx_lo(1):flx_hi(1),flx_lo(2):flx_hi(2),flx_lo(3):flx_hi(3),NVAR)

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    ! local variables

    integer :: imin(3), imax(3)

    imin(:) = lo(:)
    imax(:) = hi(:)

    imax(idir) = hi(idir)+1

#if BL_SPACEDIM <= 2
    imin(3) = 0
    imax(3) = 0
#endif

    ! Solve Riemann problem (godunov state passed back, but only (u,p) saved)
    if (level <= swe_to_comp_level) then
        call swe_HLL(qm, qp, qpd_lo, qpd_hi, &
                 qaux, qa_lo, qa_hi, &
                 flx, flx_lo, flx_hi, &
                 idir, imin, imax, &
                 domlo, domhi)
    else
        call comp_HLL(qm, qp, qpd_lo, qpd_hi, &
                        qaux, qa_lo, qa_hi, &
                        flx, flx_lo, flx_hi, &
                        idir, imin, imax, &
                        domlo, domhi)
    endif

  end subroutine cmpflx

  subroutine comp_HLL(ql, qr, qpd_lo, qpd_hi, &
                  qaux, qa_lo, qa_hi, &
                  uflx, uflx_lo, uflx_hi, &
                  idir, ilo, ihi, &
                  domlo, domhi)


    ! this is an implementation of the HLL solver described in Mignone
    ! & Bodo 05. It uses the simplest estimate of the wave speeds, setting
    ! them both to +-1.
    use prob_params_module, only : physbc_lo, physbc_hi, &
                                   Symmetry, SlipWall, NoSlipWall
    use probdata_module, only : g

    implicit none

    integer, intent(in) :: qpd_lo(3), qpd_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: uflx_lo(3), uflx_hi(3)
    integer, intent(in) :: idir, ilo(3), ihi(3)
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) :: ql(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt), intent(in) :: qr(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)

    ! note: qaux comes in dimensioned as the fully box, so use k3d to
    ! index in z
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(inout) :: uflx(uflx_lo(1):uflx_hi(1),uflx_lo(2):uflx_hi(2),uflx_lo(3):uflx_hi(3),NVAR)

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

    integer :: i, j, k

    integer :: iu, iv1, iv2, sx, sy, sz
    logical :: special_bnd_lo, special_bnd_hi, special_bnd_lo_x, special_bnd_hi_x
    integer :: bnd_fac_x, bnd_fac_y, bnd_fac_z, bnd_fac

    real(rt) :: U_hll_state(NVAR), U_state(NVAR), F_state(NVAR), Fr_state(NVAR)
    real(rt) :: S_l, S_r, S_c, Smax_l, Smax_r, Smax
    real(rt) :: rl, ul, v1l, v2l, pl, rel, cl
    real(rt) :: rr, ur, v1r, v2r, pr, rer, cr
    real(rt) :: cavg, csmall, gamcl, gamcr

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

    Smax_l = maxval(abs(ql(:,:,:,QU-1+idir))) + maxval(sqrt(g * ql(:,:,:,QRHO)))
    Smax_r = maxval(abs(qr(:,:,:,QU-1+idir))) + maxval(sqrt(g * qr(:,:,:,QRHO)))
    Smax = max(Smax_r, Smax_l)

    do k = ilo(3), ihi(3)

        bnd_fac_z = 1
        if (idir == 3) then
           if ( k == domlo(3)   .and. special_bnd_lo .or. &
                k == domhi(3)+1 .and. special_bnd_hi ) then
              bnd_fac_z = 0
           end if
        end if

        do j = ilo(2), ihi(2)

           bnd_fac_y = 1
           if (idir == 2) then
              if ( j == domlo(2)   .and. special_bnd_lo .or. &
                   j == domhi(2)+1 .and. special_bnd_hi ) then
                 bnd_fac_y = 0
              end if
           end if

           !dir$ ivdep
           do i = ilo(1), ihi(1)
               F_state(:) = 0.0d0
               U_state(:) = 0.0d0

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

              ! Enforce that the fluxes through a symmetry plane or wall are hard zero.
              if ( special_bnd_lo_x .and. i== domlo(1) .or. &
                   special_bnd_hi_x .and. i== domhi(1)+1 ) then
                 bnd_fac_x = 0
              else
                 bnd_fac_x = 1
              end if

              bnd_fac = bnd_fac_x*bnd_fac_y*bnd_fac_z

              ! signal speeds
              S_l = min(-Smax, min(ul - sqrt(gamcl*pl/rl), ur - sqrt(gamcr*pr/rr)))
              S_r = max(ul + sqrt(gamcl*pl/rl), ur + sqrt(gamcr*pr/rr))
              ! if (Smax < S_r) then
              !     write(*,*) "Smax = ", Smax, "S_r = ", S_r, "gamcl = ", gamcl, "gamcr = ", gamcr
              ! endif

               S_r = max(Smax, max(ul + sqrt(gamcl*pl/rl), ur + sqrt(gamcr*pr/rr)))

               !!!!!!!!!! HACK doesn't help
               ! signal speeds
               S_l = -1.0d0
               S_r = 1.0d0

              if (S_r <= ZERO) then
                 ! R region
                 call comp_cons_state(qr(i,j,k,:), U_state)
                 call comp_compute_flux(idir, bnd_fac, qr(i,j,k,:), U_state, F_state, pr)

             else if (S_r > ZERO .and. S_l <= ZERO) then
                 ! * region
                 call comp_cons_state(ql(i,j,k,:), U_state)
                 call comp_compute_flux(idir, bnd_fac, ql(i,j,k,:), U_state, F_state, pl)

                 call comp_cons_state(qr(i,j,k,:), U_hll_state)
                 call comp_compute_flux(idir, bnd_fac, qr(i,j,k,:), U_hll_state, Fr_state, pr)

                 ! correct the flux
                 F_state(:) = (S_r * F_state(:) - S_l * Fr_state(:) + S_r * S_l * (U_hll_state(:) - U_state(:)))/ (S_r - S_l)

              else
                 ! L region
                 call comp_cons_state(ql(i,j,k,:), U_state)
                 call comp_compute_flux(idir, bnd_fac, ql(i,j,k,:), U_state, F_state, pl)

              endif

              uflx(i,j,k,1:NVAR) = F_state(1:NVAR)
           enddo
        enddo
    enddo

  end subroutine comp_HLL

  subroutine swe_HLL(ql, qr, qpd_lo, qpd_hi, &
                  qaux, qa_lo, qa_hi, &
                  uflx, uflx_lo, uflx_hi, &
                  idir, ilo, ihi, &
                  domlo, domhi)


    ! this is an implementation of the HLL solver described in Mignone
    ! & Bodo 05. It uses the simplest estimate of the wave speeds, setting
    ! them both to +-1.
    use prob_params_module, only : physbc_lo, physbc_hi, &
                                   Symmetry, SlipWall, NoSlipWall
    use probdata_module, only : g

    implicit none

    integer, intent(in) :: qpd_lo(3), qpd_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: uflx_lo(3), uflx_hi(3)
    integer, intent(in) :: idir, ilo(3), ihi(3)
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) :: ql(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt), intent(in) :: qr(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)

    ! note: qaux comes in dimensioned as the fully box, so use k3d to
    ! index in z
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(inout) :: uflx(uflx_lo(1):uflx_hi(1),uflx_lo(2):uflx_hi(2),uflx_lo(3):uflx_hi(3),NVAR)

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

    integer :: i, j, k

    integer :: iu, iv1, iv2, sx, sy, sz
    logical :: special_bnd_lo, special_bnd_hi, special_bnd_lo_x, special_bnd_hi_x
    integer :: bnd_fac_x, bnd_fac_y, bnd_fac_z, bnd_fac

    real(rt) :: U_hll_state(NVAR), U_state(NVAR), F_state(NVAR), Fr_state(NVAR)
    real(rt) :: S_l, S_r, S_c, Smax_l, Smax_r, Smax

    ! k = ilo(3)

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

    Smax_l = maxval(abs(ql(:,:,:,QU-1+idir))) + maxval(sqrt(g * ql(:,:,:,QRHO)))
    Smax_r = maxval(abs(qr(:,:,:,QU-1+idir))) + maxval(sqrt(g * qr(:,:,:,QRHO)))
    Smax = max(Smax_r, Smax_l)

    if (Smax /= Smax .or. (Smax+1 .eq. Smax)) then
        Smax = 1.0d0
    endif

    do k = ilo(3), ihi(3)

        bnd_fac_z = 1
        if (idir == 3) then
           if ( k == domlo(3)   .and. special_bnd_lo .or. &
                k == domhi(3)+1 .and. special_bnd_hi ) then
              bnd_fac_z = 0
           end if
        end if

        do j = ilo(2), ihi(2)

           bnd_fac_y = 1
           if (idir == 2) then
              if ( j == domlo(2)   .and. special_bnd_lo .or. &
                   j == domhi(2)+1 .and. special_bnd_hi ) then
                 bnd_fac_y = 0
              end if
           end if

           !dir$ ivdep
           do i = ilo(1), ihi(1)

               F_state(:) = 0.0d0
               U_state(:) = 0.0d0

              ! Enforce that the fluxes through a symmetry plane or wall are hard zero.
              if ( special_bnd_lo_x .and. i== domlo(1) .or. &
                   special_bnd_hi_x .and. i== domhi(1)+1 ) then
                 bnd_fac_x = 0
              else
                 bnd_fac_x = 1
              end if

              bnd_fac = bnd_fac_x*bnd_fac_y*bnd_fac_z

              ! signal speeds
              S_l = -Smax
              S_r = Smax

              if (S_r <= ZERO) then
                 ! R region
                 call swe_cons_state(qr(i,j,k,:), U_state)
                 call swe_compute_flux(idir, bnd_fac, qr(i,j,k,:), U_state, F_state)

             else if (S_r > ZERO .and. S_l <= ZERO) then
                 ! * region
                 call swe_cons_state(ql(i,j,k,:), U_state)
                 call swe_compute_flux(idir, bnd_fac, ql(i,j,k,:), U_state, F_state)

                 call swe_cons_state(qr(i,j,k,:), U_hll_state)
                 call swe_compute_flux(idir, bnd_fac, qr(i,j,k,:), U_hll_state, Fr_state)

                 ! correct the flux
                 F_state(:) = (S_r * F_state(:) - S_l * Fr_state(:) + S_r * S_l * (U_hll_state(:) - U_state(:)))/ (S_r - S_l)

              else
                 ! L region
                 call swe_cons_state(ql(i,j,k,:), U_state)
                 call swe_compute_flux(idir, bnd_fac, ql(i,j,k,:), U_state, F_state)

              endif

              uflx(i,j,k,1:NVAR) = F_state(1:NVAR)
           enddo
        enddo
    enddo

end subroutine swe_HLL


subroutine swe_to_comp(swe, slo, shi, comp, clo, chi, lo, hi, dx, xlo, ignore_errors)
    use meth_params_module, only: NQ, QVAR, QRHO, QU, QV, QW, &
         NVAR, URHO, UMX, UMY, UMZ, NQAUX, QTEMP, UTEMP, UEDEN, UEINT, &
         dual_energy_eta1, QPRES, UFS, UFA
    use probdata_module, only : g
    use eos_module, only : eos
    use eos_type_module, only : eos_t, eos_input_rp
    use advection_util_module, only: swectoprim, compctoprim
    use network, only : nspec

    integer, intent(in)   :: slo(3), shi(3), clo(3), chi(3), lo(3), hi(3)
    real(rt), intent(in)  :: swe(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3), NVAR)
    real(rt), intent(out) :: comp(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3), NVAR)
    real(rt), intent(in) :: dx(3), xlo(3)
    logical, optional, intent(in) :: ignore_errors

    real(rt) :: q_swe(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3), NQ)
    real(rt) :: q_comp(NQ), U_comp(NVAR)
    real(rt) :: qaux(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3),NQAUX)
    type (eos_t)     :: eos_state

    integer i, j, k, n
    logical ignore_errs
    real(rt) :: basep, xx, rho

    if (present(ignore_errors)) then
        ignore_errs = ignore_errors
    else
        ignore_errs = .false.
    endif

    ! comp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), :) = swe(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), :)
    !
    ! return
    !
    ! write(*,*) "swe_to_comp"

    ! phi = gh
    call swectoprim(lo, hi, swe, slo, shi, q_swe, slo, shi, qaux, slo, shi, ignore_errs)

    ! write(*,*) "xlo(1) = ", xlo(1), "lo(1) = ", lo(1), "dx = ", dx(1)!, "xhi = ", xlo(1) + dx(1)*dble(hi(1) - lo(1)+HALF), "h = ", q_swe(lo(1):lo(1)+5, lo(2):lo(2)+5, lo(3), QRHO)

    ! write(*,*) "dx = ", dx(1)

    ! INCOMPRESSIBLE
    rho = 1.0d0

    do k = lo(3), hi(3)
        do j = lo(2), hi(2)
            ! basep = 0.5 * g * q_swe(lo(1),j,k,QRHO)**2 !q_swe(lo(1),j,k,QPRES)
            q_comp(1:NQ) = q_swe(lo(1),j,k,1:NQ)
            U_comp(1:NVAR) = 0.0d0
            ! NOTE: incompressible for now
            q_comp(QRHO) = rho !basep / (g * q_swe(lo(1),j,k,QRHO))

            do i = lo(1), hi(1)
                xx = xlo(1) + dx(1)*dble(i-lo(1)+HALF)
                q_comp(QPRES) = 0.5d0 * g * (q_swe(i,j,k,QRHO) - xx)**2

                eos_state % rho = q_comp(QRHO)
                eos_state % p   = q_comp(QPRES)

                call eos(eos_input_rp, eos_state)

                q_comp(QREINT) = eos_state % e * q_comp(QRHO)
                q_comp(QTEMP) = swe(i,j,k,UTEMP)

                call comp_cons_state(q_comp, U_comp)

                U_comp(UFS:UFS-1+nspec) =  U_comp(URHO) / nspec

                comp(i,j,k,1:NVAR) = U_comp(1:NVAR)
            enddo
        enddo
    enddo

end subroutine swe_to_comp


subroutine comp_to_swe(swe, slo, shi, comp, clo, chi, lo, hi, xlo, dx, ignore_errors)
    use meth_params_module, only: QVAR, QRHO, QU, QV, QW, &
         NVAR, URHO, UMX, UMY, UMZ, QTEMP, UTEMP, UFS, UEDEN, UEINT, UFA
    use probdata_module, only : g
    use advection_util_module, only: compctoprim
    use network, only : nspec

    implicit none

    integer, intent(in)   :: slo(3), shi(3), clo(3), chi(3), lo(3), hi(3)
    real(rt), intent(inout)  :: swe(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3), NVAR)
    real(rt), intent(in) :: comp(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3), NVAR)
    real(rt), intent(in) :: xlo(3), dx(3)
    logical, optional, intent(in) :: ignore_errors

    real(rt) :: q_swe(NQ), U_swe(NVAR), xx
    real(rt) :: q_comp(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3), NQ)
    real(rt) :: qaux(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3), NQAUX)

    integer i, j, k
    logical ignore_errs

    if (present(ignore_errors)) then
        ignore_errs = ignore_errors
    else
        ignore_errs = .false.
    endif

    ! write(*,*) "comp to swe"

    ! phi = gh

    ! NOTE: DON'T DO THIS IS CAUSES ALL OF COMP TO BE SET TO ZERO AS WELL
    !swe(:,:,:,:) = 0.0d0

    ! write(*,*) comp(lo(1):hi(1), lo(2), lo(3), UEINT)

    call compctoprim(lo, hi, comp, clo, chi, q_comp, clo, chi, qaux, clo, chi, xlo, dx, ignore_errs)

    do k = lo(3), hi(3)
        do j = lo(2), hi(2)
            xx = xlo(1) + dx(1)*HALF

            do i = lo(1), hi(1)
                U_swe(1:NVAR) = 0.0d0
                ! look at pressure at bottom and invert to get height
                ! using p = 0.5 * g * h**2
                q_swe(1:NQ) = q_comp(i,j,k,1:NQ)
                q_swe(QRHO) = sqrt(2.0d0 * q_comp(lo(1),j,k,QPRES) / g) + xx
                q_swe(QPRES) = 0.5d0 * g * q_swe(QRHO)**2

                q_swe(QU) = 0
#if BL_SPACEDIM == 1
                q_swe(QV) = 0
#endif
#if BL_SPACEDIM <= 2
                q_swe(QW) = 0
#endif

                call swe_cons_state(q_swe, U_swe)

                U_swe(UFS:UFS-1+nspec) = U_swe(URHO) / nspec

                swe(i,j,k,1:NVAR) = U_swe(1:NVAR)
            enddo
        enddo
    enddo

end subroutine comp_to_swe

end module riemann_module
