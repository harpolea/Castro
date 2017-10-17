module riemann_module

  use bl_types
  use bl_constants_module
  use riemann_util_module
  use meth_params_module, only : NQ, NQAUX, NVAR, QRHO, QU, QV, QW, &
                                 QPRES, QGAME, QREINT, QFS, &
                                 QFX, URHO, UMX, UMY, UMZ, UEDEN, UEINT, &
                                 UFS, UFX, &
                                 NGDNV, GDRHO, GDPRES, GDGAME, &
                                 QC, QCSML, QGAMC, &
                                 small_dens, small_pres, small_temp, &
                                 cg_maxiter, cg_tol, cg_blend, &
                                 npassive, upass_map, qpass_map, &
                                 riemann_solver, ppm_temp_fix, hybrid_riemann, &
                                 allow_negative_energy

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  private

  public cmpflx

  real(rt), parameter :: smallu = 1.e-12_rt

contains

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine cmpflx(level, qm, qp, qpd_lo, qpd_hi, &
                    flx, flx_lo, flx_hi, &
                    qint, q_lo, q_hi, &
                    qaux, qa_lo, qa_hi, &
                    shk, s_lo, s_hi, &
                    idir, ilo, ihi, jlo, jhi, kc, kflux, k3d, domlo, domhi)

    use actual_riemann_module
    use mempool_module, only : bl_allocate, bl_deallocate
    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_re
    use network, only: nspec, naux
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: level
    integer, intent(in) :: qpd_lo(3), qpd_hi(3)
    integer, intent(in) :: flx_lo(3), flx_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: s_lo(3), s_hi(3)

    integer, intent(in) :: idir,ilo,ihi,jlo,jhi,kc,kflux,k3d
    integer, intent(in) :: domlo(3),domhi(3)

    ! note: qm, qp, q come in as planes (all of x,y
    ! zones but only 2 elements in the z dir) instead of being
    ! dimensioned the same as the full box.  We index these with kc
    ! flux either comes in as planes (like qm, qp, ... above), or
    ! comes in dimensioned as the full box.  We index the flux with
    ! kflux -- this will be set correctly for the different cases.

    real(rt), intent(inout) :: qm(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt), intent(inout) :: qp(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)

    real(rt), intent(inout) ::    flx(flx_lo(1):flx_hi(1),flx_lo(2):flx_hi(2),flx_lo(3):flx_hi(3),NVAR)
    real(rt), intent(inout) ::   qint(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NGDNV)

    ! qaux come in dimensioned as the full box, so we use k3d here to
    ! index it in z

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) ::  shk(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))

    ! local variables

    integer i, j

    integer :: is_shock
    real(rt) :: cl, cr
    type (eos_t) :: eos_state

    ! Solve Riemann problem
    if (riemann_solver == 0) then
       ! Colella, Glaz, & Ferguson solver
       call riemannus(qm, qp, qpd_lo, qpd_hi, &
                      qaux, qa_lo, qa_hi, &
                      flx, flx_lo, flx_hi, &
                      qint, q_lo, q_hi, &
                      idir, ilo, ihi, jlo, jhi, kc, kflux, k3d, domlo, domhi)

    elseif (riemann_solver == 1) then
       ! Colella & Glaz solver
       call riemanncg(qm, qp, qpd_lo, qpd_hi, &
                      qaux, qa_lo, qa_hi, &
                      flx, flx_lo, flx_hi, &
                      qint, q_lo, q_hi, &
                      idir, ilo, ihi, jlo, jhi, kc, kflux, k3d, domlo, domhi)

    elseif (riemann_solver == 2) then
       ! HLLC
       call HLLC(qm, qp, qpd_lo, qpd_hi, &
                 qaux, qa_lo, qa_hi, &
                 flx, flx_lo, flx_hi, &
                 qint, q_lo, q_hi, &
                 idir, ilo, ihi, jlo, jhi, kc, kflux, k3d, domlo, domhi)
    else
       call bl_error("ERROR: invalid value of riemann_solver")
    endif

  end subroutine cmpflx

end module riemann_module
