module riemann_module

  use bl_types
  use bl_constants_module
  use riemann_util_module
  use meth_params_module, only : NQ, NQAUX, NVAR
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  private

  public cmpflx

  real(rt), parameter :: smallu = 1.e-12_rt

contains

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine cmpflx(qm, qp, qpd_lo, qpd_hi, &
                    flx, flx_lo, flx_hi, &
                    qaux, qa_lo, qa_hi, &
                    idir, lo, hi, domlo, domhi, &
                    gamma_up, glo, ghi)

    use actual_riemann_module
    use mempool_module, only : bl_allocate, bl_deallocate
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: qpd_lo(3), qpd_hi(3)
    integer, intent(in) :: flx_lo(3), flx_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: glo(3), ghi(3)

    integer, intent(in) :: idir, lo(3), hi(3)
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

    ! qaux come in dimensioned as the full box, so we use k3d here to
    ! index it in z

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(in) :: gamma_up(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),9)

    ! local variables

    integer i, j
    integer Ilo(3), Ihi(3)

    Ilo(:) = lo(:)
    Ihi(:) = hi(:)

    ! these will refer to the zone interfaces that we solve the
    ! Riemann problem across
    if (idir == 1) then
       !imin = ilo
       !imax = ihi+1
       !jmin = jlo
       !jmax = jhi
       Ihi(1) = Ihi(1) + 1
   elseif (idir == 2) then
       !imin = ilo
       !imax = ihi
       !jmin = jlo
       !jmax = jhi+1
       Ihi(2) = Ihi(2) + 1
   else
       Ihi(3) = Ihi(3) + 1
   endif

    ! Solve Riemann problem (godunov state passed back, but only (u,p) saved)
    call gr_HLL(qm, qp, qpd_lo, qpd_hi, &
             qaux, qa_lo, qa_hi, &
             flx, flx_lo, flx_hi, &
             idir, Ilo, Ihi, &
             domlo, domhi, &
             gamma_up, glo, ghi)

  end subroutine cmpflx

end module riemann_module
