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
                    idir, ilo, ihi, jlo, jhi, domlo, domhi, &
                    gamma_up, glo, ghi)

    use actual_riemann_module
    use eos_type_module, only: eos_input_re, eos_input_rt, eos_t
    use eos_module, only: eos
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: qpd_lo(3), qpd_hi(3)
    integer, intent(in) :: flx_lo(3), flx_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: glo(3), ghi(3)

    integer, intent(in) :: idir,ilo,ihi,jlo,jhi
    integer, intent(in) :: domlo(2),domhi(2)

    real(rt), intent(inout) ::  qm(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),NQ)
    real(rt), intent(inout) ::  qp(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),NQ)
    real(rt), intent(inout) :: flx(flx_lo(1):flx_hi(1),flx_lo(2):flx_hi(2),NVAR)

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),NQAUX)
    real(rt), intent(in) :: gamma_up(glo(1):ghi(1),glo(2):ghi(2),9)

    ! Local variables
    integer i, j

    ! these will refer to the zone interfaces that we solve the
    ! Riemann problem across
    integer :: imin, imax, jmin, jmax

    real(rt)         :: cl, cr
    type (eos_t) :: eos_state

    if (idir == 1) then
       imin = ilo
       imax = ihi+1
       jmin = jlo
       jmax = jhi
    else
       imin = ilo
       imax = ihi
       jmin = jlo
       jmax = jhi+1
    endif

    ! Solve Riemann problem (godunov state passed back, but only (u,p) saved)
    call gr_HLL(qm, qp, qpd_lo, qpd_hi, &
             qaux, qa_lo, qa_hi, &
             flx, flx_lo, flx_hi, &
             idir, imin, imax, jmin, jmax, 0, 0, 0, &
             [domlo(1), domlo(2), 0], [domhi(1), domhi(2), 0], &
             gamma_up, [glo(1), glo(2), 0], [ghi(1), ghi(2), 0])


  end subroutine cmpflx

end module riemann_module
