module riemann_module

  use bl_types
  use bl_constants_module
  use riemann_util_module
  use meth_params_module, only : NQ, NQAUX, NVAR, QRHO, QU, QV, QW, &
                                 QPRES, QGAME, QREINT, QFS, &
                                 QFX, URHO, UMX, UMY, UMZ, UEDEN, UEINT, &
                                 UFS, UFX, &
                                 QC, QCSML, QGAMC

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
                    qaux, qa_lo, qa_hi, &
                    idir, lo, hi, domlo, domhi)

    use actual_riemann_module
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

end module riemann_module
