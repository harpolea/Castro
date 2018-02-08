! problem-specific Fortran derive routines go here
subroutine ca_derheight(dnull,k_lo,k_hi,nk, &
                      dat,d_lo,d_hi,nc, &
                      lo,hi,domlo,domhi,delta, &
                      xlo,time,dt,bc,level,grid_no) &
                      bind(C, name="ca_derheight")
  !
  ! This routine is used by particle_count.  Yes it does nothing.
  !
  use amrex_fort_module, only : rt => amrex_real
  use meth_params_module, only : NQ, NQAUX, QRHO
  use advection_util_module, only : swectoprim
  implicit none

  integer          :: lo(3), hi(3)
  integer          :: k_lo(3), k_hi(3), nk
  integer          :: d_lo(3), d_hi(3), nc
  integer          :: domlo(3), domhi(3)
  integer          :: bc(3,2,nc)
  real(rt)         :: delta(3), xlo(3), time, dt
  real(rt)         :: dnull(k_lo(1):k_hi(1),k_lo(2):k_hi(2),k_lo(3):k_hi(3),nk)
  real(rt)         :: dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
  integer          :: level, grid_no


  real(rt)         :: q(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),NQ)
  real(rt)         :: qaux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),NQAUX)

  call swectoprim(lo, hi, dat, d_lo, d_hi, q, lo, hi, qaux, lo, hi)

  dnull(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1) = q(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),QRHO)

end subroutine ca_derheight
