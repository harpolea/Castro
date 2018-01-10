! problem-specific Fortran derive routines go here

subroutine ca_derheight(h,h_lo,h_hi,nh, &
                    u,d_lo,d_hi,nc, &
                    lo,hi,domlo,domhi,dx, &
                    xlo,time,dt,bc,level,grid_no) &
                    bind(C, name="ca_derheight")
    !
    ! Calculate height
    !
    use amrex_fort_module, only : rt => amrex_real
    use probdata_module, only: swe_to_comp_level
    use meth_params_module, only: QRHO, NQ, NQAUX, NVAR
    use advection_util_module, only: swectoprim
    use riemann_module, only: comp_to_swe
    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: h_lo(3), h_hi(3), nh
    integer, intent(in) :: d_lo(3), d_hi(3), nc
    integer, intent(in) :: domlo(3), domhi(3)
    integer, intent(in) :: bc(3,2,nc)
    real(rt), intent(in) :: dx(3), xlo(3), time, dt
    real(rt), intent(inout) :: h(h_lo(1):h_hi(1),h_lo(2):h_hi(2),h_lo(3):h_hi(3),nh)
    real(rt), intent(in) ::    u(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer, intent(in) :: level, grid_no

    real(rt) :: xx, p
    real(rt) :: q(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), NQ)
    real(rt) :: qaux(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),NQAUX)
    real(rt) :: swe(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), NVAR)
    integer :: i,j,k

    if (level > swe_to_comp_level) then
        call comp_to_swe(swe, lo, hi, u(:,:,:,1:NVAR), d_lo, d_hi, lo, hi, xlo, dx, .false.)

        call swectoprim(lo, hi, swe, lo, hi, q, lo, hi, qaux, lo, hi, .false.)
    else
        call swectoprim(lo, hi, u(:,:,:,1:NVAR), d_lo, d_hi, q, lo, hi, qaux, lo, hi, .false.)
    endif

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
              h(i,j,k,1) = q(i,j,k,QRHO)
          enddo
       enddo
    enddo

end subroutine ca_derheight

subroutine ca_derprimv(v,v_lo,v_hi,nv, &
                    u,d_lo,d_hi,nc, &
                    lo,hi,domlo,domhi,dx, &
                    xlo,time,dt,bc,level,grid_no) &
                    bind(C, name="ca_derprimv")
    !
    ! Primitive velocity
    !
    use amrex_fort_module, only : rt => amrex_real
    use probdata_module, only: swe_to_comp_level
    use meth_params_module, only: QU, QW, NQ, NQAUX, NVAR
    use advection_util_module, only: swectoprim, compctoprim
    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: v_lo(3), v_hi(3), nv
    integer, intent(in) :: d_lo(3), d_hi(3), nc
    integer, intent(in) :: domlo(3), domhi(3)
    integer, intent(in) :: bc(3,2,nc)
    real(rt), intent(in) :: dx(3), xlo(3), time, dt
    real(rt), intent(inout) :: v(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
    real(rt), intent(in) ::    u(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer, intent(in) :: level, grid_no

    real(rt) :: xx, p
    real(rt) :: q(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), NQ)
    real(rt) :: qaux(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),NQAUX)
    integer :: i,j,k

    if (level > swe_to_comp_level) then
        call compctoprim(lo, hi, u(:,:,:,1:NVAR), d_lo, d_hi, q, lo, hi, qaux, lo, hi, xlo, dx, .false.)
    else
        call swectoprim(lo, hi, u(:,:,:,1:NVAR), d_lo, d_hi, q, lo, hi, qaux, lo, hi, .false.)
    endif

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
              v(i,j,k,1:3) = q(i,j,k,QU:QW)
          enddo
       enddo
    enddo

end subroutine ca_derprimv

subroutine ca_derprimrho(rho,r_lo,r_hi,nr, &
                    u,d_lo,d_hi,nc, &
                    lo,hi,domlo,domhi,dx, &
                    xlo,time,dt,bc,level,grid_no) &
                    bind(C, name="ca_derprimrho")
    !
    ! Primitive density
    !
    use amrex_fort_module, only : rt => amrex_real
    use probdata_module, only: swe_to_comp_level
    use meth_params_module, only: QRHO, NQ, NQAUX, NVAR
    use advection_util_module, only: swectoprim, compctoprim
    use riemann_module, only: swe_to_comp
    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: r_lo(3), r_hi(3), nr
    integer, intent(in) :: d_lo(3), d_hi(3), nc
    integer, intent(in) :: domlo(3), domhi(3)
    integer, intent(in) :: bc(3,2,nc)
    real(rt), intent(in) :: dx(3), xlo(3), time, dt
    real(rt), intent(inout) :: rho(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3),nr)
    real(rt), intent(in) ::    u(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer, intent(in) :: level, grid_no

    real(rt) :: xx, p
    real(rt) :: q(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), NQ)
    real(rt) :: qaux(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),NQAUX)
    real(rt) :: comp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), NVAR)
    integer :: i,j,k

    if (level > swe_to_comp_level) then
        call compctoprim(lo, hi, u(:,:,:,1:NVAR), d_lo, d_hi, q, lo, hi, qaux, lo, hi, xlo, dx, .false.)
    else
        call swe_to_comp(u(:,:,:,1:NVAR), d_lo, d_hi, comp, lo, hi, lo, hi, dx, xlo, .false.)

        call compctoprim(lo, hi, comp(:,:,:,1:NVAR), lo, hi, q, lo, hi, qaux, lo, hi, xlo, dx, .false.)
    endif

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
              rho(i,j,k,1) = q(i,j,k,QRHO)
          enddo
       enddo
    enddo

end subroutine ca_derprimrho

subroutine ca_derW(w,w_lo,w_hi,nw, &
                    u,d_lo,d_hi,nc, &
                    lo,hi,domlo,domhi,dx, &
                    xlo,time,dt,bc,level,grid_no) &
                    bind(C, name="ca_derW")
    !
    ! Lorentz factor
    !
    use amrex_fort_module, only : rt => amrex_real
    use probdata_module, only: swe_to_comp_level
    use meth_params_module, only: QU, QW, NQ, NQAUX, NVAR
    use advection_util_module, only: swectoprim, compctoprim
    use metric_module, only: calculate_scalar_W
    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: w_lo(3), w_hi(3), nw
    integer, intent(in) :: d_lo(3), d_hi(3), nc
    integer, intent(in) :: domlo(3), domhi(3)
    integer, intent(in) :: bc(3,2,nc)
    real(rt), intent(in) :: dx(3), xlo(3), time, dt
    real(rt), intent(inout) :: w(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3),nw)
    real(rt), intent(in) ::    u(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer, intent(in) :: level, grid_no

    real(rt) :: xx, p
    real(rt) :: q(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), NQ)
    real(rt) :: qaux(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),NQAUX)
    integer :: i,j,k

    if (level > swe_to_comp_level) then
        call compctoprim(lo, hi, u(:,:,:,1:NVAR), d_lo, d_hi, q, lo, hi, qaux, lo, hi, xlo, dx, .false.)
    else
        call swectoprim(lo, hi, u(:,:,:,1:NVAR), d_lo, d_hi, q, lo, hi, qaux, lo, hi, .false.)
    endif

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
              call calculate_scalar_W(q(i,j,k,QU:QW), w(i,j,k,1))
          enddo
       enddo
    enddo

end subroutine ca_derW

subroutine ca_dereint(e,e_lo,e_hi,ncomp_e, &
                       u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                       domhi,dx,xlo,time,dt,bc,level,grid_no) &
                       bind(C, name="ca_dereint")

  use meth_params_module, only: QRHO, QREINT, NQ, NQAUX, NVAR
  use advection_util_module, only: compctoprim
  use riemann_module, only: swe_to_comp
  use probdata_module, only: swe_to_comp_level

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: e_lo(3), e_hi(3), ncomp_e
  integer, intent(in) :: u_lo(3), u_hi(3), ncomp_u
  integer, intent(in) :: domlo(3), domhi(3)
  real(rt), intent(inout) :: e(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),ncomp_e)
  real(rt), intent(in) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  real(rt), intent(in) :: dx(3), xlo(3), time, dt
  integer, intent(in) :: bc(3,2,ncomp_u), level, grid_no
  
  real(rt) :: q(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), NQ)
  real(rt) :: qaux(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),NQAUX)
  real(rt) :: comp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), NVAR)
  integer :: i,j,k

  if (level > swe_to_comp_level) then
      call compctoprim(lo, hi, u(:,:,:,1:NVAR), u_lo, u_hi, q, lo, hi, qaux, lo, hi, xlo, dx, .false.)
  else
      call swe_to_comp(u(:,:,:,1:NVAR), u_lo, u_hi, comp, lo, hi, lo, hi, dx, xlo, .false.)

      call compctoprim(lo, hi, comp(:,:,:,1:NVAR), lo, hi, q, lo, hi, qaux, lo, hi, xlo, dx, .false.)
  endif

  !
  ! Compute internal energy from (rho e).
  !
  do k = lo(3),hi(3)
     do j = lo(2),hi(2)
        do i = lo(1),hi(1)
           e(i,j,k,1) = q(i,j,k,QREINT) / q(i,j,k,QRHO)
        enddo
     enddo
  enddo

end subroutine ca_dereint
