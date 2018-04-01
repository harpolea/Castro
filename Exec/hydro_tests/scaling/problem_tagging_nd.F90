module problem_tagging_module

  use amrex_fort_module, only : rt => amrex_real
  use tagging_module, only : speederr
  implicit none

  public

contains

  ! This is a template routine for users to set their own tags based on the state.
  ! It will be overwritten by having a copy of this file in the user's problem setup.

  subroutine set_problem_tags(tag,tag_lo,tag_hi, &
                              state,state_lo,state_hi, &
                              set,clear,&
                              lo,hi,&
                              dx,problo,time,level,xlo) &
                              bind(C, name="set_problem_tags")

    use meth_params_module, only : NVAR, NQ, QU, QW, NQAUX
    use probdata_module, only: swe_to_comp_level, damn_rad
    use advection_util_module, only: compctoprim, swectoprim
    use prob_params_module, only : center
    implicit none

    integer, intent(in)          :: lo(3),hi(3)
    integer, intent(in)          :: state_lo(3),state_hi(3)
    integer, intent(in)          :: tag_lo(3),tag_hi(3)
    real(rt), intent(in)         :: state(state_lo(1):state_hi(1), &
                        state_lo(2):state_hi(2), &
                        state_lo(3):state_hi(3),NVAR)
    integer, intent(inout)         :: tag(tag_lo(1):tag_hi(1),tag_lo(2):tag_hi(2),tag_lo(3):tag_hi(3))
    real(rt), intent(in)         :: problo(3),dx(3),time,xlo(3)
    integer, intent(in)          :: level,set,clear

    integer :: i, j, k
    real(rt) :: yy, zz, r
    real(rt) :: q(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), NQ)
    real(rt) :: qaux(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),NQAUX)
    real(rt) :: speed, press

    if (level > swe_to_comp_level) then
        call compctoprim(lo, hi, state, state_lo, state_hi, q, lo, hi, qaux, lo, hi, xlo, dx, .false.)
    else
        call swectoprim(lo, hi, state, state_lo, state_hi, q, lo, hi, qaux, lo, hi, .false.)
    endif

    ! do k = lo(3), hi(3)
    !    do j = lo(2), hi(2)
    !       do i = lo(1), hi(1)
    !           speed = sqrt(sum(q(i,j,k,QU:QW)**2))
    !           ! write(*,*) "speed = ", speed
    !           if (speed .ge. speederr) then
    !               tag(i,j,k) = set
    !           endif
    !       enddo
    !    enddo
    ! enddo

    ! do k = lo(3), hi(3)
    !    do j = lo(2), hi(2)
    !       yy = xlo(2) + dx(2)*dble(j-lo(2)+0.5_rt)
    !       do i = lo(1), hi(1)
    !           if (abs(yy - center(2)) .le. speederr) then
    !               tag(i,j,k) = set
    !           endif
    !       enddo
    !    enddo
    ! enddo

    !rad dam
    ! do k = lo(3), hi(3)
    !     do j = lo(2), hi(2)
    !        yy = xlo(2) + dx(2)*dble(j-lo(2)+0.5_rt)
    !        r = sqrt((yy - center(2))**2)
    !        do i = lo(1), hi(1)
    !           if (r .le. 0.2_rt * damn_rad) then
    !               tag(i,j,k) = set
    !           endif
    !       enddo
    !    enddo
    ! enddo

    ! If it's the swe_to_comp_level or coarser, need to make sure that
    ! any refinement happens uniformly in vertical direction
    ! This only works on patches booo
    if (level <= swe_to_comp_level) then
        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              if (any(tag(lo(1):hi(1),j,k) .eq. set)) then
                  tag(lo(1):hi(1),j,k) = set
              endif
           enddo
        enddo
    endif

  end subroutine set_problem_tags

end module problem_tagging_module
