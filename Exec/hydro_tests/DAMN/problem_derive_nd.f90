! problem-specific Fortran derive routines go here

  subroutine ca_derheight(h,h_lo,h_hi,nh, &
                        u,d_lo,d_hi,nc, &
                        lo,hi,domlo,domhi,dx, &
                        xlo,time,dt,bc,level,grid_no) &
                        bind(C, name="ca_derheight")
    !
    ! This routine is used by particle_count.  Yes it does nothing.
    !
    use amrex_fort_module, only : rt => amrex_real
    use probdata_module, only: swe_to_comp_level, g
    use network, only: nspec, naux
    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_re
    use meth_params_module, only: URHO, UEINT, UTEMP, UFS, UFX
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
    integer :: i,j,k
    type (eos_t) :: eos_state

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (level <= swe_to_comp_level) then

                 h(i,j,k,1) = u(i,j,k,URHO)

             else

                 eos_state % rho  = u(i,j,k,URHO)
                 eos_state % T    = u(i,j,k,UTEMP)
                 eos_state % e    = u(i,j,k,UEINT) / u(i,j,k,URHO)
                 eos_state % xn   = u(i,j,k,UFS:UFS+nspec-1) / u(i,j,k,URHO)
                 eos_state % aux  = u(i,j,k,UFX:UFX+naux-1) / u(i,j,k,URHO)

                 call eos(eos_input_re, eos_state)

                 p = eos_state % p
                 xx = xlo(1) + dx(1)*dble(i-lo(1)+0.5d0)

                 h(i,j,k,1) = sqrt(2.0d0 * p / g) + xx
             endif
          enddo
       enddo
    enddo

end subroutine ca_derheight
