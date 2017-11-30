module timestep_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  public

contains

  ! Courant-condition limited timestep

  subroutine ca_estdt(lo,hi,u,u_lo,u_hi,dx,dt) bind(C, name="ca_estdt")

    use network, only: nspec, naux
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ
    use prob_params_module, only: dim
    use bl_constants_module
    use amrex_fort_module, only : rt => amrex_real
    use probdata_module, only : g

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: u_lo(3), u_hi(3)
    real(rt), intent(in) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
    real(rt), intent(in) :: dx(3)
    real(rt), intent(inout) :: dt

    real(rt)         :: rhoInv, ux, uy, uz, c, dt1, dt2, dt3, dt_tmp
    integer          :: i, j, k

    ! Call EOS for the purpose of computing sound speed

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             rhoInv = ONE / u(i,j,k,URHO)

             ! Compute velocity and then calculate CFL timestep.

             ux = u(i,j,k,UMX) * rhoInv
             uy = u(i,j,k,UMY) * rhoInv
             uz = u(i,j,k,UMZ) * rhoInv

             c = sqrt(g * u(i,j,k,URHO)) ! sound speed is sqrt(gh)

             dt1 = dx(1)/(c + abs(ux))
             if (dim >= 2) then
                dt2 = dx(2)/(c + abs(uy))
             else
                dt2 = dt1
             endif
             if (dim == 3) then
                dt3 = dx(3)/(c + abs(uz))
             else
                dt3 = dt1
             endif

            ! method of lines constraint is tougher
            dt_tmp = ONE/dt1
            if (dim >= 2) then
               dt_tmp = dt_tmp + ONE/dt2
            endif
            if (dim == 3) then
               dt_tmp = dt_tmp + ONE/dt3
            endif
            dt = min(dt, ONE/dt_tmp)

          enddo
       enddo
    enddo

  end subroutine ca_estdt

  ! Reactions-limited timestep


  ! Check whether the last timestep violated any of our stability criteria.
  ! If so, suggest a new timestep which would not.

  subroutine ca_check_timestep(s_old, so_lo, so_hi, &
                               s_new, sn_lo, sn_hi, &
                               lo, hi, &
                               dx, dt_old, dt_new) &
                               bind(C, name="ca_check_timestep")

    use bl_constants_module, only: HALF, ONE
    use meth_params_module, only: NVAR, URHO, UMX, UMZ, cfl
    use prob_params_module, only: dim
    use probdata_module, only : g
    use network, only: nspec, naux
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: so_lo(3), so_hi(3)
    integer, intent(in) :: sn_lo(3), sn_hi(3)
    real(rt), intent(in) :: s_old(so_lo(1):so_hi(1),so_lo(2):so_hi(2),so_lo(3):so_hi(3),NVAR)
    real(rt), intent(in) :: s_new(sn_lo(1):sn_hi(1),sn_lo(2):sn_hi(2),sn_lo(3):sn_hi(3),NVAR)
    real(rt), intent(in) :: dx(3), dt_old
    real(rt), intent(inout) :: dt_new

    integer          :: i, j, k
    real(rt)         :: tau_CFL
    real(rt)         :: h, v(3), c

    real(rt), parameter :: derivative_floor = 1.e-50_rt

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! CFL hydrodynamic stability criterion

             ! If the timestep violated (v+c) * dt / dx > 1,
             ! suggest a new timestep such that (v+c) * dt / dx <= CFL,
             ! where CFL is the user's chosen timestep constraint.
             ! We don't use the CFL choice in the test for violation
             ! because if we did, then even a small increase in velocity
             ! over the timestep would be enough to trigger a retry
             ! (and empirically we have found that this can happen
             ! basically every timestep, which can greatly increase
             ! the cost of a simulation for not much benefit);
             ! but these types of issues are exactly why a user picks a
             ! safe buffer value like CFL = 0.8 or CFL = 0.5. We only
             ! want to trigger a retry if the timestep strongly violated
             ! the stability criterion.

            v = HALF * (s_old(i,j,k,UMX:UMZ) / s_old(i,j,k,URHO) + &
                s_new(i,j,k,UMX:UMZ) / s_new(i,j,k,URHO))
            h = HALF * (s_old(i,j,k,URHO) + s_new(i,j,k,URHO))

            c = sqrt(g * h) ! sound speed is sqrt(gh)

            tau_CFL = minval(dx(1:dim) / (c + abs(v(1:dim))))

            if (dt_old > tau_CFL) then
               dt_new = min(dt_new, cfl * tau_CFL)
            endif
          enddo
       enddo
    enddo

  end subroutine ca_check_timestep

end module timestep_module
