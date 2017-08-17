module bc_ext_fill_module

  use bl_constants_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, &
                                 UEDEN, UEINT, UFS, UTEMP, const_grav, &
                                 hse_zero_vels, hse_interp_temp, hse_reflect_vels, &
                                 xl_ext, xr_ext, yl_ext, yr_ext, EXT_HSE, EXT_INTERP
  use interpolate_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  include 'AMReX_bc_types.fi'

  private

  public :: ext_fill, ext_denfill

contains

  ! this module contains different routines for filling the
  ! hydrodynamics boundary conditions

  ! NOTE: the hydrostatic boundary conditions here rely on
  ! constant gravity

  subroutine ext_fill(adv, adv_l1, adv_l2, adv_h1, adv_h2, &
                      domlo, domhi, delta, xlo, time, bc) &
                      bind(C, name="ext_fill")

    use prob_params_module, only : problo
    use network, only: nspec
    use model_parser_module

    use amrex_fort_module, only : rt => amrex_real
    integer adv_l1, adv_l2, adv_h1, adv_h2
    integer bc(2,2,*)
    integer domlo(2), domhi(2)
    real(rt)         delta(2), xlo(2), time
    real(rt)         adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)

    integer i, j, q, n, iter, m
    real(rt)         y
    real(rt)         :: dens_above, dens_base
    real(rt)         :: X_zone(nspec), dens_zone

    integer, parameter :: MAX_ITER = 100
    real(rt)        , parameter :: TOL = 1.e-8_rt

    do n = 1, NVAR

       ! XLO
       if (bc(1,1,n) == EXT_DIR .and. xl_ext == EXT_HSE .and. adv_l1 < domlo(1)) then
          call bl_error("ERROR: HSE boundaries not implemented for -X")
       end if

       ! XHI
       if (bc(1,2,n) == EXT_DIR .and. xr_ext == EXT_HSE .and. adv_h1 > domhi(1)) then
          call bl_error("ERROR: HSE boundaries not implemented for +X")
       end if

       ! YLO
       if (bc(2,1,n) == EXT_DIR .and. adv_l2 < domlo(2)) then

          if (yl_ext == EXT_HSE) then

             ! we will fill all the variables when we consider URHO
             if (n == URHO) then

                do i = adv_l1, adv_h1

                   ! we are integrating along a column at constant i.
                   ! Make sure that our starting state is well-defined
                   dens_above = adv(i,domlo(2),URHO)

                      X_zone(:) = adv(i,domlo(2),UFS:UFS-1+nspec)/dens_above

                   ! keep track of the density at the base of the domain
                   dens_base = dens_above

                   ! integrate downward
                   do j = domlo(2)-1, adv_l2, -1
                      y = problo(2) + delta(2)*(dble(j) + HALF)

                      ! HSE integration to get density, pressure

                      ! initial guesses
                      dens_zone = dens_above

                      ! velocity
                      if (hse_zero_vels == 1) then

                         ! zero normal momentum causes pi waves to pass through
                         adv(i,j,UMY) = ZERO

                         ! zero transverse momentum
                         adv(i,j,UMX) = ZERO
                         adv(i,j,UMZ) = ZERO
                      else

                         if (hse_reflect_vels == 1) then
                            adv(i,j,UMX) = -dens_zone*(adv(i,domlo(2),UMX)/dens_base)
                            adv(i,j,UMY) = -dens_zone*(adv(i,domlo(2),UMY)/dens_base)
                            adv(i,j,UMZ) = -dens_zone*(adv(i,domlo(2),UMZ)/dens_base)
                         else
                            adv(i,j,UMX) = dens_zone*(adv(i,domlo(2),UMX)/dens_base)
                            adv(i,j,UMY) = dens_zone*(adv(i,domlo(2),UMY)/dens_base)
                            adv(i,j,UMZ) = dens_zone*(adv(i,domlo(2),UMZ)/dens_base)
                         endif
                      endif

                              ! store the final state
                      adv(i,j,URHO) = dens_zone
                      adv(i,j,UFS:UFS-1+nspec) = dens_zone*X_zone(:)

                      ! for the next zone
                      dens_above = dens_zone

                   enddo
                enddo

             endif  ! n == URHO

          elseif (yl_ext == EXT_INTERP) then

             do j = domlo(2)-1, adv_l2, -1
                y = problo(2) + delta(2)*(dble(j) + HALF)

                do i = adv_l1, adv_h1

                   ! set all the variables even though we're testing on URHO
                   if (n == URHO) then

                      dens_zone = interpolate(y,npts_model,model_r, &
                                              model_state(:,idens_model))

                                      do q = 1, nspec
                         X_zone(q) = interpolate(y,npts_model,model_r, &
                                                 model_state(:,ispec_model-1+q))
                      enddo

                      ! extrap normal momentum
                      adv(i,j,UMY) = min(ZERO, adv(i,domhi(2),UMY))

                      ! zero transverse momentum
                      adv(i,j,UMX) = ZERO
                      adv(i,j,UMZ) = ZERO

                      adv(i,j,URHO) = dens_zone
                      adv(i,j,UFS:UFS-1+nspec) = dens_zone*X_zone(:)
                   endif

                enddo
             enddo
          endif  ! yl_ext check


       endif


       ! YHI
       if (bc(2,2,n) == EXT_DIR .and. adv_h2 > domhi(2)) then

          if (yr_ext == EXT_HSE) then
             call bl_error("ERROR: HSE boundaries not implemented for +Y")

          elseif (yr_ext == EXT_INTERP) then
             ! interpolate thermodynamics from initial model

             do j = domhi(2)+1, adv_h2
                y = problo(2) + delta(2)*(dble(j) + HALF)

                do i = adv_l1, adv_h1

                   ! set all the variables even though we're testing on URHO
                   if (n == URHO) then

                      dens_zone = interpolate(y,npts_model,model_r, &
                                              model_state(:,idens_model))

                                      do q = 1, nspec
                         X_zone(q) = interpolate(y,npts_model,model_r, &
                                                 model_state(:,ispec_model-1+q))
                      enddo


                      ! extrap normal momentum
                      adv(i,j,UMY) = max(ZERO, adv(i,domhi(2),UMY))

                      ! zero transverse momentum
                      adv(i,j,UMX) = ZERO
                      adv(i,j,UMZ) = ZERO

                      adv(i,j,URHO) = dens_zone
                      adv(i,j,UFS:UFS-1+nspec) = dens_zone*X_zone(:)

                   endif

                enddo
             enddo
          endif  ! yr_ext check

       endif

    enddo

  end subroutine ext_fill


  subroutine ext_denfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                         domlo,domhi,delta,xlo,time,bc) &
                         bind(C, name="ext_denfill")

    use prob_params_module, only : problo
    use interpolate_module
    use model_parser_module
    use bl_error_module

    use amrex_fort_module, only : rt => amrex_real
    integer adv_l1,adv_l2,adv_h1,adv_h2
    integer bc(2,2,*)
    integer domlo(2), domhi(2)
    real(rt)         delta(2), xlo(2), time
    real(rt)         adv(adv_l1:adv_h1,adv_l2:adv_h2)

    integer i,j
    real(rt)         y

    ! Note: this function should not be needed, technically, but is
    ! provided to filpatch because there are many times in the algorithm
    ! when just the density is needed.  We try to rig up the filling so
    ! that the same function is called here and in hypfill where all the
    ! states are filled.

    call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,bc)

    ! XLO
    if ( bc(1,1,1) == EXT_DIR .and. adv_l1 < domlo(1)) then
       call bl_error("We shoundn't be here (xlo denfill)")
    end if

    ! XHI
    if ( bc(1,2,1) == EXT_DIR .and. adv_h1 > domhi(1)) then
       call bl_error("We shoundn't be here (xlo denfill)")
    endif

    ! YLO
    if ( bc(2,1,1) == EXT_DIR .and. adv_l2 < domlo(2)) then
       do j = adv_l2, domlo(2)-1
          y = problo(2) + delta(2)*(dble(j) + HALF)
          do i = adv_l1, adv_h1
             adv(i,j) = interpolate(y,npts_model,model_r,model_state(:,idens_model))
          end do
       end do
    end if

    ! YHI
    if ( bc(2,2,1) == EXT_DIR .and. adv_h2 > domhi(2)) then
       do j = domhi(2)+1, adv_h2
          y = problo(2) + delta(2)*(dble(j)+ HALF)
          do i = adv_l1, adv_h1
             adv(i,j) = interpolate(y,npts_model,model_r,model_state(:,idens_model))
          end do
       end do
    end if

  end subroutine ext_denfill

end module bc_ext_fill_module
