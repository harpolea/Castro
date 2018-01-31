module bc_fill_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains

  subroutine ca_hypfill(adv,adv_lo,adv_hi, &
                        domlo,domhi,delta,xlo,time,bc) &
                        bind(C, name="ca_hypfill")

    use probdata_module
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UTEMP
    use interpolate_module
    use eos_module, only : eos
    use eos_type_module, only: eos_input_rt, eos_t
    use network, only: nspec
    !use model_parser_module
    use prob_params_module, only: dim

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer, intent(in) :: adv_lo(3),adv_hi(3)
    integer, intent(in) :: bc(dim,2,*)
    integer, intent(in) :: domlo(3), domhi(3)
    real(rt), intent(in)         :: delta(3), xlo(3), time
    real(rt), intent(inout)         :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NVAR)

    integer :: i,j,k,q,n,iter, npts
    real(rt)         :: x
    real(rt)         :: xs(domlo(1)-2:adv_hi(1)+1)
    real(rt)         :: pres_above,p_want,pres_zone, A
    real(rt)         :: drho,dpdr,temp_zone,eint,X_zone(nspec),dens_zone

    integer, parameter :: MAX_ITER = 100
    real(rt)        , parameter :: TOL = 1.e-8_rt
    logical :: converged_hse

    type (eos_t) :: eos_state

    do n = 1,NVAR
       call filcc_nd(adv(adv_lo(1),adv_lo(2),adv_lo(3),n),adv_lo,adv_hi, &
                  domlo,domhi,delta,xlo,bc(2,1,n))

       call filcc_nd(adv(adv_lo(1),adv_lo(2),adv_lo(3),n),adv_lo,adv_hi, &
                 domlo,domhi,delta,xlo,bc(3,1,n))
    enddo

    do n = 1, NVAR

       !        YLO
       if ( bc(2,1,n).eq.EXT_DIR .and. adv_lo(2).lt.domlo(2)) then

          ! we are FOEXTRAP in y -- we should never get here
          call bl_error("ERROR: invalid BC in Prob_nd.f90")

       end if

       !        YHI
       if ( bc(2,2,n).eq.EXT_DIR .and. adv_hi(2).gt.domhi(2)) then

          ! we are FOEXTRAP in y -- we should never get here
          call bl_error("ERROR: invalid BC in Prob_nd.f90")

       end if

       !        ZLO
       if ( bc(3,1,n).eq.EXT_DIR .and. adv_lo(3).lt.domlo(3)) then

          ! we are FOEXTRAP in z -- we should never get here
          call bl_error("ERROR: invalid BC in Prob_nd.f90")

       end if

       !        ZHI
       if ( bc(3,2,n).eq.EXT_DIR .and. adv_hi(3).gt.domhi(3)) then

          ! we are FOEXTRAP in z -- we should never get here
          call bl_error("ERROR: invalid BC in Prob_nd.f90")

       end if

       do i= domlo(1)-2, adv_hi(1)+1
           xs(i) = xlo(1) + delta(1)*dble(i-adv_lo(1))
       enddo
       npts = adv_hi(1) - adv_lo(1) + 1

       !        XLO
       if ( bc(1,1,n).eq.EXT_DIR .and. adv_lo(1).lt.domlo(1)) then

          ! this do loop counts backwards since we want to work downward
          do k = adv_lo(3), adv_hi(3)
              do j = adv_lo(2), adv_hi(2)
                 ! y = xlo(2) + delta(2)*(dble(j-adv_lo(2)) + 0.5e0_rt)

                 do i= domlo(1)-1, adv_lo(1)-1
                     x = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5e0_rt)

                    ! set all the variables even though we're testing on URHO
                    if (n .eq. URHO) then

                      ! HSE integration to get density, pressure

                      ! initial guesses
                      dens_zone = adv(i+1,j,k,URHO)

                      ! temperature and species held constant in BCs
                      temp_zone = adv(i+1,j,k,UTEMP)
                      X_zone(:) = adv(i+1,j,k,UFS:UFS-1+nspec)/adv(i+1,j,k,URHO)

                      ! get pressure in zone above
                      eos_state%rho = adv(i+1,j,k,URHO)
                      eos_state%T = adv(i+1,j,k,UTEMP)
                      eos_state%xn(:) = adv(i+1,j,k,UFS:UFS-1+nspec)/adv(i+1,j,k,URHO)

                      call eos(eos_input_rt, eos_state)

                      eint = eos_state%e
                      pres_above = eos_state%p


                      converged_hse = .FALSE.

                      do iter = 1, MAX_ITER

                         ! pressure needed from HSE
                         p_want = pres_above - &
                              delta(1)*0.5e0_rt*(dens_zone + adv(i+1,j,k,URHO))*g

                         ! pressure from EOS
                         eos_state%rho = dens_zone
                         eos_state%T = temp_zone
                         eos_state%xn(:) = X_zone

                         call eos(eos_input_rt, eos_state)

                         pres_zone = eos_state%p
                         dpdr = eos_state%dpdr
                         eint = eos_state%e

                         ! Newton-Raphson - we want to zero A = p_want - p(rho)
                         A = p_want - pres_zone
                         drho = A/(dpdr + 0.5*delta(1)*g)

                         dens_zone = max(0.9_rt*dens_zone, &
                              min(dens_zone + drho, 1.1_rt*dens_zone))


                         ! convergence?
                         if (abs(drho) < TOL*dens_zone) then
                            converged_hse = .TRUE.
                            exit
                         endif

                      enddo

                      if (.not. converged_hse) call bl_error("ERROR: failure to converge in -Y BC")

                      ! zero gradient velocity
                      adv(i,j,k,UMX) = dens_zone*(adv(domlo(1),j,k,UMX)/adv(domlo(1),j,k,URHO))
                      adv(i,j,k,UMY) = dens_zone*(adv(domlo(1),j,k,UMY)/adv(domlo(1),j,k,URHO))
                      adv(i,j,k,UMZ) = dens_zone*(adv(domlo(1),j,k,UMZ)/adv(domlo(1),j,k,URHO))

                       eos_state%rho = dens_zone
                       eos_state%T = temp_zone
                       eos_state%xn(:) = X_zone

                       call eos(eos_input_rt, eos_state)

                       pres_zone = eos_state%p
                       eint = eos_state%e

                       adv(i,j,k,URHO) = dens_zone
                       adv(i,j,k,UEINT) = dens_zone*eint
                       adv(i,j,k,UEDEN) = dens_zone*eint + &
                            0.5e0_rt*(adv(i,j,k,UMX)**2+adv(i,j,k,UMY)**2+adv(i,j,k,UMZ)**2)/dens_zone
                       adv(i,j,k,UTEMP) = temp_zone
                       adv(i,j,k,UFS:UFS-1+nspec) = dens_zone*X_zone(:)

                    end if

                 end do
              end do
           end do
       end if

       !        XHI
       if ( bc(1,2,n).eq.EXT_DIR .and. adv_hi(1).gt.domhi(1)) then

          do k=adv_lo(3),adv_hi(3)
              do j=adv_lo(2),adv_hi(2)
                 ! y = xlo(2) + delta(2)*(dble(j-adv_l2) + 0.5e0_rt)

                 do i=domhi(1)+1,adv_hi(2)
                     x = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5e0_rt)

                    ! set all the variables even though we're testing on URHO
                    if (n .eq. URHO) then

                       ! dens_zone = interpolate(x,npts_model,model_r, &
                       !      model_state(:,idens_model))
                       !
                       ! temp_zone = interpolate(y,npts_model,model_r, &
                       !      model_state(:,itemp_model))
                       !
                       ! do q = 1, nspec
                       !    X_zone(q) = interpolate(x,npts_model,model_r, &
                       !         model_state(:,ispec_model-1+q))
                       ! enddo

                       dens_zone = interpolate_noparser(x,xlo(1),adv(:,j,k,URHO),npts,i)

                       temp_zone = interpolate_noparser(x,xlo(1),adv(:,j,k,UTEMP),npts,i)

                       do q = 1, nspec
                          X_zone(q) = interpolate_noparser(x,xlo(1),adv(:,j,k,UFS-1+q),npts,i)
                       enddo


                       ! extrap normal momentum
                       adv(i,j,k,UMY) = max(0.e0_rt,adv(domhi(1),j,k,UMY))

                       ! zero transverse momentum
                       adv(i,j,k,UMX) = 0.e0_rt
                       adv(i,j,k,UMZ) = 0.e0_rt

                       eos_state%rho = dens_zone
                       eos_state%T = temp_zone
                       eos_state%xn(:) = X_zone

                       call eos(eos_input_rt, eos_state)

                       pres_zone = eos_state%p
                       eint = eos_state%e

                       adv(i,j,k,URHO) = dens_zone
                       adv(i,j,k,UEINT) = dens_zone*eint
                       adv(i,j,k,UEDEN) = dens_zone*eint + &
                            0.5e0_rt*(adv(i,j,k,UMX)**2+adv(i,j,k,UMY)**2+adv(i,j,k,UMZ)**2)/dens_zone
                       adv(i,j,k,UTEMP) = temp_zone
                       adv(i,j,k,UFS:UFS-1+nspec) = dens_zone*X_zone(:)

                    end if

                 end do
              end do
           end do
       end if

    end do

  end subroutine ca_hypfill



  subroutine ca_denfill(adv,adv_lo,adv_hi, &
                        domlo,domhi,delta,xlo,time,bc) &
                        bind(C, name="ca_denfill")

    use probdata_module
    use interpolate_module
    use model_parser_module
    use bl_error_module
    use prob_params_module, only: dim

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer, intent(in) :: adv_lo(3),adv_hi(3)
    integer, intent(in) :: bc(dim,2,*)
    integer, intent(in) :: domlo(3), domhi(3)
    real(rt), intent(in)         :: delta(3), xlo(3), time
    real(rt), intent(inout)         :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))

    integer :: i,j,k,npts
    real(rt)         :: x, xs(adv_lo(1)-1: adv_hi(1)+1)

    ! Note: this function should not be needed, technically, but is
    ! provided to filpatch because there are many times in the algorithm
    ! when just the density is needed.  We try to rig up the filling so
    ! that the same function is called here and in hypfill where all the
    ! states are filled.

    call filcc_nd(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,bc)

    !     YLO
    if ( bc(2,1,1).eq.EXT_DIR .and. adv_lo(2).lt.domlo(2)) then
       call bl_error("We shoundn't be here (ylo denfill)")
    end if

    !     YHI
    if ( bc(2,2,1).eq.EXT_DIR .and. adv_hi(2).gt.domhi(2)) then
       call bl_error("We shoundn't be here (ylo denfill)")
    endif

    !     ZLO
    if ( bc(3,1,1).eq.EXT_DIR .and. adv_lo(3).lt.domlo(3)) then
       call bl_error("We shoundn't be here (zlo denfill)")
    end if

    !     ZHI
    if ( bc(3,2,1).eq.EXT_DIR .and. adv_hi(3).gt.domhi(3)) then
       call bl_error("We shoundn't be here (zlo denfill)")
    endif

    do i= adv_lo(1)-1, adv_hi(1)+1
        xs(i) = xlo(1) + delta(1)*dble(i-adv_lo(1))
    enddo
    npts = adv_hi(1) - adv_lo(1) + 1


    !     XLO
    if ( bc(1,1,1).eq.EXT_DIR .and. adv_lo(1).lt.domlo(1)) then
       do k=adv_lo(3),adv_hi(3)
           do j=adv_lo(2),adv_hi(2)
              do i=adv_lo(1),domlo(1)-1
                 x = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5e0_rt)
                 adv(i,j,k) = interpolate_noparser(x,xlo(1),adv(:,j,k),npts,i)
                 !interpolate(x,npts_model,model_r,model_state(:,idens_model))
              end do
           end do
        end do
    end if

    !     XHI
    if ( bc(1,2,1).eq.EXT_DIR .and. adv_hi(1).gt.domhi(1)) then
       do k=adv_lo(3),adv_hi(3)
           do j=adv_lo(2),adv_hi(2)
              do i=domhi(1)+1,adv_hi(1)
                 x = xlo(1) + delta(1)*(dble(i-adv_lo(1))+ 0.5e0_rt)
                 adv(i,j,k) = interpolate_noparser(x,xlo(1),adv(:,j,k),npts,i)
                 !interpolate(x,npts_model,model_r,model_state(:,idens_model))
              end do
           end do
        end do
    end if

  end subroutine ca_denfill


  !
  ! subroutine ca_gravxfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
  !                         domlo,domhi,delta,xlo,time,bc) &
  !                         bind(C, name="ca_gravxfill")
  !
  !   use probdata_module
  !
  !   use amrex_fort_module, only : rt => amrex_real
  !   implicit none
  !
  !   include 'AMReX_bc_types.fi'
  !
  !   integer :: grav_l1,grav_l2,grav_h1,grav_h2
  !   integer :: bc(2,2,*)
  !   integer :: domlo(2), domhi(2)
  !   real(rt)         delta(2), xlo(2), time
  !   real(rt)         grav(grav_l1:grav_h1,grav_l2:grav_h2)
  !
  !   call filcc_nd(grav,grav_l1,grav_l2,grav_h1,grav_h2,domlo,domhi,delta,xlo,bc)
  !
  ! end subroutine ca_gravxfill
  !
  !
  !
  ! subroutine ca_gravyfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
  !                         domlo,domhi,delta,xlo,time,bc) &
  !                         bind(C, name="ca_gravyfill")
  !
  !   use probdata_module
  !
  !   use amrex_fort_module, only : rt => amrex_real
  !   implicit none
  !
  !   include 'AMReX_bc_types.fi'
  !
  !   integer :: grav_l1,grav_l2,grav_h1,grav_h2
  !   integer :: bc(2,2,*)
  !   integer :: domlo(2), domhi(2)
  !   real(rt)         delta(2), xlo(2), time
  !   real(rt)         grav(grav_l1:grav_h1,grav_l2:grav_h2)
  !
  !   call filcc_nd(grav,grav_l1,grav_l2,grav_h1,grav_h2,domlo,domhi,delta,xlo,bc)
  !
  ! end subroutine ca_gravyfill
  !
  !
  !
  ! subroutine ca_gravzfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
  !                         domlo,domhi,delta,xlo,time,bc) &
  !                         bind(C, name="ca_gravzfill")
  !
  !   use probdata_module
  !
  !   use amrex_fort_module, only : rt => amrex_real
  !   implicit none
  !
  !   include 'AMReX_bc_types.fi'
  !
  !   integer :: grav_l1,grav_l2,grav_h1,grav_h2
  !   integer :: bc(2,2,*)
  !   integer :: domlo(2), domhi(2)
  !   real(rt)         delta(2), xlo(2), time
  !   real(rt)         grav(grav_l1:grav_h1,grav_l2:grav_h2)
  !
  !   call filcc_nd(grav,grav_l1,grav_l2,grav_h1,grav_h2,domlo,domhi,delta,xlo,bc)
  !
  ! end subroutine ca_gravzfill
  !
  !
  !
  ! subroutine ca_reactfill(react,react_l1,react_l2, &
  !                         react_h1,react_h2,domlo,domhi,delta,xlo,time,bc) &
  !                         bind(C, name="ca_reactfill")
  !
  !   use probdata_module
  !
  !   use amrex_fort_module, only : rt => amrex_real
  !   implicit none
  !
  !   include 'AMReX_bc_types.fi'
  !
  !   integer :: react_l1,react_l2,react_h1,react_h2
  !   integer :: bc(2,2,*)
  !   integer :: domlo(2), domhi(2)
  !   real(rt)         delta(2), xlo(2), time
  !   real(rt)         react(react_l1:react_h1,react_l2:react_h2)
  !
  !   call filcc_nd(react,react_l1,react_l2,react_h1,react_h2,domlo,domhi,delta,xlo,bc)
  !
  ! end subroutine ca_reactfill
  !
  !
  ! subroutine ca_phigravfill(phi,phi_l1,phi_l2, &
  !                           phi_h1,phi_h2,domlo,domhi,delta,xlo,time,bc) &
  !                           bind(C, name="ca_phigravfill")
  !
  !   use amrex_fort_module, only : rt => amrex_real
  !   implicit none
  !
  !   include 'AMReX_bc_types.fi'
  !
  !   integer          :: phi_l1,phi_l2,phi_h1,phi_h2
  !   integer          :: bc(2,2,*)
  !   integer          :: domlo(2), domhi(2)
  !   real(rt)         :: delta(2), xlo(2), time
  !   real(rt)         :: phi(phi_l1:phi_h1,phi_l2:phi_h2)
  !
  !   call filcc_nd(phi,phi_l1,phi_l2,phi_h1,phi_h2, &
  !        domlo,domhi,delta,xlo,bc)
  !
  ! end subroutine ca_phigravfill

end module bc_fill_module
