module bc_fill_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains

  subroutine ca_hypfill(adv,adv_lo,adv_hi, &
                        domlo,domhi,delta,xlo,time,bc,level) &
                        bind(C, name="ca_hypfill")

    use probdata_module
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UTEMP, NQ, NQAUX, QRHO, QTEMP, QU, QV, QW, QFS, QPRES, QREINT
    use interpolate_module
    use eos_module, only : eos
    use eos_type_module, only: eos_input_rt, eos_t, eos_input_rp
    use network, only: nspec
    !use model_parser_module
    use prob_params_module, only: dim
    use c_interface_modules, only: ca_ctoprim
    use riemann_util_module, only: comp_cons_state

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer, intent(in) :: adv_lo(3),adv_hi(3), level
    integer, intent(in) :: bc(dim,2,*)
    integer, intent(in) :: domlo(3), domhi(3)
    real(rt), intent(in)         :: delta(3), xlo(3), time
    real(rt), intent(inout)         :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NVAR)

    integer :: i,j,k,q,n,iter
    real(rt)         :: x, h
    real(rt)         :: pres_above,p_want,pres_zone, A
    real(rt)         :: drho,dpdr,temp_zone,eint,X_zone(nspec),dens_zone
    real(rt) :: qprim(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NQ)
    real(rt) :: qaux(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NQAUX)

    integer, parameter :: MAX_ITER = 100
    real(rt)        , parameter :: TOL = 1.e-8_rt
    logical :: converged_hse

    type (eos_t) :: eos_state

    do n = 1,NVAR
       call filcc_nd(adv(:,:,:,n),adv_lo,adv_hi,domlo,domhi,delta,xlo,bc(:,:,n))
    enddo

    return

    do n = 1,NVAR
        do i = 2,3
            do j = 1,2
                call filcc_nd(adv(:,:,:,n),adv_lo,adv_hi, &
                           domlo,domhi,delta,xlo,bc(i,j,n))
           enddo
       enddo
       if (level <= swe_to_comp_level) then
           do j = 1,2
               call filcc_nd(adv(:,:,:,n),adv_lo,adv_hi, &
                          domlo,domhi,delta,xlo,bc(1,j,n))
          enddo
       endif
       ! call filcc_nd(adv(adv_lo(1),adv_lo(2),adv_lo(3),n),adv_lo,adv_hi, &
       !            domlo,domhi,delta,xlo,bc(2,1,n))
       !
       ! call filcc_nd(adv(adv_lo(1),adv_lo(2),adv_lo(3),n),adv_lo,adv_hi, &
       !           domlo,domhi,delta,xlo,bc(3,1,n))
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

       if (level > swe_to_comp_level) then

           call ca_ctoprim(adv_lo, adv_hi, adv, adv_lo, adv_hi, qprim, adv_lo, adv_hi, qaux, adv_lo, adv_hi, 1, level, xlo, delta)

           !        XLO
           if ( bc(1,1,n).eq.EXT_DIR .and. adv_lo(1).lt.domlo(1)) then

              ! this do loop counts backwards since we want to work downward
              do k = adv_lo(3), adv_hi(3)
                  do j = adv_lo(2), adv_hi(2)
                      h = height_from_p(qprim(domlo(1),j,k,QPRES), xlo(1) + delta(1)*(dble(domlo(1)-adv_lo(1)) + 0.5e0_rt))
                     do i = domlo(1)-1, adv_lo(1),-1
                         x = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5e0_rt)

                        ! set all the variables even though we're testing on URHO
                        if (n .eq. URHO) then

                          ! HSE integration to get density, pressure

                          ! initial guesses
                          dens_zone = rho_from_height(h, x)
                          pres_zone = p_from_rho(dens_zone)

                          ! get pressure in zone above
                          eos_state%rho = dens_zone
                          eos_state%p = pres_zone
                          eos_state%xn(:) = qprim(i+1,j,k,QFS:QFS-1+nspec)

                          call eos(eos_input_rp, eos_state)

                          eint = eos_state%e
                          temp_zone = eos_state%T

                          qprim(i,j,k,:) = qprim(domlo(1),j,k,:)

                          qprim(i,j,k,QRHO) = dens_zone
                          qprim(i,j,k,QPRES) = pres_zone
                          qprim(i,j,k,QTEMP) = temp_zone
                          qprim(i,j,k,QREINT) = eint * dens_zone

                          call comp_cons_state(qprim(i,j,k,:), adv(i,j,k,:))

                        end if

                     end do
                  end do
               end do
           end if

           !        XHI
           if ( bc(1,2,n).eq.EXT_DIR .and. adv_hi(1).gt.domhi(1)) then

              do k=adv_lo(3),adv_hi(3)
                  do j=adv_lo(2),adv_hi(2)
                      h = height_from_p(qprim(domhi(1),j,k,QPRES), xlo(1) + delta(1)*(dble(domhi(1)-adv_lo(1)) + 0.5e0_rt))

                     do i=domhi(1)+1,adv_hi(1)
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

                           dens_zone = interpolate_noparser(adv(:,j,k,URHO),adv_lo(1),adv_hi(1),i)

                           temp_zone = interpolate_noparser(adv(:,j,k,UTEMP),adv_lo(1),adv_hi(1),i)

                           do q = 1, nspec
                              X_zone(q) = interpolate_noparser(adv(:,j,k,UFS-1+q),adv_lo(1),adv_hi(1),i)
                           enddo

                           dens_zone = rho_from_height(h, x)
                           pres_zone = p_from_rho(dens_zone)

                           ! get pressure in zone above
                           eos_state%rho = dens_zone
                           eos_state%p = pres_zone
                           eos_state%xn(:) = qprim(i+1,j,k,QFS:QFS-1+nspec)

                           call eos(eos_input_rp, eos_state)

                           eint = eos_state%e
                           temp_zone = eos_state%T

                           qprim(i,j,k,:) = qprim(domhi(1),j,k,:)

                           qprim(i,j,k,QRHO) = dens_zone
                           qprim(i,j,k,QPRES) = pres_zone
                           qprim(i,j,k,QTEMP) = temp_zone
                           qprim(i,j,k,QREINT) = eint * dens_zone

                           call comp_cons_state(qprim(i,j,k,:), adv(i,j,k,:))


                           ! ! extrap normal momentum
                           ! adv(i,j,k,UMY) = max(0.e0_rt,adv(domhi(1),j,k,UMY))
                           !
                           ! ! zero transverse momentum
                           ! adv(i,j,k,UMX) = 0.e0_rt
                           ! adv(i,j,k,UMZ) = 0.e0_rt
                           !
                           ! eos_state%rho = dens_zone
                           ! eos_state%T = temp_zone
                           ! eos_state%xn(:) = X_zone
                           !
                           ! call eos(eos_input_rt, eos_state)
                           !
                           ! pres_zone = eos_state%p
                           ! eint = eos_state%e
                           !
                           ! adv(i,j,k,URHO) = dens_zone
                           ! adv(i,j,k,UEINT) = dens_zone*eint
                           ! adv(i,j,k,UEDEN) = dens_zone*eint + &
                           !      0.5e0_rt*(adv(i,j,k,UMX)**2+adv(i,j,k,UMY)**2+adv(i,j,k,UMZ)**2)/dens_zone
                           ! adv(i,j,k,UTEMP) = temp_zone
                           ! adv(i,j,k,UFS:UFS-1+nspec) = dens_zone*X_zone(:)

                        end if

                     end do
                  end do
               end do
           end if
       endif

    end do

  end subroutine ca_hypfill



  subroutine ca_denfill(adv,adv_lo,adv_hi, &
                        domlo,domhi,delta,xlo,time,bc,level) &
                        bind(C, name="ca_denfill")

    use probdata_module
    use interpolate_module
    use model_parser_module
    use bl_error_module
    use prob_params_module, only: dim

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer, intent(in) :: adv_lo(3),adv_hi(3),level
    integer, intent(in) :: bc(dim,2,*)
    integer, intent(in) :: domlo(3), domhi(3)
    real(rt), intent(in)         :: delta(3), xlo(3), time
    real(rt), intent(inout)         :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))

    integer :: i,j,k
    ! real(rt)         :: x

    ! Note: this function should not be needed, technically, but is
    ! provided to filpatch because there are many times in the algorithm
    ! when just the density is needed.  We try to rig up the filling so
    ! that the same function is called here and in hypfill where all the
    ! states are filled.

    call filcc_nd(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,bc)

    return 

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

    if (level > swe_to_comp_level) then
        !     XLO
        if ( bc(1,1,1).eq.EXT_DIR .and. adv_lo(1).lt.domlo(1)) then
           do k=adv_lo(3),adv_hi(3)
               do j=adv_lo(2),adv_hi(2)
                  do i=adv_lo(1),domlo(1)-1
                     ! x = xlo(1) + delta(1)*(dble(i-adv_lo(1)) + 0.5e0_rt)
                     adv(i,j,k) = interpolate_noparser(adv(:,j,k),adv_lo(1),adv_hi(1),i)
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
                     ! x = xlo(1) + delta(1)*(dble(i-adv_lo(1))+ 0.5e0_rt)
                     adv(i,j,k) = interpolate_noparser(adv(:,j,k),adv_lo(1),adv_hi(1),i)
                     !interpolate(x,npts_model,model_r,model_state(:,idens_model))
                  end do
               end do
            end do
        end if

    endif

  end subroutine ca_denfill


end module bc_fill_module
