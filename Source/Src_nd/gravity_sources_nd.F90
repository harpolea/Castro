module gravity_sources_module

  implicit none

  public

contains

  subroutine ca_gsrc(lo,hi,domlo,domhi, &
                     uold,uold_lo,uold_hi, &
                     source,src_lo,src_hi, &
                     dx,dt,time) bind(C, name="ca_gsrc")

    use amrex_fort_module, only: rt => amrex_real
    use bl_constants_module, only: ZERO, HALF, ONE
    use meth_params_module, only: NVAR, URHO, UMX, UMZ, UEDEN, grav_source_type
    use math_module, only: cross_product
    use castro_util_module, only: position
    use prob_params_module, only: center

    implicit none

    integer, intent(in)     :: lo(3), hi(3)
    integer, intent(in)     :: domlo(3), domhi(3)
    integer, intent(in)     :: uold_lo(3), uold_hi(3)
    integer, intent(in)     :: src_lo(3), src_hi(3)

    real(rt), intent(in)    :: uold(uold_lo(1):uold_hi(1),uold_lo(2):uold_hi(2),uold_lo(3):uold_hi(3),NVAR)
    real(rt), intent(inout) :: source(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),NVAR)
    real(rt), intent(in)    :: dx(3), dt, time

    real(rt) :: rho, rhoInv
    real(rt) :: Sr(3), SrE
    real(rt) :: old_ke, new_ke
    real(rt) :: loc(3)

    integer  :: i, j, k

    ! Temporary array for holding the update to the state.

    real(rt) :: src(NVAR)

    ! Temporary array for seeing what the new state would be if the update were applied here.

    real(rt) :: snew(NVAR)

    ! Initialize the update and temporary state to zero. We only need to do this once outside
    ! the loop, since the array access pattern is consistent across loop iterations.

    Sr(:) = ZERO
    src(:) = ZERO
    snew(:) = ZERO

    ! For constant gravity, we can just initialize the gravitational acceleration here.

    ! Gravitational source options for how to add the work to (rho E):
    ! grav_source_type =
    ! 1: Original version ("does work")
    ! 2: Modification of type 1 that updates the momentum before constructing the energy corrector
    ! 3: Puts all gravitational work into KE, not (rho e)
    ! 4: Conservative energy formulation

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rho    = uold(i,j,k,URHO)
             rhoInv = ONE / rho

             loc = position(i,j,k) - center

             src = ZERO
             snew = uold(i,j,k,:)

             old_ke = HALF * sum(snew(UMX:UMZ)**2) * rhoInv

             Sr(dim) = rho * const_grav

             src(UMX:UMZ) = Sr

             snew(UMX:UMZ) = snew(UMX:UMZ) + dt * src(UMX:UMZ)

             if (grav_source_type == 1 .or. grav_source_type == 2) then

                ! Src = rho u dot g, evaluated with all quantities at t^n

                SrE = dot_product(uold(i,j,k,UMX:UMZ) * rhoInv, Sr)

             else if (grav_source_type .eq. 3) then

                new_ke = HALF * sum(snew(UMX:UMZ)**2) * rhoInv
                SrE = new_ke - old_ke

             else if (grav_source_type .eq. 4) then

                ! The conservative energy formulation does not strictly require
                ! any energy source-term here, because it depends only on the
                ! fluid motions from the hydrodynamical fluxes which we will only
                ! have when we get to the 'corrector' step. Nevertheless we add a
                ! predictor energy source term in the way that the other methods
                ! do, for consistency. We will fully subtract this predictor value
                ! during the corrector step, so that the final result is correct.
                ! Here we use the same approach as grav_source_type == 2.

                SrE = dot_product(uold(i,j,k,UMX:UMZ) * rhoInv, Sr)

             else

                call bl_error("Error:: gravity_sources_nd.F90 :: invalid grav_source_type")

             end if

             src(UEDEN) = SrE

             snew(UEDEN) = snew(UEDEN) + dt * SrE

             ! Add to the outgoing source array.

             source(i,j,k,:) = src

          enddo
       enddo
    enddo

  end subroutine ca_gsrc

  ! :::
  ! ::: ------------------------------------------------------------------
  ! :::

  subroutine ca_corrgsrc(lo,hi,domlo,domhi, &
                         uold,uo_lo,uo_hi, &
                         unew,un_lo,un_hi, &
                         vol,vol_lo,vol_hi, &
                         flux1,f1_lo,f1_hi, &
                         flux2,f2_lo,f2_hi, &
                         flux3,f3_lo,f3_hi, &
                         source,sr_lo,sr_hi, &
                         dx,dt,time) bind(C, name="ca_corrgsrc")

    use amrex_fort_module, only: rt => amrex_real
    use bl_constants_module, only: ZERO, HALF, ONE, TWO
    use mempool_module, only: bl_allocate, bl_deallocate
    use meth_params_module, only: NVAR, URHO, UMX, UMZ, UEDEN, &
                                  grav_source_type, gravity_type, get_g_from_phi
    use prob_params_module, only: dg, center
    use fundamental_constants_module, only: Gconst
    use castro_util_module, only: position

    implicit none

    integer, intent(in)     :: lo(3), hi(3)
    integer, intent(in)     :: domlo(3), domhi(3)
    integer, intent(in)     :: uo_lo(3), uo_hi(3)
    integer, intent(in)     :: un_lo(3), un_hi(3)
    integer, intent(in)     :: vol_lo(3), vol_hi(3)
    integer, intent(in)     :: f1_lo(3), f1_hi(3)
    integer, intent(in)     :: f2_lo(3), f2_hi(3)
    integer, intent(in)     :: f3_lo(3), f3_hi(3)

    integer, intent(in)     :: sr_lo(3), sr_hi(3)

    ! Old and new time state data

    real(rt), intent(in)    :: uold(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3),NVAR)
    real(rt), intent(in)    :: unew(un_lo(1):un_hi(1),un_lo(2):un_hi(2),un_lo(3):un_hi(3),NVAR)

    ! Cell volume

    real(rt), intent(in)    :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))

    ! Hydrodynamics fluxes

    real(rt), intent(in)    :: flux1(f1_lo(1):f1_hi(1),f1_lo(2):f1_hi(2),f1_lo(3):f1_hi(3),NVAR)
    real(rt), intent(in)    :: flux2(f2_lo(1):f2_hi(1),f2_lo(2):f2_hi(2),f2_lo(3):f2_hi(3),NVAR)
    real(rt), intent(in)    :: flux3(f3_lo(1):f3_hi(1),f3_lo(2):f3_hi(2),f3_lo(3):f3_hi(3),NVAR)

    ! The source term to send back

    real(rt), intent(inout) :: source(sr_lo(1):sr_hi(1),sr_lo(2):sr_hi(2),sr_lo(3):sr_hi(3),NVAR)

    real(rt), intent(in)    :: dx(3), dt, time

    integer  :: i, j, k

    real(rt) :: Sr_old(3), Sr_new(3), Srcorr(3)
    real(rt) :: vold(3), vnew(3)
    real(rt) :: SrE_old, SrE_new, SrEcorr
    real(rt) :: rhoo, rhooinv, rhon, rhoninv

    real(rt) :: old_ke, new_ke
    real(rt) :: loc(3)

    real(rt) :: hdtInv

    real(rt) :: src(NVAR)

    ! Temporary array for seeing what the new state would be if the update were applied here.

    real(rt) :: snew(NVAR)

    real(rt), pointer :: phi(:,:,:)
    real(rt), pointer :: grav(:,:,:,:)
    real(rt), pointer :: gravx(:,:,:)
    real(rt), pointer :: gravy(:,:,:)
    real(rt), pointer :: gravz(:,:,:)

    Sr_old(:) = ZERO
    Sr_new(:) = ZERO
    Srcorr(:) = ZERO
    src(:) = ZERO
    snew(:) = ZERO

    hdtInv = HALF / dt

    ! Gravitational source options for how to add the work to (rho E):
    ! grav_source_type =
    ! 1: Original version ("does work")
    ! 2: Modification of type 1 that updates the U before constructing SrEcorr
    ! 3: Puts all gravitational work into KE, not (rho e)
    ! 4: Conservative gravity approach (discussed in first white dwarf merger paper).

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             loc = position(i,j,k) - center

             rhoo    = uold(i,j,k,URHO)
             rhooinv = ONE / uold(i,j,k,URHO)

             rhon    = unew(i,j,k,URHO)
             rhoninv = ONE / unew(i,j,k,URHO)

             snew = unew(i,j,k,:)

             old_ke = HALF * sum(snew(UMX:UMZ)**2) * rhoninv

             ! Define old source terms

             vold = uold(i,j,k,UMX:UMZ) * rhooinv

             Sr_old(dim) = rhoo * const_grav

             SrE_old = dot_product(vold, Sr_old)

             ! Define new source terms

             vnew = snew(UMX:UMZ) * rhoninv
             Sr_new(dim) = rhon * const_grav
             SrE_new = dot_product(vnew, Sr_new)

             ! Define corrections to source terms

             Srcorr = HALF * (Sr_new - Sr_old)

             ! Correct momenta

             src(UMX:UMZ) = Srcorr

             snew(UMX:UMZ) = snew(UMX:UMZ) + dt * src(UMX:UMZ)

             ! Correct energy

             if (grav_source_type .eq. 1) then

                ! If grav_source_type == 1, then we calculated SrEcorr before updating the velocities.

                SrEcorr = HALF * (SrE_new - SrE_old)

             else if (grav_source_type .eq. 2) then

                ! For this source type, we first update the momenta
                ! before we calculate the energy source term.

                vnew = snew(UMX:UMZ) * rhoninv
                SrE_new = dot_product(vnew, Sr_new)

                SrEcorr = HALF * (SrE_new - SrE_old)

             else if (grav_source_type .eq. 3) then

                ! Instead of calculating the energy source term explicitly,
                ! we simply update the kinetic energy.

                new_ke = HALF * sum(snew(UMX:UMZ)**2) * rhoninv
                SrEcorr = new_ke - old_ke

             else if (grav_source_type .eq. 4) then

                ! First, subtract the predictor step we applied earlier.

                SrEcorr = - SrE_old

                ! For an explanation of this approach, see wdmerger paper I.
                ! The main idea is that we are evaluating the change of the
                ! potential energy at zone edges and applying that in an equal
                ! and opposite sense to the gas energy. The physics is described
                ! in Section 2.4; the particular form of the equation we are using
                ! is found in Appendix B, as it provides the best numerical conservation
                ! properties when using AMR.


                ! For constant gravity, the only contribution is from the dimension that the gravity points in.

                if (dim .eq. 1) then
                   SrEcorr = SrEcorr + (HALF / dt) * ( flux1(i        ,j,k,URHO) * const_grav * dx(1) + &
                                                       flux1(i+1*dg(1),j,k,URHO) * const_grav * dx(1) ) / vol(i,j,k)
                else if (dim .eq. 2) then
                   SrEcorr = SrEcorr + (HALF / dt) * ( flux2(i,j        ,k,URHO) * const_grav * dx(2) + &
                                                       flux2(i,j+1*dg(2),k,URHO) * const_grav * dx(2) ) / vol(i,j,k)
                else if (dim .eq. 3) then
                   SrEcorr = SrEcorr + (HALF / dt) * ( flux3(i,j,k        ,URHO) * const_grav * dx(3) + &
                                                       flux3(i,j,k+1*dg(3),URHO) * const_grav * dx(3) ) / vol(i,j,k)
                end if

             else
                call bl_error("Error:: gravity_sources_nd.F90 :: invalid grav_source_type")
             end if

             src(UEDEN) = SrEcorr

             snew(UEDEN) = snew(UEDEN) + dt * SrEcorr

             ! Add to the outgoing source array.

             source(i,j,k,:) = src

          enddo
       enddo
    enddo

  end subroutine ca_corrgsrc

end module gravity_sources_module
