! advection routines in support of method of lines integration

subroutine ca_mol_single_stage(time, &
                               lo, hi, domlo, domhi, &
                               stage_weight, &
                               uin, uin_lo, uin_hi, &
                               uout, uout_lo, uout_hi, &
                               q, q_lo, q_hi, &
                               qaux, qa_lo, qa_hi, &
                               srcU, srU_lo, srU_hi, &
                               update, updt_lo, updt_hi, &
                               update_flux, uf_lo, uf_hi, &
                               delta, dt, &
                               flux1, flux1_lo, flux1_hi, &
                               flux2, flux2_lo, flux2_hi, &
                               pradial, p_lo, p_hi, &
                               area1, area1_lo, area1_hi, &
                               area2, area2_lo, area2_hi, &
                               dloga, dloga_lo, dloga_hi, &
                               vol, vol_lo, vol_hi, &
                               courno, verbose) bind(C, name="ca_mol_single_stage")

  use meth_params_module, only : NQ, QVAR, NVAR, UMX, UMY,&
                                 NQAUX, QFS, QFX, QRHO,&
                                 first_order_hydro, difmag, URHO
  use advection_util_2d_module, only : normalize_species_fluxes
  use advection_util_module, only : compute_cfl
  use bl_constants_module, only : ZERO, HALF, ONE
  use prob_params_module, only : coord_type
  use riemann_module, only: cmpflx
  use reconstruct_module, only : compute_reconstruction_tvd
  use amrex_fort_module, only : rt => amrex_real
  use network, only : nspec, naux

  implicit none

  integer, intent(in) :: lo(2), hi(2), verbose
  integer, intent(in) :: domlo(2), domhi(2)
  real(rt), intent(in) :: stage_weight
  integer, intent(in) :: uin_lo(3), uin_hi(3)
  integer, intent(in) :: uout_lo(3), uout_hi(3)
  integer, intent(in) :: q_lo(3), q_hi(3)
  integer, intent(in) :: qa_lo(3), qa_hi(3)
  integer, intent(in) :: srU_lo(3), srU_hi(3)
  integer, intent(in) :: updt_lo(3), updt_hi(3)
  integer, intent(in) :: uf_lo(3), uf_hi(3)
  integer, intent(in) :: flux1_lo(3), flux1_hi(3)
  integer, intent(in) :: flux2_lo(3), flux2_hi(3)
  integer, intent(in) :: p_lo(3), p_hi(3)
  integer, intent(in) :: area1_lo(3), area1_hi(3)
  integer, intent(in) :: area2_lo(3), area2_hi(3)
  integer, intent(in) :: dloga_lo(3), dloga_hi(3)
  integer, intent(in) :: vol_lo(3), vol_hi(3)

  real(rt)        , intent(in) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),NVAR)
  real(rt)        , intent(inout) :: uout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2),NVAR)
  real(rt)        , intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),NQ)
  real(rt)        , intent(inout) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),NQAUX)
  real(rt)        , intent(in) :: srcU(srU_lo(1):srU_hi(1),srU_lo(2):srU_hi(2),NVAR)
  real(rt)        , intent(inout) :: update(updt_lo(1):updt_hi(1),updt_lo(2):updt_hi(2),NVAR)
  real(rt)        , intent(inout) :: update_flux(uf_lo(1):uf_hi(1),uf_lo(2):uf_hi(2),NVAR)
  real(rt)        , intent(inout) :: flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),NVAR)
  real(rt)        , intent(inout) :: flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2),NVAR)
  real(rt)        , intent(inout) :: pradial(p_lo(1):p_hi(1),p_lo(2):p_hi(2))
  real(rt)        , intent(in) :: area1(area1_lo(1):area1_hi(1),area1_lo(2):area1_hi(2))
  real(rt)        , intent(in) :: area2(area2_lo(1):area2_hi(1),area2_lo(2):area2_hi(2))
  real(rt)        , intent(in) :: dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2))
  real(rt)        , intent(in) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2))
  real(rt)        , intent(in) :: delta(2), dt, time
  real(rt)        , intent(inout) :: courno

  ! temporary interface values of the parabola
  real(rt), allocatable :: sxm(:,:), sxp(:,:), sym(:,:), syp(:,:)

  real(rt)        , allocatable:: qxm(:,:,:), qym(:,:,:)
  real(rt)        , allocatable:: qxp(:,:,:), qyp(:,:,:)

  integer ngf
  real(rt) :: dx, dy

  integer :: lo_3D(3), hi_3D(3)
  integer :: qs_lo(3), qs_hi(3)
  real(rt) :: dx_3D(3)

  real(rt) :: div1

  integer :: i, j, n

  ngf = 1

  lo_3D  = [lo(1), lo(2), 0]
  hi_3D  = [hi(1), hi(2), 0]

  dx_3D  = [delta(1), delta(2), ZERO]

  qs_lo = [lo(1)-1, lo(2)-1, 0]
  qs_hi = [hi(1)+2, hi(2)+2, 0]

  dx = delta(1)
  dy = delta(2)

  ! Check if we have violated the CFL criterion.
  call compute_cfl(q, q_lo, q_hi, &
                   qaux, qa_lo, qa_hi, &
                   lo_3D, hi_3D, dt, dx_3D, courno)

  ! sm and sp are the minus and plus parts of the parabola -- they are
  ! defined for a single zone, so for zone i, sm is the left value of
  ! the parabola and sp is the right value of the parabola
  allocate(sxm(q_lo(1):q_hi(1), q_lo(2):q_hi(2)))
  allocate(sxp(q_lo(1):q_hi(1), q_lo(2):q_hi(2)))
  allocate(sym(q_lo(1):q_hi(1), q_lo(2):q_hi(2)))
  allocate(syp(q_lo(1):q_hi(1), q_lo(2):q_hi(2)))

  ! qm and qp are the left and right states for an interface -- they
  ! are defined for a particular interface, with the convention that
  ! qm(i) and qp(i) correspond to the i-1/2 interface
  allocate ( qxm(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),NQ) )
  allocate ( qxp(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),NQ) )
  allocate ( qym(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),NQ) )
  allocate ( qyp(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),NQ) )

  ! Do reconstruction
  do n = 1, QVAR

     call compute_reconstruction_tvd(q(:,:,n), q_lo, q_hi, &
                                    sxm, sxp, sym, syp, sxm, sxp, q_lo, q_hi, & ! extra sxm, sxp are dummy
                                    lo_3D, hi_3D, [dx, dy, ZERO], 0, 0)

     ! Construct the interface states -- this is essentially just a
     ! reshuffling of interface states from zone-center indexing to
     ! edge-centered indexing
     do j = lo(2)-1, hi(2)+1
        do i = lo(1)-1, hi(1)+1

           ! x-edges

           ! left state at i-1/2 interface
           qxm(i,j,n) = sxp(i-1,j)

           ! right state at i-1/2 interface
           qxp(i,j,n) = sxm(i,j)

           ! y-edges

           ! left state at j-1/2 interface
           qym(i,j,n) = syp(i,j-1)

           ! right state at j-1/2 interface
           qyp(i,j,n) = sym(i,j)

        enddo
     enddo
  enddo

  deallocate(sxm, sxp, sym, syp)

  ! Get the fluxes from the Riemann solver
  call cmpflx(qxm, qxp, qs_lo, qs_hi, &
              flux1, flux1_lo, flux1_hi, &
              qaux, qa_lo, qa_hi, &
              1, lo(1), hi(1), lo(2), hi(2), domlo, domhi)

  call cmpflx(qym, qyp, qs_lo, qs_hi, &
              flux2, flux2_lo, flux2_hi, &
              qaux, qa_lo, qa_hi, &
              2, lo(1), hi(1), lo(2), hi(2), domlo, domhi)

  deallocate(qxm, qxp, qym, qyp)

  ! Normalize the species fluxes
  call normalize_species_fluxes(flux1, flux1_lo, flux1_hi, &
                                flux2, flux2_lo, flux2_hi, &
                                [lo(1), lo(2), 0], [hi(1), hi(2), 0])


  ! Make the update for this state

  ! For hydro, we will create an update source term that is
  ! essentially the flux divergence.  This can be added with dt to
  ! get the update

  do n = 1, NVAR
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           update(i,j,n) = update(i,j,n) + &
                ( flux1(i,j,n) * area1(i,j) - flux1(i+1,j,n) * area1(i+1,j) + &
                  flux2(i,j,n) * area2(i,j) - flux2(i,j+1,n) * area2(i,j+1) ) / vol(i,j)

           ! for storage
           update_flux(i,j,n) = update_flux(i,j,n) + stage_weight * update(i,j,n)

           ! include source terms
           update(i,j,n) = update(i,j,n) + srcU(i,j,n)
        enddo

     enddo
  enddo

  ! Scale the fluxes for the form we expect later in refluxing.
  do n = 1, NVAR
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)+1
           flux1(i,j,n) = dt * flux1(i,j,n) * area1(i,j)
        enddo
     enddo
  enddo

  do n = 1, NVAR
     do j = lo(2), hi(2)+1
        do i = lo(1), hi(1)
           flux2(i,j,n) = dt * flux2(i,j,n) * area2(i,j)
        enddo
     enddo
  enddo

end subroutine ca_mol_single_stage
