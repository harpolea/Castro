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
                               flux, flux_lo, flux_hi, &
                               pradial, p_lo, p_hi, &
                               area, area_lo, area_hi, &
                               dloga, dloga_lo, dloga_hi, &
                               vol, vol_lo, vol_hi, &
                               courno, verbose) bind(C, name="ca_mol_single_stage")

  use meth_params_module, only : NQ, QVAR, QU, QPRES, &
                                 UMX, UMY, UMZ, UTEMP, USHK, UEINT, &
                                 NQAUX, NVAR, NHYP, use_flattening, &
                                 QTEMP, QFS, QFX, QREINT, QRHO, &
                                 NGDNV, GDU, GDPRES, first_order_hydro, difmag, &
                                 hybrid_riemann, ppm_temp_fix
  use advection_util_module, only : compute_cfl, shock, normalize_species_fluxes
  use bl_constants_module, only : ZERO, HALF, ONE
  use flatten_module, only : uflatten
  use prob_params_module, only : coord_type
  use riemann_module, only: cmpflx
  use ppm_module, only : ppm_reconstruct
  use amrex_fort_module, only : rt => amrex_real
  use eos_type_module, only : eos_t, eos_input_rt
  use eos_module, only : eos
  use network, only : nspec, naux

  implicit none

  integer, intent(in) :: lo(1), hi(1), verbose
  integer, intent(in) :: domlo(1), domhi(1)
  real(rt), intent(in) :: stage_weight
  integer, intent(in) :: uin_lo(3), uin_hi(3)
  integer, intent(in) :: uout_lo(3), uout_hi(3)
  integer, intent(in) :: q_lo(3), q_hi(3)
  integer, intent(in) :: qa_lo(3), qa_hi(3)
  integer, intent(in) :: srU_lo(3), srU_hi(3)
  integer, intent(in) :: updt_lo(3), updt_hi(3)
  integer, intent(in) :: uf_lo(3), uf_hi(3)
  integer, intent(in) :: flux_lo(3), flux_hi(3)
  integer, intent(in) :: p_lo(3), p_hi(3)
  integer, intent(in) :: area_lo(3), area_hi(3)
  integer, intent(in) :: dloga_lo(3), dloga_hi(3)
  integer, intent(in) :: vol_lo(3), vol_hi(3)

  real(rt)        , intent(in) ::      uin(  uin_lo(1):  uin_hi(1),NVAR)
  real(rt)        , intent(inout) ::  uout( uout_lo(1): uout_hi(1),NVAR)
  real(rt)        , intent(inout) ::     q(    q_lo(1):    q_hi(1),NQ)
  real(rt)        , intent(in) ::     qaux(   qa_lo(1):   qa_hi(1),NQAUX)
  real(rt)        , intent(in) ::     srcU(  srU_lo(1):  srU_hi(1),NVAR)
  real(rt)        , intent(inout) :: update(updt_lo(1): updt_hi(1),NVAR)
  real(rt)        , intent(inout) :: update_flux(uf_lo(1): uf_hi(1),NVAR)
  real(rt)        , intent(inout) ::  flux( flux_lo(1): flux_hi(1),NVAR)
  real(rt)        , intent(inout) :: pradial(  p_lo(1):   p_hi(1))
  real(rt)        , intent(in) :: area( area_lo(1): area_hi(1)     )
  real(rt)        , intent(in) :: dloga(dloga_lo(1):dloga_hi(1)     )
  real(rt)        , intent(in) ::   vol(  vol_lo(1): vol_hi(1)      )
  real(rt)        , intent(in) :: delta(1), dt, time
  real(rt)        , intent(inout) :: courno

  ! temporary interface values of the parabola
  real(rt), allocatable :: sxm(:), sxp(:)

  real(rt), allocatable :: qxm(:,:), qxp(:,:)

  real(rt)         :: dx

  integer i, n, ngf

  integer :: lo_3D(3), hi_3D(3)
  real(rt)         :: dx_3D(3)
  integer :: qp_lo(3), qp_hi(3)
  type(eos_t) :: eos_state

  ngf = 1

  lo_3D   = [lo(1), 0, 0]
  hi_3D   = [hi(1), 0, 0]

  dx_3D   = [delta(1), ZERO, ZERO]

  qp_lo = [lo(1)-1, 0, 0]
  qp_hi = [hi(1)+2, 0, 0]

  dx = delta(1)


  ! Check if we have violated the CFL criterion.
  call compute_cfl(q, q_lo, q_hi, &
                   qaux, qa_lo, qa_hi, &
                   lo_3D, hi_3D, dt, dx_3D, courno)


  ! sm and sp are the minus and plus parts of the parabola -- they are
  ! defined for a single zone, so for zone i, sm is the left value of
  ! the parabola and sp is the right value of the parabola
  allocate(sxm(q_lo(1):q_hi(1)))
  allocate(sxp(q_lo(1):q_hi(1)))

  ! qm and qp are the left and right states for an interface -- they
  ! are defined for a particular interface, with the convention that
  ! qm(i) and qp(i) correspond to the i-1/2 interface
  allocate ( qxm(qp_lo(1):qp_hi(1),NQ) )
  allocate ( qxp(qp_lo(1):qp_hi(1),NQ) )

  ! Do PPM reconstruction
  do n = 1, QVAR
     call compute_reconstruction_tvd(q(:,n), q_lo, q_hi, &
                          sxm, sxp, sxm, sxp, sxm, sxp, q_lo, q_hi, &  ! extras are dummy
                          lo(1), 0, hi(1), 0, [dx, ZERO, ZERO], 0, 0)

     ! Construct the interface states -- this is essentially just a
     ! reshuffling of interface states from zone-center indexing to
     ! edge-centered indexing
     do i = lo(1)-1, hi(1)+1
        qxm(i,n) = sxp(i-1)
        qxp(i,n) = sxm(i)
     enddo
  enddo

  deallocate(sxm, sxp)

  ! Get the fluxes from the Riemann solver
  call cmpflx(qxm, qxp, qp_lo, qp_hi, &
              flux, flux_lo, flux_hi, &
              qaux, qa_lo, qa_hi, lo(1), hi(1))

  deallocate(qxm, qxp)

  ! Normalize the species fluxes
  call normalize_species_fluxes(flux, flux_lo(1), flux_hi(1), lo_3D, hi_3D)


  ! Make the update for this state

  ! For hydro, we will create an update source term that is
  ! essentially the flux divergence.  This can be added with dt to
  ! get the update

  do n = 1, NVAR
     do i = lo(1), hi(1)
        update(i,n) = update(i,n) + ( flux(i,n) * area(i) - &
                                      flux(i+1,n) * area(i+1) ) / vol(i)

            ! for storage
        update_flux(i,n) = update_flux(i,n) + stage_weight * update(i,n)

        ! include source terms
        update(i,n) = update(i,n) + srcU(i,n)

     enddo
  enddo

  ! Scale the fluxes for the form we expect later in refluxing.

  do n = 1, NVAR
     do i = lo(1), hi(1)+1

        flux(i,n) = dt * area(i) * flux(i,n)

        ! Correct the momentum flux with the grad p part.
        if (coord_type .eq. 0 .and. n == UMX) then
           flux(i,n) = flux(i,n) + dt * area(i) * q1(i,GDPRES)
        endif

     enddo
  enddo

end subroutine ca_mol_single_stage
