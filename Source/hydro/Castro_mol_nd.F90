! advection routines in support of method of lines integration

subroutine ca_mol_single_stage(time, level, &
                               lo, hi, domlo, domhi, &
                               stage_weight, &
                               uin, uin_lo, uin_hi, &
                               uout, uout_lo, uout_hi, &
                               q, q_lo, q_hi, &
                               qaux, qa_lo, qa_hi, &
                               srcU, srU_lo, srU_hi, &
                               update, updt_lo, updt_hi, &
                               update_flux, uf_lo, uf_hi, &
                               dx, dt, &
                               flux1, flux1_lo, flux1_hi, &
#if BL_SPACEDIM >=2
                               flux2, flux2_lo, flux2_hi, &
#endif
#if BL_SPACEDIM == 3
                               flux3, flux3_lo, flux3_hi, &
#else
                               pradial, p_lo, p_hi, &
#endif
                               area1, area1_lo, area1_hi, &
#if BL_SPACEDIM >=2
                               area2, area2_lo, area2_hi, &
#endif
#if BL_SPACEDIM == 3
                               area3, area3_lo, area3_hi, &
#else
                               dloga, dloga_lo, dloga_hi, &
#endif
                               vol, vol_lo, vol_hi, &
                               courno, verbose, xlo) bind(C, name="ca_mol_single_stage")

  use meth_params_module, only : NQ, QVAR, NVAR, QTEMP, QFA,&
                                 NQAUX, QFS, QFX, QREINT, QRHO, QPRES
  use advection_util_module
  use reconstruct_module, only : compute_reconstruction_tvd
  use bl_constants_module, only : ZERO, HALF, ONE, FOURTH
  use riemann_module, only: cmpflx
  use amrex_fort_module, only : rt => amrex_real
  use eos_type_module, only : eos_t, eos_input_rt
  use network, only : nspec, naux
  use probdata_module, only: swe_to_comp_level
  use prob_params_module, only : dg

  implicit none

  integer, intent(in) :: lo(3), hi(3), verbose, level
  integer, intent(in) ::  domlo(3), domhi(3)
  real(rt), intent(in) :: stage_weight, xlo(3)
  integer, intent(in) :: uin_lo(3), uin_hi(3)
  integer, intent(in) :: uout_lo(3), uout_hi(3)
  integer, intent(in) :: q_lo(3), q_hi(3)
  integer, intent(in) :: qa_lo(3), qa_hi(3)
  integer, intent(in) :: srU_lo(3), srU_hi(3)
  integer, intent(in) :: updt_lo(3), updt_hi(3)
  integer, intent(in) :: uf_lo(3), uf_hi(3)
  integer, intent(in) :: flux1_lo(3), flux1_hi(3)
  integer, intent(in) :: area1_lo(3), area1_hi(3)
#if BL_SPACEDIM >= 2
  integer, intent(in) :: flux2_lo(3), flux2_hi(3)
  integer, intent(in) :: area2_lo(3), area2_hi(3)
#endif
#if BL_SPACEDIM == 3
  integer, intent(in) :: flux3_lo(3), flux3_hi(3)
  integer, intent(in) :: area3_lo(3), area3_hi(3)
#endif
#if BL_SPACEDIM <= 2
  integer, intent(in) :: p_lo(3), p_hi(3)
  integer, intent(in) :: dloga_lo(3), dloga_hi(3)
#endif
  integer, intent(in) :: vol_lo(3), vol_hi(3)

  real(rt)        , intent(in) :: uin(uin_lo(1):uin_hi(1), uin_lo(2):uin_hi(2), uin_lo(3):uin_hi(3), NVAR)
  real(rt)        , intent(inout) :: uout(uout_lo(1):uout_hi(1), uout_lo(2):uout_hi(2), uout_lo(3):uout_hi(3), NVAR)
  real(rt)        , intent(inout) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), NQ)
  real(rt)        , intent(inout) :: qaux(qa_lo(1):qa_hi(1), qa_lo(2):qa_hi(2), qa_lo(3):qa_hi(3), NQAUX)
  real(rt)        , intent(in) :: srcU(srU_lo(1):srU_hi(1), srU_lo(2):srU_hi(2), srU_lo(3):srU_hi(3), NVAR)
  real(rt)        , intent(inout) :: update(updt_lo(1):updt_hi(1), updt_lo(2):updt_hi(2), updt_lo(3):updt_hi(3), NVAR)
  real(rt)        , intent(inout) :: update_flux(uf_lo(1):uf_hi(1), uf_lo(2):uf_hi(2), uf_lo(3):uf_hi(3), NVAR)

  real(rt)        , intent(inout) :: flux1(flux1_lo(1):flux1_hi(1), flux1_lo(2):flux1_hi(2), flux1_lo(3):flux1_hi(3), NVAR)
  real(rt)        , intent(in) :: area1(area1_lo(1):area1_hi(1), area1_lo(2):area1_hi(2), area1_lo(3):area1_hi(3))
#if BL_SPACEDIM >= 2
  real(rt)        , intent(inout) :: flux2(flux2_lo(1):flux2_hi(1), flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3), NVAR)
  real(rt)        , intent(in) :: area2(area2_lo(1):area2_hi(1), area2_lo(2):area2_hi(2), area2_lo(3):area2_hi(3))
#endif
#if BL_SPACEDIM == 3
  real(rt)        , intent(inout) :: flux3(flux3_lo(1):flux3_hi(1), flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3), NVAR)

  real(rt)        , intent(in) :: area3(area3_lo(1):area3_hi(1), area3_lo(2):area3_hi(2), area3_lo(3):area3_hi(3))
#endif
#if BL_SPACEDIM <= 2
  real(rt)        , intent(inout) :: pradial(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
  real(rt)        , intent(in) :: dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
#endif
  real(rt)        , intent(in) :: vol(vol_lo(1):vol_hi(1), vol_lo(2):vol_hi(2), vol_lo(3):vol_hi(3))
  real(rt)        , intent(in) :: dx(3), dt, time
  real(rt)        , intent(inout) :: courno

  ! temporary interface values of the parabola
  real(rt)        , allocatable :: sxm(:,:,:), qxm(:,:,:,:)
  real(rt)        , allocatable :: sxp(:,:,:), qxp(:,:,:,:)
#if BL_SPACEDIM >= 2
    real(rt)        , allocatable :: sym(:,:,:), qym(:,:,:,:)
    real(rt)        , allocatable :: syp(:,:,:), qyp(:,:,:,:)
#endif
#if BL_SPACEDIM == 3
    real(rt)        , allocatable :: szm(:,:,:), qzm(:,:,:,:)
    real(rt)        , allocatable :: szp(:,:,:), qzp(:,:,:,:)
#endif
  real(rt) :: frac_change

  integer :: i, j, k, n
  integer :: st_lo(3), st_hi(3)
  integer :: It_lo(3), It_hi(3)

  type (eos_t) :: eos_state

  It_lo = [lo(1)-1, lo(2) - dg(2), lo(3) - dg(3)]
  It_hi = [hi(1)+1, hi(2) + dg(2), hi(3) + dg(3)]

  st_lo = [lo(1)-2, lo(2) - 2*dg(2), lo(3) - 2*dg(3)]
  st_hi = [hi(1)+2, hi(2) + 2*dg(2), hi(3) + 2*dg(3)]

#if BL_SPACEDIM <= 2
  It_lo(3) = 0
  It_hi(3) = 0
  st_lo(3) = 0
  st_hi(3) = 0
#endif

    ! write(*,*) "uin_lo, uin_hi", uin_lo, uin_hi
    ! write(*,*) "q_lo, q_hi", q_lo, q_hi
    ! write(*,*) "st_lo, st_hi", st_lo, st_hi
    ! write(*,*) "dg = ", dg

  if (level <= swe_to_comp_level) then
      call swectoprim(q_lo, q_hi, &
                   uin, uin_lo, uin_hi, &
                   q,     q_lo,   q_hi, &
                   qaux, qa_lo,  qa_hi)
   else
       call compctoprim(q_lo, q_hi, &
                    uin, uin_lo, uin_hi, &
                    q,     q_lo,   q_hi, &
                    qaux, qa_lo,  qa_hi, xlo, dx)
   endif

   ! write(*,*) "QREINT", QREINT, "QPRES", QPRES, "QFA", QFA

  ! nan check
  do n = 1, NQ
#if (BL_SPACEDIM <= 2)
      k = 0!lo(3)
#else
      do k=q_lo(3), q_hi(3)
#endif
         do j = q_lo(2), q_hi(2)
            do i = q_lo(1), q_hi(1)
                if (q(i,j,k,n) /= q(i,j,k,n)) then
                    ! write(*,*) i,j,k,n,q(i,j,k,n)
                    if (n==1) then
                        q(i,j,k,n) = 1.0e0_rt
                    else
                        q(i,j,k,n) = 0.0e0_rt
                    endif
                endif
            enddo
        enddo
#if (BL_SPACEDIM == 3)
    enddo
#endif
  enddo

  ! NOTE: don't do this
  ! call enforce_minimum_density(q, q_lo, q_hi, &
  !                              q, q_lo, q_hi, &
  !                              vol, vol_lo, vol_hi, &
  !                              lo, hi, frac_change, verbose, level)

  ! Check if we have violated the CFL criterion.
  ! call compute_cfl(q, q_lo, q_hi, &
  !                  qaux, qa_lo, qa_hi, &
  !                  lo, hi, dt, dx, courno)

  ! We come into this routine with a 3-d box of data, but we operate
  ! on it locally by considering 2 planes that encompass all of the
  ! x, y indices of the original box, but each plane corresponds to
  ! a single z index.
  !
  ! In the notation below, k3d will always been the index into the
  ! original 3-d box.  kc will be the z-index in the local "planar"
  ! data and km will be the previously used index in the local
  ! planar data.
  !
  ! With each loop in the k direction, we will overwrite the old
  ! data in the planar arrays.
  allocate(sxm(st_lo(1):st_hi(1),st_lo(2):st_hi(2),st_lo(3):st_hi(3)))
  allocate(sxp(st_lo(1):st_hi(1),st_lo(2):st_hi(2),st_lo(3):st_hi(3)))
  allocate ( qxm(It_lo(1):It_hi(1),It_lo(2):It_hi(2),It_lo(3):It_hi(3),NQ) )
  allocate ( qxp(It_lo(1):It_hi(1),It_lo(2):It_hi(2),It_lo(3):It_hi(3),NQ) )
#if BL_SPACEDIM >= 2
  allocate(sym(st_lo(1):st_hi(1),st_lo(2):st_hi(2),st_lo(3):st_hi(3)))
  allocate(syp(st_lo(1):st_hi(1),st_lo(2):st_hi(2),st_lo(3):st_hi(3)))
  allocate ( qym(It_lo(1):It_hi(1),It_lo(2):It_hi(2),It_lo(3):It_hi(3),NQ) )
  allocate ( qyp(It_lo(1):It_hi(1),It_lo(2):It_hi(2),It_lo(3):It_hi(3),NQ) )
#endif
#if BL_SPACEDIM == 3
  allocate(szm(st_lo(1):st_hi(1),st_lo(2):st_hi(2),st_lo(3):st_hi(3)))
  allocate(szp(st_lo(1):st_hi(1),st_lo(2):st_hi(2),st_lo(3):st_hi(3)))
  allocate ( qzm(It_lo(1):It_hi(1),It_lo(2):It_hi(2),It_lo(3):It_hi(3),NQ) )
  allocate ( qzp(It_lo(1):It_hi(1),It_lo(2):It_hi(2),It_lo(3):It_hi(3),NQ) )
#endif

  do n = 1, NQ

    call compute_reconstruction_tvd(q(:,:,:,n), q_lo, q_hi, &
                         sxm, sxp, &
#if BL_SPACEDIM >= 2
                         sym, syp, &
#endif
#if BL_SPACEDIM == 3
                         szm, szp, &
#endif
                         st_lo, st_hi, &
                         lo, hi, dx)

    ! Construct the interface states -- this is essentially just a
    ! reshuffling of interface states from zone-center indexing to
    ! edge-centered indexing
#if BL_SPACEDIM <= 2
    k = 0
#else
    do k = lo(3)-dg(3), hi(3)+dg(3)
#endif
        do j = lo(2)-dg(2), hi(2)+dg(2)
           do i = lo(1)-1, hi(1)+1

              ! x-edges

              ! left state at i-1/2 interface
              qxm(i,j,k,n) = sxp(i-1,j,k)

              ! right state at i-1/2 interface
              qxp(i,j,k,n) = sxm(i,j,k)
#if BL_SPACEDIM >= 2
              ! y-edges

              ! left state at j-1/2 interface
              qym(i,j,k,n) = syp(i,j-1,k)

              ! right state at j-1/2 interface
              qyp(i,j,k,n) = sym(i,j,k)
#endif
#if BL_SPACEDIM == 3
              ! z-edges

              ! left state at k3d-1/2 interface
              qzm(i,j,k,n) = szp(i,j,k-1)

              ! right state at k3d-1/2 interface
              qzp(i,j,k,n) = szm(i,j,k)
#endif
           enddo
        enddo
      enddo
#if BL_SPACEDIM == 3
    enddo
#endif

    ! Compute F^x
    call cmpflx(level, qxm, qxp, It_lo, It_hi, &
               flux1, flux1_lo, flux1_hi, &
               qaux, qa_lo, qa_hi, &
               1, lo, hi, domlo, domhi)
#if BL_SPACEDIM >= 2
    ! Compute F^y
    call cmpflx(level, qym, qyp, It_lo, It_hi, &
               flux2, flux2_lo, flux2_hi, &
               qaux, qa_lo, qa_hi, &
               2, lo, hi, domlo, domhi)
#endif
#if BL_SPACEDIM == 3
    ! Compute F^z
    call cmpflx(level, qzm, qzp, It_lo, It_hi, &
            flux3, flux3_lo, flux3_hi, &
            qaux, qa_lo, qa_hi, &
            3, lo, hi, domlo, domhi)
#endif

    deallocate(sxm, sxp, qxm, qxp)
#if BL_SPACEDIM >= 2
    deallocate(sym, syp, qym, qyp)
#endif
#if BL_SPACEDIM == 3
    deallocate(szm, szp, qzm, qzp)
#endif

  ! For hydro, we will create an update source term that is
  ! essentially the flux divergence.  This can be added with dt to
  ! get the update
  do n = 1, NVAR
#if BL_SPACEDIM <=2
     k = 0
#else
     do k = lo(3), hi(3)
#endif
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
#if BL_SPACEDIM == 1
              update(i,j,k,n) = update(i,j,k,n) + &
                   ( flux1(i,j,k,n) * area1(i,j,k) - flux1(i+1,j,k,n) * area1(i+1,j,k)) / vol(i,j,k)
#elif BL_SPACEDIM == 2
              update(i,j,k,n) = update(i,j,k,n) + &
                 ( flux1(i,j,k,n) * area1(i,j,k) - flux1(i+1,j,k,n) * area1(i+1,j,k) + &
                   flux2(i,j,k,n) * area2(i,j,k) - flux2(i,j+1,k,n) * area2(i,j+1,k) ) / vol(i,j,k)
#else
              update(i,j,k,n) = update(i,j,k,n) + &
                   (flux1(i,j,k,n) * area1(i,j,k) - flux1(i+1,j,k,n) * area1(i+1,j,k) + &
                    flux2(i,j,k,n) * area2(i,j,k) - flux2(i,j+1,k,n) * area2(i,j+1,k) + &
                    flux3(i,j,k,n) * area3(i,j,k) - flux3(i,j,k+1,n) * area3(i,j,k+1) ) / vol(i,j,k)
#endif
              ! for storage
              update_flux(i,j,k,n) = update_flux(i,j,k,n) + &
                   stage_weight * update(i,j,k,n)

              update(i,j,k,n) = update(i,j,k,n) + srcU(i,j,k,n)

              if (update(i,j,k,n) /= update(i,j,k,n)) then
                  update(i,j,k,n) = 0.0e0_rt
              endif

           enddo
        enddo
     enddo
#if BL_SPACEDIM == 3
  enddo
#endif

  ! Scale the fluxes for the form we expect later in refluxing.

  do n = 1, NVAR
#if BL_SPACEDIM <=2
     k = 0
#else
     do k = lo(3), hi(3)
#endif
        do j = lo(2), hi(2)
           do i = lo(1), hi(1) + 1
              flux1(i,j,k,n) = dt * flux1(i,j,k,n) * area1(i,j,k)
              if (flux1(i,j,k,n) /= flux1(i,j,k,n)) then
                  flux1(i,j,k,n) = 0.0e0_rt
              endif
           enddo
        enddo
     enddo
#if BL_SPACEDIM == 3
  enddo
#endif

#if BL_SPACEDIM >= 2
  do n = 1, NVAR
#if BL_SPACEDIM == 2
     k = 0
#else
     do k = lo(3), hi(3)
#endif
        do j = lo(2), hi(2) + 1
           do i = lo(1), hi(1)
              flux2(i,j,k,n) = dt * flux2(i,j,k,n) * area2(i,j,k)
              if (flux2(i,j,k,n) /= flux2(i,j,k,n)) then
                  flux2(i,j,k,n) = 0.0e0_rt
              endif
           enddo
        enddo
     enddo
#if BL_SPACEDIM == 3
  enddo
#endif
#endif

#if BL_SPACEDIM == 3
  do n = 1, NVAR
     do k = lo(3), hi(3) + 1
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              flux3(i,j,k,n) = dt * flux3(i,j,k,n) * area3(i,j,k)
              if (flux3(i,j,k,n) /= flux3(i,j,k,n)) then
                  flux3(i,j,k,n) = 0.0e0_rt
              endif
           enddo
        enddo
     enddo
  enddo
#endif

! write(*,*) flux1(lo(1), :, lo(3), :)
end subroutine ca_mol_single_stage
