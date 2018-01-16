module advection_util_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  private

  public enforce_consistent_e, enforce_minimum_density, ca_compute_cfl, swectoprim, compctoprim

contains

    subroutine enforce_consistent_e(lo,hi,state,s_lo,s_hi,level,xlo,dx)

      use meth_params_module, only: NVAR, QRHO, QU, QV, QW, UEDEN, UEINT, NQ, NQAUX
      use bl_constants_module, only: HALF, ONE
      use amrex_fort_module, only: rt => amrex_real
      use probdata_module, only: swe_to_comp_level

      implicit none

      integer,  intent(in   ) :: lo(3), hi(3), level
      integer,  intent(in   ) :: s_lo(3), s_hi(3)
      real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

      ! Local variables
      integer  :: i,j,k
      real(rt) :: u, v, w, rhoInv, xlo(3), dx(3)
      real(rt) :: q(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NQ)
      real(rt) :: qaux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),NQAUX)

      if (level <= swe_to_comp_level) then
          call swectoprim(lo, hi, &
                            state, s_lo, s_hi, &
                            q,  s_lo, s_hi, &
                            qaux, lo, hi)
      else
          call compctoprim(lo, hi, &
                            state, s_lo, s_hi, &
                            q,  s_lo, s_hi, &
                            qaux, lo, hi, xlo, dx)
      endif

      !
      ! Enforces (rho E) = (rho e) + 1/2 rho (u^2 + v^2 + w^2)
      !
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               u = q(i,j,k,QU)
               v = q(i,j,k,QV)
               w = q(i,j,k,QW)

               state(i,j,k,UEDEN) = state(i,j,k,UEINT) + &
                    HALF * q(i,j,k,QRHO) * (u*u + v*v + w*w)

            end do
         end do
      end do

  end subroutine enforce_consistent_e

  subroutine enforce_minimum_density(uin,uin_lo,uin_hi, &
                                     uout,uout_lo,uout_hi, &
                                     vol,vol_lo,vol_hi, &
                                     lo,hi,frac_change,verbose,level, &
                                     xlo, dx)

    use network, only : nspec, naux
    use meth_params_module, only : NVAR, QRHO, QREINT, UEDEN, small_dens, density_reset_method, NQ, NQAUX
    use bl_constants_module, only : ZERO
    use riemann_util_module, only : swe_cons_state, comp_cons_state
    use probdata_module, only: swe_to_comp_level

    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: lo(3), hi(3), verbose,level
    integer, intent(in) ::  uin_lo(3),  uin_hi(3)
    integer, intent(in) :: uout_lo(3), uout_hi(3)
    integer, intent(in) ::  vol_lo(3),  vol_hi(3)

    real(rt)        , intent(inout) ::  uin( uin_lo(1): uin_hi(1), uin_lo(2): uin_hi(2), uin_lo(3): uin_hi(3),NVAR)
    real(rt)        , intent(inout) :: uout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2),uout_lo(3):uout_hi(3),NVAR)
    real(rt)        , intent(in) ::  vol( vol_lo(1): vol_hi(1), vol_lo(2): vol_hi(2), vol_lo(3): vol_hi(3))
    real(rt)        , intent(inout) :: frac_change
    real(rt), intent(in) :: xlo(3), dx(3)

    ! Local variables
    integer          :: i,ii,j,jj,k,kk
    integer          :: i_set, j_set, k_set
    real(rt)         :: max_dens
    real(rt)         :: qnew(NQ)
    integer          :: num_positive_zones

    logical :: have_reset

    real(rt)     :: qin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NQ)
    real(rt)   :: qaux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),NQAUX)
    real(rt)     :: qout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2),uout_lo(3):uout_hi(3),NQ)

    !return

    max_dens = ZERO

    have_reset = .false.

    ! write(*,*) "level = ", level

    if (level <= swe_to_comp_level) then
        call swectoprim(lo, hi, &
                          uin, uin_lo, uin_hi, &
                          qin,  uin_lo, uin_hi, &
                          qaux, lo, hi)
        call swectoprim(lo, hi, &
                        uout, uout_lo, uout_hi, &
                        qout,  uout_lo, uout_hi, &
                        qaux,  lo, hi)
    else
        call compctoprim(lo, hi, &
                          uin, uin_lo, uin_hi, &
                          qin,  uin_lo, uin_hi, &
                          qaux, lo, hi, xlo, dx)
        call compctoprim(lo, hi, &
                        uout, uout_lo, uout_hi, &
                        qout,  uout_lo, uout_hi, &
                        qaux,  lo, hi, xlo, dx)
    endif

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             if (qout(i,j,k,QRHO) .eq. ZERO) then

                print *,'DENSITY EXACTLY ZERO AT CELL ',i,j,k
                print *,'  in grid ',lo(1),lo(2),lo(3),hi(1),hi(2),hi(3)
                call bl_error("Error:: advection_util_nd.f90 :: ca_enforce_minimum_density")

                qout(i,j,k,QRHO) = small_dens

            else if (qout(i,j,k,QRHO) < small_dens) then

                have_reset = .true.

                ! Store the maximum (negative) fractional change in the density

                if ( qout(i,j,k,QRHO) < ZERO .and. &
                     (qout(i,j,k,QRHO) - qin(i,j,k,QRHO)) / qin(i,j,k,QRHO) < frac_change) then

                   frac_change = (qout(i,j,k,QRHO) - qin(i,j,k,QRHO)) / qin(i,j,k,QRHO)

                endif

                if (density_reset_method == 1) then

                   ! Reset to the characteristics of the adjacent state with the highest density.

                   max_dens = qout(i,j,k,QRHO)
                   i_set = i
                   j_set = j
                   k_set = k
                   do kk = -1,1
                      do jj = -1,1
                         do ii = -1,1
                            if (i+ii.ge.lo(1) .and. j+jj.ge.lo(2) .and. k+kk.ge.lo(3) .and. &
                                 i+ii.le.hi(1) .and. j+jj.le.hi(2) .and. k+kk.le.hi(3)) then
                               if (qout(i+ii,j+jj,k+kk,QRHO) .gt. max_dens) then
                                  i_set = i+ii
                                  j_set = j+jj
                                  k_set = k+kk
                                  max_dens = qout(i_set,j_set,k_set,QRHO)
                               endif
                            endif
                         end do
                      end do
                   end do

                   if (max_dens < small_dens) then

                      ! We could not find any nearby zones with sufficient density.

                      call reset_to_small_state(qin(i,j,k,:), qout(i,j,k,:), [i, j, k], lo, hi, verbose)

                   else

                      qnew = qout(i_set,j_set,k_set,:)

                      call reset_to_zone_state(qin(i,j,k,:), qout(i,j,k,:), qnew(:), [i, j, k], lo, hi, verbose)

                   endif

                else if (density_reset_method == 2) then

                   ! Reset to the average of adjacent zones. The median is independently calculated for each variable.

                   num_positive_zones = 0
                   qnew(:) = ZERO

                   do kk = -1, 1
                      do jj = -1, 1
                         do ii = -1, 1
                            if (i+ii.ge.lo(1) .and. j+jj.ge.lo(2) .and. k+kk.ge.lo(3) .and. &
                                i+ii.le.hi(1) .and. j+jj.le.hi(2) .and. k+kk.le.hi(3)) then
                               if (qout(i+ii,j+jj,k+kk,QRHO) .ge. small_dens) then
                                  qnew(:) = qnew(:) + qout(i+ii,j+jj,k+kk,:)
                                  num_positive_zones = num_positive_zones + 1
                               endif
                            endif
                         enddo
                      enddo
                   enddo

                   if (num_positive_zones == 0) then

                      ! We could not find any nearby zones with sufficient density.

                      call reset_to_small_state(qin(i,j,k,:), qout(i,j,k,:), [i, j, k], lo, hi, verbose)

                   else

                      qnew(:) = qnew(:) / num_positive_zones

                      call reset_to_zone_state(qin(i,j,k,:), qout(i,j,k,:), qnew(:), [i, j, k], lo, hi, verbose)

                   endif

                elseif (density_reset_method == 3) then

                   ! Reset to the original zone state.

                   if (qin(i,j,k,QRHO) < small_dens) then

                      call reset_to_small_state(qin(i,j,k,:), qout(i,j,k,:), [i, j, k], lo, hi, verbose)

                   else

                      qnew(:) = qin(i,j,k,:)

                      call reset_to_zone_state(qin(i,j,k,:), qout(i,j,k,:), qnew(:), [i, j, k], lo, hi, verbose)

                   endif

                else

                   call bl_error("Unknown density_reset_method in subroutine ca_enforce_minimum_density.")

                endif

                if (level <= swe_to_comp_level) then
                    call swe_cons_state(qout(i,j,k,:), uout(i,j,k,:))
                else
                    call comp_cons_state(qout(i,j,k,:), uout(i,j,k,:))
                endif

             end if

          enddo
       enddo
    enddo

  end subroutine enforce_minimum_density



  subroutine reset_to_small_state(old_state, new_state, idx, lo, hi, verbose)

    use bl_constants_module, only: ZERO
    use network, only: nspec, naux
    use meth_params_module, only: NQ, QRHO, QU, QV, QW, QFS, small_temp, small_dens, npassive, upass_map
    use castro_util_module, only: position
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    real(rt)         :: old_state(NQ), new_state(NQ)
    integer          :: idx(3), lo(3), hi(3), verbose

    integer          :: n, ipassive

    ! If no neighboring zones are above small_dens, our only recourse
    ! is to set the density equal to small_dens, and the temperature
    ! equal to small_temp. We set the velocities to zero,
    ! though any choice here would be arbitrary.

    return

    if (verbose .gt. 0) then
       print *,'   '
       if (new_state(QRHO) < ZERO) then
          print *,'>>> RESETTING NEG.  DENSITY AT ',idx(1),idx(2),idx(3)
       else
          print *,'>>> RESETTING SMALL DENSITY AT ',idx(1),idx(2),idx(3)
       endif
       print *,'>>> FROM ',new_state(QRHO),' TO ',small_dens
       print *,'>>> IN GRID ',lo(1),lo(2),lo(3),hi(1),hi(2),hi(3)
       print *,'>>> ORIGINAL DENSITY FOR OLD STATE WAS ',old_state(QRHO)
       print *,'   '
    end if

    do ipassive = 1, npassive
       n = upass_map(ipassive)
       new_state(n) = new_state(n) * (small_dens / new_state(QRHO))
    end do

    new_state(QRHO ) = small_dens

    new_state(QU  ) = ZERO
    new_state(QV  ) = ZERO
    new_state(QW  ) = ZERO

  end subroutine reset_to_small_state


  subroutine reset_to_zone_state(old_state, new_state, input_state, idx, lo, hi, verbose)

    use bl_constants_module, only: ZERO
    use meth_params_module, only: NQ, QRHO

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    real(rt)         :: old_state(NQ), new_state(NQ), input_state(NQ)
    integer          :: idx(3), lo(3), hi(3), verbose

    return

    if (verbose .gt. 0) then
       if (new_state(QRHO) < ZERO) then
          print *,'   '
          print *,'>>> RESETTING NEG.  DENSITY AT ',idx(1),idx(2),idx(3)
          print *,'>>> FROM ',new_state(QRHO),' TO ',input_state(QRHO)
          print *,'>>> IN GRID ',lo(1),lo(2),lo(3),hi(1),hi(2),hi(3)
          print *,'>>> ORIGINAL DENSITY FOR OLD STATE WAS ',old_state(QRHO)
          print *,'   '
       else
          print *,'   '
          print *,'>>> RESETTING SMALL DENSITY AT ',idx(1),idx(2),idx(3)
          print *,'>>> FROM ',new_state(QRHO),' TO ',input_state(QRHO)
          print *,'>>> IN GRID ',lo(1),lo(2),lo(3),hi(1),hi(2),hi(3)
          print *,'>>> ORIGINAL DENSITY FOR OLD STATE WAS ',old_state(QRHO)
          print *,'   '
       end if
    end if

    new_state(:) = input_state(:)

  end subroutine reset_to_zone_state


  subroutine ca_compute_cfl(lo, hi, q, q_lo, q_hi, &
                         qaux, qa_lo, qa_hi, &
                         dt, dx, courno) &
                         bind(C, name = "ca_compute_cfl")

    use bl_constants_module, only: ZERO, ONE
    use meth_params_module, only: NQ, QRHO, QU, QV, QW, QC, NQAUX
    use prob_params_module, only: dim
    use probdata_module, only : g

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer :: lo(3), hi(3)
    integer :: q_lo(3), q_hi(3), qa_lo(3), qa_hi(3)

    real(rt)         :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt)         :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt)         :: dt, dx(3), courno

    real(rt)         :: courx, coury, courz, courmx, courmy, courmz, courtmp
    real(rt)         :: dtdx, dtdy, dtdz, c, large_number
    integer          :: i, j, k

    ! Compute running max of Courant number over grids

    courmx = courno
    courmy = courno
    courmz = courno

    dtdx = dt / dx(1)

    if (dim .ge. 2) then
       dtdy = dt / dx(2)
    else
       dtdy = ZERO
    endif

    if (dim .eq. 3) then
       dtdz = dt / dx(3)
    else
       dtdz = ZERO
    endif

    large_number = 1.0e10_rt

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

              if (q(i,j,k,QRHO) .le. ZERO) then
                 !print *, uin(:,:,:,URHO)
                 print *,'   '
                 print *,'>>> Error: advection_util_nd.F90::compute_cfl ',i, j, k
                 print *,'>>> ... negative density ', q(i,j,k,QRHO)
                 q(i,j,k,QRHO) = 1.0_rt
                 call bl_error("Error:: advection_util_nd.f90 :: compute_cfl")
             else if (q(i,j,k,QRHO) > large_number) then
                    !print *, uin(:,:,:,URHO)
                    print *,'   '
                    print *,'>>> Error: advection_util_nd.F90::compute_cfl ',i, j, k
                    print *,'>>> ... density is very large ', q(i,j,k,QRHO)
                    q(i,j,k,QRHO) = 1.0_rt
                    call bl_error("Error:: advection_util_nd.f90 :: compute_cfl")
             else if (q(i,j,k,QRHO) /= q(i,j,k,QRHO)) then
                  print *,'   '
                  print *,'>>> Error: advection_util_nd.F90::compute_cfl ',i, j, k
                  print *,'>>> ... density is nan ', q(i,j,k,QRHO)
                  write(*,*) q(i,j,k,QRHO)
                  call bl_error("Error:: advection_util_nd.f90 :: compute_cfl")
              endif

              ! stops weird nanning
              if (maxval(abs(q(i,j,k,QU:QW))) > large_number) then
                  q(i,j,k,QU:QW) = ZERO
              endif

              c = sqrt(q(i,j,k,QRHO) * g)

             courx = ( c + abs(q(i,j,k,QU)) ) * dtdx
             coury = ( c + abs(q(i,j,k,QV)) ) * dtdy
             courz = ( c + abs(q(i,j,k,QW)) ) * dtdz

             courmx = max( courmx, courx )
             courmy = max( courmy, coury )
             courmz = max( courmz, courz )

            ! method-of-lines constraint
            courtmp = courx
            if (dim >= 2) then
               courtmp = courtmp + coury
            endif
            if (dim == 3) then
               courtmp = courtmp + courz
            endif

            ! note: it might not be 1 for all RK integrators
            if (courtmp > ONE) then
               print *,'   '
               call bl_warning("Warning:: advection_util_nd.F90 :: CFL violation in compute_cfl")
               print *,'>>> ... at cell (i,j,k)   : ', i, j, k
               print *,'>>> ... u,v,w, c            ', q(i,j,k,QU), q(i,j,k,QV), q(i,j,k,QW), c
               print *,'>>> ... density             ', q(i,j,k,QRHO)
               print *,'>>> ... courtmp            ', courtmp
               call bl_error("Warning:: advection_util_nd.F90 :: CFL violation in compute_cfl")
            endif

            courno = max(min(1.0_rt,courno), min(1.0_rt,courtmp))

            if (courno > ONE) then
                print *,'   '
                call bl_warning("Warning:: advection_util_nd.F90 :: CFL violation in compute_cfl")
                print *,'>>> ... at cell (i,j,k)   : ', i, j, k
                print *,'>>> ... u,v,w, c            ', q(i,j,k,QU), q(i,j,k,QV), q(i,j,k,QW), c
                print *,'>>> ... density             ', q(i,j,k,QRHO)
                print *,'>>> ... courno            ', courno
            endif
          enddo
       enddo
    enddo

  end subroutine ca_compute_cfl

  subroutine swectoprim(lo, hi, &
                     uin, uin_lo, uin_hi, &
                     q,     q_lo,   q_hi, &
                     qaux, qa_lo,  qa_hi, ignore_errors)

    use actual_network, only : nspec, naux
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, &
                                   QRHO, QU, QV, QW, QTEMP, UTEMP, &
                                   NQ, QC, QCSML, QGAMC, QDPDR, QDPDE, NQAUX, QPRES, QREINT, UEINT, &
                                   npassive, upass_map, qpass_map, &
                                   small_dens, QFA, QFS, UFA
    use bl_constants_module, only: ZERO, HALF, ONE
    use castro_util_module, only: position
    use probdata_module, only : g, dens_incompressible

    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: uin_lo(3), uin_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)

    real(rt)        , intent(in) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
    real(rt)        , intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt)        , intent(inout) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    logical, optional, intent(in) :: ignore_errors

    real(rt)        , parameter :: small = 1.0e-8_rt

    integer          :: i, j, k, ii, jj, kk
    integer          :: n, iq, ipassive
    real(rt)         :: W

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
              if (present(ignore_errors) .and. (.not. ignore_errors)) then
                 if (uin(i,j,k,URHO) .le. ZERO) then
                     do kk = lo(3), hi(3)
                        do jj = lo(2), hi(2)
                           do ii = lo(1), hi(1)
                               print *, uin(ii,jj,kk,1:NVAR)
                           enddo
                       enddo
                   enddo
                    print *,'   '
                    print *,'>>> Error: advection_util_nd.F90::swectoprim ',i, j, k
                    print *,'>>> ... negative density ', uin(i,j,k,URHO)
                    ! uin(i,j,k,:) = 0.0_rt
                    ! uin(i,j,k,URHO) = 1.0_rt
                    call bl_error("Error:: advection_util_nd.f90 :: swectoprim")
                 else if (uin(i,j,k,URHO) /= uin(i,j,k,URHO)) then
                     print *,'   '
                     print *,'>>> Error: advection_util_nd.F90::swectoprim ',i, j, k
                     print *,'>>> ... density is nan ', uin(i,j,k,URHO)
                     write(*,*) uin(:,:,:,URHO)
                     call bl_error("Error:: advection_util_nd.f90 :: swectoprim")
                 else if (uin(i,j,k,URHO) .lt. small_dens) then
                     do kk = lo(3), hi(3)
                        do jj = lo(2), hi(2)
                           do ii = lo(1), hi(1)
                               print *, uin(ii,jj,kk,1:NVAR)
                           enddo
                       enddo
                   enddo
                    print *,'   '
                    print *,'>>> Error: advection_util_nd.F90::swectoprim ',i, j, k
                    print *,'>>> ... small density ', uin(i,j,k,URHO)
                    call bl_error("Error:: advection_util_nd.f90 :: swectoprim")
                 endif
             endif
          end do

          do i = lo(1), hi(1)
              q(i,j,k,1:NQ) = 0.0_rt

              W = sqrt(1.0_rt + sum(uin(i,j,k,UMX:UMZ)**2) / uin(i,j,k,URHO)**2)

              !q(i,j,k,:QW) = uin(i,j,k,:QW)
              if (uin(i,j,k,URHO) .le. ZERO) then
                  q(i,j,k,QRHO) = 1.0_rt / W
                  ! q(i,j,k,QU) = 0.1_rt * uin(i,j,k,UMX)
                  q(i,j,k,QV) = 0.1_rt * uin(i,j,k,UMY) / W
                  q(i,j,k,QW) = 0.1_rt * uin(i,j,k,UMZ) / W
              else
                  q(i,j,k,QRHO) = uin(i,j,k,URHO)
                  ! q(i,j,k,QU) = uin(i,j,k,UMX) / uin(i,j,k,URHO)
                  q(i,j,k,QV) = uin(i,j,k,UMY) / (uin(i,j,k,URHO) * W)
                  q(i,j,k,QW) = uin(i,j,k,UMZ) / (uin(i,j,k,URHO) * W)
              endif
#if BL_SPACEDIM == 1
              q(i,j,k,QV) = 0.0_rt
#endif
#if BL_SPACEDIM <= 2
              q(i,j,k,QW) = 0.0_rt
#endif
              q(i,j,k,QPRES) = 0.5_rt * dens_incompressible * g * q(i,j,k,QRHO)**2
              q(i,j,k,QREINT) = uin(i,j,k,UEINT)
              q(i,j,k,QTEMP) = uin(i,j,k,UTEMP)
              q(i,j,k,QFA) = uin(i,j,k,UFA) / uin(i,j,k,URHO)
              q(i,j,k,QFS:QFS-1+nspec) = q(i,j,k,QRHO) / nspec

          enddo
       enddo
    enddo

    !Load passively advected quantities into q
      ! do ipassive = 1, npassive
      !    n  = upass_map(ipassive)
      !    iq = qpass_map(ipassive)
      !    do k = lo(3),hi(3)
      !       do j = lo(2),hi(2)
      !          do i = lo(1),hi(1)
      !             q(i,j,k,iq) = uin(i,j,k,n)/q(i,j,k,QRHO)
      !          enddo
      !       enddo
      !    enddo
      ! enddo

end subroutine swectoprim

subroutine compctoprim(lo, hi, &
                   uin, uin_lo, uin_hi, &
                   q,     q_lo,   q_hi, &
                   qaux, qa_lo,  qa_hi, xlo, dx, ignore_errors)

  use actual_network, only : nspec, naux
  use eos_module, only : eos, eos_init, initialized
  use eos_type_module, only : eos_t, eos_input_rp, eos_input_re
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT,&
                                 QRHO, QU, QV, QW, QREINT, QTEMP, &
                                 NQ, QC, QCSML, QGAMC, QDPDR, QDPDE, NQAUX, QFS, QFX, QGAME, QPRES, UTEMP, &
                                 npassive, upass_map, qpass_map, &
                                 small_dens, small_temp, UFA, QFA
  use bl_constants_module, only: ZERO, HALF, ONE
  use castro_util_module, only: position
  use probdata_module, only : g, dens_incompressible
  use riemann_util_module, only : zbrent, f_of_p

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: uin_lo(3), uin_hi(3)
  integer, intent(in) :: q_lo(3), q_hi(3)
  integer, intent(in) :: qa_lo(3), qa_hi(3)

  real(rt)        , intent(in ) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
  real(rt)        , intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
  real(rt)        , intent(inout) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
  real(rt), intent(in) :: dx(3), xlo(3)

  logical, optional, intent(in) :: ignore_errors

  real(rt)        , parameter :: small = 1.0e-8_rt

  integer          :: i, j, k
  integer          :: n, iq, ipassive
  real(rt)         :: xx, pmin, pmax, sq, W2, h, ssq, p, rhoh
  real(rt)         :: fmin, fmax, eden, vel(3)
  type (eos_t)     :: eos_state

  eos_state % rho = dens_incompressible
  eos_state % e = 1.0_rt

  if (.not. initialized) call eos_init(small_dens=small_dens, small_temp=small_temp)

  call eos(eos_input_re, eos_state)

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           if (present(ignore_errors) .and. (.not. ignore_errors)) then
               if (uin(i,j,k,URHO) .le. ZERO) then
                  print *,'   '
                  print *,'>>> Error: advection_util_nd.F90::compctoprim ',i, j, k
                  print *,'>>> ... negative density ', uin(i,j,k,URHO)
                  call bl_error("Error:: advection_util_nd.f90 :: compctoprim")
               else if (uin(i,j,k,URHO) /= uin(i,j,k,URHO)) then
                   print *,'   '
                   print *,'>>> Error: advection_util_nd.F90::compctoprim ',i, j, k
                   print *,'>>> ... density is nan ', uin(i,j,k,URHO)
                   write(*,*) uin(:,:,:,URHO)
                   call bl_error("Error:: advection_util_nd.f90 :: compctoprim")
               else if (uin(i,j,k,URHO) .lt. small_dens) then
                  print *,'   '
                  print *,'>>> Error: advection_util_nd.F90::compctoprim ',i, j, k
                  print *,'>>> ... small density ', uin(i,j,k,URHO)
                  call bl_error("Error:: advection_util_nd.f90 :: compctoprim")
               endif
           endif
        end do

        do i = lo(1), hi(1)
            q(i,j,k,1:NQ) = 0.0_rt

            ! if (uin(i,j,k,URHO) .le. ZERO) then
            !     q(i,j,k,QRHO) = 1.0_rt
            !     q(i,j,k,QU) = 0.1_rt * uin(i,j,k,UMX) !/ uin(i,j,k,URHO)
            !     q(i,j,k,QV) = 0.1_rt * uin(i,j,k,UMY) !/ uin(i,j,k,URHO)
            !     q(i,j,k,QW) = 0.1_rt * uin(i,j,k,UMZ) !/ uin(i,j,k,URHO)
            ! else
            !     ! Incompressible
            !     q(i,j,k,QRHO) = 1.0_rt!uin(i,j,k,URHO)
            !     q(i,j,k,QU) = uin(i,j,k,UMX) / uin(i,j,k,URHO)
            !     q(i,j,k,QV) = uin(i,j,k,UMY) / uin(i,j,k,URHO)
            !     q(i,j,k,QW) = uin(i,j,k,UMZ) / uin(i,j,k,URHO)
            ! endif

            ! eden = uin(i,j,k,UEDEN)
            !
            ! if (eden < ZERO) then
            !     eden = abs(eden)
            ! endif
            !
            ! ssq = sum(uin(i,j,k,UMX:UMZ)**2)
            !
            ! pmin = (gamma - 1.0_rt) * (eden + uin(i,j,k,URHO) - ssq / (eden + uin(i,j,k,URHO))**2 - uin(i,j,k,URHO))
            !
            ! pmax = (gamma - 1.0_rt) * (eden + uin(i,j,k,URHO) - ssq / (eden + uin(i,j,k,URHO)))
            !
            ! if (pmin < 0.0_rt) then
            !       pmin = 0._rt
            ! end if
            !
            ! call f_of_p(fmin, pmin, uin(i,j,k,:))
            ! call f_of_p(fmax, pmax, uin(i,j,k,:))
            !
            ! if (fmin * fmax > 0.0_rt) then
            !     pmin = pmin * 0.1_rt !0._rt
            ! end if
            !
            ! call f_of_p(fmin, pmin, uin(i,j,k,:))
            !
            ! !write(*,*) "f = ", fmin, fmax, " p = ", pmin, pmax
            !
            ! if (fmin * fmax > 0.0_rt) then
            !   pmax = pmax * 10._rt
            ! end if
            !
            ! call zbrent(p, pmin, pmax, uin(i,j,k,:))
            !
            ! !write(*,*) "pressure = ", pmin, pmax
            !
            ! if (p /= p .or. p <= 0.0_rt) then! .or. p > 1.0_rt) then
            !
            !   p = abs((gamma - 1.0_rt) * ((eden + uin(i,j,k,URHO)) - ssq / (eden + uin(i,j,k,URHO))**2 - uin(i,j,k,URHO)))
            !
            !   !if (p > 1.0_rt) then
            !       !p = 1.0_rt
            !   !end if
            ! end if
            !
            ! sq = sqrt((eden + p + uin(i,j,k,URHO))**2 - ssq)
            !
            ! if (sq /= sq) then
            !   sq = eden + p + uin(i,j,k,URHO)
            ! end if
            !
            ! h = 1.0_rt + gamma * &
            ! (sq - p * (eden + p + uin(i,j,k,URHO)) / sq - uin(i,j,k,URHO)) / uin(i,j,k,URHO)
            ! W2 = 1.0_rt + ssq / (uin(i,j,k,URHO) * h)**2

            !write(*,*) "p, sq, eden, rho", p, sq, eden, uin(i,j,k,URHO)
            !return

            q(i,j,k,QRHO) = dens_incompressible! INCOMPRESSIBLE uin(i,j,k,URHO) * sq / (eden + p + uin(i,j,k,URHO))

            ! write(*,*) uin(i,j,k,URHO) * sq / (eden + p + uin(i,j,k,URHO))

            W2 = (uin(i,j,k,URHO) / dens_incompressible)**2
            ! if ((W2 - 1.0_rt) < 1.0e-8_rt) then
                ! eden = uin(i,j,k,UEDEN)
                ! p = (gamma - 1.0_rt) * (eden + uin(i,j,k,URHO) - ssq / (eden + uin(i,j,k,URHO))**2 - uin(i,j,k,URHO))
            p = (uin(i,j,k,UEDEN) - dens_incompressible * W2 + uin(i,j,k,URHO)) / (eos_state % gam1 * W2 / (eos_state % gam1 - 1.0_rt) - 1.0_rt)
            rhoh = dens_incompressible + eos_state % gam1 * p / (eos_state % gam1 - 1.0_rt)
                ! rhoh = uin(i,j,k,UEINT) * eos_state % gam1 + dens_incompressible
            ! else
            !     rhoh = ssq / (W2 - 1.0_rt)
            ! endif

            vel(1) = uin(i,j,k,UMX) / (W2 * rhoh)
            vel(2) = uin(i,j,k,UMY) / (W2 * rhoh)
            vel(3) = uin(i,j,k,UMZ) / (W2 * rhoh)
#if BL_SPACEDIM == 1
            vel(2) = 0.0_rt
#endif
#if BL_SPACEDIM <= 2
            vel(3) = 0.0_rt
#endif

            ! p = rhoh * W2 - uin(i,j,k,UEDEN) - uin(i,j,k,URHO)

            q(i,j,k,QU:QW) = vel(1:3)

            q(i,j,k,QREINT) = (rhoh - q(i,j,k,QRHO)) / eos_state % gam1!p / (gamma - 1.0_rt)

            q(i,j,k,QTEMP) = uin(i,j,k,UTEMP)

            q(i,j,k,QPRES) = p

            q(i,j,k,QFA) = uin(i,j,k,UFA) / uin(i,j,k,URHO)
            q(i,j,k,QFS:QFS-1+nspec) = q(i,j,k,QRHO) / nspec
        enddo
     enddo
  enddo

  ! ! Load passively advected quatities into q
  !   do ipassive = 1, npassive
  !      n  = upass_map(ipassive)
  !      iq = qpass_map(ipassive)
  !      do k = lo(3),hi(3)
  !         do j = lo(2),hi(2)
  !            do i = lo(1),hi(1)
  !                if (abs(q(i,j,k,QRHO)) > 1.0e-9_rt) then
  !                    q(i,j,k,iq) = uin(i,j,k,n)/q(i,j,k,QRHO)
  !                else
  !                    q(i,j,k,iq) = 0.0_rt
  !                endif
  !            enddo
  !         enddo
  !      enddo
  !   enddo

    ! get gamc, p, T, c, csml using q state
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             xx = xlo(1) + dx(1)*dble(i-lo(1)+HALF)

             eos_state % T   = q(i,j,k,QTEMP )
             eos_state % rho = q(i,j,k,QRHO  )
             eos_state % e   = q(i,j,k,QREINT)
             eos_state % xn  = q(i,j,k,QFS:QFS+nspec-1)
             eos_state % aux = q(i,j,k,QFX:QFX+naux-1)
             ! eos_state % p = 0.5_rt * g * q(i,j,k,QRHO  )**2
             eos_state % p = q(i,j,k,QPRES)
             !q(i,j,k,QPRES)  = eos_state % p

             call eos(eos_input_rp, eos_state)
             !q(i,j,k,QPRES)  = eos_state % p

             q(i,j,k,QTEMP)  = eos_state % T
             !q(i,j,k,QREINT) = eos_state % e * q(i,j,k,QRHO)
             q(i,j,k,QGAME)  = q(i,j,k,QPRES) / q(i,j,k,QREINT) + ONE

             qaux(i,j,k,QDPDR)  = eos_state % dpdr_e
             qaux(i,j,k,QDPDE)  = eos_state % dpde
             qaux(i,j,k,QGAMC)  = eos_state % gam1
             qaux(i,j,k,QC   )  = eos_state % cs
             ! qaux(i,j,k,QCSML)  = max(small, small * qaux(i,j,k,QC))
          enddo
       enddo
   enddo

end subroutine compctoprim

end module advection_util_module
