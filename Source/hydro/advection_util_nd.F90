module advection_util_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  private

  public enforce_minimum_density, compute_cfl, grctoprim, normalize_species_fluxes

contains

  subroutine enforce_minimum_density(uin,uin_lo,uin_hi, &
                                     uout,uout_lo,uout_hi, &
                                     vol,vol_lo,vol_hi, &
                                     lo,hi,frac_change,verbose)

    use network, only : nspec, naux
    use meth_params_module, only : NVAR, QRHO, QREINT, UEDEN, small_dens, density_reset_method, NQ, NQAUX
    use bl_constants_module, only : ZERO
    use riemann_util_module, only : gr_cons_state
    use metric_module, only : calculate_gamma_up, calculate_alpha

    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: lo(3), hi(3), verbose
    integer, intent(in) ::  uin_lo(3),  uin_hi(3)
    integer, intent(in) :: uout_lo(3), uout_hi(3)
    integer, intent(in) ::  vol_lo(3),  vol_hi(3)

    real(rt)        , intent(in) ::  uin( uin_lo(1): uin_hi(1), uin_lo(2): uin_hi(2), uin_lo(3): uin_hi(3),NVAR)
    real(rt)        , intent(inout) :: uout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2),uout_lo(3):uout_hi(3),NVAR)
    real(rt)        , intent(in) ::  vol( vol_lo(1): vol_hi(1), vol_lo(2): vol_hi(2), vol_lo(3): vol_hi(3))
    real(rt)        , intent(inout) :: frac_change

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
    real(rt) :: gamma_up(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),9)
    real(rt) :: alpha(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

    call calculate_gamma_up(gamma_up, lo, hi)
    call calculate_alpha(alpha, lo, hi)

    max_dens = ZERO

    have_reset = .false.

    call grctoprim(lo, hi, &
                      uin, uin_lo, uin_hi, &
                      qin,  uin_lo, uin_hi, &
                      qaux, lo, hi, &
                      gamma_up, lo, hi, &
                      alpha, lo, hi)
    call grctoprim(lo, hi, &
                    uout, uout_lo, uout_hi, &
                    qout,  uout_lo, uout_hi, &
                    qaux,  lo, hi, &
                    gamma_up, lo, hi, &
                    alpha, lo, hi)

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

                call gr_cons_state(qout(i,j,k,:), uout(i,j,k,:), gamma_up(i,j,k,:))

             end if

          enddo
       enddo
    enddo

  end subroutine enforce_minimum_density



  subroutine reset_to_small_state(old_state, new_state, idx, lo, hi, verbose)

    use bl_constants_module, only: ZERO
    use network, only: nspec, naux
    use meth_params_module, only: NQ, QRHO, QU, QV, QW, QTEMP, QREINT, QFS, small_temp, small_dens, npassive, upass_map
    use eos_type_module, only: eos_t, eos_input_rt
    use eos_module, only: eos
    use castro_util_module, only: position
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    real(rt)         :: old_state(NQ), new_state(NQ)
    integer          :: idx(3), lo(3), hi(3), verbose

    integer          :: n, ipassive
    type (eos_t)     :: eos_state

    ! If no neighboring zones are above small_dens, our only recourse
    ! is to set the density equal to small_dens, and the temperature
    ! equal to small_temp. We set the velocities to zero,
    ! though any choice here would be arbitrary.

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

    eos_state % rho = small_dens
    eos_state % T   = small_temp
    eos_state % xn  = new_state(QFS:QFS+nspec-1)
    eos_state % aux = new_state(QFS:QFS+naux-1)

    call eos(eos_input_rt, eos_state)

    new_state(QRHO ) = eos_state % rho
    new_state(QTEMP) = eos_state % T

    new_state(QU  ) = ZERO
    new_state(QV  ) = ZERO
    new_state(QW  ) = ZERO

    new_state(QREINT) = eos_state % rho * eos_state % e
  end subroutine reset_to_small_state


  subroutine reset_to_zone_state(old_state, new_state, input_state, idx, lo, hi, verbose)

    use bl_constants_module, only: ZERO
    use meth_params_module, only: NQ, QRHO

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    real(rt)         :: old_state(NQ), new_state(NQ), input_state(NQ)
    integer          :: idx(3), lo(3), hi(3), verbose

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



  subroutine compute_cfl(q, q_lo, q_hi, &
                         qaux, qa_lo, qa_hi, &
                         lo, hi, dt, dx, courno) &
                         bind(C, name = "compute_cfl")

    use bl_constants_module, only: ZERO, ONE
    use meth_params_module, only: NQ, QRHO, QU, QV, QW, QC, NQAUX
    use prob_params_module, only: dim

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer :: lo(3), hi(3)
    integer :: q_lo(3), q_hi(3), qa_lo(3), qa_hi(3)

    real(rt)         :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt)         :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt)         :: dt, dx(3), courno

    real(rt)         :: courx, coury, courz, courmx, courmy, courmz, courtmp
    real(rt)         :: dtdx, dtdy, dtdz
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

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             courx = ( qaux(i,j,k,QC) + abs(q(i,j,k,QU)) ) * dtdx
             coury = ( qaux(i,j,k,QC) + abs(q(i,j,k,QV)) ) * dtdy
             courz = ( qaux(i,j,k,QC) + abs(q(i,j,k,QW)) ) * dtdz

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
               print *,'>>> ... u,v,w, c            ', q(i,j,k,QU), q(i,j,k,QV), q(i,j,k,QW), qaux(i,j,k,QC)
               print *,'>>> ... density             ', q(i,j,k,QRHO)
            endif

            courno = max(courno, courtmp)
          enddo
       enddo
    enddo

  end subroutine compute_cfl

  subroutine grctoprim(lo, hi, &
                     uin, uin_lo, uin_hi, &
                     q,     q_lo,   q_hi, &
                     qaux, qa_lo,  qa_hi, &
                     gamma_up, g_lo, g_hi, &
                     alpha, a_lo, a_hi)

    use mempool_module, only : bl_allocate, bl_deallocate
    use actual_network, only : nspec, naux
    use eos_module, only : eos
    use eos_type_module, only : eos_t, eos_input_re
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, &
                                   UEDEN, UTEMP, &
                                   QRHO, QU, QV, QW, &
                                   QREINT, QPRES, QTEMP, QGAME, QFS, QFX, &
                                   NQ, QC, QCSML, QGAMC, QDPDR, QDPDE, NQAUX, &
                                   npassive, upass_map, qpass_map, dual_energy_eta1, &
                                   small_dens
    use bl_constants_module, only: ZERO, HALF, ONE
    use castro_util_module, only: position

    use amrex_fort_module, only : rt => amrex_real
    use riemann_util_module, only : zbrent, f_of_p
    use metric_module, only : calculate_norm

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: uin_lo(3), uin_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: g_lo(3), g_hi(3)
    integer, intent(in) :: a_lo(3), a_hi(3)

    real(rt)        , intent(in ) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)

    real(rt)        , intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt)        , intent(inout) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt)        , intent(in   ) :: gamma_up(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3),9)
    real(rt)        , intent(in   ) :: alpha(a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3))

    real(rt)        , parameter :: small = 1.d-8

    integer          :: i, j, k, g
    integer          :: n, iq, ipassive
    real(rt)         :: vel(3)
    real(rt)         :: gamma_down(9)
    real(rt)         :: pmin, pmax, ssq, p, fmin, fmax, sq, h, W2, eden, gamma

    type (eos_t) :: eos_state

    call eos(eos_input_re, eos_state)
    gamma = eos_state % gam1

    !write(*,*) "rho = ", uin(:,:,:,URHO)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (uin(i,j,k,URHO) .le. ZERO) then
                print *,'   '
                print *,'>>> Error: advection_util_nd.F90::grctoprim ',i, j, k
                print *,'>>> ... negative density ', uin(i,j,k,URHO)
                call bl_error("Error:: advection_util_nd.f90 :: grctoprim")
             else if (uin(i,j,k,URHO) /= uin(i,j,k,URHO)) then
                 print *,'   '
                 print *,'>>> Error: advection_util_nd.F90::grctoprim ',i, j, k
                 print *,'>>> ... density is nan ', uin(i,j,k,URHO)
                 write(*,*) uin(:,:,:,URHO)
                 call bl_error("Error:: advection_util_nd.f90 :: grctoprim")
             else if (uin(i,j,k,URHO) .lt. small_dens) then
                print *,'   '
                print *,'>>> Error: advection_util_nd.F90::grctoprim ',i, j, k
                print *,'>>> ... small density ', uin(i,j,k,URHO)
                call bl_error("Error:: advection_util_nd.f90 :: grctoprim")
             endif
          end do

          do i = lo(1), hi(1)

              ! calculate gamma_down. Assume for now that it's diagonal.

              gamma_down(:) = 0.0d0
              gamma_down(1) = 1.0d0 / gamma_up(i,j,k,1)
              gamma_down(5) = 1.0d0 / gamma_up(i,j,k,5)
              gamma_down(9) = 1.0d0 / gamma_up(i,j,k,9)

              eden = uin(i,j,k,UEDEN)

              if (eden < ZERO) then
                  eden = abs(eden)
              endif

              call calculate_norm(uin(i,j,k,UMX:UMZ), gamma_up(i,j,k,:), ssq)

              pmin = (gamma - 1.0d0) * (eden + uin(i,j,k,URHO) - ssq / (eden + uin(i,j,k,URHO))**2 - uin(i,j,k,URHO))

              pmax = (gamma - 1.0d0) * (eden + uin(i,j,k,URHO) - ssq / (eden + uin(i,j,k,URHO)))

              !write(*,*) "pressure = ", pmin, pmax

              if (pmin < 0.0d0) then
                  pmin = 0.d0
              end if

              !if (pmax > 1.d0 .or. pmax < pmin) then
            !      pmax = 1.0d0
             ! end

              call f_of_p(fmin, pmin, uin(i,j,k,:), gamma_up(i,j,k,:))
              call f_of_p(fmax, pmax, uin(i,j,k,:), gamma_up(i,j,k,:))

              if (fmin * fmax > 0.0d0) then
                  pmin = pmin * 0.1d0 !0.d0
              end if

              call f_of_p(fmin, pmin, uin(i,j,k,:), gamma_up(i,j,k,:))

              !write(*,*) "f = ", fmin, fmax, " p = ", pmin, pmax

              if (fmin * fmax > 0.0d0) then
                  pmax = pmax * 10.d0
              end if

              call zbrent(p, pmin, pmax, uin(i,j,k,:), gamma_up(i,j,k,:))

              !write(*,*) "pressure = ", pmin, pmax

              if (p /= p .or. p <= 0.0d0) then! .or. p > 1.0d0) then

                  p = abs((gamma - 1.0d0) * ((eden + uin(i,j,k,URHO)) - ssq / (eden + uin(i,j,k,URHO))**2 - uin(i,j,k,URHO)))

                  !if (p > 1.0d0) then
                      !p = 1.0d0
                  !end if
              end if

              sq = sqrt((eden + p + uin(i,j,k,URHO))**2 - ssq)

              if (sq /= sq) then
                  sq = eden + p + uin(i,j,k,URHO)
              end if

              h = 1.0d0 + gamma * &
                (sq - p * (eden + p + uin(i,j,k,URHO)) / sq - uin(i,j,k,URHO)) / uin(i,j,k,URHO)
              W2 = 1.0d0 + ssq / (uin(i,j,k,URHO) * h)**2

              !write(*,*) "p, sq, eden, rho", p, sq, eden, uin(i,j,k,URHO)
              !return

              q(i,j,k,QRHO) = uin(i,j,k,URHO) * sq / (eden + &
                p + uin(i,j,k,URHO))

              vel(1) = (gamma_up(i,j,k,1) * uin(i,j,k,UMX) + &
                  gamma_up(i,j,k,2) * uin(i,j,k,UMY) + &
                  gamma_up(i,j,k,3) * uin(i,j,k,UMZ)) /&
                  (W2 * h * q(i,j,k,QRHO))
              vel(2) = (gamma_up(i,j,k,2) * uin(i,j,k,UMX) + &
                  gamma_up(i,j,k,5) * uin(i,j,k,UMY) + &
                    gamma_up(i,j,k,6) * uin(i,j,k,UMZ)) /&
                  (W2 * h * q(i,j,k,QRHO))
              vel(3) = (gamma_up(i,j,k,3) * uin(i,j,k,UMX) + &
                  gamma_up(i,j,k,6) * uin(i,j,k,UMY) + &
                    gamma_up(i,j,k,9) * uin(i,j,k,UMZ)) /&
                  (W2 * h * q(i,j,k,QRHO))

              q(i,j,k,QU) = gamma_down(1) * vel(1) + &
                gamma_down(2) * vel(2) + gamma_down(3) * vel(3)
              q(i,j,k,QV) = gamma_down(2) * vel(1) + &
                gamma_down(5) * vel(2) + gamma_down(6) * vel(3)
              q(i,j,k,QW) = gamma_down(3) * vel(1) + &
                gamma_down(6) * vel(2) + gamma_down(9) * vel(3)

              q(i,j,k,QREINT) = p / (gamma - 1.0d0)

              q(i,j,k,QTEMP) = uin(i,j,k,UTEMP)

              q(i,j,k,QPRES)  = p

             !write(*,*) "pressure = ", p

          enddo
       enddo
    enddo
    !stop

    !write(*,*) "rho = ", q(:,:,:,QRHO)

    ! Load passively advected quatities into q
    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       iq = qpass_map(ipassive)
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                q(i,j,k,iq) = uin(i,j,k,n)/q(i,j,k,QRHO)
             enddo
          enddo
       enddo
    enddo

    ! get gamc, p, T, c, csml using q state
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             eos_state % T   = q(i,j,k,QTEMP )
             eos_state % rho = q(i,j,k,QRHO  )
             eos_state % e   = q(i,j,k,QREINT) / q(i,j,k,QRHO  )
             eos_state % xn  = q(i,j,k,QFS:QFS+nspec-1)
             eos_state % aux = q(i,j,k,QFX:QFX+naux-1)
             eos_state % p   = q(i,j,k,QPRES)

             call eos(eos_input_re, eos_state)

             q(i,j,k,QTEMP)  = eos_state % T
             q(i,j,k,QREINT) = eos_state % e * q(i,j,k,QRHO)

             !write(*,*) "p = ", q(i,j,k,QPRES), eos_state % p, eos_state % e, q(i,j,k,QREINT) / q(i,j,k,QRHO  )

             !q(i,j,k,QPRES)  = eos_state % p
             q(i,j,k,QGAME)  = q(i,j,k,QPRES) / q(i,j,k,QREINT) + ONE

             qaux(i,j,k,QDPDR)  = eos_state % dpdr_e
             qaux(i,j,k,QDPDE)  = eos_state % dpde
             qaux(i,j,k,QGAMC)  = eos_state % gam1
             qaux(i,j,k,QC   )  = eos_state % cs

             qaux(i,j,k,QCSML)  = max(small, small * qaux(i,j,k,QC))
          enddo
       enddo
    enddo

  end subroutine grctoprim

  subroutine normalize_species_fluxes(flux1, flux1_lo, flux1_hi, &
#if (BL_SPACEDIM >= 2)
                                      flux2, flux2_lo, flux2_hi, &
#endif
#if (BL_SPACEDIM == 3)
                                      flux3, flux3_lo, flux3_hi, &
#endif
                                      lo, hi)

    use network, only : nspec
    use meth_params_module, only : NVAR, URHO, UFS
    use bl_constants_module
    use prob_params_module, only : dg
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: lo(3),hi(3)
    integer          :: flux1_lo(3), flux1_hi(3)
    real(rt)         :: flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),flux1_lo(3):flux1_hi(3),NVAR)
#if (BL_SPACEDIM >= 2)
    integer, intent(in) :: flux2_lo(3), flux2_hi(3)
    real(rt), intent(inout) :: flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2),flux2_lo(3):flux2_hi(3),NVAR)
#endif
#if (BL_SPACEDIM == 3)
    integer, intent(in) :: flux3_lo(3), flux3_hi(3)
    real(rt), intent(inout) :: flux3(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2),flux3_lo(3):flux3_hi(3),NVAR)
#endif

    ! Local variables
    integer          :: i,j,k,n
    real(rt)         :: sum,fac

    do k = lo(3),hi(3)
        do j = lo(2),hi(2)
           do i = lo(1),hi(1)+1
              sum = ZERO
              do n = UFS, UFS+nspec-1
                 sum = sum + flux1(i,j,k,n)
              end do
              if (sum .ne. ZERO) then
                 fac = flux1(i,j,k,URHO) / sum
              else
                 fac = ONE
              end if
              do n = UFS, UFS+nspec-1
                 flux1(i,j,k,n) = flux1(i,j,k,n) * fac
              end do
           end do
        end do
#if (BL_SPACEDIM >= 2)
        do j = lo(2),hi(2)+dg(2)
           do i = lo(1),hi(1)
              sum = ZERO
              do n = UFS, UFS+nspec-1
                 sum = sum + flux2(i,j,k,n)
              end do
              if (sum .ne. ZERO) then
                 fac = flux2(i,j,k,URHO) / sum
              else
                 fac = ONE
              end if
              do n = UFS, UFS+nspec-1
                 flux2(i,j,k,n) = flux2(i,j,k,n) * fac
              end do
           end do
        end do
#endif
    end do

#if (BL_SPACEDIM == 3)
    do k = lo(3),hi(3)+dg(3)
        do j = lo(2),hi(2)
           do i = lo(1),hi(1)
               sum = ZERO
               do n = UFS, UFS+nspec-1
                  sum = sum + flux3(i,j,k,n)
               end do
               if (sum .ne. ZERO) then
                  fac = flux3(i,j,k,URHO) / sum
               else
                  fac = ONE
               end if
               do n = UFS, UFS+nspec-1
                  flux3(i,j,k,n) = flux3(i,j,k,n) * fac
               end do
            end do
         end do
     end do
#endif

  end subroutine normalize_species_fluxes

end module advection_util_module
