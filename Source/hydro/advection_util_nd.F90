module advection_util_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  private

  public enforce_minimum_density, compute_cfl, swectoprim, compctoprim

contains

  subroutine enforce_minimum_density(uin,uin_lo,uin_hi, &
                                     uout,uout_lo,uout_hi, &
                                     vol,vol_lo,vol_hi, &
                                     lo,hi,frac_change,verbose)

    use network, only : nspec, naux
    use meth_params_module, only : NVAR, QRHO, QREINT, UEDEN, small_dens, density_reset_method, NQ, NQAUX
    use bl_constants_module, only : ZERO
    use riemann_util_module, only : swe_cons_state

    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: lo(3), hi(3), verbose
    integer, intent(in) ::  uin_lo(3),  uin_hi(3)
    integer, intent(in) :: uout_lo(3), uout_hi(3)
    integer, intent(in) ::  vol_lo(3),  vol_hi(3)

    real(rt)        , intent(inout) ::  uin( uin_lo(1): uin_hi(1), uin_lo(2): uin_hi(2), uin_lo(3): uin_hi(3),NVAR)
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

    return

    max_dens = ZERO

    have_reset = .false.

    call swectoprim(lo, hi, &
                      uin, uin_lo, uin_hi, &
                      qin,  uin_lo, uin_hi, &
                      qaux, lo, hi)
    call swectoprim(lo, hi, &
                    uout, uout_lo, uout_hi, &
                    qout,  uout_lo, uout_hi, &
                    qaux,  lo, hi)

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

                call swe_cons_state(qout(i,j,k,:), uout(i,j,k,:))

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



  subroutine compute_cfl(q, q_lo, q_hi, &
                         qaux, qa_lo, qa_hi, &
                         lo, hi, dt, dx, courno) &
                         bind(C, name = "compute_cfl")

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

    large_number = 1.0d30

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

              if (q(i,j,k,QRHO) .le. ZERO) then
                 !print *, uin(:,:,:,URHO)
                 print *,'   '
                 print *,'>>> Error: advection_util_nd.F90::compute_cfl ',i, j, k
                 print *,'>>> ... negative density ', q(i,j,k,QRHO)
                 q(i,j,k,QRHO) = 1.0d0
                 !call bl_error("Error:: advection_util_nd.f90 :: compute_cfl")
             else if (q(i,j,k,QRHO) > large_number) then
                    !print *, uin(:,:,:,URHO)
                    print *,'   '
                    print *,'>>> Error: advection_util_nd.F90::compute_cfl ',i, j, k
                    print *,'>>> ... density is very large ', q(i,j,k,QRHO)
                    q(i,j,k,QRHO) = 1.0d0
                    !call bl_error("Error:: advection_util_nd.f90 :: compute_cfl")
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
            endif

            courno = max(min(1.0d0,courno), min(1.0d0,courtmp))

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

  end subroutine compute_cfl

  subroutine swectoprim(lo, hi, &
                     uin, uin_lo, uin_hi, &
                     q,     q_lo,   q_hi, &
                     qaux, qa_lo,  qa_hi, ignore_errors)

    use mempool_module, only : bl_allocate, bl_deallocate
    use actual_network, only : nspec, naux
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, &
                                   QRHO, QU, QV, QW, QTEMP, UTEMP, &
                                   NQ, QC, QCSML, QGAMC, QDPDR, QDPDE, NQAUX, QPRES, QREINT, UEINT, &
                                   npassive, upass_map, qpass_map, dual_energy_eta1, &
                                   small_dens, QFA
    use bl_constants_module, only: ZERO, HALF, ONE
    use castro_util_module, only: position
    use probdata_module, only : g

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

    real(rt)        , parameter :: small = 1.d-8

    integer          :: i, j, k, ii, jj, kk
    integer          :: n, iq, ipassive

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
                    ! uin(i,j,k,:) = 0.0d0
                    ! uin(i,j,k,URHO) = 1.0d0
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
              q(i,j,k,1:NQ) = 0.0d0

              !q(i,j,k,:QW) = uin(i,j,k,:QW)
              if (uin(i,j,k,URHO) .le. ZERO) then
                  q(i,j,k,QRHO) = 1.0d0
                  q(i,j,k,QU) = 0.1d0 * uin(i,j,k,UMX) !/ uin(i,j,k,URHO)
                  q(i,j,k,QV) = 0.1d0 * uin(i,j,k,UMY) !/ uin(i,j,k,URHO)
                  q(i,j,k,QW) = 0.1d0 * uin(i,j,k,UMZ) !/ uin(i,j,k,URHO)
              else
                  q(i,j,k,QRHO) = uin(i,j,k,URHO)
                  q(i,j,k,QU) = uin(i,j,k,UMX) / uin(i,j,k,URHO)
                  q(i,j,k,QV) = uin(i,j,k,UMY) / uin(i,j,k,URHO)
                  q(i,j,k,QW) = uin(i,j,k,UMZ) / uin(i,j,k,URHO)
              endif
              q(i,j,k,QPRES) = 0.5d0 * g * uin(i,j,k,URHO)**2

              q(i,j,k,QREINT) = uin(i,j,k,UEINT)
              q(i,j,k,QTEMP) = uin(i,j,k,UTEMP)

              q(i,j,k,QFA) = 0.0d0

          enddo
       enddo
    enddo

    !Load passively advected quatities into q
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

end subroutine swectoprim

subroutine compctoprim(lo, hi, &
                   uin, uin_lo, uin_hi, &
                   q,     q_lo,   q_hi, &
                   qaux, qa_lo,  qa_hi, ignore_errors)

  use mempool_module, only : bl_allocate, bl_deallocate
  use actual_network, only : nspec, naux
  use eos_module, only : eos
  use eos_type_module, only : eos_t, eos_input_rp
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT,&
                                 QRHO, QU, QV, QW, QREINT, QTEMP, &
                                 NQ, QC, QCSML, QGAMC, QDPDR, QDPDE, NQAUX, QFS, QFX, QGAME, QPRES, UTEMP, &
                                 npassive, upass_map, qpass_map, dual_energy_eta1, &
                                 small_dens, UFA, QFA
  use bl_constants_module, only: ZERO, HALF, ONE
  use castro_util_module, only: position
  use probdata_module, only : g

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: uin_lo(3), uin_hi(3)
  integer, intent(in) :: q_lo(3), q_hi(3)
  integer, intent(in) :: qa_lo(3), qa_hi(3)

  real(rt)        , intent(in ) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)

  real(rt)        , intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
  real(rt)        , intent(inout) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

  logical, optional, intent(in) :: ignore_errors

  real(rt)        , parameter :: small = 1.d-8

  integer          :: i, j, k
  integer          :: n, iq, ipassive
  real(rt)         :: kineng
  type (eos_t)     :: eos_state

  !q(:,:,:,:) = 0.0d0

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)

        if (present(ignore_errors) .and. (.not. ignore_errors)) then
            do i = lo(1), hi(1)
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
            end do
        endif

        do i = lo(1), hi(1)
            q(i,j,k,1:NQ) = 0.0d0

            q(i,j,k,QRHO) = uin(i,j,k,URHO)
            q(i,j,k,QU:QW) = uin(i,j,k,UMX:UMZ) / uin(i,j,k,URHO)

            kineng = HALF * q(i,j,k,QRHO) * (q(i,j,k,QU)**2 + q(i,j,k,QV)**2 + q(i,j,k,QW)**2)

            ! if ( (uin(i,j,k,UEDEN) - kineng) / uin(i,j,k,UEDEN) .gt. dual_energy_eta1) then
            !     q(i,j,k,QREINT) = (uin(i,j,k,UEDEN) - kineng) / uin(i,j,k,URHO)
            !  else
            !     q(i,j,k,QREINT) = uin(i,j,k,UEINT) / uin(i,j,k,URHO)
            ! endif

            q(i,j,k,QREINT) = uin(i,j,k,UEINT)

            q(i,j,k,QTEMP) = uin(i,j,k,UTEMP)

            q(i,j,k,QFA) = 0.0d0

        enddo
     enddo
  enddo

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
             eos_state % e   = q(i,j,k,QREINT)
             eos_state % xn  = q(i,j,k,QFS:QFS+nspec-1)
             eos_state % aux = q(i,j,k,QFX:QFX+naux-1)
             eos_state % p = 0.5d0 * g * q(i,j,k,QRHO  )**2
             q(i,j,k,QPRES)  = eos_state % p

             call eos(eos_input_rp, eos_state)

             q(i,j,k,QTEMP)  = eos_state % T
             q(i,j,k,QREINT) = eos_state % e * q(i,j,k,QRHO)
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

end subroutine compctoprim


! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine calc_pdivu(lo, hi, &
                        q1, q1_lo, q1_hi, &
                        area1, a1_lo, a1_hi, &
#if BL_SPACEDIM >= 2
                        q2, q2_lo, q2_hi, &
                        area2, a2_lo, a2_hi, &
#endif
#if BL_SPACEDIM == 3
                        q3, q3_lo, q3_hi, &
                        area3, a3_lo, a3_hi, &
#endif
                        vol, v_lo, v_hi, &
                        dx, pdivu, div_lo, div_hi)

    ! this computes the *node-centered* divergence

    use meth_params_module, only : QU, QV, QW, NQ, GDPRES, GDU, GDV, GDW
    use bl_constants_module
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) :: lo(3), hi(3)

    integer, intent(in) :: div_lo(3), div_hi(3)
    real(rt), intent(in) :: dx(3)
    real(rt), intent(inout) :: pdivu(div_lo(1):div_hi(1),div_lo(2):div_hi(2),div_lo(3):div_hi(3))

    integer, intent(in) :: q1_lo(3), q1_hi(3)
    integer, intent(in) :: a1_lo(3), a1_hi(3)
    real(rt), intent(in) :: q1(q1_lo(1):q1_hi(1),q1_lo(2):q1_hi(2),q1_lo(3):q1_hi(3),NQ)
    real(rt), intent(in) :: area1(a1_lo(1):a1_hi(1),a1_lo(2):a1_hi(2),a1_lo(3):a1_hi(3))
#if BL_SPACEDIM >= 2
    integer, intent(in) :: q2_lo(3), q2_hi(3)
    integer, intent(in) :: a2_lo(3), a2_hi(3)
    real(rt), intent(in) :: q2(q2_lo(1):q2_hi(1),q2_lo(2):q2_hi(2),q2_lo(3):q2_hi(3),NQ)
    real(rt), intent(in) :: area2(a2_lo(1):a2_hi(1),a1_lo(2):a1_hi(2),a1_lo(3):a1_hi(3))
#endif
#if BL_SPACEDIM == 3
    integer, intent(in) :: q3_lo(3), q3_hi(3)
    integer, intent(in) :: a3_lo(3), a3_hi(3)
    real(rt), intent(in) :: q3(q3_lo(1):q3_hi(1),q3_lo(2):q3_hi(2),q3_lo(3):q3_hi(3),NQ)
    real(rt), intent(in) :: area3(a3_lo(1):a3_hi(1),a1_lo(2):a1_hi(2),a1_lo(3):a1_hi(3))
#endif
    integer, intent(in) :: v_lo(3), v_hi(3)
    real(rt), intent(in) :: vol(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))

    integer  :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

#if BL_SPACEDIM == 1
             pdivu(i,j,k) = HALF * &
                  (q1(i+1,j,k,GDPRES) + q1(i,j,k,GDPRES))* &
                  (q1(i+1,j,k,GDU)*area1(i+1,j,k) - q1(i,j,k,GDU)*area1(i,j,k)) / vol(i,j,k)
#endif

#if BL_SPACEDIM == 2
             pdivu(i,j,k) = HALF*( &
                  (q1(i+1,j,k,GDPRES) + q1(i,j,k,GDPRES)) * &
                  (q1(i+1,j,k,GDU)*area1(i+1,j,k) - q1(i,j,k,GDU)*area1(i,j,k)) + &
                  (q2(i,j+1,k,GDPRES) + q2(i,j,k,GDPRES)) * &
                  (q2(i,j+1,k,GDV)*area2(i,j+1,k) - q2(i,j,k,GDV)*area2(i,j,k)) ) / vol(i,j,k)
#endif

#if BL_SPACEDIM == 3
             pdivu(i,j,k) = &
                  HALF*(q1(i+1,j,k,GDPRES) + q1(i,j,k,GDPRES)) * &
                       (q1(i+1,j,k,GDU) - q1(i,j,k,GDU))/dx(1) + &
                  HALF*(q2(i,j+1,k,GDPRES) + q2(i,j,k,GDPRES)) * &
                       (q2(i,j+1,k,GDV) - q2(i,j,k,GDV))/dx(2) + &
                  HALF*(q3(i,j,k+1,GDPRES) + q3(i,j,k,GDPRES)) * &
                       (q3(i,j,k+1,GDW) - q3(i,j,k,GDW))/dx(3)
#endif

          enddo
       enddo
    enddo

  end subroutine calc_pdivu

end module advection_util_module
