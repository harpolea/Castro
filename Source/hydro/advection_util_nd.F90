module advection_util_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  private

  public enforce_minimum_density, compute_cfl, swectoprim

contains

  subroutine enforce_minimum_density(uin,uin_lo,uin_hi, &
                                     uout,uout_lo,uout_hi, &
                                     vol,vol_lo,vol_hi, &
                                     lo,hi,frac_change,verbose)

    use network, only : nspec, naux
    use meth_params_module, only : NVAR, QRHO, QREINT, UEDEN, small_dens, density_reset_method, NQ, NQAUX
    use bl_constants_module, only : ZERO
    use riemann_util_module, only : grswe_cons_state

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

    return

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

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer :: lo(3), hi(3)
    integer :: q_lo(3), q_hi(3), qa_lo(3), qa_hi(3)

    real(rt)         :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt)         :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt)         :: dt, dx(3), courno

    real(rt)         :: courx, coury, courz, courmx, courmy, courmz, courtmp
    real(rt)         :: dtdx, dtdy, dtdz, c
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

              c = sqrt(q(i,j,k,QRHO))

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
            endif

            courno = max(courno, courtmp)
          enddo
       enddo
    enddo

  end subroutine compute_cfl

  subroutine swectoprim(lo, hi, &
                     uin, uin_lo, uin_hi, &
                     q,     q_lo,   q_hi, &
                     qaux, qa_lo,  qa_hi, &
                     gamma_up, glo, ghi)

    use mempool_module, only : bl_allocate, bl_deallocate
    use actual_network, only : nspec, naux
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, &
                                   QRHO, QU, QV, QW, &
                                   NQ, NQAUX, &
                                   small_dens
    use bl_constants_module, only: ZERO, HALF, ONE
    use castro_util_module, only: position
    use metric_module, only: calculate_norm

    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: uin_lo(3), uin_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: glo(3), ghi(3)

    real(rt)        , intent(in ) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
    real(rt)        , intent(in ) :: gamma_up(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),9)

    real(rt)        , intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt)        , intent(inout) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt)        , parameter :: small = 1.d-20

    integer          :: i, j, k
    real(rt)         :: W, ss

    q(:,:,:,:) = 0.0d0

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)

          do i = lo(1), hi(1)
             if (uin(i,j,k,URHO) .le. ZERO) then
                print *,'   '
                print *,'>>> Error: advection_util_nd.F90::swectoprim ',i, j, k
                print *,'>>> ... negative density ', uin(i,j,k,URHO)
                call bl_error("Error:: advection_util_nd.f90 :: swectoprim")
             else if (uin(i,j,k,URHO) /= uin(i,j,k,URHO)) then
                 print *,'   '
                 print *,'>>> Error: advection_util_nd.F90::swectoprim ',i, j, k
                 print *,'>>> ... density is nan ', uin(i,j,k,URHO)
                 write(*,*) uin(:,:,:,URHO)
                 call bl_error("Error:: advection_util_nd.f90 :: swectoprim")
             else if (uin(i,j,k,URHO) .lt. small_dens) then
                print *,'   '
                print *,'>>> Error: advection_util_nd.F90::swectoprim ',i, j, k
                print *,'>>> ... small density ', uin(i,j,k,URHO)
                call bl_error("Error:: advection_util_nd.f90 :: swectoprim")
             endif
          end do

          do i = lo(1), hi(1)

              call calculate_norm(uin(i,j,k,UMX:UMZ), gamma_up(i,j,k,:), ss)

              ! W^2 = 1 + S_jS^j / D^2
              W = sqrt(1.0d0 + ss / uin(i,j,k,URHO)**2)

              !if (abs(W - 1.0d0) > 0.0001d0) then
                 ! write(*,*) "W = ", W
              !end if

              ! HACK
              W = 1.0d0

              ! initialise
              q(i,j,k,:) = 0.0d0!uin(i,j,k,:QW)

              q(i,j,k,QRHO) = uin(i,j,k,URHO) / W
              q(i,j,k,QU) = uin(i,j,k,UMX) / (uin(i,j,k,URHO) * W)
              q(i,j,k,QV) = uin(i,j,k,UMY) / (uin(i,j,k,URHO) * W)
              q(i,j,k,QW) = uin(i,j,k,UMZ) / (uin(i,j,k,URHO) * W)

              if (uin(i,j,k,URHO) < 0.0d0 .or. q(i,j,k,QRHO) < 0.0d0 ) then
                  write(*,*) "D, rho", uin(i,j,k,URHO), q(i,j,k,QRHO)
                  stop
              end if

              !q(i,j,k,QW+1:) = q(i,j,k,QW+1:) / W

          enddo
       enddo
    enddo

  end subroutine swectoprim

end module advection_util_module
