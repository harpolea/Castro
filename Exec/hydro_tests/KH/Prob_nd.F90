subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

   use probdata_module
   use bl_constants_module
   use fundamental_constants_module
   use meth_params_module, only: small_temp, small_pres, small_dens

   use amrex_fort_module, only : rt => amrex_real
   implicit none

   integer :: init, namlen
   integer :: name(namlen)
   real(rt)         :: problo(3), probhi(3)

   integer :: untin
   integer :: i

   namelist /fortin/ &
        rho1, rho2, pressure, problem, bulk_velocity

   integer, parameter :: maxlen=127
   character :: probin*(maxlen)
   character :: model*(maxlen)
   integer :: ipp, ierr, ipp1

   ! Temporary storage variables in case we need to switch the primary and secondary.

   integer :: ioproc

   ! For outputting -- determine if we are the IO processor
   call bl_pd_is_ioproc(ioproc)

   ! Build "probin" filename -- the name of file containing fortin namelist.
   if (namlen .gt. maxlen) then
      call bl_error("ERROR: probin file name too long")
   end if

   do i = 1, namlen
      probin(i:i) = char(name(i))
   end do

   ! Set namelist defaults

   problem = 2

   rho1 = 1.0d-3
   rho2 = 2.0d-3
   pressure = 2.5d-3

   bulk_velocity = 0.0

   ! Read namelists -- override the defaults

   untin = 9
   open(untin,file=probin(1:namlen),form='formatted',status='old')
   read(untin,fortin)
   close(unit=untin)

   ! Force a different pressure choice for problem 5

   if (problem .eq. 5) then
      pressure = 10.0d-3
   endif

end subroutine amrex_probinit


! ::: -----------------------------------------------------------
! ::: This routine is called at problem setup time and is used
! ::: to initialize data on each grid.
! :::
! ::: NOTE:  all arrays have one cell of ghost zones surrounding
! :::        the grid interior.  Values in these cells need not
! :::        be set here.
! :::
! ::: INPUTS/OUTPUTS:
! :::
! ::: level     => amr level of grid
! ::: time      => time at which to init data
! ::: lo,hi     => index limits of grid interior (cell centered)
! ::: nstate    => number of state components.  You should know
! :::		   this already!
! ::: state     <=  Scalar array
! ::: delta     => cell size
! ::: xlo,xhi   => physical locations of lower left and upper
! :::              right hand corner of grid.  (does not include
! :::		   ghost region).
! ::: -----------------------------------------------------------
subroutine ca_initdata(level,time,lo,hi,nscal, &
                       state,state_lo,state_hi, &
                       delta,xlo,xhi)

  use probdata_module
  use eos_module, only : eos
  use eos_type_module, only: eos_t, eos_input_rp
  use meth_params_module, only : NVAR, UTEMP, &
       UFS, UFA, NQ, QRHO, QU, QV, QW, QREINT, QPRES, URHO, QTEMP, UMX, UMY, UMZ
  use network, only : nspec
  use bl_constants_module
  use prob_params_module, only: problo, center, probhi
  use riemann_util_module, only: gr_cons_state

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer :: level, nscal
  integer :: lo(3), hi(3)
  integer :: state_lo(3), state_hi(3)
  real(rt)         :: xlo(3), xhi(3), time, delta(3)
  real(rt)         :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)
  real(rt)         :: q(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NQ)

  real(rt)         :: xx, yy, zz

  type (eos_t) :: eos_state

  integer :: i, j, k, n

  real(rt)         :: dens, velx, vely, velz
  real(rt)         :: w0, sigma, ramp, delta_y
  real(rt)         :: vel1, vel2
  real(rt)         :: y1, y2
  real(rt)         :: dye, gamma_up(9), gamma

  real(rt) :: sine_n

  vel1 = -0.2d0
  vel2 =  0.2d0

  gamma_up(:) = 0.0d0
  gamma_up(1) = 1.0d0
  gamma_up(5) = 1.0d0
  gamma_up(9) = 1.0d0

  call eos(eos_input_rp, eos_state)
  gamma = eos_state % gam1

  sine_n = 4.0d0
  w0 = 0.1
  delta_y = 0.025

  y1 = center(2) - (probhi(2) - problo(2)) * 0.25e0_rt
  y2 = center(2) + (probhi(2) - problo(2)) * 0.25e0_rt

  velz = 0.0d0

  !$OMP PARALLEL DO PRIVATE(i, j, k, xx, yy, zz, dens, velx, vely, velz, ramp, eos_state)
  do k = lo(3), hi(3)
     zz = xlo(3) + delta(3)*dble(k-lo(3)+HALF)

     do j = lo(2), hi(2)
        yy = xlo(2) + delta(2)*dble(j-lo(2)+HALF)

        do i = lo(1), hi(1)
           xx = xlo(1) + delta(1)*dble(i-lo(1)+HALF)

           ! Assume the initial y-velocity represents the bulk flow
           ! which will be perturbed in the following step

           vely = bulk_velocity
           dye = ZERO

          if ( yy .lt. y1 ) then
             dens = rho1 - (rho1 - rho2) / 2 * exp( (yy-y1) / delta_y )
             velx = vel1 - (vel1 - vel2) / 2 * exp( (yy-y1) / delta_y )
          else if ( yy .le. HALF * (y1 + y2) ) then
             dens = rho2 + (rho1 - rho2) / 2 * exp( (y1-yy) / delta_y )
             velx = vel2 + (vel1 - vel2) / 2 * exp( (y1-yy) / delta_y )
          else if ( yy .lt. y2 ) then
             dens = rho2 + (rho1 - rho2) / 2 * exp( (yy-y2) / delta_y )
             velx = vel2 + (vel1 - vel2) / 2 * exp( (yy-y2) / delta_y )
          else
             dens = rho1 - (rho1 - rho2) / 2 * exp( (y2-yy) / delta_y )
             velx = vel1 - (vel1 - vel2) / 2 * exp( (y2-yy) / delta_y )
          endif

          vely = vely + w0 * sin(sine_n*M_PI*xx)

           q(i,j,k,QRHO) = dens
           q(i,j,k,QU) = velx
           q(i,j,k,QV) = vely
           q(i,j,k,QW) = velz
           q(i,j,k,QPRES) = pressure
           q(i,j,k,QREINT) = pressure / (gamma - 1.0d0)

           eos_state % rho = dens
           eos_state % p   = pressure
           eos_state % e = q(i,j,k,QREINT) / dens

           state(i,j,k,UFS:UFS-1+nspec) = ONE / nspec
           state(i,j,k,UFA)  = dye

           eos_state % xn  = state(i,j,k,UFS:UFS-1+nspec)

           call eos(eos_input_rp, eos_state)

           q(i,j,k,QTEMP) = eos_state % T

           call gr_cons_state(q(i,j,k,:), state(i,j,k,:), gamma_up)

           ! Establish the thermodynamic quantities
           state(i,j,k,UTEMP) = eos_state % T

           state(i,j,k,UFS:UFS-1+nspec) = state(i,j,k,URHO) / nspec

           !do n = 1,nspec
              !state(i,j,k,UFS+n-1) = state(i,j,k,URHO) * !state(i,j,k,UFS+n-1)
           !end do

           !write(*,*) "q: ", q
           !write(*,*) "state: ", state(i,j,k,:)
           !write(*,*) "QU, QV, QPRES, QW, QREINT, QTEMP", QU, QV, QPRES, QW, QREINT, QTEMP
           !write(*,*) "vely, qv, umy", velx, q(i,j,k,QU), state(i,j,k,UMX)

        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

  !stop

end subroutine ca_initdata
