subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

   use probdata_module
   use prob_params_module, only : center
   use bl_constants_module
   use fundamental_constants_module

   use amrex_fort_module, only : rt => amrex_real
   implicit none

   integer :: init, namlen
   integer :: name(namlen)
   real(rt)         :: problo(3), probhi(3)

   integer :: untin
   integer :: i

   namelist /fortin/ g, swe_to_comp_level, &
        nsub, eos_K, t_r

   integer, parameter :: maxlen = 127
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

   t_r = 1.0_rt
   g = 1.0_rt
   nsub = 4
   eos_K = 1._rt

   swe_to_comp_level = 0

   ! set explosion center
   center(1:3) = HALF*(problo(1:3) + probhi(1:3))

   ! Read namelists -- override the defaults

   untin = 9
   open(untin,file=probin(1:namlen),form='formatted',status='old')
   read(untin,fortin)
   close(unit=untin)

   xn_zone(:) = ZERO
   xn_zone(1) = ONE

   ! characteristic scales
   x_r = probhi(2) - problo(2)
   q_r = 0.4_rt*M_pi*x_r/t_r

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
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, &
       QFS, QFA, UEDEN, UEINT, QTEMP, NQ, QRHO, QU, QV, QW, QREINT, QPRES, UTEMP, small_dens, small_temp, nadv, UFA
  use network, only : nspec
  use bl_constants_module
  use prob_params_module, only: problo, center, probhi
  use eos_module, only : eos, initialized, eos_init
  use eos_type_module, only : eos_t, eos_input_rp, eos_input_re, eos_input_rt
  use riemann_util_module, only: comp_cons_state, swe_cons_state
  use actual_eos_module, only: gamma_const

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: level, nscal
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: state_lo(3), state_hi(3)
  real(rt), intent(in)         :: xlo(3), xhi(3), time, delta(3)
  real(rt), intent(inout)         :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)

  real(rt)         :: xx, yy, zz, r, h
  real(rt)         :: yl, zl, yc, zc, y, z

  integer :: i, j, k, jj, kk

  real(rt)         :: reint, p, u_phi, u_tot, p0, rho0
  real(rt)         :: q(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NQ)
  type(eos_t) :: eos_state

  if (.not. initialized) call eos_init(small_dens=small_dens, small_temp=small_temp)

  if (level <= swe_to_comp_level) then
      write(*,*) "Initialising level ", level, " with SWE data"
  else
      write(*,*) "Initialising level ", level, " with incompressible data"
  endif

  eos_state % rho = 1.0_rt
  eos_state % e = 1.0_rt

  if (.not. initialized) call eos_init(small_dens=small_dens, small_temp=small_temp)

  call eos(eos_input_re, eos_state)

  ! M = 0.1, uphi = 0.25
  p0 = ((gamma_const-1._rt)/gamma_const * g * eos_K * 5._rt)**(gamma_const / (gamma_const - 1._rt))
  rho0 = eos_K * p0**(1._rt / gamma_const)

  !$OMP PARALLEL DO PRIVATE(i, j, k, xx, yy, zz, r, h)
  do k = lo(3), hi(3)
      zl = problo(3) + delta(3)*dble(k)
      z = problo(3) + delta(3)*dble(k+0.5_rt)

     do j = lo(2), hi(2)
        yl = problo(2) + delta(2)*dble(j)
        y = problo(2) + delta(2)*dble(j+0.5_rt)

        reint = 0.0_rt
        u_tot = 0.0_rt

        do kk = 0, nsub-1
           zz = zl + delta(3)*dble(kk+0.5_rt)/nsub

            do jj = 0, nsub-1
               yy = yl + delta(2)*dble(jj+0.5_rt)/nsub

               r = sqrt((yy - center(2))**2 + (zz - center(3))**2)

               if (r < 0.2_rt) then
                  u_phi = 5.0_rt*r
                  p = p0 + rho0 * 12.5_rt*r**2

               else if (r < 0.4_rt) then
                  u_phi = 2.0_rt - 5.0_rt*r
                  p = p0 + rho0 * (12.5_rt*r**2 + 4.0_rt*(1.0_rt - 5.0_rt*r - log(0.2_rt) + log(r)))

               else
                  u_phi = 0.0_rt
                  p = p0 - rho0 * (2.0_rt - 4.0_rt*log(2.0_rt))
               endif

               u_tot = u_tot + u_phi
               reint = reint + p/(eos_state % gam1 - 1.0_rt)

           enddo
       enddo

       u_phi = u_tot/(nsub*nsub)
       reint = reint/(nsub*nsub)
       p = (eos_state % gam1 - 1.0_rt) * reint

       h = height_from_p(p, 0._rt)

       ! velocity is based on the reference velocity, q_r
       yc = yl + 0.5_rt*delta(2)
       zc = zl + 0.5_rt*delta(3)

       r = sqrt((z - center(3))**2 + (y - center(2))**2)

       q(lo(1):hi(1),j,k,QU) = 0.0_rt
       q(lo(1):hi(1),j,k,QV) = -rho0*q_r*u_phi*((zc-center(3))/r)  ! -sin(phi) = z/r
       q(lo(1):hi(1),j,k,QW) = rho0*q_r*u_phi*((yc-center(2))/r)

       do i = lo(1), hi(1)
           xx = xlo(1) + delta(1)*dble(i-lo(1)+HALF)

           q(i,j,k,QPRES) = p_from_height(h, xx)

           q(i,j,k,QRHO) = rho_from_p(q(i,j,k,QPRES))

           eos_state % rho = q(i,j,k,QRHO)
           eos_state % xn(:) = xn_zone(:)
           eos_state % T = 1.0e5_rt ! initial guess
           eos_state % p = q(i,j,k,QPRES)

           call eos(eos_input_rp, eos_state)

           reint = q(i,j,k,QRHO) * eos_state % e

           q(i,j,k,QREINT) = reint +  &
                0.5_rt*sum(q(i,j,k,QU:QW)**2)*q(i,j,k,QRHO)

           q(i,j,k,QTEMP) = eos_state % T

           if (level <= swe_to_comp_level) then
               ! shallow water level
               q(i,j,k,QRHO) = h
           end if

           q(i,j,k,QFA:QFA-1+nadv)  = 0.0_rt
           q(i,j,k,QFS:QFS-1+nspec) = q(i,j,k,QRHO) / nspec

           if (level <= swe_to_comp_level) then
               ! shallow water level
               call swe_cons_state(q(i,j,k,:), state(i,j,k,:))
           else
               call comp_cons_state(q(i,j,k,:), state(i,j,k,:))
           end if

        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

end subroutine ca_initdata
