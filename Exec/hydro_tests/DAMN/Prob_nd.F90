subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

   use probdata_module
   use prob_params_module, only : center
   use bl_constants_module
   use fundamental_constants_module
   use meth_params_module, only: small_temp, small_pres, small_dens
   use eos_type_module, only : eos_t, eos_input_rt, eos_input_rp
   use eos_module, only : eos

   use amrex_fort_module, only : rt => amrex_real
   implicit none

   integer :: init, namlen
   integer :: name(namlen)
   real(rt)         :: problo(3), probhi(3)

   integer :: untin
   integer :: i

   namelist /fortin/ h_in, h_out, damn_rad, g, swe_to_comp_level, p_ambient, dens_ambient, exp_energy, &
        nsub, temp_ambient, dens_incompressible

   integer, parameter :: maxlen=127
   character :: probin*(maxlen)
   character :: model*(maxlen)
   integer :: ipp, ierr, ipp1

   ! Temporary storage variables in case we need to switch the primary and secondary.

   integer :: ioproc
   type(eos_t) :: eos_state

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

   h_in = 2.e0_rt
   h_out = 1.e0_rt
   damn_rad = 2.e-1_rt

   g = 1.0e0_rt
   p_ambient = 1.e-5_rt        ! ambient pressure (in erg/cc)
   dens_ambient = 1.e0_rt      ! ambient density (in g/cc)
   exp_energy = 1.e0_rt        ! absolute energy of the explosion (in erg)
   nsub = 4
   temp_ambient = -1.e2_rt
   dens_incompressible = 1.0e0_rt

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

   ! override the pressure iwth the temperature
   if (temp_ambient > ZERO) then

      eos_state % rho = dens_incompressible
      eos_state % xn(:) = xn_zone(:)
      eos_state % T = temp_ambient

      call eos(eos_input_rt, eos_state)

      p_ambient = eos_state % p

   endif

   ! Calculate ambient state data

   eos_state % rho = dens_incompressible
   eos_state % p   = p_ambient
   eos_state % T   = 1.e5_rt ! Initial guess for iterations
   eos_state % xn  = xn_zone

   call eos(eos_input_rp, eos_state)

   e_ambient = eos_state % e

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
       QFS, QFA, UEDEN, UEINT, QTEMP, NQ, QRHO, QU, QV, QW, QREINT, QPRES, UTEMP, small_dens, small_temp
  use network, only : nspec
  use bl_constants_module
  use prob_params_module, only: problo, center, probhi
  use eos_module, only : eos, initialized, eos_init
  use eos_type_module, only : eos_t, eos_input_rp, eos_input_re, eos_input_rt
  use riemann_util_module, only: comp_cons_state, swe_cons_state

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer :: level, nscal
  integer :: lo(3), hi(3)
  integer :: state_lo(3), state_hi(3)
  real(rt)         :: xlo(3), xhi(3), time, delta(3)
  real(rt)         :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)

  real(rt)         :: xx, yy, zz, r, h

  integer :: i, j, k, n, ii, jj, kk

  real(rt)         :: dye, eint, a
  real(rt)         :: q(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NQ)
  type(eos_t) :: eos_state

  if (.not. initialized) call eos_init(small_dens=small_dens, small_temp=small_temp)

  dye = ZERO

  if (level <= swe_to_comp_level) then
      write(*,*) "Initialising level ", level, " with SWE data"
  else
      write(*,*) "Initialising level ", level, " with incompressible data"
  endif

  a = 0.0005e0_rt ! characteristic size of layer between states

  !$OMP PARALLEL DO PRIVATE(i, j, k, xx, yy, zz, r, h)
  do k = lo(3), hi(3)
     zz = xlo(3) + delta(3)*dble(k-lo(3)+HALF)

     do j = lo(2), hi(2)
        yy = xlo(2) + delta(2)*dble(j-lo(2)+HALF)


        r = yy

        ! circular dam
        !r = sqrt((yy - center(2))**2 + (zz - center(3))**2)
        
        h = h_in + 0.5e0_rt * (h_out - h_in) * (1.0e0_rt + tanh((r - damn_rad) / a))

        do i = lo(1), hi(1)

            state(i,j,k,1:NVAR) = 0.e0_rt
            q(i,j,k,1:NQ) = 0.e0_rt

            xx = xlo(1) + delta(1)*dble(i-lo(1)+HALF)

            q(i,j,k,QRHO) = dens_incompressible

            eos_state % p = 0.5e0_rt * g * (h-xx)**2

            ! if (r < damn_rad) then
            !     eos_state % p = 0.5e0_rt * g * (h_in-xx)**2
            ! else
            !     eos_state % p = 0.5e0_rt * g * (h_out-xx)**2
            ! end if

            !eos_state % e = e_zone
            eos_state % rho = q(i,j,k, QRHO)
            eos_state % xn(:) = xn_zone(:)
            eos_state%T = 100000.e0_rt ! initial guess
            !eos_state % p = 0.5e0_rt * g * state(i,j,k,URHO)**2

            call eos(eos_input_rp, eos_state)

            eint = q(i,j,k,QRHO) * eos_state % e
            q(i,j,k,QU:QW) = 0.e0_rt

            q(i,j,k,QREINT) = eint +  &
                0.5e0_rt*sum(q(i,j,k,QU:QW)**2)*q(i,j,k, QRHO)

            if (level <= swe_to_comp_level) then
                ! shallow water level
                q(i,j,k,QRHO) = h
            end if

            q(i,j,k,QFS) = q(i,j,k,QRHO)

            q(i,j,k,QTEMP) = eos_state % T

            q(i,j,k,QFA)  = dye
            q(i,j,k,QFS:QFS-1+nspec) = q(i,j,k,QRHO) / nspec

            if (level <= swe_to_comp_level) then
                ! shallow water level
                call swe_cons_state(q(i,j,k,:), state(i,j,k,:))
            else
                call comp_cons_state(q(i,j,k,:), state(i,j,k,:))
                ! state(i,j,k,URHO) = dens_incompressible
            end if

        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

end subroutine ca_initdata
