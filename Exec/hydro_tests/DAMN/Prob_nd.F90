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

   namelist /fortin/ h_in, h_out, damn_rad, g, swe_to_comp_level, dens_ambient

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

   h_in = 2.0d0
   h_out = 1.0d0
   damn_rad = 0.2d0

   g = 1.0d0
   dens_ambient = 1.0d0

   swe_to_comp_level = 0

   ! Read namelists -- override the defaults

   untin = 9
   open(untin,file=probin(1:namlen),form='formatted',status='old')
   read(untin,fortin)
   close(unit=untin)

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
       UFS, UFA
  use network, only : nspec
  use bl_constants_module
  use prob_params_module, only: problo, center, probhi
  use eos_module, only : eos
  use eos_type_module, only : eos_t, eos_input_rp, eos_input_re, eos_input_rt

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer :: level, nscal
  integer :: lo(3), hi(3)
  integer :: state_lo(3), state_hi(3)
  real(rt)         :: xlo(3), xhi(3), time, delta(3)
  real(rt)         :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)

  real(rt)         :: xx, yy, zz, r

  integer :: i, j, k, n

  real(rt)         :: dye
  type(eos_t) :: eos_state

  state(:,:,:,:) = 0.0d0
  dye = ZERO

  !$OMP PARALLEL DO PRIVATE(i, j, k, xx, yy, zz, r)
  do k = lo(3), hi(3)
     zz = xlo(3) + delta(3)*dble(k-lo(3)+HALF)

     do j = lo(2), hi(2)
        yy = xlo(2) + delta(2)*dble(j-lo(2)+HALF)

        do i = lo(1), hi(1)
           xx = xlo(1) + delta(1)*dble(i-lo(1)+HALF)

           r = sqrt((xx - center(1))**2 + (yy - center(2))**2)

           if (swe_to_comp_level <= level) then
               ! shallow water level
               if (r < damn_rad) then
                   state(i,j,k,URHO) = h_in
               else
                   state(i,j,k,URHO) = h_out
               end if

               do n = 1,nspec
                  state(i,j,k,UFS+n-1) = state(i,j,k,URHO) * state(i,j,k,UFS+n-1)
               end do

           else ! compressible level
               e_zone = exp_energy/vctr/dens_ambient

                eos_state % e = e_zone
                eos_state % rho = dens_ambient
                eos_state % xn(:) = xn_zone(:)
                eos_state % T = 1000.00 ! initial guess

                call eos(eos_input_re, eos_state)

                p_zone = (vol_pert*p_exp + vol_ambient*p_ambient)/ (vol_pert + vol_ambient)

               eos_state % p = p_zone
               eos_state % rho = dens_ambient
               eos_state % xn(:) = xn_zone(:)
               eos_state % T = 1000.0   ! initial guess

               call eos(eos_input_rp, eos_state)

               eint = dens_ambient * eos_state % e

               state(i,j,URHO) = dens_ambient
               state(i,j,UMX:UMZ) = 0.e0_rt

               state(i,j,UEDEN) = eint +  &
                    0.5e0_rt*(sum(state(i,j,UMX:UMZ)**2)/state(i,j,URHO))

               state(i,j,UEINT) = eint

               state(i,j,UFS) = state(i,j,URHO)

               state(i,j,UTEMP) = eos_state % T
           end if

           state(i,j,k,UFA)  = dye
           state(i,j,k,UFS:UFS-1+nspec) = ONE / nspec

        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

end subroutine ca_initdata
