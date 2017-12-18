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
        nsub, temp_ambient

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

   h_in = 2.0d0
   h_out = 1.0d0
   damn_rad = 0.2d0

   g = 1.0d0
   p_ambient = 1.e-5_rt        ! ambient pressure (in erg/cc)
   dens_ambient = 1.e0_rt      ! ambient density (in g/cc)
   exp_energy = 1.e0_rt        ! absolute energy of the explosion (in erg)
   nsub = 4
   temp_ambient = -1.e2_rt

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

      eos_state % rho = dens_ambient
      eos_state % xn(:) = xn_zone(:)
      eos_state % T = temp_ambient

      call eos(eos_input_rt, eos_state)

      p_ambient = eos_state % p

   endif

   ! Calculate ambient state data

   eos_state % rho = dens_ambient
   eos_state % p   = p_ambient
   eos_state % T   = 1.d5 ! Initial guess for iterations
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
       UFS, UFA, UEDEN, UEINT, UTEMP
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

  integer :: i, j, k, n, ii, jj, kk

  real(rt)         :: dye, eint, xmin, ymin
  type(eos_t) :: eos_state

  !state(:,:,:,:) = 0.0d0
  dye = ZERO

  if (level <= swe_to_comp_level) then
      write(*,*) "Initialising level ", level, " with SWE data"
  else
      write(*,*) "Initialising level ", level, " with compressible data"
  endif

  !$OMP PARALLEL DO PRIVATE(i, j, k, xx, yy, zz, r)
  do k = lo(3), hi(3)
     zz = xlo(3) + delta(3)*dble(k-lo(3)+HALF)

     do j = lo(2), hi(2)
        ymin = xlo(2) + delta(2)*dble(j-lo(2))
        yy = xlo(2) + delta(2)*dble(j-lo(2)+HALF)

        do i = lo(1), hi(1)

            state(i,j,k,1:NVAR) = 0.0d0

            xmin = xlo(1) + delta(1)*dble(i-lo(1))
            xx = xlo(1) + delta(1)*dble(i-lo(1)+HALF)

            r = yy

           ! r = sqrt((xx - center(1))**2 + (yy - center(2))**2)

           ! if (level <= swe_to_comp_level) then
           !     ! shallow water level
           !     if (r < damn_rad) then
           !         state(i,j,k,URHO) = h_in
           !     else
           !         state(i,j,k,URHO) = h_out
           !     end if
           !
           !     !do n = 1,nspec
           !      !  state(i,j,k,UFS+n-1) = state(i,j,k,URHO) * state(i,j,k,UFS+n-1)
           !     !end do
           !
           ! else ! compressible level

                state(i,j,k,URHO) = 1.0d0

                if (r < damn_rad) then
                    eos_state % p = 0.5 * g * (h_in-xx)**2
                else
                    eos_state % p = 0.5 * g * (h_out-xx)**2
                end if

                ! if (level > swe_to_comp_level) then
                !     state(i,j,k,URHO) = 2.0d0 * state(i,j,k,URHO)
                ! endif


                !eos_state % e = e_zone
                eos_state % rho = state(i,j,k,URHO)
                eos_state % xn(:) = xn_zone(:)
                !eos_state % p = 0.5d0 * g * state(i,j,k,URHO)**2

                call eos(eos_input_rp, eos_state)

                eint = state(i,j,k,URHO) * eos_state % e
                state(i,j,k,UMX:UMZ) = 0.e0_rt

                state(i,j,k,UEDEN) = eint +  &
                    0.5e0_rt*(sum(state(i,j,k,UMX:UMZ)**2)/state(i,j,k,URHO))

                if (level <= swe_to_comp_level) then
                    ! shallow water level
                    if (r < damn_rad) then
                        state(i,j,k,URHO) = h_in
                    else
                        state(i,j,k,URHO) = h_out
                    end if
                end if

                state(i,j,k,UEINT) = eint

                state(i,j,k,UFS) = state(i,j,k,URHO)

                state(i,j,k,UTEMP) = eos_state % T
          ! end if

           state(i,j,k,UFA)  = dye
           state(i,j,k,UFS:UFS-1+nspec) = state(i,j,k,URHO) / nspec

        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

end subroutine ca_initdata
