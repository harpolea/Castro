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

   namelist /fortin/ h_in, h_out, damn_rad, g

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

  state(:,:,:,:) = 0.0d0
  dye = ZERO

  !$OMP PARALLEL DO PRIVATE(i, j, k, xx, yy, zz, r)
  do k = lo(3), hi(3)
     zz = xlo(3) + delta(3)*dble(k-lo(3)+HALF)

     do j = lo(2), hi(2)
        yy = xlo(2) + delta(2)*dble(j-lo(2)+HALF)

        do i = lo(1), hi(1)
           xx = xlo(1) + delta(1)*dble(i-lo(1)+HALF)

           r = xx!sqrt((xx - center(1))**2 + (yy - center(2))**2)

           if (r < damn_rad) then
               state(i,j,k,URHO) = h_in
           else
               state(i,j,k,URHO) = h_out
           end if

           state(i,j,k,UFA)  = dye
           state(i,j,k,UFS:UFS-1+nspec) = ONE / nspec

           do n = 1,nspec
              state(i,j,k,UFS+n-1) = state(i,j,k,URHO) * state(i,j,k,UFS+n-1)
           end do

        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

end subroutine ca_initdata
