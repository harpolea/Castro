subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use amrex_constants_module, only: ZERO, HALF, ONE
  use prob_params_module, only: center
  use castro_error_module, only: castro_error
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer,  intent(in) :: init, namlen
  integer,  intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  call probdata_init(name, namlen)

  ! set explosion center
  center(:) = HALF * (problo(:) + probhi(:))

end subroutine amrex_probinit


subroutine ca_initdata(level, time, lo, hi, nscal, &
                       state, s_lo, s_hi, &
                       dx, xlo, xhi)

  use probdata_module, only: rhol, rhor, pl, pr, vl, vr
  use amrex_constants_module, only: M_PI, FOUR3RD, ZERO, HALF, ONE
  use meth_params_module , only: NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UFS
  use prob_params_module, only: center, coord_type, problo
  use amrex_fort_module, only: rt => amrex_real
  use network, only: nspec
  use extern_probin_module, only: eos_gamma

  implicit none

  integer,  intent(in   ) :: level, nscal
  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: s_lo(3), s_hi(3)
  real(rt), intent(in   ) :: xlo(3), xhi(3), time, dx(3)
  real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), NVAR)

  real(rt) :: xx, yy, zz
  real(rt) :: v, p, h, W, tau
  integer  :: i, j, k

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           xx = problo(1) + dx(1) * (dble(i) + HALF)

           state(i,j,k,UMX:UMZ) = 0.e0_rt

           if (xx <= HALF) then
              state(i,j,k,URHO) = rhol
              p = pl
              v = vl
           else
              state(i,j,k,URHO) = rhor
              p = pr
              v = vr
           endif

           ! calculate the Lorentz factor 
           W = 1 / sqrt(1 - v*v)

           ! entropy 
           h = 1.0e0_rt + eos_gamma * p / state(i,j,k,URHO)

           state(i,j,k,URHO) = state(i,j,k,URHO) * W

           state(i,j,k,UMX) = state(i,j,k,URHO) * h * W * v

           state(i,j,k,UEDEN) = state(i,j,k,URHO) * h * W - p - state(i,j,k,URHO) 
           state(i,j,k,UEINT) = state(i,j,k,UEDEN)

           state(i,j,k,UFS:UFS-1+nspec) = state(i,j,k,URHO)

        end do
     end do
  end do

end subroutine ca_initdata
