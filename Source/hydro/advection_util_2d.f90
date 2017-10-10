module advection_util_2d_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  private

  public normalize_species_fluxes

contains

  subroutine normalize_species_fluxes(  &
                    flux1,flux1_lo, flux1_hi, &
                    flux2, flux2_lo, flux2_hi, &
                    lo,hi)

    use network, only : nspec
    use meth_params_module, only : NVAR, URHO, UFS
    use bl_constants_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: lo(2),hi(2)
    integer          :: flux1_lo(3), flux1_hi(3)
    integer          :: flux2_lo(3), flux2_hi(3)
    real(rt)         :: flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),NVAR)
    real(rt)         :: flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2),NVAR)

    ! Local variables
    integer          :: i,j,n
    real(rt)         :: sum,fac

    do j = lo(2),hi(2)
       do i = lo(1),hi(1)+1
          sum = ZERO
          do n = UFS, UFS+nspec-1
             sum = sum + flux1(i,j,n)
          end do
          if (sum .ne. ZERO) then
             fac = flux1(i,j,URHO) / sum
          else
             fac = ONE
          end if
          do n = UFS, UFS+nspec-1
             flux1(i,j,n) = flux1(i,j,n) * fac
          end do
       end do
    end do
    do j = lo(2),hi(2)+1
       do i = lo(1),hi(1)
          sum = ZERO
          do n = UFS, UFS+nspec-1
             sum = sum + flux2(i,j,n)
          end do
          if (sum .ne. ZERO) then
             fac = flux2(i,j,URHO) / sum
          else
             fac = ONE
          end if
          do n = UFS, UFS+nspec-1
             flux2(i,j,n) = flux2(i,j,n) * fac
          end do
       end do
    end do

  end subroutine normalize_species_fluxes


end module advection_util_2d_module
