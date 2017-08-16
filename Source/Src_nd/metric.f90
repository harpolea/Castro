module metric_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

contains

    subroutine calculate_gamma_up(gamma_up, glo, ghi)
        implicit none

        integer, intent(in) :: glo(3), ghi(3)
        real(rt), intent(out) :: gamma_up(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), 9)

        gamma_up(:,:,:,:) = 0.0d0
        gamma_up(:,:,:,1) = 1.0d0
        gamma_up(:,:,:,5) = 1.0d0
        gamma_up(:,:,:,9) = 1.0d0

    end subroutine calculate_gamma_up

end module metric_module
