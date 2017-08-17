module metric_module
    ! Define routines used to describe the metric

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

    subroutine calculate_alpha(alpha, alo, ahi)
        implicit none

        integer, intent(in) :: alo(3), ahi(3)
        real(rt), intent(out) :: alpha(alo(1):ahi(1), alo(2):ahi(2), alo(3):ahi(3))

        alpha(:,:,:) = 1.0d0

    end subroutine calculate_alpha

    subroutine calculate_beta(beta, blo, bhi)
        implicit none

        integer, intent(in) :: blo(3), bhi(3)
        real(rt), intent(out) :: beta(blo(1):bhi(1), blo(2):bhi(2), blo(3):bhi(3), 3)

        beta(:,:,:,:) = 0.0d0

    end subroutine calculate_beta

end module metric_module
