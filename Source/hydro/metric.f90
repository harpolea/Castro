module metric_module
    ! Define routines used to describe the metric

    use amrex_fort_module, only : rt => amrex_real
    implicit none

contains

    subroutine calculate_gamma_up(lo, hi, gamma_up, glo, ghi, phi, p_lo, p_hi)

        integer, intent(in) :: lo(3), hi(3), glo(3), ghi(3), p_lo(3), p_hi(3)
        real(rt), intent(in) :: phi(p_lo(1):p_hi(1), p_lo(2):p_hi(2), p_lo(3):p_hi(3))
        real(rt), intent(out) :: gamma_up(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), 9)

        integer :: i,j,k

        gamma_up(:,:,:,:) = 0.0d0
        gamma_up(:,:,:,1) = 1.0d0
        gamma_up(:,:,:,5) = 1.0d0

        do k = lo(3), hi(3)
            do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                    gamma_up(i,j,k,9) = exp(-2.0d0 * phi(i,j,k))
                end do
            end do
        end do

    end subroutine calculate_gamma_up

    subroutine calculate_alpha(lo, hi, alpha, alo, ahi, phi, p_lo, p_hi)

        integer, intent(in) :: lo(3), hi(3), alo(3), ahi(3), p_lo(3), p_hi(3)
        real(rt), intent(in) :: phi(p_lo(1):p_hi(1), p_lo(2):p_hi(2), p_lo(3):p_hi(3))
        real(rt), intent(out) :: alpha(alo(1):ahi(1), alo(2):ahi(2), alo(3):ahi(3))

        integer :: i,j,k

        do k = lo(3), hi(3)
            do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                    alpha(i,j,k) = exp(-phi(i,j,k))
                end do
            end do
        end do

    end subroutine calculate_alpha

    subroutine calculate_beta(beta, blo, bhi)
        implicit none

        integer, intent(in) :: blo(3), bhi(3)
        real(rt), intent(out) :: beta(blo(1):bhi(1), blo(2):bhi(2), blo(3):bhi(3), 3)

        beta(:,:,:,:) = 0.0d0

    end subroutine calculate_beta

    subroutine calculate_norm(v, gamma_up, n)
        ! v is a 3-vector. Returns v_i * v^i

        real(rt), intent(in) :: v(3)
        real(rt), intent(in) :: gamma_up(9)
        real(rt), intent(out) :: n

        n = v(1)**2 * gamma_up(1) + &
            2.0d0 * v(1) * v(2) * gamma_up(2) + &
            2.0d0 * v(1) * v(3) * gamma_up(3) + &
            v(2)**2 * gamma_up(5) + &
            2.0d0 * v(2) * v(3) * gamma_up(6) + &
            v(3)**2 * gamma_up(9)

    end subroutine calculate_norm

    subroutine calculate_scalar_W(v, gamma_up, W)
        ! v is a 3-vector of velocities

        real(rt), intent(in) :: v(3)
        real(rt), intent(in) :: gamma_up(9)
        real(rt), intent(out) :: W

        integer :: i, j, k

        W = v(1)**2 * gamma_up(1) + &
            2.0d0 * v(1) * v(2) * gamma_up(2) + &
            2.0d0 * v(1) * v(3) * gamma_up(3) + &
            v(2)**2 * gamma_up(5) + &
            2.0d0 * v(2) * v(3) * gamma_up(6) + &
            v(3)**2 * gamma_up(9)

        W = 1.0d0 / sqrt(1.0d0 - W)

    end subroutine calculate_scalar_W

end module metric_module
