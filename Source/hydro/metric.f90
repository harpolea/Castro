module metric_module
    ! Define routines used to describe the metric

    use amrex_fort_module, only : rt => amrex_real
    implicit none

contains

    subroutine calculate_scalar_W(v, W)
        ! v is a 3-vector of velocities

        real(rt), intent(in) :: v(3)
        real(rt), intent(out) :: W

        integer :: i, j, k

        W = 1.0e0_rt / sqrt(1.0e0_rt - sum(v**2))

    end subroutine calculate_scalar_W

end module metric_module
