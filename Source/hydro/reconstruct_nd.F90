module reconstruct_module

    use amrex_fort_module, only : rt => amrex_real

    implicit none

    public compute_reconstruction_tvd

contains

    subroutine compute_reconstruction_tvd(s, s_lo, s_hi, &
                                   sxm, sxp, sym, syp, szm, szp, sd_lo, sd_hi, &
                                   lo, hi, dx, k3d, kc)

        use mempool_module, only : bl_allocate, bl_deallocate
        use prob_params_module, only : dg
        use bl_constants_module
        implicit none

        integer, intent(in) ::  s_lo(3),  s_hi(3)
        integer, intent(in) ::  sd_lo(3),  sd_hi(3)
        integer, intent(in) :: lo(3), hi(3)
        integer, intent(in) :: k3d, kc

        real(rt)        , intent(in) ::     s( s_lo(1): s_hi(1), s_lo(2): s_hi(2), s_lo(3): s_hi(3))
        real(rt)        , intent(inout) :: sxm( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
        real(rt)        , intent(inout) :: sxp( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
        real(rt)        , intent(inout) :: sym( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
        real(rt)        , intent(inout) :: syp( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
        real(rt)        , intent(inout) :: szm( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
        real(rt)        , intent(inout) :: szp( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
        real(rt)        , intent(in) :: dx(3)

        integer i,j,k

        real(rt), pointer :: slope_l(:,:), slope_r(:,:)

        real(rt) :: slope

        call bl_allocate(slope_l, lo(1)-2, hi(1)+2, lo(2)-2*dg(2), hi(2)+2*dg(2))
        call bl_allocate(slope_r, lo(1)-2, hi(1)+2, lo(2)-2*dg(2), hi(2)+2*dg(2))

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! x-direction
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! compute s at x-edges

        slope_l = 0.0d0
        slope_r = 0.0d0

            do j=lo(2)-dg(2),hi(2)+dg(2)
                do i=lo(1)-2,hi(1)+1

                    slope_l(i,j) = s(i,j,k3d) - s(i-1,j,k3d)
                    slope_r(i,j) = s(i+1,j,k3d) - s(i,j,k3d)

                end do
            end do

        !!$    Use Van Leer MC limiter

        do j=lo(2)-dg(2),hi(2)+dg(2)
            do i=lo(1)-2,hi(1)+1

                if (slope_l(i,j) * slope_r(i,j) < 0.0d0) then
                    slope = 0.0d0
                else
                    slope = sign( min( 2.0d0 * abs(slope_l(i,j)), &
                               2.0d0 * abs(slope_r(i,j)), &
                              0.5d0 * (abs(slope_l(i,j)) + abs(slope_r(i,j))) ), &
                           slope_l(i,j) + slope_r(i,j) )
                end if

                sxm(i,j,kc) = s(i,j,k3d) - 0.5d0 * slope
                sxp(i,j,kc) = s(i,j,k3d) + 0.5d0 * slope

            end do
        end do

#if (BL_SPACEDIM >= 2)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! y-direction
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        slope_l = 0.0d0
        slope_r = 0.0d0

        do j=lo(2)-2,hi(2)+1
            do i=lo(1)-2,hi(1)+1

                slope_l(i,j) = s(i,j,k3d) - s(i,j-1,k3d)
                slope_r(i,j) = s(i,j+1,k3d) - s(i,j,k3d)

            end do
        end do

        !!$    Use Van Leer MC limiter

        do j=lo(2)-2,hi(2)+1
            do i=lo(1)-2,hi(1)+1

                if (slope_l(i,j) * slope_r(i,j) < 0.0d0) then
                    slope = 0.0d0
                else
                    slope = sign( min( 2.0d0 * abs(slope_l(i,j)), &
                                       2.0d0 * abs(slope_r(i,j)), &
                                      0.5d0 * (abs(slope_l(i,j)) + abs(slope_r(i,j))) ), &
                                   slope_l(i,j) + slope_r(i,j) )
                end if

                sym(i,j,kc) = s(i,j,k3d) - 0.5d0 * slope
                syp(i,j,kc) = s(i,j,k3d) + 0.5d0 * slope
            end do
        end do
#endif

#if (BL_SPACEDIM == 3)
! NOTE: might need to fix this. Check out the ppm_nd one.
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! z-direction
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        slope_l = 0.0d0
        slope_r = 0.0d0

        k=k3d
        do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1

                slope_l(i,j,k) = s(i,j,k) - s(i,j,k-1)
                slope_r(i,j,k) = s(i,j,k+1) - s(i,j,k)

            end do
        end do

        !!$    Use Van Leer MC limiter

        do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1

                if (slope_l(i,j,k) * slope_r(i,j,k) < 0.0d0) then
                    slope = 0.0d0
                else
                    slope = sign( min( 2.0d0 * abs(slope_l(i,j,k)), &
                                       2.0d0 * abs(slope_r(i,j,k)), &
                                      0.5d0 * (abs(slope_l(i,j,k)) + abs(slope_r(i,j,k))) ), &
                                   slope_l(i,j,k) + slope_r(i,j,k) )
                end if

                szm(i,j,kc) = s(i,j,k) - 0.5d0 * slope
                szp(i,j,kc) = s(i,j,k) + 0.5d0 * slope
            end do
        end do
#endif
        call bl_deallocate(slope_l)
        call bl_deallocate(slope_r)

    end subroutine compute_reconstruction_tvd

end module reconstruct_module
