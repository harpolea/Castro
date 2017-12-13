module reconstruct_module

    use amrex_fort_module, only : rt => amrex_real

    implicit none

    public compute_reconstruction_tvd

contains

    subroutine compute_reconstruction_tvd(s, s_lo, s_hi, &
                                   sxm, sxp, &
#if BL_SPACEDIM >=2
                                   sym, syp, &
#endif
#if BL_SPACEDIM == 3
                                   szm, szp, &
#endif
                                   sd_lo, sd_hi, &
                                   lo, hi, dx)!, k3d, kc, km)

        use mempool_module, only : bl_allocate, bl_deallocate
        use prob_params_module, only : dg
        use bl_constants_module
        implicit none

        integer, intent(in) ::  s_lo(3),  s_hi(3)
        integer, intent(in) ::  sd_lo(3),  sd_hi(3)
        integer, intent(in) :: lo(3), hi(3)
        !integer, intent(in) :: k3d, kc, km

        real(rt)        , intent(in) ::     s( s_lo(1): s_hi(1), s_lo(2): s_hi(2), s_lo(3): s_hi(3))
        real(rt)        , intent(inout) :: sxm( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
        real(rt)        , intent(inout) :: sxp( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
#if BL_SPACEDIM >=2
        real(rt)        , intent(inout) :: sym( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
        real(rt)        , intent(inout) :: syp( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
#endif
#if BL_SPACEDIM == 3
        real(rt)        , intent(inout) :: szm( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
        real(rt)        , intent(inout) :: szp( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
#endif
        real(rt)        , intent(in) :: dx(3)

        integer i,j,k

        real(rt) :: slope_l, slope_r, slope, r

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! x-direction
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! compute s at x-edges

        slope_l = 0.0d0
        slope_r = 0.0d0
#if (BL_SPACEDIM <= 2)
        k = 0!lo(3)
#else
        do k=lo(3)-dg(3),hi(3)+dg(3)
#endif
            do j=lo(2)-dg(2),hi(2)+dg(2)
                do i=lo(1)-2,hi(1)+1

                    !write(*,*) "i, j, s(i-1, j)", i, j, s(i-1, j, k)

                    slope_l = s(i,j,k) - s(i-1,j,k)
                    slope_r = s(i+1,j,k) - s(i,j,k)

                    ! if (slope_l /= slope_l .or. slope_l-1 .eq. slope_l) then
                    !     slope_l = 0.0d0
                    ! endif
                    ! if (slope_r /= slope_r .or. slope_r-1 .eq. slope_r) then
                    !     slope_r = 0.0d0
                    ! endif

                    !!$    Use Van Leer MC limiter
                    !write(*,*) "slope_l, slope_r", slope_l, slope_r

                    if (slope_l * slope_r < 0.0d0) then
                        slope = 0.0d0
                    else
                        r = slope_r / slope_l
                        if (r /= r) then
                            r = 0.0d0
                        endif
                        ! minmod
                        ! slope = 0.5d0 * (slope_l + slope_r) * min(1.0d0, 4.0d0 / (1.0d0 + r))

                        ! van leer
                        slope = 0.5d0 * (slope_l + slope_r) * min(2.0d0 * r/(1.0d0 + r), 2.0d0 / (1.0d0 + r))


                        ! slope = sign( min( 2.0d0 * abs(slope_l), &
                        !            2.0d0 * abs(slope_r), &
                        !           0.5d0 * (abs(slope_l) + abs(slope_r))), &
                        !        slope_l + slope_r)
                    end if

                    if (slope /= slope) then
                        slope = 0.0d0
                    endif

                    sxm(i,j,k) = s(i,j,k) - 0.5d0 * slope
                    sxp(i,j,k) = s(i,j,k) + 0.5d0 * slope

                end do
            end do

#if (BL_SPACEDIM >= 2)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! y-direction
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        slope_l = 0.0d0
        slope_r = 0.0d0

        do j=lo(2)-2,hi(2)+1
            do i=lo(1)-1,hi(1)+1

                slope_l = s(i,j,k) - s(i,j-1,k)
                slope_r = s(i,j+1,k) - s(i,j,k)

                ! if (slope_l /= slope_l .or. slope_l-1 .eq. slope_l) then
                !     slope_l = 0.0d0
                ! endif
                ! if (slope_r /= slope_r .or. slope_r-1 .eq. slope_r) then
                !     slope_r = 0.0d0
                ! endif

                !!$    Use Van Leer MC limiter

                if (slope_l * slope_r < 0.0d0) then
                    slope = 0.0d0
                else
                    r = slope_r / slope_l
                    if (r /= r) then
                        r = 0.0d0
                    endif
                    ! minmod
                    ! slope = 0.5d0 * (slope_l + slope_r) * min(1.0d0, 4.0d0 / (1.0d0 + r))

                    ! van leer
                    slope = 0.5d0 * (slope_l + slope_r) * min(2.0d0 * r/(1.0d0 + r), 2.0d0 / (1.0d0 + r))


                    ! slope = sign( min( 2.0d0 * abs(slope_l), &
                    !            2.0d0 * abs(slope_r), &
                    !           0.5d0 * (abs(slope_l) + abs(slope_r))), &
                    !        slope_l + slope_r)
                end if

                if (slope /= slope) then
                    slope = 0.0d0
                endif

                sym(i,j,k) = s(i,j,k) - 0.5d0 * slope
                syp(i,j,k) = s(i,j,k) + 0.5d0 * slope
            end do
        end do
#endif
#if (BL_SPACEDIM == 3)
    end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! z-direction
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        slope_l = 0.0d0
        slope_r = 0.0d0

        do k=lo(3)-2,hi(3)+1
            do j=lo(2)-1,hi(2)+1
                do i=lo(1)-1,hi(1)+1

                    slope_l = s(i,j,k) - s(i,j,k-1)
                    slope_r = s(i,j,k+1) - s(i,j,k)

                    ! if (slope_l /= slope_l .or. slope_l-1 .eq. slope_l) then
                    !     slope_l = 0.0d0
                    ! endif
                    ! if (slope_r /= slope_r .or. slope_r-1 .eq. slope_r) then
                    !     slope_r = 0.0d0
                    ! endif

                    !!$    Use Van Leer MC limiter

                    if (slope_l * slope_r < 0.0d0) then
                        slope = 0.0d0
                    else
                        r = slope_r / slope_l
                        if (r /= r) then
                            r = 0.0d0
                        endif
                        ! minmod
                        ! slope = 0.5d0 * (slope_l + slope_r) * min(1.0d0, 4.0d0 / (1.0d0 + r))

                        ! van leer
                        slope = 0.5d0 * (slope_l + slope_r) * min(2.0d0 * r/(1.0d0 + r), 2.0d0 / (1.0d0 + r))


                        ! slope = sign( min( 2.0d0 * abs(slope_l), &
                        !            2.0d0 * abs(slope_r), &
                        !           0.5d0 * (abs(slope_l) + abs(slope_r))), &
                        !        slope_l + slope_r)
                    end if

                    szm(i,j,k) = s(i,j,k) - 0.5d0 * slope
                    szp(i,j,k) = s(i,j,k) + 0.5d0 * slope
                end do
            end do
        end do
#endif

    end subroutine compute_reconstruction_tvd

end module reconstruct_module
