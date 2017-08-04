module gr_utils_module

  implicit none

  public

contains

    subroutine W_swe(q, slo, shi, lo, hi, Ncomp, gamma_up, glo, ghi, &
                     W, wlo, whi)
        ! Calculate lorentz factor
        implicit none

        integer, intent(in) :: slo(2), shi(2), lo(2), hi(2), Ncomp, glo(2), ghi(2), wlo(2), whi(2)
        double precision, intent(in) :: q(slo(1):shi(1), slo(2):shi(2), Ncomp)
        double precision, intent(in) :: gamma_up(glo(1):ghi(1), glo(2):ghi(2), 9)
        double precision, intent(out) :: W(wlo(1):whi(1), wlo(2):whi(2))

        integer i,j

        do     j = lo(2), hi(2)
            do i = lo(1), hi(1)
                if (q(i,j,1) < 1.d-20) then
                    W(i,j) = 1.0d0
                else
                    W(i,j) = sqrt((q(i,j,2)**2 * gamma_up(i,j,1)+&
                          2.0d0 * q(i,j,2) * q(i,j,3) * &
                          gamma_up(i,j,2) + q(i,j,3)**2 * &
                          gamma_up(i,j,5)) / q(i,j,1)**2 + 1.0d0)
                end if
                ! nan check
                if (W(i,j) /= W(i,j)) then
                    W(i,j) = 1.0d0
                end if
          end do
      end do

    end subroutine W_swe

    subroutine swe_from_comp(U_prim, prlo, prhi, U_swe, slo, shi, p_comp, &
         pclo, pchi, p_swe, lo, hi, n_cons_comp, n_swe_comp, &
         alpha0, M, R, dx, prob_lo) bind(C, name="swe_from_comp")
        ! Assume nlayers = 1 as 2d
        implicit none

        integer, intent(in) :: n_cons_comp, n_swe_comp
        integer, intent(in) :: prlo(2), prhi(2), slo(2), shi(2), pclo(2), pchi(2), lo(2), hi(2)
        double precision, intent(inout)  :: U_swe(slo(1):shi(1), slo(2):shi(2), n_swe_comp)
        double precision, intent(in) :: U_prim(prlo(1):prhi(1), prlo(2):prhi(2), n_cons_comp)
        double precision, intent(in) :: p_comp(pclo(1):pchi(1), pclo(2):pchi(2))
        double precision, intent(in) :: p_swe
        double precision, intent(in) :: alpha0, M, R, dx(2), prob_lo(2)

        double precision gamma_up(lo(1):hi(1), lo(2):hi(2), 9)
        double precision h_comp, ssq
        integer neighbour, minl(1)
        double precision zfrac, W, h
        integer i, j, k

        call calc_gamma_up(gamma_up, lo, hi, lo, hi, alpha0, M, R, dx, prob_lo)

        !write(*,*) "p_swe = ", p_swe(slo(3))

        ! neighbour is the comp layer above
        ! check if this layer is above or below
        do     j = lo(2), hi(2)
            do i = lo(1), hi(1)



                !write(*,*) "zfrac, neighbour", zfrac, neighbour, h_comp(neighbour)

                ! now interpolate and stick primitive compressible variables in U_comp
                ! TODO: slope limit here?
                U_swe(i,j,1) = h_comp
                U_swe(i,j,2:3) = U_prim(i,j,2:3)
                h = h_comp
                if (n_swe_comp == 4) then
                    U_swe(i,j,4) = h
                end if

                ! interpolate W
                ! NOTE: do I interpolate the primitive velocities then calculate W
                ! or do as I've done here and calculate W then interpolate??
                ssq = U_prim(i,j,2)**2 * gamma_up(i,j,1) + &
                    2.0d0 * U_prim(i,j,2) * U_prim(i,j,3) * &
                        gamma_up(i,j,2) + &
                    2.0d0 * U_prim(i,j,2) * U_prim(i,j,4) * &
                        gamma_up(i,j,3) + &
                    U_prim(i,j,3)**2 * gamma_up(i,j,5) + &
                    2.0d0 * U_prim(i,j,3) * U_prim(i,j,4) * &
                        gamma_up(i,j,6) + &
                    U_prim(i,j,4)**2 * gamma_up(i,j,9)

                W = 1.0d0 / sqrt(1.0d0 - ssq)

                U_swe(i,j,1) = -log(alpha0 * (1.0d0 - h / R) + h / R) * W !-log(alpha0 + M * h / (R * alpha0)) * W
                U_swe(i,j,2) = U_swe(i,j,1) * W * U_swe(i,j,2)
                U_swe(i,j,3) = U_swe(i,j,1) * W * U_swe(i,j,3)
            end do
        end do

        write(*,*) "U_swe", U_swe(lo(1), lo(2), :)
        write(*,*) "U_prim", U_prim(lo(1), lo(2), :)

        ! calculate conserved variables
        !gamma_up = 0.0d0
        !gamma_up(:,:,1) = 1.0d0
        !gamma_up(:,:,5) = 1.0d0
        ! don't need this component
        !gamma_down(:,:,9) = 1.0d0 / (alpha0 + M * U_swe(:,:,4) / (R**2 * alpha0))**2

        !W(:,:) = U_swe(:,:,2)**2 * gamma_down(:,:,1) + &
        !    2.0d0 * U_swe(:,:,2) * U_swe(:,:,3) * gamma_down(:,:,2) + &
        !    U_swe(:,:,3)**2 * gamma_down(:,:,5)
        !W = 1.0d0 / sqrt (1.0d0 - W)

    end subroutine swe_from_comp

    subroutine calc_gamma_up_swe(U, slo, shi, lo, hi, n_comp, gamma_up)
        implicit none

        integer, intent(in) :: n_comp
        integer, intent(in) :: slo(2), shi(2), lo(2), hi(2)
        double precision, intent(in) :: U(slo(1):shi(1), slo(2):shi(2), n_comp)
        double precision, intent(out) :: gamma_up(lo(1):hi(1), lo(2):hi(2), 9)

        gamma_up = 0.0d0

        gamma_up(:,:,1) = 1.0d0
        gamma_up(:,:,5) = 1.0d0
        gamma_up(:,:,9) = exp(-2.0d0 * U(lo(1):hi(1), lo(2):hi(2),1))

    end subroutine calc_gamma_up_swe

    subroutine comp_from_swe(U_comp, clo, chi, U_swe, slo, shi, p, &
        rho, lo, hi, n_cons_comp, n_swe_comp, gamma, dx, alpha0, &
        M, R, nghost, prob_lo) bind(C, name="comp_from_swe")
        ! TODO: what do I do about vertical velocity component????
        !use slope_module, only: uslope
        implicit none

        integer, intent(in) :: n_cons_comp, n_swe_comp, nghost
        integer, intent(in) :: clo(2), chi(2), slo(2), shi(2), lo(2), hi(2)
        double precision, intent(in)  :: U_swe(slo(1):shi(1), slo(2):shi(2), n_swe_comp)
        double precision, intent(out) :: U_comp(clo(1):chi(1), clo(2):chi(2), n_cons_comp)
        double precision, intent(in) :: p, prob_lo(2)
        double precision, intent(in) :: rho
        double precision, intent(in)  :: gamma, dx(2), alpha0, M, R

        double precision h_swe(lo(1)-nghost:hi(1)+nghost, lo(2)-nghost:hi(2)+nghost)
        double precision v_swe(lo(1)-nghost:hi(1)+nghost, lo(2)-nghost:hi(2)+nghost, 2)
        double precision h_comp
        double precision zfrac
        integer neighbour, minl(1)
        integer i, j, k, nlo(2), nhi(2)
        double precision W(lo(1)-nghost:hi(1)+nghost, lo(2)-nghost:hi(2)+nghost)
        double precision rhoh(lo(1)-nghost:hi(1)+nghost, lo(2)-nghost:hi(2)+nghost)

        double precision gamma_up_swe(lo(1)-nghost:hi(1)+nghost, lo(2)-nghost:hi(2)+nghost, 9)
        double precision gamma_up(lo(1)-nghost:hi(1)+nghost, lo(2)-nghost:hi(2)+nghost, 9)

        nlo = lo - nghost
        nhi = hi + nghost

        write(*,*) "comp from swe"
        write(*,*) "U_swe: ", U_swe(slo(1)+nghost, slo(2)+nghost, :), slo(1)+nghost, slo(2)+nghost

        call calc_gamma_up_swe(U_swe, slo, shi, nlo, nhi, n_swe_comp, gamma_up_swe)
        call calc_gamma_up(gamma_up, nlo, nhi, nlo, nhi, alpha0, M, R, dx, prob_lo)
        call W_swe(U_swe, slo, shi, nlo, nhi, n_swe_comp, gamma_up_swe, nlo, nhi, W(nlo(1):nhi(1), nlo(2):nhi(2)), nlo, nhi)

        !write(*,*) "gamma_up", gamma_up_swe(lo(1), lo(2), lo(3), 7:9)
        !write(*,*) "W_swe", W(lo(1), lo(2), slo(3))
        !write(*,*) "p", p

        write(*,*) "hcomp: ", h_comp
        write(*,*) "U, exp(-U/w), alpha, W", U_swe(nlo(1), nlo(2),1), exp(-U_swe(nlo(1), nlo(2),1) / W(nlo(1), nlo(2))), alpha0, W(nlo(1), nlo(2))

        if (n_swe_comp > 3) then
            h_swe = U_swe(nlo(1):nhi(1),nlo(2):nhi(2),4)
        else
            ! HACK
            h_swe = R * (exp(-U_swe(nlo(1):nhi(1), nlo(2):nhi(2),1) / W(nlo(1):nhi(1), nlo(2):nhi(2))) - alpha0) / (1.0d0 - alpha0) !alpha0 * R / M * (exp(-U_swe(nlo(1):nhi(1), nlo(2):nhi(2),slo(3):shi(3),1) / W(nlo(1):nhi(1), nlo(2):nhi(2),slo(3):shi(3))) - alpha0)
        end if

            do     j = nlo(2), nhi(2)
                do i = nlo(1), nhi(1)
                    if (U_swe(i,j,1) < 1.d-20) then
                        v_swe(i,j,1) = U_swe(i,j,2)
                        v_swe(i,j,2) = U_swe(i,j,3)
                    else
                        v_swe(i,j,1) = U_swe(i,j,2) / &
                            (W(i,j) * U_swe(i,j,1))
                        v_swe(i,j,2) = U_swe(i,j,3) / &
                            (W(i,j) * U_swe(i,j,1))
                    end if
                end do
            end do

        U_comp(:,:,:) = 0.0d0

        ! calculate layer fracs and interpolate
            ! neighbour is the swe layer below
            ! check if this layer is above or below
            do     j = nlo(2), nhi(2)
                do i = nlo(1), nhi(1)

                    ! TODO: slope limit here?
                    U_comp(i,j,1) = rho
                    U_comp(i,j,2:3) = v_swe(i,j,:)

                    U_comp(i,j,5) = p
                    !write(*,*) "p: ", U_comp(i,j,5), "p(neighbour+1)", p(neighbour+1)

                    W(i,j) = U_comp(i,j,2)**2*gamma_up(i,j,1) + &
                        2.0d0 * U_comp(i,j,2) * U_comp(i,j,3) * &
                            gamma_up(i,j,2) + &
                        2.0d0 * U_comp(i,j,2) * U_comp(i,j,4) * &
                            gamma_up(i,j,3) + &
                        U_comp(i,j,3)**2 * gamma_up(i,j,5) + &
                        2.0d0 * U_comp(i,j,3) * U_comp(i,j,4) * &
                            gamma_up(i,j,6) + &
                        U_comp(i,j,4)**2 * gamma_up(i,j,9)
                    !write(*,*) "U_comp = ", U_comp(i,j,:)
                    W(i,j) = 1.0d0 / sqrt(1.0d0 - W(i,j))
                end do
            end do

        !write(*,*) "U_comp", U_comp(lo(1), lo(2), lo(3), :)
        !write(*,*) "W", W(lo(1), lo(2), lo(3))

        call rhoh_from_p(rhoh, U_comp(:,:,5), U_comp(:,:,1), gamma, clo, chi, nlo, nhi)

            do     j = nlo(2), nhi(2)
                do i = nlo(1), nhi(1)
                    U_comp(i,j,1) = U_comp(i,j,1) * W(i,j)
                    U_comp(i,j,2) = rhoh(i,j) * W(i,j)**2 * U_comp(i,j,2)
                    U_comp(i,j,3) = rhoh(i,j) * W(i,j)**2 * U_comp(i,j,3)
                    U_comp(i,j,4) = rhoh(i,j) * W(i,j)**2 * U_comp(i,j,4)
                    U_comp(i,j,5) = rhoh(i,j) * W(i,j)**2 - &
                                        U_comp(i,j,5) - U_comp(i,j,1)
                end do
            end do

        write(*,*) "U_comp", U_comp(lo(1), lo(2), :)

    end subroutine comp_from_swe

    subroutine rhoh_from_p(rhoh, p, rho, gamma, clo, chi, lo, hi)
        implicit none

        integer, intent(in) :: clo(2), chi(2), lo(2), hi(2)
        double precision, intent(in)  :: p(clo(1):chi(1), clo(2):chi(2))
        double precision, intent(in)  :: rho(clo(1):chi(1), clo(2):chi(2))
        double precision, intent(in)  :: gamma
        double precision, intent(out)  :: rhoh(lo(1):hi(1), lo(2):hi(2))

        rhoh = rho(lo(1):hi(1), lo(2):hi(2)) + &
            gamma * p(lo(1):hi(1), lo(2):hi(2)) / (gamma - 1.0d0)

    end subroutine rhoh_from_p

    subroutine p_from_rhoh(rhoh, p, rho, gamma, lo, hi)
        implicit none

        integer, intent(in) :: lo(2), hi(2)
        double precision, intent(out)  :: p(lo(1):hi(1), lo(2):hi(2))
        double precision, intent(in)  :: rho(lo(1):hi(1), lo(2):hi(2))
        double precision, intent(in)  :: gamma
        double precision, intent(in)  :: rhoh(lo(1):hi(1), lo(2):hi(2))

        p = (rhoh - rho) * (gamma - 1.0d0) / gamma

    end subroutine p_from_rhoh

    subroutine p_from_rho_eps(rho, eps, p, gamma, lo, hi)
        implicit none

        integer, intent(in) :: lo(2), hi(2)
        double precision, intent(out)  :: p(lo(1):hi(1), lo(2):hi(2))
        double precision, intent(in)  :: rho(lo(1):hi(1), lo(2):hi(2))
        double precision, intent(in)  :: gamma
        double precision, intent(in)  :: eps(lo(1):hi(1), lo(2):hi(2))

        p = (gamma - 1.0d0) * rho * eps

    end subroutine p_from_rho_eps

    subroutine calc_gamma_up(gamma_up, glo, ghi, lo, hi, alpha0, &
        M, R, dx, prob_lo)
        implicit none

        integer, intent(in) :: glo(2), ghi(2), lo(2), hi(2)
        double precision, intent(out)  :: gamma_up(glo(1):ghi(1), glo(2):ghi(2), 9)
        double precision, intent(in)  :: alpha0, M, R
        double precision, intent(in)  :: dx(2), prob_lo(2)


        gamma_up(:,:,:) = 0.0d0
        gamma_up(:,:,1) = 1.0d0
        gamma_up(:,:,5) = 1.0d0

        gamma_up(:,:,9) = alpha0
    end subroutine calc_gamma_up

    subroutine calc_gamma_down(gamma_down, glo, ghi, lo, hi, alpha0, &
        M, R, dx, prob_lo)
        implicit none

        integer, intent(in) :: glo(2), ghi(2), lo(2), hi(2)
        double precision, intent(out)  :: gamma_down(glo(1):ghi(1), glo(2):ghi(2), 9)
        double precision, intent(in)  :: alpha0, M, R
        double precision, intent(in)  :: dx(2), prob_lo(2)

        call calc_gamma_up(gamma_down, glo, ghi, lo, hi, alpha0, M, R, dx, prob_lo)

        gamma_down(:,:,9) = 1.0 / gamma_down(:,:,9)
    end subroutine calc_gamma_down

    subroutine gr_sources(S, slo, shi, U, ulo, uhi, p, plo, phi, &
        alpha, alo, ahi, gamma_up, glo, ghi, M, R, gamma, Ncomp, &
        lo, hi, dx)
        implicit none

        integer, intent(in) :: slo(2), shi(2), ulo(2), uhi(2), plo(2), phi(2), alo(2), ahi(2), glo(2), ghi(2), lo(2), hi(2), Ncomp
        double precision, intent(inout)  :: S(slo(1):shi(1), slo(2):shi(2))
        double precision, intent(in)  :: U(ulo(1):uhi(1), ulo(2):uhi(2), Ncomp)
        double precision, intent(in)  :: p(plo(1):phi(1), plo(2):phi(2))
        double precision, intent(in)  :: alpha(alo(1):ahi(1), alo(2):ahi(2))
        double precision, intent(in)  :: gamma_up(glo(1):ghi(1), glo(2):ghi(2), 9)
        double precision, intent(in)  :: M, R, gamma, dx(2)

        double precision Ssq(lo(1):hi(1), lo(2):hi(2))
        double precision h(lo(1):hi(1), lo(2):hi(2))
        double precision W2(lo(1):hi(1), lo(2):hi(2))
        double precision S_temp, alpha0
        integer i,j

        Ssq(:,:) = U(lo(1):hi(1), lo(2):hi(2), 2)**2 * &
            gamma_up(lo(1):hi(1), lo(2):hi(2),1) + &
            2.0d0 * U(lo(1):hi(1), lo(2):hi(2), 2) * &
            U(lo(1):hi(1), lo(2):hi(2), 3) * &
            gamma_up(lo(1):hi(1), lo(2):hi(2),2) + &
            2.0d0 * U(lo(1):hi(1), lo(2):hi(2), 2) * &
            U(lo(1):hi(1), lo(2):hi(2), 4) * &
            gamma_up(lo(1):hi(1), lo(2):hi(2),3) + &
            U(lo(1):hi(1), lo(2):hi(2), 3)**2 * &
            gamma_up(lo(1):hi(1), lo(2):hi(2),5) + &
            2.0d0 * U(lo(1):hi(1), lo(2):hi(2), 3) * &
            U(lo(1):hi(1), lo(2):hi(2), 4) * &
            gamma_up(lo(1):hi(1), lo(2):hi(2),6) + &
            U(lo(1):hi(1), lo(2):hi(2), 4)**2 * &
            gamma_up(lo(1):hi(1), lo(2):hi(2),9)

        h = 1.0d0 + gamma * &
            (sqrt((U(lo(1):hi(1), lo(2):hi(2), 5) + &
            p(lo(1):hi(1), lo(2):hi(2)) + &
            U(lo(1):hi(1), lo(2):hi(2), 1))**2 - Ssq) - &
            p(lo(1):hi(1), lo(2):hi(2)) * &
            (U(lo(1):hi(1), lo(2):hi(2), 5) + &
            p(lo(1):hi(1), lo(2):hi(2)) + &
            U(lo(1):hi(1), lo(2):hi(2), 1)) / &
            sqrt((U(lo(1):hi(1), lo(2):hi(2), 5) + &
            p(lo(1):hi(1), lo(2):hi(2)) + &
            U(lo(1):hi(1), lo(2):hi(2), 1))**2 - Ssq) - &
            U(lo(1):hi(1), lo(2):hi(2), 1)) / &
            U(lo(1):hi(1), lo(2):hi(2), 1)

        W2 = 1.0d0 + Ssq / &
            (U(lo(1):hi(1), lo(2):hi(2), 1) * h)**2

        alpha0 = sqrt(1.0d0 - M / R)

            do     j = lo(2), hi(2)
                do i = lo(1), hi(1)
                    S_temp = -M / R**2 * &
                        (alpha(i,j) * U(i,j, 4)**2 / W2(i,j) + &
                        (U(i,j, 5) + p(i,j) + U(i,j, 1)) / &
                        alpha(i,j))
                    S_temp = -M / R**2 / alpha0**3 * alpha(i,j)
                    if (S_temp == S_temp) then ! not nan
                        S(i,j) = S(i,j) + S_temp! / dx(3)! * dx(3) * alpha(i,j)
                    end if

                end do
            end do

    end subroutine gr_sources



    subroutine zbrent(p, x1, b, U, Ncomp, gamma, gamma_up)
        use compute_flux_module, only : f_of_p
        ! route finder using brent's method
        implicit none

        integer, intent(in) :: Ncomp
        double precision, intent(out) :: p
        double precision, intent(in)  :: U(Ncomp), gamma, gamma_up(9), x1
        double precision, intent(inout) :: b

        double precision, parameter :: TOL = 1.0d-12
        integer, parameter :: ITMAX = 100

        double precision a, c, d, fa, fb, fc, fs, s
        logical mflag, con1, con2, con3, con4, con5
        integer i

        a = x1
        c = 0.0d0
        d = 0.0d0
        call f_of_p(fa, a, U, Ncomp, gamma, gamma_up)
        call f_of_p(fb, b, U, Ncomp, gamma, gamma_up)
        fc = 0.0d0

        if (fa * fb >= 0.0d0) then
            p = b
            return
        end if

        if (abs(fa) < abs(fb)) then
            d = a
            a = b
            b = d

            d = fa
            fa = fb
            fb = d
        end if

        c = a
        fc = fa

        mflag = .true.

        do i = 1, ITMAX
            if (fa /= fc .and. fb /= fc) then
                s = a*fb*fc / ((fa-fb) * (fa-fc)) + b*fa*fc / ((fb-fa)*(fb-fc)) +&
                    c*fa*fb / ((fc-fa)*(fc-fb))
            else
                s = b - fb * (b-a) / (fb-fa)
            end if

            con1 = .false.

            if (0.25d0 * (3.0d0 * a + b) < b) then
                if ( s < 0.25d0 * (3.0d0 * a + b) .or. s > b) then
                    con1 = .true.
                end if
            else if (s < b .or. s > 0.25d0  * (3.0d0 * a + b)) then
                con1 = .true.
            end if

            con2 = mflag .and. abs(s - b) >= 0.5d0 * abs(b-c)

            con3 = (.not. mflag) .and. abs(s-b) >= 0.5d0 * abs(c-d)

            con4 = mflag .and. abs(b-c) < TOL

            con5 = (.not. mflag) .and. abs(c-d) < TOL

            if (con1 .or. con2 .or. con3 .or. con4 .or. con5) then
                s = 0.5d0 * (a + b)
                mflag = .true.
            else
                mflag = .false.
            end if

            call f_of_p(fs, s, U, Ncomp, gamma, gamma_up)

            if (abs(fa) < abs(fb)) then
                d = a
                a = b
                b = d

                d = fa
                fa = fb
                fb = d
            end if

            d = c
            c = b
            fc = fb

            if (fa * fs < 0.0d0) then
                b = s
                fb = fs
            else
                a = s
                fa = fs
            end if

            if (fb == 0.0d0 .or. fs == 0.0d0 .or. abs(b-a) < TOL) then
                p = b
                return
            end if

        end do

        p = x1

    end subroutine zbrent

end module gr_utils_module
