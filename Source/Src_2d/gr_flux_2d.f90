module compute_flux_module

  implicit none

  private

  public :: cons_to_prim, gr_swe_flux, f_of_p, gr_comp_flux

contains

  subroutine swe_flux(U, f, glo, ghi, lo, hi, Ncomp, g, x_dir)
      implicit none

      integer, intent(in) :: Ncomp
      double precision, intent(in) :: g
      integer, intent(in) :: glo(2), ghi(2), lo(2), hi(2)
      double precision, intent(in)  :: U(glo(1):ghi(1), glo(2):ghi(2), Ncomp)
      double precision, intent(out) :: f(glo(1):ghi(1), glo(2):ghi(2), Ncomp)
      integer, intent(in) :: x_dir

      integer i, j, k

      f(:,:,:) = 0.0d0

      if (x_dir == 0) then
          do    j = lo(2), hi(2)
             do i = lo(1)-1, hi(1)+1
                f(i,j,1) =  U(i,j,2)
                f(i,j,2) = U(i,j,2)**2/U(i,j,1) + 0.5d0 * g * U(i,j,1)**2
                f(i,j,3) =  U(i,j,2) * U(i,j,3) / U(i,j,1)
             end do
          end do
      else if (x_dir == 1) then
          do    j = lo(2)-1, hi(2)+1
             do i = lo(1), hi(1)
                f(i,j,1) =  U(i,j,3)
                f(i,j,2) = U(i,j,2) * U(i,j,3) / U(i,j,1)
                f(i,j,3) =  U(i,j,3)**2/U(i,j,1) + 0.5d0 * g * U(i,j,1)**2
             end do
          end do
      end if

  end subroutine swe_flux

  subroutine W_swe(q, lo, hi, Ncomp, gamma_up, glo, ghi, W)
      ! Calculate lorentz factor
      implicit none

      integer, intent(in) :: lo(2), hi(2), Ncomp, glo(2), ghi(2)
      double precision, intent(in) :: q(glo(1):ghi(1), glo(2):ghi(2), Ncomp)
      double precision, intent(in) :: gamma_up(glo(1):ghi(1), glo(2):ghi(2), 4)
      double precision, intent(out) :: W(glo(1):ghi(1), glo(2):ghi(2))

      integer i, j, k

      do j = lo(2), hi(2)
          do i = lo(1), hi(1)
              W(i,j) = sqrt((q(i,j,2)**2 * gamma_up(i,j,1)+ &
                    2.0d0 * q(i,j,2) * q(i,j,3) * gamma_up(i,j,2) + &
                    q(i,j,3)**2 * gamma_up(i,j,4)) / q(i,j,1)**2 + 1.0d0)
              ! nan check
              if (W(i,j) /= W(i,j)) then
                  W(i,j) = 1.0d0
              end if
        end do
    end do

  end subroutine W_swe

  subroutine gr_swe_flux(U, f, glo, ghi, lo, hi, Ncomp, x_dir, alpha)
      implicit none

      integer, intent(in) :: Ncomp
      integer, intent(in) :: glo(2), ghi(2), lo(2), hi(2)
      double precision, intent(in)  :: U(glo(1):ghi(1), glo(2):ghi(2), Ncomp)
      double precision, intent(out) :: f(glo(1):ghi(1), glo(2):ghi(2), Ncomp)
      integer, intent(in) :: x_dir
      double precision, intent(out) :: alpha(glo(1):ghi(1), glo(2):ghi(2))

      integer :: i, j, k
      double precision :: v(2), W(glo(1):ghi(1),glo(2):ghi(2))
      double precision :: beta(glo(1):ghi(1),glo(2):ghi(2),2)
      double precision :: gamma_up(glo(1):ghi(1),glo(2):ghi(2),4)

      alpha = exp(-U(:,:,1))
      beta = 0.0d0

      gamma_up = 0.0d0
      gamma_up(:,:,1) = 1.0d0
      gamma_up(:,:,4) = 1.0d0

      call W_swe(U, lo-1, hi+1, Ncomp, gamma_up, glo, ghi, W)

      f(:,:,:) = 0.0d0

      if (x_dir == 0) then
          do    j = lo(2), hi(2)
             do i = lo(1)-1, hi(1)+1
                if (U(i,j,1) < 1.d-20) then
                    v(1) = U(i,j,2)
                    v(2) = U(i,j,3)
                else
                    v(1) = U(i,j,2) / (U(i,j,1) * W(i,j))
                    v(2) = U(i,j,3) / (U(i,j,1) * W(i,j))
                end if

                f(i,j,1) =  U(i,j,1) * &
                      (v(1)*gamma_up(i,j,1) + v(2)*gamma_up(i,j,2) -&
                       beta(i,j,1) / alpha(i,j))
                f(i,j,2) = U(i,j,2) * &
                      (v(1)*gamma_up(i,j,1) + v(2)*gamma_up(i,j,2) -&
                       beta(i,j,1) / alpha(i,j)) + 0.5d0 * U(i,j,1)**2 / W(i,j)**2
                f(i,j,3) =  U(i,j,3) * &
                      (v(1)*gamma_up(i,j,1) + v(2)*gamma_up(i,j,2) -&
                       beta(i,j,1) / alpha(i,j))

                !f(i,j,:) = f(i,j,:)
             end do
          end do
      else if (x_dir == 1) then
          do    j = lo(2)-1, hi(2)+1
             do i = lo(1), hi(1)
                 if (U(i,j,1) < 1.d-20) then
                     v(1) = U(i,j,2)
                     v(2) = U(i,j,3)
                 else
                     v(1) = U(i,j,2) / (U(i,j,1) * W(i,j))
                     v(2) = U(i,j,3) / (U(i,j,1) * W(i,j))
                 end if

                f(i,j,1) = U(i,j,1) * &
                      (v(1)*gamma_up(i,j,2) + v(2)*gamma_up(i,j,4) -&
                       beta(i,j,2) / alpha(i,j))
                f(i,j,2) = U(i,j,2) * &
                      (v(1)*gamma_up(i,j,2) + v(2)*gamma_up(i,j,4) -&
                       beta(i,j,2) / alpha(i,j))
                f(i,j,3) =  U(i,j,3) * &
                      (v(1)*gamma_up(i,j,2) + v(2)*gamma_up(i,j,4) -&
                       beta(i,j,2) / alpha(i,j)) +  0.5d0 * U(i,j,1)**2 / W(i,j)**2

                !f(i,j,:) = f(i,j,:)
             end do
          end do
      end if

  end subroutine gr_swe_flux

  subroutine comp_flux(U, f, glo, ghi, lo, hi, Ncomp, gamma, x_dir)
      implicit none

      integer, intent(in) :: Ncomp
      double precision, intent(in) :: gamma
      integer, intent(in) :: glo(2), ghi(2), lo(2), hi(2)
      double precision, intent(in)  :: U(glo(1):ghi(1), glo(2):ghi(2), Ncomp)
      double precision, intent(out) :: f(glo(1):ghi(1), glo(2):ghi(2), Ncomp)
      integer, intent(in) :: x_dir

      integer i, j, k
      double precision :: p(glo(1):ghi(1), glo(2):ghi(2))

      p = (gamma - 1.d0) * (U(:,:,4) - 0.5d0 * (U(:,:,2)**2 + U(:,:,3)**2) / U(:,:,1))

      if (x_dir == 0) then
          do    j = lo(2), hi(2)
             do i = lo(1)-1, hi(1)+1
                f(i,j,1) =  U(i,j,2)
                f(i,j,2) = U(i,j,2)**2 / U(i,j,1) + p(i,j)
                f(i,j,3) =  U(i,j,2) * U(i,j,3) / U(i,j,1)
                f(i,j,4) = (U(i,j,4) + p(i,j)) * U(i,j,2) / U(i,j,1)
             end do
          end do
      else if (x_dir == 1) then
          do    j = lo(2)-1, hi(2)+1
             do i = lo(1), hi(1)
                f(i,j,1) =  U(i,j,3)
                f(i,j,2) = U(i,j,3) * U(i,j,2) / U(i,j,1)
                f(i,j,3) =  U(i,j,3)**2 / U(i,j,1) + p(i,j)
                f(i,j,4) = (U(i,j,4) + p(i,j)) * U(i,j,3) / U(i,j,1)
             end do
          end do
      end if

  end subroutine comp_flux

  subroutine f_of_p(f, p, U, Ncomp, gamma, gamma_up)
      implicit none

      integer, intent(in) :: Ncomp
      double precision, intent(in)  :: U(Ncomp), p, gamma, gamma_up(4)
      double precision, intent(out) :: f

      double precision :: sq

      sq = sqrt((U(4) + p + U(1))**2 - U(2)**2*gamma_up(1)- &
          2.0d0 * U(2) * U(3) * gamma_up(2) - &
          U(3)**2 * gamma_up(4))

      f = (gamma - 1.0d0) * sq / (U(4) + p + U(1)) * &
          (sq - p * (U(4) + p + U(1)) / sq - U(1)) - p

  end subroutine f_of_p

  subroutine cons_to_prim(U, clo, chi, U_prim, prlo, prhi, p, plo, &
          phi, lo, hi, Ncomp, gamma, alpha0, M, R, dx, prob_lo) &
          bind(C, name="cons_to_prim")
      use gr_utils_module, only : calc_gamma_up, calc_gamma_down, zbrent
      ! convert from conserved variables (D, Sx, Sy, tau) to primitive variables (rho, v^x, v^y, eps). Also outputs the pressure
      implicit none

      integer, intent(in) :: Ncomp
      integer, intent(in) :: clo(2), chi(2), prlo(2), prhi(2), plo(2), phi(2), lo(2), hi(2)
      double precision, intent(in)  :: U(clo(1):chi(1), clo(2):chi(2), Ncomp)
      double precision, intent(out) :: U_prim(prlo(1):prhi(1), prlo(2):prhi(2), Ncomp)
      double precision, intent(out) :: p(plo(1):phi(1), plo(2):phi(2))
      double precision, intent(in)  :: gamma, alpha0, M, R, dx(2), prob_lo(2)

      double precision gamma_up(lo(1):hi(1), lo(2):hi(2), 4)
      double precision gamma_down(lo(1):hi(1), lo(2):hi(2), 4)
      double precision :: pmin, pmax, ssq, q(Ncomp), fmin, fmax, sq, h, W2, v_up(2)
      integer :: i, j, k, l

      call calc_gamma_up(gamma_up, lo, hi, lo, hi, alpha0, M, R, dx, prob_lo)
      call calc_gamma_down(gamma_down, lo, hi, lo, hi, alpha0, M, R, dx, prob_lo)

      write(*,*) "cons_to_prim"

      do        j = lo(2), hi(2)
          do    i = lo(1), hi(1)
              q = U(i,j,:)

              !HACK
              !if (any(q /= q)) then
                  !write(*,*) "q is nan: ", q, i, j, k
                  !stop
              !end if
              if (q(1) /= q(1)) then !replace by average
                  q(1) = sum(U(:,:,1), U(:,:,1) == U(:,:,1)) / count(U(:,:,1) == U(:,:,1))
              end if
              do l = 2, 3
                  if (q(l) /= q(l)) then
                      q(l) = 0.0d0
                  end if
              end do
              if (q(4) == 1.0d0 .or. q(4) /= q(4)) then
                  q(4) = sum(U(:,:,4), U(:,:,4) == U(:,:,4)) / count(U(:,:,4) == U(:,:,4))
              end if
              ssq = q(2)**2 * gamma_up(i,j,1) + &
                  2.0d0 * q(2) * q(3) * gamma_up(i,j,2) + &
                  q(3)**2 * gamma_up(i,j,4)

              if (ssq /= ssq) then
                  ssq = 0.d0
              end if

              ! NOTE: divided pmax by 20 as it kept on finding the wrong root here
              pmin = (1.0d0 - ssq)**2 * q(4) * (gamma - 1.0d0)
              pmax = (gamma - 1.0d0) * (q(4) + q(1)) / (2.0d0 - gamma) / 20.0d0

              !if (i == lo(1)+1 .and. j == lo(2)+1 .and. k == lo(3)+1) then
                !write(*,*) "pmin: ", pmin, "pmax: ", pmax, "ssq: ", ssq, "D: ", q(1), "tau: ", q(5)
             ! end if

              if (pmin < 0.0d0) then
                  pmin = 0.d0
              end if

              if (pmax > 1.d0 .or. pmax < pmin) then
                  pmax = 1.0d0
              end if

              call f_of_p(fmin, pmin, q, Ncomp, gamma, gamma_up(i,j,:))
              call f_of_p(fmax, pmax, q, Ncomp, gamma, gamma_up(i,j,:))

              if (fmin * fmax > 0.0d0) then
                  pmin = 0.d0
              end if

              call f_of_p(fmin, pmin, q, Ncomp, gamma, gamma_up(i,j,:))

              if (fmin * fmax > 0.0d0) then
                  pmax = pmax * 10.d0
              end if

              call zbrent(p(i,j), pmin, pmax, q, Ncomp, gamma, gamma_up(i,j,:))

              if (p(i,j) /= p(i,j) .or. p(i,j) < 0.0d0 .or. p(i,j) > 1.0d0) then
                  p(i,j) = abs((gamma - 1.0d0) * (q(4) + q(1)) / (2.0d0 - gamma))

                  if (p(i,j) > 1.0d0) then
                      p(i,j) = 1.0d0
                  end if
              end if

              sq = sqrt((q(4) + p(i,j) + q(1))**2 - ssq)

              if (sq /= sq) then
                  sq = q(4) + p(i,j) + q(1)
              end if

              h = 1.0d0 + gamma * (sq - p(i,j) * (q(4) + p(i,j) + q(1)) / sq - q(1)) / q(1)
              W2 = 1.0d0 + ssq / (q(1) * h)**2

              !write(*,*) "p, sq", p(i,j), sq

              U_prim(i,j,1) = q(1) * sq / (q(4) + p(i,j) + q(1))
              v_up(1) = (gamma_up(i,j,1) * q(2) + &
                  gamma_up(i,j,2) * q(3)) /&
                  (W2 * h * U_prim(i,j,1))
              v_up(2) = (gamma_up(i,j,2) * q(2) + &
                  gamma_up(i,j,4) * q(3)) /&
                  (W2 * h * U_prim(i,j,1))

              U_prim(i,j,2) = gamma_down(i,j,1) * v_up(1) + gamma_down(i,j,2) * v_up(2)
              U_prim(i,j,3) = gamma_down(i,j,2) * v_up(1) + gamma_down(i,j,4) * v_up(2)
              U_prim(i,j,4) = (h - 1.0d0) / gamma

              if (U_prim(i,j,2) /= U_prim(i,j,2)) then
                  write(*,*) "p, q(1), sq, W, h", p(i,j), q(2), sq, W2, h
              end if

              ! HACK? fix p?
              p(i,j) = U_prim(i,j,1)**gamma

          end do
      end do

      !write(*,*) "U_comp", U(lo(1)+4, lo(2)+4, lo(3)+4, :)
      !write(*,*) "U_prim", U_prim(lo(1)+1, lo(2)+1, lo(3)+1, :)
      !write(*,*) "p: ", p(lo(1), lo(2), lo(3):hi(3))
      !write(*,*) "p", p(lo(1)+3, lo(2)+3, lo(3)+3)
      !write(*,*) "gamma_up", gamma_up(lo(1), lo(2), lo(3), 7:9)
      !write(*,*) "alpha0, R, M", alpha0, M, R

  end subroutine cons_to_prim

  subroutine gr_comp_flux(U, f, lo, hi, Ncomp, dir, &
      gamma, glo, ghi, alpha, dx, alpha0, M, R, prob_lo)
      implicit none

      integer, intent(in) :: Ncomp, dir
      integer, intent(in) :: lo(2), hi(2), glo(2), ghi(2)
      double precision, intent(inout)  :: U(glo(1):ghi(1), glo(2):ghi(2), Ncomp)
      double precision, intent(out) :: f(glo(1):ghi(1), glo(2):ghi(2), Ncomp)
      double precision, intent(in)  :: gamma, dx(2), alpha0, M, R, prob_lo(2)
      double precision, intent(out) :: alpha(glo(1):ghi(1), glo(2):ghi(2))

      integer :: i,j
      double precision :: p(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1)
      double precision :: U_prim(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, Ncomp)
      !double precision :: gamma_up(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), 9)
      double precision :: beta(glo(1):ghi(1), glo(2):ghi(2), 2)
      !double precision, parameter :: M = 1.0d0, R = 100.0d0

      beta = 0.0d0
      !gamma_up = 0.0d0
      !gamma_up(:,:,1) = 1.0d0
      !gamma_up(:,:,5) = 1.0d0
      !gamma_up(:,:,9) = alpha**2

      !alpha0 = sqrt(1.0d0 - 2.0d0 * M / R)

      write(*,*) "gr_comp_flux"
      write(*,*) "U: ", U(lo(1), lo(2), 5)

      ! HACK: while the boundaries don't work
      !do k = lo(3)-1, hi(3)+1
        !  do j = lo(2)-1, hi(2)+1
        !      do i = lo(1)-1, hi(1)+1
        !          if (U(i,j,1) == 0.0d0) then
        !              U(i,j,1) = 1.0d0
        !          end if
        !          if (U(i,j,5) == 0.0d0 .or. U(i,j,5) == 1.0d0) then
        !              U(i,j,5) = 1.5d0
        !          end if
        !      end do
         ! end do
      !end do


      call cons_to_prim(U, glo, ghi, U_prim, lo-1, hi+1, p, lo-1, hi+1, lo-1, hi+1, Ncomp, gamma, alpha0, M, R, dx, prob_lo)

      f(:,:,:) = 0.0d0

      if (dir == 0) then
          do     j = lo(2), hi(2)
              do i = lo(1)-1, hi(1)+1
                  f(i,j,1) = U(i,j,1) * (U_prim(i,j,2) - &
                       beta(i,j,1) / alpha(i,j))
                  f(i,j,2) = U(i,j,2) * (U_prim(i,j,2) - &
                       beta(i,j,1) / alpha(i,j)) + p(i,j)
                  f(i,j,3) = U(i,j,3) * (U_prim(i,j,2) - &
                       beta(i,j,1) / alpha(i,j))
                  f(i,j,4) = U(i,j,4) * (U_prim(i,j,2) - &
                       beta(i,j,1) / alpha(i,j)) + p(i,j) * U_prim(i,j,2)
               end do
           end do
      else if (dir == 1) then
          do     j = lo(2)-1, hi(2)+1
              do i = lo(1), hi(1)
                  !write(*,*) "y, alpha", alpha(lo(1),lo(2),lo(3))
                  f(i,j,1) = U(i,j,1) * (U_prim(i,j,3) - &
                       beta(i,j,2) / alpha(i,j))
                  f(i,j,2) = U(i,j,2) * (U_prim(i,j,3) - &
                      beta(i,j,2) / alpha(i,j))
                  f(i,j,3) = U(i,j,3) * (U_prim(i,j,3) - &
                       beta(i,j,2) / alpha(i,j)) + p(i,j)
                  f(i,j,4) = U(i,j,4) * (U_prim(i,j,3) - &
                       beta(i,j,2) / alpha(i,j)) + p(i,j) * U_prim(i,j,3)
               end do
           end do
      end if
  end subroutine gr_comp_flux

end module compute_flux_module
