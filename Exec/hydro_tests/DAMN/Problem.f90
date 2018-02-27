! problem-specific Fortran stuff goes here

subroutine problem_checkpoint(int_dir_name, len) bind(C, name="problem_checkpoint")

  ! called by the IO processor during checkpoint

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer :: len
  integer :: int_dir_name(len)
  character (len=len) :: dir

  integer :: i

  ! dir will be the string name of the checkpoint directory
  do i = 1, len
     dir(i:i) = char(int_dir_name(i))
  enddo



end subroutine problem_checkpoint


subroutine problem_restart(int_dir_name, len) bind(C, name="problem_restart")

  ! called by ALL processors during restart

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer :: len
  integer :: int_dir_name(len)
  character (len=len) :: dir

  integer :: i

  ! dir will be the string name of the checkpoint directory
  do i = 1, len
     dir(i:i) = char(int_dir_name(i))
  enddo

end subroutine problem_restart

subroutine smooth_initial_data(dat, dlo, dhi, lo, hi, delta, xlo) bind(C, name="smooth_initial_data")

    ! get rid of startup error
    use meth_params_module, only : NVAR
    use probdata_module, only: damn_rad

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) :: dlo(3), dhi(3), lo(3), hi(3)
    real(rt), intent(inout) :: dat(dlo(1):dhi(1), dlo(2):dhi(2), dlo(3):dhi(3), NVAR)
    real(rt), intent(in) :: delta(3), xlo(3)

    integer :: i,j,k,n
    real(rt) :: yy, avg

    integer :: smooth_zone

    ! radius of zone to smooth over
    smooth_zone = 2

    do k = lo(3), hi(3)
        do j = lo(2)+smooth_zone+1, hi(2)-smooth_zone-1
            yy = xlo(2) + delta(2)*dble(j-lo(2)+0.5_rt)

            if ((yy + delta(2)) > damn_rad) then
                ! smooth
                do i = lo(1), hi(1)
                    do n = 1, NVAR
                        avg = 0.5_rt * (dat(i,j-smooth_zone-1,k,n) + dat(i,j+smooth_zone+1,k,n))

                        dat(i,j-smooth_zone:j+smooth_zone,k,n) = avg
                    enddo
                enddo

                exit

            endif
        enddo
    enddo

end subroutine smooth_initial_data
