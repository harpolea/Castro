subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

    use bl_error_module
    use network
    use probdata_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer init, namlen
    integer name(namlen)
    real(rt)         problo(2), probhi(2)
    real(rt)         xn(nspec)

    integer untin,i

    namelist /fortin/ p_l, u_l, rho_l, p_r, u_r, rho_r, T_l, T_r, frac, idir, g

    !
    !     Build "probin" filename -- the name of file containing fortin namelist.
    !
    integer maxlen
    parameter (maxlen=256)
    character probin*(maxlen)

    if (namlen .gt. maxlen) then
        call bl_error("probin file name too long")
    end if

    do i = 1, namlen
        probin(i:i) = char(name(i))
    end do

    ! set namelist defaults

    p_l = 1.0d0               ! left pressure (erg/cc)
    u_l = 0.0d0               ! left velocity (cm/s)
    rho_l = 1.0d0             ! left density (g/cc)
    T_l = 1.0d0

    p_r = 0.1d0              ! right pressure (erg/cc)
    u_r = 0.0d0               ! right velocity (cm/s)
    rho_r = 0.125d0          ! right density (g/cc)
    T_r = 1.0d0

    idir = 1                ! direction across which to jump
    frac = 0.5              ! fraction of the domain for the interface

    g = 1.0d0

    !     Read namelists
    untin = 9
    open(untin,file=probin(1:namlen),form='formatted',status='old')
    read(untin,fortin)
    close(unit=untin)

    !     set local variable defaults -- the 'center' variables are the location of the
    !     interface
    split(1) = frac*(problo(1)+probhi(1))
    split(2) = frac*(problo(2)+probhi(2))

    !     compute the internal energy (erg/cc) for the left and right state
    xn(:) = 0.0e0_rt
    xn(1) = 1.0e0_rt


end subroutine amrex_probinit


! ::: -----------------------------------------------------------
! ::: This routine is called at problem setup time and is used
! ::: to initialize data on each grid.
! :::
! ::: NOTE:  all arrays have one cell of ghost zones surrounding
! :::        the grid interior.  Values in these cells need not
! :::        be set here.
! :::
! ::: INPUTS/OUTPUTS:
! :::
! ::: level     => amr level of grid
! ::: time      => time at which to init data
! ::: lo,hi     => index limits of grid interior (cell centered)
! ::: nstate    => number of state components.  You should know
! :::		   this already!
! ::: state     <=  Scalar array
! ::: delta     => cell size
! ::: xlo,xhi   => physical locations of lower left and upper
! :::              right hand corner of grid.  (does not include
! :::		   ghost region).
! ::: -----------------------------------------------------------
subroutine ca_initdata(level,time,lo,hi,nscal, &
                       state,state_l1,state_l2,state_h1,state_h2, &
                       delta,xlo,xhi)

  use network, only: nspec
  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UFS, QVAR, QRHO, QU, QV, QW
  use riemann_util_module, only : swe_cons_state

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer level, nscal
  integer lo(2), hi(2)
  integer state_l1,state_l2,state_h1,state_h2
  real(rt)         xlo(2), xhi(2), time, delta(2)
  real(rt)         state(state_l1:state_h1,state_l2:state_h2,NVAR)
  real(rt)         q(state_l1:state_h1,state_l2:state_h2,QVAR)

  real(rt)         xcen,ycen
  integer i,j

  q(:,:,QW) = 0.0d0

  do j = lo(2), hi(2)
     ycen = xlo(2) + delta(2)*(float(j-lo(2)) + 0.5e0_rt)

     do i = lo(1), hi(1)
        xcen = xlo(1) + delta(1)*(float(i-lo(1)) + 0.5e0_rt)

        if (idir == 1) then
           if (xcen <= split(1)) then
              q(i,j,QRHO) = rho_l
              q(i,j,QU) = rho_l*u_l
              q(i,j,QV) = 0.e0_rt
           else
              q(i,j,QRHO) = rho_r
              q(i,j,QU) = rho_r*u_r
              q(i,j,QV) = 0.e0_rt
           endif

        else if (idir == 2) then
           if (ycen <= split(2)) then
              q(i,j,QRHO) = rho_l
              q(i,j,QU) = 0.e0_rt
              q(i,j,QV) = rho_l*u_l
           else
              q(i,j,QRHO) = rho_r
              q(i,j,QU) = 0.e0_rt
              q(i,j,QV) = rho_r*u_r
           endif

        else
           call bl_abort('invalid idir')
        endif

        ! calculate conserved state from primitive initial data
        call swe_cons_state(q(i,j,:), state(i,j,:))

        state(i,j,UFS:UFS-1+nspec) = 0.0e0_rt
        state(i,j,UFS  ) = state(i,j,URHO)
     enddo
  enddo
end subroutine ca_initdata
