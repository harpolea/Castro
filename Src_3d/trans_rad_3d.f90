module transverse_rad_module

    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QFA, QFS, QFX, &
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, &
                                   UFA, UFS, UFX, &
                                   GDPRES, GDU, GDV, GDW, GDERADS, ngdnv, &
                                   nadv
    use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot, &
                                       fspace_type, comoving
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor

    implicit none

    private
    public :: transx1_rad, transx2_rad, transy1_rad, transy2_rad, transz_rad, &
              transxy_rad, transyz_rad, transxz_rad

contains  

  subroutine transy1_rad(lam, lam_lo, lam_hi, &
                         qxm, qxmo, qxp, qxpo, qd_lo, qd_hi, &
                         fy, rfy, fy_lo, fy_hi, &
                         qy, qy_lo, qy_hi, &
                         gamc, gd_lo, gd_hi, &
                         cdtdy, ilo, ihi, jlo, jhi, kc, k3d)

    integer :: lam_lo(3), lam_hi(3)
    integer :: qd_lo(3), qd_hi(3)
    integer :: fy_lo(3), fy_hi(3)
    integer :: qy_lo(3), qy_hi(3)
    integer :: gd_lo(3), gd_hi(3)
    integer ilo,ihi,jlo,jhi,kc,k3d

    double precision lam(lam_lo(1):lam_hi(1),lam_lo(2):lam_hi(2),lam_lo(3):lam_hi(3),0:ngroups-1)
    double precision  qxm(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision  qxp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision qxmo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision qxpo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision  fy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),NVAR)
    double precision rfy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),0:ngroups-1)
    double precision  qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),ngdnv)
    double precision gamc(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2),gd_lo(3):gd_hi(3))
    double precision cdtdy

    ! Local variables

    integer i, j, g
    integer n, nq
    integer iadv, ispec, iaux

    double precision rrnew, rr
    double precision compn, compu, compsn, comps
    double precision rrrx, rrlx
    double precision rurx, rulx
    double precision rvrx, rvlx
    double precision rwrx, rwlx
    double precision ekenrx, ekenlx
    double precision rerx, relx
    double precision rrnewrx, rrnewlx
    double precision runewrx, runewlx
    double precision rvnewrx, rvnewlx
    double precision rwnewrx, rwnewlx
    double precision renewrx, renewlx
    double precision pnewrx, pnewlx
    double precision rhoekenrx, rhoekenlx
    double precision ugp, ugm, dup, pav, du

    double precision :: pggp, pggm, dre, dmom
    double precision, dimension(0:ngroups-1) :: lambda, ergp, ergm, err, erl, ernewr, ernewl, &
         lamge, luge, der
    double precision eddf, f1, ugc

    do iadv = 1, nadv
       n = UFA + iadv - 1
       nq = QFA + iadv - 1

       do j = jlo, jhi 
          do i = ilo, ihi 

             compn = cdtdy*(fy(i,j+1,kc,n) - fy(i,j,kc,n))

             rr = qxp(i,j,kc,QRHO)
             rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
             compu = rr*qxp(i,j,kc,nq) - compn
             qxpo(i,j,kc,nq) = compu/rrnew

             rr = qxm(i+1,j,kc,QRHO)
             rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
             compu = rr*qxm(i+1,j,kc,nq) - compn
             qxmo(i+1,j,kc,nq) = compu/rrnew

          enddo
       enddo
    enddo

    do ispec = 1, nspec 
       n  = UFS + ispec - 1
       nq = QFS + ispec - 1

       do j = jlo, jhi 
          do i = ilo, ihi 

             compsn = cdtdy*(fy(i,j+1,kc,n) - fy(i,j,kc,n))

             rr = qxp(i,j,kc,QRHO)
             rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
             comps = rr*qxp(i,j,kc,nq) - compsn
             qxpo(i,j,kc,nq) = comps/rrnew

             rr = qxm(i+1,j,kc,QRHO)
             rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
             comps = rr*qxm(i+1,j,kc,nq) - compsn
             qxmo(i+1,j,kc,nq) = comps/rrnew

          enddo
       enddo
    enddo

    do iaux = 1, naux 
       n  = UFX + iaux - 1
       nq = QFX + iaux - 1

       do j = jlo, jhi 
          do i = ilo, ihi 

             compsn = cdtdy*(fy(i,j+1,kc,n) - fy(i,j,kc,n))

             rr = qxp(i,j,kc,QRHO)
             rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
             comps = rr*qxp(i,j,kc,nq) - compsn
             qxpo(i,j,kc,nq) = comps/rrnew

             rr = qxm(i+1,j,kc,QRHO)
             rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
             comps = rr*qxm(i+1,j,kc,nq) - compsn
             qxmo(i+1,j,kc,nq) = comps/rrnew

          enddo
       enddo
    enddo

    do j = jlo, jhi
       do i = ilo, ihi

          lambda = lam(i,j,k3d,:)
          pggp  = qy(i,j+1,kc,GDPRES)
          pggm  = qy(i,j  ,kc,GDPRES)
          ugp  =  qy(i,j+1,kc,GDV)
          ugm  =  qy(i,j  ,kc,GDV)
          ugc = 0.5d0*(ugp+ugm)
          ergp = qy(i,j+1,kc,GDERADS:GDERADS-1+ngroups)
          ergm = qy(i,j  ,kc,GDERADS:GDERADS-1+ngroups)

          ! Convert to conservation form
          rrrx = qxp(i,j,kc,QRHO)
          rurx = rrrx*qxp(i,j,kc,QU)
          rvrx = rrrx*qxp(i,j,kc,QV)
          rwrx = rrrx*qxp(i,j,kc,QW)
          ekenrx = 0.5d0*rrrx*(qxp(i,j,kc,QU)**2 + qxp(i,j,kc,QV)**2 &
               + qxp(i,j,kc,QW)**2)
          rerx = qxp(i,j,kc,QREINT) + ekenrx
          err  = qxp(i,j,kc,qrad:qradhi)

          rrlx = qxm(i+1,j,kc,QRHO)
          rulx = rrlx*qxm(i+1,j,kc,QU)
          rvlx = rrlx*qxm(i+1,j,kc,QV)
          rwlx = rrlx*qxm(i+1,j,kc,QW)
          ekenlx = 0.5d0*rrlx*(qxm(i+1,j,kc,QU)**2 + qxm(i+1,j,kc,QV)**2 &
               + qxm(i+1,j,kc,QW)**2)
          relx = qxm(i+1,j,kc,QREINT) + ekenlx
          erl  = qxm(i+1,j,kc,qrad:qradhi)

          ! Add transverse predictor
          rrnewrx = rrrx - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
          runewrx = rurx - cdtdy*(fy(i,j+1,kc,UMX) - fy(i,j,kc,UMX))
          rvnewrx = rvrx - cdtdy*(fy(i,j+1,kc,UMY) - fy(i,j,kc,UMY))
          rwnewrx = rwrx - cdtdy*(fy(i,j+1,kc,UMZ) - fy(i,j,kc,UMZ))
          lamge = lambda(:) * (ergp(:)-ergm(:))
          dmom = - cdtdy*sum(lamge(:))
          rvnewrx = rvnewrx + dmom
          luge = ugc * lamge(:)
          dre = -cdtdy*sum(luge)
          renewrx = rerx - cdtdy*(fy(i,j+1,kc,UEDEN) - fy(i,j,kc,UEDEN))  &
               + dre

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = cdtdy * ugc * f1 * (ergp(g) - ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = cdtdy * f1 * 0.5d0*(ergp(g)+ergm(g)) * (ugm-ugp)
             end do
          else ! mixed frame
             der(:) = cdtdy * luge
          end if

          ernewr = err(:) - cdtdy*(rfy(i,j+1,kc,:) - rfy(i,j,kc,:)) &
               + der(:)

          rrnewlx = rrlx - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
          runewlx = rulx - cdtdy*(fy(i,j+1,kc,UMX) - fy(i,j,kc,UMX))
          rvnewlx = rvlx - cdtdy*(fy(i,j+1,kc,UMY) - fy(i,j,kc,UMY))
          rwnewlx = rwlx - cdtdy*(fy(i,j+1,kc,UMZ) - fy(i,j,kc,UMZ))
          rvnewlx = rvnewlx + dmom
          renewlx = relx - cdtdy*(fy(i,j+1,kc,UEDEN)- fy(i,j,kc,UEDEN)) &
               + dre
          ernewl  = erl(:) - cdtdy*(rfy(i,j+1,kc,:) - rfy(i,j,kc,:)) &
               + der(:)

          dup = pggp*ugp - pggm*ugm
          pav = 0.5d0*(pggp+pggm)
          du = ugp-ugm

          pnewrx = qxp(i,j,kc,QPRES) - cdtdy*(dup + pav*du*(gamc(i,j,k3d) - 1.d0))
          pnewlx = qxm(i+1,j,kc,QPRES) - cdtdy*(dup + pav*du*(gamc(i,j,k3d) - 1.d0))

          ! Convert back to non-conservation form
          qxpo(i,j,kc,QRHO) = rrnewrx
          qxpo(i,j,kc,QU) = runewrx/qxpo(i,j,kc,QRHO)
          qxpo(i,j,kc,QV) = rvnewrx/qxpo(i,j,kc,QRHO)
          qxpo(i,j,kc,QW) = rwnewrx/qxpo(i,j,kc,QRHO)
          rhoekenrx = 0.5d0*(runewrx**2 + rvnewrx**2 + rwnewrx**2)/qxpo(i,j,kc,QRHO)
          qxpo(i,j,kc,QREINT)= renewrx - rhoekenrx
          qxpo(i,j,kc,QPRES) = pnewrx
          qxpo(i,j,kc,qrad:qradhi) = ernewr(:)
          qxpo(i,j,kc,qptot  ) = sum(lambda(:)*ernewr(:)) + qxpo(i,j,kc,QPRES)
          qxpo(i,j,kc,qreitot) = sum(qxpo(i,j,kc,qrad:qradhi)) + qxpo(i,j,kc,QREINT)

          qxmo(i+1,j,kc,QRHO) = rrnewlx
          qxmo(i+1,j,kc,QU) = runewlx/qxmo(i+1,j,kc,QRHO)
          qxmo(i+1,j,kc,QV) = rvnewlx/qxmo(i+1,j,kc,QRHO)
          qxmo(i+1,j,kc,QW) = rwnewlx/qxmo(i+1,j,kc,QRHO)
          rhoekenlx = 0.5d0*(runewlx**2 + rvnewlx**2 + rwnewlx**2)/qxmo(i+1,j,kc,QRHO)
          qxmo(i+1,j,kc,QREINT)= renewlx - rhoekenlx
          qxmo(i+1,j,kc,QPRES) = pnewlx
          qxmo(i+1,j,kc,qrad:qradhi) = ernewl(:)
          qxmo(i+1,j,kc,qptot  ) = sum(lambda(:)*ernewl(:)) + qxmo(i+1,j,kc,QPRES)
          qxmo(i+1,j,kc,qreitot) = sum(qxmo(i+1,j,kc,qrad:qradhi)) + qxmo(i+1,j,kc,QREINT)

       enddo
    enddo

  end subroutine transy1_rad

  ! ::: 
  ! ::: ------------------------------------------------------------------
  ! ::: 

  subroutine transx1_rad(lam, lam_lo, lam_hi, &
                         qym, qymo, qyp, qypo, qd_lo, qd_hi, &
                         fx, rfx, fx_lo, fx_hi, &
                         qx, qx_lo, qx_hi, &
                         gamc, gd_lo, gd_hi, &
                         cdtdx, ilo, ihi, jlo, jhi, kc, k3d)

    integer :: lam_lo(3), lam_hi(3)
    integer :: qd_lo(3), qd_hi(3)
    integer :: fx_lo(3), fx_hi(3)
    integer :: qx_lo(3), qx_hi(3)
    integer :: gd_lo(3), gd_hi(3)
    integer ilo,ihi,jlo,jhi,kc,k3d

    double precision lam(lam_lo(1):lam_hi(1),lam_lo(2):lam_hi(2),lam_lo(3):lam_hi(3),0:ngroups-1)
    double precision  qym(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision  qyp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision qymo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision qypo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision  fx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),NVAR)
    double precision rfx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),0:ngroups-1)
    double precision  qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),ngdnv)
    double precision gamc(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2),gd_lo(3):gd_hi(3))
    double precision cdtdx

    ! Local variables

    integer i, j, g
    integer n, nq
    integer iadv, ispec, iaux

    double precision rrnew, rr
    double precision rrry, rrly
    double precision rury, ruly
    double precision rvry, rvly
    double precision rwry, rwly
    double precision ekenry, ekenly
    double precision rery, rely
    double precision rrnewry, rrnewly
    double precision runewry, runewly
    double precision rvnewry, rvnewly
    double precision rwnewry, rwnewly
    double precision renewry, renewly
    double precision pnewry, pnewly
    double precision rhoekenry, rhoekenly
    double precision compn, compu, compsn, comps
    double precision ugp, ugm, dup, pav, du

    double precision :: pggp, pggm, dre, dmom
    double precision, dimension(0:ngroups-1) :: lambda, ergp, ergm, err, erl, ernewr, ernewl, &
         lamge, luge, der
    double precision eddf, f1, ugc

    ! NOTE: it is better *not* to protect against small density in this routine

    do iadv = 1, nadv
       n = UFA + iadv - 1
       nq = QFA + iadv - 1
       do j = jlo, jhi 
          do i = ilo, ihi 

             compn = cdtdx*(fx(i+1,j,kc,n) - fx(i,j,kc,n))

             rr = qyp(i,j,kc,QRHO)
             rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
             compu = rr*qyp(i,j,kc,nq) - compn
             qypo(i,j,kc,nq) = compu/rrnew

             rr = qym(i,j+1,kc,QRHO)
             rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
             compu = rr*qym(i,j+1,kc,nq) - compn
             qymo(i,j+1,kc,nq) = compu/rrnew

          enddo
       enddo
    enddo

    do ispec = 1, nspec
       n  = UFS + ispec - 1
       nq = QFS + ispec - 1
       do j = jlo, jhi 
          do i = ilo, ihi 

             compsn = cdtdx*(fx(i+1,j,kc,n) - fx(i,j,kc,n))

             rr = qyp(i,j,kc,QRHO)
             rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
             comps = rr*qyp(i,j,kc,nq) - compsn
             qypo(i,j,kc,nq) = comps/rrnew

             rr = qym(i,j+1,kc,QRHO)
             rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
             comps = rr*qym(i,j+1,kc,nq) - compsn
             qymo(i,j+1,kc,nq) = comps/rrnew

          enddo
       enddo
    enddo

    do iaux = 1, naux
       n  = UFX + iaux - 1
       nq = QFX + iaux - 1
       do j = jlo, jhi 
          do i = ilo, ihi 

             compsn = cdtdx*(fx(i+1,j,kc,n) - fx(i,j,kc,n))

             rr = qyp(i,j,kc,QRHO)
             rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
             comps = rr*qyp(i,j,kc,nq) - compsn
             qypo(i,j,kc,nq) = comps/rrnew

             rr = qym(i,j+1,kc,QRHO)
             rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
             comps = rr*qym(i,j+1,kc,nq) - compsn
             qymo(i,j+1,kc,nq) = comps/rrnew

          enddo
       enddo
    enddo

    do j = jlo, jhi 
       do i = ilo, ihi 

          lambda = lam(i,j,k3d,:)
          pggp  =  qx(i+1,j,kc,GDPRES)
          pggm  =  qx(i  ,j,kc,GDPRES)
          ugp  =  qx(i+1,j,kc,GDU)
          ugm  =  qx(i  ,j,kc,GDU)
          ugc = 0.5d0*(ugp+ugm)
          ergp = qx(i+1,j,kc,GDERADS:GDERADS-1+ngroups)
          ergm = qx(i  ,j,kc,GDERADS:GDERADS-1+ngroups)

          ! Convert to conservation form
          rrry = qyp(i,j,kc,QRHO)
          rury = rrry*qyp(i,j,kc,QU)
          rvry = rrry*qyp(i,j,kc,QV)
          rwry = rrry*qyp(i,j,kc,QW)
          ekenry = 0.5d0*rrry*(qyp(i,j,kc,QU)**2 + qyp(i,j,kc,QV)**2 + qyp(i,j,kc,QW)**2)
          rery = qyp(i,j,kc,QREINT) + ekenry
          err  = qyp(i,j,kc,qrad:qradhi)

          rrly = qym(i,j+1,kc,QRHO)
          ruly = rrly*qym(i,j+1,kc,QU)
          rvly = rrly*qym(i,j+1,kc,QV)
          rwly = rrly*qym(i,j+1,kc,QW)
          ekenly = 0.5d0*rrly* &
               (qym(i,j+1,kc,QU)**2 + qym(i,j+1,kc,QV)**2 + qym(i,j+1,kc,QW)**2)
          rely = qym(i,j+1,kc,QREINT) + ekenly
          erl  = qym(i,j+1,kc,qrad:qradhi)

          ! Add transverse predictor
          rrnewry = rrry - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
          runewry = rury - cdtdx*(fx(i+1,j,kc,UMX) - fx(i,j,kc,UMX))
          rvnewry = rvry - cdtdx*(fx(i+1,j,kc,UMY) - fx(i,j,kc,UMY))
          rwnewry = rwry - cdtdx*(fx(i+1,j,kc,UMZ) - fx(i,j,kc,UMZ))
          lamge = lambda(:) * (ergp(:)-ergm(:))
          dmom = - cdtdx*sum(lamge(:))
          runewry = runewry + dmom
          luge = ugc * lamge(:)
          dre = -cdtdx*sum(luge)
          renewry = rery - cdtdx*(fx(i+1,j,kc,UEDEN) - fx(i,j,kc,UEDEN)) &
               + dre

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = cdtdx * ugc * f1 * (ergp(g) - ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = cdtdx * f1 * 0.5d0*(ergp(g)+ergm(g)) * (ugm-ugp)
             end do
          else ! mixed frame
             der(:) = cdtdx * luge
          end if

          ernewr  = err(:) - cdtdx*(rfx(i+1,j,kc,:) - rfx(i,j,kc,:)) &
               + der(:)

          rrnewly = rrly - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
          runewly = ruly - cdtdx*(fx(i+1,j,kc,UMX) - fx(i,j,kc,UMX))
          rvnewly = rvly - cdtdx*(fx(i+1,j,kc,UMY) - fx(i,j,kc,UMY))
          rwnewly = rwly - cdtdx*(fx(i+1,j,kc,UMZ) - fx(i,j,kc,UMZ))
          runewly = runewly + dmom
          renewly = rely - cdtdx*(fx(i+1,j,kc,UEDEN) - fx(i,j,kc,UEDEN)) &
               + dre
          ernewl  = erl(:) - cdtdx*(rfx(i+1,j,kc,:) - rfx(i,j,kc,:)) &
               + der(:)

          dup = pggp*ugp - pggm*ugm
          pav = 0.5d0*(pggp+pggm)
          du = ugp-ugm

          pnewry = qyp(i,j,kc,QPRES) - cdtdx*(dup + pav*du*(gamc(i,j,k3d)-1.d0))
          pnewly = qym(i,j+1,kc,QPRES) - cdtdx*(dup + pav*du*(gamc(i,j,k3d)-1.d0))

          ! Convert back to non-conservation form
          qypo(i,j,kc,QRHO) = rrnewry
          qypo(i,j,kc,QU) = runewry/qypo(i,j,kc,QRHO)
          qypo(i,j,kc,QV) = rvnewry/qypo(i,j,kc,QRHO)
          qypo(i,j,kc,QW) = rwnewry/qypo(i,j,kc,QRHO)
          rhoekenry = 0.5d0*(runewry**2 + rvnewry**2 + rwnewry**2)/qypo(i,j,kc,QRHO)
          qypo(i,j,kc,QREINT) = renewry - rhoekenry
          qypo(i,j,kc,QPRES) = pnewry
          qypo(i,j,kc,qrad:qradhi) = ernewr(:)
          qypo(i,j,kc,qptot  ) = sum(lambda(:)*ernewr(:)) + qypo(i,j,kc,QPRES)
          qypo(i,j,kc,qreitot) = sum(qypo(i,j,kc,qrad:qradhi)) + qypo(i,j,kc,QREINT)

          qymo(i,j+1,kc,QRHO) = rrnewly
          qymo(i,j+1,kc,QU) = runewly/qymo(i,j+1,kc,QRHO)
          qymo(i,j+1,kc,QV) = rvnewly/qymo(i,j+1,kc,QRHO)
          qymo(i,j+1,kc,QW) = rwnewly/qymo(i,j+1,kc,QRHO)
          rhoekenly = 0.5d0*(runewly**2 + rvnewly**2 + rwnewly**2)/qymo(i,j+1,kc,QRHO)
          qymo(i,j+1,kc,QREINT) = renewly - rhoekenly
          qymo(i,j+1,kc,QPRES) = pnewly
          qymo(i,j+1,kc,qrad:qradhi) = ernewl(:)
          qymo(i,j+1,kc,qptot  ) = sum(lambda(:)*ernewl(:)) + qymo(i,j+1,kc,QPRES)
          qymo(i,j+1,kc,qreitot) = sum(qymo(i,j+1,kc,qrad:qradhi)) + qymo(i,j+1,kc,QREINT)

       enddo
    enddo

  end subroutine transx1_rad

  ! ::: 
  ! ::: ------------------------------------------------------------------
  ! ::: 

  subroutine transy2_rad(lam, lam_lo, lam_hi, &
                         qzm, qzmo, qzp, qzpo,qd_lo, qd_hi, &
                         fy, rfy, fy_lo, fy_hi, &
                         qy, qy_lo, qy_hi, &
                         gamc, gd_lo, gd_hi, &
                         cdtdy, ilo, ihi, jlo, jhi, kc, km, k3d)

    integer :: lam_lo(3), lam_hi(3)
    integer :: qd_lo(3), qd_hi(3)
    integer :: fy_lo(3), fy_hi(3)
    integer :: qy_lo(3), qy_hi(3)
    integer :: gd_lo(3), gd_hi(3)
    integer ilo, ihi, jlo, jhi, kc, km, k3d

    double precision lam(lam_lo(1):lam_hi(1),lam_lo(2):lam_hi(2),lam_lo(3):lam_hi(3),0:ngroups-1)
    double precision  qzm(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision  qzp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision qzmo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision qzpo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision  fy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),NVAR)
    double precision rfy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),0:ngroups-1)
    double precision  qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),ngdnv)
    double precision gamc(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2),gd_lo(3):gd_hi(3))
    double precision cdtdy

    ! Local variables

    integer i, j, g
    integer n, nq
    integer iadv, ispec, iaux

    double precision rrnew, rr
    double precision compn, compu, compsn, comps
    double precision rrrz, rrlz
    double precision rurz, rulz
    double precision rvrz, rvlz
    double precision rwrz, rwlz
    double precision ekenrz, ekenlz
    double precision rerz, relz
    double precision rrnewrz, rrnewlz
    double precision runewrz, runewlz
    double precision rvnewrz, rvnewlz
    double precision rwnewrz, rwnewlz
    double precision renewrz, renewlz
    double precision pnewrz, pnewlz
    double precision rhoekenrz, rhoekenlz
    double precision ugp, ugm, dup, pav, du

    double precision :: pggp, pggm, dre, dmom
    double precision, dimension(0:ngroups-1) :: lambda, ergp, ergm, err, erl, ernewr, ernewl, &
         lamge, luge, der
    double precision eddf, f1, ugc

    do iadv = 1, nadv
       n  = UFA + iadv - 1
       nq = QFA + iadv - 1

       do j = jlo, jhi 
          do i = ilo, ihi 

             compn = cdtdy*(fy(i,j+1,kc,n) - fy(i,j,kc,n))

             rr = qzp(i,j,kc,QRHO)
             rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
             compu = rr*qzp(i,j,kc,nq) - compn
             qzpo(i,j,kc,nq) = compu/rrnew

             compn = cdtdy*(fy(i,j+1,km,n) - fy(i,j,km,n))

             rr = qzm(i,j,kc,QRHO)
             rrnew = rr - cdtdy*(fy(i,j+1,km,URHO) - fy(i,j,km,URHO))
             compu = rr*qzm(i,j,kc,nq) - compn
             qzmo(i,j,kc,nq) = compu/rrnew

          enddo
       enddo
    enddo

    do ispec = 1, nspec 
       n  = UFS + ispec - 1
       nq = QFS + ispec - 1

       do j = jlo, jhi 
          do i = ilo, ihi 

             compsn = cdtdy*(fy(i,j+1,kc,n) - fy(i,j,kc,n))

             rr = qzp(i,j,kc,QRHO)
             rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
             comps = rr*qzp(i,j,kc,nq) - compsn
             qzpo(i,j,kc,nq) = comps/rrnew

             compsn = cdtdy*(fy(i,j+1,km,n) - fy(i,j,km,n))

             rr = qzm(i,j,kc,QRHO)
             rrnew = rr - cdtdy*(fy(i,j+1,km,URHO) - fy(i,j,km,URHO))
             comps = rr*qzm(i,j,kc,nq) - compsn
             qzmo(i,j,kc,nq) = comps/rrnew

          enddo
       enddo
    enddo

    do iaux = 1, naux 
       n  = UFX + iaux - 1
       nq = QFX + iaux - 1

       do j = jlo, jhi 
          do i = ilo, ihi 

             compsn = cdtdy*(fy(i,j+1,kc,n) - fy(i,j,kc,n))

             rr = qzp(i,j,kc,QRHO)
             rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
             comps = rr*qzp(i,j,kc,nq) - compsn
             qzpo(i,j,kc,nq) = comps/rrnew

             compsn = cdtdy*(fy(i,j+1,km,n) - fy(i,j,km,n))

             rr = qzm(i,j,kc,QRHO)
             rrnew = rr - cdtdy*(fy(i,j+1,km,URHO) - fy(i,j,km,URHO))
             comps = rr*qzm(i,j,kc,nq) - compsn
             qzmo(i,j,kc,nq) = comps/rrnew

          enddo
       enddo
    enddo

    do j = jlo, jhi
       do i = ilo, ihi

          lambda = lam(i,j,k3d,:)
          pggp  = qy(i,j+1,kc,GDPRES)
          pggm  = qy(i,j  ,kc,GDPRES)
          ugp  =  qy(i,j+1,kc,GDV)
          ugm  =  qy(i,j  ,kc,GDV)
          ugc = 0.5d0*(ugp+ugm)
          ergp = qy(i,j+1,kc,GDERADS:GDERADS-1+ngroups)
          ergm = qy(i,j  ,kc,GDERADS:GDERADS-1+ngroups)

          ! Convert to conservation form
          rrrz = qzp(i,j,kc,QRHO)
          rurz = rrrz*qzp(i,j,kc,QU)
          rvrz = rrrz*qzp(i,j,kc,QV)
          rwrz = rrrz*qzp(i,j,kc,QW)
          ekenrz = 0.5d0*rrrz*(qzp(i,j,kc,QU)**2 + qzp(i,j,kc,QV)**2 &
               + qzp(i,j,kc,QW)**2)
          rerz = qzp(i,j,kc,QREINT) + ekenrz
          err  = qzp(i,j,kc,qrad:qradhi)

          ! Add transverse predictor
          rrnewrz = rrrz - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
          runewrz = rurz - cdtdy*(fy(i,j+1,kc,UMX) - fy(i,j,kc,UMX))
          rvnewrz = rvrz - cdtdy*(fy(i,j+1,kc,UMY) - fy(i,j,kc,UMY))
          rwnewrz = rwrz - cdtdy*(fy(i,j+1,kc,UMZ) - fy(i,j,kc,UMZ))
          lamge = lambda(:) * (ergp(:)-ergm(:))
          dmom = - cdtdy*sum(lamge(:))
          rvnewrz = rvnewrz + dmom
          luge = 0.5d0*(ugp+ugm) * lamge(:)
          dre = -cdtdy*sum(luge)
          renewrz = rerz - cdtdy*(fy(i,j+1,kc,UEDEN) - fy(i,j,kc,UEDEN)) &
               + dre

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = cdtdy * ugc * f1 * (ergp(g) - ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = cdtdy * f1 * 0.5d0*(ergp(g)+ergm(g)) * (ugm-ugp)
             end do
          else ! mixed frame
             der(:) = cdtdy * luge
          end if

          ernewr  = err(:) - cdtdy*(rfy(i,j+1,kc,:) - rfy(i,j,kc,:)) &
               + der(:)

          dup = pggp*ugp - pggm*ugm
          pav = 0.5d0*(pggp+pggm)
          du = ugp-ugm

          pnewrz = qzp(i,j,kc,QPRES) - cdtdy*(dup + pav*du*(gamc(i,j,k3d) - 1.d0))

          ! Convert back to non-conservation form
          qzpo(i,j,kc,QRHO) = rrnewrz
          qzpo(i,j,kc,QU) = runewrz/qzpo(i,j,kc,QRHO)
          qzpo(i,j,kc,QV) = rvnewrz/qzpo(i,j,kc,QRHO)
          qzpo(i,j,kc,QW) = rwnewrz/qzpo(i,j,kc,QRHO)
          rhoekenrz = 0.5d0*(runewrz**2 + rvnewrz**2 + rwnewrz**2)/qzpo(i,j,kc,QRHO)
          qzpo(i,j,kc,QREINT)= renewrz - rhoekenrz
          qzpo(i,j,kc,QPRES) = pnewrz
          qzpo(i,j,kc,qrad:qradhi) = ernewr(:)
          qzpo(i,j,kc,qptot  ) = sum(lambda(:)*ernewr(:)) + qzpo(i,j,kc,QPRES)
          qzpo(i,j,kc,qreitot) = sum(qzpo(i,j,kc,qrad:qradhi)) + qzpo(i,j,kc,QREINT)

          lambda = lam(i,j,k3d-1,:)
          pggp  =  qy(i,j+1,km,GDPRES)
          pggm  =  qy(i,j  ,km,GDPRES)
          ugp  =  qy(i,j+1,km,GDV)
          ugm  =  qy(i,j  ,km,GDV)
          ugc = 0.5d0*(ugp+ugm)
          ergp = qy(i,j+1,km,GDERADS:GDERADS-1+ngroups)
          ergm = qy(i,j  ,km,GDERADS:GDERADS-1+ngroups)

          rrlz = qzm(i,j,kc,QRHO)
          rulz = rrlz*qzm(i,j,kc,QU)
          rvlz = rrlz*qzm(i,j,kc,QV)
          rwlz = rrlz*qzm(i,j,kc,QW)
          ekenlz = 0.5d0*rrlz*(qzm(i,j,kc,QU)**2 + qzm(i,j,kc,QV)**2 &
               + qzm(i,j,kc,QW)**2)
          relz = qzm(i,j,kc,QREINT) + ekenlz
          erl  = qzm(i,j,kc,qrad:qradhi)

          ! Add transverse predictor
          rrnewlz = rrlz - cdtdy*(fy(i,j+1,km,URHO) - fy(i,j,km,URHO))
          runewlz = rulz - cdtdy*(fy(i,j+1,km,UMX) - fy(i,j,km,UMX))
          rvnewlz = rvlz - cdtdy*(fy(i,j+1,km,UMY) - fy(i,j,km,UMY))
          rwnewlz = rwlz - cdtdy*(fy(i,j+1,km,UMZ) - fy(i,j,km,UMZ))
          lamge = lambda(:) * (ergp(:)-ergm(:))
          dmom = - cdtdy*sum(lamge(:))
          rvnewlz = rvnewlz + dmom
          luge = 0.5d0*(ugp+ugm) * lamge(:)
          dre = -cdtdy*sum(luge)
          renewlz = relz - cdtdy*(fy(i,j+1,km,UEDEN)- fy(i,j,km,UEDEN)) &
               + dre

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = cdtdy * ugc * f1 * (ergp(g) - ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = cdtdy * f1 * 0.5d0*(ergp(g)+ergm(g)) * (ugm-ugp)
             end do
          else ! mixed frame
             der(:) = cdtdy * luge
          end if

          ernewl  = erl(:) - cdtdy*(rfy(i,j+1,km,:)- rfy(i,j,km,:)) &
               + der

          dup = pggp*ugp - pggm*ugm
          pav = 0.5d0*(pggp+pggm)
          du = ugp-ugm

          pnewlz = qzm(i,j,kc,QPRES) - cdtdy*(dup + pav*du*(gamc(i,j,k3d-1) - 1.d0))

          qzmo(i,j,kc,QRHO) = rrnewlz
          qzmo(i,j,kc,QU) = runewlz/qzmo(i,j,kc,QRHO)
          qzmo(i,j,kc,QV) = rvnewlz/qzmo(i,j,kc,QRHO)
          qzmo(i,j,kc,QW) = rwnewlz/qzmo(i,j,kc,QRHO)
          rhoekenlz = 0.5d0*(runewlz**2 + rvnewlz**2 + rwnewlz**2)/qzmo(i,j,kc,QRHO)
          qzmo(i,j,kc,QREINT)= renewlz - rhoekenlz
          qzmo(i,j,kc,QPRES) = pnewlz
          qzmo(i,j,kc,qrad:qradhi) = ernewl(:)
          qzmo(i,j,kc,qptot  ) = sum(lambda(:)*ernewl(:)) + qzmo(i,j,kc,QPRES)
          qzmo(i,j,kc,qreitot) = sum(qzmo(i,j,kc,qrad:qradhi)) + qzmo(i,j,kc,QREINT)

       enddo
    enddo

  end subroutine transy2_rad

  ! ::: 
  ! ::: ------------------------------------------------------------------
  ! ::: 

  subroutine transx2_rad(lam, lam_lo, lam_hi, & 
                         qzm, qzmo, qzp, qzpo, qd_lo, qd_hi, &
                         fx, rfx, fx_lo, fx_hi, &
                         qx, qx_lo, qx_hi, &
                         gamc, gd_lo, gd_hi, &
                         cdtdx, ilo, ihi, jlo, jhi, kc, km, k3d)

    integer :: lam_lo(3), lam_hi(3)
    integer :: qd_lo(3), qd_hi(3)
    integer :: fx_lo(3), fx_hi(3)
    integer :: qx_lo(3), qx_hi(3)
    integer :: gd_lo(3), gd_hi(3)
    integer ilo,ihi,jlo,jhi,kc,km,k3d

    double precision lam(lam_lo(1):lam_hi(1),lam_lo(2):lam_hi(2),lam_lo(3):lam_hi(3),0:ngroups-1)
    double precision  qzm(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision  qzp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision qzmo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision qzpo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision  fx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),NVAR)
    double precision rfx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),0:ngroups-1)
    double precision  qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),ngdnv)
    double precision gamc(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2),gd_lo(3):gd_hi(3))
    double precision cdtdx

    ! Local variables

    integer i, j, g
    integer n, nq
    integer iadv, ispec, iaux

    double precision rrnew, rr
    double precision rrrz, rrlz
    double precision rurz, rulz
    double precision rvrz, rvlz
    double precision rwrz, rwlz
    double precision ekenrz, ekenlz
    double precision rerz, relz
    double precision rrnewrz, rrnewlz
    double precision runewrz, runewlz
    double precision rvnewrz, rvnewlz
    double precision rwnewrz, rwnewlz
    double precision renewrz, renewlz
    double precision pnewrz, pnewlz
    double precision rhoekenrz, rhoekenlz
    double precision compn, compu, compsn, comps
    double precision ugp, ugm, dup, pav, du

    double precision :: pggp, pggm, dre, dmom
    double precision, dimension(0:ngroups-1) :: lambda, ergp, ergm, err, erl, ernewr, ernewl, &
         lamge, luge, der
    double precision eddf, f1, ugc

    do iadv = 1, nadv
       n = UFA + iadv - 1
       nq = QFA + iadv - 1
       do j = jlo, jhi 
          do i = ilo, ihi 

             compn = cdtdx*(fx(i+1,j,kc,n) - fx(i,j,kc,n))

             rr = qzp(i,j,kc,QRHO)
             rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
             compu = rr*qzp(i,j,kc,nq) - compn
             qzpo(i,j,kc,nq) = compu/rrnew

             compn = cdtdx*(fx(i+1,j,km,n) - fx(i,j,km,n))

             rr = qzm(i,j,kc,QRHO)
             rrnew = rr - cdtdx*(fx(i+1,j,km,URHO) - fx(i,j,km,URHO))
             compu = rr*qzm(i,j,kc,nq) - compn
             qzmo(i,j,kc,nq) = compu/rrnew

          enddo
       enddo
    enddo

    do ispec = 1, nspec
       n  = UFS + ispec - 1
       nq = QFS + ispec - 1
       do j = jlo, jhi 
          do i = ilo, ihi 

             compsn = cdtdx*(fx(i+1,j,kc,n) - fx(i,j,kc,n))

             rr = qzp(i,j,kc,QRHO)
             rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
             comps = rr*qzp(i,j,kc,nq) - compsn
             qzpo(i,j,kc,nq) = comps/rrnew

             compsn = cdtdx*(fx(i+1,j,km,n) - fx(i,j,km,n))

             rr = qzm(i,j,kc,QRHO)
             rrnew = rr - cdtdx*(fx(i+1,j,km,URHO) - fx(i,j,km,URHO))
             comps = rr*qzm(i,j,kc,nq) - compsn
             qzmo(i,j,kc,nq) = comps/rrnew

          enddo
       enddo
    enddo

    do iaux = 1, naux
       n  = UFX + iaux - 1
       nq = QFX + iaux - 1
       do j = jlo, jhi 
          do i = ilo, ihi 

             compsn = cdtdx*(fx(i+1,j,kc,n) - fx(i,j,kc,n))

             rr = qzp(i,j,kc,QRHO)
             rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
             comps = rr*qzp(i,j,kc,nq) - compsn
             qzpo(i,j,kc,nq) = comps/rrnew

             compsn = cdtdx*(fx(i+1,j,km,n) - fx(i,j,km,n))

             rr = qzm(i,j,kc,QRHO)
             rrnew = rr - cdtdx*(fx(i+1,j,km,URHO) - fx(i,j,km,URHO))
             comps = rr*qzm(i,j,kc,nq) - compsn
             qzmo(i,j,kc,nq) = comps/rrnew

          enddo
       enddo
    enddo

    do j = jlo, jhi 
       do i = ilo, ihi 

          lambda = lam(i,j,k3d,:)
          pggp  =  qx(i+1,j,kc,GDPRES)
          pggm  =  qx(i  ,j,kc,GDPRES)
          ugp  =  qx(i+1,j,kc,GDU)
          ugm  =  qx(i  ,j,kc,GDU)
          ugc = 0.5d0*(ugp+ugm)
          ergp = qx(i+1,j,kc,GDERADS:GDERADS-1+ngroups)
          ergm = qx(i  ,j,kc,GDERADS:GDERADS-1+ngroups)

          dup = pggp*ugp - pggm*ugm
          pav = 0.5d0*(pggp+pggm)
          du = ugp-ugm

          ! Convert to conservation form
          rrrz = qzp(i,j,kc,QRHO)
          rurz = rrrz*qzp(i,j,kc,QU)
          rvrz = rrrz*qzp(i,j,kc,QV)
          rwrz = rrrz*qzp(i,j,kc,QW)
          ekenrz = 0.5d0*rrrz*(qzp(i,j,kc,QU)**2 + qzp(i,j,kc,QV)**2 + qzp(i,j,kc,QW)**2)
          rerz = qzp(i,j,kc,QREINT) + ekenrz
          err  = qzp(i,j,kc,qrad:qradhi)

          ! Add transverse predictor
          rrnewrz = rrrz - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
          runewrz = rurz - cdtdx*(fx(i+1,j,kc,UMX) - fx(i,j,kc,UMX))
          rvnewrz = rvrz - cdtdx*(fx(i+1,j,kc,UMY) - fx(i,j,kc,UMY))
          rwnewrz = rwrz - cdtdx*(fx(i+1,j,kc,UMZ) - fx(i,j,kc,UMZ))
          lamge = lambda(:) * (ergp(:)-ergm(:))
          dmom = - cdtdx*sum(lamge(:))
          runewrz = runewrz + dmom
          luge = 0.5d0*(ugp+ugm) * lamge(:)
          dre = -cdtdx*sum(luge)
          renewrz = rerz - cdtdx*(fx(i+1,j,kc,UEDEN) - fx(i,j,kc,UEDEN)) &
               + dre

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = cdtdx * ugc * f1 * (ergp(g) - ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = cdtdx * f1 * 0.5d0*(ergp(g)+ergm(g)) * (ugm-ugp)
             end do
          else ! mixed frame
             der(:) = cdtdx * luge
          end if

          ernewr  = err(:) - cdtdx*(rfx(i+1,j,kc,:) - rfx(i,j,kc,:)) &
               + der(:)


          pnewrz = qzp(i,j,kc,QPRES) - cdtdx*(dup + pav*du*(gamc(i,j,k3d)-1.d0))

          ! Convert back to non-conservation form
          qzpo(i,j,kc,QRHO) = rrnewrz
          qzpo(i,j,kc,QU) = runewrz/qzpo(i,j,kc,QRHO)
          qzpo(i,j,kc,QV) = rvnewrz/qzpo(i,j,kc,QRHO)
          qzpo(i,j,kc,QW) = rwnewrz/qzpo(i,j,kc,QRHO)
          rhoekenrz = 0.5d0*(runewrz**2 + rvnewrz**2 + rwnewrz**2)/qzpo(i,j,kc,QRHO)
          qzpo(i,j,kc,QREINT) = renewrz - rhoekenrz
          qzpo(i,j,kc,QPRES) = pnewrz
          qzpo(i,j,kc,qrad:qradhi) = ernewr(:)
          qzpo(i,j,kc,qptot  ) = sum(lambda(:)*ernewr(:)) + qzpo(i,j,kc,QPRES)
          qzpo(i,j,kc,qreitot) = sum(qzpo(i,j,kc,qrad:qradhi)) + qzpo(i,j,kc,QREINT)

          lambda = lam(i,j,k3d-1,:)
          pggp  =  qx(i+1,j,km,GDPRES)
          pggm  =  qx(i  ,j,km,GDPRES)
          ugp  =  qx(i+1,j,km,GDU)
          ugm  =  qx(i  ,j,km,GDU)
          ugc = 0.5d0*(ugp+ugm)
          ergp = qx(i+1,j,km,GDERADS:GDERADS-1+ngroups)
          ergm = qx(i  ,j,km,GDERADS:GDERADS-1+ngroups)

          dup = pggp*ugp - pggm*ugm
          pav = 0.5d0*(pggp+pggm)
          du = ugp-ugm

          rrlz = qzm(i,j,kc,QRHO)
          rulz = rrlz*qzm(i,j,kc,QU)
          rvlz = rrlz*qzm(i,j,kc,QV)
          rwlz = rrlz*qzm(i,j,kc,QW)
          ekenlz = 0.5d0*rrlz*(qzm(i,j,kc,QU)**2 + qzm(i,j,kc,QV)**2 + qzm(i,j,kc,QW)**2)
          relz = qzm(i,j,kc,QREINT) + ekenlz
          erl  = qzm(i,j,kc,qrad:qradhi)

          ! Add transverse predictor
          rrnewlz = rrlz - cdtdx*(fx(i+1,j,km,URHO) - fx(i,j,km,URHO))
          runewlz = rulz - cdtdx*(fx(i+1,j,km,UMX) - fx(i,j,km,UMX))
          rvnewlz = rvlz - cdtdx*(fx(i+1,j,km,UMY) - fx(i,j,km,UMY))
          rwnewlz = rwlz - cdtdx*(fx(i+1,j,km,UMZ) - fx(i,j,km,UMZ))
          lamge = lambda(:) * (ergp(:)-ergm(:))
          dmom = - cdtdx*sum(lamge(:))
          runewlz = runewlz + dmom
          luge = 0.5d0*(ugp+ugm) * lamge(:)
          dre = -cdtdx*sum(luge)
          renewlz = relz - cdtdx*(fx(i+1,j,km,UEDEN) - fx(i,j,km,UEDEN)) &
               + dre

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = cdtdx * ugc * f1 * (ergp(g) - ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = cdtdx * f1 * 0.5d0*(ergp(g)+ergm(g)) * (ugm-ugp)
             end do
          else ! mixed frame
             der(:) = cdtdx * luge
          end if

          ernewl  = erl(:) - cdtdx*(rfx(i+1,j,km,:) - rfx(i,j,km,:)) &
               +der(:)

          pnewlz = qzm(i,j,kc,QPRES) - cdtdx*(dup + pav*du*(gamc(i,j,k3d-1)-1.d0))

          qzmo(i,j,kc,QRHO) = rrnewlz
          qzmo(i,j,kc,QU) = runewlz/qzmo(i,j,kc,QRHO)
          qzmo(i,j,kc,QV) = rvnewlz/qzmo(i,j,kc,QRHO)
          qzmo(i,j,kc,QW) = rwnewlz/qzmo(i,j,kc,QRHO)
          rhoekenlz = 0.5d0*(runewlz**2 + rvnewlz**2 + rwnewlz**2)/qzmo(i,j,kc,QRHO)
          qzmo(i,j,kc,QREINT) = renewlz - rhoekenlz
          qzmo(i,j,kc,QPRES) = pnewlz
          qzmo(i,j,kc,qrad:qradhi) = ernewl(:)
          qzmo(i,j,kc,qptot  ) = sum(lambda(:)*ernewl(:)) + qzmo(i,j,kc,QPRES)
          qzmo(i,j,kc,qreitot) = sum(qzmo(i,j,kc,qrad:qradhi)) + qzmo(i,j,kc,QREINT)

       enddo
    enddo

  end subroutine transx2_rad

  ! ::: 
  ! ::: ------------------------------------------------------------------
  ! ::: 

  subroutine transxy_rad(lam, lam_lo, lam_hi, &
                         qm, qmo, qp, qpo, qd_lo, qd_hi, &
                         fxy, rfxy, fx_lo, fx_hi, &
                         fyx, rfyx, fy_lo, fy_hi, &
                         qx, qx_lo, qx_hi, &
                         qy, qy_lo, qy_hi, &
                         gamc, gd_lo, gd_hi, &
                         srcQ, src_lo, src_hi, &
                         hdt, cdtdx, cdtdy, ilo, ihi, jlo, jhi, kc, km, k3d)

    integer :: lam_lo(3), lam_hi(3)
    integer :: qd_lo(3), qd_hi(3)
    integer :: fx_lo(3), fx_hi(3)
    integer :: fy_lo(3), fy_hi(3)
    integer :: qx_lo(3), qx_hi(3)
    integer :: qy_lo(3), qy_hi(3)
    integer :: gd_lo(3), gd_hi(3)
    integer :: src_lo(3), src_hi(3)
    integer ilo, ihi, jlo, jhi, km, kc, k3d

    double precision lam(lam_lo(1):lam_hi(1),lam_lo(2):lam_hi(2),lam_lo(3):lam_hi(3),0:ngroups-1)
    double precision  qm(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision qmo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision  qp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision qpo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision  fxy(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),NVAR)
    double precision rfxy(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),0:ngroups-1)
    double precision  fyx(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),NVAR)
    double precision rfyx(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),0:ngroups-1)
    double precision  qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),ngdnv)
    double precision  qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),ngdnv)
    double precision gamc(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2),gd_lo(3):gd_hi(3))
    double precision srcQ(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),QVAR)
    double precision hdt,cdtdx,cdtdy

    ! Local variables

    integer i, j, g
    integer n , nq
    integer iadv, ispec, iaux

    double precision rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
    double precision rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
    double precision rrnewr, runewr, rvnewr, rwnewr, renewr
    double precision rrnewl, runewl, rvnewl, rwnewl, renewl
    double precision pnewr, pnewl
    double precision ugxp, ugxm, duxp, pxav, dux, pxnew
    double precision ugyp, ugym, duyp, pyav, duy, pynew
    double precision ugxpm, ugxmm, duxpm, pxavm, duxm, pxnewm
    double precision ugypm, ugymm, duypm, pyavm, duym, pynewm
    double precision compr, compl, compnr, compnl
    double precision rhotmp

    double precision :: dmx, dmy, dre, pggxp, pggyp, pggxm, pggym 
    double precision :: pggxpm, pggypm, pggxmm, pggymm 
    double precision, dimension(0:ngroups-1) :: der, lamc, lamm, lugex, lugey, lgex, lgey, &
         err, ernewr, erl, ernewl, ergxp, ergyp, ergxm, ergym, ergxpm, ergypm, ergxmm, ergymm
    double precision eddf, f1

    do iadv = 1, nadv
       n = UFA + iadv - 1
       nq = QFA + iadv - 1

       do j = jlo, jhi 
          do i = ilo, ihi 

             rrr = qp(i,j,kc,QRHO)
             rrl = qm(i,j,kc,QRHO)

             compr = rrr*qp(i,j,kc,nq)
             compl = rrl*qm(i,j,kc,nq)

             rrnewr = rrr - cdtdx*(fxy(i+1,j,kc,URHO) - fxy(i,j,kc,URHO)) &
                  - cdtdy*(fyx(i,j+1,kc,URHO) - fyx(i,j,kc,URHO))
             rrnewl = rrl - cdtdx*(fxy(i+1,j,km,URHO) - fxy(i,j,km,URHO)) &
                  - cdtdy*(fyx(i,j+1,km,URHO) - fyx(i,j,km,URHO))

             compnr = compr - cdtdx*(fxy(i+1,j,kc,n) - fxy(i,j,kc,n)) &
                  - cdtdy*(fyx(i,j+1,kc,n) - fyx(i,j,kc,n))
             compnl = compl - cdtdx*(fxy(i+1,j,km,n) - fxy(i,j,km,n)) &
                  - cdtdy*(fyx(i,j+1,km,n) - fyx(i,j,km,n))

             qpo(i,j,kc,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nq)
             qmo(i,j,kc,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d-1,nq)

          enddo
       enddo
    enddo

    do ispec = 1, nspec
       n = UFS + ispec - 1
       nq = QFS + ispec - 1

       do j = jlo, jhi 
          do i = ilo, ihi 

             rrr = qp(i,j,kc,QRHO)
             rrl = qm(i,j,kc,QRHO)

             compr = rrr*qp(i,j,kc,nq)
             compl = rrl*qm(i,j,kc,nq)

             rrnewr = rrr - cdtdx*(fxy(i+1,j,kc,URHO) - fxy(i,j,kc,URHO)) &
                  - cdtdy*(fyx(i,j+1,kc,URHO) - fyx(i,j,kc,URHO))
             rrnewl = rrl - cdtdx*(fxy(i+1,j,km,URHO) - fxy(i,j,km,URHO)) &
                  - cdtdy*(fyx(i,j+1,km,URHO) - fyx(i,j,km,URHO))

             compnr = compr - cdtdx*(fxy(i+1,j,kc,n) - fxy(i,j,kc,n)) &
                  - cdtdy*(fyx(i,j+1,kc,n) - fyx(i,j,kc,n))
             compnl = compl - cdtdx*(fxy(i+1,j,km,n) - fxy(i,j,km,n)) &
                  - cdtdy*(fyx(i,j+1,km,n) - fyx(i,j,km,n))

             qpo(i,j,kc,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nq)
             qmo(i,j,kc,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d-1,nq)

          enddo
       enddo
    enddo

    do iaux = 1, naux
       n  = UFX + iaux - 1
       nq = QFX + iaux - 1

       do j = jlo, jhi 
          do i = ilo, ihi 

             rrr = qp(i,j,kc,QRHO)
             rrl = qm(i,j,kc,QRHO)

             compr = rrr*qp(i,j,kc,nq)
             compl = rrl*qm(i,j,kc,nq)

             rrnewr = rrr - cdtdx*(fxy(i+1,j,kc,URHO) - fxy(i,j,kc,URHO)) &
                  - cdtdy*(fyx(i,j+1,kc,URHO) - fyx(i,j,kc,URHO))
             rrnewl = rrl - cdtdx*(fxy(i+1,j,km,URHO) - fxy(i,j,km,URHO)) &
                  - cdtdy*(fyx(i,j+1,km,URHO) - fyx(i,j,km,URHO))

             compnr = compr - cdtdx*(fxy(i+1,j,kc,n) - fxy(i,j,kc,n)) &
                  - cdtdy*(fyx(i,j+1,kc,n) - fyx(i,j,kc,n))
             compnl = compl - cdtdx*(fxy(i+1,j,km,n) - fxy(i,j,km,n)) &
                  - cdtdy*(fyx(i,j+1,km,n) - fyx(i,j,km,n))

             qpo(i,j,kc,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nq)
             qmo(i,j,kc,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d-1,nq)

          enddo
       enddo
    enddo

    do j = jlo, jhi 
       do i = ilo, ihi 

          lamc = lam(i,j,k3d,:)

          pggxp = qx(i+1,j,kc,GDPRES)
          pggxm = qx(i  ,j,kc,GDPRES)
          ugxp  = qx(i+1,j,kc,GDU)
          ugxm  = qx(i  ,j,kc,GDU)
          ergxp = qx(i+1,j,kc,GDERADS:GDERADS-1+ngroups)
          ergxm = qx(i  ,j,kc,GDERADS:GDERADS-1+ngroups)

          pggyp = qy(i,j+1,kc,GDPRES)
          pggym = qy(i,j  ,kc,GDPRES)
          ugyp  = qy(i,j+1,kc,GDV)
          ugym  = qy(i,j  ,kc,GDV)
          ergyp = qy(i,j+1,kc,GDERADS:GDERADS-1+ngroups)
          ergym = qy(i,j  ,kc,GDERADS:GDERADS-1+ngroups)

          lamm = lam(i,j,k3d-1,:)

          pggxpm = qx(i+1,j,km,GDPRES)
          pggxmm = qx(i  ,j,km,GDPRES)
          ugxpm  = qx(i+1,j,km,GDU)
          ugxmm  = qx(i  ,j,km,GDU)
          ergxpm = qx(i+1,j,km,GDERADS:GDERADS-1+ngroups)
          ergxmm = qx(i  ,j,km,GDERADS:GDERADS-1+ngroups)

          pggypm = qy(i,j+1,km,GDPRES)
          pggymm = qy(i,j  ,km,GDPRES)
          ugypm  = qy(i,j+1,km,GDV)
          ugymm  = qy(i,j  ,km,GDV)
          ergypm = qy(i,j+1,km,GDERADS:GDERADS-1+ngroups)
          ergymm = qy(i,j  ,km,GDERADS:GDERADS-1+ngroups)

          ! Convert to conservation form
          rrr = qp(i,j,kc,QRHO)
          rur = rrr*qp(i,j,kc,QU)
          rvr = rrr*qp(i,j,kc,QV)
          rwr = rrr*qp(i,j,kc,QW)
          ekenr = 0.5d0*rrr*(qp(i,j,kc,QU)**2 + qp(i,j,kc,QV)**2 + &
               qp(i,j,kc,QW)**2)
          rer = qp(i,j,kc,QREINT) + ekenr
          err = qp(i,j,kc,qrad:qradhi)

          rrl = qm(i,j,kc,QRHO)
          rul = rrl*qm(i,j,kc,QU)
          rvl = rrl*qm(i,j,kc,QV)
          rwl = rrl*qm(i,j,kc,QW)
          ekenl = 0.5d0*rrl*(qm(i,j,kc,QU)**2 + qm(i,j,kc,QV)**2 + &
               qm(i,j,kc,QW)**2)
          rel = qm(i,j,kc,QREINT) + ekenl
          erl = qm(i,j,kc,qrad:qradhi)

          ! Add transverse predictor
          rrnewr = rrr - cdtdx*(fxy(i+1,j  ,kc,URHO) - fxy(i,j,kc,URHO)) &
               - cdtdy*(fyx(i  ,j+1,kc,URHO) - fyx(i,j,kc,URHO))
          runewr = rur - cdtdx*(fxy(i+1,j  ,kc,UMX ) - fxy(i,j,kc,UMX)) &
               - cdtdy*(fyx(i  ,j+1,kc,UMX ) - fyx(i,j,kc,UMX))
          rvnewr = rvr - cdtdx*(fxy(i+1,j  ,kc,UMY ) - fxy(i,j,kc,UMY)) &
               - cdtdy*(fyx(i  ,j+1,kc,UMY ) - fyx(i,j,kc,UMY))
          rwnewr = rwr - cdtdx*(fxy(i+1,j  ,kc,UMZ ) - fxy(i,j,kc,UMZ)) &
               - cdtdy*(fyx(i  ,j+1,kc,UMZ ) - fyx(i,j,kc,UMZ))
          lgex = lamc(:) * (ergxp(:)-ergxm(:))
          lgey = lamc(:) * (ergyp(:)-ergym(:))
          dmx = - cdtdx*sum(lgex)
          dmy = - cdtdy*sum(lgey)
          runewr = runewr + dmx
          rvnewr = rvnewr + dmy
          lugex = 0.5d0*(ugxp+ugxm) * lgex(:)
          lugey = 0.5d0*(ugyp+ugym) * lgey(:)
          dre = -cdtdx*sum(lugex) - cdtdy*sum(lugey)
          renewr = rer - cdtdx*(fxy(i+1,j  ,kc,UEDEN) - fxy(i,j,kc,UEDEN)) &
               - cdtdy*(fyx(i  ,j+1,kc,UEDEN) - fyx(i,j,kc,UEDEN))  &
               + dre

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lamc(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = f1*(cdtdx*0.5d0*(ugxp+ugxm)*(ergxp(g)-ergxm(g)) &
                     +       cdtdy*0.5d0*(ugyp+ugym)*(ergyp(g)-ergym(g)) )
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lamc(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = f1*(cdtdx*0.5d0*(ergxp(g)+ergxm(g))*(ugxm-ugxp) &
                     +       cdtdy*0.5d0*(ergyp(g)+ergym(g))*(ugym-ugyp) )
             end do
          else ! mixed frame
             der(:) = cdtdx * lugex + cdtdy * lugey
          end if

          ernewr = err(:) - cdtdx*(rfxy(i+1,j,kc,:) - rfxy(i,j,kc,:)) &
               &          - cdtdy*(rfyx(i,j+1,kc,:) - rfyx(i,j,kc,:))  &
               &          + der(:)

          rrnewl = rrl - cdtdx*(fxy(i+1,j  ,km,URHO) - fxy(i,j,km,URHO)) &
               - cdtdy*(fyx(i  ,j+1,km,URHO) - fyx(i,j,km,URHO))
          runewl = rul - cdtdx*(fxy(i+1,j  ,km,UMX ) - fxy(i,j,km,UMX)) &
               - cdtdy*(fyx(i  ,j+1,km,UMX ) - fyx(i,j,km,UMX))
          rvnewl = rvl - cdtdx*(fxy(i+1,j  ,km,UMY ) - fxy(i,j,km,UMY)) &
               - cdtdy*(fyx(i  ,j+1,km,UMY ) - fyx(i,j,km,UMY))
          rwnewl = rwl - cdtdx*(fxy(i+1,j  ,km,UMZ ) - fxy(i,j,km,UMZ)) &
               - cdtdy*(fyx(i  ,j+1,km,UMZ ) - fyx(i,j,km,UMZ))
          lgex = lamm(:) * (ergxpm(:)-ergxmm(:))
          lgey = lamm(:) * (ergypm(:)-ergymm(:))
          dmx = - cdtdx*sum(lgex)
          dmy = - cdtdy*sum(lgey)
          runewl = runewl + dmx
          rvnewl = rvnewl + dmy 
          lugex = 0.5d0*(ugxpm+ugxmm) * lgex(:)
          lugey = 0.5d0*(ugypm+ugymm) * lgey(:)
          dre = -cdtdx*sum(lugex) - cdtdy*sum(lugey)
          renewl = rel - cdtdx*(fxy(i+1,j  ,km,UEDEN) - fxy(i,j,km,UEDEN)) &
               - cdtdy*(fyx(i  ,j+1,km,UEDEN) - fyx(i,j,km,UEDEN)) &
               + dre

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lamm(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = f1*(cdtdx*0.5d0*(ugxpm+ugxmm)*(ergxpm(g)-ergxmm(g)) &
                     +       cdtdy*0.5d0*(ugypm+ugymm)*(ergypm(g)-ergymm(g)) )
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lamm(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = f1*(cdtdx*0.5d0*(ergxpm(g)+ergxmm(g))*(ugxmm-ugxpm) &
                     +       cdtdy*0.5d0*(ergypm(g)+ergymm(g))*(ugymm-ugypm) )
             end do
          else ! mixed frame
             der(:) = cdtdx * lugex + cdtdy * lugey
          end if

          ernewl = erl(:) - cdtdx*(rfxy(i+1,j  ,km,:) - rfxy(i,j,km,:)) &
               &          - cdtdy*(rfyx(i  ,j+1,km,:) - rfyx(i,j,km,:)) &
               &          + der(:)

          duxp = pggxp*ugxp - pggxm*ugxm
          pxav = 0.5d0*(pggxp+pggxm)
          dux = ugxp-ugxm
          pxnew = cdtdx*(duxp + pxav*dux*(gamc(i,j,k3d)-1.d0))

          duxpm = pggxpm*ugxpm - pggxmm*ugxmm
          pxavm = 0.5d0*(pggxpm+pggxmm)
          duxm = ugxpm-ugxmm
          pxnewm = cdtdx*(duxpm + pxavm*duxm*(gamc(i,j,k3d-1)-1.d0))

          duyp = pggyp*ugyp - pggym*ugym
          pyav = 0.5d0*(pggyp+pggym)
          duy = ugyp-ugym
          pynew = cdtdy*(duyp + pyav*duy*(gamc(i,j,k3d)-1.d0))

          duypm = pggypm*ugypm - pggymm*ugymm
          pyavm = 0.5d0*(pggypm+pggymm)
          duym = ugypm-ugymm
          pynewm = cdtdy*(duypm + pyavm*duym*(gamc(i,j,k3d-1)-1.d0))

          pnewr = qp(i,j,kc,QPRES) - pxnew - pynew
          pnewl = qm(i,j,kc,QPRES) - pxnewm - pynewm

          ! Convert back to non-conservation form
          rhotmp = rrnewr
          qpo(i,j,kc,QRHO  ) = rhotmp        + hdt*srcQ(i,j,k3d,QRHO)
          qpo(i,j,kc,QU    ) = runewr/rhotmp + hdt*srcQ(i,j,k3d,QU)
          qpo(i,j,kc,QV    ) = rvnewr/rhotmp + hdt*srcQ(i,j,k3d,QV)
          qpo(i,j,kc,QW    ) = rwnewr/rhotmp + hdt*srcQ(i,j,k3d,QW)
          rhoekenr = 0.5d0*(runewr**2 + rvnewr**2 + rwnewr**2)/rhotmp
          qpo(i,j,kc,QREINT) = renewr - rhoekenr + hdt*srcQ(i,j,k3d,QREINT)
          qpo(i,j,kc,QPRES ) = pnewr             + hdt*srcQ(i,j,k3d,QPRES)
          qpo(i,j,kc,qrad:qradhi) = ernewr(:)
          qpo(i,j,kc,qptot  ) = sum(lamc(:)*ernewr(:)) + qpo(i,j,kc,QPRES)
          qpo(i,j,kc,qreitot) = sum(qpo(i,j,kc,qrad:qradhi)) + qpo(i,j,kc,QREINT)

          rhotmp = rrnewl
          qmo(i,j,kc,QRHO  ) = rhotmp        + hdt*srcQ(i,j,k3d-1,QRHO)
          qmo(i,j,kc,QU    ) = runewl/rhotmp + hdt*srcQ(i,j,k3d-1,QU)
          qmo(i,j,kc,QV    ) = rvnewl/rhotmp + hdt*srcQ(i,j,k3d-1,QV)
          qmo(i,j,kc,QW    ) = rwnewl/rhotmp + hdt*srcQ(i,j,k3d-1,QW)
          rhoekenl = 0.5d0*(runewl**2 + rvnewl**2 + rwnewl**2)/rhotmp
          qmo(i,j,kc,QREINT) = renewl - rhoekenl + hdt*srcQ(i,j,k3d-1,QREINT)
          qmo(i,j,kc,QPRES ) = pnewl             + hdt*srcQ(i,j,k3d-1,QPRES)
          qmo(i,j,kc,qrad:qradhi) = ernewl(:)
          qmo(i,j,kc,qptot  ) = sum(lamm(:)*ernewl(:)) + qmo(i,j,kc,QPRES)
          qmo(i,j,kc,qreitot) = sum(qmo(i,j,kc,qrad:qradhi)) + qmo(i,j,kc,QREINT)

       enddo
    enddo

  end subroutine transxy_rad

  ! ::: 
  ! ::: ------------------------------------------------------------------
  ! ::: 

  subroutine transz_rad(lam, lam_lo, lam_hi, &
                        qxm, qxmo, qxp, qxpo, qym, qymo, qyp, qypo, qd_lo, qd_hi, &
                        fz, rfz, fz_lo, fz_hi, &
                        qz, qz_lo, qz_hi, &
                        gamc, gd_lo, gd_hi, &
                        cdtdz, ilo, ihi, jlo, jhi, km, kc, k3d)

    integer :: lam_lo(3), lam_hi(3)
    integer :: qd_lo(3), qd_hi(3)
    integer :: fz_lo(3), fz_hi(3)
    integer :: qz_lo(3), qz_hi(3)
    integer :: gd_lo(3), gd_hi(3)
    integer ilo,ihi,jlo,jhi,km,kc,k3d

    double precision lam(lam_lo(1):lam_hi(1),lam_lo(2):lam_hi(2),lam_lo(3):lam_hi(3),0:ngroups-1)
    double precision  qxm(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision  qxp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision  qym(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision  qyp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision qxmo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision qxpo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision qymo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision qypo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision  fz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),NVAR)
    double precision rfz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),0:ngroups-1)
    double precision  qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),ngdnv)
    double precision gamc(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2),gd_lo(3):gd_hi(3))
    double precision cdtdz

    ! Local variables

    integer n, nq
    integer iadv, ispec, iaux
    integer i, j, g

    double precision rrnew, rr
    double precision compn, compu, compsn, comps
    double precision rrrx, rrry, rrlx, rrly
    double precision rurx, rury, rulx, ruly
    double precision rvrx, rvry, rvlx, rvly
    double precision rwrx, rwry, rwlx, rwly
    double precision ekenrx, ekenry, ekenlx, ekenly
    double precision rerx, rery, relx, rely
    double precision rrnewrx, rrnewry, rrnewlx, rrnewly
    double precision runewrx, runewry, runewlx, runewly
    double precision rvnewrx, rvnewry, rvnewlx, rvnewly
    double precision rwnewrx, rwnewry, rwnewlx, rwnewly
    double precision renewrx, renewry, renewlx, renewly
    double precision pnewrx, pnewry, pnewlx, pnewly
    double precision rhoekenrx, rhoekenry, rhoekenlx, rhoekenly
    double precision ugp, ugm, dup, pav, du

    double precision :: dmz, dre, pggp, pggm
    double precision, dimension(0:ngroups-1) :: der, lambda, luge, lamge, &
         ergp, errx, ernewrx, erry, ernewry, ergm, erlx, ernewlx, erly, ernewly
    double precision eddf, f1

    do iadv = 1, nadv
       n = UFA + iadv - 1
       nq = QFA + iadv - 1

       do j = jlo, jhi 
          do i = ilo, ihi 

             compn = cdtdz*(fz(i,j,kc,n) - fz(i,j,km,n))

             rr = qxp(i,j,km,QRHO)
             rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
             compu = rr*qxp(i,j,km,nq) - compn
             qxpo(i,j,km,nq) = compu/rrnew

             rr = qyp(i,j,km,QRHO)
             rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
             compu = rr*qyp(i,j,km,nq) - compn
             qypo(i,j,km,nq) = compu/rrnew

             rr = qxm(i+1,j,km,QRHO)
             rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
             compu = rr*qxm(i+1,j,km,nq) - compn
             qxmo(i+1,j,km,nq) = compu/rrnew

             rr = qym(i,j+1,km,QRHO)
             rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
             compu = rr*qym(i,j+1,km,nq) - compn
             qymo(i,j+1,km,nq) = compu/rrnew

          enddo
       enddo
    enddo

    do ispec = 1, nspec 
       n = UFS + ispec - 1
       nq = QFS + ispec  - 1

       do j = jlo, jhi 
          do i = ilo, ihi 

             compsn = cdtdz*(fz(i,j,kc,n) - fz(i,j,km,n))

             rr = qxp(i,j,km,QRHO)
             rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
             comps = rr*qxp(i,j,km,nq) - compsn
             qxpo(i,j,km,nq) = comps/rrnew

             rr = qyp(i,j,km,QRHO)
             rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
             comps = rr*qyp(i,j,km,nq) - compsn
             qypo(i,j,km,nq) = comps/rrnew

             rr = qxm(i+1,j,km,QRHO)
             rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
             comps = rr*qxm(i+1,j,km,nq) - compsn
             qxmo(i+1,j,km,nq) = comps/rrnew

             rr = qym(i,j+1,km,QRHO)
             rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
             comps = rr*qym(i,j+1,km,nq) - compsn
             qymo(i,j+1,km,nq) = comps/rrnew

          enddo
       enddo
    enddo

    do iaux = 1, naux 
       n  = UFX + iaux - 1
       nq = QFX + iaux  - 1

       do j = jlo, jhi 
          do i = ilo, ihi 

             compsn = cdtdz*(fz(i,j,kc,n) - fz(i,j,km,n))

             rr = qxp(i,j,km,QRHO)
             rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
             comps = rr*qxp(i,j,km,nq) - compsn
             qxpo(i,j,km,nq) = comps/rrnew

             rr = qyp(i,j,km,QRHO)
             rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
             comps = rr*qyp(i,j,km,nq) - compsn
             qypo(i,j,km,nq) = comps/rrnew

             rr = qxm(i+1,j,km,QRHO)
             rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
             comps = rr*qxm(i+1,j,km,nq) - compsn
             qxmo(i+1,j,km,nq) = comps/rrnew

             rr = qym(i,j+1,km,QRHO)
             rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
             comps = rr*qym(i,j+1,km,nq) - compsn
             qymo(i,j+1,km,nq) = comps/rrnew

          enddo
       enddo
    enddo

    do j = jlo, jhi 
       do i = ilo, ihi 

          lambda = lam(i,j,k3d-1,:)

          pggp = qz(i,j,kc,GDPRES)
          pggm = qz(i,j,km,GDPRES)
          ugp  = qz(i,j,kc,GDW)
          ugm  = qz(i,j,km,GDW)
          ergp = qz(i,j,kc,GDERADS:GDERADS-1+ngroups)
          ergm = qz(i,j,km,GDERADS:GDERADS-1+ngroups)

          ! Convert to conservation form
          rrrx = qxp(i,j,km,QRHO)
          rurx = rrrx*qxp(i,j,km,QU)
          rvrx = rrrx*qxp(i,j,km,QV)
          rwrx = rrrx*qxp(i,j,km,QW)
          ekenrx = 0.5d0*rrrx*(qxp(i,j,km,QU)**2 + qxp(i,j,km,QV)**2 &
               + qxp(i,j,km,QW)**2)
          rerx = qxp(i,j,km,QREINT) + ekenrx
          errx = qxp(i,j,km,qrad:qradhi)

          rrry = qyp(i,j,km,QRHO)
          rury = rrry*qyp(i,j,km,QU)
          rvry = rrry*qyp(i,j,km,QV)
          rwry = rrry*qyp(i,j,km,QW)
          ekenry = 0.5d0*rrry*(qyp(i,j,km,QU)**2 + qyp(i,j,km,QV)**2 &
               + qyp(i,j,km,QW)**2)
          rery = qyp(i,j,km,QREINT) + ekenry
          erry = qyp(i,j,km,qrad:qradhi)

          rrlx = qxm(i+1,j,km,QRHO)
          rulx = rrlx*qxm(i+1,j,km,QU)
          rvlx = rrlx*qxm(i+1,j,km,QV)
          rwlx = rrlx*qxm(i+1,j,km,QW)
          ekenlx = 0.5d0*rrlx*(qxm(i+1,j,km,QU)**2 + qxm(i+1,j,km,QV)**2 &
               + qxm(i+1,j,km,QW)**2)
          relx = qxm(i+1,j,km,QREINT) + ekenlx
          erlx = qxm(i+1,j,km,qrad:qradhi)

          rrly = qym(i,j+1,km,QRHO)
          ruly = rrly*qym(i,j+1,km,QU)
          rvly = rrly*qym(i,j+1,km,QV)
          rwly = rrly*qym(i,j+1,km,QW)
          ekenly = 0.5d0*rrly*(qym(i,j+1,km,QU)**2 + qym(i,j+1,km,QV)**2 &
               + qym(i,j+1,km,QW)**2)
          rely = qym(i,j+1,km,QREINT) + ekenly
          erly = qym(i,j+1,km,qrad:qradhi)

          ! Add transverse predictor
          rrnewrx = rrrx - cdtdz*(fz(i,j,kc,URHO ) - fz(i,j,km,URHO))
          runewrx = rurx - cdtdz*(fz(i,j,kc,UMX  ) - fz(i,j,km,UMX))
          rvnewrx = rvrx - cdtdz*(fz(i,j,kc,UMY  ) - fz(i,j,km,UMY))
          rwnewrx = rwrx - cdtdz*(fz(i,j,kc,UMZ  ) - fz(i,j,km,UMZ))
          lamge = lambda(:) * (ergp(:)-ergm(:))
          dmz = - cdtdz*sum(lamge)
          rwnewrx = rwnewrx + dmz
          luge = 0.5d0*(ugp+ugm) * lamge(:)
          dre = -cdtdz*sum(luge)
          renewrx = rerx - cdtdz*(fz(i,j,kc,UEDEN) - fz(i,j,km,UEDEN)) &
               + dre

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = cdtdz*0.5d0*(ugp+ugm)*(ergp(g)-ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = cdtdz*0.5d0*(ergp(g)+ergm(g))*(ugm-ugp)*f1
             end do
          else ! mixed frame
             der(:) = cdtdz * luge
          end if

          ernewrx = errx(:) - cdtdz*(rfz(i,j,kc,:) - rfz(i,j,km,:)) &
               + der(:)

          rrnewry = rrry - cdtdz*(fz(i,j,kc,URHO ) - fz(i,j,km,URHO))
          runewry = rury - cdtdz*(fz(i,j,kc,UMX  ) - fz(i,j,km,UMX))
          rvnewry = rvry - cdtdz*(fz(i,j,kc,UMY  ) - fz(i,j,km,UMY))
          rwnewry = rwry - cdtdz*(fz(i,j,kc,UMZ  ) - fz(i,j,km,UMZ))
          rwnewry = rwnewry + dmz
          renewry = rery - cdtdz*(fz(i,j,kc,UEDEN) - fz(i,j,km,UEDEN)) &
               + dre
          ernewry = erry(:) - cdtdz*(rfz(i,j,kc,:) - rfz(i,j,km,:)) &
               + der(:)

          rrnewlx = rrlx - cdtdz*(fz(i,j,kc,URHO ) - fz(i,j,km,URHO))
          runewlx = rulx - cdtdz*(fz(i,j,kc,UMX  ) - fz(i,j,km,UMX))
          rvnewlx = rvlx - cdtdz*(fz(i,j,kc,UMY  ) - fz(i,j,km,UMY))
          rwnewlx = rwlx - cdtdz*(fz(i,j,kc,UMZ  ) - fz(i,j,km,UMZ))
          rwnewlx = rwnewlx + dmz
          renewlx = relx - cdtdz*(fz(i,j,kc,UEDEN) - fz(i,j,km,UEDEN)) &
               + dre
          ernewlx = erlx(:) - cdtdz*(rfz(i,j,kc,:) - rfz(i,j,km,:)) &
               + der(:)

          rrnewly = rrly - cdtdz*(fz(i,j,kc,URHO ) - fz(i,j,km,URHO))
          runewly = ruly - cdtdz*(fz(i,j,kc,UMX  ) - fz(i,j,km,UMX))
          rvnewly = rvly - cdtdz*(fz(i,j,kc,UMY  ) - fz(i,j,km,UMY))
          rwnewly = rwly - cdtdz*(fz(i,j,kc,UMZ  ) - fz(i,j,km,UMZ))
          rwnewly = rwnewly + dmz
          renewly = rely - cdtdz*(fz(i,j,kc,UEDEN) - fz(i,j,km,UEDEN)) &
               + dre
          ernewly = erly(:) - cdtdz*(rfz(i,j,kc,:) - rfz(i,j,km,:)) &
               + der(:)

          dup = pggp*ugp - pggm*ugm
          pav = 0.5d0*(pggp+pggm)
          du = ugp-ugm

          pnewrx = qxp(i  ,j,km,QPRES) - cdtdz*(dup + pav*du*(gamc(i,j,k3d-1) - 1.d0))
          pnewlx = qxm(i+1,j,km,QPRES) - cdtdz*(dup + pav*du*(gamc(i,j,k3d-1) - 1.d0))

          pnewry = qyp(i,j  ,km,QPRES) - cdtdz*(dup + pav*du*(gamc(i,j,k3d-1) - 1.d0))
          pnewly = qym(i,j+1,km,QPRES) - cdtdz*(dup + pav*du*(gamc(i,j,k3d-1) - 1.d0))

          ! Convert back to non-conservation form
          qxpo(i,j,km,QRHO) = rrnewrx
          qxpo(i,j,km,QU) = runewrx/qxpo(i,j,km,QRHO)
          qxpo(i,j,km,QV) = rvnewrx/qxpo(i,j,km,QRHO)
          qxpo(i,j,km,QW) = rwnewrx/qxpo(i,j,km,QRHO)
          rhoekenrx = 0.5d0*(runewrx**2 + rvnewrx**2 + rwnewrx**2)/qxpo(i,j,km,QRHO)
          qxpo(i,j,km,QREINT)= renewrx - rhoekenrx
          qxpo(i,j,km,QPRES) = pnewrx
          qxpo(i,j,km,qrad:qradhi) = ernewrx(:)
          qxpo(i,j,km,qptot  ) = sum(lambda(:)*ernewrx(:)) + qxpo(i,j,km,QPRES)
          qxpo(i,j,km,qreitot) = sum(qxpo(i,j,km,qrad:qradhi)) + qxpo(i,j,km,QREINT)

          qypo(i,j,km,QRHO) = rrnewry
          qypo(i,j,km,QU) = runewry/qypo(i,j,km,QRHO)
          qypo(i,j,km,QV) = rvnewry/qypo(i,j,km,QRHO)
          qypo(i,j,km,QW) = rwnewry/qypo(i,j,km,QRHO)
          rhoekenry = 0.5d0*(runewry**2 + rvnewry**2 + rwnewry**2)/qypo(i,j,km,QRHO)
          qypo(i,j,km,QREINT)= renewry - rhoekenry
          qypo(i,j,km,QPRES) = pnewry
          qypo(i,j,km,qrad:qradhi) = ernewry(:)
          qypo(i,j,km,qptot  ) = sum(lambda(:)*ernewry(:)) + qypo(i,j,km,QPRES)
          qypo(i,j,km,qreitot) = sum(qypo(i,j,km,qrad:qradhi)) + qypo(i,j,km,QREINT)

          qxmo(i+1,j,km,QRHO) = rrnewlx
          qxmo(i+1,j,km,QU) = runewlx/qxmo(i+1,j,km,QRHO)
          qxmo(i+1,j,km,QV) = rvnewlx/qxmo(i+1,j,km,QRHO)
          qxmo(i+1,j,km,QW) = rwnewlx/qxmo(i+1,j,km,QRHO)
          rhoekenlx = 0.5d0*(runewlx**2 + rvnewlx**2 + rwnewlx**2)/qxmo(i+1,j,km,QRHO)
          qxmo(i+1,j,km,QREINT)= renewlx - rhoekenlx
          qxmo(i+1,j,km,QPRES) = pnewlx
          qxmo(i+1,j,km,qrad:qradhi) = ernewlx(:)
          qxmo(i+1,j,km,qptot  ) = sum(lambda(:)*ernewlx(:)) + qxmo(i+1,j,km,QPRES)
          qxmo(i+1,j,km,qreitot) = sum(qxmo(i+1,j,km,qrad:qradhi)) + qxmo(i+1,j,km,QREINT)

          qymo(i,j+1,km,QRHO) = rrnewly
          qymo(i,j+1,km,QU) = runewly/qymo(i,j+1,km,QRHO)
          qymo(i,j+1,km,QV) = rvnewly/qymo(i,j+1,km,QRHO)
          qymo(i,j+1,km,QW) = rwnewly/qymo(i,j+1,km,QRHO)
          rhoekenly = 0.5d0*(runewly**2 + rvnewly**2 + rwnewly**2)/qymo(i,j+1,km,QRHO)
          qymo(i,j+1,km,QREINT)= renewly - rhoekenly
          qymo(i,j+1,km,QPRES) = pnewly
          qymo(i,j+1,km,qrad:qradhi) = ernewly(:)
          qymo(i,j+1,km,qptot  ) = sum(lambda(:)*ernewly(:)) + qymo(i,j+1,km,QPRES)
          qymo(i,j+1,km,qreitot) = sum(qymo(i,j+1,km,qrad:qradhi)) + qymo(i,j+1,km,QREINT)

       enddo
    enddo

  end subroutine transz_rad

  ! ::: 
  ! ::: ------------------------------------------------------------------
  ! ::: 

  subroutine transyz_rad(lam, lam_lo, lam_hi, &
                         qm, qmo, qp, qpo, qd_lo, qd_hi, &
                         fyz, rfyz, fy_lo, fy_hi, &
                         fzy, rfzy, fz_lo, fz_hi, &
                         qy, qy_lo, qy_hi, &
                         qz, qz_lo, qz_hi, &
                         gamc, gc_lo, gc_hi, &
                         srcQ, src_lo, src_hi, &
                         hdt, cdtdy, cdtdz, ilo, ihi, jlo, jhi, km, kc, k3d)
    
    integer :: lam_lo(3), lam_hi(3)
    integer :: qd_lo(3), qd_hi(3)
    integer :: fy_lo(3), fy_hi(3)
    integer :: fz_lo(3), fz_hi(3)
    integer :: qy_lo(3), qy_hi(3)
    integer :: qz_lo(3), qz_hi(3)
    integer :: gc_lo(3), gc_hi(3)
    integer :: src_lo(3), src_hi(3)
    integer ilo,ihi,jlo,jhi,km,kc,k3d

    double precision lam(lam_lo(1):lam_hi(1),lam_lo(2):lam_hi(2),lam_lo(3):lam_hi(3),0:ngroups-1)
    double precision qm(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision qp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision qmo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision qpo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision  fyz(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),NVAR)
    double precision rfyz(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),0:ngroups-1)
    double precision  fzy(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),NVAR)
    double precision rfzy(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),0:ngroups-1)
    double precision  qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),ngdnv)
    double precision  qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),ngdnv)

    double precision gamc(gc_lo(1):gc_hi(1),gc_lo(2):gc_hi(2),gc_lo(3):gc_hi(3))
    double precision srcQ(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),QVAR)
    double precision hdt,cdtdy,cdtdz

    ! Local variables

    integer i, j, g
    integer n, nq
    integer iadv, ispec, iaux

    double precision rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
    double precision rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
    double precision rrnewr, runewr, rvnewr, rwnewr, renewr
    double precision rrnewl, runewl, rvnewl, rwnewl, renewl
    double precision pnewr, pnewl
    double precision ugyp, ugym, duyp, pyav, duy, pynew
    double precision ugzp, ugzm, duzp, pzav, duz, pznew
    double precision compr, compl, compnr, compnl
    double precision rhotmp

    double precision :: dmy, dmz, dre, pggzp, pggzm, pggyp, pggym 
    double precision, dimension(0:ngroups-1) :: der, lambda, lugey, lugez, lgey, lgez, &
         err, ernewr, erl, ernewl, ergzp, ergyp, ergzm, ergym
    double precision eddf, f1

    do iadv = 1, nadv
       n = UFA + iadv - 1
       nq = QFA + iadv - 1

       do j = jlo, jhi 
          do i = ilo, ihi 

             rrr = qp(i,j,km,QRHO)
             rrl = qm(i+1,j,km,QRHO)

             compr = rrr*qp(i,j,km,nq)
             compl = rrl*qm(i+1,j,km,nq)

             rrnewr = rrr - cdtdy*(fyz(i,j+1,km,URHO) - fyz(i,j,km,URHO)) &
                  - cdtdz*(fzy(i,j  ,kc,URHO) - fzy(i,j,km,URHO))
             rrnewl = rrl - cdtdy*(fyz(i,j+1,km,URHO) - fyz(i,j,km,URHO)) &
                  - cdtdz*(fzy(i,j  ,kc,URHO) - fzy(i,j,km,URHO))

             compnr = compr - cdtdy*(fyz(i,j+1,km,n) - fyz(i,j,km,n)) &
                  - cdtdz*(fzy(i,j  ,kc,n) - fzy(i,j,km,n))
             compnl = compl - cdtdy*(fyz(i,j+1,km,n) - fyz(i,j,km,n)) &
                  - cdtdz*(fzy(i,j  ,kc,n) - fzy(i,j,km,n))

             qpo(i,j,km,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nq)
             qmo(i+1,j,km,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d,nq)

          enddo
       enddo
    enddo

    do ispec = 1, nspec
       n = UFS + ispec - 1
       nq = QFS + ispec - 1

       do j = jlo, jhi 
          do i = ilo, ihi 

             rrr = qp(i  ,j,km,QRHO)
             rrl = qm(i+1,j,km,QRHO)

             compr = rrr*qp(i  ,j,km,nq)
             compl = rrl*qm(i+1,j,km,nq)

             rrnewr = rrr - cdtdy*(fyz(i,j+1,km,URHO) - fyz(i,j,km,URHO)) &
                  - cdtdz*(fzy(i,j  ,kc,URHO) - fzy(i,j,km,URHO))
             rrnewl = rrl - cdtdy*(fyz(i,j+1,km,URHO) - fyz(i,j,km,URHO)) &
                  - cdtdz*(fzy(i,j  ,kc,URHO) - fzy(i,j,km,URHO))

             compnr = compr - cdtdy*(fyz(i,j+1,km,n) - fyz(i,j,km,n)) &
                  - cdtdz*(fzy(i,j  ,kc,n) - fzy(i,j,km,n))
             compnl = compl - cdtdy*(fyz(i,j+1,km,n) - fyz(i,j,km,n)) &
                  - cdtdz*(fzy(i,j  ,kc,n) - fzy(i,j,km,n))

             qpo(i  ,j,km,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nq)
             qmo(i+1,j,km,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d,nq)

          enddo
       enddo
    enddo

    do iaux = 1, naux
       n  = UFX + iaux - 1
       nq = QFX + iaux - 1

       do j = jlo, jhi 
          do i = ilo, ihi 

             rrr = qp(i,j,km,QRHO)
             rrl = qm(i+1,j,km,QRHO)

             compr = rrr*qp(i  ,j,km,nq)
             compl = rrl*qm(i+1,j,km,nq)

             rrnewr = rrr - cdtdy*(fyz(i,j+1,km,URHO) - fyz(i,j,km,URHO)) &
                  - cdtdz*(fzy(i,j  ,kc,URHO) - fzy(i,j,km,URHO))
             rrnewl = rrl - cdtdy*(fyz(i,j+1,km,URHO) - fyz(i,j,km,URHO)) &
                  - cdtdz*(fzy(i,j  ,kc,URHO) - fzy(i,j,km,URHO))

             compnr = compr - cdtdy*(fyz(i,j+1,km,n) - fyz(i,j,km,n)) &
                  - cdtdz*(fzy(i,j  ,kc,n) - fzy(i,j,km,n))
             compnl = compl - cdtdy*(fyz(i,j+1,km,n) - fyz(i,j,km,n)) &
                  - cdtdz*(fzy(i,j  ,kc,n) - fzy(i,j,km,n))

             qpo(i  ,j,km,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nq)
             qmo(i+1,j,km,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d,nq)

          enddo
       enddo
    enddo

    do j = jlo, jhi 
       do i = ilo, ihi 

          lambda = lam(i,j,k3d,:)

          pggyp = qy(i,j+1,km,GDPRES)
          pggym = qy(i,j  ,km,GDPRES)
          ugyp  = qy(i,j+1,km,GDV)
          ugym  = qy(i,j  ,km,GDV)
          ergyp = qy(i,j+1,km,GDERADS:GDERADS-1+ngroups)
          ergym = qy(i,j  ,km,GDERADS:GDERADS-1+ngroups)

          pggzp = qz(i,j,kc,GDPRES)
          pggzm = qz(i,j,km,GDPRES)
          ugzp  = qz(i,j,kc,GDW)
          ugzm  = qz(i,j,km,GDW)
          ergzp = qz(i,j,kc,GDERADS:GDERADS-1+ngroups)
          ergzm = qz(i,j,km,GDERADS:GDERADS-1+ngroups)

          ! Convert to conservation form
          rrr = qp(i,j,km,QRHO)
          rur = rrr*qp(i,j,km,QU)
          rvr = rrr*qp(i,j,km,QV)
          rwr = rrr*qp(i,j,km,QW)
          ekenr = 0.5d0*rrr*(qp(i,j,km,QU)**2 + qp(i,j,km,QV)**2 + &
               qp(i,j,km,QW)**2)
          rer = qp(i,j,km,QREINT) + ekenr
          err = qp(i,j,km,qrad:qradhi)

          rrl = qm(i+1,j,km,QRHO)
          rul = rrl*qm(i+1,j,km,QU)
          rvl = rrl*qm(i+1,j,km,QV)
          rwl = rrl*qm(i+1,j,km,QW)
          ekenl = 0.5d0*rrl*(qm(i+1,j,km,QU)**2 + qm(i+1,j,km,QV)**2 + &
               qm(i+1,j,km,QW)**2)
          rel = qm(i+1,j,km,QREINT) + ekenl
          erl = qm(i+1,j,km,qrad:qradhi)

          ! Add transverse predictor
          rrnewr = rrr - cdtdy*(fyz(i,j+1,km,URHO ) - fyz(i,j,km,URHO)) &
               - cdtdz*(fzy(i,j  ,kc,URHO ) - fzy(i,j,km,URHO))
          runewr = rur - cdtdy*(fyz(i,j+1,km,UMX  ) - fyz(i,j,km,UMX)) &
               - cdtdz*(fzy(i,j  ,kc,UMX  ) - fzy(i,j,km,UMX))
          rvnewr = rvr - cdtdy*(fyz(i,j+1,km,UMY  ) - fyz(i,j,km,UMY)) &
               - cdtdz*(fzy(i,j  ,kc,UMY  ) - fzy(i,j,km,UMY))
          rwnewr = rwr - cdtdy*(fyz(i,j+1,km,UMZ  ) - fyz(i,j,km,UMZ)) &
               - cdtdz*(fzy(i,j  ,kc,UMZ  ) - fzy(i,j,km,UMZ))
          lgey = lambda(:) * (ergyp(:)-ergym(:))
          lgez = lambda(:) * (ergzp(:)-ergzm(:))
          dmy = - cdtdy*sum(lgey)
          dmz = - cdtdz*sum(lgez)
          rvnewr = rvnewr + dmy 
          rwnewr = rwnewr + dmz
          lugey = 0.5d0*(ugyp+ugym) * lgey(:)
          lugez = 0.5d0*(ugzp+ugzm) * lgez(:)
          dre = -cdtdy*sum(lugey) - cdtdz*sum(lugez)
          renewr = rer - cdtdy*(fyz(i,j+1,km,UEDEN) - fyz(i,j,km,UEDEN)) &
               - cdtdz*(fzy(i,j  ,kc,UEDEN) - fzy(i,j,km,UEDEN)) &
               + dre

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = f1*(cdtdy*0.5d0*(ugyp+ugym)*(ergyp(g)-ergym(g)) &
                     +       cdtdz*0.5d0*(ugzp+ugzm)*(ergzp(g)-ergzm(g)) )
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = f1*(cdtdy*0.5d0*(ergyp(g)+ergym(g))*(ugym-ugyp) &
                     +       cdtdz*0.5d0*(ergzp(g)+ergzm(g))*(ugzm-ugzp) )
             end do
          else ! mixed frame
             der(:) = cdtdy*lugey + cdtdz*lugez
          end if

          ernewr = err(:) - cdtdy*(rfyz(i,j+1,km,:) - rfyz(i,j,km,:)) &
               &          - cdtdz*(rfzy(i,j  ,kc,:) - rfzy(i,j,km,:)) &
               &          + der(:)

          rrnewl = rrl - cdtdy*(fyz(i,j+1,km,URHO ) - fyz(i,j,km,URHO)) &
               - cdtdz*(fzy(i,j  ,kc,URHO ) - fzy(i,j,km,URHO))
          runewl = rul - cdtdy*(fyz(i,j+1,km,UMX  ) - fyz(i,j,km,UMX)) &
               - cdtdz*(fzy(i,j  ,kc,UMX  ) - fzy(i,j,km,UMX))
          rvnewl = rvl - cdtdy*(fyz(i,j+1,km,UMY  ) - fyz(i,j,km,UMY)) &
               - cdtdz*(fzy(i,j  ,kc,UMY  ) - fzy(i,j,km,UMY))
          rwnewl = rwl - cdtdy*(fyz(i,j+1,km,UMZ  ) - fyz(i,j,km,UMZ)) &
               - cdtdz*(fzy(i,j  ,kc,UMZ  ) - fzy(i,j,km,UMZ))
          rvnewl = rvnewl + dmy
          rwnewl = rwnewl + dmz 
          renewl = rel - cdtdy*(fyz(i,j+1,km,UEDEN) - fyz(i,j,km,UEDEN)) &
               - cdtdz*(fzy(i,j  ,kc,UEDEN) - fzy(i,j,km,UEDEN)) &
               + dre
          ernewl = erl(:) - cdtdy*(rfyz(i,j+1,km,:) - rfyz(i,j,km,:)) &
               &          - cdtdz*(rfzy(i,j  ,kc,:) - rfzy(i,j,km,:)) &
               &          + der(:)

          duyp = pggyp*ugyp - pggym*ugym
          pyav = 0.5d0*(pggyp+pggym)
          duy = ugyp-ugym
          pynew = cdtdy*(duyp + pyav*duy*(gamc(i,j,k3d)-1.d0))

          duzp = pggzp*ugzp - pggzm*ugzm
          pzav = 0.5d0*(pggzp+pggzm)
          duz = ugzp-ugzm
          pznew = cdtdz*(duzp + pzav*duz*(gamc(i,j,k3d)-1.d0))

          pnewr = qp(i,j,km,QPRES) - pynew - pznew
          pnewl = qm(i+1,j,km,QPRES) - pynew - pznew

          ! Convert back to non-conservation form
          rhotmp = rrnewr
          qpo(i,j,km,QRHO  ) = rhotmp        + hdt*srcQ(i,j,k3d,QRHO)
          qpo(i,j,km,QU    ) = runewr/rhotmp + hdt*srcQ(i,j,k3d,QU)
          qpo(i,j,km,QV    ) = rvnewr/rhotmp + hdt*srcQ(i,j,k3d,QV)
          qpo(i,j,km,QW    ) = rwnewr/rhotmp + hdt*srcQ(i,j,k3d,QW)
          rhoekenr = 0.5d0*(runewr**2 + rvnewr**2 + rwnewr**2)/rhotmp
          qpo(i,j,km,QREINT)= renewr - rhoekenr + hdt*srcQ(i,j,k3d,QREINT)
          qpo(i,j,km,QPRES ) = pnewr            + hdt*srcQ(i,j,k3d,QPRES)
          qpo(i,j,km,qrad:qradhi) = ernewr(:)
          qpo(i,j,km,qptot  ) = sum(lambda(:)*ernewr(:)) + qpo(i,j,km,QPRES)
          qpo(i,j,km,qreitot) = sum(qpo(i,j,km,qrad:qradhi)) + qpo(i,j,km,QREINT)

          rhotmp = rrnewl
          qmo(i+1,j,km,QRHO  ) = rhotmp        + hdt*srcQ(i,j,k3d,QRHO)
          qmo(i+1,j,km,QU    ) = runewl/rhotmp + hdt*srcQ(i,j,k3d,QU)
          qmo(i+1,j,km,QV    ) = rvnewl/rhotmp + hdt*srcQ(i,j,k3d,QV)
          qmo(i+1,j,km,QW    ) = rwnewl/rhotmp + hdt*srcQ(i,j,k3d,QW)
          rhoekenl = 0.5d0*(runewl**2 + rvnewl**2 + rwnewl**2)/rhotmp
          qmo(i+1,j,km,QREINT)= renewl - rhoekenl + hdt*srcQ(i,j,k3d,QREINT)
          qmo(i+1,j,km,QPRES ) = pnewl            + hdt*srcQ(i,j,k3d,QPRES)
          qmo(i+1,j,km,qrad:qradhi) = ernewl(:)
          qmo(i+1,j,km,qptot) = sum(lambda(:)*ernewl(:)) + qmo(i+1,j,km,QPRES)
          qmo(i+1,j,km,qreitot) = sum(qmo(i+1,j,km,qrad:qradhi)) + qmo(i+1,j,km,QREINT)
       enddo
    enddo

  end subroutine transyz_rad

  ! ::: 
  ! ::: ------------------------------------------------------------------
  ! ::: 

  subroutine transxz_rad(lam, lam_lo, lam_hi, &
                         qm, qmo, qp, qpo, qd_lo, qd_hi, &
                         fxz, rfxz, fx_lo, fx_hi, &
                         fzx, rfzx, fz_lo, fz_hi, &
                         qx, qx_lo, qx_hi, &
                         qz, qz_lo, qz_hi, &
                         gamc, gc_lo, gc_hi, &
                         srcQ, src_lo, src_hi, &
                         hdt, cdtdx, cdtdz, ilo, ihi, jlo, jhi, km, kc, k3d)

    integer :: lam_lo(3), lam_hi(3)
    integer :: qd_lo(3), qd_hi(3)
    integer :: fx_lo(3), fx_hi(3)
    integer :: fz_lo(3), fz_hi(3)
    integer :: qx_lo(3), qx_hi(3)
    integer :: qz_lo(3), qz_hi(3)
    integer :: gc_lo(3), gc_hi(3)
    integer :: src_lo(3), src_hi(3)
    integer ilo, ihi, jlo, jhi, km, kc, k3d

    double precision lam(lam_lo(1):lam_hi(1),lam_lo(2):lam_hi(2),lam_lo(3):lam_hi(3),0:ngroups-1)
    double precision  qm(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision  qp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision qmo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision qpo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision  fxz(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),NVAR)
    double precision rfxz(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),0:ngroups-1)
    double precision  fzx(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),NVAR) 
    double precision rfzx(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),0:ngroups-1)
    double precision  qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),ngdnv)
    double precision  qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),ngdnv)
    double precision gamc(gc_lo(1):gc_hi(1),gc_lo(2):gc_hi(2),gc_lo(3):gc_hi(3))
    double precision srcQ(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),QVAR)
    double precision hdt,cdtdx,cdtdz

    ! Local variables

    integer i, j, g
    integer n, nq
    integer iadv, ispec, iaux

    double precision rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
    double precision rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
    double precision rrnewr, runewr, rvnewr, rwnewr, renewr
    double precision rrnewl, runewl, rvnewl, rwnewl, renewl
    double precision pnewr, pnewl
    double precision ugxp, ugxm, duxp, pxav, dux, pxnew
    double precision ugzp, ugzm, duzp, pzav, duz, pznew
    double precision compr, compl, compnr, compnl
    double precision rhotmp

    double precision :: dmx, dmz, dre, pggxp, pggxm, pggzp, pggzm
    double precision, dimension(0:ngroups-1) :: der, lambda, lugex, lugez, lgex, lgez, &
         err, ernewr, erl, ernewl, ergzp, ergxp, ergzm,  ergxm 
    double precision eddf, f1

    do iadv = 1, nadv
       n = UFA + iadv - 1
       nq = QFA + iadv -1 

       do j = jlo, jhi 
          do i = ilo, ihi 

             rrr = qp(i,j,km,QRHO)
             rrl = qm(i,j+1,km,QRHO)

             compr = rrr*qp(i,j,km,nq)
             compl = rrl*qm(i,j+1,km,nq)

             rrnewr = rrr - cdtdx*(fxz(i+1,j,km,URHO) - fxz(i,j,km,URHO)) &
                  - cdtdz*(fzx(i  ,j,kc,URHO) - fzx(i,j,km,URHO))
             rrnewl = rrl - cdtdx*(fxz(i+1,j,km,URHO) - fxz(i,j,km,URHO)) &
                  - cdtdz*(fzx(i  ,j,kc,URHO) - fzx(i,j,km,URHO))

             compnr = compr - cdtdx*(fxz(i+1,j,km,n) - fxz(i,j,km,n)) &
                  - cdtdz*(fzx(i  ,j,kc,n) - fzx(i,j,km,n))
             compnl = compl - cdtdx*(fxz(i+1,j,km,n) - fxz(i,j,km,n)) &
                  - cdtdz*(fzx(i  ,j,kc,n) - fzx(i,j,km,n))

             qpo(i,j,km,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nq)
             qmo(i,j+1,km,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d,nq)

          enddo
       enddo
    enddo

    do ispec = 1, nspec
       n = UFS + ispec - 1
       nq = QFS + ispec - 1

       do j = jlo, jhi 
          do i = ilo, ihi 

             rrr = qp(i,j,km,QRHO)
             rrl = qm(i,j+1,km,QRHO)

             compr = rrr*qp(i,j,km,nq)
             compl = rrl*qm(i,j+1,km,nq)

             rrnewr = rrr - cdtdx*(fxz(i+1,j,km,URHO) - fxz(i,j,km,URHO)) &
                  - cdtdz*(fzx(i  ,j,kc,URHO) - fzx(i,j,km,URHO))
             rrnewl = rrl - cdtdx*(fxz(i+1,j,km,URHO) - fxz(i,j,km,URHO)) &
                  - cdtdz*(fzx(i  ,j,kc,URHO) - fzx(i,j,km,URHO))

             compnr = compr - cdtdx*(fxz(i+1,j,km,n) - fxz(i,j,km,n)) &
                  - cdtdz*(fzx(i  ,j,kc,n) - fzx(i,j,km,n))
             compnl = compl - cdtdx*(fxz(i+1,j,km,n) - fxz(i,j,km,n)) &
                  - cdtdz*(fzx(i  ,j,kc,n) - fzx(i,j,km,n))

             qpo(i,j,km,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nq)
             qmo(i,j+1,km,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d,nq)

          enddo
       enddo
    enddo

    do iaux = 1, naux
       n  = UFX + iaux - 1
       nq = QFX + iaux - 1

       do j = jlo, jhi 
          do i = ilo, ihi 

             rrr = qp(i,j,km,QRHO)
             rrl = qm(i,j+1,km,QRHO)

             compr = rrr*qp(i,j  ,km,nq)
             compl = rrl*qm(i,j+1,km,nq)

             rrnewr = rrr - cdtdx*(fxz(i+1,j,km,URHO) - fxz(i,j,km,URHO)) &
                  - cdtdz*(fzx(i  ,j,kc,URHO) - fzx(i,j,km,URHO))
             rrnewl = rrl - cdtdx*(fxz(i+1,j,km,URHO) - fxz(i,j,km,URHO)) &
                  - cdtdz*(fzx(i  ,j,kc,URHO) - fzx(i,j,km,URHO))

             compnr = compr - cdtdx*(fxz(i+1,j,km,n) - fxz(i,j,km,n)) &
                  - cdtdz*(fzx(i  ,j,kc,n) - fzx(i,j,km,n))
             compnl = compl - cdtdx*(fxz(i+1,j,km,n) - fxz(i,j,km,n)) &
                  - cdtdz*(fzx(i  ,j,kc,n) - fzx(i,j,km,n))

             qpo(i,j  ,km,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nq)
             qmo(i,j+1,km,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d,nq)

          enddo
       enddo
    enddo

    do j = jlo, jhi 
       do i = ilo, ihi 

          lambda = lam(i,j,k3d,:)

          pggxp = qx(i+1,j,km,GDPRES)
          pggxm = qx(i  ,j,km,GDPRES)
          ugxp  = qx(i+1,j,km,GDU)
          ugxm  = qx(i  ,j,km,GDU)
          ergxp = qx(i+1,j,km,GDERADS:GDERADS-1+ngroups)
          ergxm = qx(i  ,j,km,GDERADS:GDERADS-1+ngroups)

          pggzp = qz(i,j,kc,GDPRES)
          pggzm = qz(i,j,km,GDPRES)
          ugzp  = qz(i,j,kc,GDW)
          ugzm  = qz(i,j,km,GDW)
          ergzp = qz(i,j,kc,GDERADS:GDERADS-1+ngroups)
          ergzm = qz(i,j,km,GDERADS:GDERADS-1+ngroups)

          ! Convert to conservation form
          rrr = qp(i,j,km,QRHO)
          rur = rrr*qp(i,j,km,QU)
          rvr = rrr*qp(i,j,km,QV)
          rwr = rrr*qp(i,j,km,QW)
          ekenr = 0.5d0*rrr*(qp(i,j,km,QU)**2 + qp(i,j,km,QV)**2 + qp(i,j,km,QW)**2)
          rer = qp(i,j,km,QREINT) + ekenr
          err = qp(i,j,km,qrad:qradhi)

          rrl = qm(i,j+1,km,QRHO)
          rul = rrl*qm(i,j+1,km,QU)
          rvl = rrl*qm(i,j+1,km,QV)
          rwl = rrl*qm(i,j+1,km,QW)
          ekenl = 0.5d0*rrl*(qm(i,j+1,km,QU)**2 + qm(i,j+1,km,QV)**2 + qm(i,j+1,km,QW)**2)
          rel = qm(i,j+1,km,QREINT) + ekenl
          erl = qm(i,j+1,km,qrad:qradhi)

          ! Add transverse predictor
          rrnewr = rrr - cdtdx*(fxz(i+1,j,km,URHO ) - fxz(i,j,km,URHO)) &
               - cdtdz*(fzx(i  ,j,kc,URHO ) - fzx(i,j,km,URHO))
          runewr = rur - cdtdx*(fxz(i+1,j,km,UMX  ) - fxz(i,j,km,UMX)) &
               - cdtdz*(fzx(i  ,j,kc,UMX  ) - fzx(i,j,km,UMX))
          rvnewr = rvr - cdtdx*(fxz(i+1,j,km,UMY  ) - fxz(i,j,km,UMY)) &
               - cdtdz*(fzx(i  ,j,kc,UMY  ) - fzx(i,j,km,UMY))
          rwnewr = rwr - cdtdx*(fxz(i+1,j,km,UMZ  ) - fxz(i,j,km,UMZ)) &
               - cdtdz*(fzx(i  ,j,kc,UMZ  ) - fzx(i,j,km,UMZ))
          lgex = lambda(:) * (ergxp(:)-ergxm(:))
          lgez = lambda(:) * (ergzp(:)-ergzm(:))
          dmx = - cdtdx*sum(lgex)
          dmz = - cdtdz*sum(lgez)
          runewr = runewr + dmx
          rwnewr = rwnewr + dmz
          lugex = 0.5d0*(ugxp+ugxm) * lgex(:)
          lugez = 0.5d0*(ugzp+ugzm) * lgez(:)
          dre = -cdtdx * sum(lugex) - cdtdz * sum(lugez)
          renewr = rer - cdtdx*(fxz(i+1,j,km,UEDEN) - fxz(i,j,km,UEDEN)) &
               - cdtdz*(fzx(i  ,j,kc,UEDEN) - fzx(i,j,km,UEDEN)) &
               + dre

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = f1*(cdtdx*0.5d0*(ugxp+ugxm)*(ergxp(g)-ergxm(g)) &
                     +       cdtdz*0.5d0*(ugzp+ugzm)*(ergzp(g)-ergzm(g)) )
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = f1*(cdtdx*0.5d0*(ergxp(g)+ergxm(g))*(ugxm-ugxp) &
                     +       cdtdz*0.5d0*(ergzp(g)+ergzm(g))*(ugzm-ugzp) )
             end do
          else ! mixed frame
             der(:) = cdtdx*lugex + cdtdz*lugez
          end if

          ernewr = err(:) - cdtdx*(rfxz(i+1,j,km,:) - rfxz(i,j,km,:)) &
               &          - cdtdz*(rfzx(i  ,j,kc,:) - rfzx(i,j,km,:)) &
               &          + der(:)

          rrnewl = rrl - cdtdx*(fxz(i+1,j,km,URHO ) - fxz(i,j,km,URHO)) &
               - cdtdz*(fzx(i  ,j,kc,URHO ) - fzx(i,j,km,URHO))
          runewl = rul - cdtdx*(fxz(i+1,j,km,UMX  ) - fxz(i,j,km,UMX)) &
               - cdtdz*(fzx(i  ,j,kc,UMX  ) - fzx(i,j,km,UMX))
          rvnewl = rvl - cdtdx*(fxz(i+1,j,km,UMY  ) - fxz(i,j,km,UMY)) &
               - cdtdz*(fzx(i  ,j,kc,UMY  ) - fzx(i,j,km,UMY))
          rwnewl = rwl - cdtdx*(fxz(i+1,j,km,UMZ  ) - fxz(i,j,km,UMZ)) &
               - cdtdz*(fzx(i  ,j,kc,UMZ  ) - fzx(i,j,km,UMZ))
          runewl = runewl + dmx
          rwnewl = rwnewl + dmz
          renewl = rel - cdtdx*(fxz(i+1,j,km,UEDEN) - fxz(i,j,km,UEDEN)) &
               - cdtdz*(fzx(i  ,j,kc,UEDEN) - fzx(i,j,km,UEDEN)) &
               + dre
          ernewl = erl(:) - cdtdx*(rfxz(i+1,j,km,:) - rfxz(i,j,km,:)) &
               &          - cdtdz*(rfzx(i  ,j,kc,:) - rfzx(i,j,km,:)) &
               &          + der(:)

          duxp = pggxp*ugxp - pggxm*ugxm
          pxav = 0.5d0*(pggxp+pggxm)
          dux = ugxp-ugxm
          pxnew = cdtdx*(duxp + pxav*dux*(gamc(i,j,k3d)-1.d0))

          duzp = pggzp*ugzp - pggzm*ugzm
          pzav = 0.5d0*(pggzp+pggzm)
          duz = ugzp-ugzm
          pznew = cdtdz*(duzp + pzav*duz*(gamc(i,j,k3d)-1.d0))

          pnewr = qp(i,j,km,QPRES) - pxnew - pznew
          pnewl = qm(i,j+1,km,QPRES) - pxnew - pznew

          ! Convert back to non-conservation form
          rhotmp = rrnewr
          qpo(i,j,km,QRHO  ) = rhotmp        + hdt*srcQ(i,j,k3d,QRHO)
          qpo(i,j,km,QU    ) = runewr/rhotmp + hdt*srcQ(i,j,k3d,QU)
          qpo(i,j,km,QV    ) = rvnewr/rhotmp + hdt*srcQ(i,j,k3d,QV)
          qpo(i,j,km,QW    ) = rwnewr/rhotmp + hdt*srcQ(i,j,k3d,QW)
          rhoekenr = 0.5d0*(runewr**2 + rvnewr**2 + rwnewr**2)/rhotmp
          qpo(i,j,km,QREINT)= renewr - rhoekenr + hdt*srcQ(i,j,k3d,QREINT)
          qpo(i,j,km,QPRES ) = pnewr            + hdt*srcQ(i,j,k3d,QPRES)
          qpo(i,j,km,qrad:qradhi) = ernewr(:)
          qpo(i,j,km,qptot  ) = sum(lambda(:)*ernewr(:)) + qpo(i,j,km,QPRES)
          qpo(i,j,km,qreitot) = sum(qpo(i,j,km,qrad:qradhi)) + qpo(i,j,km,QREINT)

          rhotmp = rrnewl
          qmo(i,j+1,km,QRHO  ) = rhotmp        + hdt*srcQ(i,j,k3d,QRHO)
          qmo(i,j+1,km,QU    ) = runewl/rhotmp + hdt*srcQ(i,j,k3d,QU)
          qmo(i,j+1,km,QV    ) = rvnewl/rhotmp + hdt*srcQ(i,j,k3d,QV)
          qmo(i,j+1,km,QW    ) = rwnewl/rhotmp + hdt*srcQ(i,j,k3d,QW)
          rhoekenl = 0.5d0*(runewl**2 + rvnewl**2 + rwnewl**2)/rhotmp
          qmo(i,j+1,km,QREINT)= renewl - rhoekenl + hdt*srcQ(i,j,k3d,QREINT)
          qmo(i,j+1,km,QPRES ) = pnewl            + hdt*srcQ(i,j,k3d,QPRES)
          qmo(i,j+1,km,qrad:qradhi) = ernewl(:)
          qmo(i,j+1,km,qptot  ) = sum(lambda(:)*ernewl(:)) + qmo(i,j+1,km,QPRES)
          qmo(i,j+1,km,qreitot) = sum(qmo(i,j+1,km,qrad:qradhi)) + qmo(i,j+1,km,QREINT)

       enddo
    enddo

  end subroutine transxz_rad

end module transverse_rad_module
