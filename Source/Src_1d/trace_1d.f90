module trace_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  private

  public trace

contains

  ! :::
  ! ::: ------------------------------------------------------------------
  ! :::
  subroutine trace(q, dq, flatn, q_lo, q_hi, &
                   qaux, qa_lo, qa_hi, &
                   dloga, dloga_lo, dloga_hi, &
                   srcQ, src_lo, src_hi, &
                   qxm, qxp, qpd_lo, qpd_hi, &
                   ilo, ihi, domlo, domhi, dx, dt)

    use meth_params_module, only : plm_iorder, QVAR, NQ, NQAUX, QRHO, QU, QREINT, QC, QPRES, &
         npassive, qpass_map, small_dens, ppm_type, fix_mass_flux, use_pslope
    use prob_params_module, only : physbc_lo, physbc_hi, Outflow
    use slope_module, only : uslope, pslope
    use bl_constants_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer domlo(1), domhi(1)
    integer ilo,ihi
    integer    q_lo(3), q_hi(3)
    integer    qa_lo(3), qa_hi(3)
    integer dloga_lo(3), dloga_hi(3)
    integer   qpd_lo(3),  qpd_hi(3)
    integer   src_lo(3),  src_hi(3)
    real(rt)         dx, dt
    real(rt)             q(q_lo(1):q_hi(1),NQ)
    real(rt)             qaux(qa_lo(1):qa_hi(1),NQAUX)
    real(rt)          srcQ(src_lo(1):src_hi(1),QVAR)
    real(rt)         flatn(q_lo(1):q_hi(1))
    real(rt)         dloga(dloga_lo(1):dloga_hi(1))

    real(rt)           dq( qpd_lo(1): qpd_hi(1),NQ)
    real(rt)          qxm( qpd_lo(1): qpd_hi(1),NQ)
    real(rt)          qxp( qpd_lo(1): qpd_hi(1),NQ)

    !     Local variables
    integer          :: i
    integer          :: n, ipassive

    real(rt)         :: hdt,dtdx
    real(rt)         :: cc, csq, rho, u, p, rhoe
    real(rt)         :: drho, du, dp, drhoe

    real(rt)         :: enth, alpham, alphap, alpha0r, alpha0e
    real(rt)         :: spminus, spplus, spzero
    real(rt)         :: apright, amright, azrright, azeright
    real(rt)         :: apleft, amleft, azrleft, azeleft
    real(rt)         :: acmprght, acmpleft
    real(rt)         :: sourcr,sourcp,source,courn,eta,dlogatmp

    logical :: fix_mass_flux_lo, fix_mass_flux_hi

    fix_mass_flux_lo = &
         (fix_mass_flux .eq. 1) .and. (physbc_lo(1) .eq. Outflow) .and. (ilo .eq. domlo(1))
    fix_mass_flux_hi = &
         (fix_mass_flux .eq. 1) .and. (physbc_hi(1) .eq. Outflow) .and. (ihi .eq. domhi(1))

    if (ppm_type .ne. 0) then
       print *,'Oops -- shouldnt be in trace with ppm_type != 0'
       call bl_error("Error:: Castro_1d.f90 :: trace")
    end if

    hdt = HALF * dt
    dtdx = dt/dx

    ! Compute slopes
    if (plm_iorder .eq. 1) then

       dq(ilo-1:ihi+1,1:NQ) = ZERO

    else

       call uslope(q, flatn, q_lo, q_hi, &
                   dq, dq, dq, qpd_lo, qpd_hi, &
                   ilo, 0, ihi, 0, 0, 0)

       if (use_pslope .eq. 1) &
            call pslope(q, flatn, q_lo, q_hi, &
                        dq, dq, dq, qpd_lo, qpd_hi, &
                        srcQ, src_lo, src_hi, &
                        ilo, 0, ihi, 0, 0, 0, [dx, ZERO, ZERO])

    endif

    ! Compute left and right traced states
    do i = ilo-1, ihi+1

       cc = qaux(i,QC)
       csq = cc**2
       rho = q(i,QRHO)
       u = q(i,QU)
       p = q(i,QPRES)
       rhoe = q(i,QREINT)
       enth = ( (rhoe+p)/rho )/csq

       drho  = dq(i,QRHO)
       du    = dq(i,QU)
       dp    = dq(i,QPRES)
       drhoe = dq(i,QREINT)

       alpham = HALF*(dp/(rho*cc) - du)*rho/cc
       alphap = HALF*(dp/(rho*cc) + du)*rho/cc
       alpha0r = drho - dp/csq
       alpha0e = drhoe - dp*enth

       if (u-cc .gt. ZERO) then
          spminus = -ONE
       else
          spminus = (u-cc)*dtdx
       endif
       if (u+cc .gt. ZERO) then
          spplus = -ONE
       else
          spplus = (u+cc)*dtdx
       endif
       if (u .gt. ZERO) then
          spzero = -ONE
       else
          spzero = u*dtdx
       endif

       apright = HALF*(-ONE - spplus )*alphap
       amright = HALF*(-ONE - spminus)*alpham
       azrright = HALF*(-ONE - spzero )*alpha0r
       azeright = HALF*(-ONE - spzero )*alpha0e

       if (i .ge. ilo) then
          qxp(i,QRHO  ) = rho + apright + amright + azrright
          qxp(i,QRHO  ) = max(small_dens,qxp(i,QRHO))
          qxp(i,QU    ) = u + (apright - amright)*cc/rho
          qxp(i,QPRES ) = p + (apright + amright)*csq
          qxp(i,QREINT) = rhoe + (apright + amright)*enth*csq + azeright

          ! add source term
          qxp(i  ,QRHO  ) = qxp(i,QRHO  ) + hdt*srcQ(i,QRHO)
          qxp(i  ,QRHO  ) = max(small_dens,qxp(i,QRHO))
          qxp(i  ,QU    ) = qxp(i,QU    ) + hdt*srcQ(i,QU)
          qxp(i  ,QREINT) = qxp(i,QREINT) + hdt*srcQ(i,QREINT)
          qxp(i  ,QPRES ) = qxp(i,QPRES ) + hdt*srcQ(i,QPRES)
       end if

       if (u-cc .ge. ZERO) then
          spminus = (u-cc)*dtdx
       else
          spminus = ONE
       endif
       if (u+cc .ge. ZERO) then
          spplus = (u+cc)*dtdx
       else
          spplus = ONE
       endif
       if (u .ge. ZERO) then
          spzero = u*dtdx
       else
          spzero = ONE
       endif

       apleft = HALF*(ONE - spplus )*alphap
       amleft = HALF*(ONE - spminus)*alpham
       azrleft = HALF*(ONE - spzero )*alpha0r
       azeleft = HALF*(ONE - spzero )*alpha0e

       if (i .le. ihi) then
          qxm(i+1,QRHO  ) = rho + apleft + amleft + azrleft
          qxm(i+1,QRHO  ) = max(small_dens, qxm(i+1,QRHO))
          qxm(i+1,QU    ) = u + (apleft - amleft)*cc/rho
          qxm(i+1,QPRES ) = p + (apleft + amleft)*csq
          qxm(i+1,QREINT) = rhoe + (apleft + amleft)*enth*csq + azeleft

          ! add source terms
          qxm(i+1,QRHO  ) = qxm(i+1,QRHO  ) + hdt*srcQ(i,QRHO)
          qxm(i+1,QRHO  ) = max(small_dens, qxm(i+1,QRHO))
          qxm(i+1,QU    ) = qxm(i+1,QU    ) + hdt*srcQ(i,QU)
          qxm(i+1,QREINT) = qxm(i+1,QREINT) + hdt*srcQ(i,QREINT)
          qxm(i+1,QPRES ) = qxm(i+1,QPRES ) + hdt*srcQ(i,QPRES)
       end if

       if(dloga(i).ne.0)then
          courn = dtdx*(cc+abs(u))
          eta = (ONE-courn)/(cc*dt*abs(dloga(i)))
          dlogatmp = min(eta,ONE)*dloga(i)
          sourcr = -HALF*dt*rho*dlogatmp*u
          sourcp = sourcr*csq
          source = sourcp*enth
          if (i .le. ihi) then
             qxm(i+1,QRHO  ) = qxm(i+1,QRHO  ) + sourcr
             qxm(i+1,QRHO  ) = max(small_dens, qxm(i+1,QRHO))
             qxm(i+1,QPRES ) = qxm(i+1,QPRES ) + sourcp
             qxm(i+1,QREINT) = qxm(i+1,QREINT) + source
          end if
          if (i .ge. ilo) then
             qxp(i  ,QRHO  ) = qxp(i  ,QRHO  ) + sourcr
             qxp(i  ,QRHO  ) = max(small_dens,qxp(i,QRHO))
             qxp(i  ,QPRES ) = qxp(i  ,QPRES ) + sourcp
             qxp(i  ,QREINT) = qxp(i  ,QREINT) + source
          end if
       endif
    enddo

    ! Enforce constant mass flux rate if specified
    if (fix_mass_flux_lo) then
       qxm(ilo,QRHO  ) = q(domlo(1)-1,QRHO)
       qxm(ilo,QU    ) = q(domlo(1)-1,QU  )
       qxm(ilo,QPRES ) = q(domlo(1)-1,QPRES)
       qxm(ilo,QREINT) = q(domlo(1)-1,QREINT)
    end if

    ! Enforce constant mass flux rate if specified
    if (fix_mass_flux_hi) then
       qxp(ihi+1,QRHO  ) = q(domhi(1)+1,QRHO)
       qxp(ihi+1,QU    ) = q(domhi(1)+1,QU  )
       qxp(ihi+1,QPRES ) = q(domhi(1)+1,QPRES)
       qxp(ihi+1,QREINT) = q(domhi(1)+1,QREINT)
    end if

    do ipassive = 1, npassive
       n = qpass_map(ipassive)

       ! Right state
       do i = ilo,ihi+1
          u = q(i,QU)
          if (u .gt. ZERO) then
             spzero = -ONE
          else
             spzero = u*dtdx
          endif
          acmprght = HALF*(-ONE - spzero )*dq(i,n)
          qxp(i,n) = q(i,n) + acmprght
       enddo
       if (fix_mass_flux_hi) qxp(ihi+1,n) = q(ihi+1,n)

       ! Left state
       do i = ilo-1,ihi
          u = q(i,QU)
          if (u .ge. ZERO) then
             spzero = u*dtdx
          else
             spzero = ONE
          endif
          acmpleft = HALF*(ONE - spzero )*dq(i,n)
          qxm(i+1,n) = q(i,n) + acmpleft
       enddo
       if (fix_mass_flux_lo) qxm(ilo,n) = q(ilo-1,n)
    enddo

  end subroutine trace

end module trace_module
