module bc_fill_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains

  subroutine ca_hypfill(adv,adv_l1,adv_h1, &
                        domlo,domhi,delta,xlo,time,bc) &
                        bind(C, name="ca_hypfill")

    use meth_params_module, only: NVAR

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: adv_l1,adv_h1
    integer          :: bc(1,2,*)
    integer          :: domlo(1), domhi(1)
    real(rt)         :: delta(1), xlo(1), time
    real(rt)         :: adv(adv_l1:adv_h1,NVAR)

    integer          :: n

    do n = 1,NVAR
       call filcc(adv(:,n),adv_l1,adv_h1, &
                  domlo,domhi,delta,xlo,bc(:,:,n))
    enddo

  end subroutine ca_hypfill



  subroutine ca_denfill(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="ca_denfill")

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: adv_l1,adv_h1
    integer          :: bc(1,2,*)
    integer          :: domlo(1), domhi(1)
    real(rt)         :: delta(1), xlo(1), time
    real(rt)         :: adv(adv_l1:adv_h1)

    call filcc(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,bc)

  end subroutine ca_denfill


end module bc_fill_module
