module bc_fill_module

  use bc_ext_fill_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  include 'AMReX_bc_types.fi'

  public

contains

  subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                        domlo,domhi,delta,xlo,time,bc) &
                        bind(C, name="ca_hypfill")

    use meth_params_module, only: NVAR

    use amrex_fort_module, only : rt => amrex_real
    integer          :: adv_l1,adv_l2,adv_h1,adv_h2
    integer          :: bc(2,2,*)
    integer          :: domlo(2), domhi(2)
    real(rt)         :: delta(2), xlo(2), time
    real(rt)         :: adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)

    integer          :: n

    do n = 1,NVAR
       call filcc(adv(:,:,n),adv_l1,adv_l2,adv_h1,adv_h2, &
            domlo,domhi,delta,xlo,bc(:,:,n))
    enddo

    ! process the external BCs here
    call ext_fill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                  domlo,domhi,delta,xlo,time,bc)

  end subroutine ca_hypfill



  subroutine ca_denfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                        domlo,domhi,delta,xlo,time,bc) &
                        bind(C, name="ca_denfill")

    use amrex_fort_module, only : rt => amrex_real
    integer          :: adv_l1,adv_l2,adv_h1,adv_h2
    integer          :: bc(2,2,*)
    integer          :: domlo(2), domhi(2)
    real(rt)         :: delta(2), xlo(2), time
    real(rt)         :: adv(adv_l1:adv_h1,adv_l2:adv_h2)

    call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,bc)

    ! process the external BCs here
    call ext_denfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                     domlo,domhi,delta,xlo,time,bc)

  end subroutine ca_denfill



end module bc_fill_module
