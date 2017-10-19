module c_interface_modules

  use meth_params_module, only: NVAR, NQAUX, NQ, QVAR
  use amrex_fort_module, only: rt => amrex_real

#ifdef CUDA
    use cudafor, only: cudaMemcpyAsync, cudaMemcpyHostToDevice, &
                       cudaDeviceSynchronize, dim3, cuda_stream_kind
    use cuda_module, only: threads_and_blocks, cuda_streams, max_cuda_streams
#endif

contains

    subroutine ca_enforce_consistent_e(lo,hi,state,s_lo,s_hi,idx) &
       bind(c, name='ca_enforce_consistent_e')

    use castro_util_module, only: enforce_consistent_e

    implicit none

    integer, intent(in)     :: lo(3), hi(3)
    integer, intent(in)     :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    integer, intent(in)     :: idx

    call enforce_consistent_e(lo, hi, state, s_lo, s_hi)

end subroutine ca_enforce_consistent_e

subroutine ca_compute_temp(lo, hi, state, s_lo, s_hi, idx) &
       bind(C, name="ca_compute_temp")

    use castro_util_module, only: compute_temp

    implicit none

    integer, intent(in   ) :: lo(3),hi(3)
    integer, intent(in   ) :: s_lo(3),s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    integer, intent(in)     :: idx

    call compute_temp(lo, hi, state, s_lo, s_hi)

end subroutine ca_compute_temp

subroutine ca_reset_internal_e(lo, hi, u, u_lo, u_hi, verbose, idx) &
       bind(C, name="ca_reset_internal_e")

    use castro_util_module, only: reset_internal_e

    implicit none

    integer, intent(in) :: lo(3), hi(3), verbose
    integer, intent(in) :: u_lo(3), u_hi(3)
    real(rt), intent(inout) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
    integer, intent(in)     :: idx

    call reset_internal_e(lo, hi, u, u_lo, u_hi, verbose)

end subroutine ca_reset_internal_e


  subroutine ca_normalize_species(u, u_lo, u_hi, lo, hi, idx) &
       bind(C, name="ca_normalize_species")

    use castro_util_module, only: normalize_species

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: u_lo(3), u_hi(3)
    real(rt), intent(inout) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
    integer, intent(in)     :: idx

    call normalize_species(u, u_lo, u_hi, lo, hi)

  end subroutine ca_normalize_species

  subroutine ca_swe_to_comp(swe, slo, shi, comp, clo, chi, lo, hi) &
       bind(C, name="ca_swe_to_comp")

    use actual_riemann_module, only: swe_to_comp
    use meth_params_module, only: NVAR

    implicit none

    integer, intent(in)   :: slo(3), shi(3), clo(3), chi(3), lo(3), hi(3)
    real(rt), intent(inout)  :: swe(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3), NVAR)
    real(rt), intent(inout) :: comp(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3), NVAR)

    call swe_to_comp(swe, slo, shi, comp, clo, chi, lo, hi)

  end subroutine ca_swe_to_comp

  subroutine ca_comp_to_swe(swe, slo, shi, comp, clo, chi, lo, hi) &
       bind(C, name="ca_comp_to_swe")

    use actual_riemann_module, only: comp_to_swe
    use meth_params_module, only: NVAR

    implicit none

    integer, intent(in)   :: slo(3), shi(3), clo(3), chi(3), lo(3), hi(3)
    real(rt), intent(out)  :: swe(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3), NVAR)
    real(rt), intent(in) :: comp(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3), NVAR)

    call comp_to_swe(swe, slo, shi, comp, clo, chi, lo, hi)

  end subroutine ca_comp_to_swe


  subroutine ca_enforce_minimum_density(uin, uin_lo, uin_hi, &
                                        uout, uout_lo, uout_hi, &
                                        vol, vol_lo, vol_hi, &
                                        lo, hi, frac_change, verbose, idx) &
                                        bind(C, name="ca_enforce_minimum_density")

    use advection_util_module, only: enforce_minimum_density

    implicit none

    integer, intent(in) :: lo(3), hi(3), verbose
    integer, intent(in) ::  uin_lo(3),  uin_hi(3)
    integer, intent(in) :: uout_lo(3), uout_hi(3)
    integer, intent(in) ::  vol_lo(3),  vol_hi(3)

    real(rt), intent(inout) ::  uin( uin_lo(1): uin_hi(1), uin_lo(2): uin_hi(2), uin_lo(3): uin_hi(3),NVAR)
    real(rt), intent(inout) :: uout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2),uout_lo(3):uout_hi(3),NVAR)
    real(rt), intent(in) ::  vol( vol_lo(1): vol_hi(1), vol_lo(2): vol_hi(2), vol_lo(3): vol_hi(3))
    real(rt), intent(inout) :: frac_change
    integer, intent(in)     :: idx

    call enforce_minimum_density(uin, uin_lo, uin_hi, &
                                 uout, uout_lo, uout_hi, &
                                 vol, vol_lo, vol_hi, &
                                 lo, hi, frac_change, verbose)

  end subroutine ca_enforce_minimum_density


  subroutine ca_check_initial_species(lo, hi, state, state_lo, state_hi, idx) &
                                      bind(C, name="ca_check_initial_species")

    use castro_util_module, only: check_initial_species

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: state_lo(3), state_hi(3)
    real(rt), intent(in) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)
    integer, intent(in)     :: idx

    call check_initial_species(lo, hi, state, state_lo, state_hi)

  end subroutine ca_check_initial_species


  subroutine ca_ctoprim(lo, hi, &
                        uin, uin_lo, uin_hi, &
                        q,     q_lo,   q_hi, &
                        qaux, qa_lo,  qa_hi, idx) bind(C, name = "ca_ctoprim")

    use advection_util_module, only: swectoprim

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: uin_lo(3), uin_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)

    real(rt)        , intent(inout) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)

    real(rt)        , intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt)        , intent(inout) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    integer, intent(in)     :: idx

    call swectoprim(lo, hi, &
                 uin, uin_lo, uin_hi, &
                 q,     q_lo,   q_hi, &
                 qaux, qa_lo,  qa_hi)

  end subroutine ca_ctoprim


end module c_interface_modules
