
void LorentzFac(const amrex::Box& bx, amrex::Array4<amrex::Real const> const& vel,
                amrex::Array4<amrex::Real> const& W);

AMREX_GPU_HOST_DEVICE void ConsToPrim(amrex::Real* q_zone, amrex::Real* U_zone);

AMREX_GPU_HOST_DEVICE void PrimToCons(amrex::Real* q_zone, amrex::Real* U_zone);

AMREX_GPU_HOST_DEVICE void Flux(amrex::Real* F, const amrex::Real* q_zone,
                                const amrex::Real* U_zone, const int dir);

// alias declaration for function to be passed to BrentRootFinder
using RootFindFunc = amrex::Real (*)(amrex::Real, const amrex::Real*, bool);

AMREX_GPU_HOST_DEVICE
amrex::Real BrentRootFinder(const amrex::Real lo, const amrex::Real i, RootFindFunc func,
                            const amrex::Real* args) noexcept;

///
/// this constructs the hydrodynamic source (essentially the flux
/// divergence) using the CTU framework for unsplit hydrodynamics
///
/// @param time     current time
/// @param dt       timestep
///
void construct_ctu_rhd_source(amrex::Real time, amrex::Real dt);

void plm(const amrex::Box& bx, const amrex::Box& vbx, amrex::Array4<amrex::Real const> const& q_arr,
         amrex::Array4<amrex::Real const> const& U_arr,
         amrex::Array4<amrex::Real const> const& flatn,
         amrex::Array4<amrex::Real const> const& qaux_arr,
         amrex::Array4<amrex::Real const> const& srcQ, amrex::Array4<amrex::Real> const& dU,
         amrex::Array4<amrex::Real> const& Uxm, amrex::Array4<amrex::Real> const& Uxp,
#if AMREX_SPACEDIM >= 2
         amrex::Array4<amrex::Real> const& Uym, amrex::Array4<amrex::Real> const& Uyp,
#endif
#if AMREX_SPACEDIM == 3
         amrex::Array4<amrex::Real> const& Uzm, amrex::Array4<amrex::Real> const& Uzp,
#endif
#if AMREX_SPACEDIM < 3
         amrex::Array4<amrex::Real const> const& dloga,
#endif
         const amrex::Real dt, const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx);

void consup_rhd(const amrex::Box& bx, amrex::Array4<amrex::Real> const& update,
                amrex::Array4<amrex::Real> const& flux0,
                amrex::Array4<amrex::Real const> const& area0,
#if AMREX_SPACEDIM >= 2
                amrex::Array4<amrex::Real> const& flux1,
                amrex::Array4<amrex::Real const> const& area1,
#endif
#if AMREX_SPACEDIM == 3
                amrex::Array4<amrex::Real> const& flux2,
                amrex::Array4<amrex::Real const> const& area2,
#endif
                amrex::Array4<amrex::Real const> const& vol);

void hlle(const amrex::Box& bx, amrex::Array4<amrex::Real const> const& UL,
          amrex::Array4<amrex::Real const> const& UR, amrex::Array4<amrex::Real> const& F,
          amrex::Array4<amrex::Real> const& Uint, amrex::Array4<amrex::Real> const& Ugdnv,
          const int dir, const amrex::Real dt);