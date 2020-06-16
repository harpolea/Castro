
void LorentzFac(const amrex::Box& bx, amrex::Array4<amrex::Real const> const& vel,
                amrex::Array4<amrex::Real> const& W);

AMREX_GPU_HOST_DEVICE void ConsToPrim(amrex::Real* q_zone, amrex::Real* U_zone);

AMREX_GPU_HOST_DEVICE void PrimToCons(amrex::Real* q_zone, amrex::Real* U_zone);

AMREX_GPU_HOST_DEVICE void Flux(amrex::Real* F, const amrex::Real* q_zone,
                                const amrex::Real* U_zone, const int dir);

// alias declaration for function to be passed to BrentRootFinder
using RootFindFunc = amrex::Real (*)(amrex::Real, const amrex::Real*);

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
         amrex::Array4<amrex::Real const> const& flatn,
         amrex::Array4<amrex::Real const> const& qaux_arr,
         amrex::Array4<amrex::Real const> const& srcQ, amrex::Array4<amrex::Real> const& dq,
         amrex::Array4<amrex::Real> const& qxm, amrex::Array4<amrex::Real> const& qxp,
#if AMREX_SPACEDIM >= 2
         amrex::Array4<amrex::Real> const& qym, amrex::Array4<amrex::Real> const& qyp,
#endif
#if AMREX_SPACEDIM == 3
         amrex::Array4<amrex::Real> const& qzm, amrex::Array4<amrex::Real> const& qzp,
#endif
#if AMREX_SPACEDIM < 3
         amrex::Array4<amrex::Real const> const& dloga,
#endif
         const amrex::Real dt);

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

void hlle(const amrex::Box& bx, amrex::Array4<amrex::Real const> const& qL,
          amrex::Array4<amrex::Real const> const& qR, amrex::Array4<amrex::Real> const& F,
          amrex::Array4<amrex::Real> const& qint, amrex::Array4<amrex::Real> const& qgdnv,
          const int dir, const amrex::Real dt);