void LorentzFac(const amrex::Box& bx, amrex::Array4<amrex::Real const> const& vel,
                amrex::Array4<amrex::Real> const& W);

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void ConsToPrim(amrex::Real* q_zone, amrex::Real* U_zone);

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void PrimToCons(amrex::Real* q_zone, amrex::Real* U_zone);

// alias declaration for function to be passed to BrentRootFinder
using RootFindFunc = amrex::Real (*)(amrex::Real, const amrex::Real*);

AMREX_GPU_HOST_DEVICE
amrex::Real BrentRootFinder(const amrex::Real lo, const amrex::Real i, RootFindFunc func, const amrex::Real* args) noexcept;