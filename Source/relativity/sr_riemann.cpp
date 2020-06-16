#include "Castro.H"
#include "Castro_F.H"

// #include "sr_util.H"

using namespace amrex;

void Castro::hlle(const Box& bx, Array4<Real const> const& qL, Array4<Real const> const& qR,
                  Array4<Real> const& F, Array4<Real> const& qint, Array4<Real> const& qgdnv,
                  const int dir, const Real dt) {

    const auto dx = geom.CellSizeArray();

    amrex::ParallelFor(bx, [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept {
        Real qL_cell[NQ];
        Real qR_cell[NQ];

        Real UL[NUM_STATE];
        Real UR[NUM_STATE];

        Real Uint[NUM_STATE];
        Real qint_cell[NQ];

        for (auto n = 0; n < NQ; ++n) {
            qL_cell[n] = qL(i, j, k, n);
            qR_cell[n] = qR(i, j, k, n);
        }

        // get the conserved states
        PrimToCons(qL_cell, UL);
        PrimToCons(qR_cell, UR);

        // define the signal speeds
        // assume square grids
        Real aL = -1.0_rt * amrex::min(dx[0] / dt, 1.0_rt);
        Real aR = amrex::min(dx[0] / dt, 1.0_rt);

        Real FL[NUM_STATE];
        Real FR[NUM_STATE];

        // get the fluxes
        Flux(FL, qL_cell, UL, dir);
        Flux(FR, qR_cell, UR, dir);

        // calculate numerical flux vector
        Real aL_m = amrex::min(0.0_rt, aL);
        Real aR_p = amrex::max(0.0_rt, aR);

        for (int n = 0; n < NUM_STATE; ++n) {
            F(i, j, k, n) =
                (aR_p * FL[n] - aL_m * FR[n] + aR_p * aL_m * (UR[n] - UL[n])) / (aR_p - aL_m);

            Uint[n] = (aR * UR[n] - aL * UL[n] - FR[n] + FL[n]) / (aR - aL);
        }

        ConsToPrim(qint_cell, Uint);

        for (int n = 0; n < NQ; ++n) {
            qint(i, j, k, n) = qint_cell[n];
        }
    });

    store_godunov_state(bx, qint, qgdnv);
}