#include "Castro.H"
#include "Castro_F.H"

// #include "sr_util.H"

using namespace amrex;

void Castro::hlle(const Box& bx, Array4<Real const> const& UL, Array4<Real const> const& UR,
                  Array4<Real> const& F, Array4<Real> const& Uint, Array4<Real> const& qgdnv,
                  const int dir, const Real dt) {

    const auto dx = geom.CellSizeArray();

    amrex::ParallelFor(bx, [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept {
        Real qL_cell[NQ];
        Real qR_cell[NQ];

        Real UL_cell[NUM_STATE];
        Real UR_cell[NUM_STATE];

        for (auto n = 0; n < NUM_STATE; ++n) {
            UL_cell[n] = UL(i, j, k, n);
            UR_cell[n] = UR(i, j, k, n);
        }

        // get the primitive states
        ConsToPrim(qL_cell, UL_cell);
        ConsToPrim(qR_cell, UR_cell);

        // define the signal speeds
        // assume square grids
        Real aL = -1.0_rt * amrex::min(dx[0] / dt, 1.0_rt);
        Real aR = amrex::min(dx[0] / dt, 1.0_rt);

        Real FL[NUM_STATE];
        Real FR[NUM_STATE];

        // get the fluxes
        Flux(FL, qL_cell, UL_cell, dir);
        Flux(FR, qR_cell, UR_cell, dir);

        // calculate numerical flux vector
        Real aL_m = amrex::min(0.0_rt, aL);
        Real aR_p = amrex::max(0.0_rt, aR);

        for (int n = 0; n < NUM_STATE; ++n) {
            F(i, j, k, n) =
                (aR_p * FL[n] - aL_m * FR[n] + aR_p * aL_m * (UR_cell[n] - UL_cell[n])) /
                (aR_p - aL_m);

            Uint(i, j, k, n) = (aR * UR_cell[n] - aL * UL_cell[n] - FR[n] + FL[n]) / (aR - aL);
        }
    });

    // store_godunov_state(bx, qint, qgdnv);
}