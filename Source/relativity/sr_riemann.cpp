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

        // AllPrint() << "pL = " << qL_cell[QPRES] << ", pR = " << qL_cell[QPRES] << std::endl;

        // Real pL = qL_cell(i,j,k,QPRES);
        // Real pR = qR_cell(i,j,k,QPRES);

        // define the signal speeds
        // assume square grids
        Real aL = -1.0_rt;  // * amrex::min(dx[0] / dt, 1.0_rt);
        Real aR = 1.0_rt;   //amrex::min(dx[0] / dt, 1.0_rt);

        Real FL[NUM_STATE];
        Real FR[NUM_STATE];
        Real F_HLLE[NUM_STATE];
        Real U_HLLE[NUM_STATE];

        // get the fluxes
        Flux(FL, qL_cell, UL_cell, dir);
        Flux(FR, qR_cell, UR_cell, dir);

        for (int n = 0; n < NUM_STATE; ++n) {
            F_HLLE[n] = (aR * FL[n] - aL * FR[n] + aR * aL * (UR_cell[n] - UL_cell[n])) / (aR - aL);

            if (aR <= 0.0_rt) {  // right state
                U_HLLE[n] = UR_cell[n];
            } else if (aL < 0.0_rt) {  // middle
                U_HLLE[n] = (aR * UR_cell[n] - aL * UL_cell[n] - FR[n] + FL[n]) / (aR - aL);
            } else {  // left state
                U_HLLE[n] = UL_cell[n];
            }
        }

        // shock momentum and flux
        Real S_HLLE = U_HLLE[UMX + dir];
        Real F_S = F_HLLE[UMX + dir];

        Real E_HLLE = U_HLLE[UEDEN] + U_HLLE[URHO];

        Real a_star;
        if (Math::abs(F_HLLE[UEDEN]) < 1.e-9_rt) {
            a_star = S_HLLE / (E_HLLE + F_S);
        } else {
            a_star =
                (E_HLLE + F_S -
                 std::sqrt((E_HLLE + F_S) * (E_HLLE + F_S) - S_HLLE * 2.0_rt * F_HLLE[UEDEN])) /
                (2.0_rt * F_HLLE[UEDEN]);
        }

        if (i == 65 || i == 64) {

            AllPrint() << i << " pL = " << qL_cell[QPRES] << ", pR = " << qL_cell[QPRES]
                       << ", a_star = " << a_star << std::endl;

            AllPrint() << "i = " << i << ": " << qL_cell[QRHO] << ' ' << qR_cell[QRHO]
                       << ", UEDEN = " << UL_cell[UEDEN] << ", UMX = " << UL_cell[UMX]
                       << ", FL = " << FL[UMX] << ", FR = " << FR[UMX] << std::endl;
        }

        // nan check!
        if (a_star != a_star) {
            Abort("Nan! a_star has nan'd");
        }

        // left star state
        Real A = (UL_cell[UEDEN] + UL_cell[URHO]) * aL - UL_cell[UMX + dir];
        Real B = UL_cell[UMX + dir] * (aL - qL_cell[QU + dir]) - qL_cell[QPRES];

        Real pL_star = ((A * a_star) - B) / (1.0_rt - aL * a_star);

        Real UL_star[NUM_STATE];

        for (int n = 0; n < NUM_STATE; ++n) {
            UL_star[n] = UL_cell[n] * (aL - qL_cell[QU + dir]) / (aL - a_star);
        }

        UL_star[UMX + dir] += (pL_star - qL_cell[QPRES]) / (aL - a_star);

        // right star state
        A = (UR_cell[UEDEN] + UR_cell[URHO]) * aR - UR_cell[UMX + dir];
        B = UR_cell[UMX + dir] * (aR - qR_cell[QU + dir]) - qR_cell[QPRES];

        Real pR_star = ((A * a_star) - B) / (1.0_rt - aR * a_star);

        Real UR_star[NUM_STATE];

        for (int n = 0; n < NUM_STATE; ++n) {
            UR_star[n] = UR_cell[n] * (aR - qR_cell[QU + dir]) / (aR - a_star);
        }

        UR_star[UMX + dir] += (pR_star - qR_cell[QPRES]) / (aR - a_star);

        for (int n = 0; n < NUM_STATE; ++n) {
            if (aR <= 0.0_rt) {  // right state
                F(i, j, k, n) = FR[n];
            } else if (a_star <= 0.0_rt) {  // right star state
                F(i, j, k, n) = UR_star[n] * a_star;

                if (n == UMX + dir) {
                    F(i, j, k, n) += pR_star;
                } else if (n == UEDEN) {
                    F(i, j, k, n) = UR_star[UMX + dir] - UR_star[URHO] * a_star;
                }
            } else if (a_star < 0.0_rt) {  // left star state
                F(i, j, k, n) = UL_star[n] * a_star;

                if (n == UMX + dir) {
                    F(i, j, k, n) += pL_star;
                } else if (n == UEDEN) {
                    F(i, j, k, n) = UL_star[UMX + dir] - UL_star[URHO] * a_star;
                }
            } else {  // left state
                F(i, j, k, n) = FL[n];
            }
        }

        if (i == 65 || i == 64) {
            AllPrint() << "i = " << i << ": "
                       << "F = ";
            for (int n = 0; n < NUM_STATE; ++n) {
                AllPrint() << F(i, j, k, n) << ", ";
            }
            AllPrint() << std::endl;
        }

        // calculate numerical flux vector
        // Real aL_m = amrex::min(0.0_rt, aL);
        // Real aR_p = amrex::max(0.0_rt, aR);

        // for (int n = 0; n < NUM_STATE; ++n) {
        //     F(i, j, k, n) =
        //         (aR_p * FL[n] - aL_m * FR[n] + aR_p * aL_m * (UR_cell[n] - UL_cell[n])) /
        //         (aR_p - aL_m);

        //     if (dx[dir] / dt < aL) {
        //         Uint(i, j, k, n) = UL(i, j, k, n);
        //     } else if (dx[dir] / dt < aR) {
        //         Uint(i, j, k, n) = (aR * UR_cell[n] - aL * UL_cell[n] - FR[n] + FL[n]) / (aR - aL);
        //     } else {
        //         Uint(i, j, k, n) = UR(i, j, k, n);
        //     }
        // }
    });

    // store_godunov_state(bx, qint, qgdnv);
}