#include "Castro.H"
#include "Castro_F.H"
#include "Castro_hydro.H"
#include "Castro_sr.H"
#include "Castro_util.H"
#include "sr_slope.H"

using namespace amrex;

void Castro::construct_ctu_rhd_source(Real time, Real dt) {

    BL_PROFILE("Castro::construct_ctu_rhd_source()");

    const Real strt_time = ParallelDescriptor::second();

    // this constructs the hydrodynamic source (essentially the flux
    // divergence) using the CTU framework for unsplit hydrodynamics

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... Entering construct_ctu_rhd_source()" << std::endl << std::endl;

    hydro_source.setVal(0.0);

    int coord = geom.Coord();

    const auto dx = geom.CellSizeArray();

    GpuArray<Real, 3> center;
    ca_get_center(center.begin());

    MultiFab& S_new = get_new_data(State_Type);

    Real mass_lost = 0.;
    Real xmom_lost = 0.;
    Real ymom_lost = 0.;
    Real zmom_lost = 0.;
    Real eden_lost = 0.;
    Real xang_lost = 0.;
    Real yang_lost = 0.;
    Real zang_lost = 0.;

#ifdef _OPENMP
#pragma omp parallel reduction(max:nstep_fsp) \
                     reduction(+:mass_lost,xmom_lost,ymom_lost,zmom_lost) \
                     reduction(+:eden_lost,xang_lost,yang_lost,zang_lost)
#else
#pragma omp parallel reduction(+:mass_lost,xmom_lost,ymom_lost,zmom_lost) \
                     reduction(+:eden_lost,xang_lost,yang_lost,zang_lost)
#endif
    {

        // Declare local storage now. This should be done outside the MFIter loop,
        // and then we will resize the Fabs in each MFIter loop iteration. Then,
        // we apply an Elixir to ensure that their memory is saved until it is no
        // longer needed (only relevant for the asynchronous case, usually on GPUs).

        FArrayBox q;
        FArrayBox flatn;
        FArrayBox dq;
        FArrayBox src_q;
        FArrayBox Uxm, Uxp;
        // #if AMREX_SPACEDIM >= 2
        //         FArrayBox Uym, Uyp;
        // #endif
        // #if AMREX_SPACEDIM == 3
        //         FArrayBox Uzm, Uzp;
        // #endif
        FArrayBox div;
        FArrayBox U_int;
        // #if AMREX_SPACEDIM >= 2
        //         FArrayBox ftmp1, ftmp2;
        //         FArrayBox qgdnvtmp1, qgdnvtmp2;
        //         FArrayBox ql, qr;
        // #endif
        FArrayBox flux[AMREX_SPACEDIM], qe[AMREX_SPACEDIM];
#if AMREX_SPACEDIM <= 2
        FArrayBox pradial;
#endif
        // #if AMREX_SPACEDIM == 3
        //         FArrayBox qmyx, qpyx;
        //         FArrayBox qmzx, qpzx;
        //         FArrayBox qmxy, qpxy;
        //         FArrayBox qmzy, qpzy;
        //         FArrayBox qmxz, qpxz;
        //         FArrayBox qmyz, qpyz;
        // #endif

#ifdef AMREX_USE_GPU
        size_t starting_size = MultiFab::queryMemUsage("AmrLevel_Level_" + std::to_string(level));
        size_t current_size = starting_size;
#endif

        for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

            size_t fab_size = 0;

            // the valid region box
            const Box& bx = mfi.tilebox();

            const Box& obx = amrex::grow(bx, 1);
            const Box& bx_gc = amrex::grow(bx, NUM_GROW);

            flatn.resize(bx_gc, NQ);
            Elixir elix_flatn = flatn.elixir();
            fab_size += flatn.nBytes();

            // If we are oversubscribing the GPU, performance of the hydro will be constrained
            // due to its heavy memory requirements. We can help the situation by prefetching in
            // all the data we will need, and then prefetching it out at the end. This at least
            // improves performance by mitigating the number of unified memory page faults.

            // Unfortunately in CUDA there is no easy way to see actual current memory usage when
            // using unified memory; querying CUDA for free memory usage will only tell us whether
            // we've oversubscribed at any point, not whether we're currently oversubscribing, but
            // this is still a good heuristic in most cases.

            bool oversubscribed = false;

#ifdef AMREX_USE_CUDA
            if (Gpu::Device::freeMemAvailable() < 0.005 * Gpu::Device::totalGlobalMem()) {
                oversubscribed = true;
            }
#endif

            if (oversubscribed) {
                // q[mfi].prefetchToDevice();
                qaux[mfi].prefetchToDevice();
                volume[mfi].prefetchToDevice();
                Sborder[mfi].prefetchToDevice();
                hydro_source[mfi].prefetchToDevice();
                for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                    area[i][mfi].prefetchToDevice();
                    (*fluxes[i])[mfi].prefetchToDevice();
                }
#if AMREX_SPACEDIM < 3
                dLogArea[0][mfi].prefetchToDevice();
                P_radial[mfi].prefetchToDevice();
#endif
            }

            // Array4<Real const> const q_arr = q.array();
            Array4<Real const> const U_arr = Sborder.array(mfi);
            Array4<Real const> const qaux_arr = qaux.array(mfi);

            Array4<Real const> const areax_arr = area[0].array(mfi);
#if AMREX_SPACEDIM >= 2
            Array4<Real const> const areay_arr = area[1].array(mfi);
#endif
#if AMREX_SPACEDIM == 3
            Array4<Real const> const areaz_arr = area[2].array(mfi);
#endif

            Array4<Real> const vol_arr = volume.array(mfi);

#if AMREX_SPACEDIM < 3
            Array4<Real const> const dLogArea_arr = (dLogArea[0]).array(mfi);
#endif

            // Calculate primitives based on conservatives
            q.resize(bx_gc, NQ);
            auto q_arr = q.array();
            auto elix_q = q.elixir();

            // Print() << "hello constoprim" << std::endl;

            amrex::ParallelFor(obx, [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept {
                Real q_zone[NQ];
                Real U_zone[NUM_STATE];
                for (auto n = 0; n < NUM_STATE; ++n) {
                    U_zone[n] = U_arr(i, j, k, n);
                }

                ConsToPrim(q_zone, U_zone);

                for (auto n = 0; n < NQ; ++n) {
                    q_arr(i, j, k, n) = q_zone[n];
                }
            });

            // compute the flattening coefficient

            Array4<Real> const flatn_arr = flatn.array();

            if (first_order_hydro == 1) {
                amrex::ParallelFor(obx, [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept {
                    flatn_arr(i, j, k) = 0.0;
                });
                // } else if (use_flattening == 1) {

                //     uflatten(obx, q_arr, flatn_arr, QPRES);

                // } else {
                //     amrex::ParallelFor(obx, [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept {
                //         flatn_arr(i, j, k) = 1.0;
                //     });
            }

            const Box& xbx = amrex::surroundingNodes(bx, 0);
            const Box& gxbx = amrex::grow(xbx, 1);
#if AMREX_SPACEDIM >= 2
            const Box& ybx = amrex::surroundingNodes(bx, 1);
            const Box& gybx = amrex::grow(ybx, 1);
#endif
#if AMREX_SPACEDIM == 3
            const Box& zbx = amrex::surroundingNodes(bx, 2);
            const Box& gzbx = amrex::grow(zbx, 1);
#endif

            const Box& qbx = amrex::grow(bx, NUM_GROW);

            src_q.resize(qbx, NQSRC);
            Elixir elix_src_q = src_q.elixir();
            fab_size += src_q.nBytes();
            Array4<Real> const src_q_arr = src_q.array();

            Array4<Real> const src_arr = sources_for_hydro.array(mfi);

            // src_to_prim(qbx, q_arr, src_arr, src_q_arr);

            // work on the interface states

            Uxm.resize(gxbx, NQ);
            Elixir elix_Uxm = Uxm.elixir();
            fab_size += Uxm.nBytes();

            Uxp.resize(obx, NQ);
            Elixir elix_Uxp = Uxp.elixir();
            fab_size += Uxp.nBytes();

            Array4<Real> const Uxm_arr = Uxm.array();
            Array4<Real> const Uxp_arr = Uxp.array();

            // #if AMREX_SPACEDIM >= 2
            //             Uym.resize(obx, NQ);
            //             Elixir elix_Uym = Uym.elixir();
            //             fab_size += Uym.nBytes();

            //             Uyp.resize(obx, NQ);
            //             Elixir elix_Uyp = Uyp.elixir();
            //             fab_size += Uyp.nBytes();

            //             Array4<Real> const Uym_arr = Uym.array();
            //             Array4<Real> const Uyp_arr = Uyp.array();

            // #endif

            // #if AMREX_SPACEDIM == 3
            //             Uzm.resize(obx, NQ);
            //             Elixir elix_Uzm = Uzm.elixir();
            //             fab_size += Uzm.nBytes();

            //             Uzp.resize(obx, NQ);
            //             Elixir elix_Uzp = Uzp.elixir();
            //             fab_size += Uzp.nBytes();

            //             Array4<Real> const Uzm_arr = Uzm.array();
            //             Array4<Real> const Uzp_arr = Uzp.array();

            // #endif

            // if (ppm_type == 0) {

            dq.resize(obx, NQ);
            Elixir elix_dq = dq.elixir();
            fab_size += dq.nBytes();
            auto dq_arr = dq.array();

            plm(obx, bx, q_arr, U_arr, flatn_arr, qaux_arr, src_q_arr, dq_arr, Uxm_arr, Uxp_arr,
#if AMREX_SPACEDIM >= 2
                Uym_arr, Uyp_arr,
#endif
#if AMREX_SPACEDIM == 3
                Uzm_arr, Uzp_arr,
#endif
#if (AMREX_SPACEDIM < 3)
                dLogArea_arr,
#endif
                dt, dx);

            //             } else {

            //                 ctu_ppm_states(obx, bx, q_arr, flatn_arr, qaux_arr, src_q_arr, Uxm_arr, Uxp_arr,
            // #if AMREX_SPACEDIM >= 2
            //                                Uym_arr, Uyp_arr,
            // #endif
            // #if AMREX_SPACEDIM == 3
            //                                Uzm_arr, Uzp_arr,
            // #endif
            // #if AMREX_SPACEDIM < 3
            //                                dLogArea_arr,
            // #endif
            //                                dt);
            // }

            // div.resize(obx, 1);
            // Elixir elix_div = div.elixir();
            // fab_size += div.nBytes();
            // auto div_arr = div.array();

            // // compute divu -- we'll use this later when doing the artifical viscosity
            // divu(obx, q_arr, div_arr);

            U_int.resize(obx, NQ);
            Elixir elix_U_int = U_int.elixir();
            fab_size += U_int.nBytes();
            Array4<Real> const U_int_arr = U_int.array();

            flux[0].resize(gxbx, NUM_STATE);
            Elixir elix_flux_x = flux[0].elixir();
            fab_size += flux[0].nBytes();
            Array4<Real> const flux0_arr = (flux[0]).array();

            qe[0].resize(gxbx, NGDNV);
            Elixir elix_qe_x = qe[0].elixir();
            auto qex_arr = qe[0].array();
            fab_size += qe[0].nBytes();

            // #if AMREX_SPACEDIM >= 2
            //             flux[1].resize(gybx, NUM_STATE);
            //             Elixir elix_flux_y = flux[1].elixir();
            //             fab_size += flux[1].nBytes();
            //             Array4<Real> const flux1_arr = (flux[1]).array();

            //             qe[1].resize(gybx, NGDNV);
            //             Elixir elix_qe_y = qe[1].elixir();
            //             auto qey_arr = qe[1].array();
            //             fab_size += qe[1].nBytes();
            // #endif

            // #if AMREX_SPACEDIM == 3
            //             flux[2].resize(gzbx, NUM_STATE);
            //             Elixir elix_flux_z = flux[2].elixir();
            //             fab_size += flux[2].nBytes();
            //             Array4<Real> const flux2_arr = (flux[2]).array();

            //             qe[2].resize(gzbx, NGDNV);
            //             Elixir elix_qe_z = qe[2].elixir();
            //             auto qez_arr = qe[2].array();
            //             fab_size += qe[2].nBytes();
            // #endif

#if AMREX_SPACEDIM <= 2
            if (!Geom().IsCartesian()) {
                pradial.resize(xbx, 1);
            }
            Elixir elix_pradial = pradial.elixir();
            fab_size += pradial.nBytes();
#endif

#if AMREX_SPACEDIM == 1
            // cmpflx_plus_godunov(xbx, Uxm_arr, Uxp_arr, flux0_arr, U_int_arr, qex_arr, qaux_arr,
            //                     shk_arr, 0);
            hlle(xbx, Uxm_arr, Uxp_arr, flux0_arr, U_int_arr, qex_arr, 0, dt);

            AllPrint() << "flux0_arr = " << flux0_arr(65, 0, 0, UEDEN) << std::endl;

#endif  // 1-d

            // #if AMREX_SPACEDIM >= 2
            //                 ftmp1.resize(obx, NUM_STATE);
            //             Elixir elix_ftmp1 = ftmp1.elixir();
            //             auto ftmp1_arr = ftmp1.array();
            //             fab_size += ftmp1.nBytes();

            //             ftmp2.resize(obx, NUM_STATE);
            //             Elixir elix_ftmp2 = ftmp2.elixir();
            //             auto ftmp2_arr = ftmp2.array();
            //             fab_size += ftmp2.nBytes();

            //             qgdnvtmp1.resize(obx, NGDNV);
            //             Elixir elix_qgdnvtmp1 = qgdnvtmp1.elixir();
            //             auto qgdnvtmp1_arr = qgdnvtmp1.array();
            //             fab_size += qgdnvtmp1.nBytes();

            // #if AMREX_SPACEDIM == 3
            //             qgdnvtmp2.resize(obx, NGDNV);
            //             Elixir elix_qgdnvtmp2 = qgdnvtmp2.elixir();
            //             auto qgdnvtmp2_arr = qgdnvtmp2.array();
            //             fab_size += qgdnvtmp2.nBytes();
            // #endif

            //             ql.resize(obx, NQ);
            //             Elixir elix_ql = ql.elixir();
            //             auto ql_arr = ql.array();
            //             fab_size += ql.nBytes();

            //             qr.resize(obx, NQ);
            //             Elixir elix_qr = qr.elixir();
            //             auto qr_arr = qr.array();
            //             fab_size += qr.nBytes();
            // #endif

            // #if AMREX_SPACEDIM == 2

            //             const amrex::Real hdt = 0.5 * dt;
            //             const amrex::Real hdtdx = 0.5 * dt / dx[0];
            //             const amrex::Real hdtdy = 0.5 * dt / dx[1];

            //             // compute F^x
            //             // [lo(1), lo(2)-1, 0], [hi(1)+1, hi(2)+1, 0]
            //             const Box& cxbx = amrex::grow(xbx, IntVect(AMREX_D_DECL(0, 1, 0)));

            //             // ftmp1 = fx
            //             // rftmp1 = rfx
            //             // qgdnvtmp1 = qgdnxv
            //             cmpflx_plus_godunov(cxbx, Uxm_arr, Uxp_arr, ftmp1_arr, U_int_arr, qgdnvtmp1_arr,
            //                                 qaux_arr, shk_arr, 0);

            //             // compute F^y
            //             // [lo(1)-1, lo(2), 0], [hi(1)+1, hi(2)+1, 0]
            //             const Box& cybx = amrex::grow(ybx, IntVect(AMREX_D_DECL(1, 0, 0)));

            //             // ftmp2 = fy
            //             // rftmp2 = rfy
            //             cmpflx_plus_godunov(cybx, Uym_arr, Uyp_arr, ftmp2_arr, U_int_arr, qey_arr, qaux_arr,
            //                                 shk_arr, 1);

            //             // add the transverse flux difference in y to the x states
            //             // [lo(1), lo(2), 0], [hi(1)+1, hi(2), 0]

            //             // ftmp2 = fy
            //             // rftmp2 = rfy
            //             trans_single(xbx, 1, 0, Uxm_arr, ql_arr, Uxp_arr, qr_arr, qaux_arr, ftmp2_arr, qey_arr,
            //                          areay_arr, vol_arr, hdt, hdtdy);

            //             reset_edge_state_thermo(xbx, ql.array());

            //             reset_edge_state_thermo(xbx, qr.array());

            //             // solve the final Riemann problem axross the x-interfaces

            //             cmpflx_plus_godunov(xbx, ql_arr, qr_arr, flux0_arr, U_int_arr, qex_arr, qaux_arr,
            //                                 shk_arr, 0);

            //             // add the transverse flux difference in x to the y states
            //             // [lo(1), lo(2), 0], [hi(1), hi(2)+1, 0]

            //             // ftmp1 = fx
            //             // rftmp1 = rfx
            //             // qgdnvtmp1 = qgdnvx

            //             trans_single(ybx, 0, 1, Uym_arr, ql_arr, Uyp_arr, qr_arr, qaux_arr, ftmp1_arr,
            //                          qgdnvtmp1_arr, areax_arr, vol_arr, hdt, hdtdx);

            //             reset_edge_state_thermo(ybx, ql.array());

            //             reset_edge_state_thermo(ybx, qr.array());

            //             // solve the final Riemann problem axross the y-interfaces

            //             cmpflx_plus_godunov(ybx, ql_arr, qr_arr, flux1_arr, U_int_arr, qey_arr, qaux_arr,
            //                                 shk_arr, 1);
            // #endif  // 2-d

            // #if AMREX_SPACEDIM == 3

            //             const amrex::Real hdt = 0.5 * dt;

            //             const amrex::Real hdtdx = 0.5 * dt / dx[0];
            //             const amrex::Real hdtdy = 0.5 * dt / dx[1];
            //             const amrex::Real hdtdz = 0.5 * dt / dx[2];

            //             const amrex::Real cdtdx = dt / dx[0] / 3.0;
            //             const amrex::Real cdtdy = dt / dx[1] / 3.0;
            //             const amrex::Real cdtdz = dt / dx[2] / 3.0;

            //             // compute F^x
            //             // [lo(1), lo(2)-1, lo(3)-1], [hi(1)+1, hi(2)+1, hi(3)+1]
            //             const Box& cxbx = amrex::grow(xbx, IntVect(AMREX_D_DECL(0, 1, 1)));

            //             // ftmp1 = fx
            //             // rftmp1 = rfx
            //             // qgdnvtmp1 = qgdnxv
            //             cmpflx_plus_godunov(cxbx, Uxm_arr, Uxp_arr, ftmp1_arr, U_int_arr, qgdnvtmp1_arr,
            //                                 qaux_arr, shk_arr, 0);

            //             // [lo(1), lo(2), lo(3)-1], [hi(1), hi(2)+1, hi(3)+1]
            //             const Box& tyxbx = amrex::grow(ybx, IntVect(AMREX_D_DECL(0, 0, 1)));

            //             qmyx.resize(tyxbx, NQ);
            //             Elixir elix_qmyx = qmyx.elixir();
            //             auto qmyx_arr = qmyx.array();
            //             fab_size += qmyx.nBytes();

            //             qpyx.resize(tyxbx, NQ);
            //             Elixir elix_qpyx = qpyx.elixir();
            //             auto qpyx_arr = qpyx.array();
            //             fab_size += qpyx.nBytes();

            //             // ftmp1 = fx
            //             // rftmp1 = rfx
            //             // qgdnvtmp1 = qgdnvx
            //             trans_single(tyxbx, 0, 1, Uym_arr, qmyx_arr, Uyp_arr, qpyx_arr, qaux_arr, ftmp1_arr,
            //                          qgdnvtmp1_arr, hdt, cdtdx);

            //             reset_edge_state_thermo(tyxbx, qmyx.array());

            //             reset_edge_state_thermo(tyxbx, qpyx.array());

            //             // [lo(1), lo(2)-1, lo(3)], [hi(1), hi(2)+1, hi(3)+1]
            //             const Box& tzxbx = amrex::grow(zbx, IntVect(AMREX_D_DECL(0, 1, 0)));

            //             qmzx.resize(tzxbx, NQ);
            //             Elixir elix_qmzx = qmzx.elixir();
            //             auto qmzx_arr = qmzx.array();
            //             fab_size += qmzx.nBytes();

            //             qpzx.resize(tzxbx, NQ);
            //             Elixir elix_qpzx = qpzx.elixir();
            //             auto qpzx_arr = qpzx.array();
            //             fab_size += qpzx.nBytes();

            //             trans_single(tzxbx, 0, 2, Uzm_arr, qmzx_arr, Uzp_arr, qpzx_arr, qaux_arr, ftmp1_arr,
            //                          qgdnvtmp1_arr, hdt, cdtdx);

            //             reset_edge_state_thermo(tzxbx, qmzx.array());

            //             reset_edge_state_thermo(tzxbx, qpzx.array());

            //             // compute F^y
            //             // [lo(1)-1, lo(2), lo(3)-1], [hi(1)+1, hi(2)+1, hi(3)+1]
            //             const Box& cybx = amrex::grow(ybx, IntVect(AMREX_D_DECL(1, 0, 1)));

            //             // ftmp1 = fy
            //             // rftmp1 = rfy
            //             // qgdnvtmp1 = qgdnvy
            //             cmpflx_plus_godunov(cybx, Uym_arr, Uyp_arr, ftmp1_arr, U_int_arr, qgdnvtmp1_arr,
            //                                 qaux_arr, shk_arr, 1);

            //             // [lo(1), lo(2), lo(3)-1], [hi(1)+1, hi(2), lo(3)+1]
            //             const Box& txybx = amrex::grow(xbx, IntVect(AMREX_D_DECL(0, 0, 1)));

            //             qmxy.resize(txybx, NQ);
            //             Elixir elix_qmxy = qmxy.elixir();
            //             auto qmxy_arr = qmxy.array();
            //             fab_size += qmxy.nBytes();

            //             qpxy.resize(txybx, NQ);
            //             Elixir elix_qpxy = qpxy.elixir();
            //             auto qpxy_arr = qpxy.array();
            //             fab_size += qpxy.nBytes();

            //             // ftmp1 = fy
            //             // rftmp1 = rfy
            //             // qgdnvtmp1 = qgdnvy
            //             trans_single(txybx, 1, 0, Uxm_arr, qmxy_arr, Uxp_arr, qpxy_arr, qaux_arr, ftmp1_arr,
            //                          qgdnvtmp1_arr, hdt, cdtdy);

            //             reset_edge_state_thermo(txybx, qmxy.array());

            //             reset_edge_state_thermo(txybx, qpxy.array());

            //             // [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2), lo(3)+1]
            //             const Box& tzybx = amrex::grow(zbx, IntVect(AMREX_D_DECL(1, 0, 0)));

            //             qmzy.resize(tzybx, NQ);
            //             Elixir elix_qmzy = qmzy.elixir();
            //             auto qmzy_arr = qmzy.array();
            //             fab_size += qmzy.nBytes();

            //             qpzy.resize(tzybx, NQ);
            //             Elixir elix_qpzy = qpzy.elixir();
            //             auto qpzy_arr = qpzy.array();
            //             fab_size += qpzy.nBytes();

            //             // ftmp1 = fy
            //             // rftmp1 = rfy
            //             // qgdnvtmp1 = qgdnvy
            //             trans_single(tzybx, 1, 2, Uzm_arr, qmzy_arr, Uzp_arr, qpzy_arr, qaux_arr, ftmp1_arr,
            //                          qgdnvtmp1_arr, hdt, cdtdy);

            //             reset_edge_state_thermo(tzybx, qmzy.array());

            //             reset_edge_state_thermo(tzybx, qpzy.array());

            //             // compute F^z
            //             // [lo(1)-1, lo(2)-1, lo(3)], [hi(1)+1, hi(2)+1, hi(3)+1]
            //             const Box& czbx = amrex::grow(zbx, IntVect(AMREX_D_DECL(1, 1, 0)));

            //             // ftmp1 = fz
            //             // rftmp1 = rfz
            //             // qgdnvtmp1 = qgdnvz
            //             cmpflx_plus_godunov(czbx, Uzm_arr, Uzp_arr, ftmp1_arr, U_int_arr, qgdnvtmp1_arr,
            //                                 qaux_arr, shk_arr, 2);

            //             // [lo(1)-1, lo(2)-1, lo(3)], [hi(1)+1, hi(2)+1, lo(3)]
            //             const Box& txzbx = amrex::grow(xbx, IntVect(AMREX_D_DECL(0, 1, 0)));

            //             qmxz.resize(txzbx, NQ);
            //             Elixir elix_qmxz = qmxz.elixir();
            //             auto qmxz_arr = qmxz.array();
            //             fab_size += qmxz.nBytes();

            //             qpxz.resize(txzbx, NQ);
            //             Elixir elix_qpxz = qpxz.elixir();
            //             auto qpxz_arr = qpxz.array();
            //             fab_size += qpxz.nBytes();

            //             // ftmp1 = fz
            //             // rftmp1 = rfz
            //             // qgdnvtmp1 = qgdnvz
            //             trans_single(txzbx, 2, 0, Uxm_arr, qmxz_arr, Uxp_arr, qpxz_arr, qaux_arr, ftmp1_arr,
            //                          qgdnvtmp1_arr, hdt, cdtdz);

            //             reset_edge_state_thermo(txzbx, qmxz.array());

            //             reset_edge_state_thermo(txzbx, qpxz.array());

            //             // [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2)+1, lo(3)]
            //             const Box& tyzbx = amrex::grow(ybx, IntVect(AMREX_D_DECL(1, 0, 0)));

            //             qmyz.resize(tyzbx, NQ);
            //             Elixir elix_qmyz = qmyz.elixir();
            //             auto qmyz_arr = qmyz.array();
            //             fab_size += qmyz.nBytes();

            //             qpyz.resize(tyzbx, NQ);
            //             Elixir elix_qpyz = qpyz.elixir();
            //             auto qpyz_arr = qpyz.array();
            //             fab_size += qpyz.nBytes();

            //             // ftmp1 = fz
            //             // rftmp1 = rfz
            //             // qgdnvtmp1 = qgdnvz
            //             trans_single(tyzbx, 2, 1, Uym_arr, qmyz_arr, Uyp_arr, qpyz_arr, qaux_arr, ftmp1_arr,
            //                          qgdnvtmp1_arr, hdt, cdtdz);

            //             reset_edge_state_thermo(tyzbx, qmyz.array());

            //             reset_edge_state_thermo(tyzbx, qpyz.array());

            //             // we now have q?zx, q?yx, q?zy, q?xy, q?yz, q?xz

            //             //
            //             // Use qx?, q?yz, q?zy to compute final x-flux
            //             //

            //             // compute F^{y|z}
            //             // [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2)+1, hi(3)]
            //             const Box& cyzbx = amrex::grow(ybx, IntVect(AMREX_D_DECL(1, 0, 0)));

            //             // ftmp1 = fyz
            //             // rftmp1 = rfyz
            //             // qgdnvtmp1 = qgdnvyz
            //             cmpflx_plus_godunov(cyzbx, qmyz_arr, qpyz_arr, ftmp1_arr, U_int_arr, qgdnvtmp1_arr,
            //                                 qaux_arr, shk_arr, 1);

            //             // compute F^{z|y}
            //             // [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2), hi(3)+1]
            //             const Box& czybx = amrex::grow(zbx, IntVect(AMREX_D_DECL(1, 0, 0)));

            //             // ftmp2 = fzy
            //             // rftmp2 = rfzy
            //             // qgdnvtmp2 = qgdnvzy
            //             cmpflx_plus_godunov(czybx, qmzy_arr, qpzy_arr, ftmp2_arr, U_int_arr, qgdnvtmp2_arr,
            //                                 qaux_arr, shk_arr, 2);

            //             // compute the corrected x interface states and fluxes
            //             // [lo(1), lo(2), lo(3)], [hi(1)+1, hi(2), hi(3)]

            //             trans_final(xbx, 0, 1, 2, Uxm_arr, ql_arr, Uxp_arr, qr_arr, qaux_arr, ftmp1_arr,
            //                         ftmp2_arr, qgdnvtmp1_arr, qgdnvtmp2_arr, hdt, hdtdx, hdtdy, hdtdz);

            //             reset_edge_state_thermo(xbx, ql.array());

            //             reset_edge_state_thermo(xbx, qr.array());

            //             cmpflx_plus_godunov(xbx, ql_arr, qr_arr, flux0_arr, U_int_arr, qex_arr, qaux_arr,
            //                                 shk_arr, 0);

            //             //
            //             // Use qy?, q?zx, q?xz to compute final y-flux
            //             //

            //             // compute F^{z|x}
            //             // [lo(1), lo(2)-1, lo(3)], [hi(1), hi(2)+1, hi(3)+1]
            //             const Box& czxbx = amrex::grow(zbx, IntVect(AMREX_D_DECL(0, 1, 0)));

            //             // ftmp1 = fzx
            //             // rftmp1 = rfzx
            //             // qgdnvtmp1 = qgdnvzx
            //             cmpflx_plus_godunov(czxbx, qmzx_arr, qpzx_arr, ftmp1_arr, U_int_arr, qgdnvtmp1_arr,
            //                                 qaux_arr, shk_arr, 2);

            //             // compute F^{x|z}
            //             // [lo(1), lo(2)-1, lo(3)], [hi(1)+1, hi(2)+1, hi(3)]
            //             const Box& cxzbx = amrex::grow(xbx, IntVect(AMREX_D_DECL(0, 1, 0)));

            //             // ftmp2 = fxz
            //             // rftmp2 = rfxz
            //             // qgdnvtmp2 = qgdnvxz
            //             cmpflx_plus_godunov(cxzbx, qmxz_arr, qpxz_arr, ftmp2_arr, U_int_arr, qgdnvtmp2_arr,
            //                                 qaux_arr, shk_arr, 0);

            //             // Compute the corrected y interface states and fluxes
            //             // [lo(1), lo(2), lo(3)], [hi(1), hi(2)+1, hi(3)]

            //             trans_final(ybx, 1, 0, 2, Uym_arr, ql_arr, Uyp_arr, qr_arr, qaux_arr, ftmp2_arr,
            //                         ftmp1_arr, qgdnvtmp2_arr, qgdnvtmp1_arr, hdt, hdtdx, hdtdy, hdtdz);

            //             reset_edge_state_thermo(ybx, ql.array());

            //             reset_edge_state_thermo(ybx, qr.array());

            //             // Compute the final F^y
            //             // [lo(1), lo(2), lo(3)], [hi(1), hi(2)+1, hi(3)]
            //             cmpflx_plus_godunov(ybx, ql_arr, qr_arr, flux1_arr, U_int_arr, qey_arr, qaux_arr,
            //                                 shk_arr, 1);

            //             //
            //             // Use qz?, q?xy, q?yx to compute final z-flux
            //             //

            //             // compute F^{x|y}
            //             // [lo(1), lo(2), lo(3)-1], [hi(1)+1, hi(2), hi(3)+1]
            //             const Box& cxybx = amrex::grow(xbx, IntVect(AMREX_D_DECL(0, 0, 1)));

            //             // ftmp1 = fxy
            //             // rftmp1 = rfxy
            //             // qgdnvtmp1 = qgdnvxy
            //             cmpflx_plus_godunov(cxybx, qmxy_arr, qpxy_arr, ftmp1_arr, U_int_arr, qgdnvtmp1_arr,
            //                                 qaux_arr, shk_arr, 0);

            //             // compute F^{y|x}
            //             // [lo(1), lo(2), lo(3)-1], [hi(1), hi(2)+dg(2), hi(3)+1]
            //             const Box& cyxbx = amrex::grow(ybx, IntVect(AMREX_D_DECL(0, 0, 1)));

            //             // ftmp2 = fyx
            //             // rftmp2 = rfyx
            //             // qgdnvtmp2 = qgdnvyx
            //             cmpflx_plus_godunov(cyxbx, qmyx_arr, qpyx_arr, ftmp2_arr, U_int_arr, qgdnvtmp2_arr,
            //                                 qaux_arr, shk_arr, 1);

            //             // compute the corrected z interface states and fluxes
            //             // [lo(1), lo(2), lo(3)], [hi(1), hi(2), hi(3)+1]

            //             trans_final(zbx, 2, 0, 1, Uzm_arr, ql_arr, Uzp_arr, qr_arr, qaux_arr, ftmp1_arr,
            //                         ftmp2_arr, qgdnvtmp1_arr, qgdnvtmp2_arr, hdt, hdtdx, hdtdy, hdtdz);

            //             reset_edge_state_thermo(zbx, ql.array());

            //             reset_edge_state_thermo(zbx, qr.array());

            //             // compute the final z fluxes F^z
            //             // [lo(1), lo(2), lo(3)], [hi(1), hi(2), hi(3)+1]

            //             cmpflx_plus_godunov(zbx, ql_arr, qr_arr, flux2_arr, U_int_arr, qez_arr, qaux_arr,
            //                                 shk_arr, 2);

            // #endif  // 3-d

            // clean the fluxes

            // for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

            //     const Box& nbx = amrex::surroundingNodes(bx, idir);

            //     Array4<Real> const flux_arr = (flux[idir]).array();
            //     Array4<Real const> const uin_arr = Sborder.array(mfi);

            //     // Zero out shock and temp fluxes -- these are physically meaningless here
            //     amrex::ParallelFor(nbx, [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept {
            //         flux_arr(i, j, k, UTEMP) = 0.e0_rt;
            //     });

            //     // apply_av(nbx, idir, div_arr, uin_arr, flux_arr);

            //     // if (limit_fluxes_on_small_dens == 1) {
            //     //     limit_hydro_fluxes_on_small_dens(nbx, idir, Sborder.array(mfi), q.array(mfi),
            //     //                                      volume.array(mfi), flux[idir].array(),
            //     //                                      area[idir].array(mfi), dt);
            //     // }

            //     // if (limit_fluxes_on_large_vel == 1) {
            //     //     limit_hydro_fluxes_on_large_vel(nbx, idir, Sborder.array(mfi), q.array(mfi),
            //     //                                     volume.array(mfi), flux[idir].array(),
            //     //                                     area[idir].array(mfi), dt);
            //     // }

            //     // normalize_species_fluxes(nbx, flux_arr);
            // }

            // conservative update
            Array4<Real> const update_arr = hydro_source.array(mfi);

            // Array4<Real> const flx_arr = (flux[0]).array();
            // Array4<Real> const qx_arr = (qe[0]).array();

            // #if AMREX_SPACEDIM >= 2
            //             Array4<Real> const fly_arr = (flux[1]).array();
            //             Array4<Real> const qy_arr = (qe[1]).array();
            // #endif

            // #if AMREX_SPACEDIM == 3
            //             Array4<Real> const flz_arr = (flux[2]).array();
            //             Array4<Real> const qz_arr = (qe[2]).array();
            // #endif

            consup_rhd(bx, update_arr, flux0_arr, areax_arr,
                       // #if AMREX_SPACEDIM >= 2
                       //                        fly_arr, areay_arr,
                       // #endif
                       // #if AMREX_SPACEDIM == 3
                       //                        flz_arr, areaz_arr,
                       // #endif
                       vol_arr);

            for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

                const Box& nbx = amrex::surroundingNodes(bx, idir);

                Array4<Real> const flux_arr = (flux[idir]).array();
                Array4<Real const> const area_arr = (area[idir]).array(mfi);

                //                 scale_flux(nbx,
                // #if AMREX_SPACEDIM == 1
                //                            qex_arr,
                // #endif
                //                            flux_arr, area_arr, dt);

                if (idir == 0) {
#if AMREX_SPACEDIM <= 2
                    Array4<Real> pradial_fab = pradial.array();
#endif

                    // get the scaled radial pressure -- we need to treat this specially
#if AMREX_SPACEDIM <= 2

#if AMREX_SPACEDIM == 1
                    if (!Geom().IsCartesian()) {
#elif AMREX_SPACEDIM == 2
                    if (!mom_flux_has_p(0, 0, coord)) {
#endif
                        // amrex::ParallelFor(nbx,
                        //                    [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept {
                        //                        pradial_fab(i, j, k) = qex_arr(i, j, k, GDPRES) * dt;
                        //                    });
                    }

#endif
                }

                // Store the fluxes from this advance. For simplified SDC integration we
                // only need to do this on the last iteration.

                bool add_fluxes = true;

                // if (time_integration_method == SimplifiedSpectralDeferredCorrections &&
                //     sdc_iteration != sdc_iters - 1) {
                //     add_fluxes = false;
                // }

                if (add_fluxes) {

                    Array4<Real> const flux_fab = (flux[idir]).array();
                    Array4<Real> fluxes_fab = (*fluxes[idir]).array(mfi);
                    const int numcomp = NUM_STATE;

                    AMREX_HOST_DEVICE_FOR_4D(mfi.nodaltilebox(idir), numcomp, i, j, k, n,
                                             { fluxes_fab(i, j, k, n) += flux_fab(i, j, k, n); });

#if AMREX_SPACEDIM <= 2

#if AMREX_SPACEDIM == 1
                    if (idir == 0 && !Geom().IsCartesian()) {
#elif AMREX_SPACEDIM == 2
                    if (idir == 0 && !mom_flux_has_p(0, 0, coord)) {
#endif
                        Array4<Real> pradial_fab = pradial.array();
                        Array4<Real> P_radial_fab = P_radial.array(mfi);

                        AMREX_HOST_DEVICE_FOR_4D(mfi.nodaltilebox(0), 1, i, j, k, n, {
                            P_radial_fab(i, j, k, 0) += pradial_fab(i, j, k, 0);
                        });
                    }

#endif

                }  // add_fluxes

                Array4<Real> const flux_fab = (flux[idir]).array();
                Array4<Real> mass_fluxes_fab = (*mass_fluxes[idir]).array(mfi);

                AMREX_HOST_DEVICE_FOR_4D(mfi.nodaltilebox(idir), 1, i, j, k, n, {
                    // This is a copy, not an add, since we need mass_fluxes to be
                    // only this subcycle's data when we evaluate the gravitational
                    // forces.

                    mass_fluxes_fab(i, j, k, 0) = flux_fab(i, j, k, URHO);
                });

            }  // idir loop

            if (track_grid_losses == 1) {

#pragma gpu box(bx)
                ca_track_grid_losses(
                    AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                    BL_TO_FORTRAN_ANYD(flux[0]),
#if AMREX_SPACEDIM >= 2
                    BL_TO_FORTRAN_ANYD(flux[1]),
#endif
#if AMREX_SPACEDIM == 3
                    BL_TO_FORTRAN_ANYD(flux[2]),
#endif
                    AMREX_MFITER_REDUCE_SUM(&mass_lost), AMREX_MFITER_REDUCE_SUM(&xmom_lost),
                    AMREX_MFITER_REDUCE_SUM(&ymom_lost), AMREX_MFITER_REDUCE_SUM(&zmom_lost),
                    AMREX_MFITER_REDUCE_SUM(&eden_lost), AMREX_MFITER_REDUCE_SUM(&xang_lost),
                    AMREX_MFITER_REDUCE_SUM(&yang_lost), AMREX_MFITER_REDUCE_SUM(&zang_lost));
            }

#ifdef AMREX_USE_GPU
            // Check if we're going to run out of memory in the next MFIter iteration.
            // If so, do a synchronize here so that we don't oversubscribe GPU memory.
            // Note that this will capture the case where we started with more memory
            // than what the GPU has, on the logic that even in that case, it makes
            // sense to not further pile on the oversubscription demands.

            // This could (and should) be generalized in the future to operate with
            // more granularity than the MFIter loop boundary. We would have potential
            // synchronization points prior to each of the above kernel launches, and
            // we would check whether the sum of all previously allocated fabs would
            // result in oversubscription, including any contributions from a partial
            // MFIter loop. A further optimization would be to not apply a device
            // synchronize, but rather to use CUDA events to poll on a check about
            // whether enough memory has freed up to begin the next iteration, and then
            // immediately proceed to the next kernel when there's enough space for it.

            current_size += fab_size;
            if (current_size + fab_size >= Gpu::Device::totalGlobalMem()) {
                Gpu::Device::synchronize();
                current_size = starting_size;
            }
#endif

            if (oversubscribed) {
                // q[mfi].prefetchToHost();
                qaux[mfi].prefetchToHost();
                volume[mfi].prefetchToHost();
                Sborder[mfi].prefetchToHost();
                hydro_source[mfi].prefetchToHost();
                for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                    area[i][mfi].prefetchToHost();
                    (*fluxes[i])[mfi].prefetchToHost();
                }
#if AMREX_SPACEDIM < 3
                dLogArea[0][mfi].prefetchToHost();
                P_radial[mfi].prefetchToHost();
#endif
            }

        }  // MFIter loop
    }      // OMP loop

    // Flush Fortran output

    if (verbose) flush_output();

    if (track_grid_losses) {
        material_lost_through_boundary_temp[0] += mass_lost;
        material_lost_through_boundary_temp[1] += xmom_lost;
        material_lost_through_boundary_temp[2] += ymom_lost;
        material_lost_through_boundary_temp[3] += zmom_lost;
        material_lost_through_boundary_temp[4] += eden_lost;
        material_lost_through_boundary_temp[5] += xang_lost;
        material_lost_through_boundary_temp[6] += yang_lost;
        material_lost_through_boundary_temp[7] += zang_lost;
    }

    if (print_update_diagnostics) {

        bool local = true;
        Vector<Real> hydro_update = evaluate_source_change(hydro_source, dt, local);

#ifdef BL_LAZY
        Lazy::QueueReduction([=]() mutable {
#endif
            ParallelDescriptor::ReduceRealSum(hydro_update.dataPtr(), hydro_update.size(),
                                              ParallelDescriptor::IOProcessorNumber());

            if (ParallelDescriptor::IOProcessor())
                std::cout << std::endl
                          << "  Contributions to the state from the hydro source:" << std::endl;

            print_source_change(hydro_update);

#ifdef BL_LAZY
        });
#endif
    }

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... Leaving construct_ctu_rhd_source()" << std::endl << std::endl;

    if (verbose > 0) {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
        Lazy::QueueReduction([=]() mutable {
#endif
            ParallelDescriptor::ReduceRealMax(run_time, IOProc);

            if (ParallelDescriptor::IOProcessor())
                std::cout << "Castro::construct_ctu_rhd_source() time = " << run_time << "\n"
                          << "\n";
#ifdef BL_LAZY
        });
#endif
    }
}

void Castro::plm(const Box& bx, const Box& vbx, Array4<Real const> const& q_arr,
                 Array4<Real const> const& U_arr, Array4<Real const> const& flatn_arr,
                 Array4<Real const> const& qaux_arr, Array4<Real const> const& srcQ,
                 Array4<Real> const& dU, Array4<Real> const& Uxl, Array4<Real> const& Uxr,
#if AMREX_SPACEDIM >= 2
                 Array4<Real> const& Uyl, Array4<Real> const& Uyr,
#endif
#if AMREX_SPACEDIM == 3
                 Array4<Real> const& Uzl, Array4<Real> const& Uzr,
#endif
#if AMREX_SPACEDIM < 3
                 Array4<Real const> const& dloga,
#endif
                 const Real dt, const GpuArray<Real, AMREX_SPACEDIM>& dx) {

    // Compute the normal interface states by reconstructing
    // the primitive variables using piecewise linear slopes and doing
    // characteristic tracing.  We do not apply the transverse terms here.
    //
    // .. todo::
    //    we can get rid of the the different temporary q Godunov
    //    state arrays
    //

    // Compute all slopes
    amrex::ParallelFor(bx, [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept {
        for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {

            // calculate sound speed
            eos_t eos_state;
            eos_state.rho = q_arr(i, j, k, QRHO);
            eos_state.p = q_arr(i, j, k, QPRES);
            for (int n = 0; n < NumSpec; ++n) {
                eos_state.xn[n] = q_arr(i, j, k, QFS + n);
            }

            eos(eos_input_rp, eos_state);

            const Real cs = eos_state.cs;
            const Real h = eos_state.h;

            const Real udir = q_arr(i, j, k, QU + idir);
            const Real u = q_arr(i, j, k, QU);
            const Real v = 0.0_rt;  //q_arr(i, j, k, QV);
            const Real w = 0.0_rt;  //q_arr(i, j, k, QW);
            const Real v2 = u * u + v * v + w * w;
            const Real W = 1.0_rt / std::sqrt(1.0_rt - v2);

            Real eval[5];

            eval[0] = (udir * (1.0_rt - cs * cs) -
                       cs * std::sqrt((1.0_rt - v2) *
                                      (1.0_rt - v2 * cs * cs - udir * udir * (1.0_rt - cs * cs)))) /
                      (1 - v2 * cs * cs);

            eval[4] = (udir * (1.0_rt - cs * cs) +
                       cs * std::sqrt((1.0_rt - v2) *
                                      (1.0_rt - v2 * cs * cs - udir * udir * (1.0_rt - cs * cs)))) /
                      (1 - v2 * cs * cs);

            for (int n = 1; n < 4; ++n) {
                eval[n] = udir;
            }

            Real Am = (1.0_rt - udir * udir) / (1.0_rt - udir * eval[0]);
            Real Ap = (1.0_rt - udir * udir) / (1.0_rt - udir * eval[4]);

            Real rvec[5][5] = {
                {1.0_rt, h * W * Am * eval[0], h * W * v, h * W * w, h * W * Am - 1.0_rt},
                {1.0_rt / W, u, v, w, 1.0_rt - 1.0_rt / W},
                {W * v, 2.0_rt * h * W * W * u * v, h * (1.0_rt + 2.0_rt * W * W * v * v),
                 2.0_rt * h * W * W * v * w, 2.0_rt * h * W * W * v - W * v},
                {W * w, 2.0_rt * h * W * W * u * w, 2.0_rt * h * W * W * v * w,
                 h * (1.0_rt + 2.0_rt * W * W * w * w), 2.0_rt * h * W * W * w - W * w},
                {1.0_rt, h * W * Ap * eval[4], h * W * v, h * W * w, h * W * Ap - 1.0_rt}};

            // Print() << "rvec[0] = ";
            // for (int n = 0; n < 5; n++) {
            //     Print() << rvec[0][n] << ' ';
            // }
            // Print() << std::endl;

            Real lvec[5][5] = {
                {h * W * Ap * (udir - eval[4]) - udir -
                     W * W * (v2 - udir * udir) * (2.0_rt * h - 1.0_rt) * (udir - Ap * eval[4]) +
                     h * Ap * eval[4],
                 1.0_rt + W * W * (v2 - udir * udir) * (2.0_rt * h - 1.0_rt) * (1.0_rt - Ap) -
                     h * Ap,
                 W * W * v * (2.0_rt * h - 1.0_rt) * Ap * (udir - eval[4]),
                 W * W * w * (2.0_rt * h - 1.0_rt) * Ap * (udir - eval[4]),
                 -udir -
                     W * W * (v2 - udir * udir) * (2.0_rt * h - 1.0_rt) * (udir - Ap * eval[4]) +
                     h * Ap * eval[4]},
                {h - W, W * u, W * v, W * w, -W},
                {-v, u * v, 1.0_rt - u * u, 0.0_rt, -v},
                {-w, u * w, 0.0_rt, 1.0_rt - u * u, -w},
                {h * W * Am * (udir - eval[0]) - udir -
                     W * W * (v2 - udir * udir) * (2.0_rt * h - 1.0_rt) * (udir - Am * eval[0]) +
                     h * Am * eval[0],
                 1.0_rt + W * W * (v2 - udir * udir) * (2.0_rt * h - 1.0_rt) * (1.0_rt - Am) -
                     h * Am,
                 W * W * v * (2.0_rt * h - 1.0_rt) * Am * (udir - eval[0]),
                 W * W * w * (2.0_rt * h - 1.0_rt) * Am * (udir - eval[0]),
                 -udir -
                     W * W * (v2 - udir * udir) * (2.0_rt * h - 1.0_rt) * (udir - Am * eval[0]) +
                     h * Am * eval[0]}};

            Real Delta = h * h * h * W * (h - 1.0_rt) * (1.0_rt - udir * udir) *
                         (Ap * eval[4] - Am * eval[0]);

            for (int n = 0; n < 5; ++n) {
                lvec[0][n] *= h * h / Delta;
                lvec[1][n] *= W / (h - 1.0_rt);
                lvec[2][n] *= 0.0_rt;
                lvec[3][n] *= 0.0_rt;
                // lvec[2][n] /= h * (1.0_rt - udir * udir);
                // lvec[3][n] /= h * (1.0_rt - udir * udir);
                lvec[4][n] *= -h * h / Delta;
            }

            // define the reference states
            Real factor;

            for (int n = 0; n < NUM_STATE; n++) {

                if (n == QTEMP) {
                    continue;
                }

                Real deltal, deltar;
                if (idir == 0) {
                    deltal = U_arr(i, j, k, n) - U_arr(i - 1, j, k, n);
                    deltar = U_arr(i + 1, j, k, n) - U_arr(i, j, k, n);
#if AMREX_SPACEDIM >= 2
                } else if (idir == 1) {
                    deltal = U_arr(i, j, k, n) - U_arr(i, j - 1, k, n);
                    deltar = U_arr(i, j + 1, k, n) - U_arr(i, j, k, n);
#endif
#if AMREX_SPACEDIM == 3
                } else {
                    deltal = U_arr(i, j, k, n) - U_arr(i, j, k - 1, n);
                    deltar = U_arr(i, j, k + 1, n) - U_arr(i, j, k, n);
#endif
                }

                slope(dU(i, j, k, n), deltar, deltal, flatn_arr(i, j, k, n));
            }

            for (int n = 0; n < NUM_STATE; n++) {
                if (idir == 0) {
                    // this is one the right face of the current zone,
                    // so the fastest moving eigenvalue is e_val[4] = u + c
                    factor = 0.5_rt * (1.0_rt - dt / dx[0] * amrex::max(eval[4], 0.0_rt));

                    Uxl(i + 1, j, k, n) = U_arr(i, j, k, n) + factor * dU(i, j, k, n);

                    // left face of the current zone, so the fastest moving
                    // eigenvalue is e_val[0] = u - c
                    factor = 0.5_rt * (1.0_rt + dt / dx[0] * amrex::min(eval[0], 0.0_rt));

                    Uxr(i, j, k, n) = U_arr(i, j, k, n) - factor * dU(i, j, k, n);

#if AMREX_SPACEDIM >= 2
                } else if (idir == 1) {

                    factor = 0.5_rt * (1.0_rt - dt / dx[1] * amrex::max(eval[4], 0.0_rt));

                    Uyl(i, j + 1, k, n) = U_arr(i, j, k, n) + factor * dU(i, j, k, n);

                    factor = 0.5_rt * (1.0_rt + dt / dx[1] * amrex::min(eval[0], 0.0_rt));

                    Uyr(i, j, k, n) = U_arr(i, j, k, n) - factor * dU(i, j, k, n);
#endif
#if AMREX_SPACEDIM == 3
                } else {

                    factor = 0.5_rt * (1.0_rt - dt / dx[2] * amrex::max(eval[4], 0.0_rt));

                    Uzl(i, j, k + 1, n) = U_arr(i, j, k, n) + factor * dU(i, j, k, n);

                    factor = 0.5_rt * (1.0_rt + dt / dx[2] * amrex::min(eval[0], 0.0_rt));

                    Uzr(i, j, k, n) = U_arr(i, j, k, n) - factor * dU(i, j, k, n);
#endif
                }
            }

            // compute the Vhat functions
            Real betal[5];
            Real betar[5];

            for (int n = 0; n < 5; ++n) {
                Real summ = 0.0_rt;
                for (int m = 0; m < 5; ++m) {
                    summ += lvec[n][m] * dU(i, j, k, m);
                }

                betal[n] = 0.25_rt * dt / dx[idir] * (eval[4] - eval[n]) *
                           (std::copysign(1.0_rt, eval[n]) + 1.0_rt) * summ;

                betar[n] = 0.25_rt * dt / dx[idir] * (eval[0] - eval[n]) *
                           (1.0_rt - std::copysign(1.0_rt, eval[n])) * summ;
            }

            for (int n = 0; n < 5; n++) {

                if (n == QTEMP) {
                    continue;
                }

                Real sum_l = 0.0_rt;
                Real sum_r = 0.0_rt;
                for (int m = 0; m < 5; ++m) {
                    sum_l += betal[m] * rvec[m][n];
                    sum_r += betar[m] * rvec[m][n];
                }

                // compute the interface states
                if (idir == 0) {
                    Uxl(i + 1, j, k, n) += sum_l;
                    Uxr(i, j, k, n) += sum_r;
#if AMREX_SPACEDIM >= 2
                } else if (idir == 1) {
                    Uyl(i, j + 1, k, n) += sum_l;
                    Uyr(i, j, k, n) += sum_r;
#endif
#if AMREX_SPACEDIM == 3
                } else {
                    Uzl(i, j, k + 1, n) += sum_l;
                    Uzr(i, j, k, n) += sum_r;
#endif
                }

                // Print() << "qxr(" << i << "," << j << "," << k << ") = " << qxr(i,j,k,n) << ", q_arr(" << i << "," << j << "," << k << ") = " << q_arr(i,j,k,n) << std::endl;
            }
        }
    });

    // special care for reflecting BCs
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();

    const auto domlo = geom.Domain().loVect3d();
    const auto domhi = geom.Domain().hiVect3d();

    for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {

        bool lo_bc_test = lo_bc[idir] == Symmetry;
        bool hi_bc_test = hi_bc[idir] == Symmetry;

        // we have to do this after the loops above, since here we will
        // consider interfaces, not zones

        if (idir == 0) {
            if (lo_bc_test) {

                amrex::ParallelFor(bx, [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept {
                    // reset the left state at domlo(0) if needed -- it is outside the domain
                    if (i == domlo[0]) {
                        for (int n = 0; n < NUM_STATE; n++) {
                            if (n == UMX) {
                                Uxl(i, j, k, UMX) = -Uxr(i, j, k, UMX);
                            } else {
                                Uxl(i, j, k, n) = Uxr(i, j, k, n);
                            }
                        }
                    }
                });
            }

            if (hi_bc_test) {

                amrex::ParallelFor(bx, [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept {
                    // reset the right state at domhi(0)+1 if needed -- it is outside the domain
                    if (i == domhi[0] + 1) {
                        for (int n = 0; n < NUM_STATE; n++) {
                            if (n == UMX) {
                                Uxr(i, j, k, UMX) = -Uxl(i, j, k, UMX);
                            } else {
                                Uxr(i, j, k, n) = Uxl(i, j, k, n);
                            }
                        }
                    }
                });
            }

#if AMREX_SPACEDIM >= 2
        } else if (idir == 1) {
            if (lo_bc_test) {

                amrex::ParallelFor(bx, [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept {
                    // reset the left state at domlo(0) if needed -- it is outside the domain
                    if (j == domlo[1]) {
                        for (int n = 0; n < NUM_STATE; n++) {
                            if (n == UMY) {
                                Uyl(i, j, k, UMY) = -Uyr(i, j, k, UMY);
                            } else {
                                Uyl(i, j, k, n) = Uyr(i, j, k, n);
                            }
                        }
                    }
                });
            }

            if (hi_bc_test) {

                amrex::ParallelFor(bx, [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept {
                    // reset the right state at domhi(0)+1 if needed -- it is outside the domain
                    if (j == domhi[1] + 1) {
                        for (int n = 0; n < NUM_STATE; n++) {
                            if (n == UMY) {
                                Uyr(i, j, k, UMY) = -Uyl(i, j, k, UMY);
                            } else {
                                Uyr(i, j, k, n) = Uyl(i, j, k, n);
                            }
                        }
                    }
                });
            }

#endif
#if AMREX_SPACEDIM == 3
        } else {
            if (lo_bc_test) {

                amrex::ParallelFor(bx, [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept {
                    // reset the left state at domlo(0) if needed -- it is outside the domain
                    if (k == domlo[2]) {
                        for (int n = 0; n < NUM_STATE; n++) {
                            if (n == UMZ) {
                                Uzl(i, j, k, UMZ) = -Uzr(i, j, k, UMZ);
                            } else {
                                Uzl(i, j, k, n) = Uzr(i, j, k, n);
                            }
                        }
                    }
                });
            }

            if (hi_bc_test) {

                amrex::ParallelFor(bx, [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept {
                    // reset the right state at domhi(0)+1 if needed -- it is outside the domain
                    if (k == domhi[2] + 1) {
                        for (int n = 0; n < NUM_STATE; n++) {
                            if (n == UMZ) {
                                Uzr(i, j, k, UMZ) = -Uzl(i, j, k, UMZ);
                            } else {
                                Uzr(i, j, k, n) = Uzl(i, j, k, n);
                            }
                        }
                    }
                });
            }
#endif
        }
    }
}

void Castro::consup_rhd(const Box& bx, Array4<Real> const& update, Array4<Real> const& flux0,
                        Array4<Real const> const& area0,
#if AMREX_SPACEDIM >= 2
                        Array4<Real> const& flux1, Array4<Real const> const& area1,
#endif
#if AMREX_SPACEDIM == 3
                        Array4<Real> const& flux2, Array4<Real const> const& area2,
#endif
                        Array4<Real const> const& vol) {

    // For hydro, we will create an update source term that is
    // essentially the flux divergence.  This can be added with dt to
    // get the update

    amrex::ParallelFor(bx, NUM_STATE, [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k, int n) noexcept {
        if (n != UTEMP) {
            Real volinv = 1.0 / vol(i, j, k);

            // if (n == NUM_STATE-1 && (i == 64 || i == 65)) {
            //         AllPrint() << "Flux[" << i << "] = ";
            //         for (int m = 0; m < NUM_STATE; ++m) {
            //             AllPrint() << update(i,j,k,m) + (flux0(i,j,k,m)*area0(i,j,k) - flux0(i+1,j,k,m)*area0(i+1,j,k))*volinv << ", ";
            //         }
            //         AllPrint() << std::endl;
            //     }

            update(i, j, k, n) +=
                (flux0(i, j, k, n) * area0(i, j, k) - flux0(i + 1, j, k, n) * area0(i + 1, j, k)
#if AMREX_SPACEDIM >= 2
                 + flux1(i, j, k, n) * area1(i, j, k) - flux1(i, j + 1, k, n) * area1(i, j + 1, k)
#endif
#if AMREX_SPACEDIM == 3
                 + flux2(i, j, k, n) * area2(i, j, k) - flux2(i, j, k + 1, n) * area2(i, j, k + 1)
#endif
                     ) *
                volinv;

            if (n == NUM_STATE - 1 && (i == 64 || i == 65)) {
                AllPrint() << "update[" << i << "] = ";
                for (auto m = 0; m < NUM_STATE; ++m) {
                    AllPrint() << update(i, j, k, m) << ", ";
                }
                AllPrint() << std::endl;
            }
        }
    });
}