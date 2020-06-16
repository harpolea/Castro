#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

// this can't be a member function as it doesn't play nicely with the function pointer
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real f_p(Real p, const Real* U_zone);

void Castro::LorentzFac(const Box& bx, Array4<Real const> const& vel, Array4<Real> const& W) {
    ParallelFor(bx, [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept {
        const Real v2 = vel(i, j, k, 0) * vel(i, j, k, 0) + vel(i, j, k, 1) * vel(i, j, k, 1) +
                        vel(i, j, k, 2) * vel(i, j, k, 2);
        W(i, j, k) = 1.0_rt / std::sqrt(1.0_rt - v2);
    });
}

AMREX_GPU_HOST_DEVICE void Castro::ConsToPrim(Real* q_zone, Real* U_zone) {

    for (auto n = 0; n < NUM_STATE; ++n) {
        q_zone[n] = 0.0_rt;
    }
    // find pressure using root finder
    Real a = U_zone[UMX] * U_zone[UMX] + U_zone[UMY] * U_zone[UMY] + U_zone[UMZ] * U_zone[UMZ] -
             U_zone[UEDEN] - U_zone[URHO];
    Real p_lo = amrex::max(std::sqrt(amrex::max(a, 0.0_rt)), small_pres);
    Real p_hi = amrex::max((eos_gamma - 1.0_rt) * U_zone[UEDEN] * 1.001_rt, small_pres);

    if (f_p(p_lo, U_zone) * f_p(p_hi, U_zone) > 0.0_rt) {
        p_hi *= 10.0_rt;
        p_lo *= 0.1_rt;
    }
    
    Real p = BrentRootFinder(p_lo, p_hi, f_p, U_zone);

    AllPrint() << "U_zone = (" << U_zone[URHO] << ", " << U_zone[UMX] << ", " << U_zone[UMY] << ", " << U_zone[UEDEN] << ", " << U_zone[UFS] << "), p_lo = " << p_lo << ", p_hi = " << p_hi << ", pbar = " << p << std::endl;

    q_zone[QPRES] = p;
    q_zone[QU] = U_zone[UMX] / (U_zone[UEDEN] + p);
    q_zone[QV] = U_zone[UMY] / (U_zone[UEDEN] + p);
    q_zone[QW] = U_zone[UMZ] / (U_zone[UEDEN] + p);

    Real W_star = 1.0_rt / std::sqrt(1 - q_zone[QU] * q_zone[QU] - q_zone[QV] * q_zone[QV] -
                                     q_zone[QW] * q_zone[QW]);

    q_zone[QRHO] = U_zone[URHO] / W_star;
    q_zone[QREINT] = (U_zone[UEDEN] - U_zone[URHO] * W_star + p * (1.0_rt - W_star * W_star)) /
                     (W_star * W_star);

    q_zone[QTEMP] = U_zone[UTEMP];
}

AMREX_GPU_HOST_DEVICE void Castro::PrimToCons(Real* q_zone, Real* U_zone) {
    for (auto n = 0; n < NUM_STATE; ++n) {
        U_zone[n] = 0.0_rt;
    }

    Real v2 = q_zone[QU] * q_zone[QU] + q_zone[QV] * q_zone[QV] + q_zone[QW] * q_zone[QW];
    Real W = 1.0_rt / std::sqrt(1.0_rt - v2);
    Real h = 1.0_rt + q_zone[QREINT] / q_zone[QRHO] + q_zone[QPRES] / q_zone[QRHO];

    U_zone[URHO] = q_zone[QRHO] * W;

    U_zone[UMX] = U_zone[URHO] * h * W * q_zone[QU];
    U_zone[UMY] = U_zone[URHO] * h * W * q_zone[QV];
    U_zone[UMZ] = U_zone[URHO] * h * W * q_zone[QW];

    U_zone[UEDEN] = U_zone[URHO] * h * W - q_zone[QPRES] - U_zone[URHO];
    U_zone[UEINT] = U_zone[UEDEN];
    U_zone[UTEMP] = q_zone[QTEMP];
}

AMREX_GPU_HOST_DEVICE void Castro::Flux(Real* F, const Real* q_zone, const Real* U_zone,
                                        const int dir) {
    for (auto n = 0; n < NUM_STATE; ++n) {
        F[n] = 0.0_rt;
    }

    F[URHO] = U_zone[URHO] * q_zone[QU + dir];
    F[UMX] = U_zone[UMX] * q_zone[QU + dir];
    F[UMY] = U_zone[UMY] * q_zone[QU + dir];
    F[UMZ] = U_zone[UMZ] * q_zone[QU + dir];
    F[UMX + dir] += q_zone[QPRES];
    F[UEDEN] = U_zone[UMX + dir] - U_zone[URHO] * q_zone[QU + dir];
    F[UEINT] = F[UEDEN];
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real f_p(Real pbar, const Real* U_zone) {
    Real u_star = U_zone[UMX] / (U_zone[UEDEN] + pbar + U_zone[URHO]);
    Real v_star = U_zone[UMY] / (U_zone[UEDEN] + pbar + U_zone[URHO]);
    Real w_star = U_zone[UMZ] / (U_zone[UEDEN] + pbar + U_zone[URHO]);

    Real W_star = 1.0_rt / std::sqrt(1 - u_star * u_star - v_star * v_star - w_star * w_star);

    Real eps_star = (U_zone[UEDEN] + U_zone[URHO] * (1.0_rt - W_star) + pbar * (1.0_rt - W_star * W_star)) /
                    (U_zone[URHO] * W_star);

    Real rho_star = U_zone[URHO] / W_star;

    Real p = (eos_gamma - 1.0_rt) * rho_star * eps_star;

    // return p - pbar;

    eos_t eos_state;
    eos_state.rho = rho_star;
    eos_state.e = eps_star;
    eos_state.xn[0] = 1.0_rt;
    for (auto n = 1; n < NumSpec; ++n) {
        eos_state.xn[n] = 0.0_rt;
    }

    eos(eos_input_re, eos_state);

    return eos_state.p - pbar;
}

AMREX_GPU_HOST_DEVICE
Real Castro::BrentRootFinder(const Real lo, const Real hi, RootFindFunc func,
                             const Real* args) noexcept {
    const Real tol = 1.e-12_rt;
    const int MAXITER = 100;
    const Real EPS = 3.0e-15_rt;

    Real p, q, r, s;

    Real a = lo;
    Real b = hi;

    Real fa = func(a, args);
    Real fb = func(b, args);

    if (fb * fa > 0.0_rt) {
        //        AllPrint() << "fa " << fa << " fb " << fb << "\n";
        Error(
            "BrentRootFinder. Root must be bracketed, but instead the supplied end points have the "
            "same sign.");
        return 0.0_rt;
    } else if (fa == 0.0_rt) {
        return a;
    } else if (fb == 0.0_rt) {
        return b;
    } 

    Real c = a;
    Real fc = fa;

    if (Math::abs(fa) < Math::abs(fb)) {
        a = b;
        b = c;
        c = a;
        fa = fb;
        fb = fc;
        fc = fa;
    }

    Real d = 0.0_rt, e = 0.0_rt;
    int i;
    bool mflag = true;
    for (i = 0; i < MAXITER; ++i) {

        //  Convergence check
        Real tol1 = 2.0_rt * EPS * Math::abs(b) + 0.5_rt * tol;

        if (Math::abs(0.5_rt * (b - a)) <= tol1 || fb == 0.0_rt) {
            break;
        }

        if (fa != fb && fb != fc) {
            // inverse quadratic interpolation
            s = a * fb * fc / (fa - fb) / (fa - fc) + b * fa * fc / (fb - fa) / (fb - fc) + c * fa * fb / (fc - fa) / (fc - fb);
        } else {
            // secant method
            s = b - fb * (b - a) / (fb - fa);
        }

        bool condition1 = !((s > (3*a+b)/4.0_rt && s < b) || (s > b && s < (3*a+b)/4.0_rt));
        bool condition2 = mflag && (Math::abs(s-b) >= Math::abs(b-c)/2.0_rt);
        bool condition3 = !mflag && (Math::abs(s-b) >= Math::abs(c-d)/2.0_rt);
        bool condition4 = mflag && (Math::abs(b-c) < EPS);
        bool condition5 = !mflag && (Math::abs(c-d) < EPS);

        if (condition1 || condition2 || condition3 || condition4 || condition5) {
            // bisection
            s = (a + b) / 2.0_rt;
            mflag = true;
        } else {
            mflag = false;
        }

        Real fs = func(s, args);

        d = c;
        c = b;

        if (fa * fs < 0.0_rt) {
            b = s;
        } else {
            a = s;
        }

        // swap a, b
        if (Math::abs(fa) < Math::abs(fb)) {
            s = a;
            a = b;
            b = s;
            fs = fa;
            fa = fb;
            fb = fs;
        }
        

        // //  Convergence check
        // Real tol1 = 2.0_rt * EPS * Math::abs(b) + 0.5_rt * tol;
        // Real xm = 0.5_rt * (c - b);

        // if (Math::abs(xm) <= tol1 || fb == 0.0_rt) {
        //     break;
        // }

        // if (Math::abs(e) >= tol1 && Math::abs(fa) > Math::abs(fb)) {
        //     //  Attempt inverse quadratic interpolation
        //     s = fb / fa;
        //     if (a == c) {
        //         p = 2.0_rt * xm * s;
        //         q = 1.0_rt - s;
        //     } else {
        //         q = fa / fc;
        //         r = fb / fc;
        //         p = s * (2.0_rt * xm * q * (q - r) - (b - a) * (r - 1.0_rt));
        //         q = (q - 1.0_rt) * (r - 1.0_rt) * (s - 1.0_rt);
        //     }

        //     //  Check whether in bounds
        //     if (p > 0.0_rt) q = -q;

        //     p = Math::abs(p);

        //     if (2.0 * p < min(3.0 * xm * q - Math::abs(tol1 * q), 1.0 * Math::abs(e * q))) {
        //         //  Accept interpolation
        //         e = d;
        //         d = p / q;
        //     } else {
        //         //  Interpolation failed, use bisection
        //         d = xm;
        //         e = d;
        //     }
        // } else {
        //     //  Bounds decreasing too slowly, use bisection
        //     d = xm;
        //     e = d;
        // }

        // //  Move last best guess to a
        // a = b;
        // fa = fb;

        // //  Evaluate new trial root
        // if (Math::abs(d) > tol1) {
        //     b += d;
        // } else {
        //     if (xm < 0.0_rt) {
        //         b -= tol1;
        //     } else {
        //         b += tol1;
        //     }
        // }

        // fb = func(b, args);
    }

    if (i >= MAXITER) {
        Error("BrentRootFinder: exceeding maximum iterations.");
    }

    return b;
}