#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real f_p(amrex::Real p, const amrex::Real* U_zone);

void Castro::LorentzFac(const Box& bx, Array4<Real const> const& vel, Array4<Real> const& W) {
    amrex::ParallelFor(bx, [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept {
        const Real v2 = vel(i, j, k, 0) * vel(i, j, k, 0) + vel(i, j, k, 1) * vel(i, j, k, 1) +
                        vel(i, j, k, 2) * vel(i, j, k, 2);
        W(i, j, k) = 1.0_rt / std::sqrt(1.0_rt - v2);
    });
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Castro::ConsToPrim(Real* q_zone, Real* U_zone) 
{
    // find pressure using root finder
    Real p_lo = amrex::max(amrex::max(std::abs(U_zone[UMX]), amrex::max(std::abs(U_zone[UMY]), std::abs(U_zone[UMZ]))) - U_zone[UEDEN], 0.0_rt);
    Real p_hi = amrex::max(std::abs(U_zone[UMX]), amrex::max(std::abs(U_zone[UMY]), std::abs(U_zone[UMZ]))); // FIXME: what should this be???

    if (f_p(p_lo, U_zone) * f_p(p_hi, U_zone) > 0.0_rt) {
        p_hi *= 10.0_rt;
        p_lo *= 0.1_rt;
    }
    Real p = BrentRootFinder(p_lo, p_hi, f_p, U_zone);

    q_zone[QPRES] = p;
    q_zone[QU] = U_zone[UMX] / (U_zone[UEDEN] + p);
    q_zone[QV] = U_zone[UMY] / (U_zone[UEDEN] + p);
    q_zone[QW] = U_zone[UMZ] / (U_zone[UEDEN] + p);
    
    Real W_star = 1.0_rt / std::sqrt(1 - q_zone[QU]*q_zone[QU] - 
        q_zone[QV]*q_zone[QV] - q_zone[QW]*q_zone[QW]);
    
    q_zone[QRHO] = U_zone[URHO] / W_star;
    q_zone[QREINT] = (U_zone[UEDEN] - U_zone[URHO] * W_star + p * (1.0_rt - W_star*W_star)) / (W_star*W_star);
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Castro::PrimToCons(Real* q_zone, Real* U_zone) 
{
    Real v2 = q_zone[QU] * q_zone[QU] + q_zone[QV] * q_zone[QV] + q_zone[QW] * q_zone[QW];
    Real W = 1.0_rt / std::sqrt(1.0_rt - v2);
    Real h = (q_zone[QREINT] / q_zone[QRHO] + q_zone[QPRES]) / q_zone[QRHO];

    U_zone[URHO] = q_zone[QRHO] * W;

    U_zone[UMX] = U_zone[URHO] * h * W * q_zone[QU];
    U_zone[UMY] = U_zone[URHO] * h * W * q_zone[QV];
    U_zone[UMZ] = U_zone[URHO] * h * W * q_zone[QW];

    U_zone[UEDEN] = U_zone[URHO] * h * W - q_zone[QPRES] - U_zone[URHO];
    U_zone[UEINT] = U_zone[UEDEN];
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real f_p(amrex::Real pbar, const amrex::Real* U_zone)
{
    Real u_star = U_zone[UMX] / (U_zone[UEDEN] + pbar);
    Real v_star = U_zone[UMY] / (U_zone[UEDEN] + pbar);
    Real w_star = U_zone[UMZ] / (U_zone[UEDEN] + pbar);

    Real W_star = 1.0_rt / std::sqrt(1 - u_star*u_star - v_star*v_star - w_star*w_star);

    Real eps_star = (U_zone[UEDEN] - U_zone[URHO] * W_star + 
        pbar * (1.0_rt - W_star*W_star)) / (U_zone[URHO] * W_star);

    Real rho_bar = U_zone[URHO] / W_star;

    eos_t eos_state;
    eos_state.rho = rho_bar;
    eos_state.e = eps_star;
    eos_state.xn[0] = 1.0_rt;
    for (auto n = 1; n < NumSpec; ++n) {
        eos_state.xn[n] = 0.0_rt;
    }

    eos(eos_input_rt, eos_state);

    return eos_state.p - pbar;
}


AMREX_GPU_HOST_DEVICE
Real Castro::BrentRootFinder(const amrex::Real lo, const amrex::Real hi,
                             RootFindFunc func, const amrex::Real* args) noexcept 
{
    const Real tol = 1.e-12_rt;
    const int MAXITER = 100;
    const Real EPS = 3.0e-15_rt;

    Real p, q, r, s;

    Real a = lo;
    Real b = hi;

    Real fa = func(a, args);
    Real fb = func(b, args);
    Real c = b;
    Real fc = fb;

    if (fb * fa > 0.0_rt) {
        //        amrex::AllPrint() << "fa " << fa << " fb " << fb << "\n";
        amrex::Error(
            "BrentRootFinder. Root must be bracketed, but instead the supplied end points have the "
            "same sign.");
        return 0.0_rt;
    } else if (fa == 0.0_rt) {
        return a;
    } else if (fb == 0.0_rt) {
        return b;
    }

    Real d = 0.0_rt, e = 0.0_rt;
    int i;
    for (i = 0; i < MAXITER; ++i) {
        if (fb * fc > 0) {
            //  Rename a, b, c and adjust bounding interval d
            c = a;
            fc = fa;
            d = b - a;
            e = d;
        }

        if (amrex::Math::abs(fc) < amrex::Math::abs(fb)) {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }

        //  Convergence check
        Real tol1 = 2.0_rt * EPS * amrex::Math::abs(b) + 0.5_rt * tol;
        Real xm = 0.5_rt * (c - b);

        if (amrex::Math::abs(xm) <= tol1 || fb == 0.0_rt) {
            break;
        }

        if (amrex::Math::abs(e) >= tol1 && amrex::Math::abs(fa) > amrex::Math::abs(fb)) {
            //  Attempt inverse quadratic interpolation
            s = fb / fa;
            if (a == c) {
                p = 2.0_rt * xm * s;
                q = 1.0_rt - s;
            } else {
                q = fa / fc;
                r = fb / fc;
                p = s * (2.0_rt * xm * q * (q - r) - (b - a) * (r - 1.0_rt));
                q = (q - 1.0_rt) * (r - 1.0_rt) * (s - 1.0_rt);
            }

            //  Check whether in bounds
            if (p > 0.0_rt) q = -q;

            p = amrex::Math::abs(p);

            if (2.0 * p < amrex::min(3.0 * xm * q - amrex::Math::abs(tol1 * q),
                                     1.0 * amrex::Math::abs(e * q))) {
                //  Accept interpolation
                e = d;
                d = p / q;
            } else {
                //  Interpolation failed, use bisection
                d = xm;
                e = d;
            }
        } else {
            //  Bounds decreasing too slowly, use bisection
            d = xm;
            e = d;
        }

        //  Move last best guess to a
        a = b;
        fa = fb;

        //  Evaluate new trial root
        if (amrex::Math::abs(d) > tol1) {
            b += d;
        } else {
            if (xm < 0.0_rt) {
                b -= tol1;
            } else {
                b += tol1;
            }
        }

        fb = func(b, args);
    }

    if (i >= MAXITER) {
        amrex::Error("BrentRootFinder: exceeding maximum iterations.");
    }

    return b;
}