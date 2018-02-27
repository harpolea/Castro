/* Implementations of functions in Problem.H go here */

#include "Castro.H"
#include "Castro_F.H"
#include "Problem_F.H"

using namespace amrex;

void Castro::problem_post_timestep(MultiFab& S_new) {

    if (startup_error_removed) return;

    if (level == 0)
    {
        Real dtlev = parent->dtLevel(0);
        Real cumtime = parent->cumTime() + dtlev;

    	if (cumtime > smooth_time) {
            // do smoothing
            const Real* dx        = geom.CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi) {

                const Box& bx = mfi.tilebox();

            	const int* lo = bx.loVect();
            	const int* hi = bx.hiVect();
                const RealBox& pbx  = RealBox(bx,geom.CellSize(),geom.ProbLo());
                const Real* xlo     = pbx.lo();

            	smooth_initial_data(BL_TO_FORTRAN_3D(S_new[mfi]),
            			  ARLIM_3D(lo), ARLIM_3D(hi), ZFILL(dx),
            			  ZFILL(xlo));
            }

            avgDown();

            startup_error_removed = true;
        }
    }
}
