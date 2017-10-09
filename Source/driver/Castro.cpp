
#ifndef WIN32
#include <unistd.h>
#endif

#include <iomanip>

#include <algorithm>
#include <cstdio>
#include <vector>
#include <iostream>
#include <string>
#include <ctime>

#include <AMReX_Utility.H>
#include <AMReX_CONSTANTS.H>
#include <Castro.H>
#include <Castro_F.H>
#include <Derive_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_ParmParse.H>
#include <Castro_error_F.H>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace amrex;

bool         Castro::signalStopJob = false;

bool         Castro::dump_old      = false;

int          Castro::verbose       = 0;
ErrorList    Castro::err_list;
int          Castro::radius_grow   = 1;
BCRec        Castro::phys_bc;
int          Castro::NUM_STATE     = -1;
int          Castro::NUM_GROW      = -1;

Real         Castro::frac_change   = 1.e200;

int          Castro::Density       = -1;
int          Castro::Eden          = -1;
int          Castro::Eint          = -1;
int          Castro::Temp          = -1;
int          Castro::Xmom          = -1;
int          Castro::Ymom          = -1;
int          Castro::Zmom          = -1;

int          Castro::NumSpec       = 0;
int          Castro::FirstSpec     = -1;

int          Castro::NumAux        = 0;
int          Castro::FirstAux      = -1;

int          Castro::NumAdv        = 0;
int          Castro::FirstAdv      = -1;

int          Castro::QVAR          = -1;
int          Castro::QRADVAR       = 0;
int          Castro::NQAUX         = -1;
int          Castro::NQ            = -1;

Array<std::string> Castro::source_names;

int          Castro::MOL_STAGES;
Array< Array<Real> > Castro::a_mol;
Array<Real> Castro::b_mol;
Array<Real> Castro::c_mol;

#include <castro_defaults.H>

std::string  Castro::probin_file = "probin";


#if BL_SPACEDIM == 1
IntVect      Castro::hydro_tile_size(1024);
#elif BL_SPACEDIM == 2
IntVect      Castro::hydro_tile_size(1024,16);
#else
IntVect      Castro::hydro_tile_size(1024,16,16);
#endif

// this will be reset upon restart
Real         Castro::previousCPUTimeUsed = 0.0;

Real         Castro::startCPUTime = 0.0;

int          Castro::Knapsack_Weight_Type = -1;
int          Castro::num_state_type = 0;

// Note: Castro::variableSetUp is in Castro_setup.cpp
// variableCleanUp is called once at the end of a simulation
void
Castro::variableCleanUp ()
{
    desc_lst.clear();

    ca_finalize_meth_params();

    network_finalize();
}

void
Castro::read_params ()
{
    static bool done = false;

    if (done) return;

    done = true;

    ParmParse pp("castro");

#include <castro_queries.H>

    pp.query("v",verbose);
    pp.query("sum_interval",sum_interval);

    pp.query("dump_old",dump_old);

    // Get boundary conditions
    Array<int> lo_bc(BL_SPACEDIM), hi_bc(BL_SPACEDIM);
    pp.getarr("lo_bc",lo_bc,0,BL_SPACEDIM);
    pp.getarr("hi_bc",hi_bc,0,BL_SPACEDIM);
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        phys_bc.setLo(i,lo_bc[i]);
        phys_bc.setHi(i,hi_bc[i]);
    }

    //
    // Check phys_bc against possible periodic geometry
    // if periodic, must have internal BC marked.
    //
    if (Geometry::isAnyPeriodic())
    {
        //
        // Do idiot check.  Periodic means interior in those directions.
        //
        for (int dir = 0; dir<BL_SPACEDIM; dir++)
        {
            if (Geometry::isPeriodic(dir))
            {
                if (lo_bc[dir] != Interior)
                {
                    std::cerr << "Castro::read_params:periodic in direction "
                              << dir
                              << " but low BC is not Interior\n";
                    amrex::Error();
                }
                if (hi_bc[dir] != Interior)
                {
                    std::cerr << "Castro::read_params:periodic in direction "
                              << dir
                              << " but high BC is not Interior\n";
                    amrex::Error();
                }
            }
        }
    }
    else
    {
        //
        // Do idiot check.  If not periodic, should be no interior.
        //
        for (int dir=0; dir<BL_SPACEDIM; dir++)
        {
            if (lo_bc[dir] == Interior)
            {
                std::cerr << "Castro::read_params:interior bc in direction "
                          << dir
                          << " but not periodic\n";
                amrex::Error();
            }
            if (hi_bc[dir] == Interior)
            {
                std::cerr << "Castro::read_params:interior bc in direction "
                          << dir
                          << " but not periodic\n";
                amrex::Error();
            }
        }
    }

    if ( Geometry::IsRZ() && (lo_bc[0] != Symmetry) ) {
        std::cerr << "ERROR:Castro::read_params: must set r=0 boundary condition to Symmetry for r-z\n";
        amrex::Error();
    }

#if (BL_SPACEDIM == 1)
    if ( Geometry::IsSPHERICAL() )
    {
      if ( (lo_bc[0] != Symmetry) && (Geometry::ProbLo(0) == 0.0) )
      {
        std::cerr << "ERROR:Castro::read_params: must set r=0 boundary condition to Symmetry for spherical\n";
        amrex::Error();
      }
    }
#elif (BL_SPACEDIM == 2)
    if ( Geometry::IsSPHERICAL() )
      {
	         amrex::Abort("We don't support spherical coordinate systems in 2D");
      }
#elif (BL_SPACEDIM == 3)
    if ( Geometry::IsRZ() )
      {
	         amrex::Abort("We don't support cylindrical coordinate systems in 3D");
      }
    else if ( Geometry::IsSPHERICAL() )
      {
	         amrex::Abort("We don't support spherical coordinate systems in 3D");
      }
#endif

    // sanity checks

    if (grown_factor < 1)
       amrex::Error("grown_factor must be integer >= 1");

    if (cfl <= 0.0 || cfl > 1.0)
      amrex::Error("Invalid CFL factor; must be between zero and one.");

        if (use_colglaz >= 0)
          {
    	         std::cerr << "ERROR:: use_colglaz is deprecated.  Use riemann_solver instead\n";
    	            amrex::Error();
          }


        // Make sure not to call refluxing if we're not actually doing any hydro.
        if (do_hydro == 0) do_reflux = 0;

        if (max_dt < fixed_dt)
          {
    	         std::cerr << "cannot have max_dt < fixed_dt\n";
    	            amrex::Error();
          }


   StateDescriptor::setBndryFuncThreadSafety(bndry_func_thread_safe);

   ParmParse ppa("amr");
   ppa.query("probin_file",probin_file);

    Array<int> tilesize(BL_SPACEDIM);
    if (pp.queryarr("hydro_tile_size", tilesize, 0, BL_SPACEDIM))
    {
	       for (int i=0; i<BL_SPACEDIM; i++) hydro_tile_size[i] = tilesize[i];
    }

}

Castro::Castro ()
    :
    old_sources(num_src),
    new_sources(num_src),
    prev_state(num_state_type)
{
}

Castro::Castro (Amr&            papa,
                int             lev,
                const Geometry& level_geom,
                const BoxArray& bl,
		        const DistributionMapping& dm,
                Real            time)
    :
    AmrLevel(papa,lev,level_geom,bl,dm,time),
    old_sources(num_src),
    new_sources(num_src),
    prev_state(num_state_type)
{
    buildMetrics();

    initMFs();

    for (int i = 0; i < n_lost; i++) {
      material_lost_through_boundary_cumulative[i] = 0.0;
      material_lost_through_boundary_temp[i] = 0.0;
    }


   // Initialize source term data to zero.

   MultiFab& dSdt_new = get_new_data(Source_Type);
   dSdt_new.setVal(0.0);


    // initialize the Godunov state array used in hydro -- we wait
    // until here so that ngroups is defined (if needed) in
    // rad_params_module
    ca_init_godunov_indices();

    // NQ will be used to dimension the primitive variable state
    // vector it will include the "pure" hydrodynamical variables +
    // any radiation variables
    NQ = QVAR + QRADVAR;

}

Castro::~Castro ()
{
}

void
Castro::buildMetrics ()
{
    const int ngrd = grids.size();

    radius.resize(ngrd);

    const Real* dx = geom.CellSize();

    for (int i = 0; i < ngrd; i++)
    {
        const Box& b = grids[i];
        int ilo      = b.smallEnd(0)-radius_grow;
        int ihi      = b.bigEnd(0)+radius_grow;
        int len      = ihi - ilo + 1;

        radius[i].resize(len);

        Real* rad = radius[i].dataPtr();

        if (Geometry::IsCartesian())
        {
            for (int j = 0; j < len; j++)
            {
                rad[j] = 1.0;
            }
        }
        else
        {
            RealBox gridloc = RealBox(grids[i],geom.CellSize(),geom.ProbLo());

            const Real xlo = gridloc.lo(0) + (0.5 - radius_grow)*dx[0];

            for (int j = 0; j < len; j++)
            {
                rad[j] = xlo + j*dx[0];
            }
        }
    }

    volume.clear();
    volume.define(grids,dmap,1,NUM_GROW);
    geom.GetVolume(volume);

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        area[dir].clear();
    	area[dir].define(getEdgeBoxArray(dir),dmap,1,NUM_GROW);
        geom.GetFaceArea(area[dir],dir);
    }

    dLogArea[0].clear();
#if (BL_SPACEDIM <= 2)
    geom.GetDLogA(dLogArea[0],grids,dmap,0,NUM_GROW);
#endif

    if (level == 0) setGridInfo();
}

// Initialize the MultiFabs and flux registers that live as class members.

void
Castro::initMFs()
{
    fluxes.resize(3);

    for (int dir = 0; dir < BL_SPACEDIM; ++dir)
	   fluxes[dir].reset(new MultiFab(getEdgeBoxArray(dir), dmap, NUM_STATE, 0));

    for (int dir = BL_SPACEDIM; dir < 3; ++dir)
	   fluxes[dir].reset(new MultiFab(get_new_data(State_Type).boxArray(), dmap, NUM_STATE, 0));

#if (BL_SPACEDIM <= 2)
    if (!Geometry::IsCartesian())
	   P_radial.define(getEdgeBoxArray(0), dmap, 1, 0);
#endif

    if (do_reflux && level > 0) {

    	flux_reg.define(grids, dmap, crse_ratio, level, NUM_STATE);
    	flux_reg.setVal(0.0);

#if (BL_SPACEDIM < 3)
    	if (!Geometry::IsCartesian()) {
    	    pres_reg.define(grids, dmap, crse_ratio, level, 1);
    	    pres_reg.setVal(0.0);
    	}
#endif
    }

    // Set the flux register scalings.

    if (do_reflux) {

    	flux_crse_scale = -1.0;
    	flux_fine_scale = 1.0;

    	// The fine pressure scaling depends on dimensionality,
    	// as the dimensionality determines the number of
    	// adjacent zones. In 1D the face is a point so
    	// there's only one fine neighbor for a given coarse
    	// face; in 2D there's crse_ratio[1] faces adjacent
    	// to a face perpendicular to the radial dimension;
    	// and in 3D there would be crse_ratio**2, though
    	// we do not separate the pressure out in 3D. Note
    	// that the scaling by dt has already been handled
    	// in the construction of the P_radial array.

    	// The coarse pressure scaling is the same as for the
    	// fluxes, we want the total refluxing contribution
    	// over the full set of fine timesteps to equal P_radial.

#if (BL_SPACEDIM == 1)
    	pres_crse_scale = 1.0;
    	pres_fine_scale = 1.0;
#elif (BL_SPACEDIM == 2)
    	pres_crse_scale = 1.0;
    	pres_fine_scale = 1.0 / crse_ratio[1];
#endif

    }

    if (keep_sources_until_end || (do_reflux && update_sources_after_reflux)) {

    	// These arrays hold all source terms that update the state.

    	for (int n = 0; n < num_src; ++n) {
    	    old_sources[n].reset(new MultiFab(grids, dmap, NUM_STATE, NUM_GROW));
    	    new_sources[n].reset(new MultiFab(grids, dmap, NUM_STATE, get_new_data(State_Type).nGrow()));
    	}

    	// This array holds the hydrodynamics update.

    	hydro_source.define(grids,dmap,NUM_STATE,0);

    }

    use_post_step_regrid = 0;
}

void
Castro::setTimeLevel (Real time,
                      Real dt_old,
                      Real dt_new)
{
    AmrLevel::setTimeLevel(time,dt_old,dt_new);
}

void
Castro::setGridInfo ()
{

    // Send refinement data to Fortran. We do it here
    // because now the grids have been initialized and
    // we need this data for setting up the problem.
    // Note that this routine will always get called
    // on level 0, even if we are doing a restart,
    // so it is safe to put this here.

    if (level == 0) {

      int max_level = parent->maxLevel();
      int nlevs = max_level + 1;

      Real dx_level[3*nlevs];
      int domlo_level[3*nlevs];
      int domhi_level[3*nlevs];
      int ref_ratio_to_f[3*nlevs];
      int n_error_buf_to_f[nlevs];
      int blocking_factor_to_f[nlevs];

      const Real* dx_coarse = geom.CellSize();

      const int* domlo_coarse = geom.Domain().loVect();
      const int* domhi_coarse = geom.Domain().hiVect();

      for (int dir = 0; dir < 3; dir++) {
    	dx_level[dir] = (ZFILL(dx_coarse))[dir];

    	domlo_level[dir] = (ARLIM_3D(domlo_coarse))[dir];
    	domhi_level[dir] = (ARLIM_3D(domhi_coarse))[dir];

    	// Refinement ratio and error buffer on finest level are meaningless,
    	// and we want them to be zero on the finest level because some
    	// of the algorithms depend on this feature.

    	ref_ratio_to_f[dir + 3 * (nlevs - 1)] = 0;
    	n_error_buf_to_f[nlevs-1] = 0;
      }

      for (int lev = 0; lev <= max_level; lev++)
	     blocking_factor_to_f[lev] = parent->blockingFactor(lev)[0];

      for (int lev = 1; lev <= max_level; lev++) {
	         IntVect ref_ratio = parent->refRatio(lev-1);

    	// Note that we are explicitly calculating here what the
    	// data would be on refined levels rather than getting the
    	// data directly from those levels, because some potential
    	// refined levels may not exist at the beginning of the simulation.

    	for (int dir = 0; dir < 3; dir++)
    	  if (dir < BL_SPACEDIM) {
    	    dx_level[3 * lev + dir] = dx_level[3 * (lev - 1) + dir] / ref_ratio[dir];
    	    int ncell = (domhi_level[3 * (lev - 1) + dir] - domlo_level[3 * (lev - 1) + dir] + 1) * ref_ratio[dir];
    	    domlo_level[3 * lev + dir] = domlo_level[dir];
    	    domhi_level[3 * lev + dir] = domlo_level[3 * lev + dir] + ncell - 1;
    	    ref_ratio_to_f[3 * (lev - 1) + dir] = ref_ratio[dir];
    	  } else {
    	    dx_level[3 * lev + dir] = 0.0;
    	    domlo_level[3 * lev + dir] = 0;
    	    domhi_level[3 * lev + dir] = 0;
    	    ref_ratio_to_f[3 * (lev - 1) + dir] = 0;
    	  }

    	n_error_buf_to_f[lev - 1] = parent->nErrorBuf(lev - 1);
      }

      ca_set_grid_info(max_level, dx_level, domlo_level, domhi_level,
		       ref_ratio_to_f, n_error_buf_to_f, blocking_factor_to_f);

    }
}

void
Castro::initData ()
{
    BL_PROFILE("Castro::initData()");

    //
    // Loop over grids, call FORTRAN function to init with data.
    //
    int ns          = NUM_STATE;
    const Real* dx  = geom.CellSize();
    MultiFab& S_new = get_new_data(State_Type);
    Real cur_time   = state[State_Type].curTime();

    S_new.setVal(0.);

    // make sure dx = dy = dz -- that's all we guarantee to support
#if (BL_SPACEDIM == 2)
    const Real SMALL = 1.e-13;
    if (fabs(dx[0] - dx[1]) > SMALL*dx[0])
      {
	         amrex::Abort("We don't support dx != dy");
      }
#elif (BL_SPACEDIM == 3)
    const Real SMALL = 1.e-13;
    if ( (fabs(dx[0] - dx[1]) > SMALL*dx[0]) || (fabs(dx[0] - dx[2]) > SMALL*dx[0]) )
      {
	         amrex::Abort("We don't support dx != dy != dz");
      }
#endif

    ca_set_amr_info(level, -1, -1, -1.0, -1.0);

    if (verbose && ParallelDescriptor::IOProcessor())
       std::cout << "Initializing the data at level " << level << std::endl;

   if (Knapsack_Weight_Type > 0) {
       get_new_data(Knapsack_Weight_Type).setVal(1.0);
   }

#ifdef MAESTRO_INIT
    MAESTRO_init();
#else
    {
       for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
       {
    	  RealBox gridloc = RealBox(grids[mfi.index()],geom.CellSize(),geom.ProbLo());
              const Box& box     = mfi.validbox();
              const int* lo      = box.loVect();
              const int* hi      = box.hiVect();

              int NVAR;
              ca_get_nvar(&NVAR);

              // NOTE: UNCOMMENT FOR PYTHON
              ca_initdata(level, cur_time, ARLIM_3D(lo), ARLIM_3D(hi), ns,
              S_new[mfi].dataPtr(), ARLIM_3D((S_new[mfi]).loVect()), ARLIM_3D((S_new[mfi]).hiVect()), ZFILL(dx),
              ZFILL(gridloc.lo()), ZFILL(gridloc.hi()));

/*
#ifdef DIMENSION_AGNOSTIC
          BL_FORT_PROC_CALL(CA_INITDATA,ca_initdata)
          (level, cur_time, ARLIM_3D(lo), ARLIM_3D(hi), ns,
      	   BL_TO_FORTRAN_3D(S_new[mfi]), ZFILL(dx),
      	   ZFILL(gridloc.lo()), ZFILL(gridloc.hi()));
#else
          BL_FORT_PROC_CALL(CA_INITDATA,ca_initdata)
      	  (level, cur_time, lo, hi, ns,
      	   BL_TO_FORTRAN(S_new[mfi]), dx,
      	   gridloc.lo(), gridloc.hi());
#endif
*/
              // Verify that the sum of (rho X)_i = rho at every cell
    	  const int idx = mfi.tileIndex();

              ca_check_initial_species(ARLIM_3D(lo), ARLIM_3D(hi),
    				   BL_TO_FORTRAN_3D(S_new[mfi]), &idx);
       }

       enforce_consistent_e(S_new);

       // Do a FillPatch so that we can get the ghost zones filled.

       int ng = S_new.nGrow();

       if (ng > 0)
	      AmrLevel::FillPatch(*this, S_new, ng, cur_time, State_Type, 0, S_new.nComp());
    }

#endif // MAESTRO_INIT

    set_special_tagging_flag(cur_time);

    MultiFab& dSdt_new = get_new_data(Source_Type);
    dSdt_new.setVal(0.);

    if (verbose && ParallelDescriptor::IOProcessor())
       std::cout << "Done initializing the level " << level << " data " << std::endl;
}

void
Castro::init (AmrLevel &old)
{
    BL_PROFILE("Castro::init(old)");

    Castro* oldlev = (Castro*) &old;

    //
    // Create new grid data by fillpatching from old.
    //
    Real dt_new    = parent->dtLevel(level);
    Real cur_time  = oldlev->state[State_Type].curTime();
    Real prev_time = oldlev->state[State_Type].prevTime();
    Real dt_old    = cur_time - prev_time;
    setTimeLevel(cur_time,dt_old,dt_new);

    for (int s = 0; s < num_state_type; ++s) {
    	MultiFab& state_MF = get_new_data(s);
    	FillPatch(old, state_MF, state_MF.nGrow(), cur_time, s, 0, state_MF.nComp());
    }

}

//
// This version inits the data on a new level that did not
// exist before regridding.
//
void
Castro::init ()
{
    BL_PROFILE("Castro::init()");

    Real dt        = parent->dtLevel(level);
    Real cur_time  = getLevel(level-1).state[State_Type].curTime();
    Real prev_time = getLevel(level-1).state[State_Type].prevTime();

    Real dt_old = (cur_time - prev_time)/(Real)parent->MaxRefRatio(level-1);

    Real time = cur_time;

    // If we just triggered a regrid, we need to account for the fact that
    // the data on the coarse level has already been advanced.

    if (getLevel(level-1).use_post_step_regrid)
	   time = prev_time;

    setTimeLevel(time,dt_old,dt);

    for (int s = 0; s < num_state_type; ++s) {
    	MultiFab& state_MF = get_new_data(s);
    	FillCoarsePatch(state_MF, 0, time, s, 0, state_MF.nComp());
    }
}

Real
Castro::initialTimeStep ()
{
    Real dummy_dt = 0.0;
    Real init_dt  = 0.0;

    if (initial_dt > 0.0)
    {
       init_dt = initial_dt;
    }
    else
    {
       init_dt = init_shrink*estTimeStep(dummy_dt);
    }

    return init_dt;
}

Real
Castro::estTimeStep (Real dt_old)
{
    BL_PROFILE("Castro::estTimeStep()");

    if (fixed_dt > 0.0)
        return fixed_dt;

    ca_set_amr_info(level, -1, -1, -1.0, -1.0);

    Real estdt = max_dt;

    const MultiFab& stateMF = get_new_data(State_Type);

    const Real* dx = geom.CellSize();

    std::string limiter = "castro.max_dt";

    // Start the hydro with the max_dt value, but divide by CFL
    // to account for the fact that we multiply by it at the end.
    // This ensures that if max_dt is more restrictive than the hydro
    // criterion, we will get exactly max_dt for a timestep.

    Real estdt_hydro = max_dt / cfl;
    if (do_hydro)
    {

    	  // Compute hydro-limited timestep.
    	if (do_hydro)
    	  {

#ifdef _OPENMP
#pragma omp parallel reduction(min:estdt_hydro)
#endif
	    {
	      Real dt = max_dt / cfl;

	      for (MFIter mfi(stateMF,true); mfi.isValid(); ++mfi)
    		{
    		  const Box& box = mfi.tilebox();

    		  ca_estdt(ARLIM_3D(box.loVect()),
                   ARLIM_3D(box.hiVect()),
    			   BL_TO_FORTRAN_3D(stateMF[mfi]),
    			   ZFILL(dx),&dt);
    		}
#ifdef _OPENMP
#pragma omp critical (castro_estdt)
#endif
	      {
		      estdt_hydro = std::min(estdt_hydro,dt);
	      }
	    }
	  }

       ParallelDescriptor::ReduceRealMin(estdt_hydro);
       estdt_hydro *= cfl;
       if (verbose && ParallelDescriptor::IOProcessor())
           std::cout << "...estimated hydro-limited timestep at level " << level << ": " << estdt_hydro << std::endl;

       // Determine if this is more restrictive than the maximum timestep limiting

       if (estdt_hydro < estdt) {
    	 limiter = "hydro";
    	 estdt = estdt_hydro;
       }
    }

    if (verbose && ParallelDescriptor::IOProcessor())
      std::cout << "Castro::estTimeStep (" << limiter << "-limited) at level " << level << ":  estdt = " << estdt << '\n';

    return estdt;
}

void
Castro::computeNewDt (int                   finest_level,
                      int                   sub_cycle,
                      Array<int>&           n_cycle,
                      const Array<IntVect>& ref_ratio,
                      Array<Real>&          dt_min,
                      Array<Real>&          dt_level,
                      Real                  stop_time,
                      int                   post_regrid_flag)
{
    BL_PROFILE("Castro::computeNewDt()");

    //
    // We are at the start of a coarse grid timecycle.
    // Compute the timesteps for the next iteration.
    //
    if (level > 0)
        return;

    int i;

    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (i = 0; i <= finest_level; i++)
    {
        Castro& adv_level = getLevel(i);
        dt_min[i] = adv_level.estTimeStep(dt_level[i]);
    }

    if (fixed_dt <= 0.0)
    {
       if (post_regrid_flag == 1)
       {
          //
          // Limit dt's by pre-regrid dt
          //
          for (i = 0; i <= finest_level; i++)
          {
              dt_min[i] = std::min(dt_min[i],dt_level[i]);
          }
       }
       else
       {
          //
          // Limit dt's by change_max * old dt
          //
          for (i = 0; i <= finest_level; i++)
          {
             if (verbose && ParallelDescriptor::IOProcessor())
                 if (dt_min[i] > change_max*dt_level[i])
                 {
                        std::cout << "Castro::compute_new_dt : limiting dt at level "
                             << i << '\n';
                        std::cout << " ... new dt computed: " << dt_min[i]
                             << '\n';
                        std::cout << " ... but limiting to: "
                             << change_max * dt_level[i] << " = " << change_max
                             << " * " << dt_level[i] << '\n';
                 }
              dt_min[i] = std::min(dt_min[i],change_max*dt_level[i]);
          }
       }
    }

    //
    // Find the minimum over all levels
    //
    for (i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_min[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[State_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

void
Castro::computeInitialDt (int                   finest_level,
                          int                   sub_cycle,
                          Array<int>&           n_cycle,
                          const Array<IntVect>& ref_ratio,
                          Array<Real>&          dt_level,
                          Real                  stop_time)
{
    BL_PROFILE("Castro::computeInitialDt()");

    //
    // Grids have been constructed, compute dt for all levels.
    //
    if (level > 0)
        return;

    int i;

    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    ///TODO/DEBUG: This will need to change for optimal subcycling.
    for (i = 0; i <= finest_level; i++)
    {
        dt_level[i] = getLevel(i).initialTimeStep();
        n_factor   *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_level[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;

    Real cur_time  = state[State_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

void
Castro::post_timestep (int iteration)
{
    BL_PROFILE("Castro::post_timestep()");

    // Pass some information about the state of the simulation to a Fortran module.

    ca_set_amr_info(level, iteration, -1, -1.0, -1.0);

    //
    // Integration cycle on fine level grids is complete
    // do post_timestep stuff here.
    //
    int finest_level = parent->finestLevel();

    // Now do the refluxing. If we're using gravity it
    // will also do the sync solve associated with the reflux.

    if (do_reflux && level < parent->finestLevel())
	   reflux(level, level+1);

    // Ensure consistency with finer grids.

    if (level < finest_level)
	   avgDown();

    MultiFab& S_new = get_new_data(State_Type);

    // Clean up any aberrant state data generated by the reflux and average-down,
    // and then update quantities like temperature to be consistent.

    clean_state(S_new);

    // Flush Fortran output

    if (verbose)
	   flush_output();

#ifdef DO_PROBLEM_POST_TIMESTEP

    // Provide a hook for the user to do things after all of
    // the normal updates have been applied. The user is
    // responsible for any actions after this point, like
    // doing a computeTemp call if they change the state data.

    problem_post_timestep();

#endif

    if (level == 0)
    {
        int nstep = parent->levelSteps(0);
    	Real dtlev = parent->dtLevel(0);
    	Real cumtime = parent->cumTime() + dtlev;

    	bool sum_int_test = false;

    	if (sum_interval > 0) {

    	  if (nstep%sum_interval == 0)
    	    sum_int_test = true;

    	}

    	bool sum_per_test = false;

    	if (sum_per > 0.0) {

    	  const int num_per_old = floor((cumtime - dtlev) / sum_per);
    	  const int num_per_new = floor((cumtime        ) / sum_per);

    	  if (num_per_old != num_per_new)
    	    sum_per_test = true;

    	}
            if (sum_int_test || sum_per_test)
    	       sum_integrated_quantities();
    }
}

void
Castro::post_restart ()
{
   BL_PROFILE("Castro::post_restart()");

   Real cur_time = state[State_Type].curTime();

    set_special_tagging_flag(cur_time);

    // initialize the Godunov state array used in hydro -- we wait
    // until here so that ngroups is defined (if needed) in
    // rad_params_module
    ca_init_godunov_indices();

    // NQ will be used to dimension the primitive variable state
    // vector it will include the "pure" hydrodynamical variables +
    // any radiation variables
    NQ = QVAR + QRADVAR;

#ifdef DO_PROBLEM_POST_RESTART
    problem_post_restart();
#endif

}

void
Castro::postCoarseTimeStep (Real cumtime)
{
    // postCoarseTimeStep() is only called by level 0.
    BL_ASSERT(level == 0);
    AmrLevel::postCoarseTimeStep(cumtime);
}

void
Castro::check_for_post_regrid (Real time)
{

}

void
Castro::post_regrid (int lbase,
                     int new_finest)
{
    fine_mask.clear();
}

void
Castro::post_init (Real stop_time)
{

    BL_PROFILE("Castro::post_init()");

    if (level > 0)
        return;

    //
    // Average data down from finer levels
    // so that conserved data is consistent between levels.
    //
    int finest_level = parent->finestLevel();
    for (int k = finest_level-1; k>= 0; k--)
        getLevel(k).avgDown();

// Allow the user to define their own post_init functions.

#ifdef DO_PROBLEM_POST_INIT

    problem_post_init();

#endif

    int nstep = parent->levelSteps(0);
	Real dtlev = parent->dtLevel(0);
	Real cumtime = parent->cumTime();
	if (cumtime != 0.0) cumtime += dtlev;

	bool sum_int_test = false;

	if (sum_interval > 0) {

	  if (nstep%sum_interval == 0)
	    sum_int_test = true;

	}

	bool sum_per_test = false;

	if (sum_per > 0.0) {

	  const int num_per_old = floor((cumtime - dtlev) / sum_per);
	  const int num_per_new = floor((cumtime        ) / sum_per);

	  if (num_per_old != num_per_new)
	    sum_per_test = true;

	}

        if (sum_int_test || sum_per_test)
	       sum_integrated_quantities();
}

void
Castro::post_grown_restart ()
{
    if (level > 0)
        return;
}

int
Castro::okToContinue ()
{
    if (level > 0)
        return 1;

    int test = 1;

    if (signalStopJob) {
      test = 0;
      if (ParallelDescriptor::IOProcessor())
	     std::cout << " Signalling a stop of the run due to signalStopJob = true." << std::endl;
    }
    else if (parent->dtLevel(0) < dt_cutoff) {
      test = 0;
      if (ParallelDescriptor::IOProcessor())
	     std::cout << " Signalling a stop of the run because dt < dt_cutoff." << std::endl;
    }

    return test;
}

#ifdef AUX_UPDATE
void
Castro::advance_aux(Real time, Real dt)
{
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... special update for auxiliary variables \n";

    MultiFab&  S_old = get_old_data(State_Type);
    MultiFab&  S_new = get_new_data(State_Type);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_old,true); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.tilebox();
        FArrayBox& old_fab = S_old[mfi];
        FArrayBox& new_fab = S_new[mfi];
	    void ca_auxupdate(BL_TO_FORTRAN(old_fab),
			  BL_TO_FORTRAN(new_fab),
			  box.loVect(), box.hiVect(),
			  &dt);
    }
}
#endif


void
Castro::FluxRegCrseInit() {

    if (level == parent->finestLevel()) return;

    Castro& fine_level = getLevel(level+1);

    for (int i = 0; i < BL_SPACEDIM; ++i)
	fine_level.flux_reg.CrseInit(*fluxes[i], i, 0, 0, NUM_STATE, flux_crse_scale);

#if (BL_SPACEDIM <= 2)
    if (!Geometry::IsCartesian())
	fine_level.pres_reg.CrseInit(P_radial, 0, 0, 0, 1, pres_crse_scale);
#endif

}


void
Castro::FluxRegFineAdd() {

    if (level == 0) return;

    for (int i = 0; i < BL_SPACEDIM; ++i)
	flux_reg.FineAdd(*fluxes[i], i, 0, 0, NUM_STATE, flux_fine_scale);

#if (BL_SPACEDIM <= 2)
    if (!Geometry::IsCartesian())
	getLevel(level).pres_reg.FineAdd(P_radial, 0, 0, 0, 1, pres_fine_scale);
#endif

}


void
Castro::reflux(int crse_level, int fine_level)
{
    BL_PROFILE("Castro::reflux()");

    BL_ASSERT(fine_level > crse_level);

    const Real strt = ParallelDescriptor::second();

    FluxRegister* reg;

    for (int lev = fine_level; lev > crse_level; --lev) {

	reg = &getLevel(lev).flux_reg;

	Castro& crse_lev = getLevel(lev-1);
	Castro& fine_lev = getLevel(lev);

	MultiFab& state = crse_lev.get_new_data(State_Type);

	// Clear out the data that's not on coarse-fine boundaries so that this register only
	// modifies the fluxes on coarse-fine interfaces.

	reg->ClearInternalBorders(crse_lev.geom);

	// Trigger the actual reflux on the coarse level now.

	reg->Reflux(state, crse_lev.volume, 1.0, 0, 0, NUM_STATE, crse_lev.geom);

	// Store the density change, for the gravity sync.

	// Also update the coarse fluxes MultiFabs using the reflux data. This should only make
	// a difference if we re-evaluate the source terms later.

	Array<std::unique_ptr<MultiFab> > temp_fluxes(3);

	if (update_sources_after_reflux) {

	    for (int i = 0; i < BL_SPACEDIM; ++i) {
    		temp_fluxes[i].reset(new MultiFab(crse_lev.fluxes[i]->boxArray(),
    						  crse_lev.fluxes[i]->DistributionMap(),
    						  crse_lev.fluxes[i]->nComp(), crse_lev.fluxes[i]->nGrow()));
    		temp_fluxes[i]->setVal(0.0);
	    }
	    for (OrientationIter fi; fi; ++fi) {
    		const FabSet& fs = (*reg)[fi()];
    		int idir = fi().coordDir();
    		fs.copyTo(*temp_fluxes[idir], 0, 0, 0, temp_fluxes[idir]->nComp());
	    }
	    for (int i = 0; i < BL_SPACEDIM; ++i) {
    		MultiFab::Add(*crse_lev.fluxes[i], *temp_fluxes[i], 0, 0, crse_lev.fluxes[i]->nComp(), 0);
    		temp_fluxes[i].reset();
	    }

	    // Reflux into the hydro_source array so that we have the most up-to-date version of it.

	    reg->Reflux(crse_lev.hydro_source, crse_lev.volume, 1.0, 0, 0, NUM_STATE, crse_lev.geom);

	}

	// We no longer need the flux register data, so clear it out.

	reg->setVal(0.0);

#if (BL_SPACEDIM <= 2)
	if (!Geometry::IsCartesian()) {

	    reg = &getLevel(lev).pres_reg;

	    MultiFab dr(crse_lev.grids, crse_lev.dmap, 1, 0);
	    dr.setVal(crse_lev.geom.CellSize(0));

	    reg->ClearInternalBorders(crse_lev.geom);

	    reg->Reflux(state, dr, 1.0, 0, Xmom, 1, crse_lev.geom);

	    if (update_sources_after_reflux) {

    		temp_fluxes[0].reset(new MultiFab(crse_lev.P_radial.boxArray(),
    						  crse_lev.P_radial.DistributionMap(),
    						  crse_lev.P_radial.nComp(), crse_lev.P_radial.nGrow()));
    		temp_fluxes[0]->setVal(0.0);

                    for (OrientationIter fi; fi; ++fi)
    		{
    		    const FabSet& fs = (*reg)[fi()];
    		    int idir = fi().coordDir();
    		    if (idir == 0) {
        			fs.copyTo(*temp_fluxes[idir], 0, 0, 0, temp_fluxes[idir]->nComp());
    		    }
            }

    		MultiFab::Add(crse_lev.P_radial, *temp_fluxes[0], 0, 0, crse_lev.P_radial.nComp(), 0);
    		temp_fluxes[0].reset();

                    reg->Reflux(crse_lev.hydro_source, dr, 1.0, 0, Xmom, 1, crse_lev.geom);

	    }

	    reg->setVal(0.0);

	}
#endif

    }

    // Do the sync solve across all levels.

    // Now subtract the new-time updates to the state data,
    // recompute it, and add it back. This corrects for the fact
    // that the new-time data was calculated the first time around
    // using a state that hadn't yet been refluxed. Note that this
    // needs to come after the gravity sync solve because the gravity
    // source depends on an up-to-date value of phi. We'll do this
    // on the fine level in addition to the coarser levels, because
    // global sources like gravity or source terms that rely on
    // ghost zone fills like diffusion depend on the data in the
    // coarser levels.

    if (update_sources_after_reflux) {

    	for (int lev = fine_level; lev >= crse_level; --lev) {

    	    MultiFab& S_new = getLevel(lev).get_new_data(State_Type);
    	    Real time = getLevel(lev).state[State_Type].curTime();
    	    Real dt = parent->dtLevel(lev);

    	    for (int n = 0; n < num_src; ++n)
                    if (source_flag(n))
    		    getLevel(lev).apply_source_to_state(S_new, *getLevel(lev).new_sources[n], -dt);

    	    // Make the state data consistent with this earlier version before
    	    // recalculating the new-time source terms.

    	    getLevel(lev).clean_state(S_new);

    	    getLevel(lev).do_new_sources(time, dt);

    	}

    }

    if (verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real      end    = ParallelDescriptor::second() - strt;

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(end,IOProc);
        if (ParallelDescriptor::IOProcessor())
            std::cout << "Castro::reflux() at level " << level << " : time = " << end << std::endl;
#ifdef BL_LAZY
	});
#endif
    }
}

void
Castro::avgDown ()
{
    BL_PROFILE("Castro::avgDown()");

  if (level == parent->finestLevel()) return;

  avgDown(State_Type);
  avgDown(Source_Type);
}

void
Castro::normalize_species (MultiFab& S_new)
{
    int ng = S_new.nGrow();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
    {
       const Box& bx = mfi.growntilebox(ng);
       const int idx = mfi.tileIndex();

       ca_normalize_species(BL_TO_FORTRAN_3D(S_new[mfi]),
			    ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()), &idx);
    }
}

void
Castro::enforce_consistent_e (MultiFab& S)
{

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S,true); mfi.isValid(); ++mfi)
    {
        const Box& box     = mfi.tilebox();
        const int* lo      = box.loVect();
        const int* hi      = box.hiVect();

    	const int idx      = mfi.tileIndex();
            ca_enforce_consistent_e(ARLIM_3D(lo), ARLIM_3D(hi), BL_TO_FORTRAN_3D(S[mfi]), &idx);
    }
}

Real
Castro::enforce_min_density (MultiFab& S_old, MultiFab& S_new)
{

    // This routine sets the density in S_new to be larger than the density floor.
    // Note that it will operate everywhere on S_new, including ghost zones.
    // S_old is present so that, after the hydro call, we know what the old density
    // was so that we have a reference for comparison. If you are calling it elsewhere
    // and there's no meaningful reference state, just pass in the same MultiFab twice.

    // The return value is the the negative fractional change in the state that has the
    // largest magnitude. If there is no reference state, this is meaningless.

    Real dens_change = 1.e0;

    MultiFab reset_source;

    if (print_update_diagnostics)
    {

    	// Before we do anything, make a copy of the state.

    	reset_source.define(S_new.boxArray(), S_new.DistributionMap(), S_new.nComp(), 0);

    	MultiFab::Copy(reset_source, S_new, 0, 0, S_new.nComp(), 0);

    }

#ifdef _OPENMP
#pragma omp parallel reduction(min:dens_change)
#endif
    for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi) {

    	const Box& bx = mfi.growntilebox();

    	FArrayBox& stateold = S_old[mfi];
    	FArrayBox& statenew = S_new[mfi];
    	FArrayBox& vol      = volume[mfi];
    	const int idx = mfi.tileIndex();

    	ca_enforce_minimum_density(stateold.dataPtr(), ARLIM_3D(stateold.loVect()), ARLIM_3D(stateold.hiVect()),
    				   statenew.dataPtr(), ARLIM_3D(statenew.loVect()), ARLIM_3D(statenew.hiVect()),
    				   vol.dataPtr(), ARLIM_3D(vol.loVect()), ARLIM_3D(vol.hiVect()),
    				   ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
    				   &dens_change, &verbose, &idx);

    }

    if (print_update_diagnostics)
    {

	// Evaluate what the effective reset source was.

	MultiFab::Subtract(reset_source, S_new, 0, 0, S_old.nComp(), 0);

	bool local = true;
	Array<Real> reset_update = evaluate_source_change(reset_source, 1.0, local);

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
	    ParallelDescriptor::ReduceRealSum(reset_update.dataPtr(), reset_update.size(), ParallelDescriptor::IOProcessorNumber());

	    if (ParallelDescriptor::IOProcessor()) {
    		if (std::abs(reset_update[0]) != 0.0) {
    		    std::cout << std::endl << "  Contributions to the state from negative density resets:" << std::endl;

    		    print_source_change(reset_update);
    		}
	    }

#ifdef BL_LAZY
        });
#endif

    }

    return dens_change;

}

void
Castro::avgDown (int state_indx)
{
    BL_PROFILE("Castro::avgDown(state_indx)");

    if (level == parent->finestLevel()) return;

    Castro& fine_lev = getLevel(level+1);

    const Geometry& fgeom = fine_lev.geom;
    const Geometry& cgeom =          geom;

    MultiFab&  S_crse   = get_new_data(state_indx);
    MultiFab&  S_fine   = fine_lev.get_new_data(state_indx);

    amrex::average_down(S_fine, S_crse,
			 fgeom, cgeom,
			 0, S_fine.nComp(), fine_ratio);
}

void
Castro::allocOldData ()
{
    for (int k = 0; k < num_state_type; k++)
        state[k].allocOldData();
}

void
Castro::removeOldData()
{
    AmrLevel::removeOldData();
}

void
Castro::errorEst (TagBoxArray& tags,
                  int          clearval,
                  int          tagval,
                  Real         time,
                  int          n_error_buf,
                  int          ngrow)
{
    BL_PROFILE("Castro::errorEst()");

    ca_set_amr_info(level, -1, -1, -1.0, -1.0);

    Real t = time;

    // If we are forcing a post-timestep regrid,
    // note that we need to use the new time here,
    // not the old time.

    if (use_post_step_regrid)
	   t = get_state_data(State_Type).curTime();

    // Apply each of the built-in tagging functions.

    for (int j = 0; j < err_list.size(); j++)
	   apply_tagging_func(tags, clearval, tagval, t, j);

    // Now we'll tag any user-specified zones using the full state array.

    const int*  domain_lo = geom.Domain().loVect();
    const int*  domain_hi = geom.Domain().hiVect();
    const Real* dx        = geom.CellSize();
    const Real* prob_lo   = geom.ProbLo();

    MultiFab& S_new = get_new_data(State_Type);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Array<int>  itags;

	for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
	{
	    // tile box
	    const Box&  tilebx  = mfi.tilebox();

            TagBox&     tagfab  = tags[mfi];

	    // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
	    // So we are going to get a temporary integer array.
	    tagfab.get_itags(itags, tilebx);

            // data pointer and index space
	    int*        tptr    = itags.dataPtr();
	    const int*  tlo     = tilebx.loVect();
	    const int*  thi     = tilebx.hiVect();

#ifdef DIMENSION_AGNOSTIC
	    set_problem_tags(tptr,  ARLIM_3D(tlo), ARLIM_3D(thi),
			     BL_TO_FORTRAN_3D(S_new[mfi]),
			     &tagval, &clearval,
			     ARLIM_3D(tilebx.loVect()), ARLIM_3D(tilebx.hiVect()),
			     ZFILL(dx), ZFILL(prob_lo), &time, &level);
#else
	    set_problem_tags(tptr,  ARLIM(tlo), ARLIM(thi),
			     BL_TO_FORTRAN(S_new[mfi]),
			     &tagval, &clearval,
			     tilebx.loVect(), tilebx.hiVect(),
			     dx, prob_lo, &time, &level);
#endif

	    //
	    // Now update the tags in the TagBox.
	    //
            tagfab.tags_and_untags(itags, tilebx);
	}
    }
}



void
Castro::apply_tagging_func(TagBoxArray& tags, int clearval, int tagval, Real time, int j)
{

    const int*  domain_lo = geom.Domain().loVect();
    const int*  domain_hi = geom.Domain().hiVect();
    const Real* dx        = geom.CellSize();
    const Real* prob_lo   = geom.ProbLo();

    for (int j = 0; j < err_list.size(); j++)
    {
        auto mf = derive(err_list[j].name(), time, err_list[j].nGrow());

        BL_ASSERT(mf);

#ifdef _OPENMP
#pragma omp parallel
#endif
	{
	    Array<int>  itags;

	    for (MFIter mfi(*mf,true); mfi.isValid(); ++mfi)
	    {
		// FABs
		FArrayBox&  datfab  = (*mf)[mfi];
		TagBox&     tagfab  = tags[mfi];

		// tile box
		const Box&  tilebx  = mfi.tilebox();

		// physical tile box
		const RealBox& pbx  = RealBox(tilebx,geom.CellSize(),geom.ProbLo());

		//fab box
		const Box&  datbox  = datfab.box();

		// We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
		// So we are going to get a temporary integer array.
		tagfab.get_itags(itags, tilebx);

		// data pointer and index space
		int*        tptr    = itags.dataPtr();
		const int*  tlo     = tilebx.loVect();
		const int*  thi     = tilebx.hiVect();
		//
		const int*  lo      = tlo;
		const int*  hi      = thi;
		//
		const Real* xlo     = pbx.lo();
		//
		Real*       dat     = datfab.dataPtr();
		const int*  dlo     = datbox.loVect();
		const int*  dhi     = datbox.hiVect();
		const int   ncomp   = datfab.nComp();

		err_list[j].errFunc()(tptr, tlo, thi, &tagval,
				      &clearval, dat, dlo, dhi,
				      lo,hi, &ncomp, domain_lo, domain_hi,
				      dx, xlo, prob_lo, &time, &level);
		//
		// Now update the tags in the TagBox.
		//
                tagfab.tags_and_untags(itags, tilebx);
	    }
	}

    }
}



std::unique_ptr<MultiFab>
Castro::derive (const std::string& name,
                Real           time,
                int            ngrow)
{
   return AmrLevel::derive(name,time,ngrow);
}

void
Castro::derive (const std::string& name,
                Real           time,
                MultiFab&      mf,
                int            dcomp)
{

    AmrLevel::derive(name,time,mf,dcomp);
}

void
Castro::network_init ()
{
   ca_network_init();
}

void
Castro::network_finalize ()
{
   ca_network_finalize();
}

void
Castro::extern_init ()
{
  // initialize the external runtime parameters -- these will
  // live in the probin

  if (ParallelDescriptor::IOProcessor()) {
    std::cout << "reading extern runtime parameters ..." << std::endl;
  }

  int probin_file_length = probin_file.length();
  Array<int> probin_file_name(probin_file_length);

  for (int i = 0; i < probin_file_length; i++)
    probin_file_name[i] = probin_file[i];

  ca_extern_init(probin_file_name.dataPtr(),&probin_file_length);
}

void
Castro::reset_internal_energy(MultiFab& S_new)
{

    MultiFab old_state;

    // Make a copy of the state so we can evaluate how much changed.

    if (print_update_diagnostics)
    {
	old_state.define(S_new.boxArray(), S_new.DistributionMap(), S_new.nComp(), 0);
        MultiFab::Copy(old_state, S_new, 0, 0, S_new.nComp(), 0);
    }

    int ng = S_new.nGrow();

    // Ensure (rho e) isn't too small or negative
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(ng);
	    const int idx = mfi.tileIndex();

        ca_reset_internal_e(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
			    BL_TO_FORTRAN_3D(S_new[mfi]),
			    print_fortran_warnings, &idx);
    }

    // Flush Fortran output

    if (verbose)
      flush_output();

    if (print_update_diagnostics)
    {
    	// Evaluate what the effective reset source was.

    	MultiFab reset_source(S_new.boxArray(), S_new.DistributionMap(), S_new.nComp(), 0);

    	MultiFab::Copy(reset_source, S_new, 0, 0, S_new.nComp(), 0);

    	MultiFab::Subtract(reset_source, old_state, 0, 0, old_state.nComp(), 0);

    	bool local = true;
    	Array<Real> reset_update = evaluate_source_change(reset_source, 1.0, local);

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
	    ParallelDescriptor::ReduceRealSum(reset_update.dataPtr(), reset_update.size(), ParallelDescriptor::IOProcessorNumber());

	    if (ParallelDescriptor::IOProcessor()) {
    		if (std::abs(reset_update[Eint]) != 0.0) {
    		    std::cout << std::endl << "  Contributions to the state from negative energy resets:" << std::endl;

    		    print_source_change(reset_update);
    		}
	    }

#ifdef BL_LAZY
	});
#endif
    }
}

void
Castro::computeTemp(MultiFab& State)
{

  reset_internal_energy(State);

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(State,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox();
        const int idx = mfi.tileIndex();
        ca_compute_temp(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
        		BL_TO_FORTRAN_3D(State[mfi]), &idx);
    }
}

void
Castro::set_special_tagging_flag(Real time)
{
   if (!do_special_tagging) return;

   MultiFab& S_new = get_new_data(State_Type);
   Real max_den = S_new.norm0(Density);

   int flag_was_changed = 0;
   ca_set_special_tagging_flag(max_den,&flag_was_changed);
   if (ParallelDescriptor::IOProcessor()) {
      if (flag_was_changed == 1) {
        std::ofstream os("Bounce_time",std::ios::out);
        os << "T_Bounce " << time << std::endl;
        os.close();
      }
   }
}


Real
Castro::getCPUTime()
{

  int numCores = ParallelDescriptor::NProcs();
#ifdef _OPENMP
  numCores = numCores*omp_get_max_threads();
#endif

  Real T = numCores*(ParallelDescriptor::second() - startCPUTime) +
    previousCPUTimeUsed;

  return T;
}


MultiFab&
Castro::build_fine_mask()
{
    BL_ASSERT(level > 0); // because we are building a mask for the coarser level

    if (!fine_mask.empty()) return fine_mask;

    BoxArray baf = parent->boxArray(level);
    baf.coarsen(crse_ratio);

    const BoxArray& bac = parent->boxArray(level-1);
    const DistributionMapping& dmc = parent->DistributionMap(level-1);
    fine_mask.define(bac,dmc,1,0);
    fine_mask.setVal(1.0);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(fine_mask); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = fine_mask[mfi];

    	const std::vector< std::pair<int,Box> >& isects = baf.intersections(fab.box());

    	for (int ii = 0; ii < isects.size(); ++ii)
    	{
    	    fab.setVal(0.0,isects[ii].second,0);
    	}
    }

    return fine_mask;
}

iMultiFab&
Castro::build_interior_boundary_mask (int ng)
{
    for (int i = 0; i < ib_mask.size(); ++i)
    {
    	if (ib_mask[i]->nGrow() == ng) {
    	    return *ib_mask[i];
    	}
    }

    //  If we got here, we need to build a new one
    ib_mask.push_back(std::unique_ptr<iMultiFab>(new iMultiFab(grids, dmap, 1, ng)));

    iMultiFab& imf = *ib_mask.back();

    int ghost_covered_by_valid = 0;
    int other_cells = 1; // uncovered ghost, valid, and outside domain cells are set to 1

    imf.BuildMask(geom.Domain(), geom.periodicity(),
		  ghost_covered_by_valid, other_cells, other_cells, other_cells);

    return imf;
}

// Fill a version of the state with ng ghost zones from the state data.

void
Castro::expand_state(MultiFab& S, Real time, int ng)
{
    BL_ASSERT(S.nGrow() >= ng);

    AmrLevel::FillPatch(*this,S,ng,time,State_Type,0,NUM_STATE);

    clean_state(S);

}


void
Castro::check_for_nan(MultiFab& state, int check_ghost)
{

  int ng = 0;
  if (check_ghost == 1) {
    ng = state.nComp();
  }

  if (state.contains_nan(Density,state.nComp(),ng,true))
    {
      for (int i = 0; i < state.nComp(); i++)
        {
    	  if (state.contains_nan(Density + i, 1, ng, true)) {
    	      std::string abort_string = std::string("State has NaNs in the ") + desc_lst[State_Type].name(i) + std::string(" component::check_for_nan()");
    	      amrex::Abort(abort_string.c_str());
            }
        }
    }
}


// Given State_Type state data, perform a number of cleaning steps to make
// sure the data is sensible. The return value is the same as the return
// value of enforce_min_density.

Real
Castro::clean_state(MultiFab& state) {

    // Enforce a minimum density.

    MultiFab temp_state(state.boxArray(), state.DistributionMap(), state.nComp(), state.nGrow());

    MultiFab::Copy(temp_state, state, 0, 0, state.nComp(), state.nGrow());

    Real frac_change = enforce_min_density(temp_state, state);

    // Ensure all species are normalized.

    normalize_species(state);

    // Compute the temperature (note that this will also reset
    // the internal energy for consistency with the total energy).

    computeTemp(state);

    return frac_change;

}



Real
Castro::clean_state(MultiFab& state, MultiFab& state_old) {

    // Enforce a minimum density.

    Real frac_change = enforce_min_density(state_old, state);

    // Ensure all species are normalized.

    normalize_species(state);

    // Compute the temperature (note that this will also reset
    // the internal energy for consistency with the total energy).

    computeTemp(state);

    return frac_change;

}
