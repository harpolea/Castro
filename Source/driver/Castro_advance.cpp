
#include "Castro.H"
#include "Castro_F.H"

#include <cmath>
#include <climits>

using std::string;
using namespace amrex;

Real
Castro::advance (Real time,
                 Real dt,
                 int  amr_iteration,
                 int  amr_ncycle)

  // the main driver for a single level.  This will do either the SDC
  // algorithm or the Strang-split reactions algorithm.
  //
  // arguments:
  //    time          : the current simulation time
  //    dt            : the timestep to advance (e.g., go from time to
  //                    time + dt)
  //    amr_iteration : where we are in the current AMR subcycle.  Each
  //                    level will take a number of steps to reach the
  //                    final time of the coarser level below it.  This
  //                    counter starts at 1
  //    amr_ncycle    : the number of subcycles at this level

{
    BL_PROFILE("Castro::advance()");

    if (do_ctu) {
        std::cout << "Use method of lines, not ctu!";
        exit(EXIT_FAILURE);
    }

    Real dt_new = dt;

    initialize_advance(time, dt, amr_iteration, amr_ncycle);

    // Do the advance.

      for (int iter = 0; iter < MOL_STAGES; ++iter) {
    	mol_iteration = iter;
    	dt_new = do_advance(time, dt, amr_iteration, amr_ncycle);
      }

      if (cfl_violation && hard_cfl_limit)
        amrex::Abort("CFL is too high at this level -- go back to a checkpoint and restart with lower cfl number");

        // If we didn't kill the job, reset the violation counter.

    cfl_violation = 0;

    // Check to see if this advance violated certain stability criteria.
    // If so, get a new timestep and do subcycled advances until we reach
    // t = time + dt.

    if (use_retry)
        dt_new = std::min(dt_new, retry_advance(time, dt, amr_iteration, amr_ncycle));

    if (use_post_step_regrid)
	   check_for_post_regrid(time + dt);

#ifdef AUX_UPDATE
    advance_aux(time, dt);
#endif

    finalize_advance(time, dt, amr_iteration, amr_ncycle);

    return dt_new;
}

Real
Castro::do_advance (Real time,
                    Real dt,
                    int  amr_iteration,
                    int  amr_ncycle)
{

  // this routine will advance the old state data (called S_old here)
  // to the new time, for a single level.  The new data is called
  // S_new here.  The update includes reactions (if we are not doing
  // SDC), hydro, and the source terms.

    BL_PROFILE("Castro::do_advance()");

    const Real prev_time = state[State_Type].prevTime();
    const Real  cur_time = state[State_Type].curTime();

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

    // Perform initialization steps.
    initialize_do_advance(time, dt, amr_iteration, amr_ncycle);

    // Check for NaN's.
    check_for_nan(S_old);

    // Since we are Strang splitting the reactions, do them now (only
    // for first stage of MOL)
    if (mol_iteration == 0) {

      // Initialize the new-time data. This copy needs to come after the
      // reactions.

      MultiFab::Copy(S_new, Sborder, 0, 0, NUM_STATE, S_new.nGrow());

      // Construct the old-time sources from Sborder.  For both CTU
      // and integration, this will already be applied to S_new (with
      // full dt weighting), to be correctly later.  Note: this
      // implies that we do not use the sources as integrate by the
      // MOL integrator.  Also note -- this does not affect the
      // prediction of the interface state, an explict source will be
      // traced there as needed.
      do_old_sources(prev_time, dt, amr_iteration, amr_ncycle);

      // store the result of the burn and old-time sources in Sburn for later stages
      MultiFab::Copy(Sburn, S_new, 0, 0, NUM_STATE, S_new.nGrow());
    }

    cons_to_prim(time);

    check_for_cfl_violation(dt);

    if (cfl_violation)
        return dt;

    // Do the hydro update.  We build directly off of Sborder, which
    // is the state that has already seen the burn
    construct_mol_hydro_source(time, dt);

    // For MOL integration, we are done with this stage, unless it is
    // the last stage

    if (mol_iteration == MOL_STAGES-1) {
      // we just finished the last stage of the MOL integration.
      // Construct S_new now using the weighted sum of the k_mol
      // updates

      // Apply the update -- we need to build on Sburn, so
      // start with that state
      MultiFab::Copy(S_new, Sburn, 0, 0, S_new.nComp(), 0);
      MultiFab::Saxpy(S_new, dt, hydro_source, 0, 0, S_new.nComp(), 0);

      enforce_consistent_e(S_new); // not sure does anything
      clean_state(S_new);
    }

    finalize_do_advance(time, dt, amr_iteration, amr_ncycle);

    return dt;
}



void
Castro::initialize_do_advance(Real time, Real dt, int amr_iteration, int amr_ncycle)
{

    // Reset the change from density resets
    frac_change = 1.e0;

    int finest_level = parent->finestLevel();

    // Reset the grid loss tracking.

    if (track_grid_losses)
      for (int i = 0; i < n_lost; i++)
	     material_lost_through_boundary_temp[i] = 0.0;

    // For the hydrodynamics update we need to have NUM_GROW ghost zones available,
    // but the state data does not carry ghost zones. So we use a FillPatch
    // using the state data to give us Sborder, which does have ghost zones.

    // for Method of lines, our initialization of Sborder depends on
    // which stage in the RK update we are working on

    if (mol_iteration == 0) {

        //   MultiFab& S_old = get_old_data(State_Type);
        //   std::cout << "nan checking\n";
        //   // Check for NaN's.
        //   check_for_nan(S_old);
        //   std::cout << "checked nans\n";

    	// first MOL stage
    	Sborder.define(grids, dmap, NUM_STATE, NUM_GROW);
    	const Real prev_time = state[State_Type].prevTime();
    	expand_state(Sborder, prev_time, NUM_GROW);
    } else {

    	// the initial state for the kth stage follows the Butcher
    	// tableau.  We need to create the proper state starting with
    	// the result after the first dt/2 burn (which we copied into
    	// Sburn) and we need to fill ghost cells.

    	// We'll overwrite S_new with this information, since we don't
    	// need it anymorebuild this state temporarily in S_new (which
    	// is State_Data) to allow for ghost filling.
    	MultiFab& S_new = get_new_data(State_Type);

    	MultiFab::Copy(S_new, Sburn, 0, 0, S_new.nComp(), 0);
    	for (int i = 0; i < mol_iteration; ++i) {
    	  MultiFab::Saxpy(S_new, dt*a_mol[mol_iteration][i], *k_mol[i], 0, 0, S_new.nComp(), 0);
        }

    	Sborder.define(grids, dmap, NUM_STATE, NUM_GROW);
    	const Real new_time = state[State_Type].curTime();
    	expand_state(Sborder, new_time, NUM_GROW);


        // NOTE: this helps
        enforce_min_density(S_new, S_new);
        enforce_consistent_e(S_new);

    }
    // NOTE: this does not
    // enforce_min_density(Sborder, Sborder);
    // enforce_consistent_e(Sborder);
}



void
Castro::finalize_do_advance(Real time, Real dt, int amr_iteration, int amr_ncycle)
{

    // Update the dSdt MultiFab. Since we want (S^{n+1} - S^{n}) / dt,
    // we only need to take twice the new-time source term, since in
    // the predictor-corrector approach, the new-time source term is
    // 1/2 * S^{n+1} - 1/2 * S^{n}. This is untrue in general for the
    // non-momentum sources, but those don't appear in the hydro
    // anyway, and for safety we'll only do this on the momentum
    // terms.

    if (source_term_predictor == 1) {

        MultiFab& dSdt_new = get_new_data(Source_Type);

    	dSdt_new.setVal(0.0, NUM_GROW);

    	for (int n = 0; n < num_src; ++n) {
    	    MultiFab::Add(dSdt_new, *new_sources[n], Xmom, Xmom, 3, 0);
    	}

    	dSdt_new.mult(2.0 / dt);
    }

    Sborder.clear();
    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

    // NOTE: this helps stop it blowing up
    enforce_min_density(S_old, S_new);
    enforce_consistent_e(S_new);
}



void
Castro::initialize_advance(Real time, Real dt, int amr_iteration, int amr_ncycle)
{
    // Pass some information about the state of the simulation to a Fortran module.
    ca_set_amr_info(level, amr_iteration, amr_ncycle, time, dt);

    // Save the current iteration.
    iteration = amr_iteration;

    do_subcycle = false;
    sub_iteration = 0;
    sub_ncycle = 0;
    dt_subcycle = 1.e200;

    cfl_violation = 0;

    if (use_post_step_regrid && level > 0) {

	if (getLevel(level-1).post_step_regrid && amr_iteration == 1) {

            // If the level below this just triggered a special regrid,
            // the coarse contribution to this level's FluxRegister
            // is no longer valid because the grids have, in general, changed.
            // Zero it out, and add them back using the saved copy of the fluxes.

            if (use_post_step_regrid && level > 0)
        	   if (getLevel(level-1).post_step_regrid && amr_iteration == 1)
        	      getLevel(level-1).FluxRegCrseInit();

            // If we're coming off a new regrid at the end of the last coarse
            // timestep, then we want to subcycle this timestep at the timestep
            // suggested by this level, since the data on this level will not
            // have been taken into account when calculating the timestep
            // constraint using the coarser data. This is true even if the level
            // previously existed, because in general there can be new data at this
            // level as a result of the regrid.

            // This step MUST be done before the time level swap because estTimeStep
            // looks at the "new" time data for calculating the timestep constraint.
            // It should also be done before the call to ca_set_amr_info since estTimeStep
            // temporarily resets the level data.

            dt_subcycle = estTimeStep(dt);

            if (dt_subcycle < dt) {

                sub_ncycle = ceil(dt / dt_subcycle);

                if (ParallelDescriptor::IOProcessor()) {
                    std::cout << std::endl;
                    std::cout << "  Subcycling with maximum dt = " << dt_subcycle << " at level " << level
                              << " to avoid timestep constraint violations after a post-timestep regrid."
                              << std::endl << std::endl;
                }
                do_subcycle = true;
            }
        }
    }

    // Pass some information about the state of the simulation to a Fortran module.

    ca_set_amr_info(level, amr_iteration, amr_ncycle, time, dt);

    // The option of whether to do a multilevel initialization is
    // controlled within the radiation class.  This step belongs
    // before the swap.

    // Swap the new data from the last timestep into the old state data.

    for (int k = 0; k < num_state_type; k++) {

    	// The following is a hack to make sure that we only
    	// ever have new data for a few state types that only
    	// ever need new time data; by doing a swap now, we'll
    	// guarantee that allocOldData() does nothing. We do
    	// this because we never need the old data, so we
    	// don't want to allocate memory for it.

    	if (k == Source_Type)
    	    state[k].swapTimeLevels(0.0);

    	state[k].allocOldData();
    	state[k].swapTimeLevels(dt);
    }

    // Ensure data is valid before beginning advance. This addresses
    // the fact that we may have new data on this level that was interpolated
    // from a coarser level, and the interpolation in general cannot be
    // trusted to respect the consistency between certain state variables
    // (e.g. UEINT and UEDEN) that we demand in every zone.

    // NOTE: don't do this
    //enforce_consistent_e(get_old_data(State_Type));
    clean_state(get_old_data(State_Type));

    // Make a copy of the MultiFabs in the old and new state data in case we may do a retry.

    if (use_retry) {
      // Store the old and new time levels.

      for (int k = 0; k < num_state_type; k++) {

    	prev_state[k].reset(new StateData());

    	StateData::Initialize(*prev_state[k], state[k]);
      }
    }

    if (!(keep_sources_until_end || (do_reflux && update_sources_after_reflux))) {

      // These arrays hold all source terms that update the state.  We
      // create and destroy them in the advance to save memory outside
      // of the advance, unless we keep_sources_until_end (for
      // diagnostics) or update the sources after reflux.  In those
      // cases, we initialize these sources in Castro::initMFs and
      // keep them for the life of the simulation.

      for (int n = 0; n < num_src; ++n) {
    	old_sources[n].reset(new MultiFab(grids, dmap, NUM_STATE, NUM_GROW));
    	new_sources[n].reset(new MultiFab(grids, dmap, NUM_STATE, get_new_data(State_Type).nGrow()));
      }

      // This array holds the hydrodynamics update.

      hydro_source.define(grids,dmap,NUM_STATE,0);
    }

    // This array holds the sum of all source terms that affect the
    // hydrodynamics.  If we are doing the source term predictor,
    // we'll also use this after the hydro update to store the sum of
    // the new-time sources, so that we can compute the time
    // derivative of the source terms.

    sources_for_hydro.define(grids,dmap,NUM_STATE,NUM_GROW);

    q.define(grids, dmap, QVAR, NUM_GROW);

    qaux.define(grids, dmap, NQAUX, NUM_GROW);

    src_q.define(grids, dmap, QVAR, NUM_GROW);

    // if we are not doing CTU advection, then we are doing a method
    // of lines, and need storage for hte intermediate stages
    k_mol.resize(MOL_STAGES);
    for (int n = 0; n < MOL_STAGES; ++n) {
        k_mol[n].reset(new MultiFab(grids, dmap, NUM_STATE, 0));
        k_mol[n]->setVal(0.0);
    }

    // for the post-burn state
    Sburn.define(grids, dmap, NUM_STATE, 0);

    // Zero out the current fluxes.

    for (int dir = 0; dir < 3; ++dir)
	   fluxes[dir]->setVal(0.0);

#if (BL_SPACEDIM <= 2)
    if (!Geometry::IsCartesian())
	   P_radial.setVal(0.0);
#endif

    mass_fluxes.resize(3);

    for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
	    mass_fluxes[dir].reset(new MultiFab(getEdgeBoxArray(dir), dmap, 1, 0));
        mass_fluxes[dir]->setVal(0.0);
    }

    for (int dir = BL_SPACEDIM; dir < 3; ++dir) {
    	mass_fluxes[dir].reset(new MultiFab(get_new_data(State_Type).boxArray(), dmap, 1, 0));
        mass_fluxes[dir]->setVal(0.0);
    }
}



void
Castro::finalize_advance(Real time, Real dt, int amr_iteration, int amr_ncycle)
{

    // Add the material lost in this timestep to the cumulative losses.
    if (track_grid_losses) {

      ParallelDescriptor::ReduceRealSum(material_lost_through_boundary_temp, n_lost);

      for (int i = 0; i < n_lost; i++)
    	material_lost_through_boundary_cumulative[i] += material_lost_through_boundary_temp[i];
    }

    if (do_reflux) {
    	FluxRegCrseInit();
    	FluxRegFineAdd();
    }

    Real cur_time = state[State_Type].curTime();
    set_special_tagging_flag(cur_time);

    if (!(keep_sources_until_end || (do_reflux && update_sources_after_reflux))) {

    	amrex::FillNull(old_sources);
    	amrex::FillNull(new_sources);
    	hydro_source.clear();
    }

    sources_for_hydro.clear();
    q.clear();
    qaux.clear();
    src_q.clear();

    amrex::FillNull(prev_state);

    k_mol.clear();
    Sburn.clear();

    // NOTE: doesn't help
    // MultiFab& S_old = get_old_data(State_Type);
    // MultiFab& S_new = get_new_data(State_Type);
    //
    // enforce_min_density(S_old, S_new);
    // enforce_consistent_e(S_new);
}



Real
Castro::retry_advance(Real time, Real dt, int amr_iteration, int amr_ncycle)
{

    Real dt_new = 1.e200;
    Real dt_sub = 1.e200;

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

    const Real* dx = geom.CellSize();

#ifdef _OPENMP
#pragma omp parallel reduction(min:dt_sub)
#endif
    for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

    	const int* lo = bx.loVect();
    	const int* hi = bx.hiVect();

    	ca_check_timestep(BL_TO_FORTRAN_3D(S_old[mfi]),
    			  BL_TO_FORTRAN_3D(S_new[mfi]),
    			  ARLIM_3D(lo), ARLIM_3D(hi), ZFILL(dx),
    			  &dt, &dt_subcycle);
    }

    if (retry_neg_dens_factor > 0.0) {

        // Negative density criterion
    	// Reset so that the desired maximum fractional change in density
    	// is not larger than retry_neg_dens_factor.

        ParallelDescriptor::ReduceRealMin(frac_change);

    	if (frac_change < 0.0)
    	  dt_subcycle = std::min(dt_subcycle, dt * -(retry_neg_dens_factor / frac_change));
    }

    ParallelDescriptor::ReduceRealMin(dt_sub);

    // Do the retry if the suggested timestep is smaller than the actual one.
    // A user-specified tolerance parameter can be used here to prevent
    // retries that are caused by small differences.

    if (dt_sub * (1.0 + retry_tolerance) < std::min(dt, dt_subcycle)) {

        dt_subcycle = dt_sub;

    	if (verbose && ParallelDescriptor::IOProcessor()) {
    	  std::cout << std::endl;
    	  std::cout << "  Timestep " << dt << " rejected at level " << level << "." << std::endl;
    	  std::cout << "  Performing a retry, with subcycled timesteps of maximum length dt = " << dt_subcycle << std::endl;
    	  std::cout << std::endl;
    	}

    	// Restore the original values of the state data.

    	for (int k = 0; k < num_state_type; k++) {

    	  if (prev_state[k]->hasOldData())
    	      state[k].copyOld(*prev_state[k]);

    	  state[k].setTimeLevel(time + dt, dt, 0.0);
    	}
	}

    return dt_new;
}
