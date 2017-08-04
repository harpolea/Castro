#include "Castro.H"
#include "Castro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

using namespace amrex;

void
Castro::construct_hydro_source(Real time, Real dt)
{

  // this constructs the hydrodynamic source (essentially the flux
  // divergence) using the CTU framework for unsplit hydrodynamics

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... Entering hydro advance" << std::endl << std::endl;

    hydro_source.setVal(0.0);

    // Set up the source terms to go into the hydro.

    sources_for_hydro.setVal(0.0);

    for (int n = 0; n < num_src; ++n)
	MultiFab::Add(sources_for_hydro, *old_sources[n], 0, 0, NUM_STATE, NUM_GROW);

    sources_for_hydro.FillBoundary(geom.periodicity());

#ifndef SDC
    // Optionally we can predict the source terms to t + dt/2,
    // which is the time-level n+1/2 value, To do this we use a
    // lagged predictor estimate: dS/dt_n = (S_n - S_{n-1}) / dt, so
    // S_{n+1/2} = S_n + (dt / 2) * dS/dt_n.

    if (source_term_predictor == 1) {

      MultiFab& dSdt_new = get_new_data(Source_Type);

      dSdt_new.FillBoundary(geom.periodicity());

      MultiFab::Saxpy(sources_for_hydro, 0.5 * dt, dSdt_new, 0, 0, NUM_STATE, NUM_GROW);

    }
#else
    // If we're doing SDC, time-center the source term (using the
    // current iteration's old sources and the last iteration's new
    // sources).

    MultiFab& SDC_source = get_new_data(SDC_Source_Type);

    MultiFab::Add(sources_for_hydro, SDC_source, 0, 0, NUM_STATE, 0);

    sources_for_hydro.FillBoundary(geom.periodicity());

#endif

    int finest_level = parent->finestLevel();

    const Real *dx = geom.CellSize();
    Real courno    = -1.0e+200;

    MultiFab& S_new = get_new_data(State_Type);


    // note: the radiation consup currently does not fill these
    Real mass_lost       = 0.;
    Real xmom_lost       = 0.;
    Real ymom_lost       = 0.;
    Real zmom_lost       = 0.;
    Real eden_lost       = 0.;
    Real xang_lost       = 0.;
    Real yang_lost       = 0.;
    Real zang_lost       = 0.;

    BL_PROFILE_VAR("Castro::advance_hydro_ca_umdrv()", CA_UMDRV);

#ifdef _OPENMP
#pragma omp parallel reduction(+:mass_lost,xmom_lost,ymom_lost,zmom_lost) \
		     reduction(+:eden_lost,xang_lost,yang_lost,zang_lost)
#endif
    {

      FArrayBox flux[BL_SPACEDIM];
#if (BL_SPACEDIM <= 2)
      FArrayBox pradial(Box::TheUnitBox(),1);
#endif
      FArrayBox q, qaux, src_q;

      int priv_nstep_fsp = -1;

      Real cflLoc = -1.0e+200;
      int is_finest_level = (level == finest_level) ? 1 : 0;
      const int*  domain_lo = geom.Domain().loVect();
      const int*  domain_hi = geom.Domain().hiVect();

      for (MFIter mfi(S_new,hydro_tile_size); mfi.isValid(); ++mfi)
      {
	  const Box& bx    = mfi.tilebox();
	  const Box& qbx = amrex::grow(bx, NUM_GROW);

	  const int* lo = bx.loVect();
	  const int* hi = bx.hiVect();

	  FArrayBox &statein  = Sborder[mfi];
	  FArrayBox &stateout = S_new[mfi];

	  FArrayBox &source_in  = sources_for_hydro[mfi];
	  FArrayBox &source_out = hydro_source[mfi];

	  q.resize(qbx, QVAR);
	  qaux.resize(qbx, NQAUX);
	  src_q.resize(qbx, QVAR);

	  // convert the conservative state to the primitive variable state.
	  // this fills both q and qaux.

	  const int idx = mfi.tileIndex();

	  ca_ctoprim(ARLIM_3D(qbx.loVect()), ARLIM_3D(qbx.hiVect()),
		     statein.dataPtr(), ARLIM_3D(statein.loVect()), ARLIM_3D(statein.hiVect()),
		     q.dataPtr(), ARLIM_3D(q.loVect()), ARLIM_3D(q.hiVect()),
		     qaux.dataPtr(), ARLIM_3D(qaux.loVect()), ARLIM_3D(qaux.hiVect()), &idx);

	  // convert the source terms expressed as sources to the conserved state to those
	  // expressed as sources for the primitive state.

	  ca_srctoprim(ARLIM_3D(qbx.loVect()), ARLIM_3D(qbx.hiVect()),
		       q.dataPtr(), ARLIM_3D(q.loVect()), ARLIM_3D(q.hiVect()),
		       qaux.dataPtr(), ARLIM_3D(qaux.loVect()), ARLIM_3D(qaux.hiVect()),
		       source_in.dataPtr(), ARLIM_3D(source_in.loVect()), ARLIM_3D(source_in.hiVect()),
		       src_q.dataPtr(), ARLIM_3D(src_q.loVect()), ARLIM_3D(src_q.hiVect()), &idx);

#ifndef RADIATION

	  // Add in the reactions source term; only done in SDC.

#endif
	  // Allocate fabs for fluxes
	  for (int i = 0; i < BL_SPACEDIM ; i++)  {
	    const Box& bxtmp = amrex::surroundingNodes(bx,i);
	    flux[i].resize(bxtmp,NUM_STATE);

	  }

#if (BL_SPACEDIM <= 2)
	  if (!Geometry::IsCartesian()) {
	    pradial.resize(amrex::surroundingNodes(bx,0),1);
	  }
#endif

	  ca_ctu_update
	    (&is_finest_level, &time,
	     lo, hi, domain_lo, domain_hi,
	     BL_TO_FORTRAN_3D(statein),
	     BL_TO_FORTRAN_3D(stateout),

	     BL_TO_FORTRAN_3D(q),
	     BL_TO_FORTRAN_3D(qaux),
	     BL_TO_FORTRAN_3D(src_q),
	     BL_TO_FORTRAN_3D(source_out),
	     dx, &dt,
	     D_DECL(BL_TO_FORTRAN_3D(flux[0]),
		    BL_TO_FORTRAN_3D(flux[1]),
		    BL_TO_FORTRAN_3D(flux[2])),

#if (BL_SPACEDIM < 3)
	     BL_TO_FORTRAN_3D(pradial),
#endif
	     D_DECL(BL_TO_FORTRAN_3D(area[0][mfi]),
		    BL_TO_FORTRAN_3D(area[1][mfi]),
		    BL_TO_FORTRAN_3D(area[2][mfi])),
#if (BL_SPACEDIM < 3)
	     BL_TO_FORTRAN_3D(dLogArea[0][mfi]),
#endif
	     BL_TO_FORTRAN_3D(volume[mfi]),
	     &cflLoc, verbose,
	     mass_lost, xmom_lost, ymom_lost, zmom_lost,
	     eden_lost, xang_lost, yang_lost, zang_lost);

	  // Store the fluxes from this advance.
	  // For normal integration we want to add the fluxes from this advance
	  // since we may be subcycling the timestep. But for SDC integration
	  // we want to copy the fluxes since we expect that there will not be
	  // subcycling and we only want the last iteration's fluxes.

	  for (int i = 0; i < BL_SPACEDIM ; i++) {
#ifndef SDC
	    (*fluxes    [i])[mfi].plus(    flux[i],mfi.nodaltilebox(i),0,0,NUM_STATE);
#else
	    (*fluxes    [i])[mfi].copy(    flux[i],mfi.nodaltilebox(i),0,mfi.nodaltilebox(i),0,NUM_STATE);
#endif
	  }

#if (BL_SPACEDIM <= 2)
	  if (!Geometry::IsCartesian()) {
#ifndef SDC
	    P_radial[mfi].plus(pradial,mfi.nodaltilebox(0),0,0,1);
#else
	    P_radial[mfi].copy(pradial,mfi.nodaltilebox(0),0,mfi.nodaltilebox(0),0,1);
#endif
	  }
#endif
      } // MFIter loop

#ifdef _OPENMP
#pragma omp critical (hydro_courno)
#endif
      {
	courno = std::max(courno,cflLoc);

      }
    }  // end of omp parallel region

    BL_PROFILE_VAR_STOP(CA_UMDRV);

#ifdef RADIATION
    if (radiation->verbose>=1) {
#ifdef BL_LAZY
      Lazy::QueueReduction( [=] () mutable {
#endif
	  ParallelDescriptor::ReduceIntMax(nstep_fsp, ParallelDescriptor::IOProcessorNumber());
	  if (ParallelDescriptor::IOProcessor() && nstep_fsp > 0) {
	    std::cout << "Radiation f-space advection on level " << level
		      << " takes as many as " << nstep_fsp;
	    if (nstep_fsp == 1) {
	      std::cout<< " substep.\n";
	    }
	    else {
	      std::cout<< " substeps.\n";
	    }
	  }
#ifdef BL_LAZY
	});
#endif
    }
#else
    // Flush Fortran output

    if (verbose)
      flush_output();

    if (track_grid_losses)
    {
	material_lost_through_boundary_temp[0] += mass_lost;
	material_lost_through_boundary_temp[1] += xmom_lost;
	material_lost_through_boundary_temp[2] += ymom_lost;
	material_lost_through_boundary_temp[3] += zmom_lost;
	material_lost_through_boundary_temp[4] += eden_lost;
	material_lost_through_boundary_temp[5] += xang_lost;
	material_lost_through_boundary_temp[6] += yang_lost;
	material_lost_through_boundary_temp[7] += zang_lost;
    }

    if (print_update_diagnostics)
    {

	bool local = true;
	Array<Real> hydro_update = evaluate_source_change(hydro_source, dt, local);

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
	    ParallelDescriptor::ReduceRealSum(hydro_update.dataPtr(), hydro_update.size(), ParallelDescriptor::IOProcessorNumber());

	    if (ParallelDescriptor::IOProcessor())
		std::cout << std::endl << "  Contributions to the state from the hydro source:" << std::endl;

	    print_source_change(hydro_update);

#ifdef BL_LAZY
	});
#endif
      }
#endif

    if (courno > 1.0) {
	std::cout << "WARNING -- EFFECTIVE CFL AT THIS LEVEL " << level << " IS " << courno << '\n';
	if (hard_cfl_limit == 1)
	  amrex::Abort("CFL is too high at this level -- go back to a checkpoint and restart with lower cfl number");
    }

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << std::endl << "... Leaving hydro advance" << std::endl << std::endl;

}



void
Castro::construct_mol_hydro_source(Real time, Real dt)
{

  // this constructs the hydrodynamic source (essentially the flux
  // divergence) using method of lines integration.  The output, as a
  // update to the state, is stored in the k_mol array of multifabs.

  if (verbose && ParallelDescriptor::IOProcessor())
    std::cout << "... hydro MOL stage " << mol_iteration << std::endl;


  // we'll add each stage's contribution to -div{F(U)} as we compute them
  if (mol_iteration == 0) {
    hydro_source.setVal(0.0);
  }

  // Set up the source terms to go into the hydro -- note: the
  // sources_for_hydro MF has ghost zones, but we don't need them
  // here, since sources don't explicitly enter into the prediction
  // for MOL integration

  sources_for_hydro.setVal(0.0);

  for (int n = 0; n < num_src; ++n)
    MultiFab::Add(sources_for_hydro, *old_sources[n], 0, 0, NUM_STATE, 0);

  int finest_level = parent->finestLevel();

  const Real *dx = geom.CellSize();
  Real courno    = -1.0e+200;

  MultiFab& S_new = get_new_data(State_Type);

  MultiFab& k_stage = *k_mol[mol_iteration];

  BL_PROFILE_VAR("Castro::advance_hydro_ca_umdrv()", CA_UMDRV);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {

    FArrayBox flux[BL_SPACEDIM];
#if (BL_SPACEDIM <= 2)
    FArrayBox pradial(Box::TheUnitBox(),1);
#endif
    FArrayBox q, qaux;

    int priv_nstep_fsp = -1;

    Real cflLoc = -1.0e+200;

    const int*  domain_lo = geom.Domain().loVect();
    const int*  domain_hi = geom.Domain().hiVect();

    for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi)
      {
	const Box& bx  = mfi.tilebox();
	const Box& qbx = amrex::grow(bx, NUM_GROW);

	const int* lo = bx.loVect();
	const int* hi = bx.hiVect();

	FArrayBox &statein  = Sborder[mfi];
	FArrayBox &stateout = S_new[mfi];

	FArrayBox &source_in  = sources_for_hydro[mfi];

	// the output of this will be stored in the correct stage MF
	FArrayBox &source_out = k_stage[mfi];
	FArrayBox &source_hydro_only = hydro_source[mfi];


	FArrayBox& vol = volume[mfi];

	q.resize(qbx, QVAR);
	qaux.resize(qbx, NQAUX);

	// convert the conservative state to the primitive variable state.
	// this fills both q and qaux.

	const int idx = mfi.tileIndex();

	ca_ctoprim(ARLIM_3D(qbx.loVect()), ARLIM_3D(qbx.hiVect()),
		   statein.dataPtr(), ARLIM_3D(statein.loVect()), ARLIM_3D(statein.hiVect()),
		   q.dataPtr(), ARLIM_3D(q.loVect()), ARLIM_3D(q.hiVect()),
		   qaux.dataPtr(), ARLIM_3D(qaux.loVect()), ARLIM_3D(qaux.hiVect()), &idx);


	// Allocate fabs for fluxes
	for (int i = 0; i < BL_SPACEDIM ; i++)  {
	  const Box& bxtmp = amrex::surroundingNodes(bx,i);
	  flux[i].resize(bxtmp,NUM_STATE);
	}

#if (BL_SPACEDIM <= 2)
	if (!Geometry::IsCartesian()) {
	  pradial.resize(amrex::surroundingNodes(bx,0),1);
	}
#endif

	ca_mol_single_stage
	  (&time,
	   lo, hi, domain_lo, domain_hi,
	   &(b_mol[mol_iteration]),
	   BL_TO_FORTRAN_3D(statein),
	   BL_TO_FORTRAN_3D(stateout),
	   BL_TO_FORTRAN_3D(q),
	   BL_TO_FORTRAN_3D(qaux),
	   BL_TO_FORTRAN_3D(source_in),
	   BL_TO_FORTRAN_3D(source_out),
	   BL_TO_FORTRAN_3D(source_hydro_only),
	   dx, &dt,
	   D_DECL(BL_TO_FORTRAN_3D(flux[0]),
		  BL_TO_FORTRAN_3D(flux[1]),
		  BL_TO_FORTRAN_3D(flux[2])),
#if (BL_SPACEDIM < 3)
	   BL_TO_FORTRAN_3D(pradial),
#endif
	   D_DECL(BL_TO_FORTRAN_3D(area[0][mfi]),
		  BL_TO_FORTRAN_3D(area[1][mfi]),
		  BL_TO_FORTRAN_3D(area[2][mfi])),
#if (BL_SPACEDIM < 3)
	   BL_TO_FORTRAN_3D(dLogArea[0][mfi]),
#endif
	   BL_TO_FORTRAN_3D(volume[mfi]),
	   &cflLoc, verbose);

	// Store the fluxes from this advance -- we weight them by the
	// integrator weight for this stage
	for (int i = 0; i < BL_SPACEDIM ; i++) {
	  (*fluxes    [i])[mfi].saxpy(b_mol[mol_iteration], flux[i],
				      mfi.nodaltilebox(i), mfi.nodaltilebox(i), 0, 0, NUM_STATE);
	}

#if (BL_SPACEDIM <= 2)
	if (!Geometry::IsCartesian()) {
	  P_radial[mfi].plus(pradial,mfi.nodaltilebox(0),0,0,1);
	}
#endif
      } // MFIter loop

#ifdef _OPENMP
#pragma omp critical (hydro_courno)
#endif
    {
      courno = std::max(courno,cflLoc);
    }
  }  // end of omp parallel region

  BL_PROFILE_VAR_STOP(CA_UMDRV);

  // Flush Fortran output

  if (verbose)
    flush_output();


  if (print_update_diagnostics)
    {

      bool local = true;
      Array<Real> hydro_update = evaluate_source_change(k_stage, dt, local);

#ifdef BL_LAZY
      Lazy::QueueReduction( [=] () mutable {
#endif
	  ParallelDescriptor::ReduceRealSum(hydro_update.dataPtr(), hydro_update.size(), ParallelDescriptor::IOProcessorNumber());

	  if (ParallelDescriptor::IOProcessor())
	    std::cout << std::endl << "  Contributions to the state from the hydro source:" << std::endl;

	  print_source_change(hydro_update);

#ifdef BL_LAZY
	});
#endif
    }


  if (courno > 1.0) {
    std::cout << "WARNING -- EFFECTIVE CFL AT THIS LEVEL " << level << " IS " << courno << '\n';
    if (hard_cfl_limit == 1)
      amrex::Abort("CFL is too high at this level -- go back to a checkpoint and restart with lower cfl number");
  }

}
