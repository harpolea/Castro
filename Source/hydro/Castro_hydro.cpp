#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

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

  //for (int n = 0; n < num_src; ++n)
//    MultiFab::Add(sources_for_hydro, *old_sources[n], 0, 0, NUM_STATE, 0);

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

    for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {
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
