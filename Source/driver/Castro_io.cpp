
#ifndef WIN32
#include <unistd.h>
#endif

#include <iomanip>
#include <iostream>
#include <string>
#include <ctime>

#include <AMReX_Utility.H>
#include "Castro.H"
#include "Castro_F.H"
#include "Castro_io.H"
#include <AMReX_ParmParse.H>

#ifdef _OPENMP
#include <omp.h>
#endif


#include "AMReX_buildInfo.H"

using std::string;
using namespace amrex;

// Castro maintains an internal checkpoint version numbering system.
// This allows us to maintain backwards compatibility with checkpoints
// generated by old versions of the code, so that new versions can
// restart from them. The version number is stored in the CastroHeader
// file inside a checkpoint. The following were the changes that were made
// in updating version numbers:
// 0: all checkpoints before we began the numbering system (July 27, 2015)
// 1: PhiGrav_Type was added to the checkpoint
// 2: Source_Type was added to the checkpoint
// 3: A ReactHeader file was generated and the maximum de/dt was stored there
// 4: Reactions_Type added to checkpoint; ReactHeader functionality deprecated
// 5: SDC_Source_Type and SDC_React_Type added to checkpoint

namespace
{
    int input_version = -1;
    int current_version = 5;
}

// I/O routines for Castro

void
Castro::restart (Amr&     papa,
                 istream& is,
                 bool     bReadSpecial)
{
    // Let's check Castro checkpoint version first;
    // trying to read from checkpoint; if nonexisting, set it to 0.
    if (input_version == -1) {
       	if (ParallelDescriptor::IOProcessor()) {
       	    std::ifstream CastroHeaderFile;
       	    std::string FullPathCastroHeaderFile = papa.theRestartFile();
       	    FullPathCastroHeaderFile += "/CastroHeader";
       	    CastroHeaderFile.open(FullPathCastroHeaderFile.c_str(), std::ios::in);
       	    if (CastroHeaderFile.good()) {
        		char foo[256];
        		// first line: Checkpoint version: ?
        		CastroHeaderFile.getline(foo, 256, ':');
        		CastroHeaderFile >> input_version;
           		CastroHeaderFile.close();
      	    } else {
       		       input_version = 0;
       	    }
       	}
      	ParallelDescriptor::Bcast(&input_version, 1, ParallelDescriptor::IOProcessorNumber());
    }

    BL_ASSERT(input_version >= 0);

    // also need to mod checkPoint function to store the new version in a text file

    AmrLevel::restart(papa,is,bReadSpecial);

    if (input_version == 0) { // old checkpoint without PhiGrav_Type
    }

    if (input_version < 3) { // old checkpoint without Source_Type
      state[Source_Type].restart(desc_lst[Source_Type], state[State_Type]);
    }

    // For versions < 2, we didn't store all three components
    // of the momenta in the checkpoint when doing 1D or 2D simulations.
    // So the state data that was read in will be a MultiFab with a
    // number of components that doesn't include the extra momenta,
    // which is incompatible with what we want. Our strategy is therefore
    // to create a new MultiFab with the right number of components, and
    // copy the data from the old MultiFab into the new one in the correct
    // slots. Then we'll swap pointers so that the new MultiFab is where
    // the new state data lives, and delete the old data as we no longer need it.

#if (BL_SPACEDIM < 3)

    if (input_version < 2) {

      int ns = desc_lst[State_Type].nComp();
      int ng = desc_lst[State_Type].nExtra();
      MultiFab* new_data = new MultiFab(grids,dmap,ns,ng);
      MultiFab& chk_data = get_state_data(State_Type).newData();

#if (BL_SPACEDIM == 1)

      // In 1D, we can copy everything below the y-momentum as normal,
      // and everything above the z-momentum as normal but shifted by
      // two components. The y- and z-momentum are zeroed out.

      for (int n = 0; n < ns; n++) {
    	if (n < Ymom)
    	  MultiFab::Copy(*new_data, chk_data, n,   n, 1, ng);
    	else if (n == Ymom || n == Zmom)
    	  new_data->setVal(0.0, n, 1, ng);
    	else
    	  MultiFab::Copy(*new_data, chk_data, n-2, n, 1, ng);
      }

#elif (BL_SPACEDIM == 2)

      // Strategy is the same in 2D but we only need to worry about
      // shifting by one component.

      for (int n = 0; n < ns; n++) {
    	if (n < Zmom)
    	  MultiFab::Copy(*new_data, chk_data, n,   n, 1, ng);
    	else if (n == Zmom)
    	  new_data->setVal(0.0, n, 1, ng);
    	else
    	  MultiFab::Copy(*new_data, chk_data, n-1, n, 1, ng);
      }

#endif

      // Now swap the pointers.

      get_state_data(State_Type).replaceNewData(new_data);

    }

#endif

    buildMetrics();

    initMFs();

    // get the elapsed CPU time to now;
    if (level == 0 && ParallelDescriptor::IOProcessor())
    {
      // get ellapsed CPU time
      std::ifstream CPUFile;
      std::string FullPathCPUFile = parent->theRestartFile();
      FullPathCPUFile += "/CPUtime";
      CPUFile.open(FullPathCPUFile.c_str(), std::ios::in);

      CPUFile >> previousCPUTimeUsed;
      CPUFile.close();

      std::cout << "read CPU time: " << previousCPUTimeUsed << "\n";

    }

    if (track_grid_losses && level == 0)
    {

      // get the current value of the diagnostic quantities
      std::ifstream DiagFile;
      std::string FullPathDiagFile = parent->theRestartFile();
      FullPathDiagFile += "/Diagnostics";
      DiagFile.open(FullPathDiagFile.c_str(), std::ios::in);

      for (int i = 0; i < n_lost; i++) {
    	DiagFile >> material_lost_through_boundary_cumulative[i];
    	material_lost_through_boundary_temp[i] = 0.0;
      }

      DiagFile.close();

    }

    if (level == 0)
    {
    	// get problem-specific stuff -- note all processors do this,
    	// eliminating the need for a broadcast
    	std::string dir = parent->theRestartFile();

    	char * dir_for_pass = new char[dir.size() + 1];
    	std::copy(dir.begin(), dir.end(), dir_for_pass);
    	dir_for_pass[dir.size()] = '\0';

    	int len = dir.size();

    	Array<int> int_dir_name(len);
    	for (int j = 0; j < len; j++)
    	  int_dir_name[j] = (int) dir_for_pass[j];

    	problem_restart(int_dir_name.dataPtr(), &len);

    	delete [] dir_for_pass;

    }

    const Real* dx  = geom.CellSize();

    if ( (grown_factor > 1) && (parent->maxLevel() < 1) )
    {
       std::cout << "grown_factor is " << grown_factor << std::endl;
       std::cout << "max_level is " << parent->maxLevel() << std::endl;
       amrex::Error("Must have max_level > 0 if doing special restart with grown_factor");
    }

    if (grown_factor > 1 && level == 0)
    {
       if (verbose && ParallelDescriptor::IOProcessor())
          std::cout << "Doing special restart with grown_factor " << grown_factor << std::endl;

       MultiFab& S_new = get_new_data(State_Type);
       Real cur_time   = state[State_Type].curTime();

       Box orig_domain;
       if (star_at_center == 0) {
          orig_domain = amrex::coarsen(geom.Domain(),grown_factor);
       } else if (star_at_center == 1) {

          Box domain(geom.Domain());
          int d,lo=0,hi=0;
          if (Geometry::IsRZ()) {
             if (grown_factor != 2)
                amrex::Abort("Must have grown_factor = 2");

             d = 0;
             int dlen =  domain.size()[d];
             lo = 0;
             hi = dlen/2;
             orig_domain.setSmall(d,lo);
             orig_domain.setBig(d,hi);

             d = 1;
             dlen =  domain.size()[d];
             lo =   dlen/4    ;
             hi = 3*dlen/4 - 1;
             orig_domain.setSmall(d,lo);
             orig_domain.setBig(d,hi);

          } else {
             for (int d = 0; d < BL_SPACEDIM; d++)
             {
                int dlen =  domain.size()[d];
                if (grown_factor == 2) {
                   lo =   dlen/4    ;
                   hi = 3*dlen/4 - 1;
                } else if (grown_factor == 3) {
                   lo =   (dlen)/3    ;
                   hi = 2*(dlen)/3 - 1;
                } else {
                   amrex::Abort("Must have grown_factor = 2 or 3");
                }
                orig_domain.setSmall(d,lo);
                orig_domain.setBig(d,hi);
             }
          }
       } else {
          if (ParallelDescriptor::IOProcessor())
             std::cout << "... invalid value of star_at_center: " << star_at_center << std::endl;
          amrex::Abort();
       }

       int ns = NUM_STATE;

       for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
       {
           const Box& bx      = mfi.validbox();
           const int* lo      = bx.loVect();
           const int* hi      = bx.hiVect();

           if (! orig_domain.contains(bx)) {

#ifdef DIMENSION_AGNOSTIC
              BL_FORT_PROC_CALL(CA_INITDATA,ca_initdata)
                (level, cur_time, ARLIM_3D(lo), ARLIM_3D(hi), ns,
        		 BL_TO_FORTRAN_3D(S_new[mfi]), ZFILL(dx),
        		 ZFILL(geom.ProbLo()), ZFILL(geom.ProbHi()));
#else
	      BL_FORT_PROC_CALL(CA_INITDATA,ca_initdata)
    		(level, cur_time, lo, hi, ns,
    		 BL_TO_FORTRAN(S_new[mfi]), dx,
    		 geom.ProbLo(), geom.ProbHi());
#endif
           }
       }
    }

    if (grown_factor > 1 && level == 1)
        getLevel(0).avgDown();

    // If we want, we can restart the checkpoint at a new time.

    if (reset_checkpoint_time > -1.e199) {

        if (!parent->RegridOnRestart())
            amrex::Error("reset_checkpoint_time only makes sense when amr.regrid_on_restart=1");

        const Real cur_time = get_state_data(State_Type).curTime();
        const Real prev_time = get_state_data(State_Type).prevTime();
        const Real dt = cur_time - prev_time;

        parent->setStartTime(reset_checkpoint_time);
        parent->setCumTime(reset_checkpoint_time);

        for (int n = 0; n < num_state_type; ++n) {
            StateData& state = get_state_data(n);
            state.setOldTimeLevel(reset_checkpoint_time-dt);
            state.setNewTimeLevel(reset_checkpoint_time   );
        }
    }

    if (reset_checkpoint_step > -1) {

        if (!parent->RegridOnRestart())
            amrex::Error("reset_checkpoint_step only makes sense when amr.regrid_on_restart=1");

        parent->setLevelSteps(level, reset_checkpoint_step);
        parent->setLevelCount(level, reset_checkpoint_step);

    }
}

void
Castro::set_state_in_checkpoint (Array<int>& state_in_checkpoint)
{
  for (int i=0; i<num_state_type; ++i)
    state_in_checkpoint[i] = 1;

  for (int i=0; i<num_state_type; ++i) {
    if (input_version < 3 && i == Source_Type) {
      // We are reading an old checkpoint with no Source_Type
      state_in_checkpoint[i] = 0;
    }
  }
}

void
Castro::checkPoint(const std::string& dir,
                   std::ostream&  os,
                   VisMF::How     how,
                   bool dump_old_default)
{
  AmrLevel::checkPoint(dir, os, how, dump_old);

  if (level == 0 && ParallelDescriptor::IOProcessor())
  {
	{
	    std::ofstream CastroHeaderFile;
	    std::string FullPathCastroHeaderFile = dir;
	    FullPathCastroHeaderFile += "/CastroHeader";
	    CastroHeaderFile.open(FullPathCastroHeaderFile.c_str(), std::ios::out);

	    CastroHeaderFile << "Checkpoint version: " << current_version << std::endl;
	    CastroHeaderFile.close();
	}

	{
	    // store elapsed CPU time
	    std::ofstream CPUFile;
	    std::string FullPathCPUFile = dir;
	    FullPathCPUFile += "/CPUtime";
	    CPUFile.open(FullPathCPUFile.c_str(), std::ios::out);

	    CPUFile << std::setprecision(15) << getCPUTime();
	    CPUFile.close();
	}

	if (track_grid_losses) {

	    // store diagnostic quantities
            std::ofstream DiagFile;
	    std::string FullPathDiagFile = dir;
	    FullPathDiagFile += "/Diagnostics";
	    DiagFile.open(FullPathDiagFile.c_str(), std::ios::out);

	    for (int i = 0; i < n_lost; i++)
	      DiagFile << std::setprecision(15) << material_lost_through_boundary_cumulative[i] << std::endl;

	    DiagFile.close();

	}

	{
	    // store any problem-specific stuff
	    char * dir_for_pass = new char[dir.size() + 1];
	    std::copy(dir.begin(), dir.end(), dir_for_pass);
	    dir_for_pass[dir.size()] = '\0';

	    int len = dir.size();

	    Array<int> int_dir_name(len);
	    for (int j = 0; j < len; j++)
		int_dir_name[j] = (int) dir_for_pass[j];

	    problem_checkpoint(int_dir_name.dataPtr(), &len);

	    delete [] dir_for_pass;
	}
    }

}

std::string
Castro::thePlotFileType () const
{
    //
    // Increment this whenever the writePlotFile() format changes.
    //
    static const std::string the_plot_file_type("HyperCLaw-V1.1");

    return the_plot_file_type;
}

void
Castro::setPlotVariables ()
{
  AmrLevel::setPlotVariables();

  // Don't add the Source_Type data to the plotfile, we only
  // want to store it in the checkpoints.

  for (int i = 0; i < desc_lst[Source_Type].nComp(); i++)
      parent->deleteStatePlotVar(desc_lst[Source_Type].name(i));

  ParmParse pp("castro");

  bool plot_X;

  if (pp.query("plot_X",plot_X))
  {
      if (plot_X)
      {
          //
    	  // Get the number of species from the network model.
              //
    	  ca_get_num_spec(&NumSpec);
              //
    	  // Get the species names from the network model.
              //
    	  for (int i = 0; i < NumSpec; i++)
          {
              int len = 20;
              Array<int> int_spec_names(len);
              //
              // This call return the actual length of each string in "len"
              //
              ca_get_spec_names(int_spec_names.dataPtr(),&i,&len);
              char* spec_name = new char[len+1];
              for (int j = 0; j < len; j++)
                  spec_name[j] = int_spec_names[j];
              spec_name[len] = '\0';
	          string spec_string = "X(";
              spec_string += spec_name;
              spec_string += ')';
	          parent->addDerivePlotVar(spec_string);
              delete [] spec_name;
    	  }
      }
  }
}


void
Castro::writeJobInfo (const std::string& dir)
{

  // job_info file with details about the run
  std::ofstream jobInfoFile;
  std::string FullPathJobInfoFile = dir;
  FullPathJobInfoFile += "/job_info";
  jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);

  std::string PrettyLine = std::string(78, '=') + "\n";
  std::string OtherLine = std::string(78, '-') + "\n";
  std::string SkipSpace = std::string(8, ' ');

  // job information
  jobInfoFile << PrettyLine;
  jobInfoFile << " Castro Job Information\n";
  jobInfoFile << PrettyLine;

  jobInfoFile << "job name: " << job_name << "\n\n";
  jobInfoFile << "inputs file: " << inputs_name << "\n\n";

  jobInfoFile << "number of MPI processes: " << ParallelDescriptor::NProcs() << "\n";
#ifdef _OPENMP
  jobInfoFile << "number of threads:       " << omp_get_max_threads() << "\n";
#endif

  jobInfoFile << "\n";
  jobInfoFile << "CPU time used since start of simulation (CPU-hours): " <<
    getCPUTime()/3600.0;

  jobInfoFile << "\n\n";

  // plotfile information
  jobInfoFile << PrettyLine;
  jobInfoFile << " Plotfile Information\n";
  jobInfoFile << PrettyLine;

  time_t now = time(0);

  // Convert now to tm struct for local timezone
  tm* localtm = localtime(&now);
  jobInfoFile   << "output data / time: " << asctime(localtm);

  char currentDir[FILENAME_MAX];
  if (getcwd(currentDir, FILENAME_MAX)) {
    jobInfoFile << "output dir:         " << currentDir << "\n";
  }

  jobInfoFile << "\n\n";


  // build information
  jobInfoFile << PrettyLine;
  jobInfoFile << " Build Information\n";
  jobInfoFile << PrettyLine;

  jobInfoFile << "build date:    " << buildInfoGetBuildDate() << "\n";
  jobInfoFile << "build machine: " << buildInfoGetBuildMachine() << "\n";
  jobInfoFile << "build dir:     " << buildInfoGetBuildDir() << "\n";
  jobInfoFile << "AMReX dir:     " << buildInfoGetAMReXDir() << "\n";

  jobInfoFile << "\n";

  jobInfoFile << "COMP:          " << buildInfoGetComp() << "\n";
  jobInfoFile << "COMP version:  " << buildInfoGetCompVersion() << "\n";

  jobInfoFile << "\n";

  jobInfoFile << "C++ compiler:  " << buildInfoGetCXXName() << "\n";
  jobInfoFile << "C++ flags:     " << buildInfoGetCXXFlags() << "\n";

  jobInfoFile << "\n";

  jobInfoFile << "Fortran comp:  " << buildInfoGetFName() << "\n";
  jobInfoFile << "Fortran flags: " << buildInfoGetFFlags() << "\n";

  jobInfoFile << "\n";

  jobInfoFile << "Link flags:    " << buildInfoGetLinkFlags() << "\n";
  jobInfoFile << "Libraries:     " << buildInfoGetLibraries() << "\n";

  jobInfoFile << "\n";

  for (int n = 1; n <= buildInfoGetNumModules(); n++) {
    jobInfoFile << buildInfoGetModuleName(n) << ": " << buildInfoGetModuleVal(n) << "\n";
  }

  jobInfoFile << "\n";

  const char* githash1 = buildInfoGetGitHash(1);
  const char* githash2 = buildInfoGetGitHash(2);
  const char* githash3 = buildInfoGetGitHash(3);
  if (strlen(githash1) > 0) {
    jobInfoFile << "Castro       git describe: " << githash1 << "\n";
  }
  if (strlen(githash2) > 0) {
    jobInfoFile << "AMReX        git describe: " << githash2 << "\n";
  }
  if (strlen(githash3) > 0) {
    jobInfoFile << "Microphysics git describe: " << githash3 << "\n";
  }

  const char* buildgithash = buildInfoGetBuildGitHash();
  const char* buildgitname = buildInfoGetBuildGitName();
  if (strlen(buildgithash) > 0){
    jobInfoFile << buildgitname << " git describe: " << buildgithash << "\n";
  }

  jobInfoFile << "\n\n";


  // grid information
  jobInfoFile << PrettyLine;
  jobInfoFile << " Grid Information\n";
  jobInfoFile << PrettyLine;

  int f_lev = parent->finestLevel();

  for (int i = 0; i <= f_lev; i++)
    {
      jobInfoFile << " level: " << i << "\n";
      jobInfoFile << "   number of boxes = " << parent->numGrids(i) << "\n";
      jobInfoFile << "   maximum zones   = ";
      for (int n = 0; n < BL_SPACEDIM; n++)
	{
	  jobInfoFile << parent->Geom(i).Domain().length(n) << " ";
	  //jobInfoFile << parent->Geom(i).ProbHi(n) << " ";
	}
      jobInfoFile << "\n\n";
    }

  jobInfoFile << " Boundary conditions\n";
  Array<int> lo_bc_out(BL_SPACEDIM), hi_bc_out(BL_SPACEDIM);
  ParmParse pp("castro");
  pp.getarr("lo_bc",lo_bc_out,0,BL_SPACEDIM);
  pp.getarr("hi_bc",hi_bc_out,0,BL_SPACEDIM);


  // these names correspond to the integer flags setup in the
  // Castro_setup.cpp
  const char* names_bc[] =
    { "interior", "inflow", "outflow",
      "symmetry", "slipwall", "noslipwall" };


  jobInfoFile << "   -x: " << names_bc[lo_bc_out[0]] << "\n";
  jobInfoFile << "   +x: " << names_bc[hi_bc_out[0]] << "\n";
  if (BL_SPACEDIM >= 2) {
    jobInfoFile << "   -y: " << names_bc[lo_bc_out[1]] << "\n";
    jobInfoFile << "   +y: " << names_bc[hi_bc_out[1]] << "\n";
  }
  if (BL_SPACEDIM == 3) {
    jobInfoFile << "   -z: " << names_bc[lo_bc_out[2]] << "\n";
    jobInfoFile << "   +z: " << names_bc[hi_bc_out[2]] << "\n";
  }

  jobInfoFile << "\n\n";


  // species info
  Real Aion = 0.0;
  Real Zion = 0.0;

  int mlen = 20;

  jobInfoFile << PrettyLine;
  jobInfoFile << " Species Information\n";
  jobInfoFile << PrettyLine;

  jobInfoFile <<
    std::setw(6) << "index" << SkipSpace <<
    std::setw(mlen+1) << "name" << SkipSpace <<
    std::setw(7) << "A" << SkipSpace <<
    std::setw(7) << "Z" << "\n";
  jobInfoFile << OtherLine;

  for (int i = 0; i < NumSpec; i++) {

    int len = mlen;
    Array<int> int_spec_names(len);
    //
    // This call return the actual length of each string in "len"
    //
    ca_get_spec_names(int_spec_names.dataPtr(),&i,&len);
    char* spec_name = new char[len+1];
    for (int j = 0; j < len; j++)
    spec_name[j] = int_spec_names[j];
    spec_name[len] = '\0';

    // get A and Z
    ca_get_spec_az(&i, &Aion, &Zion);

    jobInfoFile <<
    std::setw(6) << i << SkipSpace <<
    std::setw(mlen+1) << std::setfill(' ') << spec_name << SkipSpace <<
    std::setw(7) << Aion << SkipSpace <<
    std::setw(7) << Zion << "\n";
    delete [] spec_name;
  }
  jobInfoFile << "\n\n";


  // runtime parameters
  jobInfoFile << PrettyLine;
  jobInfoFile << " Inputs File Parameters\n";
  jobInfoFile << PrettyLine;

  ParmParse::dumpTable(jobInfoFile, true);

  jobInfoFile.close();

}

void
Castro::writePlotFile(const std::string& dir,
                      ostream& os,
                      VisMF::How how)
{
  plotFileOutput(dir, os, how, 0);
}


void
Castro::writeSmallPlotFile (const std::string& dir,
			    ostream&       os,
			    VisMF::How     how)
{
  plotFileOutput(dir, os, how, 1);
}


void
Castro::plotFileOutput(const std::string& dir,
                       ostream& os,
                       VisMF::How how,
                       const int is_small)
{

    int i, n;
    //
    // The list of indices of State to write to plotfile.
    // first component of pair is state_type,
    // second component of pair is component # within the state_type
    //
    std::vector<std::pair<int,int> > plot_var_map;
    for (int typ = 0; typ < desc_lst.size(); typ++)
        for (int comp = 0; comp < desc_lst[typ].nComp(); comp++)
            if (((parent->isStatePlotVar(desc_lst[typ].name(comp)) && is_small == 0) ||
                 (parent->isStateSmallPlotVar(desc_lst[typ].name(comp)) && is_small == 1)) &&
                desc_lst[typ].getType() == IndexType::TheCellType())
                plot_var_map.push_back(std::pair<int,int>(typ,comp));

    int num_derive = 0;
    std::list<std::string> derive_names;
    const std::list<DeriveRec>& dlist = derive_lst.dlist();

    for (std::list<DeriveRec>::const_iterator it = dlist.begin();
	 it != dlist.end();
	 ++it)
    {
        if ((parent->isDerivePlotVar(it->name()) && is_small == 0) ||
            (parent->isStateSmallPlotVar(it->name()) && is_small == 1))
        {
    		derive_names.push_back(it->name());
    		num_derive++;
	    }
    }

    int n_data_items = plot_var_map.size() + num_derive;

    Real cur_time = state[State_Type].curTime();

    if (level == 0 && ParallelDescriptor::IOProcessor())
    {
        //
        // The first thing we write out is the plotfile type.
        //
        os << thePlotFileType() << '\n';

        if (n_data_items == 0)
            amrex::Error("Must specify at least one valid data item to plot");

        os << n_data_items << '\n';

    	//
    	// Names of variables -- first state, then derived
    	//
    	for (i =0; i < plot_var_map.size(); i++)
        {
    	    int typ = plot_var_map[i].first;
    	    int comp = plot_var_map[i].second;
    	    os << desc_lst[typ].name(comp) << '\n';
        }

    	for ( std::list<std::string>::iterator it = derive_names.begin();
    	      it != derive_names.end(); ++it)
        {
    	    const DeriveRec* rec = derive_lst.get(*it);
                os << rec->variableName(0) << '\n';
        }

        os << BL_SPACEDIM << '\n';
        os << parent->cumTime() << '\n';
        int f_lev = parent->finestLevel();
        os << f_lev << '\n';
        for (i = 0; i < BL_SPACEDIM; i++)
            os << Geometry::ProbLo(i) << ' ';
        os << '\n';
        for (i = 0; i < BL_SPACEDIM; i++)
            os << Geometry::ProbHi(i) << ' ';
        os << '\n';
        for (i = 0; i < f_lev; i++)
            os << parent->refRatio(i)[0] << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
            os << parent->Geom(i).Domain() << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
            os << parent->levelSteps(i) << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
        {
            for (int k = 0; k < BL_SPACEDIM; k++)
                os << parent->Geom(i).CellSize()[k] << ' ';
            os << '\n';
        }
        os << (int) Geometry::Coord() << '\n';
        os << "0\n"; // Write bndry data.

	    writeJobInfo(dir);

    }
    // Build the directory to hold the MultiFab at this level.
    // The name is relative to the directory containing the Header file.
    //
    static const std::string BaseName = "/Cell";
    char buf[64];
    sprintf(buf, "Level_%d", level);
    std::string Level = buf;
    //
    // Now for the full pathname of that directory.
    //
    std::string FullPath = dir;
    if (!FullPath.empty() && FullPath[FullPath.size()-1] != '/')
        FullPath += '/';
    FullPath += Level;
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!amrex::UtilCreateDirectory(FullPath, 0755))
            amrex::CreateDirectoryFailed(FullPath);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor())
    {
        os << level << ' ' << grids.size() << ' ' << cur_time << '\n';
        os << parent->levelSteps(level) << '\n';

        for (i = 0; i < grids.size(); ++i)
        {
            RealBox gridloc = RealBox(grids[i],geom.CellSize(),geom.ProbLo());
            for (n = 0; n < BL_SPACEDIM; n++)
                os << gridloc.lo(n) << ' ' << gridloc.hi(n) << '\n';
        }
        //
        // The full relative pathname of the MultiFabs at this level.
        // The name is relative to the Header file containing this name.
        // It's the name that gets written into the Header.
        //
        if (n_data_items > 0)
        {
            std::string PathNameInHeader = Level;
            PathNameInHeader += BaseName;
            os << PathNameInHeader << '\n';
        }
    }
    //
    // We combine all of the multifabs -- state, derived, etc -- into one
    // multifab -- plotMF.
    // NOTE: we are assuming that each state variable has one component,
    // but a derived variable is allowed to have multiple components.
    int       cnt   = 0;
    const int nGrow = 0;
    MultiFab  plotMF(grids,dmap,n_data_items,nGrow);
    MultiFab* this_dat = 0;

    //
    // Cull data from state variables -- use no ghost cells.
    //
    for (i = 0; i < plot_var_map.size(); i++)
    {
    	int typ  = plot_var_map[i].first;
    	int comp = plot_var_map[i].second;
    	this_dat = &state[typ].newData();
    	MultiFab::Copy(plotMF,*this_dat,comp,cnt,1,nGrow);
    	cnt++;
    }
    //
    // Cull data from derived variables.
    //
    if (derive_names.size() > 0)
    {
    	for (std::list<std::string>::iterator it = derive_names.begin();
    	     it != derive_names.end(); ++it)
    	{
    	    auto derive_dat = derive(*it,cur_time,nGrow);
    	    MultiFab::Copy(plotMF,*derive_dat,0,cnt,1,nGrow);
    	    cnt++;
    	}
    }

    // int swe_to_comp_level;
    // ca_get_swe_to_comp_level(&swe_to_comp_level);
    // const Real* dx        = geom.CellSize();

    // convert to compressible

    // if ((level <= swe_to_comp_level) && (State_Type == 0)) {
    //     for (MFIter mfi(plotMF); mfi.isValid(); ++mfi)
    //     {
    //         const Box& bx = mfi.tilebox();
    //         // do some conversion stuff
    //         bool ignore_errors = false;
    //         RealBox gridloc = RealBox(grids[mfi.index()],geom.CellSize(),geom.ProbLo());
    //
    //         ca_swe_to_comp_self(BL_TO_FORTRAN_3D(plotMF[mfi]),
    //             ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()), ZFILL(dx), ZFILL(gridloc.lo()), &ignore_errors);
    //     }
    // }

    //
    // Use the Full pathname when naming the MultiFab.
    //
    std::string TheFullPath = FullPath;
    TheFullPath += BaseName;
    VisMF::Write(plotMF,TheFullPath,how,true);
}
