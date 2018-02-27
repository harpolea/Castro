/* problem-specific stuff goes here */

#ifdef __cplusplus
extern "C" {
#endif

  void problem_checkpoint(int * int_dir_name, int * len);

  void problem_restart(int * int_dir_name, int * len);

  void smooth_initial_data(BL_FORT_FAB_ARG_3D(dat), const int* lo, const int* hi, const amrex::Real* delta, const amrex::Real* xlo);

#ifdef __cplusplus
}
#endif
