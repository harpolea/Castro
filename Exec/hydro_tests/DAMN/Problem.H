/* problem-specific Castro:: declarations go here */

void problem_post_timestep(amrex::MultiFab& S_new);

#ifdef DO_PROBLEM_POST_RESTART
void problem_post_restart();
#endif

#ifdef DO_PROBLEM_POST_INIT
void problem_post_init();
#endif

bool startup_error_removed = false;
const amrex::Real smooth_time = 0.1;
