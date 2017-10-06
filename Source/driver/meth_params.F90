
! This file is automatically created by parse_castro_params.py.  To update
! or add runtime parameters, please edit _cpp_parameters and then run
! mk_params.sh

! This module stores the runtime parameters and integer names for
! indexing arrays.
!
! The Fortran-specific parameters are initialized in set_method_params(),
! and the ones that we are mirroring from C++ and obtaining through the
! ParmParse module are initialized in ca_set_castro_method_params().

module meth_params_module

  use bl_error_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  ! number of ghost cells for the hyperbolic solver
  integer, parameter     :: NHYP    = 4

  ! NTHERM: number of thermodynamic variables
  integer, save :: NTHERM, NVAR
  integer, save :: URHO, UMX, UMY, UMZ, UMR, UML, UMP, UEDEN, UEINT, UTEMP, UFA, UFS, UFX
  integer, save :: USHK

  ! QTHERM: number of primitive variables
  integer, save :: QTHERM, QVAR
  integer, save :: QRHO, QU, QV, QW, QPRES, QREINT, QTEMP, QGAME
  integer, save :: NQAUX, QGAMC, QC, QCSML, QDPDR, QDPDE
  integer, save :: QFA, QFS, QFX

  integer, save :: nadv

  ! NQ will be the total number of primitive variables, hydro + radiation
  integer, save :: NQ

  integer, save :: npassive
  integer, save, allocatable :: qpass_map(:), upass_map(:)

  ! These are used for the Godunov state
  ! Note that the velocity indices here are picked to be the same value
  ! as in the primitive variable array
  integer, save :: NGDNV, GDRHO, GDU, GDV, GDW, GDPRES, GDGAME

  integer         , save :: numpts_1d

  real(rt)        , save, allocatable :: outflow_data_old(:,:)
  real(rt)        , save, allocatable :: outflow_data_new(:,:)
  real(rt)        , save :: outflow_data_old_time
  real(rt)        , save :: outflow_data_new_time
  logical         , save :: outflow_data_allocated
  real(rt)        , save :: max_dist

  character(len=:), allocatable :: gravity_type

  ! these flags are for interpreting the EXT_DIR BCs
  integer, parameter :: EXT_UNDEFINED = -1
  integer, parameter :: EXT_HSE = 1
  integer, parameter :: EXT_INTERP = 2

  integer, save :: xl_ext, yl_ext, zl_ext, xr_ext, yr_ext, zr_ext

  ! Create versions of these variables on the GPU
  ! the device update is then done in Castro_nd.f90

  !$acc declare &
  !$acc create(NTHERM, NVAR) &
  !$acc create(URHO, UMX, UMY, UMZ, UMR, UML, UMP, UEDEN, UEINT, UTEMP, UFA, UFS,UFX) &
  !$acc create(USHK) &
  !$acc create(QTHERM, QVAR) &
  !$acc create(QRHO, QU, QV, QW, QPRES, QREINT, QTEMP) &
  !$acc create(QGAMC, QGAME) &
  !$acc create(NQ) &
  !$acc create(QFA, QFS, QFX) &
  !$acc create(xl_ext, yl_ext, zl_ext, xr_ext, yr_ext, zr_ext)

  ! Begin the declarations of the ParmParse parameters

  real(rt), save :: difmag
  real(rt), save :: small_dens
  real(rt), save :: small_temp
  real(rt), save :: small_pres
  real(rt), save :: small_ener
  integer         , save :: do_hydro
  integer         , save :: do_ctu
  integer         , save :: hybrid_hydro
  integer         , save :: ppm_type
  integer         , save :: ppm_trace_sources
  integer         , save :: ppm_temp_fix
  integer         , save :: ppm_predict_gammae
  integer         , save :: ppm_reference_eigenvectors
  integer         , save :: plm_iorder
  integer         , save :: hybrid_riemann
  integer         , save :: riemann_solver
  integer         , save :: cg_maxiter
  real(rt), save :: cg_tol
  integer         , save :: cg_blend
  integer         , save :: use_eos_in_riemann
  integer         , save :: use_flattening
  integer         , save :: transverse_use_eos
  integer         , save :: transverse_reset_density
  integer         , save :: transverse_reset_rhoe
  integer         , save :: dual_energy_update_E_from_e
  real(rt), save :: dual_energy_eta1
  real(rt), save :: dual_energy_eta2
  real(rt), save :: dual_energy_eta3
  integer         , save :: use_pslope
  integer         , save :: fix_mass_flux
  integer         , save :: limit_fluxes_on_small_dens
  integer         , save :: density_reset_method
  integer         , save :: allow_negative_energy
  integer         , save :: allow_small_energy
  integer         , save :: do_sponge
  integer         , save :: sponge_implicit
  integer         , save :: first_order_hydro
  character (len=:), allocatable, save :: xl_ext_bc_type
  character (len=:), allocatable, save :: xr_ext_bc_type
  character (len=:), allocatable, save :: yl_ext_bc_type
  character (len=:), allocatable, save :: yr_ext_bc_type
  character (len=:), allocatable, save :: zl_ext_bc_type
  character (len=:), allocatable, save :: zr_ext_bc_type
  integer         , save :: hse_zero_vels
  integer         , save :: hse_interp_temp
  integer         , save :: hse_reflect_vels
  real(rt), save :: cfl
  real(rt), save :: dtnuc_e
  real(rt), save :: dtnuc_X
  integer         , save :: dtnuc_mode
  real(rt), save :: dxnuc
  integer         , save :: do_react
  real(rt), save :: react_T_min
  real(rt), save :: react_T_max
  real(rt), save :: react_rho_min
  real(rt), save :: react_rho_max
  integer         , save :: disable_shock_burning
  real(rt), save :: diffuse_cutoff_density
  real(rt), save :: diffuse_cond_scale_fac
  integer         , save :: do_grav
  integer         , save :: grav_source_type
  integer         , save :: do_rotation
  real(rt), save :: rot_period
  real(rt), save :: rot_period_dot
  integer         , save :: rotation_include_centrifugal
  integer         , save :: rotation_include_coriolis
  integer         , save :: rotation_include_domegadt
  integer         , save :: state_in_rotating_frame
  integer         , save :: rot_source_type
  integer         , save :: implicit_rotation_update
  integer         , save :: rot_axis
  real(rt), save :: point_mass
  integer         , save :: point_mass_fix_solution
  integer         , save :: do_acc
  integer         , save :: grown_factor
  integer         , save :: track_grid_losses
  real(rt), save :: const_grav

  !$acc declare &
  !$acc create(difmag, small_dens, small_temp) &
  !$acc create(small_pres, small_ener, do_hydro) &
  !$acc create(do_ctu, hybrid_hydro, ppm_type) &
  !$acc create(ppm_trace_sources, ppm_temp_fix, ppm_predict_gammae) &
  !$acc create(ppm_reference_eigenvectors, plm_iorder, hybrid_riemann) &
  !$acc create(riemann_solver, cg_maxiter, cg_tol) &
  !$acc create(cg_blend, use_eos_in_riemann, use_flattening) &
  !$acc create(transverse_use_eos, transverse_reset_density, transverse_reset_rhoe) &
  !$acc create(dual_energy_update_E_from_e, dual_energy_eta1, dual_energy_eta2) &
  !$acc create(dual_energy_eta3, use_pslope, fix_mass_flux) &
  !$acc create(limit_fluxes_on_small_dens, density_reset_method, allow_negative_energy) &
  !$acc create(allow_small_energy, do_sponge, sponge_implicit) &
  !$acc create(first_order_hydro, hse_zero_vels, hse_interp_temp) &
  !$acc create(hse_reflect_vels, cfl, dtnuc_e) &
  !$acc create(dtnuc_X, dtnuc_mode, dxnuc) &
  !$acc create(do_react, react_T_min, react_T_max) &
  !$acc create(react_rho_min, react_rho_max, disable_shock_burning) &
  !$acc create(diffuse_cutoff_density, diffuse_cond_scale_fac, do_grav) &
  !$acc create(grav_source_type, do_rotation, rot_period) &
  !$acc create(rot_period_dot, rotation_include_centrifugal, rotation_include_coriolis) &
  !$acc create(rotation_include_domegadt, state_in_rotating_frame, rot_source_type) &
  !$acc create(implicit_rotation_update, rot_axis, point_mass) &
  !$acc create(point_mass_fix_solution, do_acc, grown_factor)

  ! End the declarations of the ParmParse parameters

  real(rt)        , save :: rot_vec(3)

contains

  subroutine ca_set_castro_method_params() bind(C, name="ca_set_castro_method_params")

    use amrex_parmparse_module, only: amrex_parmparse_build, amrex_parmparse_destroy, amrex_parmparse

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type (amrex_parmparse) :: pp


    const_grav = 0.0d0;

    call amrex_parmparse_destroy(pp)


    difmag = 0.1d0;
    small_dens = -1.d200;
    small_temp = -1.d200;
    small_pres = -1.d200;
    small_ener = -1.d200;
    do_hydro = -1;
    do_ctu = 1;
    hybrid_hydro = 0;
    ppm_type = 1;
    ppm_trace_sources = 1;
    ppm_temp_fix = 0;
    ppm_predict_gammae = 0;
    ppm_reference_eigenvectors = 0;
    plm_iorder = 2;
    hybrid_riemann = 0;
    riemann_solver = 0;
    cg_maxiter = 12;
    cg_tol = 1.0d-5;
    cg_blend = 2;
    use_eos_in_riemann = 0;
    use_flattening = 1;
    transverse_use_eos = 0;
    transverse_reset_density = 1;
    transverse_reset_rhoe = 0;
    dual_energy_update_E_from_e = 1;
    dual_energy_eta1 = 1.0d0;
    dual_energy_eta2 = 1.0d-4;
    dual_energy_eta3 = 1.0d0;
    use_pslope = 1;
    fix_mass_flux = 0;
    limit_fluxes_on_small_dens = 0;
    density_reset_method = 1;
    allow_negative_energy = 0;
    allow_small_energy = 1;
    do_sponge = 0;
    sponge_implicit = 1;
    first_order_hydro = 0;
    allocate(character(len=1)::xl_ext_bc_type)
    xl_ext_bc_type = "";
    allocate(character(len=1)::xr_ext_bc_type)
    xr_ext_bc_type = "";
    allocate(character(len=1)::yl_ext_bc_type)
    yl_ext_bc_type = "";
    allocate(character(len=1)::yr_ext_bc_type)
    yr_ext_bc_type = "";
    allocate(character(len=1)::zl_ext_bc_type)
    zl_ext_bc_type = "";
    allocate(character(len=1)::zr_ext_bc_type)
    zr_ext_bc_type = "";
    hse_zero_vels = 0;
    hse_interp_temp = 0;
    hse_reflect_vels = 0;
    cfl = 0.8d0;
    dtnuc_e = 1.d200;
    dtnuc_X = 1.d200;
    dtnuc_mode = 1;
    dxnuc = 1.d200;
    do_react = -1;
    react_T_min = 0.0d0;
    react_T_max = 1.d200;
    react_rho_min = 0.0d0;
    react_rho_max = 1.d200;
    disable_shock_burning = 0;
    diffuse_cutoff_density = -1.d200;
    diffuse_cond_scale_fac = 1.0d0;
    do_grav = -1;
    grav_source_type = 4;
    do_rotation = -1;
    rot_period = -1.d200;
    rot_period_dot = 0.0d0;
    rotation_include_centrifugal = 1;
    rotation_include_coriolis = 1;
    rotation_include_domegadt = 1;
    state_in_rotating_frame = 1;
    rot_source_type = 4;
    implicit_rotation_update = 1;
    rot_axis = 3;
    point_mass = 0.0d0;
    point_mass_fix_solution = 0;
    do_acc = -1;
    grown_factor = 1;
    track_grid_losses = 0;

    call amrex_parmparse_build(pp, "castro")
    call pp%query("difmag", difmag)
    call pp%query("small_dens", small_dens)
    call pp%query("small_temp", small_temp)
    call pp%query("small_pres", small_pres)
    call pp%query("small_ener", small_ener)
    call pp%query("do_hydro", do_hydro)
    call pp%query("do_ctu", do_ctu)
    call pp%query("hybrid_hydro", hybrid_hydro)
    call pp%query("ppm_type", ppm_type)
    call pp%query("ppm_trace_sources", ppm_trace_sources)
    call pp%query("ppm_temp_fix", ppm_temp_fix)
    call pp%query("ppm_predict_gammae", ppm_predict_gammae)
    call pp%query("ppm_reference_eigenvectors", ppm_reference_eigenvectors)
    call pp%query("plm_iorder", plm_iorder)
    call pp%query("hybrid_riemann", hybrid_riemann)
    call pp%query("riemann_solver", riemann_solver)
    call pp%query("cg_maxiter", cg_maxiter)
    call pp%query("cg_tol", cg_tol)
    call pp%query("cg_blend", cg_blend)
    call pp%query("use_eos_in_riemann", use_eos_in_riemann)
    call pp%query("use_flattening", use_flattening)
    call pp%query("transverse_use_eos", transverse_use_eos)
    call pp%query("transverse_reset_density", transverse_reset_density)
    call pp%query("transverse_reset_rhoe", transverse_reset_rhoe)
    call pp%query("dual_energy_update_E_from_e", dual_energy_update_E_from_e)
    call pp%query("dual_energy_eta1", dual_energy_eta1)
    call pp%query("dual_energy_eta2", dual_energy_eta2)
    call pp%query("dual_energy_eta3", dual_energy_eta3)
    call pp%query("use_pslope", use_pslope)
    call pp%query("fix_mass_flux", fix_mass_flux)
    call pp%query("limit_fluxes_on_small_dens", limit_fluxes_on_small_dens)
    call pp%query("density_reset_method", density_reset_method)
    call pp%query("allow_negative_energy", allow_negative_energy)
    call pp%query("allow_small_energy", allow_small_energy)
    call pp%query("do_sponge", do_sponge)
    call pp%query("sponge_implicit", sponge_implicit)
    call pp%query("first_order_hydro", first_order_hydro)
    call pp%query("xl_ext_bc_type", xl_ext_bc_type)
    call pp%query("xr_ext_bc_type", xr_ext_bc_type)
    call pp%query("yl_ext_bc_type", yl_ext_bc_type)
    call pp%query("yr_ext_bc_type", yr_ext_bc_type)
    call pp%query("zl_ext_bc_type", zl_ext_bc_type)
    call pp%query("zr_ext_bc_type", zr_ext_bc_type)
    call pp%query("hse_zero_vels", hse_zero_vels)
    call pp%query("hse_interp_temp", hse_interp_temp)
    call pp%query("hse_reflect_vels", hse_reflect_vels)
    call pp%query("cfl", cfl)
    call pp%query("dtnuc_e", dtnuc_e)
    call pp%query("dtnuc_X", dtnuc_X)
    call pp%query("dtnuc_mode", dtnuc_mode)
    call pp%query("dxnuc", dxnuc)
    call pp%query("do_react", do_react)
    call pp%query("react_T_min", react_T_min)
    call pp%query("react_T_max", react_T_max)
    call pp%query("react_rho_min", react_rho_min)
    call pp%query("react_rho_max", react_rho_max)
    call pp%query("disable_shock_burning", disable_shock_burning)
    call pp%query("do_grav", do_grav)
    call pp%query("grav_source_type", grav_source_type)
    call pp%query("do_rotation", do_rotation)
    call pp%query("do_acc", do_acc)
    call pp%query("grown_factor", grown_factor)
    call pp%query("track_grid_losses", track_grid_losses)
    call amrex_parmparse_destroy(pp)



    !$acc update &
    !$acc device(difmag, small_dens, small_temp) &
    !$acc device(small_pres, small_ener, do_hydro) &
    !$acc device(do_ctu, hybrid_hydro, ppm_type) &
    !$acc device(ppm_trace_sources, ppm_temp_fix, ppm_predict_gammae) &
    !$acc device(ppm_reference_eigenvectors, plm_iorder, hybrid_riemann) &
    !$acc device(riemann_solver, cg_maxiter, cg_tol) &
    !$acc device(cg_blend, use_eos_in_riemann, use_flattening) &
    !$acc device(transverse_use_eos, transverse_reset_density, transverse_reset_rhoe) &
    !$acc device(dual_energy_update_E_from_e, dual_energy_eta1, dual_energy_eta2) &
    !$acc device(dual_energy_eta3, use_pslope, fix_mass_flux) &
    !$acc device(limit_fluxes_on_small_dens, density_reset_method, allow_negative_energy) &
    !$acc device(allow_small_energy, do_sponge, sponge_implicit) &
    !$acc device(first_order_hydro, hse_zero_vels, hse_interp_temp) &
    !$acc device(hse_reflect_vels, cfl, dtnuc_e) &
    !$acc device(dtnuc_X, dtnuc_mode, dxnuc) &
    !$acc device(do_react, react_T_min, react_T_max) &
    !$acc device(react_rho_min, react_rho_max, disable_shock_burning) &
    !$acc device(diffuse_cutoff_density, diffuse_cond_scale_fac, do_grav) &
    !$acc device(grav_source_type, do_rotation, rot_period) &
    !$acc device(rot_period_dot, rotation_include_centrifugal, rotation_include_coriolis) &
    !$acc device(rotation_include_domegadt, state_in_rotating_frame, rot_source_type) &
    !$acc device(implicit_rotation_update, rot_axis, point_mass) &
    !$acc device(point_mass_fix_solution, do_acc, grown_factor)


    ! now set the external BC flags
    select case (xl_ext_bc_type)
    case ("hse", "HSE")
       xl_ext = EXT_HSE
    case ("interp", "INTERP")
       xl_ext = EXT_INTERP
    case default
       xl_ext = EXT_UNDEFINED
    end select

    select case (yl_ext_bc_type)
    case ("hse", "HSE")
       yl_ext = EXT_HSE
    case ("interp", "INTERP")
       yl_ext = EXT_INTERP
    case default
       yl_ext = EXT_UNDEFINED
    end select

    select case (zl_ext_bc_type)
    case ("hse", "HSE")
       zl_ext = EXT_HSE
    case ("interp", "INTERP")
       zl_ext = EXT_INTERP
    case default
       zl_ext = EXT_UNDEFINED
    end select

    select case (xr_ext_bc_type)
    case ("hse", "HSE")
       xr_ext = EXT_HSE
    case ("interp", "INTERP")
       xr_ext = EXT_INTERP
    case default
       xr_ext = EXT_UNDEFINED
    end select

    select case (yr_ext_bc_type)
    case ("hse", "HSE")
       yr_ext = EXT_HSE
    case ("interp", "INTERP")
       yr_ext = EXT_INTERP
    case default
       yr_ext = EXT_UNDEFINED
    end select

    select case (zr_ext_bc_type)
    case ("hse", "HSE")
       zr_ext = EXT_HSE
    case ("interp", "INTERP")
       zr_ext = EXT_INTERP
    case default
       zr_ext = EXT_UNDEFINED
    end select

    !$acc update device(xl_ext, yl_ext, zl_ext, xr_ext, yr_ext, zr_ext)


end subroutine ca_set_castro_method_params

subroutine f_set_castro_method_params()
    use network, only : nspec, naux
    implicit none

    integer :: dm, Density, Xmom, Ymom, Zmom, Eden, Eint, Temp, numadv, QLAST, FirstSpec, FirstAux, FirstAdv, cnt

  call ca_set_castro_method_params()

  numadv = 0

  NTHERM = 7
  NVAR = NTHERM + nspec + naux + numadv
  nadv = numadv

  Density = 0
  Xmom = 1
  Ymom = 2
  Zmom = 3
  Eden = 4
  Eint = 5
  Temp = 6

  cnt = Temp + 1

  if (nadv > 0) then
      FirstAdv = cnt
      cnt = cnt + nadv
  else
      FirstAdv = 0
  endif

  if (nspec > 0) then
      FirstSpec = cnt
      cnt = cnt + nspec
  else
      FirstSpec = 0
  endif

  if (naux > 0) then
      FirstAux = cnt
      cnt = cnt + naux
  else
      FirstAux = 0
  endif

  URHO  = Density   + 1
  UMX   = Xmom      + 1
  UMY   = Xmom      + 2
  UMZ   = Xmom      + 3
  UEDEN = Eden      + 1
  UEINT = Eint      + 1
  UTEMP = Temp      + 1

  if (numadv .ge. 1) then
     UFA   = FirstAdv  + 1
  else
     UFA = 1
  end if

  UFS   = FirstSpec + 1

  if (naux .ge. 1) then
     UFX = FirstAux  + 1
  else
     UFX = 1
  end if

  QTHERM = NTHERM + 1 ! the + 1 is for QGAME which is always defined in primitive mode

  QVAR = QTHERM + nspec + naux + numadv

  ! NQ will be the number of hydro + radiation variables in the primitive
  ! state.  Initialize it just for hydro here
  NQ = QVAR

  ! We use these to index into the state "Q"
  QRHO  = 1

  QU    = 2
  QV    = 3
  QW    = 4

  QGAME = 5

  QLAST   = QGAME

  QPRES   = QLAST + 1
  QREINT  = QLAST + 2

  QTEMP   = QTHERM ! = QLAST + 3

  if (numadv >= 1) then
     QFA = QTHERM + 1
     QFS = QFA + numadv

  else
     QFA = 1   ! density
     QFS = QTHERM + 1

  end if

  if (naux >= 1) then
     QFX = QFS + nspec

  else
     QFX = 1

  end if


end subroutine f_set_castro_method_params


  subroutine ca_finalize_meth_params() bind(C, name="ca_finalize_meth_params")
    implicit none

    deallocate(xl_ext_bc_type)
    deallocate(xr_ext_bc_type)
    deallocate(yl_ext_bc_type)
    deallocate(yr_ext_bc_type)
    deallocate(zl_ext_bc_type)
    deallocate(zr_ext_bc_type)

end subroutine ca_finalize_meth_params

subroutine f_finalize_meth_params()

  call ca_finalize_meth_params()

end subroutine f_finalize_meth_params


end module meth_params_module
