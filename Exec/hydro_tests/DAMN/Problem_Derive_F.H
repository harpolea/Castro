/* problem-specific stuff goes here */

#ifndef _Problem_Derive_F_H_
#define _Problem_Derive_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif


void ca_derheight
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
     const int* level, const int* grid_no);

void ca_derprimv
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
     const int* level, const int* grid_no);

void ca_derprimrho
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
    const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
    const int* lo, const int* hi,
    const int* domain_lo, const int* domain_hi,
    const amrex::Real* delta, const amrex::Real* xlo,
    const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
    const int* level, const int* grid_no);

void ca_derW
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
    const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
    const int* lo, const int* hi,
    const int* domain_lo, const int* domain_hi,
    const amrex::Real* delta, const amrex::Real* xlo,
    const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
    const int* level, const int* grid_no);

void ca_dereint
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
    const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
    const int* lo, const int* hi,
    const int* domain_lo, const int* domain_hi,
    const amrex::Real* delta, const amrex::Real* xlo,
    const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
    const int* level, const int* grid_no);

void ca_derprimmom
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
    const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
    const int* lo, const int* hi,
    const int* domain_lo, const int* domain_hi,
    const amrex::Real* delta, const amrex::Real* xlo,
    const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
    const int* level, const int* grid_no);

#ifdef __cplusplus
}
#endif
#endif
