#include<cuda.h>

__global__ void  flux_continuity_CUDA(
                                     double* const e_flux,
                                     double* const f_flux,
                                     double* const g_flux,
                                     double* const um,
                                     unsigned const int Nx,
                                     unsigned const int Ny,
                                     unsigned const int Nz,
                                     double const c1 
                                     );

__global__ void  flux_momentumX_CUDA(
                                     double* const e_flux,
                                     double* const f_flux,
                                     double* const g_flux,
                                     double* const u,
                                     double* const um,
                                     double* const press,
                                     double* const dcvel,
                                     double* const ddvelp,
                                     double* const vis,
                                     double* const jbn,
                                     unsigned const int Nx,
                                     unsigned const int Ny,
                                     unsigned const int Nz,
                                     double const cdiv 
                                     );
__global__ void  flux_momentumY_CUDA(
                                     double* const e_flux,
                                     double* const f_flux,
                                     double* const g_flux,
                                     double* const u,
                                     double* const um,
                                     double* const press,
                                     double* const dcvel,
                                     double* const ddvelp,
                                     double* const vis,
                                     double* const jbn,
                                     unsigned const int Nx,
                                     unsigned const int Ny,
                                     unsigned const int Nz,
                                     double const cdiv 
                                     );

__global__ void  flux_momentumZ_CUDA(
                                     double* const e_flux,
                                     double* const f_flux,
                                     double* const g_flux,
                                     double* const u,
                                     double* const um,
                                     double* const press,
                                     double* const dcvel,
                                     double* const ddvelp,
                                     double* const vis,
                                     double* const jbn,
                                     double* const xmesh,
                                     unsigned const int Nx,
                                     unsigned const int Ny,
                                     unsigned const int Nz,
                                     double const cdiv,
                                     double const froude
                                     );

__global__ void  flux_Energy_CUDA(
                                  double* const e_flux,
                                  double* const f_flux,
                                  double* const g_flux,
                                  double* const u,
                                  double* const um,
                                  double* const dTemp,
                                  double* const vis5,
                                  double* const jbn,
                                  unsigned const int Nx,
                                  unsigned const int Ny,
                                  unsigned const int Nz
                                  );

