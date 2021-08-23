#include<cuda.h>


__global__ void divDevXPlus(
                           double* const d_func,
                           const double* const func,
                           const unsigned int Nx, 	
                           const unsigned int Ny, 	
                           const unsigned int Nz,  	
                           const double deltax,
                           const double c1,
                           const double c2,
                           const double c3
                           );

__global__ void divDevXMin(
                           double* const d_func,
                           const double* const func,
                           const unsigned int Nx, 	
                           const unsigned int Ny, 	
                           const unsigned int Nz,  	
                           const double deltax,
                           const double c1,
                           const double c2,
                           const double c3
                           ) ;

__global__ void divDevYPlus(
                           double* const d_func,
                           const double* const func,
                           const unsigned int Nx, 	
                           const unsigned int Ny, 	
                           const unsigned int Nz,  	
                           const double deltaY,
                           const double c1,
                           const double c2,
                           const double c3
                           );

__global__ void divDevYMin(
                           double* const d_func,
                           const double* const func,
                           const unsigned int Nx, 	
                           const unsigned int Ny, 	
                           const unsigned int Nz,  	
                           const double deltaY,
                           const double c1,
                           const double c2,
                           const double c3
                           );

__global__ void divDevZPlus(
                           double* const d_func,
                           const double* const func,
                           const unsigned int Nx, 	
                           const unsigned int Ny, 	
                           const unsigned int Nz,  	
                           const double deltaZ,
                           const double c1,
                           const double c2,
                           const double c3
                           );

__global__ void divDevZMin(
                           double* const d_func,
                           const double* const func,
                           const unsigned int Nx, 	
                           const unsigned int Ny, 	
                           const unsigned int Nz,  	
                           const double deltaZ,
                           const double c1,
                           const double c2,
                           const double c3
                           );

__global__ void get_flux_e_CUDA(
                               double* const flux,
                               const double* const e_f,
                               const double* const f_f,
                               const double* const g_f,
                               const double* const jbn,
                               const unsigned int Nx, 
                               const unsigned int Ny, 
                               const unsigned int Nz 
                               );

__global__ void get_flux_f_CUDA(
                               double* const flux,
                               const double* const e_f,
                               const double* const f_f,
                               const double* const g_f,
                               const double* const jbn,
                               const unsigned int Nx, 
                               const unsigned int Ny, 
                               const unsigned int Nz 
                               );

__global__ void get_flux_g_CUDA(
                               double* const flux,
                               const double* const e_f,
                               const double* const f_f,
                               const double* const g_f,
                               const double* const jbn,
                               const unsigned int Nx, 
                               const unsigned int Ny, 
                               const unsigned int Nz 
                               );

__global__ void rs_divergence_CUDA(
                                   double* const rs,
                                   const double* const e_x,
                                   const double* const f_y,
                                   const double* const g_z,
                                   const double* const jbn,
                                   unsigned int Nx, 
                                   unsigned int Ny, 
                                   unsigned int Nz
                                   );
