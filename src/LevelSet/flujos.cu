#include<cuda.h>

__device__ double devDiv(
                         const double c1,
                         const double e1,
                         const double c2,
                         const double e2,
                         const double c3,
                         const double e3,
                         const double delta
                         )
{ 
   return delta*(c1*e1+c2*e2+c3*e3);
}
__global__ void  flux_continuity_CUDA(
                                     double* const e_flux,
                                     double* const f_flux,
                                     double* const g_flux,
                                     double* const um,
                                     unsigned const int Nx,
                                     unsigned const int Ny,
                                     unsigned const int Nz,
                                     double const c1 
                                     )
{

   unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x,
                idy = blockDim.y*blockIdx.y + threadIdx.y,
                idz = blockDim.z*blockIdx.z + threadIdx.z;

   unsigned int id = idx + idy*Nx + idz*Nx*Ny,
                Offset = Nx*Ny*Nz;

   e_flux[id] = -c1*um[id + 1*Offset];
   f_flux[id] = -c1*um[id + 2*Offset];
   g_flux[id] = -c1*um[id + 3*Offset];

   //Viscous terms
   e_flux[id + Offset] = 0.0;
   f_flux[id + Offset] = 0.0;
   g_flux[id + Offset] = 0.0;

   return;
}

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
                                     )
{

   unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x,
                idy = blockDim.y*blockIdx.y + threadIdx.y,
                idz = blockDim.z*blockIdx.z + threadIdx.z;

   unsigned int id = idx + idy*Nx + idz*Nx*Ny,
                Offset = Nx*Ny*Nz;

   double param0, param1, param2;

   e_flux[id] = -um[id + 1*Offset]*u[id] - press[id];
   f_flux[id] = -um[id + 2*Offset]*u[id];
   g_flux[id] = -um[id + 3*Offset]*u[id];


   //Viscous terms
   param0 = jbn[id           ]*ddvelp[id           ]   
          + jbn[id + 3*Offset]*ddvelp[id + 3*Offset]
          + jbn[id + 6*Offset]*ddvelp[id + 6*Offset];

   param1 = jbn[id + 1*Offset]*dcvel[id + 1*Offset]   
          + jbn[id + 4*Offset]*dcvel[id + 4*Offset]
          + jbn[id + 7*Offset]*dcvel[id + 7*Offset];

   param2 = jbn[id + 2*Offset]*dcvel[id + 2*Offset]   
          + jbn[id + 5*Offset]*dcvel[id + 5*Offset]
          + jbn[id + 8*Offset]*dcvel[id + 8*Offset];

   e_flux[id + Offset] = vis[id]*cdiv*(2.0*param0 - param1 - param2);

   param1 = jbn[id           ]*ddvelp[id + 1*Offset]   
          + jbn[id + 1*Offset]*ddvelp[id + 4*Offset]
          + jbn[id + 2*Offset]*ddvelp[id + 7*Offset];

   param2 = jbn[id + 1*Offset]*dcvel[id           ]   
          + jbn[id + 4*Offset]*dcvel[id + 3*Offset]
          + jbn[id + 7*Offset]*dcvel[id + 6*Offset];

   f_flux[id + Offset] = vis[id]*(param1 + param2);

   param1 = jbn[id           ]*dcvel[id + 2*Offset]   
          + jbn[id + 1*Offset]*dcvel[id + 5*Offset]
          + jbn[id + 2*Offset]*dcvel[id + 8*Offset];

   param2 = jbn[id + 2*Offset]*ddvelp[id           ]   
          + jbn[id + 5*Offset]*ddvelp[id + 3*Offset]
          + jbn[id + 8*Offset]*ddvelp[id + 6*Offset];

   g_flux[id + Offset] = vis[id]*(param1 + param2);

   return;
}


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
                                     )
{

   unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x,
                idy = blockDim.y*blockIdx.y + threadIdx.y,
                idz = blockDim.z*blockIdx.z + threadIdx.z;

   unsigned int id = idx + idy*Nx + idz*Nx*Ny,
                Offset = Nx*Ny*Nz;

   double param0, param1, param2;

   e_flux[id] = -um[id + 1*Offset]*u[id + 1*Offset];
   f_flux[id] = -um[id + 2*Offset]*u[id + 1*Offset] -  press[id];
   g_flux[id] = -um[id + 3*Offset]*u[id + 1*Offset];


   //Viscous terms

   param1 = jbn[id           ]*ddvelp[id + 1*Offset]   
          + jbn[id + 3*Offset]*ddvelp[id + 4*Offset]
          + jbn[id + 6*Offset]*ddvelp[id + 7*Offset];

   param2 = jbn[id + 1*Offset]*dcvel[id + 0*Offset]   
          + jbn[id + 4*Offset]*dcvel[id + 3*Offset]
          + jbn[id + 7*Offset]*dcvel[id + 6*Offset];

   e_flux[id + Offset] = vis[id]*(param1 + param2);


   param0 = jbn[id + 0*Offset]*dcvel[id + 0*Offset]   
          + jbn[id + 3*Offset]*dcvel[id + 3*Offset]
          + jbn[id + 6*Offset]*dcvel[id + 6*Offset];

   param1 = jbn[id + 1*Offset]*ddvelp[id + 1*Offset]   
          + jbn[id + 4*Offset]*ddvelp[id + 4*Offset]
          + jbn[id + 7*Offset]*ddvelp[id + 7*Offset];

   param2 = jbn[id + 2*Offset]*dcvel[id + 2*Offset]   
          + jbn[id + 5*Offset]*dcvel[id + 5*Offset]
          + jbn[id + 8*Offset]*dcvel[id + 8*Offset];

   f_flux[id + Offset] = vis[id]*cdiv*(2.0*param1 - param0 - param2);

   param1 = jbn[id + 1*Offset]*dcvel[id + 2*Offset]   
          + jbn[id + 4*Offset]*dcvel[id + 5*Offset]
          + jbn[id + 7*Offset]*dcvel[id + 8*Offset];

   param2 = jbn[id + 2*Offset]*ddvelp[id + 1*Offset]   
          + jbn[id + 5*Offset]*ddvelp[id + 4*Offset]
          + jbn[id + 8*Offset]*ddvelp[id + 7*Offset];

   g_flux[id + Offset] = vis[id]*(param1 + param2);

   return;
}


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
                                     )
{

   unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x,
                idy = blockDim.y*blockIdx.y + threadIdx.y,
                idz = blockDim.z*blockIdx.z + threadIdx.z;

   unsigned int id = idx + idy*Nx + idz*Nx*Ny,
                Offset = Nx*Ny*Nz;

   double param0, param1, param2;

   e_flux[id] = -um[id + 1*Offset]*u[id + 2*Offset];
   f_flux[id] = -um[id + 2*Offset]*u[id + 2*Offset];
   g_flux[id] = -um[id + 3*Offset]*u[id + 2*Offset] - press[id] 
                -froude*xmesh[id + 2*Offset];


   //Viscous terms
   param1 = jbn[id           ]*ddvelp[id + 2*Offset]   
          + jbn[id + 3*Offset]*ddvelp[id + 5*Offset]
          + jbn[id + 6*Offset]*ddvelp[id + 8*Offset];

   param2 = jbn[id + 2*Offset]*dcvel[id + 0*Offset]   
          + jbn[id + 5*Offset]*dcvel[id + 3*Offset]
          + jbn[id + 8*Offset]*dcvel[id + 6*Offset];

   e_flux[id + Offset] = vis[id]*(param1 + param2);


   param1 = jbn[id + 2*Offset]*dcvel[id + 1*Offset]   
          + jbn[id + 5*Offset]*dcvel[id + 4*Offset]
          + jbn[id + 8*Offset]*dcvel[id + 7*Offset];

   param2 = jbn[id + 1*Offset]*ddvelp[id + 2*Offset]   
          + jbn[id + 4*Offset]*ddvelp[id + 5*Offset]
          + jbn[id + 7*Offset]*ddvelp[id + 8*Offset];


   f_flux[id + Offset] = vis[id]*(param1 + param2);

   param0 = jbn[id + 0*Offset]*dcvel[id + 0*Offset]   
          + jbn[id + 3*Offset]*dcvel[id + 3*Offset]
          + jbn[id + 6*Offset]*dcvel[id + 6*Offset];

   param1 = jbn[id + 1*Offset]*dcvel[id + 1*Offset]   
          + jbn[id + 4*Offset]*dcvel[id + 4*Offset]
          + jbn[id + 7*Offset]*dcvel[id + 7*Offset];

   param2 = jbn[id + 2*Offset]*ddvelp[id + 2*Offset]   
          + jbn[id + 5*Offset]*ddvelp[id + 5*Offset]
          + jbn[id + 8*Offset]*ddvelp[id + 8*Offset];

   g_flux[id + Offset] = vis[id]*cdiv*(2.0*param2 - param1 - param0);

   return;

}


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
                                  )
{

   unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x,
                idy = blockDim.y*blockIdx.y + threadIdx.y,
                idz = blockDim.z*blockIdx.z + threadIdx.z;

   unsigned int id = idx + idy*Nx + idz*Nx*Ny,
                Offset = Nx*Ny*Nz;

   e_flux[id] = -um[id + 4*Offset]*u[id + 0*Offset];
   f_flux[id] = -um[id + 4*Offset]*u[id + 1*Offset];
   g_flux[id] = -um[id + 4*Offset]*u[id + 2*Offset];


   //Viscous terms
   e_flux[id + Offset] = vis5[id]*(
                         jbn[id           ]*dTemp[id + 0*Offset]   
                       + jbn[id + 3*Offset]*dTemp[id + 1*Offset]
                       + jbn[id + 6*Offset]*dTemp[id + 2*Offset]);

   f_flux[id + Offset] = vis5[id]*(
                         jbn[id + 1*Offset]*dTemp[id + 0*Offset]   
                       + jbn[id + 4*Offset]*dTemp[id + 1*Offset]
                       + jbn[id + 7*Offset]*dTemp[id + 2*Offset]);

   g_flux[id + Offset] = vis5[id]*(
                         jbn[id + 2*Offset]*dTemp[id + 0*Offset]   
                       + jbn[id + 5*Offset]*dTemp[id + 1*Offset]
                       + jbn[id + 8*Offset]*dTemp[id + 2*Offset]);


   return;

}

