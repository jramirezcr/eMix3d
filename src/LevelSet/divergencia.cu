#include<cuda.h>
//only kernels

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
                           ) 
{

   double e1, e2, e3;
   unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x,
                idy = blockDim.y*blockIdx.y + threadIdx.y,
                idz = blockDim.z*blockIdx.z + threadIdx.z;


   unsigned int id = idx + Nx*idy + Nx*Ny*idz,
                id_ip1 = (idx + 1) + Nx*idy + Nx*Ny*idz,
                id_ip2 = (idx + 2) + Nx*idy + Nx*Ny*idz;

  
   e1 = func[id];
   e2 = func[id_ip1];
   e3 = func[id_ip2];

   if(idx == Nx-2){

      unsigned int id_ipNxm1 = (Nx - 1) + Nx*idy + Nx*Ny*idz,
                   id_ipNxm2 = (Nx - 2) + Nx*idy + Nx*Ny*idz,
                   id_ipNxm3 = (Nx - 3) + Nx*idy + Nx*Ny*idz,
                   id_ipNxm4 = (Nx - 4) + Nx*idy + Nx*Ny*idz;

      e1 = func[id_ipNxm2];
      e2 = func[id_ipNxm1];
      e3 = 4.0*func[id_ipNxm1] - 6.0*func[id_ipNxm2] + 4.0*func[id_ipNxm3] 
         - func[id_ipNxm4];
   }

   if(idx == Nx-1){
      e1 = 0.0;
      e2 = 0.0;
      e3 = 0.0;
   }


   d_func[id] = devDiv(c1, e1, c2, e2, c3, e3, deltax);

   return;
}

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
                           ) 
{

   double e1, e2, e3;
   unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x,
                idy = blockDim.y*blockIdx.y + threadIdx.y,
                idz = blockDim.z*blockIdx.z + threadIdx.z;


   unsigned int id = idx + Nx*idy + Nx*Ny*idz,
                id_im1 = (idx - 1) + Nx*idy + Nx*Ny*idz,
                id_im2 = (idx - 2) + Nx*idy + Nx*Ny*idz;
  
   e1 = func[id];
   e2 = func[id_im1];
   e3 = func[id_im2];

   if(idx == 1){

   unsigned int id_im0 = Nx*idy + Nx*Ny*idz,
                id_im0p1 = 1 + Nx*idy + Nx*Ny*idz,
                id_im0p2 = 2 + Nx*idy + Nx*Ny*idz,
                id_im0p3 = 3 + Nx*idy + Nx*Ny*idz;

      e1 = func[id_im0p1];
      e2 = func[id_im0];
      e3 = 4.0*func[id_im0] - 6.0*func[id_im0p1] + 4.0*func[id_im0p2] 
         - func[id_im0p3];
   }

   if(idx == 0){
      e1 = 0.0;
      e2 = 0.0;
      e3 = 0.0;
   }


   d_func[id] = devDiv(c1, e1, c2, e2, c3, e3, deltax);

   return;
}



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
                           ) 
{

   double e1, e2, e3;
   unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x,
                idy = blockDim.y*blockIdx.y + threadIdx.y,
                idz = blockDim.z*blockIdx.z + threadIdx.z;


   unsigned int id = idx + Nx*idy + Nx*Ny*idz,
                id_jp1 = idx + Nx*(idy + 1) + Nx*Ny*idz,
                id_jp2 = idx + Nx*(idy + 2) + Nx*Ny*idz;

  
   e1 = func[id];
   e2 = func[id_jp1];
   e3 = func[id_jp2];

   if(idy == Ny-2){

      unsigned int id_jpNym1 = idx + Nx*(Ny -1) + Nx*Ny*idz,
                   id_jpNym2 = idx + Nx*(Ny -2) + Nx*Ny*idz,
                   id_jpNym3 = idx + Nx*(Ny -3) + Nx*Ny*idz,
                   id_jpNym4 = idx + Nx*(Ny -4) + Nx*Ny*idz;

      e1 = func[id_jpNym2];
      e2 = func[id_jpNym1];
      e3 = 4.0*func[id_jpNym1] - 6.0*func[id_jpNym2] + 4.0*func[id_jpNym3] 
         - func[id_jpNym4];
   }

   if(idy == Ny - 1){
      e1 = 0.0;
      e2 = 0.0;
      e3 = 0.0;
   }


   d_func[id] = devDiv(c1, e1, c2, e2, c3, e3, deltaY);

   return;


}

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
                           ) 
{

   double e1, e2, e3;
   unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x,
                idy = blockDim.y*blockIdx.y + threadIdx.y,
                idz = blockDim.z*blockIdx.z + threadIdx.z;


   unsigned int id = idx + Nx*idy + Nx*Ny*idz,
                id_jm1 = idx + Nx*(idy - 1) + Nx*Ny*idz,
                id_jm2 = idx + Nx*(idy - 2) + Nx*Ny*idz;
  
   e1 = func[id];
   e2 = func[id_jm1];
   e3 = func[id_jm2];

   if(idy == 1){

   unsigned int id_jm0 = idx + Nx*Ny*idz,
                id_jm0p1 = idx + Nx*1 + Nx*Ny*idz,
                id_jm0p2 = idx + Nx*2 + Nx*Ny*idz,
                id_jm0p3 = idx + Nx*3 + Nx*Ny*idz;

      e1 = func[id_jm0p1];
      e2 = func[id_jm0];
      e3 = 4.0*func[id_jm0] - 6.0*func[id_jm0p1] + 4.0*func[id_jm0p2] 
         - func[id_jm0p3];
   }

   if(idy == 0){
      e1 = 0.0;
      e2 = 0.0;
      e3 = 0.0;
   }


   d_func[id] = devDiv(c1, e1, c2, e2, c3, e3, deltaY);

   return;
}


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
                           ) 
{

   double e1, e2, e3;
   unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x,
                idy = blockDim.y*blockIdx.y + threadIdx.y,
                idz = blockDim.z*blockIdx.z + threadIdx.z;


   unsigned int id = idx + Nx*idy + Nx*Ny*idz,
                id_kp1 = idx + Nx*idy + Nx*Ny*(idz + 1),
                id_kp2 = idx + Nx*idy + Nx*Ny*(idz + 2);

  
   e1 = func[id];
   e2 = func[id_kp1];
   e3 = func[id_kp2];

   if(idz == Nz - 2){

      unsigned int id_kpNzm1 = idx + Nx*idy + Nx*Ny*(Nz - 1),
                   id_kpNzm2 = idx + Nx*idy + Nx*Ny*(Nz - 2),
                   id_kpNzm3 = idx + Nx*idy + Nx*Ny*(Nz - 3),
                   id_kpNzm4 = idx + Nx*idy + Nx*Ny*(Nz - 4);

      e1 = func[id_kpNzm2];
      e2 = func[id_kpNzm1];
      e3 = 4.0*func[id_kpNzm1] - 6.0*func[id_kpNzm2] + 4.0*func[id_kpNzm3] 
         - func[id_kpNzm4];
   }

   if(idz == Nz - 1){
      e1 = 0.0;
      e2 = 0.0;
      e3 = 0.0;
   }


   d_func[id] = devDiv(c1, e1, c2, e2, c3, e3, deltaZ);

   return;
}

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
                           ) 
{

   double e1, e2, e3;
   unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x,
                idy = blockDim.y*blockIdx.y + threadIdx.y,
                idz = blockDim.z*blockIdx.z + threadIdx.z;


   unsigned int id = idx + Nx*idy + Nx*Ny*idz,
                id_km1 = idx + Nx*idy + Nx*Ny*(idz - 1),
                id_km2 = idx + Nx*idy + Nx*Ny*(idz - 2);
  
   e1 = func[id];
   e2 = func[id_km1];
   e3 = func[id_km2];

   if(idz == 1){

   unsigned int id_km0 = idx + Nx*Ny*idz,
                id_km0p1 = idx + idy + Nx*Ny*1,
                id_km0p2 = idx + idy + Nx*Ny*2,
                id_km0p3 = idx + idy + Nx*Ny*3;

      e1 = func[id_km0p1];
      e2 = func[id_km0];
      e3 = 4.0*func[id_km0] - 6.0*func[id_km0p1] + 4.0*func[id_km0p2] 
         - func[id_km0p3];
   }

   if(idz == 0){
      e1 = 0.0;
      e2 = 0.0;
      e3 = 0.0;
   }


   d_func[id] = devDiv(c1, e1, c2, e2, c3, e3, deltaZ);

   return;
}




__global__ void get_flux_e_CUDA(
                               double* const flux,
                               const double* const e_f,
                               const double* const f_f,
                               const double* const g_f,
                               const double* const jbn,
                               const unsigned int Nx, 
                               const unsigned int Ny, 
                               const unsigned int Nz 
                               )
{
   unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x,
                idy = blockDim.y*blockIdx.y + threadIdx.y,
                idz = blockDim.z*blockIdx.z + threadIdx.z;

   unsigned int id = idx + Nx*idy + Nx*Ny*idz, 
                Offset = Nx*Ny*Nz;

   flux[id] = jbn[id + 9*Offset]*(
                                 jbn[id + 0*Offset]*e_f[id] + 
                                 jbn[id + 1*Offset]*f_f[id] + 
                                 jbn[id + 2*Offset]*g_f[id] 
                                 );


   return;
}

__global__ void get_flux_f_CUDA(
                               double* const flux,
                               const double* const e_f,
                               const double* const f_f,
                               const double* const g_f,
                               const double* const jbn,
                               const unsigned int Nx, 
                               const unsigned int Ny, 
                               const unsigned int Nz 
                               )
{
   unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x,
                idy = blockDim.y*blockIdx.y + threadIdx.y,
                idz = blockDim.z*blockIdx.z + threadIdx.z;

   unsigned int id = idx + Nx*idy + Nx*Ny*idz, 
                Offset = Nx*Ny*Nz;

   flux[id] = jbn[id + 9*Offset]*(
                                 jbn[id + 3*Offset]*e_f[id] + 
                                 jbn[id + 4*Offset]*f_f[id] + 
                                 jbn[id + 5*Offset]*g_f[id] 
                                 );


   return;
}

__global__ void get_flux_g_CUDA(
                               double* const flux,
                               const double* const e_f,
                               const double* const f_f,
                               const double* const g_f,
                               const double* const jbn,
                               const unsigned int Nx, 
                               const unsigned int Ny, 
                               const unsigned int Nz 
                               )
{
   unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x,
                idy = blockDim.y*blockIdx.y + threadIdx.y,
                idz = blockDim.z*blockIdx.z + threadIdx.z;

   unsigned int id = idx + Nx*idy + Nx*Ny*idz, 
                Offset = Nx*Ny*Nz;

   flux[id] = jbn[id + 9*Offset]*(
                                 jbn[id + 6*Offset]*e_f[id] + 
                                 jbn[id + 7*Offset]*f_f[id] + 
                                 jbn[id + 8*Offset]*g_f[id]  
                                 );

   return;
}

__global__ void rs_divergence_CUDA(
                                   double* const rs,
                                   const double* const e_x,
                                   const double* const f_y,
                                   const double* const g_z,
                                   const double* const jbn,
                                   unsigned int Nx, 
                                   unsigned int Ny, 
                                   unsigned int Nz
                                   )
{
   unsigned int idx = threadIdx.x + blockIdx.x*blockDim.x,
                idy = threadIdx.y + blockIdx.y*blockDim.y,
                idz = threadIdx.z + blockIdx.z*blockDim.z;
 
   unsigned int id = idx + Nx*idy + Nx*Ny*idz,
                Offset = Nx*Ny*Nz;

   rs[id] = jbn[id + 10*Offset]*(
                                 e_x[id] + e_x[id + Offset] 
                               + f_y[id] + f_y[id + Offset] 
                               + g_z[id] + g_z[id + Offset] 
                                );

   return;
}
