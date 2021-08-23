#include<cuda.h>
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

__global__ void extrapolKernel(
                               double* const rs,          //RS
                               const double* const extVal,//Var extrapol
                               const double* const phiS,  //Level Set F
                               const double* const jbn,   //Jacobian
                               const double* const d_Phi, //Phi Der
                               const double deltaX,
                               const double deltaY,
                               const double deltaZ,
                               const unsigned int Nx,
                               const unsigned int Ny,
                               const unsigned int Nz,
                               const int extFlag
                              )
{

   const int Offset = Nx*Ny*Nz;
   int id2;
   double so;
   double phiDeltaX, phiDeltaY, phiDeltaZ,
         d_ext_xu, d_ext_xd,
         d_ext_yu, d_ext_yd,
         d_ext_zu, d_ext_zd;

   int idx = blockIdx.x*blockDim.x + threadIdx.x,
                      idy = blockIdx.y*blockDim.y + threadIdx.y,
                      idz = blockIdx.z*blockDim.z + threadIdx.z;  


   //Offsets example (id_ip) EQ (i+1,j,k) 
   int id = Nx*Ny*idz + Nx*idy + idx,
                id_ip = Nx*Ny*idz + Nx*idy + idx + 1, 
                id_im = Nx*Ny*idz + Nx*idy + idx - 1, 
                id_jp = Nx*Ny*idz + Nx*(idy + 1) + idx, 
                id_jm = Nx*Ny*idz + Nx*(idy - 1) + idx, 
                id_kp = Nx*Ny*(idz + 1) + Nx*idy + idx, 
                id_km = Nx*Ny*(idz - 1) + Nx*idy + idx; 


   //Dealing with boundaries
   id2 = id;
   if(idx==0){id2 = id_ip; id_im = id;}
   if(idy==0){id2 = id_jp; id_jm = id;} 
   if(idz==0){id2 = id_kp; id_km = id;} 
   if(idx==Nx-1){id2 = id_im ; id_ip = id;} 
   if(idy==Ny-1){id2 = id_jm ; id_jp = id;} 
   if(idz==Nz-1){id2 = id_km ; id_kp = id;} 

  

// pick up the side to extrapol 
   if(extFlag>0){

      so = (double)(phiS[id]>0.0);       
   }
   else{
      so = -1.0*(double)(phiS[id]<=0.0);       
   }

  
   phiDeltaX = so*d_Phi[id           ];
   phiDeltaY = so*d_Phi[id + 1*Offset];
   phiDeltaZ = so*d_Phi[id + 2*Offset];


// Downwind derivatives of ext 
   d_ext_xd = deltaX*jbn[id           ]*(extVal[id2] - extVal[id_im]); 
   d_ext_yd = deltaY*jbn[id + 4*Offset]*(extVal[id2] - extVal[id_jm]); 
   d_ext_zd = deltaZ*jbn[id + 8*Offset]*(extVal[id2] - extVal[id_km]); 

// Upwind derivatives of ext 
   d_ext_xu = deltaX*jbn[id           ]*(extVal[id_ip] - extVal[id2]); 
   d_ext_yu = deltaY*jbn[id + 4*Offset]*(extVal[id_jp] - extVal[id2]); 
   d_ext_zu = deltaZ*jbn[id + 8*Offset]*(extVal[id_kp] - extVal[id2]); 
   

   double xMax = (double)(phiDeltaX > 0.0) 
               - (double)(phiDeltaX < 0.0);

   rs[id] = (0.5*(xMax + 1.0)*d_ext_xd 
          + 0.5*abs(xMax - 1.0)*d_ext_xu)*phiDeltaX;

   xMax = (double)(phiDeltaY > 0.0) 
        - (double)(phiDeltaY < 0.0);

   rs[id] += (0.5*(xMax + 1.0)*d_ext_yd 
          + 0.5*abs(xMax - 1.0)*d_ext_yu)*phiDeltaY;


   xMax = (double)(phiDeltaZ > 0.0) 
        - (double)(phiDeltaZ < 0.0);

   rs[id] += (0.5*(xMax + 1.0)*d_ext_zd 
          + 0.5*abs(xMax - 1.0)*d_ext_zu)*phiDeltaZ;

          
   return; 

}

__global__ void DevFirstOrder_LS(
                                double* const d_Phi,
                                const double* const phiS,
                                const double* const jbn,
                                const double deltaX,
                                const double deltaY,
                                const double deltaZ,
                                const unsigned int Nx, 
                                const unsigned int Ny, 
                                const unsigned int Nz
                                )
{

   const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x,
                      idy = blockIdx.y*blockDim.y + threadIdx.y,
                      idz = blockIdx.z*blockDim.z + threadIdx.z;  

   //Offsets example (id_ip) EQ (i+1,j,k) 
   unsigned int id = Nx*Ny*idz + Nx*idy + idx,
                id_ip = Nx*Ny*idz + Nx*idy + idx + 1, 
                id_im = Nx*Ny*idz + Nx*idy + idx - 1, 
                id_jp = Nx*Ny*idz + Nx*(idy + 1) + idx, 
                id_jm = Nx*Ny*idz + Nx*(idy - 1) + idx, 
                id_kp = Nx*Ny*(idz + 1) + Nx*idy + idx, 
                id_km = Nx*Ny*(idz - 1) + Nx*idy + idx; 
   
   double factor = 0.5;

   //Dealing with boundaries
   if(idx==0){id_im = id; factor = 1.0;}
   if(idy==0){id_jm = id; factor = 1.0;}
   if(idz==0){id_km = id; factor = 1.0;}
   if(idx==Nx-1){id_ip = id; factor = 1.0;}
   if(idy==Ny-1){id_jp = id; factor = 1.0;}
   if(idz==Nz-1){id_kp = id; factor = 1.0;}

   const unsigned int Offset = Nx*Ny*Nz;

   d_Phi[           id] = factor*deltaX*jbn[id           ]
                        * (phiS[id_ip] - phiS[id_im]);

   d_Phi[1*Offset + id] = factor*deltaY*jbn[id + 4*Offset]
                        * (phiS[id_jp] - phiS[id_jm]);

   d_Phi[2*Offset + id] = factor*deltaZ*jbn[id + 8*Offset]
                        * (phiS[id_kp] - phiS[id_km]);
	
   return;
}

__global__ void RunGK_FirstS(
                            double* d,
                            double* d0,
                            double  dt, 
                            double* rs,
                            const int Nx, const int Ny, const int Nz
                            )
{
   const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x,
                      idy = blockIdx.y*blockDim.y + threadIdx.y,
                      idz = blockIdx.z*blockDim.z + threadIdx.z;

   const unsigned int id = idx + idy*Nx + idz*Nx*Ny;
   
   d[id] = d0[id] - dt*rs[id];

   return;
}

__global__ void RunGK_SecondS(
                             double* d,
                             double* d0,
                             double* d1,
                             double  dt, 
                             double* rs,
                             const int Nx, const int Ny, const int Nz
                            )
{
   const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x,
                      idy = blockIdx.y*blockDim.y + threadIdx.y,
                      idz = blockIdx.z*blockDim.z + threadIdx.z;

   const unsigned int id = idx + idy*Nx + idz*Nx*Ny;
   
   d[id] = 0.75*d0[id] +0.25*( d1[id] - dt*rs[id]);

   return;
}

__global__ void RunGK_ThirdS(
                             double* d,
                             double* d0,
                             double* d1,
                             const double dt, 
                             double* rs,
                             const int Nx, const int Ny, const int Nz
                            )
{
   const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x,
                      idy = blockIdx.y*blockDim.y + threadIdx.y,
                      idz = blockIdx.z*blockDim.z + threadIdx.z;

   const unsigned int id = idx + idy*Nx + idz*Nx*Ny;
   
   d[id] = (d0[id] + 2.0*( d1[id] - dt*rs[id])) / 3.0 ;

   return;
}

__global__ void copyLSGas(
                          double* const value,
                          const double* const copyVal, 
                          const double* const phiS,
                          int Nx, int Ny, int Nz
                         )
{

   const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x,
                      idy = blockIdx.y*blockDim.y + threadIdx.y,
                      idz = blockIdx.z*blockDim.z + threadIdx.z;

   const unsigned int id = idx + idy*Nx + idz*Nx*Ny;
   
   value[id] = (phiS[id] > 0.0) ? copyVal[id] : value[id];
   
   return;
}


__global__ void copyLSLiquid(
                             double* const value,
                             const double* const copyVal, 
                             const double* const phiS,
                             int Nx, int Ny, int Nz,
                             int disp
                            )
{

   const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x,
                      idy = blockIdx.y*blockDim.y + threadIdx.y,
                      idz = blockIdx.z*blockDim.z + threadIdx.z;

   const unsigned int id = idx + idy*Nx + idz*Nx*Ny;

   
   const unsigned int offset = Nx*Ny*Nz;
   
   value[id + disp*offset] = (phiS[id] < 0.0) ? copyVal[id] :  
                                     value[id + disp*offset];
   
   return;
}
