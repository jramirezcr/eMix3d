#include<cuda.h>
#include<math.h>
#include"LevelSet/lsTools.h"
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#define PI 3.14159265359

__device__ double Phi_x_WENO(
                          double beta1,
                          double beta2,
                          double beta3,
                          double beta4,
                          double beta5
                          )
{
   
   double  s_b1, s_b2, s_b3,
          alpha_1, alpha_2, alpha_3,
          omega_1, omega_2, omega_3, result;

   s_b1 = (13.0/12.0)*(beta1 - 2.0*beta2 + beta3)
                     *(beta1 - 2.0*beta2 + beta3)
        + (0.25)*(beta1 - 4.0*beta2 + 3.0*beta3)
                *(beta1 - 4.0*beta2 + 3.0*beta3);

   s_b2 = (13.0/12.0)*(beta2 - 2.0*beta3 + beta4)
                     *(beta2 - 2.0*beta3 + beta4)
        + (0.25)*(beta2 - beta4)*(beta2 - beta4);

   s_b3 = (13.0/12.0)*(beta3 - 2.0*beta4 + beta5)
                     *(beta3 - 2.0*beta4 + beta5)
        + (0.25)*(3.0*beta3 - 4.0*beta4 + beta5)
                *(3.0*beta3 - 4.0*beta4 + beta5);


   alpha_1 = 0.1 /((s_b1 + 1.0e-6)*(s_b1 + 1.0e-6));
   alpha_2 = 0.6 /((s_b2 + 1.0e-6)*(s_b2 + 1.0e-6));
   alpha_3 = 0.3 /((s_b3 + 1.0e-6)*(s_b3 + 1.0e-6));

   omega_1 = alpha_1 / (alpha_1 + alpha_2 + alpha_3);
   omega_2 = alpha_2 / (alpha_1 + alpha_2 + alpha_3);
   omega_3 = alpha_3 / (alpha_1 + alpha_2 + alpha_3);
  
   result = ((omega_1*(2.0*beta1 - 7.0*beta2 + 11.0*beta3) 
       + omega_2*(-1.0*beta2 + 5.0*beta3 + 2.0*beta4)
       + omega_3*(2.0*beta3 + 5.0*beta4 - beta5))*(1.0/6.0));

   return result;
}


__global__ void Dev1thO_Downwind(
                                double* const       d_Phi,
                                const double* const phiS,
                                const double        deltaX,
                                const double        deltaY,
                                const double        deltaZ,
                                const unsigned int  Nx, 
                                const unsigned int  Ny, 
                                const unsigned int  Nz
                                )
{

   const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x,
                      idy = blockIdx.y*blockDim.y + threadIdx.y,
                      idz = blockIdx.z*blockDim.z + threadIdx.z;  

   //Offsets sample (id_ip) EQ (i+1,j,k) 
   unsigned int id = Nx*Ny*idz + Nx*idy + idx,
                id_im = Nx*Ny*idz + Nx*idy + idx - 1, 
                id_jm = Nx*Ny*idz + Nx*(idy - 1) + idx, 
                id_km = Nx*Ny*(idz - 1) + Nx*idy + idx; 
                   
   unsigned int ix = id, 
                iy = id, 
                iz = id;


   //Dealing with boundaries
   if(idx==0){id_im = id; ix = Nx*Ny*idz + Nx*idy + 1;} 
   if(idy==0){id_jm = id; iy = Nx*Ny*idz + Nx*1 + idx;} 
   if(idz==0){id_km = id; iz = Nx*Ny*1 + Nx*idy + idx;} 

   const unsigned int Offset = Nx*Ny*Nz;

   d_Phi[           id] = deltaX*(phiS[ix] - phiS[id_im]);

   d_Phi[1*Offset + id] = deltaY*(phiS[iy] - phiS[id_jm]);

   d_Phi[2*Offset + id] = deltaZ*(phiS[iz] - phiS[id_km]);
	
   return;

}

__global__ void PhiDevPlusParameter(
                                    double* const       phi_xyz,
                                    const double* const d_Phi,
                                    unsigned const int  Nx,
                                    unsigned const int  Ny,
                                    unsigned const int  Nz
                                    )
{
   unsigned const int Offset = Nx*Ny*Nz; 

   unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x,
                idy = blockDim.y*blockIdx.y + threadIdx.y, 
                idz = blockDim.z*blockIdx.z + threadIdx.z;

   unsigned int id = idx + idy*Nx + idz*Nx*Ny,
                id_im1 = (idx - 1) + idy*Nx + idz*Nx*Ny,
                id_ip1 = (idx + 1) + idy*Nx + idz*Nx*Ny,

                id_jm1 = idx + (idy - 1)*Nx + idz*Nx*Ny,
                id_jp1 = idx + (idy + 1)*Nx + idz*Nx*Ny,

                id_km1 = idx + idy*Nx + (idz - 1)*Nx*Ny,
                id_kp1 = idx + idy*Nx + (idz + 1)*Nx*Ny,

                id_im2 = (idx - 2) + idy*Nx + idz*Nx*Ny,
                id_ip2 = (idx + 2) + idy*Nx + idz*Nx*Ny,

                id_jm2 = idx + (idy - 2)*Nx + idz*Nx*Ny,
                id_jp2 = idx + (idy + 2)*Nx + idz*Nx*Ny,

                id_km2 = idx + idy*Nx + (idz - 2)*Nx*Ny,
                id_kp2 = idx + idy*Nx + (idz + 2)*Nx*Ny;


   //Dealing with boundaries

   if(idx == 0    ){id_im1 = id; id_im2 = id_ip1;} 
   if(idx == 1    ){id_im2 = id_im1;} 
   if(idx == Nx -1){id_ip1 = id; id_ip2 = id_im1;} 
   if(idx == Nx -2){id_ip2 = id_ip1;} 

   if(idy == 0    ){id_jm1 = id; id_jm2 = id_jp1;}
   if(idy == 1    ){id_jm2 = id_jm1;} 
   if(idy == Ny -1){id_jp1 = id; id_jp2 = id_jm1;} 
   if(idy == Ny -2){id_jp2 = id_jp1;} 

   if(idz == 0    ){id_km1 = id; id_km2 = id_kp1;} 
   if(idz == 1    ){id_km2 = id_km1;} 
   if(idz == Nz -1){id_kp1 = id; id_kp2 = id_jm1;} 
   if(idz == Nz -2){id_kp2 = id_kp1;} 

   double beta1, beta2, beta3, beta4, beta5;

  
   //Axis X

   beta1 = d_Phi[id_im2];
   beta2 = d_Phi[id_im1]; 
   beta3 = d_Phi[id];
   beta4 = d_Phi[id_ip1];
   beta5 = d_Phi[id_ip2];

   phi_xyz[id] = Phi_x_WENO(beta1, beta2, beta3, beta4, beta5);
   //Axis Y

   beta1 = d_Phi[id_jm2 + 1*Offset];
   beta2 = d_Phi[id_jm1 + 1*Offset]; 
   beta3 = d_Phi[id     + 1*Offset];
   beta4 = d_Phi[id_jp1 + 1*Offset];
   beta5 = d_Phi[id_jp2 + 1*Offset];


   phi_xyz[id + 1*Offset] = Phi_x_WENO(beta1, beta2, beta3, beta4, beta5);

 
   //Axis Z

   beta1 = d_Phi[id_km2 + 2*Offset];
   beta2 = d_Phi[id_km1 + 2*Offset]; 
   beta3 = d_Phi[id     + 2*Offset];
   beta4 = d_Phi[id_kp1 + 2*Offset];
   beta5 = d_Phi[id_kp2 + 2*Offset];

   phi_xyz[id + 2*Offset] = Phi_x_WENO(beta1, beta2, beta3, beta4, beta5);
   return;
}


__global__ void PhiDevMinusParameter(
                                     double* const       phi_xyz,
                                     const double* const d_Phi,
                                     unsigned const int  Nx,
                                     unsigned const int  Ny,
                                     unsigned const int  Nz
                                    )
{
   unsigned const int Offset = Nx*Ny*Nz; 

   unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x,
                idy = blockDim.y*blockIdx.y + threadIdx.y, 
                idz = blockDim.z*blockIdx.z + threadIdx.z;

   unsigned int id = idx + idy*Nx + idz*Nx*Ny,
                id_im1 = (idx - 1) + idy*Nx + idz*Nx*Ny,
                id_im2 = (idx - 2) + idy*Nx + idz*Nx*Ny,
                id_ip1 = (idx + 1) + idy*Nx + idz*Nx*Ny,

                id_jm1 = idx + (idy - 1)*Nx + idz*Nx*Ny,
                id_jm2 = idx + (idy - 2)*Nx + idz*Nx*Ny,
                id_jp1 = idx + (idy + 1)*Nx + idz*Nx*Ny,

                id_km1 = idx + idy*Nx + (idz - 1)*Nx*Ny,
                id_km2 = idx + idy*Nx + (idz - 2)*Nx*Ny,
                id_kp1 = idx + idy*Nx + (idz + 1)*Nx*Ny,

                id_ip2 = (idx + 2) + idy*Nx + idz*Nx*Ny,

                id_jp2 = idx + (idy + 2)*Nx + idz*Nx*Ny,

                id_kp2 = idx + idy*Nx + (idz + 2)*Nx*Ny,

                id_ip3 = (idx + 3) + idy*Nx + idz*Nx*Ny,

                id_jp3 = idx + (idy + 3)*Nx + idz*Nx*Ny,

                id_kp3 = idx + idy*Nx + (idz + 3)*Nx*Ny;

   //Dealing with boundaries

   if(idx == 0    ){id_im1 = id;} 
   if(idx == Nx -1){id_ip1 = id; id_ip2 = id_im1; id_ip3 = id_im2;} 
   if(idx == Nx -2){id_ip2 = id_ip1; id_ip3 = id;} 
   if(idx == Nx -3){id_ip3 = id_ip2;} 


   if(idy == 0    ){id_jm1 = id;} 
   if(idy == Ny -1){id_jp1 = id; id_jp2 = id_jm1; id_jp3 = id_jm2;} 
   if(idy == Ny -2){id_jp2 = id_jp1; id_jp3 = id;} 
   if(idy == Ny -3){id_jp3 = id_jp2;} 

   if(idz == 0    ){id_km1 = id;} 
   if(idz == Nz -1){id_kp1 = id; id_kp2 = id_km1;id_kp3 = id_km2;} 
   if(idz == Nz -2){id_kp2 = id_kp1; id_kp3 = id;} 
   if(idz == Nz -3){id_kp3 = id_kp2;} 


   double beta1, beta2, beta3, beta4, beta5;
  
   //Axis X

   beta1 = d_Phi[id_ip3];
   beta2 = d_Phi[id_ip2]; 
   beta3 = d_Phi[id_ip1];
   beta4 = d_Phi[id    ];
   beta5 = d_Phi[id_im1];


   phi_xyz[id           ] = Phi_x_WENO(beta1, beta2, beta3, beta4, beta5);
   //Axis Y

   beta1 = d_Phi[id_jp3 + 1*Offset];
   beta2 = d_Phi[id_jp2 + 1*Offset]; 
   beta3 = d_Phi[id_jp1 + 1*Offset];
   beta4 = d_Phi[id     + 1*Offset];
   beta5 = d_Phi[id_jm1 + 1*Offset];

   phi_xyz[id + 1*Offset] = Phi_x_WENO(beta1, beta2, beta3, beta4, beta5);

   //Axis Z


   beta1 = d_Phi[id_kp3 + 2*Offset];
   beta2 = d_Phi[id_kp2 + 2*Offset]; 
   beta3 = d_Phi[id_kp1 + 2*Offset];
   beta4 = d_Phi[id     + 2*Offset];
   beta5 = d_Phi[id_km1 + 2*Offset];

   phi_xyz[id + 2*Offset] = Phi_x_WENO(beta1, beta2, beta3, beta4, beta5);
   return;

}


__global__ void reini_RS_WENO(
                             double* const       rs,
                             const double* const phiS,                    
                             const double        deltaXYZ,
                             const double* const d_phiP,
                             const double* const d_phiM,
                             const double* const phiS0,                    
                             unsigned int        Nx,
                             unsigned int        Ny,
                             unsigned int        Nz 
                             )
{

   unsigned int Offset = Nx*Ny*Nz;

   unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x,
                idy = blockDim.y*blockIdx.y + threadIdx.y,
                idz = blockDim.z*blockIdx.z + threadIdx.z;

   unsigned int id = idx + Nx*idy + Nx*Ny*idz;
   double       so, rs_x, rs_y, rs_z, ta, grad_mod;
   double       phiMax, phiMin;

   ta = (double)(phiS0[id] > 0.0) - (double)(phiS0[id] < 0.0);

   //Getting gradient axis X
   phiMax = MAX(d_phiP[id   ], 0.0)*MAX(d_phiP[id   ], 0.0);  
   phiMin = MIN(d_phiM[id   ], 0.0)*MIN(d_phiM[id   ], 0.0);  

   rs_x   = 0.5*(ta + 1.0)*MAX(phiMax, phiMin);

   phiMax = MAX(d_phiM[id   ], 0.0)*MAX(d_phiM[id   ], 0.0);  
   phiMin = MIN(d_phiP[id   ], 0.0)*MIN(d_phiP[id   ], 0.0);  

   rs_x   += 0.5*abs(ta - 1.0)*MAX(phiMax, phiMin);

   //Getting gradient axis Y
   phiMax = MAX(d_phiP[id + 1*Offset], 0.0)
           *MAX(d_phiP[id + 1*Offset], 0.0);  

   phiMin = MIN(d_phiM[id + 1*Offset], 0.0)
           *MIN(d_phiM[id + 1*Offset], 0.0);  

   rs_y   = 0.5*(ta + 1.0)*MAX(phiMax, phiMin);

   phiMax = MAX(d_phiM[id + 1*Offset], 0.0)
           *MAX(d_phiM[id + 1*Offset], 0.0);  

   phiMin = MIN(d_phiP[id + 1*Offset], 0.0)
           *MIN(d_phiP[id + 1*Offset], 0.0);  

   rs_y   += 0.5*abs(ta - 1.0)*MAX(phiMax, phiMin);

   //Getting gradient axis Z
   phiMax = MAX(d_phiP[id + 2*Offset], 0.0)
           *MAX(d_phiP[id + 2*Offset], 0.0);  

   phiMin = MIN(d_phiM[id + 2*Offset], 0.0)
           *MIN(d_phiM[id + 2*Offset], 0.0);  

   rs_z   = 0.5*(ta + 1.0)*MAX(phiMax, phiMin);

   phiMax = MAX(d_phiM[id + 2*Offset], 0.0)
           *MAX(d_phiM[id + 2*Offset], 0.0);  

   phiMin = MIN(d_phiP[id + 2*Offset], 0.0)
           *MIN(d_phiP[id + 2*Offset], 0.0);  

   rs_z   += 0.5*abs(ta - 1.0)*MAX(phiMax, phiMin);

   grad_mod = sqrt(rs_x + rs_y + rs_z);

   so = phiS0[id] 
      / sqrt(phiS0[id]*phiS0[id] + deltaXYZ*deltaXYZ );

   rs[id] = 1.0*so*(grad_mod - 1.0);

   return;
}


__global__ void advect_RS_WENO(
                              double* const       rs,         //RHS 
                              const double* const velocity,   
                              const double* const d_phiP_d,
                              const double* const d_phiM_d,
                              unsigned int        Nx,
                              unsigned int        Ny,
                              unsigned int        Nz
                              )
{
   unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x,
                idy = blockDim.y*blockIdx.y + threadIdx.y,
                idz = blockDim.z*blockIdx.z + threadIdx.z;

   unsigned int id = idx + Nx*idy + Nx*Ny*idz,
                Offset = Nx*Ny*Nz;

   double rs_x, rs_y, rs_z;
   double grad_x, grad_y, grad_z;
   double rsign;

   rsign  = (double)(velocity[id] > 0.0) 
          - (double)(velocity[id] < 0.0); 
            
   rs_x   = 0.5*   (rsign + 1.0)*velocity[id]*d_phiP_d[id] 
          + 0.5*abs(rsign - 1.0)*velocity[id]*d_phiM_d[id];

   grad_x = 0.5*   (rsign + 1.0)*d_phiP_d[id] 
          + 0.5*abs(rsign - 1.0)*d_phiM_d[id];

   rsign  = (double)(velocity[id + 1*Offset] > 0.0) 
          - (double)(velocity[id + 1*Offset] < 0.0); 

   rs_y   = 0.5*(rsign + 1.0)*velocity[id + 1*Offset]
               *d_phiP_d[id + 1*Offset] 
          + 0.5*abs(rsign - 1.0)*velocity[id + 1*Offset]
               *d_phiM_d[id + 1*Offset];

   grad_y = 0.5*   (rsign + 1.0)*d_phiP_d[id + 1*Offset] 
          + 0.5*abs(rsign - 1.0)*d_phiM_d[id + 1*Offset];

   rsign  = (double)(velocity[id + 2*Offset] > 0.0) 
          - (double)(velocity[id + 2*Offset] < 0.0); 

   rs_z   = 0.5*(rsign + 1.0)*velocity[id + 2*Offset]
               *d_phiP_d[id + 2*Offset] 
          + 0.5*abs(rsign - 1.0)*velocity[id + 2*Offset]
               *d_phiM_d[id + 2*Offset];

   grad_z = 0.5*   (rsign + 1.0)*d_phiP_d[id + 2*Offset] 
          + 0.5*abs(rsign - 1.0)*d_phiM_d[id + 2*Offset];


   rs[id] = rs_x + rs_y + rs_z;
   
   return;
}


__global__
void enrightVelocityProfile(
                            double       *vel,    //Velocity Array
                            double       *xMesh,  //Mesh values 
                            double       *yMesh,
                            double       *zMesh,
                            const int     Nx,     //Mesh dimensions
                            const int     Ny,
                            const int     Nz,
                            const double  time,   //current time
                            const double  period 
                           )
{

   unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x,
                idy = blockDim.y*blockIdx.y + threadIdx.y,
                idz = blockDim.z*blockIdx.z + threadIdx.z;

   unsigned int id = idx + Nx*idy + Nx*Ny*idz,
                offset = Nx*Ny*Nz;
   
   vel[id           ] = 2.0*sin(PI*xMesh[id])*sin(PI*xMesh[id])
                           *sin(2.0*PI*yMesh[id])
                           *sin(2.0*PI*zMesh[id])*cos(PI*time/period);

   vel[id + 1*offset] =   -sin(PI*yMesh[id])*sin(PI*yMesh[id])
                          *sin(2.0*PI*xMesh[id])
                          *sin(2.0*PI*zMesh[id])*cos(PI*time/period);

   vel[id + 2*offset] =   -sin(PI*zMesh[id])*sin(PI*zMesh[id])
                          *sin(2.0*PI*yMesh[id])
                          *sin(2.0*PI*xMesh[id])*cos(PI*time/period);

      
}

__global__
void meshRegularStructured(
                           double       *xMesh,  //Mesh values 
                           double       *yMesh,
                           double       *zMesh,
                           double        deltaX,
                           double        deltaY,
                           double        deltaZ,
                           const int     Nx,     //Mesh dimensions
                           const int     Ny,
                           const int     Nz
                          )
{

   unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x,
                idy = blockDim.y*blockIdx.y + threadIdx.y,
                idz = blockDim.z*blockIdx.z + threadIdx.z;

   unsigned int id   = idx + Nx*idy + Nx*Ny*idz;

   xMesh[id] = (double)(idx - 5.0)*deltaX;   
   yMesh[id] = (double)(idy - 5.0)*deltaY;   
   zMesh[id] = (double)(idz - 5.0)*deltaZ;   

}

__global__
void cuGhostCellsMirror3dZ( 
                          double    *ghostArray,
                          const int  ncells,
                          const int  Nx,
                          const int  Ny,
                          const int  Nz,
                          double     direction
                          )
{

   int NxG = Nx + 2*ncells,
       NyG = Ny + 2*ncells;

   unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x,
                idy = blockDim.y*blockIdx.y + threadIdx.y;

   //Left boundary
   for(unsigned int idz = 0; idz < ncells; idz++){
      unsigned int id  = idx + NxG*idy + NxG*NyG*(ncells + idz);
      unsigned int idG = idx + NxG*idy + NxG*NyG*(ncells - idz - 1);

      ghostArray[idG] = ghostArray[id]*direction;    
   }

   //right boundary
   for(unsigned int idz = 0; idz < ncells; idz++){
      unsigned int id  = idx + NxG*idy + NxG*NyG*(Nz - idz - 1 + ncells);
      unsigned int idG = idx + NxG*idy + NxG*NyG*(Nz + idz + ncells);

      ghostArray[idG] = ghostArray[id]*direction;    
   }
}

__global__
void cuGhostCellsMirror3dY( 
                          double    *ghostArray,
                          const int  ncells,
                          const int  Nx,
                          const int  Ny,
                          const int  Nz,
                          double     direction
                          )
{

   int NxG = Nx + 2*ncells,
       NyG = Ny + 2*ncells,
       NzG = Nz + 2*ncells;

   unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x,
                idz = blockDim.y*blockIdx.y + threadIdx.y;

   //Left boundary
   for(unsigned int idy = 0; idy < ncells; idy++){
      unsigned int id  = idx + NxG*(ncells + idy) + NxG*NyG*idz;
      unsigned int idG = idx + NxG*(ncells - idy - 1) + NxG*NyG*idz;

      ghostArray[idG] = ghostArray[id]*direction;    
   }

   //right boundary
   for(unsigned int idy = 0; idy < ncells; idy++){
      unsigned int id  = idx + NxG*(Ny - idy - 1 + ncells) + NxG*NyG*idz;
      unsigned int idG = idx + NxG*(Ny + idy + ncells) + NxG*NyG*idz;

      ghostArray[idG] = ghostArray[id]*direction;    
   }
}

__global__
void cuGhostCellsMirror3dX( 
                          double    *ghostArray,
                          const int  ncells,
                          const int  Nx,
                          const int  Ny,
                          const int  Nz,
                          double     direction
                          )
{


   int NxG = Nx + 2*ncells,
       NyG = Ny + 2*ncells,
       NzG = Nz + 2*ncells;

   unsigned int idy = blockDim.x*blockIdx.x + threadIdx.x,
                idz = blockDim.y*blockIdx.y + threadIdx.y;

   //Left boundary
   for(unsigned int idx = 0; idx < ncells; idx++){
      unsigned int id  = (ncells + idx) + NxG*idy + NxG*NyG*idz;
      unsigned int idG = (ncells - idx - 1) + NxG*idy + NxG*NyG*idz;

      ghostArray[idG] = ghostArray[id]*direction;    
   }

   //right boundary
   for(unsigned int idx = 0; idx < ncells; idx++){
      unsigned int id  = (Nx - idx - 1 + ncells) + NxG*idy + NxG*NyG*idz;
      unsigned int idG = (Nx + idx + ncells) + NxG*idy + NxG*NyG*idz;

      ghostArray[idG] = ghostArray[id]*direction;    
   }
}

__global__ void PhiDevPlusParameterJB(
                                    double* const       phi_xyz,
                                    const double* const d_Phi,
                                    const double* const d_jbn,
                                    unsigned const int  Nx,
                                    unsigned const int  Ny,
                                    unsigned const int  Nz
                                    )
{
   unsigned const int Offset = Nx*Ny*Nz; 

   unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x,
                idy = blockDim.y*blockIdx.y + threadIdx.y, 
                idz = blockDim.z*blockIdx.z + threadIdx.z;

   unsigned int id = idx + idy*Nx + idz*Nx*Ny,
                id_im1 = (idx - 1) + idy*Nx + idz*Nx*Ny,
                id_ip1 = (idx + 1) + idy*Nx + idz*Nx*Ny,

                id_jm1 = idx + (idy - 1)*Nx + idz*Nx*Ny,
                id_jp1 = idx + (idy + 1)*Nx + idz*Nx*Ny,

                id_km1 = idx + idy*Nx + (idz - 1)*Nx*Ny,
                id_kp1 = idx + idy*Nx + (idz + 1)*Nx*Ny,

                id_im2 = (idx - 2) + idy*Nx + idz*Nx*Ny,
                id_ip2 = (idx + 2) + idy*Nx + idz*Nx*Ny,

                id_jm2 = idx + (idy - 2)*Nx + idz*Nx*Ny,
                id_jp2 = idx + (idy + 2)*Nx + idz*Nx*Ny,

                id_km2 = idx + idy*Nx + (idz - 2)*Nx*Ny,
                id_kp2 = idx + idy*Nx + (idz + 2)*Nx*Ny;


   //Dealing with boundaries

   if(idx == 0    ){id_im1 = id; id_im2 = id_ip1;} 
   if(idx == 1    ){id_im2 = id_im1;} 
   if(idx == Nx -1){id_ip1 = id; id_ip2 = id_im1;} 
   if(idx == Nx -2){id_ip2 = id_ip1;} 

   if(idy == 0    ){id_jm1 = id; id_jm2 = id_jp1;}
   if(idy == 1    ){id_jm2 = id_jm1;} 
   if(idy == Ny -1){id_jp1 = id; id_jp2 = id_jm1;} 
   if(idy == Ny -2){id_jp2 = id_jp1;} 

   if(idz == 0    ){id_km1 = id; id_km2 = id_kp1;} 
   if(idz == 1    ){id_km2 = id_km1;} 
   if(idz == Nz -1){id_kp1 = id; id_kp2 = id_jm1;} 
   if(idz == Nz -2){id_kp2 = id_kp1;} 

   double beta1, beta2, beta3, beta4, beta5;

  
   //Axis X

   beta1 = d_Phi[id_im2]*d_jbn[id];
   beta2 = d_Phi[id_im1]*d_jbn[id]; 
   beta3 = d_Phi[id    ]*d_jbn[id];
   beta4 = d_Phi[id_ip1]*d_jbn[id];
   beta5 = d_Phi[id_ip2]*d_jbn[id];

   phi_xyz[id] = Phi_x_WENO(beta1, beta2, beta3, beta4, beta5);
   //Axis Y

   beta1 = d_Phi[id_jm2 + 1*Offset]*d_jbn[id + 4*Offset];
   beta2 = d_Phi[id_jm1 + 1*Offset]*d_jbn[id + 4*Offset]; 
   beta3 = d_Phi[id     + 1*Offset]*d_jbn[id + 4*Offset];
   beta4 = d_Phi[id_jp1 + 1*Offset]*d_jbn[id + 4*Offset];
   beta5 = d_Phi[id_jp2 + 1*Offset]*d_jbn[id + 4*Offset];


   phi_xyz[id + 1*Offset] = Phi_x_WENO(beta1, beta2, beta3, beta4, beta5);

 
   //Axis Z

   beta1 = d_Phi[id_km2 + 2*Offset]*d_jbn[id + 8*Offset];
   beta2 = d_Phi[id_km1 + 2*Offset]*d_jbn[id + 8*Offset]; 
   beta3 = d_Phi[id     + 2*Offset]*d_jbn[id + 8*Offset];
   beta4 = d_Phi[id_kp1 + 2*Offset]*d_jbn[id + 8*Offset];
   beta5 = d_Phi[id_kp2 + 2*Offset]*d_jbn[id + 8*Offset];

   phi_xyz[id + 2*Offset] = Phi_x_WENO(beta1, beta2, beta3, beta4, beta5);
   return;
}


__global__ void PhiDevMinusParameterJB(
                                     double* const       phi_xyz,
                                     const double* const d_Phi,
                                     const double* const d_jbn,
                                     unsigned const int  Nx,
                                     unsigned const int  Ny,
                                     unsigned const int  Nz
                                    )
{
   unsigned const int Offset = Nx*Ny*Nz; 

   unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x,
                idy = blockDim.y*blockIdx.y + threadIdx.y, 
                idz = blockDim.z*blockIdx.z + threadIdx.z;

   unsigned int id = idx + idy*Nx + idz*Nx*Ny,
                id_im1 = (idx - 1) + idy*Nx + idz*Nx*Ny,
                id_im2 = (idx - 2) + idy*Nx + idz*Nx*Ny,
                id_ip1 = (idx + 1) + idy*Nx + idz*Nx*Ny,

                id_jm1 = idx + (idy - 1)*Nx + idz*Nx*Ny,
                id_jm2 = idx + (idy - 2)*Nx + idz*Nx*Ny,
                id_jp1 = idx + (idy + 1)*Nx + idz*Nx*Ny,

                id_km1 = idx + idy*Nx + (idz - 1)*Nx*Ny,
                id_km2 = idx + idy*Nx + (idz - 2)*Nx*Ny,
                id_kp1 = idx + idy*Nx + (idz + 1)*Nx*Ny,

                id_ip2 = (idx + 2) + idy*Nx + idz*Nx*Ny,

                id_jp2 = idx + (idy + 2)*Nx + idz*Nx*Ny,

                id_kp2 = idx + idy*Nx + (idz + 2)*Nx*Ny,

                id_ip3 = (idx + 3) + idy*Nx + idz*Nx*Ny,

                id_jp3 = idx + (idy + 3)*Nx + idz*Nx*Ny,

                id_kp3 = idx + idy*Nx + (idz + 3)*Nx*Ny;

   //Dealing with boundaries

   if(idx == 0    ){id_im1 = id;} 
   if(idx == Nx -1){id_ip1 = id; id_ip2 = id_im1; id_ip3 = id_im2;} 
   if(idx == Nx -2){id_ip2 = id_ip1; id_ip3 = id;} 
   if(idx == Nx -3){id_ip3 = id_ip2;} 


   if(idy == 0    ){id_jm1 = id;} 
   if(idy == Ny -1){id_jp1 = id; id_jp2 = id_jm1; id_jp3 = id_jm2;} 
   if(idy == Ny -2){id_jp2 = id_jp1; id_jp3 = id;} 
   if(idy == Ny -3){id_jp3 = id_jp2;} 

   if(idz == 0    ){id_km1 = id;} 
   if(idz == Nz -1){id_kp1 = id; id_kp2 = id_km1;id_kp3 = id_km2;} 
   if(idz == Nz -2){id_kp2 = id_kp1; id_kp3 = id;} 
   if(idz == Nz -3){id_kp3 = id_kp2;} 


   double beta1, beta2, beta3, beta4, beta5;
  
   //Axis X

   beta1 = d_Phi[id_ip3]*d_jbn[id];
   beta2 = d_Phi[id_ip2]*d_jbn[id]; 
   beta3 = d_Phi[id_ip1]*d_jbn[id];
   beta4 = d_Phi[id    ]*d_jbn[id];
   beta5 = d_Phi[id_im1]*d_jbn[id];


   phi_xyz[id           ] = Phi_x_WENO(beta1, beta2, beta3, beta4, beta5);
   //Axis Y

   beta1 = d_Phi[id_jp3 + 1*Offset]*d_jbn[id + 4*Offset];
   beta2 = d_Phi[id_jp2 + 1*Offset]*d_jbn[id + 4*Offset]; 
   beta3 = d_Phi[id_jp1 + 1*Offset]*d_jbn[id + 4*Offset];
   beta4 = d_Phi[id     + 1*Offset]*d_jbn[id + 4*Offset];
   beta5 = d_Phi[id_jm1 + 1*Offset]*d_jbn[id + 4*Offset];

   phi_xyz[id + 1*Offset] = Phi_x_WENO(beta1, beta2, beta3, beta4, beta5);

   //Axis Z


   beta1 = d_Phi[id_kp3 + 2*Offset]*d_jbn[id + 8*Offset];
   beta2 = d_Phi[id_kp2 + 2*Offset]*d_jbn[id + 8*Offset]; 
   beta3 = d_Phi[id_kp1 + 2*Offset]*d_jbn[id + 8*Offset];
   beta4 = d_Phi[id     + 2*Offset]*d_jbn[id + 8*Offset];
   beta5 = d_Phi[id_km1 + 2*Offset]*d_jbn[id + 8*Offset];

   phi_xyz[id + 2*Offset] = Phi_x_WENO(beta1, beta2, beta3, beta4, beta5);

   return;

}

__global__
void cuSwapToGhost( 
                   double    *ghostArray,
                   double    *valueArray,
                   const int  gcells,
                   const int  Nx,
                   const int  Ny,
                   const int  Nz
                  )
{

   int NxG = Nx + 2*gcells; 
   int NyG = Ny + 2*gcells; 

   unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x,
                idy = blockDim.y*blockIdx.y + threadIdx.y, 
                idz = blockDim.z*blockIdx.z + threadIdx.z;

   unsigned int id  = idx + idy*Nx + idz*Nx*Ny;
   unsigned int idG = idx + gcells + (idy + gcells)*NxG 
                    + (idz + gcells)*NxG*NyG;

   ghostArray[idG] = valueArray[id];

}

__global__
void cuSwapFromGhost( 
                   double    *valueArray,
                   double    *ghostArray,
                   const int  gcells,
                   const int  Nx,
                   const int  Ny,
                   const int  Nz
                  )
{

   int NxG = Nx + 2*gcells; 
   int NyG = Ny + 2*gcells; 

   unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x,
                idy = blockDim.y*blockIdx.y + threadIdx.y, 
                idz = blockDim.z*blockIdx.z + threadIdx.z;

   unsigned int id  = idx + idy*Nx + idz*Nx*Ny;
   unsigned int idG = idx + gcells + (idy + gcells)*NxG 
                    + (idz + gcells)*NxG*NyG;

   valueArray[id] = ghostArray[idG];

}

__global__ void reini_RS_WENOJB(
                             double* const       rs,
                             const double* const phiS,                    
                             const double* const deltaXYZ,
                             const double* const d_phiP,
                             const double* const d_phiM,
                             const double* const phiS0,                    
                             unsigned int        Nx,
                             unsigned int        Ny,
                             unsigned int        Nz 
                             )
{

   unsigned int Offset = Nx*Ny*Nz;

   unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x,
                idy = blockDim.y*blockIdx.y + threadIdx.y,
                idz = blockDim.z*blockIdx.z + threadIdx.z;

   unsigned int id = idx + Nx*idy + Nx*Ny*idz;
   double       so, rs_x, rs_y, rs_z, ta, grad_mod;
   double       phiMax, phiMin;

   ta = (double)(phiS0[id] > 0.0) - (double)(phiS0[id] < 0.0);

   //Getting gradient axis X
   phiMax = MAX(d_phiP[id   ], 0.0)*MAX(d_phiP[id   ], 0.0);  
   phiMin = MIN(d_phiM[id   ], 0.0)*MIN(d_phiM[id   ], 0.0);  

   rs_x   = 0.5*(ta + 1.0)*MAX(phiMax, phiMin);

   phiMax = MAX(d_phiM[id   ], 0.0)*MAX(d_phiM[id   ], 0.0);  
   phiMin = MIN(d_phiP[id   ], 0.0)*MIN(d_phiP[id   ], 0.0);  

   rs_x   += 0.5*abs(ta - 1.0)*MAX(phiMax, phiMin);

   //Getting gradient axis Y
   phiMax = MAX(d_phiP[id + 1*Offset], 0.0)
           *MAX(d_phiP[id + 1*Offset], 0.0);  

   phiMin = MIN(d_phiM[id + 1*Offset], 0.0)
           *MIN(d_phiM[id + 1*Offset], 0.0);  

   rs_y   = 0.5*(ta + 1.0)*MAX(phiMax, phiMin);

   phiMax = MAX(d_phiM[id + 1*Offset], 0.0)
           *MAX(d_phiM[id + 1*Offset], 0.0);  

   phiMin = MIN(d_phiP[id + 1*Offset], 0.0)
           *MIN(d_phiP[id + 1*Offset], 0.0);  

   rs_y   += 0.5*abs(ta - 1.0)*MAX(phiMax, phiMin);

   //Getting gradient axis Z
   phiMax = MAX(d_phiP[id + 2*Offset], 0.0)
           *MAX(d_phiP[id + 2*Offset], 0.0);  

   phiMin = MIN(d_phiM[id + 2*Offset], 0.0)
           *MIN(d_phiM[id + 2*Offset], 0.0);  

   rs_z   = 0.5*(ta + 1.0)*MAX(phiMax, phiMin);

   phiMax = MAX(d_phiM[id + 2*Offset], 0.0)
           *MAX(d_phiM[id + 2*Offset], 0.0);  

   phiMin = MIN(d_phiP[id + 2*Offset], 0.0)
           *MIN(d_phiP[id + 2*Offset], 0.0);  

   rs_z   += 0.5*abs(ta - 1.0)*MAX(phiMax, phiMin);

   grad_mod = sqrt(rs_x + rs_y + rs_z);

   so = phiS0[id] 
      / sqrt(phiS0[id]*phiS0[id] + deltaXYZ[id]*deltaXYZ[id] );

   rs[id] = 1.0*so*(grad_mod - 1.0);

   return;
}

