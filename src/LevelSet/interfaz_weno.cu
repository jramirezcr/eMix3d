#include<cuda.h>
#include<stdio.h>
#include"extrapol.h"
#include"womegas.h"
#include"dimdef.h"


extern "C"{
void reiniz_CUDA(
                double* const phiS,
                const double* const jbn,
                const double* const deltaXYZ,
                const double deltaX,
                const double deltaY,
                const double deltaZ,
                int Nx,
                int Ny,
                int Nz,
                const double dtin
                )
{
   double *phiS_d, *phiS0_d,*jbn_d, *deltaXYZ_d,
          *d_phiP_d, *d_phiM_d, *rs_d, *d_Phi_d, *d_d;   

   unsigned int Offset = Nx*Ny*Nz;
   cudaMalloc((void**)&phiS_d,     sizeof(double)*Offset);
   cudaMalloc((void**)&phiS0_d,    sizeof(double)*Offset);
   cudaMalloc((void**)&jbn_d,      sizeof(double)*11*Offset);
   cudaMalloc((void**)&deltaXYZ_d, sizeof(double)*Offset);
   cudaMalloc((void**)&d_Phi_d,    sizeof(double)*3*Offset);
   cudaMalloc((void**)&d_phiP_d,   sizeof(double)*3*Offset);
   cudaMalloc((void**)&d_phiM_d,   sizeof(double)*3*Offset);
   cudaMalloc((void**)&rs_d,       sizeof(double)*Offset);
   cudaMalloc((void**)&d_d,        sizeof(double)*Offset);

   cudaMemcpy(phiS_d, phiS, sizeof(double)*Offset, 
              cudaMemcpyHostToDevice);
   cudaMemcpy(phiS0_d, phiS, sizeof(double)*Offset, 
              cudaMemcpyHostToDevice);
   cudaMemcpy(jbn_d, jbn, sizeof(double)*11*Offset, 
              cudaMemcpyHostToDevice);
   cudaMemcpy(deltaXYZ_d, deltaXYZ, sizeof(double)*Offset, 
              cudaMemcpyHostToDevice);

   dim3 DimBlock(BLOCKDMX,BLOCKDMY,BLOCKDMZ);   
   dim3 DimGrid(GRIDMX,GRIDMY,GRIDMZ);   

   for(int itera = 1 ; itera <= 20; itera++){
       //First Step
       
       Dev1thO_Downwind<<<DimGrid, DimBlock>>>( d_Phi_d, phiS_d,
                                    deltaX, deltaY, deltaZ,
                                    Nx, Ny, Nz);
       PhiDevPlusParameter<<<DimGrid, DimBlock>>>( d_phiP_d, d_Phi_d,
                                                   jbn_d, Nx, Ny, Nz);
     
       PhiDevMinusParameter<<<DimGrid, DimBlock>>>( d_phiM_d, d_Phi_d,
                                                    jbn_d, Nx, Ny, Nz);
     
       reini_RS_WENO<<<DimGrid, DimBlock>>>(rs_d, phiS_d, phiS0_d,deltaXYZ_d,
                                            d_phiP_d, d_phiM_d, Nx, Ny, Nz);

       RunGK_FirstS<<<DimGrid, DimBlock>>>(d_d, phiS_d, dtin, 
                                           rs_d, Nx, Ny, Nz);  

       //Second Step
     
       Dev1thO_Downwind<<<DimGrid, DimBlock>>>( d_Phi_d, d_d,
                                    deltaX, deltaY, deltaZ,
                                    Nx, Ny, Nz);
     
       PhiDevPlusParameter<<<DimGrid, DimBlock>>>( d_phiP_d, d_Phi_d,
                                                   jbn_d, Nx, Ny, Nz);
     
       PhiDevMinusParameter<<<DimGrid, DimBlock>>>( d_phiM_d, d_Phi_d,
                                                    jbn_d, Nx, Ny, Nz);
     
       reini_RS_WENO<<<DimGrid, DimBlock>>>(rs_d, d_d,phiS0_d, deltaXYZ_d,
                                            d_phiP_d, d_phiM_d, Nx, Ny, Nz );
     
       RunGK_SecondS<<<DimGrid, DimBlock>>>(d_d, phiS_d, d_d, dtin,
                                            rs_d, Nx, Ny, Nz);  
     
     
       //Third Step
       Dev1thO_Downwind<<<DimGrid, DimBlock>>>( d_Phi_d, d_d,
                                    deltaX, deltaY, deltaZ,
                                    Nx, Ny, Nz);
     
       PhiDevPlusParameter<<<DimGrid, DimBlock>>>( d_phiP_d, d_Phi_d,
                                                   jbn_d, Nx, Ny, Nz);
     
       PhiDevMinusParameter<<<DimGrid, DimBlock>>>( d_phiM_d, d_Phi_d,
                                                    jbn_d, Nx, Ny, Nz);
     
       reini_RS_WENO<<<DimGrid, DimBlock>>>(rs_d, d_d,phiS0_d, deltaXYZ_d,
                                          d_phiP_d, d_phiM_d, Nx, Ny, Nz );
     
       RunGK_ThirdS<<<DimGrid, DimBlock>>>(phiS_d, phiS_d, d_d, dtin,
                                            rs_d, Nx, Ny, Nz);  
   } 


   cudaMemcpy(phiS, phiS_d, sizeof(double)*Offset, 
              cudaMemcpyDeviceToHost);

 // check for error
  cudaError_t error = cudaGetLastError();
  if(error != cudaSuccess)
  {
    // print the CUDA error message and exit
    printf("CUDA error: %s\n", cudaGetErrorString(error));
    exit(-1);
  }


   cudaFree(phiS_d);
   cudaFree(phiS0_d);
   cudaFree(jbn_d);
   cudaFree(deltaXYZ_d);
   cudaFree(d_phiP_d);
   cudaFree(d_phiM_d);
   cudaFree(d_Phi_d);
   cudaFree(rs_d);
   cudaFree(d_d);
   return;
}
}

extern "C"{
void advect_CUDA(
                double* const phiS,
                const double* const velocity,
                const double* const jbn,
                double deltaX,
                double deltaY,
                double deltaZ,
                unsigned int Nx, unsigned int Ny, unsigned int Nz,
                const double dt
                )
{
   double *phiS_d, *velocity_d, *jbn_d,
          *rs_d, *d_d, *d_Phi_d, *d_phiP_d, *d_phiM_d;

   unsigned int Offset = Nx*Ny*Nz;

   cudaMalloc((void**)&phiS_d, sizeof(double)*Offset);
   cudaMalloc((void**)&velocity_d, sizeof(double)*3*Offset);
   cudaMalloc((void**)&jbn_d, sizeof(double)*11*Offset);
   cudaMalloc((void**)&rs_d, sizeof(double)*Offset);
   cudaMalloc((void**)&d_d, sizeof(double)*Offset);
   cudaMalloc((void**)&d_Phi_d, sizeof(double)*3*Offset);
   cudaMalloc((void**)&d_phiP_d, sizeof(double)*3*Offset);
   cudaMalloc((void**)&d_phiM_d, sizeof(double)*3*Offset);

   cudaMemcpy(phiS_d, phiS, sizeof(double)*Offset, 
              cudaMemcpyHostToDevice);
   cudaMemcpy(velocity_d, velocity, sizeof(double)*3*Offset, 
              cudaMemcpyHostToDevice);
   cudaMemcpy(jbn_d, jbn, sizeof(double)*11*Offset, 
              cudaMemcpyHostToDevice);


   dim3 DimBlock(BLOCKDMX,BLOCKDMY,BLOCKDMZ);   
   dim3 DimGrid(GRIDMX,GRIDMY,GRIDMZ);   

   Dev1thO_Downwind<<<DimGrid, DimBlock>>>( d_Phi_d, phiS_d,
                                deltaX, deltaY, deltaZ,
                                Nx, Ny, Nz);
   PhiDevPlusParameter<<<DimGrid, DimBlock>>>( d_phiP_d, d_Phi_d,
                                               jbn_d, Nx, Ny, Nz);
 
   PhiDevMinusParameter<<<DimGrid, DimBlock>>>( d_phiM_d, d_Phi_d,
                                                jbn_d, Nx, Ny, Nz);
 
   advect_RS_WENO<<<DimGrid, DimBlock>>>(rs_d, velocity_d,
                                        d_phiP_d, d_phiM_d, Nx, Ny, Nz);

   RunGK_FirstS<<<DimGrid, DimBlock>>>(d_d, phiS_d, dt, 
                                       rs_d, Nx, Ny, Nz);  

   //Second Step
 
   Dev1thO_Downwind<<<DimGrid, DimBlock>>>( d_Phi_d, d_d,
                                deltaX, deltaY, deltaZ,
                                Nx, Ny, Nz);
 
   PhiDevPlusParameter<<<DimGrid, DimBlock>>>( d_phiP_d, d_Phi_d,
                                               jbn_d, Nx, Ny, Nz);
 
   PhiDevMinusParameter<<<DimGrid, DimBlock>>>( d_phiM_d, d_Phi_d,
                                                jbn_d, Nx, Ny, Nz);
 

   advect_RS_WENO<<<DimGrid, DimBlock>>>(rs_d, velocity_d,
                                        d_phiP_d, d_phiM_d, Nx, Ny, Nz);
 
   RunGK_SecondS<<<DimGrid, DimBlock>>>(d_d, phiS_d, d_d, dt,
                                        rs_d, Nx, Ny, Nz);  
 
 
   //Third Step
   Dev1thO_Downwind<<<DimGrid, DimBlock>>>( d_Phi_d, d_d,
                                deltaX, deltaY, deltaZ,
                                Nx, Ny, Nz);
 
   PhiDevPlusParameter<<<DimGrid, DimBlock>>>( d_phiP_d, d_Phi_d,
                                               jbn_d, Nx, Ny, Nz);
 
   PhiDevMinusParameter<<<DimGrid, DimBlock>>>( d_phiM_d, d_Phi_d,
                                                jbn_d, Nx, Ny, Nz);
 

   advect_RS_WENO<<<DimGrid, DimBlock>>>(rs_d, velocity_d,
                                        d_phiP_d, d_phiM_d, Nx, Ny, Nz);
 
   RunGK_ThirdS<<<DimGrid, DimBlock>>>(phiS_d, phiS_d, d_d, dt,
                                           rs_d, Nx, Ny, Nz);  


   cudaMemcpy(phiS, phiS_d, sizeof(double)*Offset, 
              cudaMemcpyDeviceToHost);
   
   
   // check for error
   cudaError_t error = cudaGetLastError();
   if(error != cudaSuccess)
   {
     // print the CUDA error message and exit
      printf("CUDA error: %s\n", cudaGetErrorString(error));
      exit(-1);
   }
   
   cudaFree(phiS_d);
   cudaFree(jbn_d);
   cudaFree(velocity_d);
   cudaFree(d_phiP_d);
   cudaFree(d_phiM_d);
   cudaFree(d_Phi_d);
   cudaFree(rs_d);
   cudaFree(d_d);

   return;
}
}

