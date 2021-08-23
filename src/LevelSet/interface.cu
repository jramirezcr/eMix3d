#include<cuda.h>
#include<stdio.h>
#include "LevelSet/extrapol.h"
#include "LevelSet/dimdef.h"


void exta(
         double* extVal_d,
         double* phiS_d,
         double* jbn_d,
         double* d_Phi_d,
         double* rs_d,
         double* d_d,
         double deltaX, double deltaY, double deltaZ,
         int Nx, int Ny, int Nz,
         double dtext,
         int Flag);

extern "C"{


void callCUDA(
             double* pressg,
             double* velocity,
             double* phiS,
             double* jbn,
             double deltaX,
             double deltaY,
             double deltaZ,
             unsigned int Nx,
             unsigned int Ny,
             unsigned int Nz,
             double dtext
             )
{

   unsigned int Offset = Nx*Ny*Nz;
   double *pressg_d, *velocity_d, *phiS_d, *jbn_d,
          *rs_d, *d_Phi_d, *extVal_d,
          *d_d;

   cudaMalloc((void**)&pressg_d,sizeof(double)*Offset);
   cudaMalloc((void**)&velocity_d,sizeof(double)*3*Offset);
   cudaMalloc((void**)&phiS_d,sizeof(double)*Offset);
   cudaMalloc((void**)&jbn_d,sizeof(double)*11*Offset);

   cudaMalloc((void**)&d_d,sizeof(double)*Offset);


   cudaMalloc((void**)&extVal_d,sizeof(double)*Offset);
   cudaMalloc((void**)&d_Phi_d,sizeof(double)*3*Offset);
   cudaMalloc((void**)&rs_d,sizeof(double)*Offset);

   cudaMemcpy(pressg_d, pressg, sizeof(double)*Offset, 
              cudaMemcpyHostToDevice );

   cudaMemcpy(velocity_d,velocity,sizeof(double)*3*Offset,
              cudaMemcpyHostToDevice );

   cudaMemcpy(phiS_d, phiS, sizeof(double)*Offset, 
              cudaMemcpyHostToDevice );

   cudaMemcpy(jbn_d, jbn, sizeof(double)*11*Offset, 
              cudaMemcpyHostToDevice );

   dim3 DimBlock(BLOCKDMX,BLOCKDMY,BLOCKDMZ);   
   dim3 DimGrid(GRIDMX,GRIDMY,GRIDMZ);   

   DevFirstOrder_LS<<<DimGrid, DimBlock>>>(
                                          d_Phi_d,
                                          phiS_d,
                                          jbn_d,
                                          deltaX,
                                          deltaY,
                                          deltaZ,
                                          Nx,
                                          Ny,
                                          Nz 
                                          );

  


// Extrapolating Velocity liquid variables

   exta(velocity_d, phiS_d, jbn_d, d_Phi_d, rs_d, d_d,
        deltaX, deltaY, deltaZ, Nx, Ny, Nz,
        dtext, -1);

       printf("  U Velocity Liquid \n");

   exta(&(velocity_d[1*Offset]), phiS_d, jbn_d, d_Phi_d, rs_d, d_d,
        deltaX, deltaY, deltaZ, Nx, Ny, Nz,
        dtext, -1);

       printf("  V Velocity Liquid \n");

   exta(&(velocity_d[2*Offset]), phiS_d, jbn_d, d_Phi_d, rs_d, d_d,
        deltaX, deltaY, deltaZ, Nx, Ny, Nz,
        dtext, -1);

       printf("  W Velocity Liquid \n");

// Extrapolating Gas Pressure Variable

   exta(pressg_d, phiS_d, jbn_d, d_Phi_d, rs_d, d_d,
        deltaX, deltaY, deltaZ, Nx, Ny, Nz,
        dtext, 1);

       printf("  Pressure Gas \n");

// Returning values from Device to Host

   cudaMemcpy(velocity,velocity_d,sizeof(double)*3*Offset,
              cudaMemcpyDeviceToHost );

   cudaMemcpy(pressg,pressg_d,sizeof(double)*Offset,
              cudaMemcpyDeviceToHost );
    
   cudaFree(pressg_d);
   cudaFree(velocity_d);
   cudaFree(jbn_d);
   cudaFree(phiS_d);
   cudaFree(d_Phi_d);
   cudaFree(d_d);
   cudaFree(extVal_d);
   cudaFree(rs_d);

   return;
}

}

extern "C"{

void normalLEvelSetCUDA(
                       double* nV,
                       double* phiS,
                       double* jbn,
                       double deltaX,
                       double deltaY,
                       double deltaZ,
                       unsigned int Nx,
                       unsigned int Ny,
                       unsigned int Nz
                       )
{

   unsigned int Offset = Nx*Ny*Nz;

   double *nV_d, *phiS_d, *jbn_d;

   cudaMalloc((void**)&nV_d,sizeof(double)*3*Offset);
   cudaMalloc((void**)&phiS_d,sizeof(double)*Offset);
   cudaMalloc((void**)&jbn_d,sizeof(double)*11*Offset);

   cudaMemcpy(phiS_d, phiS, sizeof(double)*Offset, 
              cudaMemcpyHostToDevice );

   cudaMemcpy(jbn_d, jbn, sizeof(double)*11*Offset, 
              cudaMemcpyHostToDevice );

   dim3 DimBlock(BLOCKDMX,BLOCKDMY,BLOCKDMZ);   
   dim3 DimGrid(GRIDMX,GRIDMY,GRIDMZ);   

   DevFirstOrder_LS<<<DimGrid, DimBlock>>>(
                                          nV_d,
                                          phiS_d,
                                          jbn_d,
                                          deltaX,
                                          deltaY,
                                          deltaZ,
                                          Nx,
                                          Ny,
                                          Nz 
                                          );

  

// Returning values from Device to Host

   cudaMemcpy(nV,nV_d,sizeof(double)*3*Offset,
              cudaMemcpyDeviceToHost );
    
   cudaFree(jbn_d);
   cudaFree(phiS_d);
   cudaFree(nV_d);

   return;
}

}


extern "C"{

void extrapolVarCUDA(
                    double* valToExt,
                    double* phiS,
                    double* d_Phi,
                    double* jbn,
                    double deltaX,
                    double deltaY,
                    double deltaZ,
                    unsigned int Nx,
                    unsigned int Ny,
                    unsigned int Nz,
                    double dtext
                    )
{

   unsigned int Offset = Nx*Ny*Nz;
   double *valToExt_d, *phiS_d, *jbn_d,
          *rs_d, *d_Phi_d, *extVal_d,
          *d_d;

   cudaMalloc((void**)&valToExt_d,sizeof(double)*Offset);
   cudaMalloc((void**)&phiS_d,sizeof(double)*Offset);
   cudaMalloc((void**)&d_Phi_d,sizeof(double)*3*Offset);
   cudaMalloc((void**)&jbn_d,sizeof(double)*11*Offset);

   cudaMalloc((void**)&extVal_d,sizeof(double)*Offset);
   cudaMalloc((void**)&rs_d,sizeof(double)*Offset);
   cudaMalloc((void**)&d_d,sizeof(double)*Offset);

   cudaMemcpy(valToExt_d,valToExt,sizeof(double)*Offset,
              cudaMemcpyHostToDevice );

   cudaMemcpy(phiS_d, phiS, sizeof(double)*Offset, 
              cudaMemcpyHostToDevice );

   cudaMemcpy(jbn_d, jbn, sizeof(double)*11*Offset, 
              cudaMemcpyHostToDevice );

   cudaMemcpy(d_Phi_d,d_Phi,sizeof(double)*3*Offset,
              cudaMemcpyHostToDevice );

   dim3 DimBlock(BLOCKDMX,BLOCKDMY,BLOCKDMZ);   
   dim3 DimGrid(GRIDMX,GRIDMY,GRIDMZ);   


// Extrapolating Velocity liquid variables

   exta(valToExt_d, phiS_d, jbn_d, d_Phi_d, rs_d, d_d,
        deltaX, deltaY, deltaZ, Nx, Ny, Nz,
        dtext, 1);

   printf("  Ext-Some Val \n");

// Returning values from Device to Host

   cudaMemcpy(valToExt,valToExt_d,sizeof(double)*Offset,
              cudaMemcpyDeviceToHost );
    
   cudaFree(valToExt_d);
   cudaFree(jbn_d);
   cudaFree(phiS_d);
   cudaFree(d_Phi_d);
   cudaFree(d_d);
   cudaFree(extVal_d);
   cudaFree(rs_d);

   return;
}

}



void exta(
         double* extVal_d,
         double* phiS_d,
         double* jbn_d,
         double* d_Phi_d,
         double* rs_d,
         double* d_d,
         double deltaX, double deltaY, double deltaZ,
         int Nx, int Ny, int Nz,
         double dtext,
         int Flag
         )
{


   dim3 DimBlock(BLOCKDMX,BLOCKDMY,BLOCKDMZ);   
   dim3 DimGrid(GRIDMX,GRIDMY,GRIDMZ);   
   printf("\n\n  Extrapolating on CUDA Device: \n");

   for(int itera = 1 ; itera <=10 ; itera++){

       extrapolKernel<<<DimGrid, DimBlock>>>(
                                          rs_d,         
                                          extVal_d, phiS_d, jbn_d, d_Phi_d,
                                            deltaX, deltaY, deltaZ,
                                            Nx, Ny, Nz, 
                                            Flag
                                            );  
     
       RunGK_FirstS<<<DimGrid, DimBlock>>>( d_d, extVal_d, 
                                           dtext, rs_d, Nx, Ny, Nz);
     
     
       extrapolKernel<<<DimGrid, DimBlock>>>(
                                            rs_d, d_d, 
                                            phiS_d, jbn_d, d_Phi_d,
                                            deltaX, deltaY, deltaZ,
                                            Nx, Ny, Nz, 
                                            Flag);
     
       RunGK_SecondS<<<DimGrid, DimBlock>>>( d_d, 
                                             extVal_d, d_d, 
                                             dtext, rs_d, 
                                             Nx, Ny, Nz);
     
     
       extrapolKernel<<<DimGrid, DimBlock>>>(
                                            rs_d, d_d, 
                                            phiS_d, jbn_d, d_Phi_d,
                                            deltaX, deltaY, deltaZ,
                                            Nx, Ny, Nz, 
                                            Flag);
     
       RunGK_ThirdS<<<DimGrid, DimBlock>>>(  extVal_d, 
                                             extVal_d, d_d, 
                                             dtext, rs_d, 
                                             Nx, Ny, Nz);
       
   }
 // check for error
  cudaError_t error = cudaGetLastError();
  if(error != cudaSuccess)
  {
    // print the CUDA error message and exit
    printf("CUDA error: %s\n", cudaGetErrorString(error));
    exit(-1);
  }
 
 }

 
void cuExtrapolation(
                    double* extVal_d,
                    double* phiS_d,
                    double* jbn_d,
                    double deltaX, double deltaY, double deltaZ,
                    int Nx, int Ny, int Nz,
                    double dtext,
                    int Flag
                    )
{

   double *d_dPhi;
   double *rs_d;
   double *d_d;

   int Offset = Nx*Ny*Nz;

   cudaMalloc((void**)&d_dPhi, 3*sizeof(double)*Offset);
   cudaMalloc((void**)&rs_d,     sizeof(double)*Offset);
   cudaMalloc((void**)&d_d,      sizeof(double)*Offset);

   int numGBX, numGBY,numGBZ;

   dim3 dimBlock(10,10,5);   

   numGBX = Nx / 10;
   numGBY = Ny / 10;
   numGBZ = Nz / 5;

   dim3 dimGrid(numGBX,numGBY,numGBZ);   
   

   DevFirstOrder_LS<<<dimGrid, dimBlock>>>(
                                          d_dPhi,
                                          phiS_d,
                                          jbn_d,
                                          1.0/deltaX,
                                          1.0/deltaY,
                                          1.0/deltaZ,
                                          Nx,
                                          Ny,
                                          Nz 
                                          );
   printf("\n\n  Extrapolating on CUDA Device: \n");

   for(int itera = 1 ; itera <=10 ; itera++){
       extrapolKernel<<<dimGrid, dimBlock>>>(
                                          rs_d,         
                                          extVal_d, phiS_d, jbn_d, d_dPhi,
                                       1.0/deltaX, 1.0/deltaY, 1.0/deltaZ,
                                            Nx, Ny, Nz, 
                                            Flag
                                            );  
     
       RunGK_FirstS<<<dimGrid, dimBlock>>>( d_d, extVal_d, 
                                           dtext, rs_d, Nx, Ny, Nz);
     
       extrapolKernel<<<dimGrid, dimBlock>>>(
                                            rs_d, d_d, 
                                            phiS_d, jbn_d, d_dPhi,
                                       1.0/deltaX, 1.0/deltaY, 1.0/deltaZ,
                                            Nx, Ny, Nz, 
                                            Flag);
     
       RunGK_SecondS<<<dimGrid, dimBlock>>>( d_d, 
                                             extVal_d, d_d, 
                                             dtext, rs_d, 
                                             Nx, Ny, Nz);
     
     
       extrapolKernel<<<dimGrid, dimBlock>>>(
                                            rs_d, d_d, 
                                            phiS_d, jbn_d, d_dPhi,
                                       1.0/deltaX, 1.0/deltaY, 1.0/deltaZ,
                                            Nx, Ny, Nz, 
                                            Flag);
     
       RunGK_ThirdS<<<dimGrid, dimBlock>>>(  extVal_d, 
                                             extVal_d, d_d, 
                                             dtext, rs_d, 
                                             Nx, Ny, Nz);
   }
 // check for error
  cudaError_t error = cudaGetLastError();
  if(error != cudaSuccess)
  {
    // print the CUDA error message and exit
    printf("CUDA error: %s\n", cudaGetErrorString(error));
    exit(-1);
  }

   cudaFree(d_d);
   cudaFree(d_dPhi);
   cudaFree(rs_d);
 
 }
 
 
 
