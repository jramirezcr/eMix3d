#include<cuda.h>
#include<iostream>
#include <fstream>
#include<stdio.h>
#include<stdlib.h>
#include"LevelSet/lsTools.h"
#include"LevelSet/extrapol.h"
#include"LevelSet/lsTransportEc.h"


extern "C"{
/*
    dim3 dimBlock(5,5,5);
    dim3 dimGrid(32,32,32);

    dim3 dimBlockG(5,5,5);
    dim3 dimGridG(30,30,30);

    dim3 dimBlockB(5,5);
    dim3 dimGridB(32,32);
*/


void advect_CUDA(
                 double *press,
                 double *jbn,
                 double *sobject,
                 double  deltaX,
                 double  deltaY,
                 double  deltaZ,
                 int     Nx, 
                 int     Ny, 
                 int     Nz,
                 double  deltamin
               )
{


   int gcells = 5;

   double *d_jbn;
   double *d_jbnW;
   double *d_sobject;
   double *d_sobjectW;
   double *d_press;
   double *d_pressW;
   

   int NxG = Nx + 2*gcells,
       NyG = Ny + 2*gcells,
       NzG = Nz + 2*gcells; 

   int offset  = Nx*Ny*Nz;
   int offsetG = NxG*NyG*NzG;

    int numGBX, numGBY, numGBZ;

    dim3 dimBlock(10,10,5);

    numGBX = NxG / 10;
    numGBY = NyG / 10;
    numGBZ = NzG / 5;

    dim3 dimGrid(numGBX,numGBY,numGBZ);

    dim3 dimBlockG(10,10,5);

    numGBX = Nx / 10;
    numGBY = Ny / 10;
    numGBZ = Nz / 5;

    dim3 dimGridG(numGBX,numGBY,numGBZ);

    dim3 dimBlockB(10,10);

    numGBX = NxG / 10;
    numGBY = NyG / 10;

    dim3 dimGridB(numGBX,numGBY);

    cudaSetDevice(1);

   cudaMalloc((void**)&d_jbnW,      11*sizeof(double)*offset);
   cudaMalloc((void**)&d_sobjectW,     sizeof(double)*offset);
   cudaMalloc((void**)&d_pressW,       sizeof(double)*offset);

   cudaMalloc((void**)&d_jbn,     11*sizeof(double)*offsetG);
   cudaMalloc((void**)&d_sobject,    sizeof(double)*offsetG);
   cudaMalloc((void**)&d_press,      sizeof(double)*offsetG);

   cudaMemcpy(d_jbnW,     jbn,    11*sizeof(double)*offset,
              cudaMemcpyHostToDevice);
   cudaMemcpy(d_sobjectW, sobject,   sizeof(double)*offset,
              cudaMemcpyHostToDevice);
   cudaMemcpy(d_pressW,press,        sizeof(double)*offset,
              cudaMemcpyHostToDevice);



   cuSwapToGhost<<<dimGridG, dimBlockG>>>
                 (d_press, d_pressW, gcells, Nx, Ny,Nz);

   cuSwapToGhost<<<dimGridG, dimBlockG>>>
            (&(d_jbn[0*offsetG]), &(d_jbnW[0*offset]), gcells, Nx, Ny,Nz);

   cuSwapToGhost<<<dimGridG, dimBlockG>>>
            (&(d_jbn[4*offsetG]), &(d_jbnW[4*offset]), gcells, Nx, Ny,Nz);

   cuSwapToGhost<<<dimGridG, dimBlockG>>>
            (&(d_jbn[8*offsetG]), &(d_jbnW[8*offset]), gcells, Nx, Ny,Nz);

   cuSwapToGhost<<<dimGridG, dimBlockG>>>
                 (d_sobject, d_sobjectW, gcells, Nx, Ny,Nz);


   //Boundary Conditions Ghost Cells press

   cuGhostCellsMirror3dZ<<<dimGridB, dimBlockB>>>
                        (d_press, gcells, Nx, Ny, Nz, 1.0);

   cuGhostCellsMirror3dY<<<dimGridB, dimBlockB>>>
                        (d_press, gcells, Nx, Ny, Nz, 1.0);

   cuGhostCellsMirror3dX<<<dimGridB, dimBlockB>>>
                        (d_press, gcells, Nx, Ny, Nz, 1.0);

   //Boundary Conditions Ghost Cells geometry

   cuGhostCellsMirror3dZ<<<dimGridB, dimBlockB>>>
                        (d_sobject, gcells, Nx, Ny, Nz, 1.0);

   cuGhostCellsMirror3dY<<<dimGridB, dimBlockB>>>
                        (d_sobject, gcells, Nx, Ny, Nz, 1.0);

   cuGhostCellsMirror3dX<<<dimGridB, dimBlockB>>>
                        (d_sobject, gcells, Nx, Ny, Nz, 1.0);


   //Ghost cells Jacobean

   cuGhostCellsMirror3dX<<<dimGridB, dimBlockB>>>
                        (&(d_jbn[0*offsetG]), gcells, Nx, Ny, Nz,  1.0);

   cuGhostCellsMirror3dY<<<dimGridB, dimBlockB>>>
                        (&(d_jbn[0*offsetG]), gcells, Nx, Ny, Nz,  1.0);

   cuGhostCellsMirror3dZ<<<dimGridB, dimBlockB>>>
                        (&(d_jbn[0*offsetG]), gcells, Nx, Ny, Nz,  1.0);

   cuGhostCellsMirror3dX<<<dimGridB, dimBlockB>>>
                        (&(d_jbn[4*offsetG]), gcells, Nx, Ny, Nz,  1.0);

   cuGhostCellsMirror3dY<<<dimGridB, dimBlockB>>>
                        (&(d_jbn[4*offsetG]), gcells, Nx, Ny, Nz,  1.0);

   cuGhostCellsMirror3dZ<<<dimGridB, dimBlockB>>>
                        (&(d_jbn[4*offsetG]), gcells, Nx, Ny, Nz,  1.0);

   cuGhostCellsMirror3dX<<<dimGridB, dimBlockB>>>
                        (&(d_jbn[8*offsetG]), gcells, Nx, Ny, Nz,  1.0);

   cuGhostCellsMirror3dY<<<dimGridB, dimBlockB>>>
                        (&(d_jbn[8*offsetG]), gcells, Nx, Ny, Nz,  1.0);

   cuGhostCellsMirror3dZ<<<dimGridB, dimBlockB>>>
                        (&(d_jbn[8*offsetG]), gcells, Nx, Ny, Nz,  1.0);


   //Extrapolating pressure and ls (geometry boundary condition)
   cuExtrapolation(d_press,d_sobject,d_jbn,deltaX,deltaY,deltaZ, 
                   NxG, NyG, NzG, deltamin*0.5, 1.0);

   cuSwapFromGhost<<<dimGridG, dimBlockG>>>
                  (d_pressW, d_press, gcells, Nx, Ny,Nz);

   //Returning values from gpu to cpu (phi, vel, press)
   cudaMemcpy(press,d_pressW,sizeof(double)*offset, cudaMemcpyDeviceToHost);

   cudaFree(d_jbn);
   cudaFree(d_jbnW);
   cudaFree(d_sobject);
   cudaFree(d_sobjectW);
   cudaFree(d_press);
   cudaFree(d_pressW);

   // check for error
   cudaError_t error = cudaGetLastError();
   if(error != cudaSuccess)
   {
     // print the CUDA error message and exit
     std::cout << "CUDA error extrapol:\n" <<  cudaGetErrorString(error) 
               << std::endl;
     exit(-1);
   }


} 

}// end extern C

void cuAdvectLsJB(
                 double *d_phi, 
                 double *d_vel,
                 double *d_jbn,
                 double  deltaX,
                 double  deltaY,
                 double  deltaZ,
                 double  dt,
                 int     NxG, 
                 int     NyG, 
                 int     NzG
                 )
{

    cudaSetDevice(1);

    //Derived of Level-Set function
    double *d_dPhi;
    double *d_rsPhi;
    double *d_phiTemp;
    double *d_dPhiPlus;
    double *d_dPhiMinus;

    int offsetG = NxG*NyG*NzG; 

    int numGBX,numGBY,numGBZ;
   
    dim3 dimBlock(10,10,5);

    numGBX = NxG / 10;
    numGBY = NyG / 10;
    numGBZ = NxG / 5;

    dim3 dimGrid(numGBX,numGBY,numGBZ);

    cudaMalloc((void**)&d_dPhi,      3*sizeof(double)*offsetG);
    cudaMalloc((void**)&d_dPhiPlus,  3*sizeof(double)*offsetG);
    cudaMalloc((void**)&d_dPhiMinus, 3*sizeof(double)*offsetG);
    cudaMalloc((void**)&d_rsPhi,       sizeof(double)*offsetG);
    cudaMalloc((void**)&d_phiTemp,     sizeof(double)*offsetG);


    //Runge First Step
    Dev1thO_Downwind<<<dimGrid,dimBlock>>>(
                                          d_dPhi,
                                          d_phi,
                                          1.0/deltaX,
                                          1.0/deltaY,
                                          1.0/deltaZ,
                                          NxG, 
                                          NyG, 
                                          NzG
                                          );


    PhiDevPlusParameterJB<<<dimGrid,dimBlock>>>(
                                              d_dPhiPlus,
                                              d_dPhi,
                                              d_jbn,
                                              NxG, 
                                              NyG, 
                                              NzG 
                                             );  


    PhiDevMinusParameterJB<<<dimGrid,dimBlock>>>( d_dPhiMinus,
                                                  d_dPhi,
                                                  d_jbn,
                                               NxG, 
                                               NyG, 
                                               NzG 
                                              );  

    advect_RS_WENO<<<dimGrid,dimBlock>>>(
                                         d_rsPhi, 
                                         d_vel,
                                         d_dPhiPlus,
                                         d_dPhiMinus,
                                         NxG, 
                                         NyG, 
                                         NzG
                                        );

     RunGK_FirstS<<<dimGrid,dimBlock>>>(
                                        d_phiTemp,
                                        d_phi, 
                                        dt, 
                                        d_rsPhi, 
                                        NxG, 
                                        NyG, 
                                        NzG
                                       );  

    //Second step Runge-Kutta
    Dev1thO_Downwind<<<dimGrid,dimBlock>>>(
                                          d_dPhi,
                                          d_phiTemp,
                                          1.0/deltaX,
                                          1.0/deltaY,
                                          1.0/deltaZ,
                                          NxG, 
                                          NyG, 
                                          NzG
                                          );


    PhiDevMinusParameterJB<<<dimGrid,dimBlock>>>(
                                               d_dPhiMinus,
                                               d_dPhi,
                                               d_jbn,
                                               NxG, 
                                               NyG, 
                                               NzG 
                                              );  

    PhiDevPlusParameterJB<<<dimGrid,dimBlock>>>(
                                              d_dPhiPlus,
                                              d_dPhi,
                                              d_jbn,
                                              NxG, 
                                              NyG, 
                                              NzG 
                                             );  


    advect_RS_WENO<<<dimGrid,dimBlock>>>(
                                         d_rsPhi, 
                                         d_vel,
                                         d_dPhiPlus,
                                         d_dPhiMinus,
                                         NxG, 
                                         NyG, 
                                         NzG
                                        );


    RunGK_SecondS<<<dimGrid, dimBlock>>>(
                                         d_phiTemp, 
                                         d_phi, 
                                         d_phiTemp, 
                                         dt,
                                         d_rsPhi, 
                                         NxG, 
                                         NyG, 
                                         NzG
                                        );

    //Third step Runge-Kutta
    Dev1thO_Downwind<<<dimGrid,dimBlock>>>(
                                          d_dPhi,
                                          d_phiTemp,
                                          1.0/deltaX,
                                          1.0/deltaY,
                                          1.0/deltaZ,
                                          NxG, 
                                          NyG, 
                                          NzG
                                          );


    PhiDevMinusParameterJB<<<dimGrid,dimBlock>>>(
                                               d_dPhiMinus,
                                               d_dPhi,
                                               d_jbn,
                                               NxG, 
                                               NyG, 
                                               NzG 
                                              );  


    PhiDevPlusParameterJB<<<dimGrid,dimBlock>>>(
                                              d_dPhiPlus,
                                              d_dPhi,
                                              d_jbn,
                                              NxG, 
                                              NyG, 
                                              NzG 
                                             );  

    advect_RS_WENO<<<dimGrid,dimBlock>>>(
                                         d_rsPhi, 
                                         d_vel,
                                         d_dPhiPlus,
                                         d_dPhiMinus,
                                         NxG, 
                                         NyG, 
                                         NzG
                                        );

    RunGK_ThirdS<<<dimGrid, dimBlock>>>(
                                        d_phi, 
                                        d_phi, 
                                        d_phiTemp, 
                                        dt,
                                        d_rsPhi, 
                                        NxG, 
                                        NyG, 
                                        NzG
                                        );


   // check for error
   cudaError_t error = cudaGetLastError();
   if(error != cudaSuccess)
   {
     // print the CUDA error message and exit
     std::cout << "CUDA error cuAdvectLsJB:\n" <<  cudaGetErrorString(error) 
               << std::endl;
     exit(-1);
   }

   cudaFree(d_dPhi);
   cudaFree(d_phiTemp);
   cudaFree(d_rsPhi);
   cudaFree(d_dPhiPlus);
   cudaFree(d_dPhiMinus);

}

     

void cuReinitLsJB(
                 double *d_phi,    // Level set function on DEVICE
                 double *d_jbn,    // Level set function on DEVICE
                 double  deltaX,
                 double  deltaY,
                 double  deltaZ,
                 double *d_dmins,    // Level set function on DEVICE
                 int     Nx, 
                 int     Ny, 
                 int     Nz,
                 int     gcells,
                 double  deltamin
                 )
{


    int numGBX,numGBY,numGBZ;
   
    dim3 dimBlock(10,10,5);

    numGBX = Nx / 10;
    numGBY = Ny / 10;
    numGBZ = Nz / 5;

    dim3 dimGrid(numGBX,numGBY,numGBZ);

    dim3 dimBlockB(10,10);

    numGBX = Nx / 10;
    numGBY = Ny / 10;

    dim3 dimGridB(numGBX,numGBY);

    //Derived of Level-Set function
    double *d_dPhi;
    double *d_phi0;
    double *d_rsPhi;
    double *d_phiTemp;
    double *d_dPhiPlus;
    double *d_dPhiMinus;

    cudaMalloc((void**)&d_dPhi,      3*sizeof(double)*Nx*Ny*Nz);
    cudaMalloc((void**)&d_phi0,        sizeof(double)*Nx*Ny*Nz);
    cudaMalloc((void**)&d_dPhiPlus,  3*sizeof(double)*Nx*Ny*Nz);
    cudaMalloc((void**)&d_dPhiMinus, 3*sizeof(double)*Nx*Ny*Nz);
    cudaMalloc((void**)&d_rsPhi,       sizeof(double)*Nx*Ny*Nz);
    cudaMalloc((void**)&d_phiTemp,     sizeof(double)*Nx*Ny*Nz);
  

    for(int itera = 1; itera <= 10; itera++){

        //Boundary Conditions Ghost Cells Phi

        cuGhostCellsMirror3dZ<<<dimGridB, dimBlockB>>>
        (d_phi, gcells, Nx-2*gcells, Ny-2*gcells, Nz-2*gcells, 1.0);

        cuGhostCellsMirror3dY<<<dimGridB, dimBlockB>>>
        (d_phi, gcells, Nx-2*gcells, Ny-2*gcells, Nz-2*gcells, 1.0);

        cuGhostCellsMirror3dX<<<dimGridB, dimBlockB>>>
        (d_phi, gcells, Nx-2*gcells, Ny-2*gcells, Nz-2*gcells, 1.0);


        cudaMemcpy(d_phi0,d_phi,sizeof(double)*Nx*Ny*Nz,
                  cudaMemcpyDeviceToDevice );


        //Runge First Step
        Dev1thO_Downwind<<<dimGrid,dimBlock>>>(
                                              d_dPhi,
                                              d_phi,
                                              1.0/deltaX,
                                              1.0/deltaY,
                                              1.0/deltaZ,
                                              Nx, 
                                              Ny, 
                                              Nz
                                              );
       
       
        PhiDevPlusParameterJB<<<dimGrid,dimBlock>>>(
                                                  d_dPhiPlus,
                                                  d_dPhi,
                                                  d_jbn,
                                                  Nx, 
                                                  Ny, 
                                                  Nz 
                                                 );  
       
       
        PhiDevMinusParameterJB<<<dimGrid,dimBlock>>>(
                                                   d_dPhiMinus,
                                                   d_dPhi,
                                                   d_jbn,
                                                   Nx, 
                                                   Ny, 
                                                   Nz 
                                                  );  
       
        reini_RS_WENOJB<<<dimGrid,dimBlock>>>(
                                             d_rsPhi, 
                                             d_phi, 
                                             d_dmins,
                                             d_dPhiPlus,
                                             d_dPhiMinus,
                                             d_phi0, 
                                             Nx, 
                                             Ny, 
                                             Nz
                                            );
       
         RunGK_FirstS<<<dimGrid,dimBlock>>>(
                                            d_phiTemp,
                                            d_phi, 
                                            0.5*deltamin,
                                            d_rsPhi, 
                                            Nx, 
                                            Ny, 
                                            Nz
                                           );  
       
        //Second step Runge-Kutta
        Dev1thO_Downwind<<<dimGrid,dimBlock>>>(
                                              d_dPhi,
                                              d_phiTemp,
                                              1.0/deltaX,
                                              1.0/deltaY,
                                              1.0/deltaZ,
                                              Nx, 
                                              Ny, 
                                              Nz
                                              );
       
       
        PhiDevMinusParameterJB<<<dimGrid,dimBlock>>>(
                                                   d_dPhiMinus,
                                                   d_dPhi,
                                                   d_jbn,
                                                   Nx, 
                                                   Ny, 
                                                   Nz 
                                                  );  
       
        PhiDevPlusParameterJB<<<dimGrid,dimBlock>>>(
                                                  d_dPhiPlus,
                                                  d_dPhi,
                                                  d_jbn,
                                                  Nx, 
                                                  Ny, 
                                                  Nz 
                                                 );  
       
       
        reini_RS_WENOJB<<<dimGrid,dimBlock>>>(
                                             d_rsPhi, 
                                             d_phiTemp,
                                             d_dmins,
                                             d_dPhiPlus,
                                             d_dPhiMinus,
                                             d_phi0, 
                                             Nx, 
                                             Ny, 
                                             Nz
                                            );
       
       
        RunGK_SecondS<<<dimGrid, dimBlock>>>(
                                             d_phiTemp, 
                                             d_phi, 
                                             d_phiTemp, 
                                             0.5*deltamin,
                                             d_rsPhi, 
                                             Nx, 
                                             Ny, 
                                             Nz
                                            );
       
        //Third step Runge-Kutta
        Dev1thO_Downwind<<<dimGrid,dimBlock>>>(
                                              d_dPhi,
                                              d_phiTemp,
                                              1.0/deltaX,
                                              1.0/deltaY,
                                              1.0/deltaZ,
                                              Nx, 
                                              Ny, 
                                              Nz
                                              );
       
       
        PhiDevMinusParameterJB<<<dimGrid,dimBlock>>>(
                                                   d_dPhiMinus,
                                                   d_dPhi,
                                                   d_jbn,
                                                   Nx, 
                                                   Ny, 
                                                   Nz 
                                                  );  
       
       
        PhiDevPlusParameterJB<<<dimGrid,dimBlock>>>(
                                                  d_dPhiPlus,
                                                  d_dPhi,
                                                  d_jbn,
                                                  Nx, 
                                                  Ny, 
                                                  Nz 
                                                 );  

       
        reini_RS_WENOJB<<<dimGrid,dimBlock>>>(
                                             d_rsPhi, 
                                             d_phiTemp,
                                             d_dmins,
                                             d_dPhiPlus,
                                             d_dPhiMinus,
                                             d_phi0, 
                                             Nx, 
                                             Ny, 
                                             Nz
                                            );
       

        RunGK_ThirdS<<<dimGrid, dimBlock>>>(
                                            d_phi, 
                                            d_phi, 
                                            d_phiTemp, 
                                            0.5*deltamin,
                                            d_rsPhi, 
                                            Nx, 
                                            Ny, 
                                            Nz
                                        );

   }

   // check for error
   cudaError_t error = cudaGetLastError();
   if(error != cudaSuccess)
   {
     // print the CUDA error message and exit
     std::cout << "CUDA error cuReinitLsJB: \n" <<  cudaGetErrorString(error) 
               << std::endl;
     exit(-1);
   }

   cudaFree(d_dPhi);
   cudaFree(d_phi0);
   cudaFree(d_phiTemp);
   cudaFree(d_rsPhi);
   cudaFree(d_dPhiPlus);
   cudaFree(d_dPhiMinus);
 

}


