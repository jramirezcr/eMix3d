#include<cuda.h>
#include<stdio.h>
#include"flujos.h"
#include"divergencia.h"
#include"extrapol.h"

void rsDivergence(
                 double* const rs_d,
                 double* const e_d,
                 double* const f_d,
                 double* const g_d,
                 double* const ex_d,
                 double* const fy_d,
                 double* const gz_d,
                 const double* const jbn_d,
                 unsigned int Nx,
                 unsigned int Ny,
                 unsigned int Nz,
                 const double deltaX,
                 const double deltaY,
                 const double deltaZ,
                 unsigned int ncp,
                 unsigned int itera
                 );


extern "C"{
void flujos(
            double* const um1,
            double* const um,
            double* const um02,
            const double* const u,
            const double* const ug,
            const double* const dcvel,
            const double* const dcvelg,
            const double* const ddvelp,
            const double* const ddvelpg,
            const double* const press,
            const double* const pressg,
            const double* const dTemp,
            const double* const dTempg,
            const double* const vis,
            const double* const visg,
            const double* const vis5,
            const double* const vis5g,
            const double* const jbn_p,
            const double* const jbn_n,
            const double* const xmesh,
            const double c1, 
            const double dt, 
            const double froude, 
            const double deltaX, 
            const double deltaY, 
            const double deltaZ, 
            const unsigned int Nx, const unsigned int Ny, const unsigned int Nz,
            int itera, int ncp
            )
{


   double *um_d, *um1_d, *um02_d, *u_d, *ug_d, *dcvel_d, *dcvelg_d,
          *ddvelp_d, *ddvelpg_d, *press_d, *pressg_d, *dTemp_d, *dTempg_d, 
          *vis_d, *visg_d, *vis5_d, *vis5g_d, *jbn_d, *xmesh_d;
  
   double *e_d, *f_d, *g_d; 
   double *ex_d, *fy_d, *gz_d, *rs_d; 

   unsigned int Offset = Nx*Ny*Nz;

   cudaMalloc((void**)&um_d, sizeof(double)*Offset*10);
   cudaMalloc((void**)&um1_d, sizeof(double)*Offset*10);
   cudaMalloc((void**)&um02_d, sizeof(double)*Offset*10);
   cudaMalloc((void**)&u_d, sizeof(double)*Offset*3);
   cudaMalloc((void**)&ug_d, sizeof(double)*Offset*3);
   cudaMalloc((void**)&dcvel_d, sizeof(double)*Offset*9);
   cudaMalloc((void**)&dcvelg_d, sizeof(double)*Offset*9);
   cudaMalloc((void**)&ddvelp_d, sizeof(double)*Offset*9);
   cudaMalloc((void**)&ddvelpg_d, sizeof(double)*Offset*9);
   cudaMalloc((void**)&dTemp_d, sizeof(double)*Offset*3);
   cudaMalloc((void**)&dTempg_d, sizeof(double)*Offset*3);
   cudaMalloc((void**)&press_d, sizeof(double)*Offset);
   cudaMalloc((void**)&pressg_d, sizeof(double)*Offset);
   cudaMalloc((void**)&vis_d, sizeof(double)*Offset);
   cudaMalloc((void**)&visg_d, sizeof(double)*Offset);
   cudaMalloc((void**)&vis5_d, sizeof(double)*Offset);
   cudaMalloc((void**)&vis5g_d, sizeof(double)*Offset);
   cudaMalloc((void**)&xmesh_d, sizeof(double)*3*Offset);

   cudaMalloc((void**)&e_d, sizeof(double)*2*Offset);
   cudaMalloc((void**)&f_d, sizeof(double)*2*Offset);
   cudaMalloc((void**)&g_d, sizeof(double)*2*Offset);

   cudaMalloc((void**)&ex_d, sizeof(double)*2*Offset);
   cudaMalloc((void**)&fy_d, sizeof(double)*2*Offset);
   cudaMalloc((void**)&gz_d, sizeof(double)*2*Offset);
   cudaMalloc((void**)&rs_d, sizeof(double)*Offset);

   cudaMalloc((void**)&jbn_d, sizeof(double)*Offset*11);

   cudaMemcpy(um_d, um,  sizeof(double)*Offset*10, cudaMemcpyHostToDevice);
   cudaMemcpy(um02_d, um02,  sizeof(double)*Offset*10, cudaMemcpyHostToDevice);
   cudaMemcpy(u_d, u,  sizeof(double)*Offset*3, cudaMemcpyHostToDevice);
   cudaMemcpy(ug_d, ug,  sizeof(double)*Offset*3, cudaMemcpyHostToDevice);
   cudaMemcpy(dcvel_d, dcvel,  sizeof(double)*Offset*9, cudaMemcpyHostToDevice);
   cudaMemcpy(dcvelg_d, dcvelg,  sizeof(double)*Offset*9, cudaMemcpyHostToDevice);
   cudaMemcpy(ddvelp_d, ddvelp,  sizeof(double)*Offset*9, cudaMemcpyHostToDevice);
   cudaMemcpy(ddvelpg_d, ddvelpg,  sizeof(double)*Offset*9, cudaMemcpyHostToDevice);
   cudaMemcpy(press_d, press,  sizeof(double)*Offset, cudaMemcpyHostToDevice);
   cudaMemcpy(pressg_d, pressg,  sizeof(double)*Offset, cudaMemcpyHostToDevice);
   cudaMemcpy(dTemp_d, dTemp,  sizeof(double)*Offset*3, cudaMemcpyHostToDevice);
   cudaMemcpy(dTempg_d, dTempg,  sizeof(double)*Offset*3, cudaMemcpyHostToDevice);
   cudaMemcpy(vis_d, vis,  sizeof(double)*Offset, cudaMemcpyHostToDevice);
   cudaMemcpy(visg_d, visg,  sizeof(double)*Offset, cudaMemcpyHostToDevice);
   cudaMemcpy(vis5_d, vis5,  sizeof(double)*Offset, cudaMemcpyHostToDevice);
   cudaMemcpy(vis5g_d, vis5g,  sizeof(double)*Offset, cudaMemcpyHostToDevice);
   cudaMemcpy(xmesh_d, xmesh,  sizeof(double)*3*Offset, cudaMemcpyHostToDevice);


   const double cdiv = 2.0/3.0;
   
   if((ncp + itera)%2){
     cudaMemcpy(jbn_d, jbn_n,  sizeof(double)*Offset*11, cudaMemcpyHostToDevice);
   }
   else{
     cudaMemcpy(jbn_d, jbn_p,  sizeof(double)*Offset*11, cudaMemcpyHostToDevice);
   }

   dim3 DimGrid(12,12,24);
   dim3 DimBlock(8,8,4);


   for(int num_flujo = 0; num_flujo <= 9 ; num_flujo++){

      if(num_flujo == 0){//continuity liquid
        flux_continuity_CUDA<<<DimGrid, DimBlock>>>( e_d, f_d , g_d, um_d, Nx, Ny, Nz, c1);
      }

      
      else if(num_flujo == 1){//momentum X liquid
         flux_momentumX_CUDA<<<DimGrid, DimBlock>>>(e_d, f_d, g_d, 
                                                    u_d, um_d, press_d, dcvel_d, ddvelp_d,
                                                    vis_d, jbn_d, Nx, Ny, Nz, cdiv);
      }

      else if(num_flujo == 2){//momentum Y liquid
         flux_momentumY_CUDA<<<DimGrid, DimBlock>>>(e_d, f_d, g_d, 
                                                    u_d, um_d, press_d, dcvel_d, ddvelp_d,
                                                    vis_d, jbn_d, Nx, Ny, Nz, cdiv);
      }


      else if(num_flujo == 3){//momentum Z liquid
         flux_momentumZ_CUDA<<<DimGrid, DimBlock>>>(e_d, f_d, g_d,
                                                    u_d, um_d, press_d, dcvel_d,
                                                    ddvelp_d, vis_d, jbn_d, xmesh_d, 
                                                    Nx, Ny, Nz, cdiv, froude);
      }


      else if(num_flujo == 4){//Energy liquid
           flux_Energy_CUDA<<<DimGrid, DimBlock>>>(e_d, f_d, g_d,
                                                   u_d, um_d, dTemp_d, vis5_d, jbn_d,
                                                   Nx, Ny, Nz);
      }

      else if(num_flujo == 5){//continuity gas
        flux_continuity_CUDA<<<DimGrid, DimBlock>>>( e_d, f_d , g_d, &(um_d[5*Offset]), Nx, Ny, Nz, c1);
      }

      else if(num_flujo == 6){//momentum X gas
         flux_momentumX_CUDA<<<DimGrid, DimBlock>>>(e_d, f_d, g_d, 
                                             ug_d, &(um_d[5*Offset]), pressg_d, dcvelg_d, ddvelpg_d,
                                              visg_d, jbn_d, Nx, Ny, Nz, cdiv);
      }

      else if(num_flujo == 7){//momentum Y gas
         flux_momentumX_CUDA<<<DimGrid, DimBlock>>>(e_d, f_d, g_d, 
                             ug_d, &(um_d[5*Offset]), pressg_d, dcvelg_d, ddvelpg_d,
                             visg_d, jbn_d, Nx, Ny, Nz, cdiv);
      }

      else if(num_flujo == 8){//momentum Z liquid
         flux_momentumZ_CUDA<<<DimGrid, DimBlock>>>(e_d, f_d, g_d,
                                                    ug_d, &(um_d[5*Offset]), pressg_d, dcvelg_d,
                                                    ddvelpg_d, visg_d, jbn_d, xmesh_d, 
                                                    Nx, Ny, Nz, cdiv, 0.0);
      }

      else if(num_flujo == 9){//Energy liquid
           flux_Energy_CUDA<<<DimGrid, DimBlock>>>(e_d, f_d, g_d,
                                                   ug_d, &(um_d[5*Offset]), dTempg_d, vis5g_d, jbn_d,
                                                   Nx, Ny, Nz);
      }
    
      rsDivergence(rs_d, e_d, f_d, g_d, ex_d, fy_d, gz_d, jbn_d, Nx, Ny, Nz, 
                   deltaX, deltaY, deltaZ, ncp, itera);
   
      if(ncp==1){
          RunGK_FirstS<<<DimGrid,DimBlock>>>(&(um1_d[num_flujo*Offset]), &(um_d[num_flujo*Offset]), 
                                             dt, rs_d, Nx, Ny, Nz);
      }
      else if(ncp==2){
          RunGK_SecondS<<<DimGrid,DimBlock>>>(&(um1_d[num_flujo*Offset]), &(um02_d[num_flujo*Offset]), 
                                              &(um_d[num_flujo*Offset]), dt, rs_d, Nx, Ny, Nz);
      }

 
   }

   cudaMemcpy(um1, um1_d,  sizeof(double)*Offset*10, cudaMemcpyDeviceToHost);

   // check for error
   cudaError_t error = cudaGetLastError();
   if(error != cudaSuccess)
   {   
   // print the CUDA error message and exit
      printf("CUDA error: %s\n", cudaGetErrorString(error));
      exit(-1);
   }   
    
 
   return;
}
}

void rsDivergence(
                 double* const rs_d,
                 double* const e_d,
                 double* const f_d,
                 double* const g_d,
                 double* const ex_d,
                 double* const fy_d,
                 double* const gz_d,
                 const double* const jbn_d,
                 unsigned int Nx,
                 unsigned int Ny,
                 unsigned int Nz,
                 const double deltaX,
                 const double deltaY,
                 const double deltaZ,
                 unsigned int ncp,
                 unsigned int itera
                 )
{

   unsigned int Offset  = Nx*Ny*Nz;
 
   dim3 DimGrid(12,12,24);
   dim3 DimBlock(8,8,4);

   //Convective fluxes 
   get_flux_e_CUDA<<<DimGrid, DimBlock>>>(ex_d, e_d, f_d, g_d, jbn_d, Nx, Ny, Nz);       
   get_flux_f_CUDA<<<DimGrid, DimBlock>>>(fy_d, e_d, f_d, g_d, jbn_d, Nx, Ny, Nz);       
   get_flux_g_CUDA<<<DimGrid, DimBlock>>>(gz_d, e_d, f_d, g_d, jbn_d, Nx, Ny, Nz);       

   //Diffusive fluxes
   get_flux_e_CUDA<<<DimGrid, DimBlock>>>(&(ex_d[Offset]), &(e_d[Offset]), &(f_d[Offset]), 
                                          &(g_d[Offset]), jbn_d, Nx, Ny, Nz);       
   get_flux_f_CUDA<<<DimGrid, DimBlock>>>(&(fy_d[Offset]), &(e_d[Offset]), &(f_d[Offset]), 
                                          &(g_d[Offset]), jbn_d, Nx, Ny, Nz);       
   get_flux_g_CUDA<<<DimGrid, DimBlock>>>(&(gz_d[Offset]), &(e_d[Offset]), &(f_d[Offset]), 
                                          &(g_d[Offset]), jbn_d, Nx, Ny, Nz);       
                   

   double cons1 = 7.0 / 6.0;
   double cons2 = -8.0 / 6.0;
   double cons3 = 1.0 / 6.0;

   if((ncp + itera)%2){

      divDevXMin<<<DimGrid, DimBlock>>>( e_d, ex_d, Nx, Ny, Nz, deltaX, cons1, cons2, cons3);
      divDevXMin<<<DimGrid, DimBlock>>>( &(e_d[Offset]), &(ex_d[Offset]), Nx, Ny, Nz, deltaX, cons1, cons2, cons3);
 
      divDevYMin<<<DimGrid, DimBlock>>>( f_d, fy_d, Nx, Ny, Nz, deltaY, cons1, cons2, cons3);
      divDevYMin<<<DimGrid, DimBlock>>>( &(f_d[Offset]), &(fy_d[Offset]), Nx, Ny, Nz, deltaY, cons1, cons2, cons3);
 
      divDevZMin<<<DimGrid, DimBlock>>>( g_d, gz_d, Nx, Ny, Nz, deltaZ, cons1, cons2, cons3);
      divDevZMin<<<DimGrid, DimBlock>>>( &(g_d[Offset]), &(gz_d[Offset]), Nx, Ny, Nz, deltaZ, cons1, cons2, cons3);

   }
   else{

      divDevXPlus<<<DimGrid, DimBlock>>>( e_d, ex_d, Nx, Ny, Nz, deltaX, cons1, cons2, cons3);
      divDevXPlus<<<DimGrid, DimBlock>>>( &(e_d[Offset]), &(ex_d[Offset]), Nx, Ny, Nz, deltaX, cons1, cons2, cons3);
 
      divDevYPlus<<<DimGrid, DimBlock>>>( f_d, fy_d, Nx, Ny, Nz, deltaY, cons1, cons2, cons3);
      divDevYPlus<<<DimGrid, DimBlock>>>( &(f_d[Offset]), &(fy_d[Offset]), Nx, Ny, Nz, deltaY, cons1, cons2, cons3);
 
      divDevZPlus<<<DimGrid, DimBlock>>>( g_d, gz_d, Nx, Ny, Nz, deltaZ, cons1, cons2, cons3);
      divDevZPlus<<<DimGrid, DimBlock>>>( &(g_d[Offset]), &(gz_d[Offset]), Nx, Ny, Nz, deltaZ, cons1, cons2, cons3);

   }

   rs_divergence_CUDA<<<DimGrid, DimBlock>>>(rs_d, e_d, f_d, g_d, jbn_d, Nx, Ny, Nz); 

   // check for error
   cudaError_t error = cudaGetLastError();
   if(error != cudaSuccess)
   {   
   // print the CUDA error message and exit
      printf("CUDA error: %s\n", cudaGetErrorString(error));
      exit(-1);
   }   
   
   return;
}

