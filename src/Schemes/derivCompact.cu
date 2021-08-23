#include<iostream>
#include<stdio.h>
#include"Engine/cuArray.hpp"

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

template<typename Tprec>
__global__
void transposeArrayJ(ArrayDev3D<Tprec> a, ArrayDev3D<Tprec> b){

   int i = blockDim.x*blockIdx.x + threadIdx.x;
   int j = blockDim.y*blockIdx.y + threadIdx.y;
   int k = blockDim.z*blockIdx.z + threadIdx.z;

   if(i < a.getDim(1) && j < a.getDim(2) && k < a.getDim(3)){

      a(i,j,k) = b(j,i,k);

   }

}

template<typename Tprec>
__global__
void transposeArrayK(ArrayDev3D<Tprec> a, ArrayDev3D<Tprec> b){

   int i = blockDim.x*blockIdx.x + threadIdx.x;
   int j = blockDim.y*blockIdx.y + threadIdx.y;
   int k = blockDim.z*blockIdx.z + threadIdx.z;

   if(i < a.getDim(1) && j < a.getDim(2) && k < a.getDim(3)){

      a(i,j,k) = b(k,j,i);

   }
}



template<typename Tprec>
__global__
void getC4RHS(ArrayDev3D<Tprec> f,  ArrayDev3D<Tprec> r, Tprec delta){

   int idx = blockDim.x*blockIdx.x + threadIdx.x,
       idy = blockDim.y*blockIdx.y + threadIdx.y,
       idz = blockDim.z*blockIdx.z + threadIdx.z;

   Tprec ars  =   3.0/4.0,
         ars1 = -17.0/6.0,
         brs1 =   3.0/2.0,
         crs1 =   3.0/2.0,
         drs1 =  -1.0/6.0;

   if(idx > 0 && idx < r.getDim(1) - 1 && idy < r.getDim(2) && 
      idz < r.getDim(3)){
     
      f(idx,idy,idz) =   (1.0/delta)*(ars*(r(idx+1,idy,idz) 
                                          -r(idx-1,idy,idz))); 
   }

   if(idx == 0 && idy < r.getDim(2) && idz < r.getDim(3)){
        f(idx,idy,idz) = (1.0/delta)*(ars1*r(0,idy,idz)
                                     +brs1*r(1,idy,idz)
                                     +crs1*r(2,idy,idz) 
                                     +drs1*r(3,idy,idz)
                                     );
   }

   if(idx == r.getDim(1) - 1 && idy < r.getDim(2) && idz < r.getDim(3)){
        f(idx,idy,idz) = (1.0/delta)*(-ars1*r(r.getDim(1) - 1,idy,idz)
                                      -brs1*r(r.getDim(1) - 2,idy,idz)
                                      -crs1*r(r.getDim(1) - 3,idy,idz) 
                                      -drs1*r(r.getDim(1) - 4,idy,idz)
                                     );
   }

}


template<typename Tprec>
__global__
void getC6RHS(ArrayDev3D<Tprec> f,  ArrayDev3D<Tprec> r, Tprec delta){

   int idx = blockDim.x*blockIdx.x + threadIdx.x,
       idy = blockDim.y*blockIdx.y + threadIdx.y,
       idz = blockDim.z*blockIdx.z + threadIdx.z;


   Tprec ars  =   14.0/(9.0*2.0),
         brs  =   1.0/(9.0*4.0),
         ars1 = -17.0/6.0,
         brs1 =   3.0/2.0,
         crs1 =   3.0/2.0,
         drs1 =  -1.0/6.0,
         ars2 =  -3.0/4.0,
         brs2 =       0.0,
         crs2 =   3.0/4.0,
         drs2 =       0.0;


   if(idx > 1 && idx < r.getDim(1) - 2 && idy < r.getDim(2) && 
      idz < r.getDim(3)){
     
      f(idx,idy,idz) =   (1.0/delta)*(ars*(r(idx+1,idy,idz)-r(idx-1,idy,idz)) 
                          + (brs*(r(idx+2,idy,idz)-r(idx-2,idy,idz)))); 
   }

   if(idx == 0 && idy < r.getDim(2) && idz < r.getDim(3)){
        f(idx,idy,idz) = (1.0/delta)*(ars1*r(0,idy,idz)
                                     +brs1*r(1,idy,idz)
                                     +crs1*r(2,idy,idz) 
                                     +drs1*r(3,idy,idz)
                                     );
   }


   if(idx == 1 && idy < r.getDim(2) && idz < r.getDim(3)){
        f(idx,idy,idz) = (1.0/delta)*(ars2*r(0,idy,idz)
                                     +brs2*r(1,idy,idz)
                                     +crs2*r(2,idy,idz) 
                                     +drs2*r(3,idy,idz)
                                     );
   }

   if(idx == r.getDim(1) - 1 && idy < r.getDim(2) && idz < r.getDim(3)){
        f(idx,idy,idz) = (1.0/delta)*(-ars1*r(r.getDim(1) - 1,idy,idz)
                                      -brs1*r(r.getDim(1) - 2,idy,idz)
                                      -crs1*r(r.getDim(1) - 3,idy,idz) 
                                      -drs1*r(r.getDim(1) - 4,idy,idz)
                                     );
   }

   if(idx == r.getDim(1) - 2 && idy < r.getDim(2) && idz < r.getDim(3)){
        f(idx,idy,idz) = (1.0/delta)*(-ars2*r(r.getDim(1) - 1,idy,idz)
                                      -brs2*r(r.getDim(1) - 2,idy,idz)
                                      -crs2*r(r.getDim(1) - 3,idy,idz) 
                                      -drs2*r(r.getDim(1) - 4,idy,idz)
                                     );
   }

}



template<typename Tprec>
__global__
void cuDerivCS(ArrayDev3D<Tprec> f, 
             ArrayDev3D<Tprec> r,
             ArrayDev3D<Tprec> gam
            ){

   extern __shared__ Tprec diag[]; 

   int block     = blockDim.x*blockDim.y;
   int blockloop = 3*f.getDim(1) / block; 

   int id_inside;
   
   if(3*f.getDim(1) % block) blockloop += 1;

   for(int it = 0; it < blockloop; it++){
     id_inside = it*block + blockDim.x*threadIdx.y + threadIdx.x;

     //Diagonal a matrix
     if(id_inside < f.getDim(1)){
       diag[id_inside] = 1.0/3.0; 
     } 
     //Diagonal b matrix
     if(id_inside >= f.getDim(1) && id_inside < 2*f.getDim(1)){
       diag[id_inside] = 1.0; 
     } 

     //Diagonal c matrix
     if(id_inside >= 2*f.getDim(1) && id_inside < 3*f.getDim(1)){
       diag[id_inside] = 1.0/3.0; 
     } 

//   A cof
     if(id_inside == 0            ) diag[id_inside] = 0.0;
     if(id_inside == 1            ) diag[id_inside] = 1.0 / 4.0;
     if(id_inside == f.getDim(1) - 1) diag[id_inside] = 3.0;
     if(id_inside == f.getDim(1) - 2) diag[id_inside] = 1.0/4.0;

//   C cof
     if(id_inside == 2*f.getDim(1)    ) diag[id_inside] = 3.0;                
     if(id_inside == 2*f.getDim(1) + 1) diag[id_inside] = 1.0/4.0;
     if(id_inside == 3*f.getDim(1) - 1) diag[id_inside] = 0.0;
     if(id_inside == 3*f.getDim(1) - 2) diag[id_inside] = 1.0/4.0;
   }
   __syncthreads();


   int idx = blockDim.x*blockIdx.x + threadIdx.x;  
   int idy = blockDim.y*blockIdx.y + threadIdx.y;  

   if(idx < f.getDim(2) && idy < f.getDim(3)){

      Tprec bet = diag[f.getDim(1)];
    
      f(0,idx,idy) = r(0,idx,idy) / bet; 
      for(int i = 1; i <f.getDim(1); i++){  
           gam(i,idx,idy) = diag[2*f.getDim(1) + i -1] / bet;
           bet = diag[f.getDim(1) + i] - diag[i]*gam(i,idx,idy);
    
         f(i,idx,idy) =  (r(i,idx,idy) - diag[i]*f(i-1,idx,idy))/bet;
      }
    
      for(int i = f.getDim(1) - 2; i >= 0 ; i--){  
         f(i,idx,idy) = f(i,idx,idy) - gam(i+1,idx,idy)*f(i+1,idx,idy);
    
      }

   }
}


//Deriv functions
namespace cu{
namespace cs{

template<typename Tprec>
void dFdX(
            ArrayDev3D<Tprec>  &fprima, 
            ArrayDev3D<Tprec>  &f, 
            Tprec              delta
           ){

    Tprec* rhs_d;

    cudaMalloc((void **)&rhs_d, f.getDim(1)*f.getDim(2)*f.getDim(3)*sizeof(Tprec));

    ArrayDev3D<Tprec> rhs(f.getDim(1), f.getDim(2), f.getDim(3)),
                      trash(f.getDim(1), f.getDim(2), f.getDim(3));


//--configuring launch options (fermi card)


    int sharedSize = 3*rhs.getDim(1)*sizeof(Tprec);

    int gridDimX, gridDimY, gridDimZ;
    int gridDimYDer, gridDimZDer;
   
    gridDimX = f.getDim(1) / 8; 
    if(f.getDim(1) % 8) gridDimX++;

    gridDimY = f.getDim(2) / 8; 
    if(f.getDim(2) % 8) gridDimY++;

    gridDimYDer = f.getDim(2) / 16; 
    if(f.getDim(2) % 16) gridDimYDer++;

    gridDimZ = f.getDim(3) / 4; 
    if(f.getDim(3) % 4) gridDimZ++;

    gridDimZDer = f.getDim(3) / 16; 
    if(f.getDim(3) % 16) gridDimZDer++;

    dim3 BlockDim(8,8,4);
    dim3 GridDim(gridDimX,gridDimY,gridDimZ);

    dim3 BlockDimDer(16,16);
    dim3 GridDimDer(gridDimYDer,gridDimZDer);


//--end configuring launch options (fermi card)

      getC6RHS<<<GridDim,BlockDim>>>(rhs, f, delta);
      cuDerivCS<<<GridDimDer,BlockDimDer,sharedSize>>>(fprima,
                                                       rhs, 
                                                       trash);

     cudaError_t error = cudaGetLastError();
     if(error != cudaSuccess)
     {   
        // print the CUDA error message and exit
        printf("CUDA dFdX error: %s\n", cudaGetErrorString(error));
        exit(-1);
     }  

     cudaFree(rhs_d);
 
}

template<typename Tprec>
void dFdY(
            ArrayDev3D<Tprec>  &fprima, 
            ArrayDev3D<Tprec>  &f, 
            Tprec               delta
           ){

    ArrayDev3D<Tprec> rhs(f.getDim(2), f.getDim(1), f.getDim(3)),
                      ftrans(f.getDim(2), f.getDim(1), f.getDim(3)),
                      trash(f.getDim(2), f.getDim(1), f.getDim(3));
    
//--configuring launch options (fermi card)
   
    int gridDimX, gridDimY, gridDimZ;
    int gridDimYDer, gridDimZDer;

    int sharedSize = 3*rhs.getDim(1)*sizeof(Tprec);

    gridDimX = ftrans.getDim(1) / 8; 
    if(ftrans.getDim(1) % 8) gridDimX++;

    gridDimY = ftrans.getDim(2) / 8; 
    if(ftrans.getDim(2) % 8) gridDimY++;

    gridDimYDer = ftrans.getDim(2) / 16; 
    if(ftrans.getDim(2) % 16) gridDimYDer++;

    gridDimZ = ftrans.getDim(3) / 4; 
    if(ftrans.getDim(3) % 4) gridDimZ++;

    gridDimZDer = ftrans.getDim(3) / 16; 
    if(ftrans.getDim(3) % 16) gridDimZDer++;

    dim3 BlockDim(8,8,4);
    dim3 GridDim(gridDimX,gridDimY,gridDimZ);

    dim3 BlockDimDer(16,16);
    dim3 GridDimDer(gridDimYDer,gridDimZDer);

//--end configuation

    transposeArrayJ<<<GridDim,BlockDim>>>(ftrans,f);

    getC6RHS<<<GridDim,BlockDim>>>(rhs, ftrans, delta);
    cuDerivCS<<<GridDimDer,BlockDimDer,sharedSize>>>(ftrans,
                                                     rhs, 
                                                     trash);

    GridDim = dim3(gridDimY,gridDimX,gridDimZ);
    transposeArrayJ<<<GridDim,BlockDim>>>(fprima, ftrans);

     cudaError_t error = cudaGetLastError();
     if(error != cudaSuccess)
     {   
        // print the CUDA error message and exit
        printf("CUDA dFdY error: %s\n", cudaGetErrorString(error));
        exit(-1);
     }  

}

template<typename Tprec>
void dFdZ(
          ArrayDev3D<Tprec>  &fprima, 
          ArrayDev3D<Tprec>  &f, 
          Tprec               delta
          ){

    ArrayDev3D<Tprec> rhs(f.getDim(3), f.getDim(2), f.getDim(1)),
                      ftrans(f.getDim(3), f.getDim(2), f.getDim(1)),
                      trash(f.getDim(3), f.getDim(2), f.getDim(1));
    
//--configuring launch options (fermi card)

   
    int gridDimX, gridDimY, gridDimZ;
    int gridDimYDer, gridDimZDer;

    int sharedSize = 3*rhs.getDim(1)*sizeof(Tprec);

    gridDimX = ftrans.getDim(1) / 8; 
    if(ftrans.getDim(1) % 8) gridDimX++;

    gridDimY = ftrans.getDim(2) / 8; 
    if(ftrans.getDim(2) % 8) gridDimY++;

    gridDimYDer = ftrans.getDim(2) / 16; 
    if(ftrans.getDim(2) % 16) gridDimYDer++;

    gridDimZ = ftrans.getDim(3) / 4; 
    if(ftrans.getDim(3) % 4) gridDimZ++;

    gridDimZDer = ftrans.getDim(3) / 16; 
    if(ftrans.getDim(3) % 16) gridDimZDer++;

    dim3 BlockDim(8,8,4);
    dim3 GridDim(gridDimX,gridDimY,gridDimZ);

    dim3 BlockDimDer(16,16);
    dim3 GridDimDer(gridDimYDer,gridDimZDer);

//--end configuation

    transposeArrayK<<<GridDim,BlockDim>>>(ftrans,f);

    getC6RHS<<<GridDim,BlockDim>>>(rhs, ftrans, delta);
    cuDerivCS<<<GridDimDer,BlockDimDer,sharedSize>>>(ftrans,
                                                     rhs, 
                                                     trash);

    BlockDim = dim3(4,8,8);
    GridDim  = dim3(gridDimZ,gridDimY,gridDimX);
    transposeArrayK<<<GridDim,BlockDim>>>(fprima, ftrans);

     cudaError_t error = cudaGetLastError();
     if(error != cudaSuccess)
     {   
        // print the CUDA error message and exit
        printf("CUDA dFdZ error: %s\n", cudaGetErrorString(error));
        exit(-1);
     }  


}

//end namespaces
}
}

