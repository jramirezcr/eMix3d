#include"Schemes/derivCompact.hpp"
#include"Engine/cuArray.hpp"

template<typename Tprec>
__global__ void kernProyectDataJBN(
                          ArrayDev3D<Tprec> FLUX,
                          ArrayDev3D<Tprec> jbndet,
                          ArrayDev3D<Tprec> jbn
                          ){


   int i = blockIdx.x*blockDim.x + threadIdx.x;
   int j = blockIdx.y*blockDim.y + threadIdx.y;
   int k = blockIdx.z*blockDim.z + threadIdx.z;

   if(i < FLUX.getDim(1) && j < FLUX.getDim(2) && k < FLUX.getDim(3)){
       FLUX(i,j,k) = jbndet(i,j,k)*jbn(i,j,k)*FLUX(i,j,k);
   }
}

template<typename Tprec>
__global__ void kernGetDiv(
                          ArrayDev3D<Tprec> DIV,
                          ArrayDev3D<Tprec> jbn,
                          ArrayDev3D<Tprec> FLUXE,
                          ArrayDev3D<Tprec> FLUXF,
                          ArrayDev3D<Tprec> FLUXG
                          ){

   int i = blockIdx.x*blockDim.x + threadIdx.x;
   int j = blockIdx.y*blockDim.y + threadIdx.y;
   int k = blockIdx.z*blockDim.z + threadIdx.z;

   if(i < DIV.getDim(1) && j < DIV.getDim(2) && k < DIV.getDim(3)){
       DIV(i,j,k) = (FLUXE(i,j,k) + FLUXF(i,j,k) + FLUXG(i,j,k))*
                    jbn(i,j,k);
   }
}

template<typename Tprec>
void getDiv(
            ArrayDev3D<Tprec> DIV,
            ArrayDev3D<Tprec> jbn,
            ArrayDev3D<Tprec> FLUXE,
            ArrayDev3D<Tprec> FLUXF,
            ArrayDev3D<Tprec> FLUXG
            ){

    int gridDimX, gridDimY, gridDimZ;
     
    gridDimX = DIV.getDim(1) / 8; 
    if(DIV.getDim(1) % 8) gridDimX++;

    gridDimY = DIV.getDim(2) / 8; 
    if(DIV.getDim(2) % 8) gridDimY++;

    gridDimZ = DIV.getDim(3) / 4; 
    if(DIV.getDim(3) % 4) gridDimZ++;

    dim3 BlockDim(8,8,4);
    dim3 GridDim(gridDimX,gridDimY,gridDimZ);

    kernGetDiv<<<GridDim,BlockDim>>>(DIV,jbn,FLUXE,FLUXF,FLUXG);

    cudaError_t error; 
    error = cudaGetLastError();
    if(error!=cudaSuccess){
      std::cout<< "cuda error getDiv: " << cudaGetErrorString(error) <<std::endl;
      exit(-1);
    
    }
}

template<typename Tprec>
void proyectDataJBN(
                    ArrayDev3D<Tprec> FLUX,
                    ArrayDev3D<Tprec> jbndet,
                    ArrayDev3D<Tprec> jbn
                   ){
    int gridDimX, gridDimY, gridDimZ;
     
    gridDimX = FLUX.getDim(1) / 8; 
    if(FLUX.getDim(1) % 8) gridDimX++;

    gridDimY = FLUX.getDim(2) / 8; 
    if(FLUX.getDim(2) % 8) gridDimY++;

    gridDimZ = FLUX.getDim(3) / 4; 
    if(FLUX.getDim(3) % 4) gridDimZ++;

    dim3 BlockDim(8,8,4);
    dim3 GridDim(gridDimX,gridDimY,gridDimZ);

    kernProyectDataJBN<<<GridDim, BlockDim>>>(FLUX, jbndet, jbn);

    cudaError_t error; 
    error = cudaGetLastError();
    if(error!=cudaSuccess){
      std::cout<< "cuda error proyectDataJBN: " << cudaGetErrorString(error) <<std::endl;
      exit(-1);
    }

}

extern "C"{

void cuDivergence(
                 double *rs,
                 double *rsv,
                 double *e,
                 double *f, 
                 double *g,
                 double *ev,
                 double *fv, 
                 double *gv,
                 double *jbn,
                 double deltaX, double deltaY, double deltaZ,
                 int Nx, int Ny, int Nz 
                 ){

   int offset = Nx*Ny*Nz; 


   ArrayDev3D<double> FLUXE(Nx,Ny,Nz);
   ArrayDev3D<double> FLUXF(Nx,Ny,Nz);
   ArrayDev3D<double> FLUXG(Nx,Ny,Nz);

   ArrayDev3D<double> FLUXEv(Nx,Ny,Nz);
   ArrayDev3D<double> FLUXFv(Nx,Ny,Nz);
   ArrayDev3D<double> FLUXGv(Nx,Ny,Nz);

   ArrayDev3D<double> dEdx(Nx,Ny,Nz);
   ArrayDev3D<double> dFdy(Nx,Ny,Nz);
   ArrayDev3D<double> dGdz(Nx,Ny,Nz);

   ArrayDev3D<double> dEvdx(Nx,Ny,Nz);
   ArrayDev3D<double> dFvdy(Nx,Ny,Nz);
   ArrayDev3D<double> dGvdz(Nx,Ny,Nz);

   ArrayDev3D<double> jbn0(Nx,Ny,Nz);
   ArrayDev3D<double> jbn4(Nx,Ny,Nz);
   ArrayDev3D<double> jbn8(Nx,Ny,Nz);

   ArrayDev3D<double> jbndet(Nx,Ny,Nz);
   ArrayDev3D<double> jbn10(Nx,Ny,Nz);

   ArrayDev3D<double> rhs(Nx,Ny,Nz);
   ArrayDev3D<double> rhsv(Nx,Ny,Nz);

   FLUXE.copyFromHost(e);
   FLUXF.copyFromHost(f);
   FLUXG.copyFromHost(g);

   FLUXEv.copyFromHost(ev);
   FLUXFv.copyFromHost(fv);
   FLUXGv.copyFromHost(gv);

   jbn0.copyFromHost(&(jbn[0*offset]));
   jbn4.copyFromHost(&(jbn[4*offset]));
   jbn8.copyFromHost(&(jbn[8*offset]));

   jbndet.copyFromHost(&(jbn[9*offset]));
   jbn10.copyFromHost(&(jbn[10*offset]));


   proyectDataJBN(FLUXE, jbndet, jbn0);
   proyectDataJBN(FLUXF, jbndet, jbn4);
   proyectDataJBN(FLUXG, jbndet, jbn8);

   proyectDataJBN(FLUXEv, jbndet, jbn0);
   proyectDataJBN(FLUXFv, jbndet, jbn4);
   proyectDataJBN(FLUXGv, jbndet, jbn8);

   cu::cs::dFdX(dEdx,  FLUXE ,  deltaX);
   cu::cs::dFdX(dEvdx, FLUXEv , deltaX);

   cu::cs::dFdY(dFdy,  FLUXF ,  deltaY);
   cu::cs::dFdY(dFvdy, FLUXFv , deltaY);

   cu::cs::dFdZ(dGdz,  FLUXG ,  deltaZ);
   cu::cs::dFdZ(dGvdz, FLUXGv , deltaZ);

   getDiv(rhs, jbn10, dEdx, dFdy, dGdz);
   getDiv(rhsv, jbn10, dEvdx, dFvdy, dGvdz);

   rhs.copyToHost(rs);
   rhsv.copyToHost(rsv);

}

}
