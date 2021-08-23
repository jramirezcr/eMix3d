#include<cuda.h>
#include<stdio.h>

#include"Schemes/derivCompact.hpp"
#include"Engine/cuArray.hpp"

template<typename Tprec>
struct TensorDev2{

   TensorDev2(int i, int j, int k){
      resize(i,j,k);
   }

   void resize(int i, int j, int k){
       e11.resize(i,j,k);
       e21.resize(i,j,k);
       e31.resize(i,j,k);
       e12.resize(i,j,k);
       e22.resize(i,j,k);
       e32.resize(i,j,k);
       e13.resize(i,j,k);
       e23.resize(i,j,k);
       e33.resize(i,j,k);
    
   }

   ArrayDev3D<Tprec> e11, e21, e31; 
   ArrayDev3D<Tprec> e12, e22, e32; 
   ArrayDev3D<Tprec> e13, e23, e33; 

};

template<typename Tprec>
struct PseudoConservative{

   PseudoConservative(int i, int j, int k){
      resize(i,j,k);
   }

   void resize(int i, int j, int k){
       press.resize(i,j,k);
       rhoU.resize(i,j,k);
       rhoV.resize(i,j,k);
       rhoW.resize(i,j,k);
   }

   ArrayDev3D<Tprec> press, rhoU, rhoV, rhoW; 
};


template<typename Tprec>
struct Fluxes{

   Fluxes(int i, int j, int k){
      resize(i,j,k);
   }

   void resize(int i, int j, int k){
       e.resize(i,j,k);
       f.resize(i,j,k);
       g.resize(i,j,k);

       ev.resize(i,j,k);
       fv.resize(i,j,k);
       gv.resize(i,j,k);
   }

   ArrayDev3D<Tprec> e, f, g; 
   ArrayDev3D<Tprec> ev, fv, gv; 
};



template<typename Tprec>
struct VectorDev3D{

   VectorDev3D(int i, int j, int k){
      resize(i,j,k);
   }

   void resize(int i, int j, int k){
       x.resize(i,j,k);
       y.resize(i,j,k);
       z.resize(i,j,k);

   }

   ArrayDev3D<Tprec> x, y, z; 
};



template<typename Tprec>
struct FlowVariables{

   Tprec coefMa; 
   Tprec coefRe; 
   Tprec coefPr; 
   Tprec coefFr; 

};


template<typename Tprec>
__global__
void cuGetGradientJacob(TensorDev2<Tprec> Tensor,
                        ArrayDev3D<Tprec> jbnX,
                        ArrayDev3D<Tprec> jbnY,
                        ArrayDev3D<Tprec> jbnZ
                        ){

     int i = blockDim.x*blockIdx.x + threadIdx.x;
     int j = blockDim.y*blockIdx.y + threadIdx.y;
     int k = blockDim.z*blockIdx.z + threadIdx.z;

     if(i < Tensor.e11.getDim(1) && j < Tensor.e11.getDim(2) 
                                 && k < Tensor.e11.getDim(3)){

        Tensor.e11(i,j,k) = Tensor.e11(i,j,k)*jbnX(i,j,k);
        Tensor.e12(i,j,k) = Tensor.e12(i,j,k)*jbnX(i,j,k);
        Tensor.e13(i,j,k) = Tensor.e13(i,j,k)*jbnX(i,j,k);

        Tensor.e21(i,j,k) = Tensor.e21(i,j,k)*jbnY(i,j,k);
        Tensor.e22(i,j,k) = Tensor.e22(i,j,k)*jbnY(i,j,k);
        Tensor.e23(i,j,k) = Tensor.e23(i,j,k)*jbnY(i,j,k);

        Tensor.e31(i,j,k) = Tensor.e31(i,j,k)*jbnZ(i,j,k);
        Tensor.e32(i,j,k) = Tensor.e32(i,j,k)*jbnZ(i,j,k);
        Tensor.e33(i,j,k) = Tensor.e33(i,j,k)*jbnZ(i,j,k);
     }
}

template<typename Tprec>
__global__
void cuGetViscousSGWALLE(ArrayDev3D<Tprec> vt, 
                         ArrayDev3D<Tprec> Delta, 
                         ArrayDev3D<Tprec> SdSd, 
                         ArrayDev3D<Tprec> SS
                         ){

     int i = blockDim.x*blockIdx.x + threadIdx.x;
     int j = blockDim.y*blockIdx.y + threadIdx.y;
     int k = blockDim.z*blockIdx.z + threadIdx.z;

     if(i < SdSd.getDim(1) && j < SdSd.getDim(2) 
                             && k < SdSd.getDim(3)){

       vt(i,j,k)  = pow(0.5*Delta(i,j,k),2.0)*pow(SdSd(i,j,k),(3.0/2.0))  ;
       vt(i,j,k) /=(pow(SS(i,j,k),(5.0/2.0)) + pow(SdSd(i,j,k),(5.0/4.0)) + 
                    0.00000000000000001);
     }
}

template<typename Tprec>
__global__
void cuGetSdTensorWALLE(TensorDev2<Tprec> Sd,
                        TensorDev2<Tprec> S2,
                        TensorDev2<Tprec> R2,
                        ArrayDev3D<Tprec> SS, 
                        ArrayDev3D<Tprec> RR
                        ){

     int i = blockDim.x*blockIdx.x + threadIdx.x;
     int j = blockDim.y*blockIdx.y + threadIdx.y;
     int k = blockDim.z*blockIdx.z + threadIdx.z;

     if(i < Sd.e11.getDim(1) && j < Sd.e11.getDim(2) 
                             && k < Sd.e11.getDim(3)){

          Sd.e11(i,j,k) = S2.e11(i,j,k) + R2.e11(i,j,k) -
                          (SS(i,j,k) + RR(i,j,k))/3.0;

          Sd.e22(i,j,k) = S2.e22(i,j,k) + R2.e22(i,j,k) -
                          (SS(i,j,k) + RR(i,j,k))/3.0;

          Sd.e33(i,j,k) = S2.e33(i,j,k) + R2.e33(i,j,k) -
                          (SS(i,j,k) + RR(i,j,k))/3.0;

          Sd.e12(i,j,k) = S2.e12(i,j,k) + R2.e12(i,j,k);

          Sd.e13(i,j,k) = S2.e13(i,j,k) + R2.e13(i,j,k);

          Sd.e21(i,j,k) = S2.e21(i,j,k) + R2.e21(i,j,k);

          Sd.e23(i,j,k) = S2.e23(i,j,k) + R2.e23(i,j,k);

          Sd.e31(i,j,k) = S2.e31(i,j,k) + R2.e31(i,j,k);

          Sd.e32(i,j,k) = S2.e32(i,j,k) + R2.e32(i,j,k);

     }
}

template<typename Tprec>
__global__
void cuGetDoubleContraction(ArrayDev3D<Tprec> result, 
                           TensorDev2<Tprec>  Tensor
                           ){

     int i = blockDim.x*blockIdx.x + threadIdx.x;
     int j = blockDim.y*blockIdx.y + threadIdx.y;
     int k = blockDim.z*blockIdx.z + threadIdx.z;

     if(i < Tensor.e11.getDim(1) && j < Tensor.e11.getDim(2) 
                                 && k < Tensor.e11.getDim(3)){

        result(i,j,k) = Tensor.e11(i,j,k)*Tensor.e11(i,j,k) +
                        Tensor.e22(i,j,k)*Tensor.e22(i,j,k) +
                        Tensor.e33(i,j,k)*Tensor.e33(i,j,k) +
                        Tensor.e12(i,j,k)*Tensor.e12(i,j,k) +
                        Tensor.e13(i,j,k)*Tensor.e13(i,j,k) +
                        Tensor.e21(i,j,k)*Tensor.e21(i,j,k) +
                        Tensor.e23(i,j,k)*Tensor.e23(i,j,k) +
                        Tensor.e31(i,j,k)*Tensor.e31(i,j,k) +
                        Tensor.e32(i,j,k)*Tensor.e32(i,j,k);
                       
     }
}
                         

template<typename Tprec>
__global__
void cuGetStrainTensor(TensorDev2<Tprec> Strain,
                       TensorDev2<Tprec> VelocityGradient
                      ){

     int i = blockDim.x*blockIdx.x + threadIdx.x;
     int j = blockDim.y*blockIdx.y + threadIdx.y;
     int k = blockDim.z*blockIdx.z + threadIdx.z;

     if(i < Strain.e11.getDim(1) && j < Strain.e11.getDim(2) 
                                 && k < Strain.e11.getDim(3)){

         Strain.e11(i,j,k) = 2.0*VelocityGradient.e11(i,j,k);
         Strain.e21(i,j,k) = (VelocityGradient.e21(i,j,k) + 
                              VelocityGradient.e12(i,j,k));
         Strain.e31(i,j,k) = (VelocityGradient.e31(i,j,k) + 
                              VelocityGradient.e13(i,j,k));

         Strain.e12(i,j,k) = (VelocityGradient.e12(i,j,k) + 
                              VelocityGradient.e21(i,j,k));
         Strain.e22(i,j,k) = 2.0*VelocityGradient.e22(i,j,k);
         Strain.e32(i,j,k) = (VelocityGradient.e32(i,j,k) + 
                              VelocityGradient.e23(i,j,k));


         Strain.e13(i,j,k) = (VelocityGradient.e13(i,j,k) + 
                              VelocityGradient.e31(i,j,k));
         Strain.e23(i,j,k) = (VelocityGradient.e23(i,j,k) + 
                              VelocityGradient.e32(i,j,k));
         Strain.e33(i,j,k) = 2.0*VelocityGradient.e33(i,j,k);


     }


}


template<typename Tprec>
__global__
void cuGetRotationTensor(TensorDev2<Tprec> Rotation,
                         TensorDev2<Tprec> VelocityGradient
                      ){

     int i = blockDim.x*blockIdx.x + threadIdx.x;
     int j = blockDim.y*blockIdx.y + threadIdx.y;
     int k = blockDim.z*blockIdx.z + threadIdx.z;

     if(i < Rotation.e11.getDim(1) && j < Rotation.e11.getDim(2) 
                                   && k < Rotation.e11.getDim(3)){

         Rotation.e11(i,j,k) = 0.0;
         Rotation.e21(i,j,k) = (VelocityGradient.e21(i,j,k) -
                                VelocityGradient.e12(i,j,k));
         Rotation.e31(i,j,k) = (VelocityGradient.e31(i,j,k) - 
                                VelocityGradient.e13(i,j,k));

         Rotation.e12(i,j,k) = (VelocityGradient.e12(i,j,k) -
                                VelocityGradient.e21(i,j,k));
         Rotation.e22(i,j,k) = 0.0;
         Rotation.e32(i,j,k) = (VelocityGradient.e32(i,j,k) - 
                                VelocityGradient.e23(i,j,k));


         Rotation.e13(i,j,k) = (VelocityGradient.e13(i,j,k) - 
                                VelocityGradient.e31(i,j,k));
         Rotation.e23(i,j,k) = (VelocityGradient.e23(i,j,k) - 
                                VelocityGradient.e32(i,j,k));
         Rotation.e33(i,j,k) = 0.0;


     }
}


template<typename Tprec>
__global__
void cuGetDyadicProduct(TensorDev2<Tprec> T2,
                        TensorDev2<Tprec> Tensor
                      ){

     int i = blockDim.x*blockIdx.x + threadIdx.x;
     int j = blockDim.y*blockIdx.y + threadIdx.y;
     int k = blockDim.z*blockIdx.z + threadIdx.z;

     if(i < T2.e11.getDim(1) && j < T2.e11.getDim(2) 
                                   && k < T2.e11.getDim(3)){
         // i=1, j=1
         T2.e11(i,j,k) = Tensor.e11(i,j,k)*Tensor.e11(i,j,k) +
                         Tensor.e12(i,j,k)*Tensor.e21(i,j,k) +
                         Tensor.e13(i,j,k)*Tensor.e31(i,j,k);

         // i=2, j=2
         T2.e22(i,j,k) = Tensor.e21(i,j,k)*Tensor.e12(i,j,k) +
                         Tensor.e22(i,j,k)*Tensor.e22(i,j,k) +
                         Tensor.e23(i,j,k)*Tensor.e32(i,j,k);

         // i=3, j=3
         T2.e33(i,j,k) = Tensor.e31(i,j,k)*Tensor.e13(i,j,k) +
                         Tensor.e32(i,j,k)*Tensor.e23(i,j,k) +
                         Tensor.e33(i,j,k)*Tensor.e33(i,j,k);

         // i=1, j=2
         T2.e12(i,j,k) = Tensor.e11(i,j,k)*Tensor.e12(i,j,k) +
                         Tensor.e12(i,j,k)*Tensor.e22(i,j,k) +
                         Tensor.e13(i,j,k)*Tensor.e32(i,j,k);

         // i=1, j=3
         T2.e13(i,j,k) = Tensor.e11(i,j,k)*Tensor.e13(i,j,k) +
                         Tensor.e12(i,j,k)*Tensor.e23(i,j,k) +
                         Tensor.e13(i,j,k)*Tensor.e33(i,j,k);

         // i=2, j=1
         T2.e21(i,j,k) = Tensor.e21(i,j,k)*Tensor.e11(i,j,k) +
                         Tensor.e22(i,j,k)*Tensor.e21(i,j,k) +
                         Tensor.e23(i,j,k)*Tensor.e31(i,j,k);

         // i=2, j=3
         T2.e23(i,j,k) = Tensor.e21(i,j,k)*Tensor.e13(i,j,k) +
                         Tensor.e22(i,j,k)*Tensor.e23(i,j,k) +
                         Tensor.e23(i,j,k)*Tensor.e33(i,j,k);

         // i=3, j=1
         T2.e31(i,j,k) = Tensor.e31(i,j,k)*Tensor.e11(i,j,k) +
                         Tensor.e32(i,j,k)*Tensor.e21(i,j,k) +
                         Tensor.e33(i,j,k)*Tensor.e31(i,j,k);

         // i=3, j=2
         T2.e32(i,j,k) = Tensor.e31(i,j,k)*Tensor.e12(i,j,k) +
                         Tensor.e32(i,j,k)*Tensor.e22(i,j,k) +
                         Tensor.e33(i,j,k)*Tensor.e32(i,j,k);
     }

}
///
                       
                        
template<typename Tprec>
void getStrainTensor(TensorDev2<Tprec>& Strain, 
                     TensorDev2<Tprec>& VelocityGradient
                    ){

    int gridDimX, gridDimY, gridDimZ;

    gridDimX = Strain.e11.getDim(1) / 8; 
    if(Strain.e11.getDim(1) % 8) gridDimX++;

    gridDimY = Strain.e11.getDim(2) / 8; 
    if(Strain.e11.getDim(2) % 8) gridDimY++;

    gridDimZ = Strain.e11.getDim(3) / 4; 
    if(Strain.e11.getDim(3) % 4) gridDimZ++;

    dim3 BlockDim(8,8,4);
    dim3 GridDim(gridDimX,gridDimY,gridDimZ);

    cuGetStrainTensor<<<GridDim, BlockDim>>>(Strain, VelocityGradient);

}

template<typename Tprec>
void getRotationTensor(TensorDev2<Tprec>& Rotation, 
                       TensorDev2<Tprec>& VelocityGradient
                    ){

    int gridDimX, gridDimY, gridDimZ;

    gridDimX = Rotation.e11.getDim(1) / 8; 
    if(Rotation.e11.getDim(1) % 8) gridDimX++;

    gridDimY = Rotation.e11.getDim(2) / 8; 
    if(Rotation.e11.getDim(2) % 8) gridDimY++;

    gridDimZ = Rotation.e11.getDim(3) / 4; 
    if(Rotation.e11.getDim(3) % 4) gridDimZ++;

    dim3 BlockDim(8,8,4);
    dim3 GridDim(gridDimX,gridDimY,gridDimZ);

    cuGetRotationTensor<<<GridDim, BlockDim>>>(Rotation, VelocityGradient);

}

template<typename Tprec>
void getDyadicProduct(TensorDev2<Tprec>& T2, 
                      TensorDev2<Tprec>& Tensor
                    ){

    int gridDimX, gridDimY, gridDimZ;

    gridDimX = Tensor.e11.getDim(1) / 8; 
    if(Tensor.e11.getDim(1) % 8) gridDimX++;

    gridDimY = Tensor.e11.getDim(2) / 8; 
    if(Tensor.e11.getDim(2) % 8) gridDimY++;

    gridDimZ = Tensor.e11.getDim(3) / 4; 
    if(Tensor.e11.getDim(3) % 4) gridDimZ++;

    dim3 BlockDim(8,8,4);
    dim3 GridDim(gridDimX,gridDimY,gridDimZ);

    cuGetStrainTensor<<<GridDim, BlockDim>>>(T2, Tensor);

}

template<typename Tprec>
void getDoubleContraction(ArrayDev3D<Tprec>& result, 
                          TensorDev2<Tprec>& Tensor
                         ){

    int gridDimX, gridDimY, gridDimZ;

    gridDimX = Tensor.e11.getDim(1) / 8; 
    if(Tensor.e11.getDim(1) % 8) gridDimX++;

    gridDimY = Tensor.e11.getDim(2) / 8; 
    if(Tensor.e11.getDim(2) % 8) gridDimY++;

    gridDimZ = Tensor.e11.getDim(3) / 4; 
    if(Tensor.e11.getDim(3) % 4) gridDimZ++;

    dim3 BlockDim(8,8,4);
    dim3 GridDim(gridDimX,gridDimY,gridDimZ);

    cuGetDoubleContraction<<<GridDim, BlockDim>>>(result, Tensor);

}

template<typename Tprec>
void getSdTensorWALLE(TensorDev2<Tprec>& Sd,
                      TensorDev2<Tprec>& S2,
                      TensorDev2<Tprec>& R2,
                      ArrayDev3D<Tprec>& SS, 
                      ArrayDev3D<Tprec>& RR
                     ){
    
    int gridDimX, gridDimY, gridDimZ;

    gridDimX = Sd.e11.getDim(1) / 8; 
    if(Sd.e11.getDim(1) % 8) gridDimX++;

    gridDimY = Sd.e11.getDim(2) / 8; 
    if(Sd.e11.getDim(2) % 8) gridDimY++;

    gridDimZ = Sd.e11.getDim(3) / 4; 
    if(Sd.e11.getDim(3) % 4) gridDimZ++;

    dim3 BlockDim(8,8,4);
    dim3 GridDim(gridDimX,gridDimY,gridDimZ);

    cuGetSdTensorWALLE<<<GridDim, BlockDim>>>(Sd, S2, R2, SS, RR);

}

template<typename Tprec>
void getViscousSGWALLE(ArrayDev3D<Tprec>& vt,
                       ArrayDev3D<Tprec>& Delta,
                       ArrayDev3D<Tprec>& SdSd,
                       ArrayDev3D<Tprec>& SS
                      ){

    int gridDimX, gridDimY, gridDimZ;

    gridDimX = SdSd.getDim(1) / 8; 
    if(SdSd.getDim(1) % 8) gridDimX++;

    gridDimY = SdSd.getDim(2) / 8; 
    if(SdSd.getDim(2) % 8) gridDimY++;

    gridDimZ = SdSd.getDim(3) / 4; 
    if(SdSd.getDim(3) % 4) gridDimZ++;

    dim3 BlockDim(8,8,4);
    dim3 GridDim(gridDimX,gridDimY,gridDimZ);

    cuGetViscousSGWALLE<<<GridDim, BlockDim>>>(vt, Delta, SdSd, SS);
     
}

template<typename Tprec>
void getGradientJacob(TensorDev2<Tprec>& Tensor,
                      ArrayDev3D<Tprec>& jbnX,
                      ArrayDev3D<Tprec>& jbnY,
                      ArrayDev3D<Tprec>& jbnZ
                      ){

    int gridDimX, gridDimY, gridDimZ;

    gridDimX = Tensor.e11.getDim(1) / 8; 
    if(Tensor.e11.getDim(1) % 8) gridDimX++;

    gridDimY = Tensor.e11.getDim(2) / 8; 
    if(Tensor.e11.getDim(2) % 8) gridDimY++;

    gridDimZ = Tensor.e11.getDim(3) / 4; 
    if(Tensor.e11.getDim(3) % 4) gridDimZ++;

    dim3 BlockDim(8,8,4);
    dim3 GridDim(gridDimX,gridDimY,gridDimZ);

    cuGetGradientJacob<<<GridDim, BlockDim>>>(Tensor, jbnX, jbnY, jbnZ);

}

extern "C"{

void cudaDeriv(
               double *dVel,
               double *dTemp,
               double *vel,
               double *temp,
               double *amut,
               double *dxyzsgd,
               double *jbn,
               double deltaX, double deltaY, double deltaZ,
               int Nx, int Ny, int Nz
              )
{

     int offset = Nx*Ny*Nz;

     //Cinematic variables on GPU device (Velocity)
     ArrayDev3D<double> u(Nx,Ny,Nz), v(Nx,Ny,Nz), w(Nx,Ny,Nz);

     //Thermodinamic variable on GPU device (Temperture)
     ArrayDev3D<double> T(Nx,Ny,Nz);

     //Gradient by element on GPU device (Velocity)
     ArrayDev3D<double> dudx(Nx,Ny,Nz), dudy(Nx,Ny,Nz), dudz(Nx,Ny,Nz);
     ArrayDev3D<double> dvdx(Nx,Ny,Nz), dvdy(Nx,Ny,Nz), dvdz(Nx,Ny,Nz);
     ArrayDev3D<double> dwdx(Nx,Ny,Nz), dwdy(Nx,Ny,Nz), dwdz(Nx,Ny,Nz);

     //Gradient by element on GPU device (Temperature and Pressure)
     ArrayDev3D<double> dTdx(Nx,Ny,Nz), dTdy(Nx,Ny,Nz), dTdz(Nx,Ny,Nz);
 //    ArrayDev3D<double> sinA(Nx,Ny,Nz), cosA(Nx,Ny,Nz);

     //std::cout<<offset << std::endl; 
/*
     for(int k = 0; k < Nz; k++){
        for(int j = 0; j < Ny; j++){
           for(int i = 0; i < Nx; i++){
              int id = i + j*Nx + k*Nx*Ny;
              temp[id] = sin(2.0*3.1415926*(double)k/((double)(Nz-1)));
           } 
        }
     }

      
     sinA.copyFromHost(temp);

     cu::cs::dFdZ(cosA, sinA, (2.0*3.1415926/((double)(Nz-1))));

     cosA.copyToHost(temp);    
 
     for(int k = 0; k < Nz; k++){
              int j = 1, i = 1;
              int id = i + j*Nx + k*Nx*Ny;
              std::cout<< temp[id]
              <<" "<< temp[id] << std::endl;
     } 

*/
/*
     //Copy from host
*/
     u.copyFromHost(&(vel[0*offset]));
     v.copyFromHost(&(vel[1*offset]));
     w.copyFromHost(&(vel[2*offset]));

     T.copyFromHost(temp);
     
/*
     //Gradients
*/

    //Velocity
      cu::cs::dFdX(dudx, u, deltaX);
      cu::cs::dFdY(dudy, u, deltaY);
      cu::cs::dFdZ(dudz, u, deltaZ);

       cu::cs::dFdX(dvdx, v, deltaX);
      cu::cs::dFdY(dvdy, v, deltaY);
      cu::cs::dFdZ(dvdz, v, deltaZ);

      cu::cs::dFdX(dwdx, w, deltaX);
      cu::cs::dFdY(dwdy, w, deltaY);
      cu::cs::dFdZ(dwdz, w, deltaZ);

     //Temperute
      cu::cs::dFdX(dTdx, T, deltaX);
      cu::cs::dFdY(dTdy, T, deltaY);
      cu::cs::dFdZ(dTdz, T, deltaZ);

     cudaError_t error = cudaGetLastError();
     if(error != cudaSuccess)
     {
        // print the CUDA error message and exit
        printf("CUDA Derivadas error: %s\n", cudaGetErrorString(error));
        exit(-1);
     }


/*
     //Copy values to Host
*/
     //Velocity Gradient without jacobian
     dudx.copyToHost(&(dVel[0*offset]));
     dvdx.copyToHost(&(dVel[1*offset]));
     dwdx.copyToHost(&(dVel[2*offset]));

     dudy.copyToHost(&(dVel[3*offset]));
     dvdy.copyToHost(&(dVel[4*offset]));
     dwdy.copyToHost(&(dVel[5*offset]));

     dudz.copyToHost(&(dVel[6*offset]));
     dvdz.copyToHost(&(dVel[7*offset]));
     dwdz.copyToHost(&(dVel[8*offset]));

     //Temperature Gradient without jacobian
     dTdx.copyToHost(&(dTemp[0*offset]));
     dTdy.copyToHost(&(dTemp[1*offset]));
     dTdz.copyToHost(&(dTemp[2*offset]));

     error = cudaGetLastError();
     if(error != cudaSuccess)
     {
        // print the CUDA error message and exit
        printf("CUDA Todo error: %s\n", cudaGetErrorString(error));
        exit(-1);
     }

}

}
