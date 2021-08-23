#ifndef __TYPEDATAMSV__
#define __TYPEDATAMSV__

#include "Engine/array.hpp"

template<typename Tprec>
class MSVMeshData{

public:


   MSVMeshData(){
   }

   MSVMeshData(int _nx, int _ny, int _nz ){

     nx = _nx > 0 ? _nx : 1;
     ny = _ny > 0 ? _ny : 1;
     nz = _nz > 0 ? _nz : 1;

     x.resize(_nx, _ny, _nz);  
     y.resize(_nx, _ny, _nz);  
     z.resize(_nx, _ny, _nz);  

     jbnX.resize(_nx, _ny, _nz);
     jbnY.resize(_nx, _ny, _nz);
     jbnZ.resize(_nx, _ny, _nz);
     jbnDET.resize(_nx, _ny, _nz);

   }

   void setMeshValues(Tprec* _x, Tprec* _y, Tprec* _z){
     for(int k = 0; k < nz ; k++){
        for(int j = 0; j < ny ; j++){
           for(int i = 0; i < nx ; i++){

              int id = i + j*nx + k*nx*ny;

              x(i,j,k) = _x[id]; 
              y(i,j,k) = _y[id]; 
              z(i,j,k) = _z[id]; 
           } 
        } 
     } 

     deltaX = 1.0 / ((float)nx - 1.0);
     deltaY = 1.0 / ((float)ny - 1.0);
     deltaZ = 1.0 / ((float)nz - 1.0);

     invDeltaX = 1.0 / deltaX;
     invDeltaY = 1.0 / deltaY;
     invDeltaZ = 1.0 / deltaZ;

   }

   void setMeshJacobean(Tprec* _jbnX, Tprec* _jbnY, Tprec* _jbnZ, Tprec* _jbnDET){
     for(int k = 0; k < nz ; k++){
        for(int j = 0; j < ny ; j++){
           for(int i = 0; i < nx ; i++){

              int id = i + j*nx + k*nx*ny;

              jbnX(i,j,k) = _jbnX[id]; 
              jbnY(i,j,k) = _jbnY[id]; 
              jbnZ(i,j,k) = _jbnZ[id]; 
              jbnX(i,j,k) = _jbnDET[id];

           } 
        } 
     } 

   }

    void setMeshMaxAndMin(Tprec _deltaMin, Tprec _deltaMax){
      deltaMin = _deltaMin;
      deltaMax = _deltaMax;

    }


   
   Array<Tprec, 3> x, y, z;
   Tprec nx, ny, nz;

   Array<Tprec, 3> jbnX, jbnY, jbnZ, jbnDET;

   Tprec deltaX, deltaY, deltaZ;
   Tprec invDeltaX, invDeltaY, invDeltaZ;

   Tprec deltaMin, deltaMax;  
   

};


template<typename Tprec>
class MSVFlowData{

   public:

   MSVFlowData(){
   }

   MSVFlowData(int _nx, int _ny, int _nz ){

     nx = _nx > 0 ? _nx : 1;
     ny = _ny > 0 ? _ny : 1;
     nz = _nz > 0 ? _nz : 1;

     u.resize(_nx, _ny, _nz);  
     v.resize(_nx, _ny, _nz);  
     w.resize(_nx, _ny, _nz);  

     press.resize(_nx, _ny, _nz);
     temp.resize(_nx, _ny, _nz);
   }

   void setFlowData(Tprec* _u, Tprec* _v, Tprec* _w, Tprec* _press, Tprec* _temp){
     for(int k = 0; k < nz ; k++){
        for(int j = 0; j < ny ; j++){
           for(int i = 0; i < nx ; i++){

              int id = i + j*nx + k*nx*ny;

              u(i,j,k) = _u[id]; 
              v(i,j,k) = _v[id]; 
              w(i,j,k) = _w[id]; 

              press(i,j,k) = _press[id];
              temp(i,j,k) = _temp[id];

           } 
        } 
     } 
   }

   Array<Tprec, 3> u, v, w;
   Array<Tprec, 3> press, temp;
   Tprec nx, ny, nz;

};   

#endif
