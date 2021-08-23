#include "Mesh/typeDataMSV.hpp"
#include <iostream>

int main(){

   int nx, ny, nz;

   nx = ny = nz = 100;

   double *x = new double[nx*ny*nz];
   double *y = new double[nx*ny*nz];
   double *z = new double[nx*ny*nz];

   double *jbnX = new double[nx*ny*nz];
   double *jbnY = new double[nx*ny*nz];
   double *jbnZ = new double[nx*ny*nz];

   double *jbnDET = new double[nx*ny*nz];

   double *u = new double[nx*ny*nz];
   double *v = new double[nx*ny*nz];
   double *w = new double[nx*ny*nz];

   double *press = new double[nx*ny*nz];
   double *temp = new double[nx*ny*nz];

   MSVMeshData<double> mesh(nx,ny,nz);
   MSVFlowData<double> flow(nx,ny,nz);

   mesh.setMeshValues(x,y,z);
   mesh.setMeshJacobean(jbnX, jbnY, jbnZ, jbnDET);

   flow.setFlowData(u,v,w,press,temp);

   flow.u(10,10,10) = 10.0;

   return 0;
}
