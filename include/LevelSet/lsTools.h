#ifndef _LSTOOLS_H
#define _LSTOOLS_H

__global__ 
void Dev1thO_Downwind(
                      double* const d_Phi,
                      const double* const phiS,
                      const double deltaX,
                      const double deltaY,
                      const double deltaZ,
                      const unsigned int Nx, 
                      const unsigned int Ny, 
                      const unsigned int Nz
                      );  


__global__ 
void PhiDevPlusParameter(
                         double* const       phi_xyz,
                         const double* const d_Phi,
                         unsigned const int  Nx,
                         unsigned const int  Ny,
                         unsigned const int  Nz
                         );

__global__ 
void PhiDevMinusParameter(
                          double* const       phi_xyz,
                          const double* const d_Phi,
                          unsigned const int  Nx,
                          unsigned const int  Ny,
                          unsigned const int  Nz
                          );

__global__ 
void PhiDevPlusParameterJB(
                           double* const       phi_xyz,
                           const double* const d_Phi,
                           const double* const d_jbn,
                           unsigned const int  Nx,
                           unsigned const int  Ny,
                           unsigned const int  Nz
                           ); 

__global__ 
void PhiDevMinusParameterJB(
                          double* const       phi_xyz,
                          const double* const d_Phi,
                          const double* const d_jbn,
                          unsigned const int  Nx,
                          unsigned const int  Ny,
                          unsigned const int  Nz
                          );

__global__ 
void reini_RS_WENOJB(
                   double* const       rs,
                   const double* const phiS,
                   const double* const deltaXYZ,
                   const double* const d_phiP,
                   const double* const d_phiM,
                   const double* const phiS0,
                   unsigned int        Nx,
                   unsigned int        Ny,
                   unsigned int        Nz
                   );

__global__ 
void advect_RS_WENO(
                    double* const       rs,
                    const double* const velocity,
                    const double* const d_phiP_d,
                    const double* const d_phiM_d,
                    unsigned int        Nx, 
                    unsigned int        Ny, 
                    unsigned int        Nz
                    );

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
                          );

__global__
void enrightVelocityProfile(
                            double       *vel,    //Velocity Array
                            double       *xMesh,  //Mesh values 
                            double       *yMesh,
                            double       *zMesh,
                            const int     Nx,     //Mesh dimensions
                            const int     Ny,
                            const int     Nz,
                            const double  time,    //current time
                            const double  period
                           );

__global__
void cuGhostCellsMirror3dZ( 
                          double    *ghostArray,
                          const int  ncells,
                          const int  Nx, 
                          const int  Ny, 
                          const int  Nz, 
                          double     direction
                          ); 

__global__
void cuGhostCellsMirror3dY( 
                          double    *ghostArray,
                          const int  ncells,
                          const int  Nx, 
                          const int  Ny, 
                          const int  Nz,
                          double     direction
                          ); 


__global__
void cuGhostCellsMirror3dX( 
                          double    *ghostArray,
                          const int  ncells,
                          const int  Nx, 
                          const int  Ny, 
                          const int  Nz,
                          double     direction
                          ); 

__global__
void cuSwapToGhost( 
                    double    *ghostArray,
                    double    *valueArray,
                    const int  gcells,
                    const int  Nx, 
                    const int  Ny, 
                    const int  Nz  
                   );

__global__
void cuSwapFromGhost( 
                   double    *valueArray,
                   double    *ghostArray,
                   const int  gcells,
                   const int  Nx, 
                   const int  Ny, 
                   const int  Nz  
                   );

__global__ 
void reini_RS_WENO(
                   double* const       rs,
                   const double* const phiS,
                   const double        deltaXYZ,
                   const double* const d_phiP,
                   const double* const d_phiM,
                   const double* const phiS0,
                   unsigned int        Nx,
                   unsigned int        Ny,
                   unsigned int        Nz
                   );

#endif
