#ifndef EXTRAPOL_KERNEL_H
#define EXTRAPOL_KERNEL_H


__global__ void extrapolKernel(
                               double* const rs,          //RS
                               const double* const extVal,//Var extrapol
                               const double* const phiS,  //Level Set F
                               const double* const jbn,   //Jacobian
                               const double* const d_Phi, //Phi Der
                               const double deltaX,
                               const double deltaY,
                               const double deltaZ,
                               const unsigned int Nx, 
                               const unsigned int Ny, 
                               const unsigned int Nz, 
                               const int extFlag
                              );

__global__ void DevFirstOrder_LS(
                                double* const d_Phi,
                                const double* const phiS,
                                const double* const jbn,
                                const double deltaX,
                                const double deltaY,
                                const double deltaZ,
                                const unsigned int Nx, 
                                const unsigned int Ny, 
                                const unsigned int Nz
                                );

__global__ void RunGK_FirstS(
                            double*  d,
                            double*  d0, 
                            double dt, 
                            double*  rs,
                            const int Nx, const int Ny, const int Nz
                            );

__global__ void RunGK_SecondS(
                             double* const d,
                             double* const d0, 
                             double* const d1, 
                             double dt, 
                             double* const rs,  
                             const int Nx, const int Ny, const int Nz
                            );

__global__ void RunGK_ThirdS(
                             double* d,
                             double* d0,
                             double* d1,
                             double dt,
                             double* rs,
                             const int Nx, const int Ny, const int Nz
                            );
 
	
__global__ void copyLSGas(
                          double* const value,
                          const double* const copyVal, 
                          const double* const phiS,
                          int Nx, int Ny, int Nz
                         );

__global__ void copyLSLiquid(
                             double* const value,
                             const double* const copyVal, 
                             const double* const phiS,
                             int Nx, int Ny, int Nz, 
                             int disp
                            );


void cuExtrapolation(
                    double* extVal_d,
                    double* phiS_d,
                    double* jbn_d,
                    double deltaX, double deltaY, double deltaZ,
                    int Nx, int Ny, int Nz, 
                    double dtext,
                    int Flag
                    );


#endif
