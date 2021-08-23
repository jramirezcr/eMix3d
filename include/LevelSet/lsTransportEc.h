#ifndef LEVELSETEQ_H
#define LEVELSETEQ_H

void cuAdvectLS(
               double *d_phi, 
               double *d_vel,
               double  deltaX,
               double  deltaY,
               double  deltaZ,
               double  dt, 
               int     Nx, 
               int     Ny, 
               int     Nz  
               );  

void cuReinitLS(
               double *d_phi,    // Level set function on DEVICE
               double  deltaX,
               double  deltaY,
               double  deltaZ,
               int     Nx, 
               int     Ny, 
               int     Nz  
               );  

void cuAdvectLsJB(
               double *d_phi, 
               double *d_vel,
               double *d_jbn,
               double  deltaX,
               double  deltaY,
               double  deltaZ,
               double  dt, 
               int     Nx, 
               int     Ny, 
               int     Nz  
               );  

void cuReinitLsJB(
               double *d_phi,    // Level set function on DEVICE
               double *d_jbn,
               double  deltaX,
               double  deltaY,
               double  deltaZ,
               double *d_mins,
               int     Nx, 
               int     Ny, 
               int     Nz,
               int     gcells,
               double  deltamin
               );  


#endif

