      !-----Iterface to an external C
      !function--------------------------------

      interface
         subroutine reiniz_CUDA(s, jbn, deltaxyz, deltax, deltay, &
         deltaz, Nx, Ny, Nz, dtin)bind(C, name='reiniz_CUDA')
         use iso_c_binding
         IMPLICIT NONE    
         real(c_double) :: s(Nx, Ny, Nz) 
         real(c_double) :: jbn(Nx, Ny, Nz,11)
         real(c_double) :: deltaxyz(Nx, Ny, Nz) 
         real(c_double), value :: deltax, deltay, deltaz  
         integer(c_int), value :: Nx, Ny, Nz
         real(c_double), value :: dtin 
         end subroutine
      end interface
 
!-----

