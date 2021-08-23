!-----Interfaz to an external C function

      interface
         subroutine advectLS_CUDA(pres, jbn, sobject,    &
                  deltax, deltay, deltaz, Nx, Ny, Nz,    &
                  deltamin)bind(C, name = 'advect_CUDA')
         use iso_c_binding
         IMPLICIT NONE
         real(c_double) :: pres(Nx, Ny, Nz)
         real(c_double) :: jbn(Nx, Ny, Nz,11)
         real(c_double) :: sobject(Nx, Ny, Nz)
         real(c_double), value :: deltax, deltay, deltaz 
         integer(c_int), value :: Nx, Ny, Nz
         real(c_double), value :: deltamin
         end subroutine 
      end interface


!______________________________________________________________________
