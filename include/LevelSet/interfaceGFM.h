
!-----Interface to an external C function-------------------------------

      interface
          subroutine callCUDA(pressg, u,phiS, jbn, deltaX, &
                     deltaY, deltaZ,Nx,Ny,Nz, dtext ) &
                     bind (C,name='callCUDA')
          use iso_c_binding    
          IMPLICIT NONE
          real(c_double) :: pressg(Nx,Ny,Nz)
          real(c_double) :: u(Nx,Ny,Nz,3)
          real(c_double) :: phiS(Nx,Ny,Nz)
          real(c_double) :: jbn(Nx,Ny,Nz,11)
          real(c_double),value :: deltaX, deltaY, deltaZ
          integer(c_int),value :: Nx, Ny, Nz
          real(c_double),value :: dtext
    
          end subroutine callCUDA
    
      end interface

      interface
          subroutine extrapolVarCUDA(press, phiS, d_phi, &
                     jbn, deltaX, &
                     deltaY, deltaZ,Nx,Ny,Nz, dtext ) &
                     bind (C,name='extrapolVarCUDA')
          use iso_c_binding    
          IMPLICIT NONE
          real(c_double) :: press(Nx,Ny,Nz)
          real(c_double) :: phiS(Nx,Ny,Nz)
          real(c_double) :: d_phi(Nx,Ny,Nz,3)
          real(c_double) :: jbn(Nx,Ny,Nz,11)
          real(c_double),value :: deltaX, deltaY, deltaZ
          integer(c_int),value :: Nx, Ny, Nz
          real(c_double),value :: dtext
    
          end subroutine extrapolVarCUDA
      end interface
    
