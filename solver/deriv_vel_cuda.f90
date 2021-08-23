      SUBROUTINE DERIV_VEL()
      use dimensiones
      use velocidades
      use derivvel
      use deltas
      use jacobtools
      use consderper
      use consdernper
      use derivtools
      use sgdmodel
      IMPLICIT NONE
      integer i,j,k,l,m

      real, allocatable, dimension(:,:,:) :: visWALLE


!-----------------------------------------------

      interface
        subroutine cudaDeriv(dvel, dtemp,            &
                             u, temp,                &
                             vt, dxyzsgd,            &
                             jbn,                    &
                             deltax, deltay, deltaz, & 
                             Nx, Ny, Nz)bind(C, name='cudaDeriv')
        use iso_c_binding
        IMPLICIT NONE
        real(c_double) :: dvel(Nx,Ny,Nz,9)
        real(c_double) :: dtemp(Nx,Ny,Nz,3)
        real(c_double) :: u(Nx,Ny,Nz,3)
        real(c_double) :: temp(Nx,Ny,Nz)
        real(c_double) :: vt(Nx,Ny,Nz)
        real(c_double) :: dxyzsgd(Nx,Ny,Nz)
        real(c_double) :: jbn(Nx,Ny,Nz,11)
        real(c_double),value :: deltax, deltay, deltaz
        integer(c_int),value :: Nx, Ny, Nz
        end subroutine
 
      end interface


!----------------------------------------------

        call cudaDeriv(dvel, dtemp,            &
                       u, temp,                &
                       amut, dxyzsgd, jbn, &
                       1.0/ deltax, 1.0 / deltay, 1.0 / deltaz, & 
                       Nx, Ny, Nz)

       write(6,*)'WALE MAX VAL', maxval(amut)

      return
      end subroutine deriv_vel
