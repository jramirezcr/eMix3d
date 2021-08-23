!__________________________________________________________________
       SUBROUTINE rushton()
!__________________________________________________________________

      use dimensiones
      use variables
      use consadim
      use velocidades
      use mallagrid
      use tiempo

      implicit none
      integer neq, irk, i, j, itc, nmcit, k, no_imp
      real :: ue, ve, ancho, movTheta, Ro
      real, allocatable, dimension(:,:,:,:) ::  xmp
      real, allocatable, dimension(:,:,:) ::  maskRush, maskDye
!      thetas=atheta*dto
!      phis=aphi*dto
    
      ancho = 2.d0*dmin
 
      write(6,*) dto

      Ro = 1.0d0 / (2.0d0*acos(-1.0d0))

      movTheta = (1.0d0/Ro)*(dto) 

      allocate(xmp(Nx,Ny,Nz,3))
      allocate(maskRush(Nx,Ny,Nz))
      allocate(maskDye(Nx,Ny,Nz))

      DO i=1, Nx
         DO j=1, Ny
           DO k=1, Nz
                maskRush(i,j,k) = 1.d0
                maskDye(i,j,k) = 1.d0

           ENDDO
         ENDDO
      ENDDO

      !Disco

      CALL addMAskCil(maskRush, x, 0.0d0, 0.375d0, (2.d0/3.d0 - ancho*0.5d0), &
                      (2.d0/3.d0 + ancho*0.5d0), Nx, Ny, Nz)

      CALL addMAskCil(maskRush, x, 0.0d0, 0.1d0, 2.d0/3.d0, &
                      2.d0, Nx, Ny, Nz)
!
!     CALL addMAskCil(maskRush, x, 0.0, 0.1891, &
!                     (1.0 + ancho*0.5) , &
!                     (1.0 + ancho*0.5) + 0.27, Nx, Ny, Nz)
!
!     CALL addMAskCil(maskRush, x, 0.0, 0.1351, &
!                     (1.0 + ancho*0.5) + 0.27, &
!                     (1.0 + ancho*0.5) + 0.75, Nx, Ny, Nz)

      !Primera rotacion que dara movimiento al primer grupo de 
      !paletas

      CALL rotateOnAxisZ(xmp,x, movTheta, Nx, Ny, Nz)


      CALL addMaskCube(maskRush,xmp, -1.d0*ancho*0.5d0,              &
                       1.d0*ancho*0.5d0,&
                       0.25d0, 0.5d0, 2.d0/3.d0 - 0.1d0, 2.d0/3.d0 + 0.1d0, Nx, Ny, Nz) 
      
      CALL addMaskCube(maskRush,xmp, -1.d0*ancho*0.5,              &
                       1.d0*ancho*0.5d0,&
                     -0.5d0, -0.25d0, 2.d0/3.d0 - 0.1d0, 2.d0/3.d0 + 0.1d0, Nx, Ny, Nz) 


     !------------------

      CALL rotateOnAxisZ(xmp,x, movTheta + 1.04719d0,  Nx, Ny, Nz)


      CALL addMaskCube(maskRush,xmp, -1.d0*ancho*0.5d0,              &
                       1.d0*ancho*0.5,&
                       0.25d0, 0.5d0, 2.d0/3.d0 - 0.1d0, 2.d0/3.d0 + 0.1d0, Nx, Ny, Nz) 
      
      CALL addMaskCube(maskRush,xmp,  -1.0*ancho*0.5,              &
                       1.0*ancho*0.5,&
                      -0.5, -0.25, 2.d0/3.d0 - 0.1, 2.d0/3.d0 + 0.1, Nx, Ny, Nz) 


      !-----

      CALL rotateOnAxisZ(xmp,x, movTheta +2.0*1.04719,  Nx, Ny, Nz)


      CALL addMaskCube(maskRush,xmp, -1.0*ancho*0.5,              &
                       1.0*ancho*0.5,&
                       0.25, 0.5, 2.d0/3.d0 - 0.1, 2.d0/3.d0 + 0.1, Nx, Ny, Nz) 
      
      CALL addMaskCube(maskRush,xmp, -1.0*ancho*0.5,              &
                       1.0*ancho*0.5,&
                      -0.5, -0.25, 2.d0/3.d0 - 0.1, 2.d0/3.d0 + 0.1, Nx, Ny, Nz) 

       CALL addAngularVel(um1, x, maskRush, -1.0*Ro, Nx, Ny, Nz)
      

      RETURN
      END SUBROUTINE rushton

