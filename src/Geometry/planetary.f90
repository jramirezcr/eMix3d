!__________________________________________________________________
       SUBROUTINE planetary()
!__________________________________________________________________

      use dimensiones
      use variables
      use consadim
      use velocidades
      use mallagrid
      use tiempo
      use maskVar

      implicit none
      integer neq, irk, i, j, itc, nmcit, k, no_imp
      real :: ue, ve, ancho, movTheta, Ro, PI
      real, allocatable, dimension(:,:,:,:) ::  xmp
      real, allocatable, dimension(:,:,:,:) ::  xmp1
      real, allocatable, dimension(:,:,:) ::  maskRush
!      thetas=atheta*dto
!      phis=aphi*dto

      PI = 3.1415926d0
    
      ancho = 3.0*dmin
 
      write(6,*) dto

      Ro = 1.0 / (2.0*3.1415926)

      movTheta = (1.0/Ro)*(dto) 

      allocate(xmp(Nx,Ny,Nz,3))
      allocate(xmp1(Nx,Ny,Nz,3))
      allocate(maskRush(Nx,Ny,Nz))

      DO i=1, Nx
         DO j=1, Ny
           DO k=1, Nz
                maskRush(i,j,k) = 1.0
           ENDDO
         ENDDO
      ENDDO

      

      CALL rotateOnAxisZ(xmp,x, -movTheta, Nx, Ny, Nz)
      CALL traslateOnAxisX(xmp1,xmp, 0.2086 , Nx, Ny, Nz)
      CALL rotateOnAxisZ(xmp,xmp1, 1.0*movTheta, Nx, Ny, Nz)
      CALL rotateOnAxisZ(xmp1,xmp,-2.0*movTheta, Nx, Ny, Nz)

      CALL addMaskCube(maskRush,xmp1, -0.1579, -0.09, &
                    -0.5*ancho, 0.5*ancho, 0.0684, 0.93579, Nx, Ny, Nz) 

      CALL addMaskCube(maskRush,xmp1, 0.09, 0.1579, &
                    -0.5*ancho, 0.5*ancho, 0.0684, 0.93579, Nx, Ny, Nz) 

      CALL addMaskCube(maskRush,xmp1, -0.09 ,0.09, &
                    -0.5*ancho, 0.5*ancho, 0.0684, 0.1368, Nx, Ny, Nz) 

      CALL addPlanetaryVel(um1, x, xmp, maskRush, 1.0*Ro, Nx, Ny, Nz)
      
      DO i=1, Nx
         DO j=1, Ny
           DO k=1, Nz
               impellermask(i,j,k) =  maskRush(i,j,k) 
               maskRush(i,j,k) = 1.0
           ENDDO
         ENDDO
      ENDDO
      

      !Disco

      
     !------------------


      movTheta = movTheta + 3.1415926

      CALL rotateOnAxisZ(xmp,x, -movTheta, Nx, Ny, Nz)
      CALL traslateOnAxisX(xmp1,xmp, 0.2086 , Nx, Ny, Nz)
      CALL rotateOnAxisZ(xmp,xmp1, 1.0*movTheta, Nx, Ny, Nz)
      CALL rotateOnAxisZ(xmp1,xmp,-2.0*movTheta - PI*0.5, Nx, Ny, Nz)

      CALL addMaskCube(maskRush,xmp1, -0.1579, -0.09, &
                    -0.5*ancho, 0.5*ancho, 0.0684, 0.93579, Nx, Ny, Nz) 

      CALL addMaskCube(maskRush,xmp1, 0.09, 0.1579, &
                    -0.5*ancho, 0.5*ancho, 0.0684, 0.93579, Nx, Ny, Nz) 

      CALL addMaskCube(maskRush,xmp1, -0.09 ,0.09, &
                    -0.5*ancho, 0.5*ancho, 0.0684, 0.1368, Nx, Ny, Nz) 

      CALL addPlanetaryVel(um1, x, xmp, maskRush, 1.0*Ro, Nx, Ny, Nz)

      DO i=1, Nx
         DO j=1, Ny
           DO k=1, Nz
               impellermask(i,j,k) =  impellermask(i,j,k)*maskRush(i,j,k) 
           ENDDO
         ENDDO
      ENDDO
      

      RETURN
      END SUBROUTINE planetary

