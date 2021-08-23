      SUBROUTINE rotateOnAxisZ(meshp, mesh, theta, Nx, Ny, Nz)
!------------------------------------------------------------
!     Jorge Ramirez Cruz
!     This subroutine rotates the mesh entry in a theta angle
!
!-----Parameters---------------------------------------------
      IMPLICIT NONE
      real    :: meshp(Nx,Ny,Nz,3)    !Array result
      real    :: mesh(Nx,Ny,Nz,3)     !Array to rotate
      real    :: theta                !Angle 
      integer :: Nx                   !Grid dimension axis X
      integer :: Ny                   !Grid dimension axis Y    
      integer :: Nz                   !Grid dimension axis z
!------------------------------------------------------------
      integer :: i, j, k

     

      DO k = 1, Nz
         DO j = 1, Ny
            DO i = 1, Nx
                meshp(i,j,k,1) = mesh(i,j,k,1)*cos(theta) -            &
                                 mesh(i,j,k,2)*sin(theta)

                meshp(i,j,k,2) = mesh(i,j,k,1)*sin(theta) +            &
                                 mesh(i,j,k,2)*cos(theta)

                meshp(i,j,k,3) = mesh(i,j,k,3) 
               
            ENDDO
         ENDDO
      ENDDO
      
      END SUBROUTINE rotateOnAxisZ

      SUBROUTINE traslateOnAxisX(meshp, mesh, value, Nx, Ny, Nz)
!------------------------------------------------------------
!     Jorge Ramirez Cruz
!     This subroutine translate the mesh entry some offset value
!
!-----Parameters---------------------------------------------
      IMPLICIT NONE
      real    :: meshp(Nx,Ny,Nz,3)    !Array result
      real    :: mesh(Nx,Ny,Nz,3)     !Array to rotate
      real    :: value                !offset value 
      integer :: Nx                   !Grid dimension axis X
      integer :: Ny                   !Grid dimension axis Y    
      integer :: Nz                   !Grid dimension axis z
!------------------------------------------------------------
      integer :: i, j, k

     

      DO k = 1, Nz
         DO j = 1, Ny
            DO i = 1, Nx
                meshp(i,j,k,1) = mesh(i,j,k,1) - value

                meshp(i,j,k,2) = mesh(i,j,k,2)

                meshp(i,j,k,3) = mesh(i,j,k,3) 
               
            ENDDO
         ENDDO
      ENDDO
      
      END SUBROUTINE traslateOnAxisX


!--------------------------------------------------------------


      SUBROUTINE addMaskCube(maskVar, mesh, Lx1, Lx2, Ly1, Ly2,        &
                                            Lz1, Lz2, Nx, Ny, Nz)
!------------------------------------------------------------
!     Jorge Ramirez Cruz
!     This subroutine creates a cube mas a cube maskk 
!
!-----Parameters---------------------------------------------
      IMPLICIT NONE
      real :: maskVar(Nx,Ny,Nz)       !Mask variable
      real :: mesh(Nx,Ny,Nz,3)        !Mesh dimentions
      real :: Lx1, Lx2                !Space X delimeters
      real :: Ly1, Ly2                !Space Y delimeters
      real :: Lz1, Lz2                !Space Z delimeters
      integer :: Nx                   !Grid dimension axis X
      integer :: Ny                   !Grid dimension axis Y    
      integer :: Nz                   !Grid dimension axis z
!------------------------------------------------------------
      integer :: i, j, k
      
      
      DO k = 1, Nz
         DO j = 1, Ny
            DO i = 1, Nx
                maskVar(i,j,k) = maskVar(i,j,k)*merge(0.d0, 1.d0,      &
                                         ((mesh(i,j,k,1).GE.Lx1).AND.  &
                                          (mesh(i,j,k,1).LE.Lx2).AND.  &
                                          (mesh(i,j,k,2).GE.Ly1).AND.  &
                                          (mesh(i,j,k,2).LE.Ly2).AND.  &
                                          (mesh(i,j,k,3).GE.Lz1).AND.  &
                                          (mesh(i,j,k,3).LE.Lz2))      &
                                               )
            ENDDO
         ENDDO
      ENDDO
    
      END SUBROUTINE addMaskCube

      SUBROUTINE addMaskCil(maskVar, mesh, r1, r2, Lz1, Lz2, Nx, Ny, Nz)
!------------------------------------------------------------
!     Jorge Ramirez Cruz
!     This subroutine creates a cilynder mask 
!
!-----Parameters---------------------------------------------
      IMPLICIT NONE
      real :: maskVar(Nx,Ny,Nz)       !Mask variable
      real :: mesh(Nx,Ny,Nz,3)        !Mesh dimentions
      real :: r1, r2                  !Space r delimeters
      real :: Lz1, Lz2                !Space Z delimeters
      integer :: Nx                   !Grid dimension axis X
      integer :: Ny                   !Grid dimension axis Y    
      integer :: Nz                   !Grid dimension axis z
!------------------------------------------------------------
      integer :: i, j, k
      real    :: radio

      DO k = 1, Nz
         DO j = 1, Ny
            DO i = 1, Nx
                radio = sqrt(mesh(i,j,k,1)*mesh(i,j,k,1) +             &
                            mesh(i,j,k,2)*mesh(i,j,k,2))
                maskVar(i,j,k) = maskVar(i,j,k)*merge(0.d0, 1.d0,      &
                                 ((radio.GE.r1).AND.(radio.LE.r2).AND. &
                                 (mesh(i,j,k,3).GE.Lz1).AND.           &
                                 (mesh(i,j,k,3).LE.Lz2))) 

            ENDDO
         ENDDO
      ENDDO

       
      END SUBROUTINE addMaskCil


      SUBROUTINE addAngularVel(um, mesh, mask, Ro, Nx, Ny, Nz)
!------------------------------------------------------------
!     Jorge Ramirez Cruz
!
!-----Parameters---------------------------------------------
      IMPLICIT NONE
      real    :: um(Nx,Ny,Nz,10)         !Conservative variables
      real    :: mesh(Nx,Ny,Nz,3)        !Mesh dimentions
      real    :: mask(Nx,Ny,Nz)          !Mask variable
      real    :: Ro                      !Rossby number
      integer :: Nx                      !Grid dimension axis X
      integer :: Ny                      !Grid dimension axis Y    
      integer :: Nz                      !Grid dimension axis z
!------------------------------------------------------------
      integer :: i, j, k
      real    :: radio, thetas, ue, ve

      Ro = 1.0 / Ro

      DO k = 1, Nz
         DO j = 1, Ny
            DO i = 1, Nx
               radio = sqrt(mesh(i,j,k,1)*mesh(i,j,k,1) +             &
                            mesh(i,j,k,2)*mesh(i,j,k,2))

               IF(mesh(i,j,k,1).NE.0.0)then
                 thetas = atan(abs(mesh(i,j,k,2)/mesh(i,j,k,1)))
               ELSEIF(mesh(i,j,k,1).EQ.0.0)then
                 thetas = 3.1415926/2.
               ENDIF 

               IF((mesh(i,j,k,1).LT.0.0).AND.(mesh(i,j,k,2).GT.0.0))then
                 thetas = 3.141592 - thetas
               ELSEIF((mesh(i,j,k,1).LT.0.0).AND.                      &
                     (mesh(i,j,k,2).LE.0.0))then
                 thetas = thetas + 3.141592
               ELSEIF((mesh(i,j,k,1).GE.0.0).AND.                      &
                     (mesh(i,j,k,2).LT.0.0))then
                 thetas = 2.0*3.141592 - thetas
               ENDIF 

               ue =-1.0*Ro*radio*sin(thetas)
               ve = Ro*radio*cos(thetas)

               um(i,j,k,2) = ((1.0 - mask(i,j,k))*ue                   &
                           + mask(i,j,k)*(um(i,j,k,2)/um(i,j,k,1)))*   &
                             um(i,j,k,1) 

               um(i,j,k,3) = ((1.0 - mask(i,j,k))*ve                   &
                           + mask(i,j,k)*(um(i,j,k,3)/um(i,j,k,1)))*   &
                             um(i,j,k,1)

               um(i,j,k,4) = mask(i,j,k)*um(i,j,k,4)


            ENDDO
         ENDDO
      ENDDO

      END SUBROUTINE addAngularVel

      SUBROUTINE addPlanetaryVel(um, mesh1, mesh2, mask, Ro, Nx, Ny, Nz)
!------------------------------------------------------------
!     Jorge Ramirez Cruz
!
!-----Parameters---------------------------------------------
      IMPLICIT NONE
      real    :: um(Nx,Ny,Nz,10)         !Conservative variables
      real    :: mesh1(Nx,Ny,Nz,3)       !Mesh1 dimentions
      real    :: mesh2(Nx,Ny,Nz,3)       !Mesh2 dimentions
      real    :: mask(Nx,Ny,Nz)          !Mask variable
      real    :: Ro                      !Rossby number
      integer :: Nx                      !Grid dimension axis X
      integer :: Ny                      !Grid dimension axis Y    
      integer :: Nz                      !Grid dimension axis z
!------------------------------------------------------------
      integer :: i, j, k
      real    :: radio, thetas, ue, ve, ue2, ve2

      Ro = 1.0 / Ro

      DO k = 1, Nz
         DO j = 1, Ny
            DO i = 1, Nx
               radio = sqrt(mesh1(i,j,k,1)*mesh1(i,j,k,1) +             &
                            mesh1(i,j,k,2)*mesh1(i,j,k,2))

               IF(mesh1(i,j,k,1).NE.0.0)then
                 thetas = atan(abs(mesh1(i,j,k,2)/mesh1(i,j,k,1)))
               ELSEIF(mesh1(i,j,k,1).EQ.0.0)then
                 thetas = 3.1415926/2.
               ENDIF 

               IF((mesh1(i,j,k,1).LT.0.0).AND.(mesh1(i,j,k,2).GT.0.0))then
                 thetas = 3.141592 - thetas
               ELSEIF((mesh1(i,j,k,1).LT.0.0).AND.                      &
                     (mesh1(i,j,k,2).LE.0.0))then
                 thetas = thetas + 3.141592
               ELSEIF((mesh1(i,j,k,1).GE.0.0).AND.                      &
                     (mesh1(i,j,k,2).LT.0.0))then
                 thetas = 2.0*3.141592 - thetas
               ENDIF 

               ue =-1.0*Ro*radio*sin(thetas)
               ve = Ro*radio*cos(thetas)



               radio = sqrt(mesh2(i,j,k,1)*mesh2(i,j,k,1) +             &
                            mesh2(i,j,k,2)*mesh2(i,j,k,2))

               IF(mesh2(i,j,k,1).NE.0.0)then
                 thetas = atan(abs(mesh2(i,j,k,2)/mesh2(i,j,k,1)))
               ELSEIF(mesh2(i,j,k,1).EQ.0.0)then
                 thetas = 3.1415926/2.
               ENDIF 

               IF((mesh2(i,j,k,1).LT.0.0).AND.(mesh2(i,j,k,2).GT.0.0))then
                 thetas = 3.141592 - thetas
               ELSEIF((mesh2(i,j,k,1).LT.0.0).AND.                      &
                     (mesh2(i,j,k,2).LE.0.0))then
                 thetas = thetas + 3.141592
               ELSEIF((mesh2(i,j,k,1).GE.0.0).AND.                      &
                     (mesh2(i,j,k,2).LT.0.0))then
                 thetas = 2.0*3.141592 - thetas
               ENDIF 

               ue2 =-2.0*Ro*radio*sin(thetas)
               ve2 = 2.0*Ro*radio*cos(thetas)


               um(i,j,k,2) = (1.0 - mask(i,j,k))*ue                    &
                           + (1.0 - mask(i,j,k))*ue2                   &
                           + mask(i,j,k)*um(i,j,k,2)

               um(i,j,k,3) = (1.0 - mask(i,j,k))*ve                    &
                           + (1.0 - mask(i,j,k))*ve2                   &
                           + mask(i,j,k)*um(i,j,k,3)

               um(i,j,k,4) = mask(i,j,k)*um(i,j,k,4)

            ENDDO
         ENDDO
      ENDDO

      END SUBROUTINE addPlanetaryVel


