!__________________________________________________________________
       SUBROUTINE rushton()
!__________________________________________________________________

      use dimensiones
      use variables
      use consadim
      use velocidades
      use mallagrid
      use tiempo
      use maskvar

      implicit none
      integer neq, irk, i, j, itc, nmcit, k, no_imp
      real :: ue, ve,  theta, Ro
      real, allocatable, dimension(:,:,:,:) ::  xmp
      real, allocatable, dimension(:,:,:) ::  maskRush, maskDye

      integer,allocatable, dimension(:,:):: cflag
      real,allocatable, dimension(:,:)   :: point1_x, point1_y, point1_z
      real,allocatable, dimension(:,:)   :: point2_x, point2_y, point2_z
      real,allocatable, dimension(:,:)   :: point3_x, point3_y, point3_z
      real,allocatable, dimension(:,:)   :: n_x, n_y
      real, allocatable,dimension(:,:,:)  :: matrix, matrixpress
      integer,allocatable, dimension(:,:) :: index1i, index1j
      integer,allocatable, dimension(:,:) :: index2i, index2j
      real :: ux, vy, ancho
      real coefB1, coefB2, coefB3, pi
      integer :: it_im

      allocate(cflag(nx,ny))
      allocate(point1_x(nx,ny))
      allocate(point1_y(nx,ny))
      allocate(point1_z(nx,ny))
      allocate(point2_x(nx,ny))
      allocate(point2_y(nx,ny))
      allocate(point2_z(nx,ny))
      allocate(index1i(nx,ny))
      allocate(index1j(nx,ny))
      allocate(index2i(nx,ny))
      allocate(index2j(nx,ny))
      allocate(point3_x(nx,ny))
      allocate(point3_y(nx,ny))
      allocate(point3_z(nx,ny))
      allocate(matrix(nx,ny,9))
      allocate(matrixpress(nx,ny,9))
      allocate(n_x(nx,ny))
      allocate(n_y(nx,ny))


      ancho = 1.3d0*dmin
 
      write(6,*) dto

      Ro = 1.0d0 / (2.0d0*acos(-1.0d0)) 
      pi = 2.d0*acos(0.d0)

      theta = -((1.0d0/Ro)*(dto) + pi/6.d0)

      allocate(xmp(Nx,Ny,Nz,3))
      allocate(maskRush(Nx,Ny,Nz))
      allocate(maskDye(Nx,Ny,Nz))

      DO i=1, Nx
         DO j=1, Ny
           DO k=1, Nz
                maskRush(i,j,k) = 1.d0
                maskDye(i,j,k)  = 1.d0

                cflag(i,j)    = 1

           ENDDO
         ENDDO
      ENDDO

      !ancho = 0.0085

      DO it_im = 1, 6
         !R1, R2, Z1, Z2, tickness
         call makeAspa(x, maskRush, theta, 0.25, 0.5, 1.5d0 - 0.1, & 
                        1.5d0 + 0.1,ancho, nx, ny, nz)

         call makeAspa(x, maskRush, theta, 0.25, 0.5, 2.5d0 - 0.1, & 
                        2.5d0 + 0.1,ancho, nx, ny, nz)

         call makeAspa(x, maskRush, theta, 0.25, 0.5, 3.5d0 - 0.1, & 
                        3.5d0 + 0.1,ancho, nx, ny, nz)

         theta = theta + pi/3.d0 
      ENDDO

      CALL addMAskCil(maskRush, x, 0.0, 0.375, (1.5d0 - 1.0*ancho), &
                     (1.5d0 + 1.0*ancho), Nx, Ny, Nz) 

      CALL addMAskCil(maskRush, x, 0.0, 0.375, (2.5d0 - 1.0*ancho), &
                     (2.5d0 + 1.0*ancho), Nx, Ny, Nz) 

      CALL addMAskCil(maskRush, x, 0.0, 0.375, (3.5d0 - 1.0*ancho), &
                     (3.5d0 + 1.0*ancho), Nx, Ny, Nz) 


      CALL addMAskCil(maskRush, x, 0.0d0, 0.07895d0, 1.5d0, &
                      5.d0, Nx, Ny, Nz) 

      CALL addAngularVel(um1, x, maskRush, -1.0*Ro, Nx, Ny, Nz)

!----------------------------------------------------------------------------
      DO it_im = 1, 6
 
      call getInterpolData(x, matrix, matrixpress,               &
                           point1_x, point1_y,                   &
                           point2_x, point2_y,point3_x, point3_y,&
                           n_x, n_y,                             &
                           index1i, index1j,                     &
                           index2i, index2j,                     &
                           cflag, theta, 0.25, 0.5, ancho, nx, ny, nz)

     
      DO k =1,nz
         DO j =1,ny
            DO i =1,nx
               if((cflag(i,j).EQ.-1).AND.(&
       ((x(i,j,k,3).GE.(1.5d0 - 0.1)).AND.(x(i,j,k,3).LE.(1.5d0 + 0.1))).OR.&
       ((x(i,j,k,3).GE.(2.5d0 - 0.1)).AND.(x(i,j,k,3).LE.(2.5d0 + 0.1))).OR.&
       ((x(i,j,k,3).GE.(3.5d0 - 0.1)).AND.(x(i,j,k,3).LE.(3.5d0 + 0.1)))&
                 )&
                 )then

                  coefB1 = matrix(i,j,1)*0.0 &
                         + matrix(i,j,2)*maskRush(index1i(i,j),index1j(i,j),k) &
                         + matrix(i,j,3)*maskRush(index2i(i,j),index2j(i,j),k) 

                  coefB2 = matrix(i,j,4)*0.0 &
                         + matrix(i,j,5)*maskRush(index1i(i,j),index1j(i,j),k) &
                         + matrix(i,j,6)*maskRush(index2i(i,j),index2j(i,j),k)

                  coefB3 = matrix(i,j,7)*0.0 &
                         + matrix(i,j,8)*maskRush(index1i(i,j),index1j(i,j),k) &
                         + matrix(i,j,9)*maskRush(index2i(i,j),index2j(i,j),k) 

                  maskRush(i,j,k) = coefB1 + x(i,j,k,1)*coefB2 + x(i,j,k,2)*coefB3

                  call addAngularVelPoint(ux, vy, point1_x(i,j),point1_y(i,j), -Ro)

                  !Extrapol velocidad X
                  coefB1 = matrix(i,j,1)*ux &
                  + matrix(i,j,2)*(um1(index1i(i,j),index1j(i,j),k,2)/  &
                                   um1(index1i(i,j),index1j(i,j),k,1))  &
                  + matrix(i,j,3)*(um1(index2i(i,j),index2j(i,j),k,2)/  &
                                   um1(index2i(i,j),index2j(i,j),k,1))

                  coefB2 = matrix(i,j,4)*ux &
                  + matrix(i,j,5)*(um1(index1i(i,j),index1j(i,j),k,2)/  &
                                   um1(index1i(i,j),index1j(i,j),k,1))  &
                  + matrix(i,j,6)*(um1(index2i(i,j),index2j(i,j),k,2)/  &
                                   um1(index2i(i,j),index2j(i,j),k,1))


                  coefB3 = matrix(i,j,7)*ux &
                  + matrix(i,j,8)*(um1(index1i(i,j),index1j(i,j),k,2)/  &
                                   um1(index1i(i,j),index1j(i,j),k,1))  &
                  + matrix(i,j,9)*(um1(index2i(i,j),index2j(i,j),k,2)/  &
                                   um1(index2i(i,j),index2j(i,j),k,1))

                  um1(i,j,k,2) = (coefB1+x(i,j,k,1)*coefB2+x(i,j,k,2)*coefB3)&
                                 *um1(i,j,k,1)

                  !Extrapol velocidad Y

                  coefB1 = matrix(i,j,1)*vy &
                  + matrix(i,j,2)*(um1(index1i(i,j),index1j(i,j),k,3)/  &
                                   um1(index1i(i,j),index1j(i,j),k,1))  &
                  + matrix(i,j,3)*(um1(index2i(i,j),index2j(i,j),k,3)/  &
                                   um1(index2i(i,j),index2j(i,j),k,1))

                  coefB2 = matrix(i,j,4)*vy &
                  + matrix(i,j,5)*(um1(index1i(i,j),index1j(i,j),k,3)/  &
                                   um1(index1i(i,j),index1j(i,j),k,1))  &
                  + matrix(i,j,6)*(um1(index2i(i,j),index2j(i,j),k,3)/  &
                                   um1(index2i(i,j),index2j(i,j),k,1))

                  coefB3 = matrix(i,j,7)*vy &
                  + matrix(i,j,8)*(um1(index1i(i,j),index1j(i,j),k,3)/  &
                                   um1(index1i(i,j),index1j(i,j),k,1))  &
                  + matrix(i,j,9)*(um1(index2i(i,j),index2j(i,j),k,3)/  &
                                   um1(index2i(i,j),index2j(i,j),k,1))


                  um1(i,j,k,3) =(coefB1+x(i,j,k,1)*coefB2+x(i,j,k,2)*coefB3) &
                                *um1(i,j,k,1)

                   
               endif
            ENDDO
         ENDDO
      ENDDO

      theta = theta + pi/3.d0 
      ENDDO


      DO i=1, Nx
         DO j=1, Ny
           DO k=1, Nz
                impellermask(i,j,k) =  maskRush(i,j,k) 
           ENDDO
         ENDDO
      ENDDO

      END SUBROUTINE rushton


      SUBROUTINE addAngularVelPoint(ux, vy, pointx, pointy, Ro)
!------------------------------------------------------------
!     Jorge Ramirez Cruz
!
!-----Parameters---------------------------------------------
      IMPLICIT NONE
      real    :: ux                      
      real    :: vy                      
      real    :: Ro
      real    :: pointx                  
      real    :: pointy
!------------------------------------------------------------
      real    :: radio, thetas, ue, ve

      Ro = 1.d0 / Ro

      radio = sqrt(pointx*pointx +   &
                   pointy*pointy)

      IF(pointx.NE.0.d0)then
        thetas = atan(abs(pointy/pointx))
      ELSEIF(pointx.EQ.0.d0)then
        thetas = 3.1415926d0/2.d0
      ENDIF

      IF((pointx.LT.0.d0).AND.(pointy.GT.0.0))then
        thetas = 3.141592d0 - thetas
      ELSEIF((pointx.LT.0.d0).AND.                      &
            (pointy.LE.0.d0))then
        thetas = thetas + 3.141592d0
      ELSEIF((pointx.GE.0.d0).AND.                      &
            (pointy.LT.0.d0))then
        thetas = 2.d0*3.141592d0 - thetas
      ENDIF

      ux =-1.d0*Ro*radio*sin(thetas)
      vy = Ro*radio*cos(thetas)

      END SUBROUTINE addAngularVelPoint

      SUBROUTINE getMatrix(matrix,                      &
                           point1_x, point1_y,          &
                           point2_x, point2_y,          &
                           point3_x, point3_y,          &
                           cflag,                       &
                           nx,ny)
        IMPLICIT NONE
        real    :: matrix(nx, ny, 9)
        real    :: point1_x(nx,ny),  point1_y(nx,ny)
        real    :: point2_x(nx,ny),  point2_y(nx,ny)
        real    :: point3_x(nx,ny),  point3_y(nx,ny)
        integer :: nx, ny
        real    :: deter
        integer :: i,j,k
        integer :: cflag(nx,ny)

        DO j =2,ny-1
           DO i =2,nx-1

           IF(cflag(i,j).EQ.-1)then
        
                deter = point2_x(i,j)*point3_y(i,j) &
                      - point3_x(i,j)*point2_y(i,j) &
                      + point1_x(i,j)*point2_y(i,j) &
                      - point1_x(i,j)*point3_y(i,j) &
                      + point3_x(i,j)*point1_y(i,j) &
                      - point2_x(i,j)*point1_y(i,j) 


                matrix(i,j,1) = (point2_x(i,j)*point3_y(i,j)           &
                              -  point3_x(i,j)*point2_y(i,j))/deter

                matrix(i,j,2) = (point3_x(i,j)*point1_y(i,j)           &
                              -  point1_x(i,j)*point3_y(i,j))/deter

                matrix(i,j,3) = (point1_x(i,j)*point2_y(i,j)          &
                              -  point2_x(i,j)*point1_y(i,j))/deter
         
         
                matrix(i,j,4) = (point2_y(i,j)-point3_y(i,j))/deter

                matrix(i,j,5) = (point3_y(i,j)-point1_y(i,j))/deter

                matrix(i,j,6) = (point1_y(i,j)-point2_y(i,j))/deter
         
         
                matrix(i,j,7) = (point3_x(i,j)-point2_x(i,j))/deter

                matrix(i,j,8) = (point1_x(i,j)-point3_x(i,j))/deter

                matrix(i,j,9) = (point2_x(i,j)-point1_x(i,j))/deter

            ENDIF
          ENDDO
        ENDDO
      
      ENDSUBROUTINE


      SUBROUTINE getInterpolData(                                      &
                                 mesh,                                 &
                                 matrix,                               &
                                 matrixpress,                          &
                                 point1_x, point1_y,                   &
                                 point2_x, point2_y,                   &
                                 point3_x, point3_y,                   &
                                 n_x, n_y,                             &
                                 index1i, index1j,                     &
                                 index2i, index2j,                     &
                                 cflag,                                &
                                 theta,                                &
                                 c1, c2, thickness,                    &
                                 nx, ny, nz                            &
                                )

      IMPLICIT NONE

      real    :: mesh(nx,ny,nz,3)                                 
      real    :: matrix(nx,ny,9)                          
      real    :: matrixpress(nx,ny,9)                          
      real    :: point1_x(nx,ny), point1_y(nx,ny)                
      real    :: point2_x(nx,ny), point2_y(nx,ny)               
      real    :: point3_x(nx,ny), point3_y(nx,ny)                   
      real    :: n_x(nx,ny), n_y(nx,ny)                   
      real    :: theta, c1, c2
      integer :: cflag(nx,ny)                           
      integer :: i,j,k
      integer :: nx,ny,nz
      integer :: index1, index2
      integer :: index1i(nx,ny),   index1j(nx,ny)                   
      integer :: index2i(nx,ny),   index2j(nx,ny)                   

      real    :: A1, B1, D1
      real    :: A2, B2, D2
      real    :: A3, B3, D3
      real    :: A4, B4, D4

      real    :: A1p, B1p, D1p
      real    :: A2p, B2p, D2p
      real    :: A3p, B3p, D3p
      real    :: A4p, B4p, D4p

      real    :: Agen, Bgen, Dgen
      real    :: rango, deltax, deltay, deltaz

      real    :: angle, pi

      real plane1, plane2, plane3, plane4, t
      real, allocatable, dimension(:,:) :: mask2d
      real    :: thickness, magnormal

      allocate(mask2d(nx,ny))

      deltax = (mesh(nx,1,1,1) - mesh(1,1,1,1))  / float(nx-1)
      deltay = (mesh(1,ny,1,2) - mesh(1,1,1,2))  / float(ny-1)
      deltaz = (mesh(1,1,nz,3) - mesh(1,1,1,3))  / float(nz-1)

      rango = sqrt(deltax**2.d0 + deltay**2.d0 + deltaz**2.d0)


      A1 = 0.d0
      B1 = 1.d0
      
      A2 =  0.d0
      B2 = -1.d0

      A3 =  1.d0
      B3 =  0.d0

      A4 =  -1.d0
      B4 =   0.d0

      D1 = -thickness
      D2 = -thickness

      D3 = -c2
      D4 =  c1

      A1p = A1*cos(theta)  - B1*sin(theta)     
      B1p = A1*sin(theta)  + B1*cos(theta)     

      A2p = A2*cos(theta)  - B2*sin(theta)     
      B2p = A2*sin(theta)  + B2*cos(theta)     

      A3p = A3*cos(theta)  - B3*sin(theta)     
      B3p = A3*sin(theta)  + B3*cos(theta)     

      A4p = A4*cos(theta)  - B4*sin(theta)     
      B4p = A4*sin(theta)  + B4*cos(theta)     

      pi = 2.d0*acos(0.d0)

      DO j =1,ny
         DO i =1,nx
            k = 1

            mask2d(i,j) = 1.d0
            cflag(i,j)  = 1
            plane1 = A1p*mesh(i,j,k,1) + B1p*mesh(i,j,k,2) + D1
            plane2 = A2p*mesh(i,j,k,1) + B2p*mesh(i,j,k,2) + D2
            plane3 = A3p*mesh(i,j,k,1) + B3p*mesh(i,j,k,2) + D3
            plane4 = A4p*mesh(i,j,k,1) + B4p*mesh(i,j,k,2) + D4
            if((plane1.LE.0.d0).AND.(plane2.LE.0.d0).AND.(plane3.LE.0.d0)&
            .AND.(plane4.LE.0.d0))then 
                  mask2d(i,j) = 0.d0

            endif
         ENDDO
      ENDDO


      DO j =2,ny-1
         DO i =2,nx-1
              IF(mask2d(i,j).EQ.1.d0)then
                 IF((mask2d(i+1,j).EQ.0.0).OR.(mask2d(i-1,j).EQ.0.0).OR. &
                    (mask2d(i,j+1).EQ.0.0).OR.(mask2d(i,j-1).EQ.0.0))THEN
                   cflag(i,j) = -1  
                ENDIF
              ENDIF
            
         ENDDO
      ENDDO

      DO j =2,ny-1
         DO i =2,nx-1
            k = 1
            plane1 = A1p*mesh(i,j,k,1) + B1p*mesh(i,j,k,2) + D1
            plane2 = A2p*mesh(i,j,k,1) + B2p*mesh(i,j,k,2) + D2
            plane3 = A3p*mesh(i,j,k,1) + B3p*mesh(i,j,k,2) + D3
            plane4 = A4p*mesh(i,j,k,1) + B4p*mesh(i,j,k,2) + D4
            IF(cflag(i,j).EQ.-1)then
               IF((plane1.GE.0.d0).AND.(plane1.LE.rango))then
                 Agen = A1p 
                 Bgen = B1p 
                 Dgen = D1 

               ELSEIF((plane2.GE.0.d0).AND.(plane2.LE.rango))then
                 Agen = A2p 
                 Bgen = B2p 
                 Dgen = D2

               ELSEIF((plane3.GE.0.d0).AND.(plane3.LE.rango))then
                 Agen = A3p 
                 Bgen = B3p 
                 Dgen = D3
               ELSEIF((plane4.GE.0.d0).AND.(plane4.LE.rango))then
                 Agen = A4p 
                 Bgen = B4p 
                 Dgen = D4

               ELSE
                 write(6,*) 'problemas'

               ENDIF
            
               t = -(Dgen + mesh(i,j,k,1)*Agen + mesh(i,j,k,2)*Bgen) / &
                    (Agen*Agen + Bgen*Bgen) 

               angle = acos(Agen)/ sqrt(Agen**2.d0 + Bgen**2.d0)
 
               if(((angle.GE.pi/6.d0).AND.(angle.LE.pi/3.d0)).OR.&
                 ((angle.GE.2.d0*pi/3.d0).AND.(angle.LE.5.d0*pi/6.d0)))then

                  index1i(i,j) = i + int(sign(1.0, Agen))
                  index1j(i,j) = j

                  index2i(i,j) = i
                  index2j(i,j) = j + int(sign(1.0, Bgen))

               else if((angle.LE.pi/6.d0).OR.(angle.GE.5.d0*pi/6.d0))then

                  index1i(i,j) = i + int(sign(1.0, Agen)) 
                  index1j(i,j) = j

                  index2i(i,j) = i + int(sign(1.0, Agen))
                  index2j(i,j) = j + int(sign(1.0, Bgen))
               else

                  index1i(i,j) = i
                  index1j(i,j) = j + int(sign(1.0, Bgen))

                  index2i(i,j) = i + int(sign(1.0, Agen))
                  index2j(i,j) = j + int(sign(1.0, Bgen))
               endif

               point1_x(i,j) = mesh(i,j,k,1) + t*Agen
               point1_y(i,j) = mesh(i,j,k,2) + t*Bgen

               point2_x(i,j) = mesh(index1i(i,j),index1j(i,j),k,1) 
               point2_y(i,j) = mesh(index1i(i,j),index1j(i,j),k,2)

               point3_x(i,j) = mesh(index2i(i,j),index2j(i,j),k,1) 
               point3_y(i,j) = mesh(index2i(i,j),index2j(i,j),k,2)


               n_x(i,j) = mesh(i,j,k,1) - point1_x(i,j) 
               n_y(i,j) = mesh(i,j,k,2) - point1_y(i,j) 

               magnormal = sqrt(n_x(i,j)**2.d0 + n_y(i,j)**2.d0)

               n_x(i,j) = n_x(i,j) /magnormal 
               n_y(i,j) = n_y(i,j) /magnormal 

            ENDIF

         ENDDO
      ENDDO

      call getMatrix(matrix, point1_x, point1_y,point2_x, point2_y,    &
                     point3_x, point3_y,cflag, nx, ny)

      call getPressureMatrix(matrixpress, n_x, n_y,point2_x, point2_y, &
                     point3_x, point3_y,cflag, nx, ny)

      ENDSUBROUTINE

      SUBROUTINE makeAspa(                                             &
                          mesh,                                        &
                          mask,                                        &
                          theta,                                       &
                          c1, c2, z1, z2, thickness,                   &
                          nx, ny, nz                                   &
                         )

      IMPLICIT NONE
      real    :: mesh(nx,ny,nz,3)                                 
      real    :: mask(nx,ny,nz)                          
      real    :: theta
      integer :: nx,ny,nz, i,j,k

      real    :: A1, B1, D1
      real    :: A2, B2, D2
      real    :: A3, B3, D3
      real    :: A4, B4, D4

      real    :: A1p, B1p, D1p
      real    :: A2p, B2p, D2p
      real    :: A3p, B3p, D3p
      real    :: A4p, B4p, D4p
      real    :: plane1, plane2, plane3, plane4

      real    :: thickness, c1, c2, z1, z2

      A1 = 0.d0
      B1 = 1.d0
      
      A2 =  0.d0
      B2 = -1.d0

      A3 =  1.d0
      B3 =  0.d0

      A4 =  -1.d0
      B4 =   0.d0


      D1 = -thickness
      D2 = -thickness

      D3 = -c2
      D4 =  c1

      A1p = A1*cos(theta)  - B1*sin(theta)     
      B1p = A1*sin(theta)  + B1*cos(theta)     

      A2p = A2*cos(theta)  - B2*sin(theta)     
      B2p = A2*sin(theta)  + B2*cos(theta)     

      A3p = A3*cos(theta)  - B3*sin(theta)     
      B3p = A3*sin(theta)  + B3*cos(theta)     

      A4p = A4*cos(theta)  - B4*sin(theta)     
      B4p = A4*sin(theta)  + B4*cos(theta)     

      DO k =1,nz
         DO j =1,ny
            DO i =1,nx

              plane1 = A1p*mesh(i,j,k,1) + B1p*mesh(i,j,k,2) + D1
              plane2 = A2p*mesh(i,j,k,1) + B2p*mesh(i,j,k,2) + D2
              plane3 = A3p*mesh(i,j,k,1) + B3p*mesh(i,j,k,2) + D3
              plane4 = A4p*mesh(i,j,k,1) + B4p*mesh(i,j,k,2) + D4
              if((plane1.LE.0.d0).AND.(plane2.LE.0.d0).AND.(plane3.LE.0.d0)&
              .AND.(plane4.LE.0.d0).AND.(mesh(i,j,k,3).GE.z1).AND.&
               (mesh(i,j,k,3).LE.z2))mask(i,j,k) = 0.d0
 
            ENDDO
         ENDDO
      ENDDO

      ENDSUBROUTINE


      SUBROUTINE getPressureMatrix(matrix,                      &
                                   n_x, n_y,                    &
                                   point2_x, point2_y,          &
                                   point3_x, point3_y,          &
                                   cflag,                       &
                                   nx,ny)
        IMPLICIT NONE
        real    :: matrix(nx, ny, 9)
        real    :: n_x(nx,ny),  n_y(nx,ny)
        real    :: point2_x(nx,ny),  point2_y(nx,ny)
        real    :: point3_x(nx,ny),  point3_y(nx,ny)
        integer :: nx, ny
        real    :: deter
        integer :: i,j,k
        integer :: cflag(nx,ny)

        DO j =2,ny-1
           DO i =2,nx-1

           IF(cflag(i,j).EQ.-1)then
        
                deter = n_x(i,j)*point2_y(i,j) &
                      - n_x(i,j)*point3_y(i,j) &
                      + point3_x(i,j)*n_y(i,j) &
                      - point2_x(i,j)*n_y(i,j) 


                matrix(i,j,1) = (point2_x(i,j)*point3_y(i,j)           &
                              -  point3_x(i,j)*point2_y(i,j))/deter

                matrix(i,j,2) = (point3_x(i,j)*n_y(i,j)           &
                              -  n_x(i,j)*point3_y(i,j))/deter

                matrix(i,j,3) = (n_x(i,j)*point2_y(i,j)          &
                              -  point2_x(i,j)*n_y(i,j))/deter
         
                matrix(i,j,4) = (point2_y(i,j)-point3_y(i,j))/deter

                matrix(i,j,5) = (-n_y(i,j))/deter

                matrix(i,j,6) = (n_y(i,j))/deter
         
         
                matrix(i,j,7) = (point3_x(i,j)-point2_x(i,j))/deter

                matrix(i,j,8) = (n_x(i,j))/deter

                matrix(i,j,9) = (-n_x(i,j))/deter

            ENDIF
          ENDDO
        ENDDO


      
      ENDSUBROUTINE

