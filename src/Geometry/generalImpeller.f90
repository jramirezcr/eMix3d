      module stlTools
         real, allocatable, dimension(:)     ::  normals, vertex1
         real, allocatable, dimension(:)     ::  vertex2, vertex3
         integer, allocatable, dimension(:)  ::  neigh
         integer, allocatable, dimension(:)  ::  index1, index2, index3
         integer                             ::  numOfTriangles
         integer                             ::  numOfPoints
         real                                ::  Ro
         integer, allocatable, dimension(:,:,:)  ::  tubeMask
      end module stlTools

!_____________________________________________________________
       SUBROUTINE rushton()
!__________________________________________________________________

      use dimensiones
      use variables
      use consadim
      use velocidades
      use mallagrid
      use tiempo
      use iso_c_binding
      use impellermotion
      use maskvar
      use stlTools
      implicit none
      include 'Geometry/interfaceImpeller.h'
 
      integer neq, irk, i, j, itc, nmcit, k, no_imp, n_p
   
      call impellerRead(numOfTriangles)

      no_imp = numOfTriangles
      n_p = no_imp*3
   
      allocate(neigh(no_imp*3))
      allocate(normals(no_imp*3))
      allocate(vertex1(no_imp*3))
      allocate(vertex2(no_imp*3))
      allocate(vertex3(no_imp*3))
      allocate(index1(Nx*Ny*Nz))
      allocate(index2(Nx*Ny*Nz))
      allocate(index3(Nx*Ny*Nz))
      allocate(tubeMask(Nx,Ny,Nz))

      numOfPoints = 0

      DO k = 2, Nz-1
        DO j = 2, Ny-1
          DO i = 2, Nx-1
             numOfPoints = numOfPoints + 1
             index1(numOfPoints) = i
             index2(numOfPoints) = j
             index3(numOfPoints) = k
             
          ENDDO
        ENDDO
      ENDDO

      write(6,*) "vcvvccvvcvc", numOfPoints

      call impellerSource(neigh,                              &
                          normals, vertex1, vertex2, vertex3, n_p)
      do i = 1, no_imp
         write(6,*) neigh(i),neigh(i+1),neigh(i+2) 
      enddo


      call impellerFindMask()


      RETURN
      END SUBROUTINE rushton

!-----------------------------------------------------------------

       SUBROUTINE moveImpeller()
!__________________________________________________________________

      use dimensiones
      use variables
      use consadim
      use velocidades
      use mallagrid
      use tiempo
      use iso_c_binding
      use impellermotion
      use maskvar
      use stlTools
      implicit none
      include 'Geometry/interfaceImpeller.h'
 
      integer neq, irk, i, j, itc, nmcit, k, no_imp, n_p

      no_imp = numOfTriangles
      n_p = no_imp*3

      call impellerFindMask()

      CALL addAngularVel(um1, x, impellermask, -1.0*Ro, Nx, Ny, Nz)

      RETURN
      END SUBROUTINE moveImpeller

      SUBROUTINE rotateSTLonAxisZ(theta, normals,           &
                                  vertex1, vertex2, vertex3,&
                                  numOfTriangles)
      REAL    :: theta
      REAL    :: normals(3*numOfTriangles) 
      REAL    :: vertex1(3*numOfTriangles) 
      REAL    :: vertex2(3*numOfTriangles) 
      REAL    :: vertex3(3*numOfTriangles) 
      REAL    :: valXPrima
      REAL    :: valYPrima
      INTEGER :: numOfTriangles
      INTEGER :: triagleId

      DO triangleId = 0, numOfTriangles
         valXPrima= &
                     cos(theta)*normals(3*triangleID + 1)&
                    -sin(theta)*normals(3*triangleID + 2)

         valYPrima= &
                    sin(theta)*normals(3*triangleId + 1)&
                   +cos(theta)*normals(3*triangleId + 2)

         normals(3*triangleId + 1) = valXprima 
         normals(3*triangleId + 2) = valYprima 

         valXPrima= &
                     cos(theta)*vertex1(3*triangleID + 1)&
                    -sin(theta)*vertex1(3*triangleID + 2)

         valYPrima= &
                    sin(theta)*vertex1(3*triangleId + 1)&
                   +cos(theta)*vertex1(3*triangleId + 2)

         vertex1(3*triangleId + 1) = valXprima 
         vertex1(3*triangleId + 2) = valYprima 

         valXPrima= &
                     cos(theta)*vertex2(3*triangleID + 1)&
                    -sin(theta)*vertex2(3*triangleID + 2)

         valYPrima= &
                    sin(theta)*vertex2(3*triangleId + 1)&
                   +cos(theta)*vertex2(3*triangleId + 2)

         vertex2(3*triangleId + 1) = valXprima 
         vertex2(3*triangleId + 2) = valYprima 

         valXPrima= &
                     cos(theta)*vertex3(3*triangleID + 1)&
                    -sin(theta)*vertex3(3*triangleID + 2)

         valYPrima= &
                    sin(theta)*vertex3(3*triangleId + 1)&
                   +cos(theta)*vertex3(3*triangleId + 2)

         vertex3(3*triangleId + 1) = valXprima 
         vertex3(3*triangleId + 2) = valYprima 

      ENDDO

      END SUBROUTINE


