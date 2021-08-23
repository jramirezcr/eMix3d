
!____________________________________________________________________
      SUBROUTINE advectLevelSet(evolve)
!____________________________________________________________________

      use dimensiones
      use variables
      use tiempo
      use consadim
      use velocidades
      use mallagrid
      use jacobtools
      use source_calor
      use flow
      use multiphase
      use deltas
      IMPLICIT NONE
      integer evolve
      include 'LevelSet/interfaceAdvect.h'

      call advectLS_CUDA(pres, jbn, sobject, 1.d0 / deltax, &
                         1.d0/ deltay,  1.d0 / deltaz, &
                         Nx, Ny ,Nz, deltamin)

      RETURN
      END SUBROUTINE advectLevelSet


      SUBROUTINE heviside()
      
      use variables
      use dimensiones
      use multiphase
      IMPLICIT NONE
      integer :: i,j,k
      real    :: eps, PI
 
      PI = 2.0*acos(0.0)
      eps = 1.5*deltamax

      DO k=1,nz
         DO j=1,ny
            DO i=1,nx
               if(phi(i,j,k).LT.-eps) hevi(i,j,k) = 0.0
            
               if((phi(i,j,k).GE.-eps).AND.(phi(i,j,k).LE.eps))then
                  hevi(i,j,k) = 0.5 + phi(i,j,k) / ( 2.0* eps) &
                              + sin(PI*phi(i,j,k)/eps) / (2.0*PI)
               endif

               if(phi(i,j,k).GT.eps) hevi(i,j,k)  = 1.0

            ENDDO
          ENDDO
       ENDDO
     


      ENDSUBROUTINE
