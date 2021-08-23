      SUBROUTINE viscosidad()
      use dimensiones
      use velocidades
      use viscosidades
      use sgdmodel
      use flow
      use acoustic
      use consadim
      IMPLICIT NONE
      integer :: i,j,k,l,m

      DO k=1,nz
       DO j=1,ny
        DO i=1,nx
          vis(i,j,k)=1.0
          vis5(i,j,k)=amut(i,j,k)/0.7d0
          vis(i,j,k)=c4*vis(i,j,k)+amut(i,j,k)
        END DO
       END DO
      END DO

      RETURN
      END SUBROUTINE VISCOSIDAD
