      SUBROUTINE DERIV_VEL()
      use dimensiones
      use velocidades
      use derivvel
      use deltas
!      use cons_mac
!      use combustion
      IMPLICIT NONE
      integer i,j,k,l,m,ni



!     VELOCIDADES
       DO k=1,nz
        DO j=1,ny
         DO i=2,nx-1
          dcvel(i,j,k,1)=deltax*(u(i+1,j,k,1)-u(i-1,j,k,1))*0.5
          dcvel(i,j,k,2)=deltax*(u(i+1,j,k,2)-u(i-1,j,k,2))*0.5
          dcvel(i,j,k,3)=deltax*(u(i+1,j,k,3)-u(i-1,j,k,3))*0.5
         END DO
         dcvel(nx,j,k,1)=deltax*(u(nx,j,k,1)-u(nx-1,j,k,1))
         dcvel(1,j,k,1)=deltax*(u(2,j,k,1)-u(1,j,k,1))
         dcvel(nx,j,k,2)=deltax*(u(nx,j,k,2)-u(nx-1,j,k,2))
         dcvel(1,j,k,2)=deltax*(u(2,j,k,2)-u(1,j,k,2))
         dcvel(nx,j,k,3)=deltax*(u(nx,j,k,3)-u(nx-1,j,k,3))
         dcvel(1,j,k,3)=deltax*(u(2,j,k,3)-u(1,j,k,3))
        END DO
       END DO

       DO k=1,nz
        DO j=2,ny-1
         DO i=1,nx
          dcvel(i,j,k,4)=deltay*(u(i,j+1,k,1)-u(i,j-1,k,1))*0.5
          dcvel(i,j,k,5)=deltay*(u(i,j+1,k,2)-u(i,j-1,k,2))*0.5
          dcvel(i,j,k,6)=deltay*(u(i,j+1,k,3)-u(i,j-1,k,3))*0.5
         END DO
        END DO
       END DO
       DO k=1,nz
        DO i=1,nx
          dcvel(i,1,k,4)=deltay*(u(i,2,k,1)-u(i,1,k,1))
          dcvel(i,1,k,5)=deltay*(u(i,2,k,2)-u(i,1,k,2))
          dcvel(i,1,k,6)=deltay*(u(i,2,k,3)-u(i,1,k,3))
          dcvel(i,ny,k,4)=deltay*(u(i,ny,k,1)-u(i,ny-1,k,1))
          dcvel(i,ny,k,5)=deltay*(u(i,ny,k,2)-u(i,ny-1,k,2))
          dcvel(i,ny,k,6)=deltay*(u(i,ny,k,3)-u(i,ny-1,k,3))
        END DO
       END DO

       DO k=2,nz-1
        DO j=1,ny
         DO i=1,nx
          dcvel(i,j,k,7)=deltaz*(u(i,j,k+1,1)-u(i,j,k-1,1))*0.5
          dcvel(i,j,k,8)=deltaz*(u(i,j,k+1,2)-u(i,j,k-1,2))*0.5
          dcvel(i,j,k,9)=deltaz*(u(i,j,k+1,3)-u(i,j,k-1,3))*0.5
         END DO
        END DO
       END DO
       DO j=1,ny
        DO i=1,nx
         dcvel(i,j,1,7)=deltaz*(u(i,j,2,1)-u(i,j,nz,1))*0.5
         dcvel(i,j,1,8)=deltaz*(u(i,j,2,2)-u(i,j,nz,2))*0.5
         dcvel(i,j,1,9)=deltaz*(u(i,j,2,3)-u(i,j,nz,3))*0.5
         dcvel(i,j,nz,7)=deltaz*(u(i,j,1,1)-u(i,j,nz-1,1))*0.5
         dcvel(i,j,nz,8)=deltaz*(u(i,j,1,2)-u(i,j,nz-1,2))*0.5
         dcvel(i,j,nz,9)=deltaz*(u(i,j,1,3)-u(i,j,nz-1,3))*0.5
        END DO
       END DO

       DO k=1,nz
        DO j=1,ny
         DO i=1,nx-1
          ddvelp(i,j,k,1)=deltax*(u(i+1,j,k,1)-u(i,j,k,1))
          ddvelp(i,j,k,2)=deltax*(u(i+1,j,k,2)-u(i,j,k,2))
          ddvelp(i,j,k,3)=deltax*(u(i+1,j,k,3)-u(i,j,k,3))
         END DO
         ddvelp(nx,j,k,1)=ddvelp(nx-1,j,k,1)
         ddvelp(nx,j,k,2)=ddvelp(nx-1,j,k,2)
         ddvelp(nx,j,k,3)=ddvelp(nx-1,j,k,3)
        END DO
       END DO

       DO k=1,nz
        DO j=1,ny-1
         DO i=1,nx
          ddvelp(i,j,k,4)=deltay*(u(i,j+1,k,1)-u(i,j,k,1))
          ddvelp(i,j,k,5)=deltay*(u(i,j+1,k,2)-u(i,j,k,2))
          ddvelp(i,j,k,6)=deltay*(u(i,j+1,k,3)-u(i,j,k,3))
         END DO
        END DO
       END DO
       DO k=1,nz
         DO i=1,nx
          ddvelp(i,ny,k,4)=ddvelp(i,ny-1,k,4)
          ddvelp(i,ny,k,5)=ddvelp(i,ny-1,k,5)
          ddvelp(i,ny,k,6)=ddvelp(i,ny-1,k,6)
        END DO
       END DO

       DO k=1,nz-1
        DO j=1,ny
         DO i=1,nx
          ddvelp(i,j,k,7)=deltaz*(u(i,j,k+1,1)-u(i,j,k,1))
          ddvelp(i,j,k,8)=deltaz*(u(i,j,k+1,2)-u(i,j,k,2))
          ddvelp(i,j,k,9)=deltaz*(u(i,j,k+1,3)-u(i,j,k,3))
         END DO
        END DO
       END DO
       DO j=1,ny
        DO i=1,nx
         ddvelp(i,j,nz,7)=deltaz*(u(i,j,1,1)-u(i,j,nz,1))
         ddvelp(i,j,nz,8)=deltaz*(u(i,j,1,2)-u(i,j,nz,2))
         ddvelp(i,j,nz,9)=deltaz*(u(i,j,1,3)-u(i,j,nz,3))
        END DO
       END DO


       DO k=1,nz
        DO j=1,ny
         DO i=2,nx
          ddvelm(i,j,k,1)=deltax*(u(i,j,k,1)-u(i-1,j,k,1))
          ddvelm(i,j,k,2)=deltax*(u(i,j,k,2)-u(i-1,j,k,2))
          ddvelm(i,j,k,3)=deltax*(u(i,j,k,3)-u(i-1,j,k,3))
         END DO
         ddvelm(1,j,k,1)=ddvelm(2,j,k,1)
         ddvelm(1,j,k,2)=ddvelm(2,j,k,2)
         ddvelm(1,j,k,3)=ddvelm(2,j,k,3)
        END DO
       END DO

       DO k=1,nz
        DO j=2,ny
         DO i=1,nx
          ddvelm(i,j,k,4)=deltay*(u(i,j,k,1)-u(i,j-1,k,1))
          ddvelm(i,j,k,5)=deltay*(u(i,j,k,2)-u(i,j-1,k,2))
          ddvelm(i,j,k,6)=deltay*(u(i,j,k,3)-u(i,j-1,k,3))
         END DO
        END DO
       END DO
       DO k=1,nz
         DO i=1,nx
         ddvelm(i,1,k,4)=ddvelm(i,2,k,4)
         ddvelm(i,1,k,5)=ddvelm(i,2,k,5)
         ddvelm(i,1,k,6)=ddvelm(i,2,k,6)
        END DO
       END DO

       DO k=2,nz
        DO j=1,ny
         DO i=1,nx
          ddvelm(i,j,k,7)=deltaz*(u(i,j,k,1)-u(i,j,k-1,1))
          ddvelm(i,j,k,8)=deltaz*(u(i,j,k,2)-u(i,j,k-1,2))
          ddvelm(i,j,k,9)=deltaz*(u(i,j,k,3)-u(i,j,k-1,3))
         END DO
        END DO
       END DO
       DO j=1,ny
        DO i=1,nx
         ddvelm(i,j,1,7)=deltaz*(u(i,j,1,1)-u(i,j,nz,1))
         ddvelm(i,j,1,8)=deltaz*(u(i,j,1,2)-u(i,j,nz,2))
         ddvelm(i,j,1,9)=deltaz*(u(i,j,1,3)-u(i,j,nz,3))
        END DO
       END DO


!    CONCENTRACION

       DO k=1,nz
        DO j=1,ny
         DO i=2,nx-1
          dconcm(i,j,k,1)=deltax*(conc(i+1,j,k)-conc(i-1,j,k))*0.5
          dconcp(i,j,k,1)=deltax*(conc(i+1,j,k)-conc(i-1,j,k))*0.5
         END DO
        END DO
       END DO
       DO k=1,nz
        DO j=1,ny
         dconcm(1,j,k,1)=deltax*(conc(2,j,k)-conc(1,j,k))
         dconcm(nx,j,k,1)=deltax*(conc(nx,j,k)-conc(nx-1,j,k))
         dconcp(1,j,k,1)=deltax*(conc(2,j,k)-conc(1,j,k))
         dconcp(nx,j,k,1)=deltax*(conc(nx,j,k)-conc(nx-1,j,k))
         END DO
        END DO

       DO k=1,nz
        DO j=2,ny-1
         DO i=1,nx
          dconcm(i,j,k,2)=deltay*(conc(i,j+1,k)-conc(i,j-1,k))*0.5
          dconcp(i,j,k,2)=deltay*(conc(i,j+1,k)-conc(i,j-1,k))*0.5
         END DO
        END DO
       END DO
       DO k=1,nz
         DO i=1,nx
         dconcm(i,1,k,2)=deltay*(conc(i,2,k)-conc(i,1,k))
         dconcm(i,ny,k,2)=deltay*(conc(i,ny,k)-conc(i,ny-1,k))
         dconcp(i,1,k,2)=deltay*(conc(i,2,k)-conc(i,1,k))
         dconcp(i,ny,k,2)=deltay*(conc(i,ny,k)-conc(i,ny-1,k))
        END DO
       END DO

       DO k=2,nz-1
        DO j=1,ny
         DO i=1,nx
          dconcm(i,j,k,3)=deltaz*(conc(i,j,k+1)-conc(i,j,k-1))*0.5
          dconcp(i,j,k,3)=deltaz*(conc(i,j,k+1)-conc(i,j,k-1))*0.5
         END DO
        END DO
       END DO
       DO j=1,ny
        DO i=1,nx
         dconcm(i,j,1,3)=deltaz*(conc(i,j,2)-conc(i,j,nz))*0.5
         dconcm(i,j,nz,3)=deltaz*(conc(i,j,1)-conc(i,j,nz-1))*0.5
         dconcp(i,j,1,3)=deltaz*(conc(i,j,2)-conc(i,j,nz))*0.5
         dconcp(i,j,nz,3)=deltaz*(conc(i,j,1)-conc(i,j,nz-1))*0.5
        END DO
       END DO






!    TEMPERATURA

       DO k=1,nz
        DO j=1,ny
         DO i=2,nx-1
          dtempm(i,j,k,1)=deltax*(temp(i,j,k)-temp(i-1,j,k))
          dtempp(i,j,k,1)=deltax*(temp(i+1,j,k)-temp(i,j,k))
          dtemp(i,j,k,1)=deltax*(temp(i+1,j,k)-temp(i-1,j,k))*0.5d0
         END DO
         dtempm(1,j,k,1)=deltax*(temp(2,j,k)-temp(1,j,k))
         dtempm(nx,j,k,1)=deltax*(temp(nx,j,k)-temp(nx-1,j,k))
         dtemp(1,j,k,1)=deltax*(temp(2,j,k)-temp(1,j,k))
         dtempp(1,j,k,1)=deltax*(temp(2,j,k)-temp(1,j,k))
         dtempp(nx,j,k,1)=deltax*(temp(nx,j,k)-temp(nx-1,j,k))
         dtemp(nx,j,k,1)=deltax*(temp(nx,j,k)-temp(nx-1,j,k))
        END DO
       END DO

       DO k=1,nz
        DO j=2,ny-1
         DO i=1,nx
          dtempm(i,j,k,2)=deltay*(temp(i,j,k)-temp(i,j-1,k))
          dtempp(i,j,k,2)=deltay*(temp(i,j+1,k)-temp(i,j,k))
          dtempp(i,j,k,2)=deltay*(temp(i,j+1,k)-temp(i,j-1,k))*0.5d0
         END DO
        END DO
       END DO

       DO k=1,nz
         DO i=1,nx
         dtempm(i,1,k,2)=deltay*(temp(i,2,k)-temp(i,1,k))
         dtempm(i,ny,k,2)=deltay*(temp(i,ny,k)-temp(i,ny-1,k))
         dtempp(i,1,k,2)=deltay*(temp(i,2,k)-temp(i,1,k))
         dtempp(i,ny,k,2)=deltay*(temp(i,ny,k)-temp(i,ny-1,k))
         dtemp(i,1,k,2)=deltay*(temp(i,2,k)-temp(i,1,k))
         dtemp(i,ny,k,2)=deltay*(temp(i,ny,k)-temp(i,ny-1,k))
        END DO
       END DO

       DO k=2,nz-1
        DO j=1,ny
         DO i=1,nx
          dtempm(i,j,k,3)=deltaz*(temp(i,j,k)-temp(i,j,k-1))
          dtempp(i,j,k,3)=deltaz*(temp(i,j,k+1)-temp(i,j,k))
          dtemp(i,j,k,3)=deltaz*(temp(i,j,k+1)-temp(i-1,j,k))*0.5d0
         END DO
        END DO
       END DO

       DO j=1,ny
        DO i=1,nx
         dtempm(i,j,1,3)=deltaz*(temp(i,j,1)-temp(i,j,nz))
         dtempm(i,j,nz,3)=deltaz*(temp(i,j,nz)-temp(i,j,nz-1))
         dtempp(i,j,1,3)=deltaz*(temp(i,j,2)-temp(i,j,1))
         dtempp(i,j,nz,3)=deltaz*(temp(i,j,1)-temp(i,j,nz))
         dtemp(i,j,1,3)=deltaz*(temp(i,j,2)-temp(i,j,1))
         dtemp(i,j,nz,3)=deltaz*(temp(i,j,1)-temp(i,j,nz))
        END DO
       END DO

      return
      end subroutine deriv_vel
