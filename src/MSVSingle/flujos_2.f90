      SUBROUTINE flujos(neq)


      use mflujos
      use viscosidades
      use consadim
      use dimensiones
      use variables
      use velocidades
      use jacobtools
      use derivvel
      use acoustic
      use flow
      use combustion
      IMPLICIT NONE
      integer :: i,j,k,l,m
      real :: p0,p1,p2,cdiv
      integer :: neq,mi,nit,ni
      real :: gen

      gen = 1.0/ 1700.0

      cdiv=2./3.
      IF (neq.eq.1) THEN
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
            e(i,j,k)=-c1*um(i,j,k,2)
            f(i,j,k)=-c1*um(i,j,k,3)
            g(i,j,k)=-c1*um(i,j,k,4)
         END DO
        END DO
       END DO

       DO k=1,nz
        DO j=1,ny
         DO i=1,nx
            evp(i,j,k)=0.0
            fvp(i,j,k)=0.0
            gvp(i,j,k)=0.0
         END DO
        END DO
       END DO

       DO k=1,nz
        DO j=1,ny
         DO i=1,nx
            evm(i,j,k)=0.0
            fvm(i,j,k)=0.0
            gvm(i,j,k)=0.0
         END DO
        END DO
       END DO

      ELSEIF (neq.eq.2) THEN
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
            e(i,j,k)=-um(i,j,k,2)*u(i,j,k,1) - pres(i,j,k)
            f(i,j,k)=-um(i,j,k,3)*u(i,j,k,1)
            g(i,j,k)=-um(i,j,k,4)*u(i,j,k,1)
          END DO
         END DO
        END DO

        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p0=jbnm(i,j,k,1)*ddvelp(i,j,k,1)+  &
              jbnm(i,j,k,4)*ddvelp(i,j,k,4)+  &
              jbnm(i,j,k,7)*ddvelp(i,j,k,7)
           p1=jbnm(i,j,k,2)*dcvel(i,j,k,2)+  &
              jbnm(i,j,k,5)*dcvel(i,j,k,5)+  &
              jbnm(i,j,k,8)*dcvel(i,j,k,8)
           p2=jbnm(i,j,k,3)*dcvel(i,j,k,3)+  &
              jbnm(i,j,k,6)*dcvel(i,j,k,6)+  &
              jbnm(i,j,k,9)*dcvel(i,j,k,9)
           evm(i,j,k)=vis(i,j,k)*cdiv*(2.*p0-p1-p2)
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbnm(i,j,k,1)*ddvelp(i,j,k,2)+  &
              jbnm(i,j,k,2)*ddvelp(i,j,k,5)+  &
              jbnm(i,j,k,3)*ddvelp(i,j,k,8)
           p2=jbnm(i,j,k,2)*dcvel(i,j,k,1)+  &
              jbnm(i,j,k,5)*dcvel(i,j,k,4)+  &
              jbnm(i,j,k,8)*dcvel(i,j,k,7)
           fvm(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbnm(i,j,k,1)*dcvel(i,j,k,3)+  &
              jbnm(i,j,k,2)*dcvel(i,j,k,6)+  &
              jbnm(i,j,k,3)*dcvel(i,j,k,9)
           p2=jbnm(i,j,k,3)*ddvelp(i,j,k,1)+  &
              jbnm(i,j,k,6)*ddvelp(i,j,k,4)+  &
              jbnm(i,j,k,9)*ddvelp(i,j,k,7)
           gvm(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO

        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p0=jbnp(i,j,k,1)*ddvelm(i,j,k,1)+  &
              jbnp(i,j,k,4)*ddvelm(i,j,k,4)+  &
              jbnp(i,j,k,7)*ddvelm(i,j,k,7)
           p1=jbnp(i,j,k,2)*dcvel(i,j,k,2)+  &
              jbnp(i,j,k,5)*dcvel(i,j,k,5)+  &
              jbnp(i,j,k,8)*dcvel(i,j,k,8)
           p2=jbnp(i,j,k,3)*dcvel(i,j,k,3)+  &
              jbnp(i,j,k,6)*dcvel(i,j,k,6)+  &
              jbnp(i,j,k,9)*dcvel(i,j,k,9)
           evp(i,j,k)=vis(i,j,k)*cdiv*(2.*p0-p1-p2)
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbnp(i,j,k,1)*dcvel(i,j,k,2)+  &
              jbnp(i,j,k,4)*dcvel(i,j,k,5)+  &
              jbnp(i,j,k,7)*dcvel(i,j,k,8)
           p2=jbnp(i,j,k,2)*ddvelm(i,j,k,1)+  &
              jbnp(i,j,k,5)*ddvelm(i,j,k,4)+  &
              jbnp(i,j,k,8)*ddvelm(i,j,k,7)
           fvp(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbnp(i,j,k,1)*dcvel(i,j,k,3)+  &
              jbnp(i,j,k,4)*dcvel(i,j,k,6)+  &
              jbnp(i,j,k,7)*dcvel(i,j,k,9)
           p2=jbnp(i,j,k,3)*ddvelm(i,j,k,1)+  &
              jbnp(i,j,k,6)*ddvelm(i,j,k,4)+  &
              jbnp(i,j,k,9)*ddvelm(i,j,k,7)
           gvp(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO



      ELSEIF (neq.eq.3) THEN
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
            e(i,j,k)=-um(i,j,k,2)*u(i,j,k,2)
            f(i,j,k)=-um(i,j,k,3)*u(i,j,k,2) - pres(i,j,k)
            g(i,j,k)=-um(i,j,k,4)*u(i,j,k,2)
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbnm(i,j,k,1)*ddvelp(i,j,k,2)+  &
              jbnm(i,j,k,4)*ddvelp(i,j,k,5)+  &
              jbnm(i,j,k,7)*ddvelp(i,j,k,8)
           p2=jbnm(i,j,k,2)*dcvel(i,j,k,1)+  &
              jbnm(i,j,k,5)*dcvel(i,j,k,4)+  &
              jbnm(i,j,k,8)*dcvel(i,j,k,7)
           evm(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p0=jbnm(i,j,k,1)*dcvel(i,j,k,1)+  &
              jbnm(i,j,k,4)*dcvel(i,j,k,4)+  &
              jbnm(i,j,k,7)*dcvel(i,j,k,7)
           p1=jbnm(i,j,k,2)*ddvelp(i,j,k,2)+  &
              jbnm(i,j,k,5)*ddvelp(i,j,k,5)+  &
              jbnm(i,j,k,8)*ddvelp(i,j,k,8)
           p2=jbnm(i,j,k,3)*dcvel(i,j,k,3)+  &
              jbnm(i,j,k,6)*dcvel(i,j,k,6)+  &
              jbnm(i,j,k,9)*dcvel(i,j,k,9)
           fvm(i,j,k)=vis(i,j,k)*cdiv*(2.*p1-p0-p2)
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbnm(i,j,k,2)*dcvel(i,j,k,3)+  &
              jbnm(i,j,k,5)*dcvel(i,j,k,6)+  &
              jbnm(i,j,k,8)*dcvel(i,j,k,9)
           p2=jbnm(i,j,k,3)*ddvelp(i,j,k,2)+  &
              jbnm(i,j,k,6)*ddvelp(i,j,k,5)+  &
              jbnm(i,j,k,9)*ddvelp(i,j,k,8)
           gvm(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO

        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbnp(i,j,k,1)*ddvelm(i,j,k,2)+  &
              jbnp(i,j,k,4)*ddvelm(i,j,k,5)+  &
              jbnp(i,j,k,7)*ddvelm(i,j,k,8)
           p2=jbnp(i,j,k,2)*dcvel(i,j,k,1)+  &
              jbnp(i,j,k,5)*dcvel(i,j,k,4)+  &
              jbnp(i,j,k,8)*dcvel(i,j,k,7)
           evp(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p0=jbnp(i,j,k,1)*dcvel(i,j,k,1)+  &
              jbnp(i,j,k,4)*dcvel(i,j,k,4)+  &
              jbnp(i,j,k,7)*dcvel(i,j,k,7)
           p1=jbnp(i,j,k,2)*ddvelm(i,j,k,2)+  &
              jbnp(i,j,k,5)*ddvelm(i,j,k,5)+  &
              jbnp(i,j,k,8)*ddvelm(i,j,k,8)
           p2=jbnp(i,j,k,3)*dcvel(i,j,k,3)+  &
              jbnp(i,j,k,6)*dcvel(i,j,k,6)+  &
              jbnp(i,j,k,9)*dcvel(i,j,k,9)
           fvp(i,j,k)=vis(i,j,k)*cdiv*(2.*p1-p0-p2)
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbnp(i,j,k,3)*ddvelm(i,j,k,2)+  &
              jbnp(i,j,k,6)*ddvelm(i,j,k,5)+  &
              jbnp(i,j,k,9)*ddvelm(i,j,k,8)
           p2=jbnp(i,j,k,2)*dcvel(i,j,k,3)+  &
              jbnp(i,j,k,5)*dcvel(i,j,k,6)+  &
              jbnp(i,j,k,8)*dcvel(i,j,k,9)
           gvp(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO




      ELSEIF (neq.eq.4) THEN
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           e(i,j,k)=-um(i,j,k,2)*u(i,j,k,3)
           f(i,j,k)=-um(i,j,k,3)*u(i,j,k,3)
           g(i,j,k)=-um(i,j,k,4)*u(i,j,k,3)  -pres(i,j,k)
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbnm(i,j,k,1)*ddvelp(i,j,k,3)+  &
              jbnm(i,j,k,4)*ddvelp(i,j,k,6)+  &
              jbnm(i,j,k,7)*ddvelp(i,j,k,9)
           p2=jbnm(i,j,k,3)*dcvel(i,j,k,1)+  &
              jbnm(i,j,k,6)*dcvel(i,j,k,4)+  &
              jbnm(i,j,k,9)*dcvel(i,j,k,7)
           evm(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbnm(i,j,k,3)*dcvel(i,j,k,2)+  &
              jbnm(i,j,k,6)*dcvel(i,j,k,5)+  &
              jbnm(i,j,k,9)*dcvel(i,j,k,8)
           p2=jbnm(i,j,k,2)*ddvelp(i,j,k,3)+  &
              jbnm(i,j,k,5)*ddvelp(i,j,k,6)+  &
              jbnm(i,j,k,8)*ddvelp(i,j,k,9)
           fvm(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p0=jbnm(i,j,k,1)*dcvel(i,j,k,1)+  &
              jbnm(i,j,k,4)*dcvel(i,j,k,4)+  &
              jbnm(i,j,k,7)*dcvel(i,j,k,7)
           p1=jbnm(i,j,k,2)*dcvel(i,j,k,2)+  &
              jbnm(i,j,k,5)*dcvel(i,j,k,5)+  &
              jbnm(i,j,k,8)*dcvel(i,j,k,8)
           p2=jbnm(i,j,k,3)*ddvelp(i,j,k,3)+  &
              jbnm(i,j,k,6)*ddvelp(i,j,k,6)+  &
              jbnm(i,j,k,9)*ddvelp(i,j,k,9)
           gvm(i,j,k)=vis(i,j,k)*cdiv*(2.*p2-p1-p0)
          END DO
         END DO
        END DO

        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbnp(i,j,k,1)*ddvelm(i,j,k,3)+  &
              jbnp(i,j,k,4)*ddvelm(i,j,k,6)+  &
              jbnp(i,j,k,7)*ddvelm(i,j,k,9)
           p2=jbnp(i,j,k,3)*dcvel(i,j,k,1)+  &
              jbnp(i,j,k,6)*dcvel(i,j,k,4)+  &
              jbnp(i,j,k,9)*dcvel(i,j,k,7)
           evp(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO

        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbnp(i,j,k,3)*dcvel(i,j,k,2)+  &
              jbnp(i,j,k,6)*dcvel(i,j,k,5)+  &
              jbnp(i,j,k,9)*dcvel(i,j,k,8)
           p2=jbnp(i,j,k,2)*ddvelm(i,j,k,3)+  &
              jbnp(i,j,k,5)*ddvelm(i,j,k,6)+  &
              jbnp(i,j,k,8)*ddvelm(i,j,k,9)
           fvp(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p0=jbnp(i,j,k,1)*dcvel(i,j,k,1)+  &
              jbnp(i,j,k,4)*dcvel(i,j,k,4)+  &
              jbnp(i,j,k,7)*dcvel(i,j,k,7)
           p1=jbnp(i,j,k,2)*dcvel(i,j,k,2)+  &
              jbnp(i,j,k,5)*dcvel(i,j,k,5)+  &
              jbnp(i,j,k,8)*dcvel(i,j,k,8)
           p2=jbnp(i,j,k,3)*ddvelm(i,j,k,3)+  &
              jbnp(i,j,k,6)*ddvelm(i,j,k,6)+  &
              jbnp(i,j,k,9)*ddvelm(i,j,k,9)
           gvp(i,j,k)=vis(i,j,k)*cdiv*(2.*p2-p1-p0)
          END DO
         END DO
        END DO


      ELSEIF (neq.eq.5) THEN
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
            e(i,j,k)=-u(i,j,k,1)*(um(i,j,k,5))! + Eckert*pres(i,j,k))
            f(i,j,k)=-u(i,j,k,2)*(um(i,j,k,5))! + Eckert*pres(i,j,k))
            g(i,j,k)=-u(i,j,k,3)*(um(i,j,k,5))! + Eckert*pres(i,j,k))
          END DO
         END DO
        END DO

       !difusive part 

       DO k=1,nz
         DO j=1,ny
          DO i=1,nx
              
             evm(i,j,k)=-vis5(i,j,k)*(jbnm(i,j,k,1)*dtemp(i,j,k,1))
             evp(i,j,k)=-vis5(i,j,k)*(jbnp(i,j,k,1)*dtemp(i,j,k,1))
          END DO
         END DO
        END DO
         
         
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
             fvm(i,j,k)=-vis5(i,j,k)*(jbnm(i,j,k,5)*dtemp(i,j,k,2))
             fvp(i,j,k)=-vis5(i,j,k)*(jbnp(i,j,k,5)*dtemp(i,j,k,2))
          END DO
         END DO
        END DO
         
       
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
             gvm(i,j,k)=-vis5(i,j,k)*(jbnm(i,j,k,9)*dtemp(i,j,k,3))
             gvp(i,j,k)=-vis5(i,j,k)*(jbnp(i,j,k,9)*dtemp(i,j,k,3))
          END DO
         END DO
        END DO

!     !viscous work part 
!      DO k=1,nz
!       DO j=1,ny
!        DO i=1,nx
!          p0=jbnm(i,j,k,1)*dcvel(i,j,k,1)
!          p1=jbnm(i,j,k,5)*dcvel(i,j,k,5)
!          p2=jbnm(i,j,k,9)*dcvel(i,j,k,9)
!          evm(i,j,k)=evm(i,j,k)+u(i,j,k,1)*vis2(i,j,k)*(2.*p0)
!
!          p0=jbnp(i,j,k,1)*dcvel(i,j,k,1)
!          p1=jbnp(i,j,k,5)*dcvel(i,j,k,5)
!          p2=jbnp(i,j,k,9)*dcvel(i,j,k,9)
!          evp(i,j,k)=evp(i,j,k)+u(i,j,k,1)*vis2(i,j,k)*(2.*p0)
!        END DO
!       END DO
!      END DO
!      
!      
!      DO k=1,nz
!       DO j=1,ny
!        DO i=1,nx
!          p1=jbnm(i,j,k,5)*dcvel(i,j,k,5)
!          p2=jbnm(i,j,k,1)*dcvel(i,j,k,1)
!          evm(i,j,k)=evm(i,j,k) + u(i,j,k,2)*vis2(i,j,k)*(p1+p2)
!
!          p1=jbnp(i,j,k,5)*dcvel(i,j,k,5)
!          p2=jbnp(i,j,k,1)*dcvel(i,j,k,1)
!          evp(i,j,k)=evp(i,j,k) + u(i,j,k,2)*vis2(i,j,k)*(p1+p2)
!        END DO
!       END DO
!      END DO
!      
!      
!      DO k=1,nz
!       DO j=1,ny
!        DO i=1,nx
!          p1=jbnm(i,j,k,9)*dcvel(i,j,k,9)
!          p2=jbnm(i,j,k,1)*dcvel(i,j,k,1)
!          evm(i,j,k)=evm(i,j,k)+u(i,j,k,3)*vis2(i,j,k)*(p1+p2)
!
!          p1=jbnp(i,j,k,9)*dcvel(i,j,k,9)
!          p2=jbnp(i,j,k,1)*dcvel(i,j,k,1)
!          evp(i,j,k)=evp(i,j,k)+u(i,j,k,3)*vis2(i,j,k)*(p1+p2)
!        END DO
!       END DO
!      END DO
!
!
!
!     DO k=1,nz
!       DO j=1,ny
!        DO i=1,nx
!          p1=jbnm(i,j,k,5)*dcvel(i,j,k,5)
!          p2=jbnm(i,j,k,1)*dcvel(i,j,k,1)
!          fvm(i,j,k)=fvm(i,j,k) + u(i,j,k,1)*vis2(i,j,k)*(p1+p2)
!
!          p1=jbnp(i,j,k,5)*dcvel(i,j,k,5)
!          p2=jbnp(i,j,k,1)*dcvel(i,j,k,1)
!          fvp(i,j,k)=fvp(i,j,k) + u(i,j,k,1)*vis2(i,j,k)*(p1+p2)
!        END DO
!       END DO
!      END DO
!      
!      
!      DO k=1,nz
!       DO j=1,ny
!        DO i=1,nx
!          p0=jbnm(i,j,k,1)*dcvel(i,j,k,1)
!          p1=jbnm(i,j,k,5)*dcvel(i,j,k,5)
!          p2=jbnm(i,j,k,9)*dcvel(i,j,k,9)
!          fvm(i,j,k)=fvm(i,j,k)+u(i,j,k,2)*vis2(i,j,k)*(2.e0*p1)
!
!          p0=jbnp(i,j,k,1)*dcvel(i,j,k,1)
!          p1=jbnp(i,j,k,5)*dcvel(i,j,k,5)
!          p2=jbnp(i,j,k,9)*dcvel(i,j,k,9)
!          fvp(i,j,k)=fvp(i,j,k)+u(i,j,k,2)*vis2(i,j,k)*(2.e0*p1)
!        END DO
!       END DO
!      END DO
!      
!
!      DO k=1,nz
!       DO j=1,ny
!        DO i=1,nx
!          p1=jbnm(i,j,k,5)*dcvel(i,j,k,5)
!          p2=jbnm(i,j,k,9)*dcvel(i,j,k,9)
!          fvm(i,j,k)=fvm(i,j,k) + u(i,j,k,3)*vis2(i,j,k)*(p1+p2)
!
!          p1=jbnp(i,j,k,5)*dcvel(i,j,k,5)
!          p2=jbnp(i,j,k,9)*dcvel(i,j,k,9)
!          fvp(i,j,k)=fvp(i,j,k) + u(i,j,k,3)*vis2(i,j,k)*(p1+p2)
!        END DO
!       END DO
!      END DO
!
!   
!      DO k=1,nz
!       DO j=1,ny
!        DO i=1,nx
!          p1=jbnm(i,j,k,9)*dcvel(i,j,k,9)
!          p2=jbnm(i,j,k,1)*dcvel(i,j,k,1)
!          gvm(i,j,k)=gvm(i,j,k) + u(i,j,k,1)*vis2(i,j,k)*(p1+p2)
!
!          p1=jbnp(i,j,k,9)*dcvel(i,j,k,9)
!          p2=jbnp(i,j,k,1)*dcvel(i,j,k,1)
!          gvp(i,j,k)=gvp(i,j,k) + u(i,j,k,1)*vis2(i,j,k)*(p1+p2)
!        END DO
!       END DO
!      END DO
!      
!      
!      DO k=1,nz
!       DO j=1,ny
!        DO i=1,nx
!          p1=jbnm(i,j,k,5)*dcvel(i,j,k,5)
!          p2=jbnm(i,j,k,9)*dcvel(i,j,k,9)
!          gvm(i,j,k)=gvm(i,j,k) + u(i,j,k,2)*vis2(i,j,k)*(p1+p2)
!
!          p1=jbnp(i,j,k,5)*dcvel(i,j,k,5)
!          p2=jbnp(i,j,k,9)*dcvel(i,j,k,9)
!          gvp(i,j,k)=gvp(i,j,k) + u(i,j,k,2)*vis2(i,j,k)*(p1+p2)
!        END DO
!       END DO
!      END DO
!      
!
!      DO k=1,nz
!       DO j=1,ny
!        DO i=1,nx
!          p0=jbnm(i,j,k,1)*dcvel(i,j,k,1)
!          p1=jbnm(i,j,k,5)*dcvel(i,j,k,5)
!          p2=jbnm(i,j,k,9)*dcvel(i,j,k,9)
!          gvm(i,j,k)=gvm(i,j,k) + u(i,j,k,3)*vis2(i,j,k)*(2.e0*p2)
!
!          p0=jbnp(i,j,k,1)*dcvel(i,j,k,1)
!          p1=jbnp(i,j,k,5)*dcvel(i,j,k,5)
!          p2=jbnp(i,j,k,9)*dcvel(i,j,k,9)
!          gvp(i,j,k)=gvp(i,j,k) + u(i,j,k,3)*vis2(i,j,k)*(2.e0*p2)
!        END DO
!       END DO
!      END DO
!       

      ELSEIF (neq.eq.6) THEN
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
            e(i,j,k)=-um(i,j,k,6)*u(i,j,k,1)
            f(i,j,k)=-um(i,j,k,6)*u(i,j,k,2)
            g(i,j,k)=-um(i,j,k,6)*u(i,j,k,3)
         END DO
        END DO
       END DO

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            evm(i,j,k)=vis5(i,j,k)*(          &
                     jbnm(i,j,k,1)*dconcp(i,j,k,1)+    &
                     jbnm(i,j,k,4)*dconcp(i,j,k,2)+    &
                     jbnm(i,j,k,7)*dconcp(i,j,k,3))
           END DO
          END DO
         END DO
         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            fvm(i,j,k)=vis5(i,j,k)*(          &
                     jbnm(i,j,k,2)*dconcp(i,j,k,1)+    &
                     jbnm(i,j,k,5)*dconcp(i,j,k,2)+    &
                     jbnm(i,j,k,8)*dconcp(i,j,k,3))
           END DO
          END DO
         END DO
         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            gvm(i,j,k)=vis5(i,j,k)*(          &
                     jbnm(i,j,k,3)*dconcp(i,j,k,1)+    &
                     jbnm(i,j,k,6)*dconcp(i,j,k,2)+    &
                     jbnm(i,j,k,9)*dconcp(i,j,k,3))
           END DO
          END DO
         END DO
!
         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            evp(i,j,k)=vis5(i,j,k)*(          &
                     jbnp(i,j,k,1)*dconcm(i,j,k,1)+    &
                     jbnp(i,j,k,4)*dconcm(i,j,k,2)+    &
                     jbnp(i,j,k,7)*dconcm(i,j,k,3))
           END DO
          END DO
         END DO
         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            fvp(i,j,k)=vis5(i,j,k)*(          &
                     jbnp(i,j,k,2)*dconcm(i,j,k,1)+    &
                     jbnp(i,j,k,5)*dconcm(i,j,k,2)+    &
                     jbnp(i,j,k,8)*dconcm(i,j,k,3))
           END DO
          END DO
         END DO
         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            gvp(i,j,k)=vis5(i,j,k)*(          &
                     jbnp(i,j,k,3)*dconcm(i,j,k,1)+    &
                     jbnp(i,j,k,6)*dconcm(i,j,k,2)+    &
                     jbnp(i,j,k,9)*dconcm(i,j,k,3))
          END DO
          END DO
         END DO


      ENDIF
      return
      end subroutine flujos
