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
      use multiphase
      use mallagrid
      use deltas
      IMPLICIT NONE
      integer :: i,j,k,l,m
      real :: p0,p1,p2,cdiv
      integer :: neq,mi,nit,ni

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
            e(i,j,k)=-um(i,j,k,2)*u(i,j,k,1)-pres(i,j,k)
            f(i,j,k)=-um(i,j,k,3)*u(i,j,k,1)
            g(i,j,k)=-um(i,j,k,4)*u(i,j,k,1)
          END DO
         END DO
        END DO

        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p0=jbnm(i,j,k,1)*ddvelp(i,j,k,1)
           p1=jbnm(i,j,k,5)*dcvel(i,j,k,5)
           p2=jbnm(i,j,k,9)*dcvel(i,j,k,9)
           evm(i,j,k)=vis(i,j,k)*cdiv*(2.*p0-p1-p2)
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbnm(i,j,k,1)*ddvelp(i,j,k,2)
           p2=jbnm(i,j,k,5)*dcvel(i,j,k,4)
           fvm(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbnm(i,j,k,1)*dcvel(i,j,k,3)
           p2=jbnm(i,j,k,9)*ddvelp(i,j,k,7)
           gvm(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO

        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p0=jbnp(i,j,k,1)*ddvelm(i,j,k,1)
           p1=jbnp(i,j,k,5)*dcvel(i,j,k,5)
           p2=jbnp(i,j,k,9)*dcvel(i,j,k,9)
           evp(i,j,k)=vis(i,j,k)*cdiv*(2.*p0-p1-p2)
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbnp(i,j,k,1)*dcvel(i,j,k,2)
           p2=jbnp(i,j,k,5)*ddvelm(i,j,k,4)
           fvp(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbnp(i,j,k,1)*dcvel(i,j,k,3)
           p2=jbnp(i,j,k,9)*ddvelm(i,j,k,7)
           gvp(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO



      ELSEIF (neq.eq.3) THEN
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
            e(i,j,k)=-um(i,j,k,2)*u(i,j,k,2)
            f(i,j,k)=-um(i,j,k,3)*u(i,j,k,2)-pres(i,j,k)
            g(i,j,k)=-um(i,j,k,4)*u(i,j,k,2)
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbnm(i,j,k,1)*ddvelp(i,j,k,2)
           p2=jbnm(i,j,k,5)*dcvel(i,j,k,4)
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
           p1=jbnm(i,j,k,5)*dcvel(i,j,k,6)
           p2=jbnm(i,j,k,9)*ddvelp(i,j,k,8)
           gvm(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO

        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbnp(i,j,k,1)*ddvelm(i,j,k,2)
           p2=jbnp(i,j,k,5)*dcvel(i,j,k,4)
           evp(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p0=jbnp(i,j,k,1)*dcvel(i,j,k,1)
           p1=jbnp(i,j,k,5)*ddvelm(i,j,k,5)
           p2=jbnp(i,j,k,9)*dcvel(i,j,k,9)
           fvp(i,j,k)=vis(i,j,k)*cdiv*(2.*p1-p0-p2)
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbnp(i,j,k,9)*ddvelm(i,j,k,8)
           p2=jbnp(i,j,k,5)*dcvel(i,j,k,6)
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
           g(i,j,k)=-um(i,j,k,4)*u(i,j,k,3) -pres(i,j,k) -  &
                     froude*rhotwophase(i,j,k)/deltaz
                                
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbnm(i,j,k,1)*ddvelp(i,j,k,3)
           p2=jbnm(i,j,k,9)*dcvel(i,j,k,7)
           evm(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbnm(i,j,k,9)*dcvel(i,j,k,8)
           p2=jbnm(i,j,k,5)*ddvelp(i,j,k,6)
           fvm(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p0=jbnm(i,j,k,1)*dcvel(i,j,k,1)
           p1=jbnm(i,j,k,5)*dcvel(i,j,k,5)
           p2=jbnm(i,j,k,9)*ddvelp(i,j,k,9)
           gvm(i,j,k)=vis(i,j,k)*cdiv*(2.*p2-p1-p0)
          END DO
         END DO
        END DO

        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbnp(i,j,k,1)*ddvelm(i,j,k,3)
           p2=jbnp(i,j,k,9)*dcvel(i,j,k,7)
           evp(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO

        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbnp(i,j,k,9)*dcvel(i,j,k,8)
           p2=jbnp(i,j,k,5)*ddvelm(i,j,k,6)
           fvp(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p0=jbnp(i,j,k,1)*dcvel(i,j,k,1)
           p1=jbnp(i,j,k,5)*dcvel(i,j,k,5)
           p2= jbnp(i,j,k,9)*ddvelm(i,j,k,9)
           gvp(i,j,k)=vis(i,j,k)*cdiv*(2.*p2-p1-p0)
          END DO
         END DO
        END DO


      ELSEIF (neq.eq.5) THEN
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
            e(i,j,k)=-u(i,j,k,1)*um(i,j,k,5)
            f(i,j,k)=-u(i,j,k,2)*um(i,j,k,5)
            g(i,j,k)=-u(i,j,k,3)*um(i,j,k,5)
          END DO
         END DO
        END DO

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            evm(i,j,k)=vis5(i,j,k)*(          &
                     jbnm(i,j,k,1)*dtempp(i,j,k,1))
           END DO
          END DO
         END DO
         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            fvm(i,j,k)=vis5(i,j,k)*(          &
                     jbnm(i,j,k,5)*dtempp(i,j,k,2))
           END DO
          END DO
         END DO
         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            gvm(i,j,k)=vis5(i,j,k)*(          &
                     jbnm(i,j,k,9)*dtempp(i,j,k,3))
           END DO
          END DO
         END DO
!
         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            evp(i,j,k)=vis5(i,j,k)*(          &
                     jbnp(i,j,k,1)*dtempm(i,j,k,1))
           END DO
          END DO
         END DO
         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            fvp(i,j,k)=vis5(i,j,k)*(          &
                     jbnp(i,j,k,5)*dtempm(i,j,k,2))
           END DO
          END DO
         END DO
         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            gvp(i,j,k)=vis5(i,j,k)*(          &
                     jbnp(i,j,k,9)*dtempm(i,j,k,3))
          END DO
          END DO
         END DO

      ENDIF
      return
      end subroutine flujos
