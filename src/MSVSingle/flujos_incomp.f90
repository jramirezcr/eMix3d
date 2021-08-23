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
            ev(i,j,k)=0.0
            fv(i,j,k)=0.0
            gv(i,j,k)=0.0
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
           p0=jbn(i,j,k,1)*dvel(i,j,k,1)
           p1=jbn(i,j,k,5)*dvel(i,j,k,5)
           p2=jbn(i,j,k,9)*dvel(i,j,k,9)
           ev(i,j,k)=vis(i,j,k)*cdiv*(2.*p0-p1-p2)
          END DO
         END DO
        END DO

        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbn(i,j,k,1)*dvel(i,j,k,2)
           p2=jbn(i,j,k,5)*dvel(i,j,k,4)
           fv(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO


        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbn(i,j,k,1)*dvel(i,j,k,3)
           p2=jbn(i,j,k,9)*dvel(i,j,k,7)
           gv(i,j,k)=vis(i,j,k)*(p1+p2)
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
           p1=jbn(i,j,k,1)*dvel(i,j,k,2)
           p2=jbn(i,j,k,5)*dvel(i,j,k,4)
           ev(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO

        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p0=jbn(i,j,k,1)*dvel(i,j,k,1)
           p1=jbn(i,j,k,5)*dvel(i,j,k,5)
           p2=jbn(i,j,k,9)*dvel(i,j,k,9)
           fv(i,j,k)=vis(i,j,k)*cdiv*(2.*p1-p0-p2)
          END DO
         END DO
        END DO

        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbn(i,j,k,9)*dvel(i,j,k,8)
           p2=jbn(i,j,k,5)*dvel(i,j,k,6)
           gv(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO




      ELSEIF (neq.eq.4) THEN
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           e(i,j,k)=-um(i,j,k,2)*u(i,j,k,3)
           f(i,j,k)=-um(i,j,k,3)*u(i,j,k,3)
           g(i,j,k)=-um(i,j,k,4)*u(i,j,k,3)-pres(i,j,k)
          END DO
         END DO
        END DO

        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbn(i,j,k,1)*dvel(i,j,k,3)
           p2=jbn(i,j,k,9)*dvel(i,j,k,7)
           ev(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO

        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p1=jbn(i,j,k,9)*dvel(i,j,k,8)
           p2=jbn(i,j,k,5)*dvel(i,j,k,6)
           fv(i,j,k)=vis(i,j,k)*(p1+p2)
          END DO
         END DO
        END DO

        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           p0=jbn(i,j,k,1)*dvel(i,j,k,1)
           p1=jbn(i,j,k,5)*dvel(i,j,k,5)
           p2=jbn(i,j,k,9)*dvel(i,j,k,9)
           gv(i,j,k)=vis(i,j,k)*cdiv*(2.*p2-p1-p0)
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

!
         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            ev(i,j,k)=vis5(i,j,k)*(          &
                     jbn(i,j,k,1)*dtemp(i,j,k,1))
           END DO
          END DO
         END DO
         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            fv(i,j,k)=vis5(i,j,k)*(          &
                     jbn(i,j,k,5)*dtemp(i,j,k,2))
           END DO
          END DO
         END DO

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            gv(i,j,k)=vis5(i,j,k)*(          &
                     jbn(i,j,k,9)*dtemp(i,j,k,3))
          END DO
          END DO
         END DO

      ENDIF
      return
      end subroutine flujos
