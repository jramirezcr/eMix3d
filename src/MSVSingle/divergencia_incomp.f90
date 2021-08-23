      SUBROUTINE divergencia(neq,ncp,nmcit)
      use cons_mac
      use derivvel
      use jacobtools
      use mflujos
      use bob
      use bob_yz
      use variables
      use right
      use dmflujos
      use dimensiones
      use fronterapl
      use velocidades
      use tiempo
      use afluid
      use mallagrid
      use sgdmodel
      IMPLICIT NONE
      integer i,j,k
      integer neq,ncp,nmcit
      real sy,ue,masa1,masa2,masa3



       IF(ncp.eq.1) THEN
       IF((nmcit.eq.1).or.(nmcit.eq.3).or.(nmcit.eq.5).or.   &
           (nmcit.eq.7)) THEN


         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            ex(i,j,k)=jbnp(i,j,k,10)*(jbnp(i,j,k,1)*e(i,j,k))
           END DO
          END DO
         END DO

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            evx(i,j,k)=jbnp(i,j,k,10)*(jbnp(i,j,k,1)*evp(i,j,k))
           END DO
          END DO
         END DO

         CALL DIV_x(cp1,cp2,cp3,evx,ex,1)

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            dere(i,j,k)=jbnp(i,j,k,11)*dere(i,j,k)
            derev(i,j,k)=jbnp(i,j,k,11)*derev(i,j,k)
           END DO
          END DO
         END DO

        ELSE

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            ex(i,j,k)=jbnm(i,j,k,10)*(jbnm(i,j,k,1)*e(i,j,k))
           END DO
          END DO
         END DO

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            evx(i,j,k)=jbnm(i,j,k,10)*(jbnm(i,j,k,1)*evm(i,j,k))
           END DO
          END DO
         END DO


         CALL DIV_x(cm1,cm2,cm3,evx,ex,-1)

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            dere(i,j,k)=jbnm(i,j,k,11)*dere(i,j,k)
            derev(i,j,k)=jbnm(i,j,k,11)*derev(i,j,k)
           END DO
          END DO
         END DO

        ENDIF



 



        IF((nmcit.eq.1).or.(nmcit.eq.3).or.(nmcit.eq.5).or.   &
           (nmcit.eq.7)) THEN

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            fx(i,j,k)=jbnp(i,j,k,10)*(jbnp(i,j,k,5)*f(i,j,k))
           END DO
          END DO
         END DO

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            fvx(i,j,k)=jbnp(i,j,k,10)*(jbnp(i,j,k,5)*fvp(i,j,k))
           END DO
          END DO
         END DO

         CALL DIV_y(cp1,cp2,cp3,fvx,fx,1)

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            derf(i,j,k)=jbnp(i,j,k,11)*derf(i,j,k)
            derfv(i,j,k)=jbnp(i,j,k,11)*derfv(i,j,k)
           END DO
          END DO
         END DO
        
        ELSE

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            fx(i,j,k)=jbnm(i,j,k,10)*(jbnm(i,j,k,5)*f(i,j,k))
           END DO
          END DO
         END DO

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            fvx(i,j,k)=jbnm(i,j,k,10)*(jbnm(i,j,k,5)*fvm(i,j,k))
           END DO
          END DO
         END DO

         CALL DIV_y(cm1,cm2,cm3,fvx,fx,-1)

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            derf(i,j,k)=jbnm(i,j,k,11)*derf(i,j,k)
            derfv(i,j,k)=jbnm(i,j,k,11)*derfv(i,j,k)
           END DO
          END DO
         END DO

        ENDIF



        IF((nmcit.eq.1).or.(nmcit.eq.3).or.(nmcit.eq.5).or.   &
           (nmcit.eq.7)) THEN

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            gx(i,j,k)=jbnp(i,j,k,10)*(jbnp(i,j,k,9)*g(i,j,k))
           END DO
          END DO
         END DO

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            gvx(i,j,k)=jbnp(i,j,k,10)*(jbnp(i,j,k,9)*gvp(i,j,k))
           END DO
          END DO
         END DO
         CALL DIVnp_z(cp1,cp2,cp3,gvx,gx,1)

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            derg(i,j,k)=jbnp(i,j,k,11)*derg(i,j,k)
            dergv(i,j,k)=jbnp(i,j,k,11)*dergv(i,j,k)
           END DO
          END DO
         END DO

        ELSE

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            gx(i,j,k)=jbnm(i,j,k,10)*(jbnm(i,j,k,9)*g(i,j,k))

           END DO
          END DO
         END DO

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            gvx(i,j,k)=jbnm(i,j,k,10)*(jbnm(i,j,k,9)*gvm(i,j,k))
           END DO
          END DO
         END DO

         CALL DIVnp_z(cm1,cm2,cm3,gvx,gx,-1)

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            derg(i,j,k)=jbnm(i,j,k,11)*derg(i,j,k)
            dergv(i,j,k)=jbnm(i,j,k,11)*dergv(i,j,k)
           END DO
          END DO
         END DO

        ENDIF


       ELSE

     
        IF((nmcit.eq.1).or.(nmcit.eq.3).or.(nmcit.eq.5).or.   &
           (nmcit.eq.7)) THEN

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            ex(i,j,k)=jbnm(i,j,k,10)*(jbnm(i,j,k,1)*e(i,j,k))
           END DO
          END DO
         END DO

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            evx(i,j,k)=jbnm(i,j,k,10)*(jbnm(i,j,k,1)*evm(i,j,k))
           END DO
          END DO
         END DO


         CALL DIV_x(cm1,cm2,cm3,evx,ex,-1)
         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            dere(i,j,k)=jbnm(i,j,k,11)*dere(i,j,k)
            derev(i,j,k)=jbnm(i,j,k,11)*derev(i,j,k)
           END DO
          END DO
         END DO

        ELSE

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            ex(i,j,k)=jbnp(i,j,k,10)*(jbnp(i,j,k,1)*e(i,j,k))
           END DO
          END DO
         END DO

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
           evx(i,j,k)=jbnp(i,j,k,10)*(jbnp(i,j,k,1)*evp(i,j,k))
           END DO
          END DO
         END DO

         CALL DIV_x(cp1,cp2,cp3,evx,ex,1)

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            dere(i,j,k)=jbnp(i,j,k,11)*dere(i,j,k)
            derev(i,j,k)=jbnp(i,j,k,11)*derev(i,j,k)
           END DO
          END DO
         END DO

        ENDIF


        IF((nmcit.eq.1).or.(nmcit.eq.3).or.(nmcit.eq.5).or.   &
           (nmcit.eq.7)) THEN

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            fx(i,j,k)=jbnm(i,j,k,10)*(jbnm(i,j,k,5)*f(i,j,k))
           END DO
          END DO
         END DO

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            fvx(i,j,k)=jbnm(i,j,k,10)*(jbnm(i,j,k,5)*fvm(i,j,k))
           END DO
          END DO
         END DO

         CALL DIV_y(cm1,cm2,cm3,fvx,fx,-1)
         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            derf(i,j,k)=jbnm(i,j,k,11)*derf(i,j,k)
            derfv(i,j,k)=jbnm(i,j,k,11)*derfv(i,j,k)
           END DO
          END DO
         END DO

        ELSE

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            fx(i,j,k)=jbnp(i,j,k,10)*(jbnp(i,j,k,5)*f(i,j,k))
           END DO
          END DO
         END DO

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            fvx(i,j,k)=jbnp(i,j,k,10)*(jbnp(i,j,k,5)*fvp(i,j,k))
           END DO
          END DO
         END DO

         CALL DIV_y(cp1,cp2,cp3,fvx,fx,1)

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            derf(i,j,k)=jbnp(i,j,k,11)*derf(i,j,k)
            derfv(i,j,k)=jbnp(i,j,k,11)*derfv(i,j,k)
           END DO
          END DO
         END DO


        ENDIF


        IF((nmcit.eq.1).or.(nmcit.eq.3).or.(nmcit.eq.5).or.   &
           (nmcit.eq.7)) THEN

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            gx(i,j,k)=jbnm(i,j,k,10)*(jbnm(i,j,k,9)*g(i,j,k))
           END DO
          END DO
         END DO

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            gvx(i,j,k)=jbnm(i,j,k,10)*(jbnm(i,j,k,9)*gvm(i,j,k))
           END DO
          END DO
         END DO
         CALL DIVnp_z(cm1,cm2,cm3,gvx,gx,-1)
         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            derg(i,j,k)=jbnm(i,j,k,11)*derg(i,j,k)
            dergv(i,j,k)=jbnm(i,j,k,11)*dergv(i,j,k)
           END DO
          END DO
         END DO

        ELSE

        DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            gx(i,j,k)=jbnp(i,j,k,10)*(jbnp(i,j,k,9)*g(i,j,k))
           END DO
          END DO
         END DO

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            gvx(i,j,k)=jbnp(i,j,k,10)*(jbnp(i,j,k,9)*gvp(i,j,k))
           END DO
          END DO
         END DO
         CALL DIVnp_z(cp1,cp2,cp3,gvx,gx,1)

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
            derg(i,j,k)=jbnp(i,j,k,11)*derg(i,j,k)
            dergv(i,j,k)=jbnp(i,j,k,11)*dergv(i,j,k)
           END DO
          END DO
         END DO

        ENDIF

       ENDIF



      DO k=1,nz
       DO j=1,ny
        DO i=1,nx
         rs(i,j,k)=(dere(i,j,k)+derf(i,j,k)+derg(i,j,k))  
         rsv(i,j,k)=(derev(i,j,k)+derfv(i,j,k)+dergv(i,j,k))
        END DO
       END DO
      END DO

      return
      end SUBROUTINE divergencia



      SUBROUTINE div_x(c1,c2,c3,ev,e,ndir)
      use dimensiones
      use jacobtools
      use dmflujos
      use deltas
      use variables
      USE cons_mac

      IMPLICIT NONE
      INTEGER :: ndir
      integer i,j,k
      real :: c1,c2,c3,e1
      real,dimension (nx,ny,nz) :: ev,e




      IF(ndir.eq.-1) then
      DO k=1,nz
       DO j=1,ny
        DO i=3,nx
         dere(i,j,k)=deltax*(c1*e(i,j,k)+c2*e(i-1,j,k)  &
                            +c3*e(i-2,j,k))
        END DO
       END DO
      END DO
      DO k=1,nz
       DO j=1,ny
        e1=4.*e(1,j,k)-6.*e(2,j,k)+4.*e(3,j,k)-e(4,j,k)
        dere(2,j,k)=deltax*(c1*e(2,j,k)+c2*e(1,j,k)+c3*e1)
       END DO
      END DO

      DO k=1,nz
       DO j=1,ny
        DO i=3,nx
         derev(i,j,k)=deltax*(c1*ev(i,j,k)+c2*ev(i-1,j,k)  &
                             +c3*ev(i-2,j,k))
        END DO
       END DO
      END DO
      DO k=1,nz
       DO j=1,ny
        e1=4.*ev(1,j,k)-6.*ev(2,j,k)+4.*ev(3,j,k)-ev(4,j,k)
        derev(2,j,k)=deltax*(c1*ev(2,j,k)+c2*ev(1,j,k)+c3*e1)
       END DO
      END DO

      ELSE

      DO k=1,nz
       DO j=1,ny
        DO i=1,nx-2
         dere(i,j,k)=deltax*(c1*e(i,j,k)+c2*e(i+1,j,k)   &
                           +c3*e(i+2,j,k))
        END DO
       END DO
      END DO
      DO k=1,nz
       DO j=1,ny
        e1=4.*e(nx,j,k)-6.*e(nx-1,j,k)+4.*e(nx-2,j,k)-e(nx-3,j,k)
        dere(nx-1,j,k)=deltax*(c1*e(nx-1,j,k)+c2*e(nx,j,k)+c3*e1)
       END DO
       END DO
      DO k=1,nz
       DO j=1,ny
        DO i=1,nx-2
         derev(i,j,k)=deltax*(c1*ev(i,j,k)+c2*ev(i+1,j,k)   &
                             +c3*ev(i+2,j,k))
        END DO
       END DO
      END DO
       DO k=1,nz
        DO j=1,ny
         e1=4.*ev(nx,j,k)-6.*ev(nx-1,j,k)+4.*ev(nx-2,j,k)-ev(nx-3,j,k)
         derev(nx-1,j,k)=deltax*(c1*ev(nx-1,j,k)+c2*ev(nx,j,k)+c3*e1)
        END DO
       END DO

       ENDIF

      DO k=1,nz
       DO j=1,ny
        dere(nx,j,k)=0.0
        dere(1,j,k)=0.0
        derev(nx,j,k)=0.0
        derev(1,j,k)=0.0
       ENDDO
      ENDDO
 
      return

      end subroutine div_x



      SUBROUTINE div_y(c1,c2,c3,fv,f,ndir)
      use dimensiones
      use jacobtools
      use dmflujos
      use deltas
      use variables

      IMPLICIT NONE
      INTEGER :: ndir
      integer i,j,k
      real :: c1,c2,c3,f1
      real,dimension (nx,ny,nz) :: fv,f



      IF(ndir.eq.-1) then

      DO k=1,nz
       DO j=3,ny
        DO i=1,nx
        derf(i,j,k)=deltay*(c1*f(i,j,k)+c2*f(i,j-1,k)   &
                           +c3*f(i,j-2,k))
        END DO
       END DO
      END DO
      DO k=1,nz
       DO i=1,nx
        f1=4.*f(i,1,k)-6.*f(i,2,k)+4.*f(i,3,k)-f(i,4,k)
        derf(i,2,k)=deltay*(c1*f(i,2,k)+c2*f(i,1,k)+c3*f1)
       END DO
      END DO

      DO k=1,nz
       DO j=3,ny
        DO i=1,nx
        derfv(i,j,k)=deltay*(c1*fv(i,j,k)+c2*fv(i,j-1,k)   &
                            +c3*fv(i,j-2,k))
        END DO
       END DO
      END DO
      DO k=1,nz
       DO i=1,nx
        f1=4.*fv(i,1,k)-6.*fv(i,2,k)+4.*fv(i,3,k)-fv(i,4,k)
        derfv(i,2,k)=deltay*(c1*fv(i,2,k)+c2*fv(i,1,k)+c3*f1)
       END DO
      END DO

      ELSE

       DO k=1,nz
        DO j=1,ny-2
         DO i=1,nx
          derf(i,j,k)=deltay*(c1*f(i,j,k)+c2*f(i,j+1,k)  &
                          +c3*f(i,j+2,k))
         END DO
        END DO
       END DO
       DO k=1,nz
        DO i=1,nx
         f1=4.*f(i,ny,k)-6.*f(i,ny-1,k)+4.*f(i,ny-2,k)-f(i,ny-3,k)
         derf(i,ny-1,k)=deltay*(c1*f(i,ny-1,k)+c2*f(i,ny,k)+c3*f1)
        END DO
       END DO

       DO k=1,nz
        DO j=1,ny-2
         DO i=1,nx
          derfv(i,j,k)=deltay*(c1*fv(i,j,k)+c2*fv(i,j+1,k)  &
                          +c3*fv(i,j+2,k))
         END DO
        END DO
       END DO
       DO k=1,nz
        DO i=1,nx
         f1=4.*fv(i,ny,k)-6.*fv(i,ny-1,k)+4.*fv(i,ny-2,k)-fv(i,ny-3,k)
         derfv(i,ny-1,k)=deltay*(c1*fv(i,ny-1,k)+c2*fv(i,ny,k)+c3*f1)
        END DO
       END DO

      ENDIF

      DO k=1,nz
       DO i=1,nx
        derf(i,ny,k)=0.0
        derf(i,1,k)=0.0
        derfv(i,ny,k)=0.0
        derfv(i,1,k)=0.0
       ENDDO
      ENDDO

      return
      END SUBROUTINE div_y


      SUBROUTINE div_z(c1,c2,c3,gv,g,ndir)
      use dimensiones
      use jacobtools
      use dmflujos
      use deltas
      use variables

      IMPLICIT NONE
      INTEGER :: ndir
      integer i,j,k
      real :: c1,c2,c3,g1
      real,dimension (nx,ny,nz) :: gv,g


      IF(ndir.eq.-1) then

      DO k=3,nz
       DO j=1,ny
        DO i=1,nx
         derg(i,j,k)=deltaz*(c1*g(i,j,k)+c2*g(i,j,k-1)   &
                            +c3*g(i,j,k-2))
        END DO
       END DO
      END DO
      DO j=1,ny
       DO i=1,nx
        derg(i,j,2)=deltaz*(c1*g(i,j,2)+c2*g(i,j,1)+c3*g(i,j,nz))
        derg(i,j,1)=deltaz*(c1*g(i,j,1)+c2*g(i,j,nz)+c3*g(i,j,nz-1))
       END DO
      END DO

      DO k=3,nz
       DO j=1,ny
        DO i=1,nx
         dergv(i,j,k)=deltaz*(c1*gv(i,j,k)+c2*gv(i,j,k-1)  &
                            +c3*gv(i,j,k-2))
        END DO
       END DO
      END DO
      DO j=1,ny
       DO i=1,nx
        dergv(i,j,2)=deltaz*(c1*gv(i,j,2)+c2*gv(i,j,1)+c3*gv(i,j,nz))
        dergv(i,j,1)=deltaz*(c1*gv(i,j,1)+c2*gv(i,j,nz)+c3*gv(i,j,nz-1))
       END DO
      END DO

      ELSE

      DO k=1,nz-2
       DO j=1,ny
        DO i=1,nx
         derg(i,j,k)=deltaz*(c1*g(i,j,k)+c2*g(i,j,k+1)   &
                         +c3*g(i,j,k+2))
        END DO
       END DO
      END DO
      DO j=1,ny
       DO i=1,nx
        derg(i,j,nz-1)=deltaz*(c1*g(i,j,nz-1)+c2*g(i,j,nz)+c3*g(i,j,1))
        derg(i,j,nz)=deltaz*(c1*g(i,j,nz)+c2*g(i,j,1)+c3*g(i,j,2))
       END DO
      END DO

      DO k=1,nz-2
       DO j=1,ny
        DO i=1,nx
         dergv(i,j,k)=deltaz*(c1*gv(i,j,k)+c2*gv(i,j,k+1)   &
                         +c3*gv(i,j,k+2))
        END DO
       END DO
      END DO
      DO j=1,ny
       DO i=1,nx
        dergv(i,j,nz-1)=deltaz*(c1*gv(i,j,nz-1)+c2*gv(i,j,nz)+c3*gv(i,j,1))
        dergv(i,j,nz)=deltaz*(c1*gv(i,j,nz)+c2*gv(i,j,1)+c3*gv(i,j,2))
       END DO
      END DO

      ENDIF

  
      RETURN
      END SUBROUTINE div_z

      SUBROUTINE divnp_z(c1,c2,c3,gv,g,ndir)
      use dimensiones
      use jacobtools
      use dmflujos
      use deltas
      use variables

      IMPLICIT NONE
      INTEGER :: ndir
      integer i,j,k
      real :: c1,c2,c3,g1
      real,dimension (nx,ny,nz) :: gv,g


      IF(ndir.eq.-1) then

      DO k=3,nz
       DO j=1,ny
        DO i=1,nx
         derg(i,j,k)=deltaz*(c1*g(i,j,k)+c2*g(i,j,k-1)   &
                            +c3*g(i,j,k-2))
        END DO
       END DO
      END DO
      DO j=1,ny
       DO i=1,nx
        g1=4.*g(i,j,1)-6.*g(i,j,2)+4.*g(i,j,3)-g(i,j,4)
        derg(i,j,2)=deltaz*(c1*g(i,j,2)+c2*g(i,j,1)+c3*g1)
       END DO
      END DO

      DO k=3,nz
       DO j=1,ny
        DO i=1,nx
         dergv(i,j,k)=deltaz*(c1*gv(i,j,k)+c2*gv(i,j,k-1)  &
                            +c3*gv(i,j,k-2))
        END DO
       END DO
      END DO
      DO j=1,ny
       DO i=1,nx
        g1=4.*gv(i,j,1)-6.*gv(i,j,2)+4.*gv(i,j,3)-gv(i,j,4)
        dergv(i,j,2)=deltaz*(c1*gv(i,j,2)+c2*gv(i,j,1)+c3*g1)
       END DO
      END DO

      ELSE

      DO k=1,nz-2
       DO j=1,ny
        DO i=1,nx
         derg(i,j,k)=deltaz*(c1*g(i,j,k)+c2*g(i,j,k+1)   &
                         +c3*g(i,j,k+2))
        END DO
       END DO
      END DO
      DO j=1,ny
       DO i=1,nx
        g1=4.*g(i,j,nz)-6.*g(i,j,nz-1)+4.*g(i,j,nz-2)-g(i,j,nz-3)
        derg(i,j,nz-1)=deltaz*(c1*g(i,j,nz-1)+c2*g(i,j,nz)+c3*g1)
       END DO
      END DO

      DO k=1,nz-2
       DO j=1,ny
        DO i=1,nx
         dergv(i,j,k)=deltaz*(c1*gv(i,j,k)+c2*gv(i,j,k+1)   &
                         +c3*gv(i,j,k+2))
        END DO
       END DO
      END DO
      DO j=1,ny
       DO i=1,nx
        g1=4.*gv(i,j,nz)-6.*gv(i,j,nz-1)+4.*gv(i,j,nz-2)-gv(i,j,nz-3)
        dergv(i,j,nz-1)=deltaz*(c1*gv(i,j,nz-1)+c2*gv(i,j,nz)+c3*g1)
       END DO
      END DO

      ENDIF

      DO j=1,ny
       DO i=1,nx
        derg(i,j,nz)=0.0
        derg(i,j,1)=0.0
        dergv(i,j,nz)=0.0
        dergv(i,j,1)=0.0
       END DO
      END DO
  
      RETURN
      END SUBROUTINE divnp_z

