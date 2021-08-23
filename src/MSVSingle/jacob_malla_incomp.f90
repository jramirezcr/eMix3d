      SUBROUTINE JACOBEANO()
      use dimensiones
      use jacobtools
      use mallagrid
      use deltas
      use cons_mac
      IMPLICIT NONE
      integer :: i,j,k,l,m
      real :: a


      DO k=1,nz
       DO j=1,ny
        DO i=1,nx-2
         jbnp(i,j,k,1)=deltax*(cp1*x(i,j,k,1)+cp2*x(i+1,j,k,1)  &
                              +cp3*x(i+2,j,k,1))
         jbnp(i,j,k,2)=deltax*(cp1*x(i,j,k,2)+cp2*x(i+1,j,k,2)  &
                             +cp3*x(i+2,j,k,2))
         jbnp(i,j,k,3)=deltax*(cp1*x(i,j,k,3)+cp2*x(i+1,j,k,3)  &
                              +cp3*x(i+2,j,k,3))
        END DO
       END DO
      END DO
      DO k=1,nz
       DO j=1,ny
        jbnp(nx-1,j,k,1)=deltax*(2.*x(nx,j,k,1)+3.*x(nx-1,j,k,1)-6.*  &
                      x(nx-2,j,k,1)+x(nx-3,j,k,1))/6.
        jbnp(nx-1,j,k,2)=deltax*(2.*x(nx,j,k,2)+3.*x(nx-1,j,k,2)-6.*  &
                      x(nx-2,j,k,2)+x(nx-3,j,k,2))/6.
        jbnp(nx-1,j,k,3)=deltax*(2.*x(nx,j,k,3)+3.*x(nx-1,j,k,3)-6.*  &
                      x(nx-2,j,k,3)+x(nx-3,j,k,3))/6.
        jbnp(nx,j,k,1)=deltax*(11.*x(nx,j,k,1)-18.*x(nx-1,j,k,1)+9.*  &
                      x(nx-2,j,k,1)-2.*x(nx-3,j,k,1))/6.
        jbnp(nx,j,k,2)=deltax*(11.*x(nx,j,k,2)-18.*x(nx-1,j,k,2)+9.*  &
                      x(nx-2,j,k,2)-2.*x(nx-3,j,k,2))/6.
        jbnp(nx,j,k,3)=deltax*(11.*x(nx,j,k,3)-18.*x(nx-1,j,k,3)+9.*  &
                      x(nx-2,j,k,3)-2.*x(nx-3,j,k,3))/6.
       END DO
      END DO

      DO k=1,nz
       DO j=1,ny-2
        DO i=1,nx
         jbnp(i,j,k,4)=deltay*(cp1*x(i,j,k,1)+cp2*x(i,j+1,k,1)  &
                              +cp3*x(i,j+2,k,1))
         jbnp(i,j,k,5)=deltay*(cp1*x(i,j,k,2)+cp2*x(i,j+1,k,2)  &
                              +cp3*x(i,j+2,k,2))
         jbnp(i,j,k,6)=deltay*(cp1*x(i,j,k,3)+cp2*x(i,j+1,k,3)  &
                              +cp3*x(i,j+2,k,3))
        END DO
       END DO
      END DO

      DO k=1,nz
       DO i=1,nx
        jbnp(i,ny-1,k,4)=deltay*(2.*x(i,ny,k,1)+3.*x(i,ny-1,k,1)-6.* &
                      x(i,ny-2,k,1)+x(i,ny-3,k,1))/6.
        jbnp(i,ny-1,k,5)=deltay*(2.*x(i,ny,k,2)+3.*x(i,ny-1,k,2)-6.* &
                      x(i,ny-2,k,2)+x(i,ny-3,k,2))/6.
        jbnp(i,ny-1,k,6)=deltay*(2.*x(i,ny,k,3)+3.*x(i,ny-1,k,3)-6.* &
                      x(i,ny-2,k,3)+x(i,ny-3,k,3))/6.
        jbnp(i,ny,k,4)=deltay*(11.*x(i,ny,k,1)-18.*x(i,ny-1,k,1)+9.* &
                      x(i,ny-2,k,1)-2.*x(i,ny-3,k,1))/6.
        jbnp(i,ny,k,5)=deltay*(11.*x(i,ny,k,2)-18.*x(i,ny-1,k,2)+9.* &
                      x(i,ny-2,k,2)-2.*x(i,ny-3,k,2))/6.
        jbnp(i,ny,k,6)=deltay*(11.*x(i,ny,k,3)-18.*x(i,ny-1,k,3)+9.* &
                      x(i,ny-2,k,3)-2.*x(i,ny-3,k,3))/6.
       END DO
      END DO

      DO k=1,nz-2
       DO j=1,ny
        DO i=1,nx
         jbnp(i,j,k,7)=deltaz*(cp1*x(i,j,k,1)+cp2*x(i,j,k+1,1)    &
                              +cp3*x(i,j,k+2,1))
         jbnp(i,j,k,8)=deltaz*(cp1*x(i,j,k,2)+cp2*x(i,j,k+1,2)    &
                              +cp3*x(i,j,k+2,2))
         jbnp(i,j,k,9)=deltaz*(cp1*x(i,j,k,3)+cp2*x(i,j,k+1,3)    &
                              +cp3*x(i,j,k+2,3))
        END DO
       END DO
      END DO

      DO j=1,ny
       DO i=1,nx
        jbnp(i,j,nz-1,7)=deltaz*(2.*x(i,j,nz,1)+3.*x(i,j,nz-1,1)-6.* &
                      x(i,j,nz-2,1)+x(i,j,nz-3,1))/6.
        jbnp(i,j,nz-1,8)=deltaz*(2.*x(i,j,nz,2)+3.*x(i,j,nz-1,2)-6.* &
                      x(i,j,nz-2,2)+x(i,j,nz-3,2))/6.
        jbnp(i,j,nz-1,9)=deltaz*(2.*x(i,j,nz,3)+3.*x(i,j,nz-1,3)-6.* &
                      x(i,j,nz-2,3)+x(i,j,nz-3,3))/6.
        jbnp(i,j,nz,7)=deltaz*(11.*x(i,j,nz,1)-18.*x(i,j,nz-1,1)+9.* &
                      x(i,j,nz-2,1)-2*x(i,j,nz-3,1))/6.
        jbnp(i,j,nz,8)=deltaz*(11.*x(i,j,nz,2)-18.*x(i,j,nz-1,2)+9.* &
                      x(i,j,nz-2,2)-2*x(i,j,nz-3,2))/6.
        jbnp(i,j,nz,9)=deltaz*(11.*x(i,j,nz,3)-18.*x(i,j,nz-1,3)+9.* &
                      x(i,j,nz-2,3)-2*x(i,j,nz-3,3))/6.
       END DO
      END DO



      DO k=1,nz
       DO j=1,ny
        DO i=1,nx
         jbnp(i,j,k,10)=jbnp(i,j,k,1)*jbnp(i,j,k,5)*jbnp(i,j,k,9)   &
                       +jbnp(i,j,k,2)*jbnp(i,j,k,5)*jbnp(i,j,k,7)   &
                       +jbnp(i,j,k,2)*jbnp(i,j,k,4)*jbnp(i,j,k,8)   &
                       -jbnp(i,j,k,1)*jbnp(i,j,k,6)*jbnp(i,j,k,8)   &
                       -jbnp(i,j,k,2)*jbnp(i,j,k,4)*jbnp(i,j,k,9)   &
                       -jbnp(i,j,k,3)*jbnp(i,j,k,5)*jbnp(i,j,k,7) 

        END DO
       END DO
      END DO

      DO k=1,nz
       DO j=1,ny
        DO i=1,nx
         IF (abs(jbnp(i,j,k,10)).le.1.e-12) THEN
          print *,'ERROR EN EL JACOBIANO P  ',i,j,k
          print *,'  VALOR CERO DE TRAZA :',jbnp(i,j,k,10)
          print *,'    ',jbnp(i,j,k,1),jbnp(i,j,k,5),jbnp(i,j,k,9)
         ENDIF
        END DO
       END DO
      END DO
      DO k=1,nz
       DO j=1,ny
        DO i=1,nx
         jbnp(i,j,k,11)=1./jbnp(i,j,k,10)
         t1=jbnp(i,j,k,5)*jbnp(i,j,k,9)-jbnp(i,j,k,6)*jbnp(i,j,k,8)
         t2=jbnp(i,j,k,4)*jbnp(i,j,k,9)-jbnp(i,j,k,6)*jbnp(i,j,k,7)
         t3=jbnp(i,j,k,4)*jbnp(i,j,k,8)-jbnp(i,j,k,5)*jbnp(i,j,k,7)
         t4=jbnp(i,j,k,2)*jbnp(i,j,k,9)-jbnp(i,j,k,3)*jbnp(i,j,k,8)
         t5=jbnp(i,j,k,1)*jbnp(i,j,k,9)-jbnp(i,j,k,3)*jbnp(i,j,k,7)
         t6=jbnp(i,j,k,1)*jbnp(i,j,k,8)-jbnp(i,j,k,2)*jbnp(i,j,k,7)
         t7=jbnp(i,j,k,2)*jbnp(i,j,k,6)-jbnp(i,j,k,3)*jbnp(i,j,k,5)
         t8=jbnp(i,j,k,1)*jbnp(i,j,k,6)-jbnp(i,j,k,3)*jbnp(i,j,k,4)
         t9=jbnp(i,j,k,1)*jbnp(i,j,k,5)-jbnp(i,j,k,2)*jbnp(i,j,k,4)
         jbnp(i,j,k,1)= t1*jbnp(i,j,k,11)
         jbnp(i,j,k,2)=-t2*jbnp(i,j,k,11)
         jbnp(i,j,k,3)= t3*jbnp(i,j,k,11)
         jbnp(i,j,k,4)=-t4*jbnp(i,j,k,11)
         jbnp(i,j,k,5)= t5*jbnp(i,j,k,11)
         jbnp(i,j,k,6)=-t6*jbnp(i,j,k,11)
         jbnp(i,j,k,7)= t7*jbnp(i,j,k,11)
         jbnp(i,j,k,8)=-t8*jbnp(i,j,k,11)
         jbnp(i,j,k,9)= t9*jbnp(i,j,k,11)
        END DO
       END DO
      END DO



      DO k=1,nz
       DO j=1,ny
        DO i=3,nx
         jbnm(i,j,k,1)=deltax*(cm1*x(i,j,k,1)+cm2*x(i-1,j,k,1)  &
                           +cm3*x(i-2,j,k,1))
         jbnm(i,j,k,2)=deltax*(cm1*x(i,j,k,2)+cm2*x(i-1,j,k,2)  &
                           +cm3*x(i-2,j,k,2))
         jbnm(i,j,k,3)=deltax*(cm1*x(i,j,k,3)+cm2*x(i-1,j,k,3)  &
                           +cm3*x(i-2,j,k,3))
        END DO
       END DO
      END DO

      DO k=1,nz
       DO j=1,ny
        jbnm(2,j,k,1)=deltax*(-2.*x(1,j,k,1)-3.*x(2,j,k,1)+6.*  &
                      x(3,j,k,1)-x(4,j,k,1))/6.
        jbnm(2,j,k,2)=deltax*(-2.*x(1,j,k,2)-3.*x(2,j,k,2)+6.*  &
                      x(3,j,k,2)-x(4,j,k,2))/6.
        jbnm(2,j,k,3)=deltax*(-2.*x(1,j,k,3)-3.*x(2,j,k,3)+6.*  &
                      x(3,j,k,3)-x(4,j,k,3))/6.
        jbnm(1,j,k,1)=deltax*(-11.*x(1,j,k,1)+18.*x(2,j,k,1)-9.*&
                      x(3,j,k,1)+2.*x(4,j,k,1))/6.
        jbnm(1,j,k,2)=deltax*(-11.*x(1,j,k,2)+18.*x(2,j,k,2)-9.*&
                      x(3,j,k,2)+2.*x(4,j,k,2))/6.
        jbnm(1,j,k,3)=deltax*(-11.*x(1,j,k,3)+18.*x(2,j,k,3)-9.*&
                      x(3,j,k,3)+2.*x(4,j,k,3))/6.
       END DO
      END DO


      DO k=1,nz
       DO j=3,ny
        DO i=1,nx
         jbnm(i,j,k,4)=deltay*(cm1*x(i,j,k,1)+cm2*x(i,j-1,k,1)  &
                            +cm3*x(i,j-2,k,1))
         jbnm(i,j,k,5)=deltay*(cm1*x(i,j,k,2)+cm2*x(i,j-1,k,2)  &
                            +cm3*x(i,j-2,k,2))
         jbnm(i,j,k,6)=deltay*(cm1*x(i,j,k,3)+cm2*x(i,j-1,k,3)  &
                            +cm3*x(i,j-2,k,3))
        END DO
       END DO
      END DO

      DO k=1,nz
       DO i=1,nx
        jbnm(i,2,k,4)=deltay*(-2.*x(i,1,k,1)-3.*x(i,2,k,1)+6.*   &
                      x(i,3,k,1)-x(i,4,k,1))/6.
        jbnm(i,2,k,5)=deltay*(-2.*x(i,1,k,2)-3.*x(i,2,k,2)+6.*   &
                      x(i,3,k,2)-x(i,4,k,2))/6.
        jbnm(i,2,k,6)=deltay*(-2.*x(i,1,k,3)-3.*x(i,2,k,3)+6.*   &
                      x(i,3,k,3)-x(i,4,k,3))/6.
        jbnm(i,1,k,4)=deltay*(-11.*x(i,1,k,1)+18.*x(i,2,k,1)-9.* &
                      x(i,3,k,1)+2.*x(i,4,k,1))/6.
        jbnm(i,1,k,5)=deltay*(-11.*x(i,1,k,2)+18.*x(i,2,k,2)-9.* &
                      x(i,3,k,2)+2.*x(i,4,k,2))/6.
        jbnm(i,1,k,6)=deltay*(-11.*x(i,1,k,3)+18.*x(i,2,k,3)-9.* &
                      x(i,3,k,3)+2.*x(i,4,k,3))/6.
       END DO
      END DO


      DO k=3,nz
       DO j=1,ny
        DO i=1,nx
         jbnm(i,j,k,7)=deltaz*(cm1*x(i,j,k,1)+cm2*x(i,j,k-1,1) &
                            +cm3*x(i,j,k-2,1))
         jbnm(i,j,k,8)=deltaz*(cm1*x(i,j,k,2)+cm2*x(i,j,k-1,2)  &
                            +cm3*x(i,j,k-2,2))
         jbnm(i,j,k,9)=deltaz*(cm1*x(i,j,k,3)+cm2*x(i,j,k-1,3)  &
                            +cm3*x(i,j,k-2,3))
        END DO
       END DO
      END DO

      DO j=1,ny
       DO i=1,nx
        jbnm(i,j,2,7)=deltaz*(-2.*x(i,j,1,1)-3.*x(i,j,2,1)+6.*  &
                      x(i,j,3,1)-x(i,j,4,1))/6.
        jbnm(i,j,2,8)=deltaz*(-2.*x(i,j,1,2)-3.*x(i,j,2,2)+6.*  &
                      x(i,j,3,2)-x(i,j,4,2))/6.
        jbnm(i,j,2,9)=deltaz*(-2.*x(i,j,1,3)-3.*x(i,j,2,3)+6.*  &
                      x(i,j,3,3)-x(i,j,4,3))/6.
        jbnm(i,j,1,7)=deltaz*(-11.*x(i,j,1,1)+18.*x(i,j,2,1)-9.*&
                      x(i,j,3,1)+2*x(i,j,4,1))/6.
        jbnm(i,j,1,8)=deltaz*(-11.*x(i,j,1,2)+18.*x(i,j,2,2)-9.*&
                      x(i,j,3,2)+2*x(i,j,4,2))/6.
        jbnm(i,j,1,9)=deltaz*(-11.*x(i,j,1,3)+18.*x(i,j,2,3)-9.*&
                      x(i,j,3,3)+2*x(i,j,4,3))/6.
       END DO
      END DO


      DO k=1,nz
       DO j=1,ny
        DO i=1,nx
         jbnm(i,j,k,10)=jbnm(i,j,k,1)*jbnm(i,j,k,5)*jbnm(i,j,k,9)   &
                       +jbnm(i,j,k,2)*jbnm(i,j,k,5)*jbnm(i,j,k,7)   &
                       +jbnm(i,j,k,2)*jbnm(i,j,k,4)*jbnm(i,j,k,8)   &
                       -jbnm(i,j,k,1)*jbnm(i,j,k,6)*jbnm(i,j,k,8)   &
                       -jbnm(i,j,k,2)*jbnm(i,j,k,4)*jbnm(i,j,k,9)   &
                       -jbnm(i,j,k,3)*jbnm(i,j,k,5)*jbnm(i,j,k,7)

        END DO
       END DO
      END DO
      DO k=1,nz
       DO j=1,ny
        DO i=1,nx
         IF (abs(jbnm(i,j,k,10)).le.1.e-12) THEN
          print *,'ERROR EN EL JACOBIANO M  ',i,j,k
          print *,'  VALOR CERO DE TRAZA :',jbnm(i,j,k,10)
          print *,'    ',jbnm(i,j,k,1),jbnm(i,j,k,5),jbnm(i,j,k,9)
         ENDIF
        END DO
       END DO
      END DO
      DO k=1,nz
       DO j=1,ny
        DO i=1,nx
         jbnm(i,j,k,11)=1./jbnm(i,j,k,10)
         t1=jbnm(i,j,k,5)*jbnm(i,j,k,9)-jbnm(i,j,k,6)*jbnm(i,j,k,8)
         t2=jbnm(i,j,k,4)*jbnm(i,j,k,9)-jbnm(i,j,k,6)*jbnm(i,j,k,7)
         t3=jbnm(i,j,k,4)*jbnm(i,j,k,8)-jbnm(i,j,k,5)*jbnm(i,j,k,7)
         t4=jbnm(i,j,k,2)*jbnm(i,j,k,9)-jbnm(i,j,k,3)*jbnm(i,j,k,8)
         t5=jbnm(i,j,k,1)*jbnm(i,j,k,9)-jbnm(i,j,k,3)*jbnm(i,j,k,7)
         t6=jbnm(i,j,k,1)*jbnm(i,j,k,8)-jbnm(i,j,k,2)*jbnm(i,j,k,7)
         t7=jbnm(i,j,k,2)*jbnm(i,j,k,6)-jbnm(i,j,k,3)*jbnm(i,j,k,5)
         t8=jbnm(i,j,k,1)*jbnm(i,j,k,6)-jbnm(i,j,k,3)*jbnm(i,j,k,4)
         t9=jbnm(i,j,k,1)*jbnm(i,j,k,5)-jbnm(i,j,k,2)*jbnm(i,j,k,4)
         jbnm(i,j,k,1)= t1*jbnm(i,j,k,11)
         jbnm(i,j,k,2)=-t2*jbnm(i,j,k,11)
         jbnm(i,j,k,3)= t3*jbnm(i,j,k,11)
         jbnm(i,j,k,4)=-t4*jbnm(i,j,k,11)
         jbnm(i,j,k,5)= t5*jbnm(i,j,k,11)
         jbnm(i,j,k,6)=-t6*jbnm(i,j,k,11)
         jbnm(i,j,k,7)= t7*jbnm(i,j,k,11)
         jbnm(i,j,k,8)=-t8*jbnm(i,j,k,11)
         jbnm(i,j,k,9)= t9*jbnm(i,j,k,11)
        END DO
       END DO
      END DO



 118  FORMAT(f16.8,f16.8,f16.8,f16.8)

      return
      end subroutine jacobeano
