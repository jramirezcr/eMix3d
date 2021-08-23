
       SUBROUTINE INITSTEP()
       use tiempo
       use dimensiones
       use mallagrid
       use afluid
       IMPLICIT NONE
      integer :: i,j,k,l,m

       up=0.5
       down=1.2
       do k=1,nz
        do j=1,ny
         do i=2,nx-1
          dx(i,j,k)=MIN(x(i,j,k,1)-x(i-1,j,k,1),   &
                        x(i+1,j,k,1)-x(i,j,k,1))
         ENDDO
        ENDDO
       ENDDO
       do k=1,nz
        do j=1,ny
         dx(1,j,k)=(x(2,j,k,1)-x(1,j,k,1))
         dx(nx,j,k)=(x(nx,j,k,1)-x(nx-1,j,k,1))
        ENDDO
       ENDDO

       do k=1,nz
        do j=2,ny-1
         do i=1,nx
          dy(i,j,k)=MIN(x(i,j,k,2)-x(i,j-1,k,2),   &
                        x(i,j+1,k,2)-x(i,j,k,2))
         ENDDO
        ENDDO
       ENDDO
       do k=1,nz
        do i=1,nx
         dy(i,1,k)=(x(i,2,k,2)-x(i,1,k,2))
         dy(i,ny,k)=(x(i,ny,k,2)-x(i,ny-1,k,2))
        ENDDO
       ENDDO

       do k=2,nz-1
        do j=1,ny
         do i=1,nx
          dz(i,j,k)=MIN(x(i,j,k,3)-x(i,j,k-1,3),   &
                        x(i,j,k+1,3)-x(i,j,k,3))
         ENDDO
        ENDDO
       ENDDO
       do j=1,ny
        do i=1,nx
         dz(i,j,1)=(x(i,j,2,3)-x(i,j,1,3))
         dz(i,j,nz)=(x(i,j,nz,3)-x(i,j,nz-1,3))
        ENDDO
       ENDDO
       dxyz_min=1000000000.0

       do k=1,nz
        do j=1,ny
         do i=1,nx
          dmin=min(dy(i,j,k),dx(i,j,k),dz(i,j,k))
          if(dmin.lt.dxyz_min)dxyz_min=dmin
         ENDDO
        ENDDO
       ENDDO
       write(6,*)dmin

       do k=1,nz
        do j=1,ny
         do i=1,nx
          dy(i,j,k)=1./dy(i,j,k)
          dx(i,j,k)=1./dx(i,j,k)
          dz(i,j,k)=1./dz(i,j,k)
         ENDDO
        ENDDO
       ENDDO



       RETURN
       END SUBROUTINE INITSTEP
     

       SUBROUTINE TSTEP()
       use tiempo
       use dimensiones
       use velocidades
       use viscosidades
       use consadim
       use flow
       use afluid
       IMPLICIT NONE
      integer :: i,j,k,l,m
       real :: c,uv,vv,wv,umaxi,cx,cy,cz, rey

       real :: diftstep, dxdy, dxdz, dydz, dxd, dyd, dzd,cfldif


       diftstep = 100000000000.0
       olddt = dt
       cfla=0.
       rey=1./reynolds
      DO k=1,nz
       DO j=1,ny
        DO i=1,nx
!          uv=max(ABS(u(i,j,k,1)),dx(i,j,k)*vis(i,j,k))
!          vv=max(ABS(u(i,j,k,2)),dy(i,j,k)*vis(i,j,k))
!          wv=max(ABS(u(i,j,k,3)),dz(i,j,k)*vis(i,j,k))
!          cx=sqrt(cs*cs+uv*uv)
!          cy=sqrt(cs*cs+vv*vv)
!          cz=sqrt(cs*cs+wv*wv)
!          cfla = MAX(cfla,dx(i,j,k)*(uv+cx), &
!                          dy(i,j,k)*(vv+cy), &
!                          dz(i,j,k)*(wv+cz))


          uv=ABS(u(i,j,k,1))
          vv=ABS(u(i,j,k,2))
          wv=ABS(u(i,j,k,3))
          uv=(uv+sqrt(cs*cs+uv*uv))
          vv=(vv+sqrt(cs*cs+vv*vv))
          wv=(wv+sqrt(cs*cs+wv*wv))
          cfla = MAX(cfla,dx(i,j,k)*(uv), &
                          dy(i,j,k)*(vv), &
                          dz(i,j,k)*(wv))

          dxdy = 1.0 / (dx(i,j,k)*dy(i,j,k))**2.0
          dxdz = 1.0 / (dx(i,j,k)*dz(i,j,k))**2.0
          dydz = 1.0 / (dy(i,j,k)*dz(i,j,k))**2.0

          dxd = 1.0/dx(i,j,k)**2.0
          dyd = 1.0/dy(i,j,k)**2.0     
          dzd = 1.0/dz(i,j,k)**2.0     

          cfldif = 0.25d0*dxd*dyd*dzd*reynolds / (dxdy + dxdz + dydz)

          diftstep = min(diftstep, cfldif) 

        END DO
       END DO
      END DO

      write(6,*) 'El valor del dt difusivo es:', diftstep
      cflm = olddt*cfla
!
!       ----------------------------------------------
!       RNUEVO PASO DE TIEMPO
!       ----------------------------------------------
!
      dt = olddt/cflm*cfl
      diffdt = dt - olddt
      dt = olddt + 0.5*(up+down)*diffdt &
                + 0.5*(up-down)*ABS(diffdt)

      !dt = min(dt, diftstep)


       umax=0.0

       do k=1,nz
        do j=1,ny
         do i=1,nx
!          umaxi=max(abs(u(i,j,k,1)),abs(u(i,j,k,2)),abs(u(i,j,k,3)))
          umaxi=sqrt(u(i,j,k,1)*u(i,j,k,1)+u(i,j,k,2)*u(i,j,k,2)+u(i,j,k,3)*u(i,j,k,3))
          if(umaxi.gt.umax)umax=umaxi
         ENDDO
        ENDDO
       ENDDO
 

      nc=max(1,int((umax/cs)/mach))
      nc=1
       RETURN
       END SUBROUTINE TSTEP

