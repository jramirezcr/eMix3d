
      Subroutine TVD()
       use variables
       use dimensiones
       use velocidades
       use tiempo
       use tvdvar
       use flow

       IMPLICIT NONE
       integer i,j,k
       real c,cr,cons

        do k=1,nz
        do j=1,ny
        do i=2,nx-1
        do m=1,nd
         dmas(i,j,k,m)=um1(i+1,j,k,m)-um1(i,j,k,m)
         dmen(i,j,k,m)=um1(i,j,k,m)-um1(i-1,j,k,m)
        enddo
        enddo
        enddo
        enddo
        do k=1,nz
        do j=1,ny
        do m=1,nd
         dmas(1,j,k,m)=um1(2,j,k,m)-um1(1,j,k,m)
         dmas(nx,j,k,m)=dmas(nx-1,j,k,m)
         dmen(1,j,k,m)=dmen(2,j,k,m)
         dmen(nx,j,k,m)=um1(nx,j,k,m)-um1(nx-1,j,k,m)
        enddo
        enddo
        enddo
         
        do k=1,nz
        do j=1,ny
        do i=1,nx
         rmas(i,j,k)=(dmen(i,j,k,1)*dmas(i,j,k,1)+dmen(i,j,k,2)*dmas(i,j,k,2)+dmen(i,j,k,3)*dmas(i,j,k,3))/ &
              (dmas(i,j,k,1)*dmas(i,j,k,1)+dmas(i,j,k,2)*dmas(i,j,k,2)+dmas(i,j,k,3)*dmas(i,j,k,3))
         rmen(i,j,k)=(dmen(i,j,k,1)*dmas(i,j,k,1)+dmen(i,j,k,2)*dmas(i,j,k,2)+dmen(i,j,k,3)*dmas(i,j,k,3))/ &
              (dmen(i,j,k,1)*dmen(i,j,k,1)+dmen(i,j,k,2)*dmen(i,j,k,2)+dmen(i,j,k,3)*dmen(i,j,k,3))
        enddo
        enddo
        do k=1,nz
        do j=1,ny
        do i=1,nx
          c=SQRT(grav*h(i,j,k))
          cr = dx(i,j,k)*(ABS(u(i,j,k,1))+c)*dt
          cons=cr*(1.0-cr)
          if(cr.gt.0.25)cons=0.25

         grmas(i,j,k)=0.5*cons*(1.0-max(0.0,min(2.0*rmas(i,j,k),1.)))
         grmen(i,j,k)=0.5*cons*(1.0-max(0.0,min(2.0*rmen(i,j,k),1.)))
        enddo
        enddo
        enddo

        do k=1,nz
        do j=1,ny
        do i=2,nx-1
        do m=1,nd
         tvd_x(i,j,k,m)=(grmas(i,j,k)+grmen(i+1,j,k))*dmas(i,j,k,m)-(grmas(i-1,j,k)+grmen(i,j,k))*dmen(i,j,k,m)
        enddo
        enddo
        enddo
        enddo
        
        do k=1,nz
        do j=1,ny
        do m=1,nd
         tvd_x(1,j,k,m)=(grmas(i,j,k)+grmen(i+1,j,k))*dmas(i,j,k,m)-(grmas(i,j,k)+grmen(i,j,k))*dmen(i,j,k,m)
         tvd_x(nx,j,k,m)=(grmas(i,j,k)+grmen(i,j,k))*dmas(i,j,k,m)-(grmas(i-1,j,k)+grmen(i,j,k))*dmen(i,j,k,m)
        enddo
        enddo
        enddo

        do k=1,nz
        do j=2,ny-1
        do i=1,nx
        do m=1,nd
         dmas(i,j,k,m)=um1(i,j+1,k,m)-um1(i,j,k,m)
         dmen(i,j,k,m)=um1(i,j,k,m)-um1(i,j-1,k,m)
        enddo
        enddo
        enddo
        do k=1,nz
        do i=1,nx
        do m=1,nd
         dmas(i,1,k,m)=um1(i,2,k,m)-um1(i,1,k,m)
         dmas(i,ny,k,m)=dmas(i,ny-1,k,m)
         dmen(i,1,k,m)=dmen(i,2,k,m)
         dmen(i,ny,k,m)=um1(i,ny,k,m)-um1(i,ny-1,k,m)
        enddo
        enddo
        enddo

        do k=1,nz
        do j=1,ny
        do i=1,nx
         rmas(i,j)=(dmen(i,j,k,1)*dmas(i,j,k,1)+dmen(i,j,k,2)*dmas(i,j,k,2)+dmen(i,j,k,3)*dmas(i,j,k,3))/ &
              (dmas(i,j,k,1)*dmas(i,j,k,1)+dmas(i,j,k,2)*dmas(i,j,k,2)+dmas(i,j,k,3)*dmas(i,j,k,3))
         rmen(i,j)=(dmen(i,j,k,1)*dmas(i,j,k,1)+dmen(i,j,k,2)*dmas(i,j,k,2)+dmen(i,j,k,3)*dmas(i,j,k,3))/ &
              (dmen(i,j,k,1)*dmen(i,j,k,1)+dmen(i,j,k,2)*dmen(i,j,k,2)+dmen(i,j,k,3)*dmen(i,j,k,3))
        enddo
        enddo
        do i=1,nx
        do j=1,ny
          c=SQRT(grav*h(i,j))
          cr = dy(i,j)*(ABS(u(i,j,2))+c)*dt
          cons=cr*(1.0-cr)
          if(cr.gt.0.25)cons=0.25

         grmas(i,j)=0.5*cons*(1.0-max(0.0,min(2.0*rmas(i,j),1.)))
         grmen(i,j)=0.5*cons*(1.0-max(0.0,min(2.0*rmen(i,j),1.)))
        enddo
        enddo

        do i=1,nx
        do j=2,ny-1
        do k=1,3
         tvd_y(i,j,k)=(grmas(i,j)+grmen(i,j+1))*dmas(i,j,k)-(grmas(i,j-1)+grmen(i,j))*dmen(i,j,k)
        enddo
        enddo
        enddo

        do i=1,nx
        do k=1,3
         tvd_y(i,1,k)=(grmas(i,j)+grmen(i,j+1))*dmas(i,j,k)-(grmas(i,j)+grmen(i,j))*dmen(i,j,k)
         tvd_y(i,ny,k)=(grmas(i,j)+grmen(i,j))*dmas(i,j,k)-(grmas(i,j-1)+grmen(i,j))*dmen(i,j,k)
        enddo
        enddo


        do i=1,nx
        do j=1,ny
        do k=1,3
         um1(i,j,k)=um1(i,j,k)+tvd_x(i,j,k)+tvd_y(i,j,k)
        enddo
        enddo
        enddo

        RETURN
        END SUBROUTINE TVD
