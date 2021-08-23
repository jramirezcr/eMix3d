!________________________________________________________
      SUBROUTINE INIESTADISTICAS
!________________________________________________________

      use dimensiones
      use variables
      use velocidades
      use stat
      IMPLICIT NONE
      integer i,j,k,l,m

       do k=1,nz
       do j=1,ny
       do i=1,nx
       do l=1,17
        st(i,j,k,l)=0.0
       enddo
       enddo
       enddo
       enddo

       RETURN
       END SUBROUTINE INIESTADISTICAS


!________________________________________________________
      SUBROUTINE ESTADISTICAS
!________________________________________________________

      use dimensiones
      use variables
      use velocidades
      use stat
      use tiempo
      use sgdmodel
      IMPLICIT NONE

      integer  i,j,k,l,m,icor


      do k=1,nz
      do j=1,ny
      do i=1,nx
       st(i,j,k,1)=st(i,j,k,1)+dt*u(i,j,k,1)
       st(i,j,k,2)=st(i,j,k,2)+dt*u(i,j,k,1)**2.
       st(i,j,k,3)=st(i,j,k,3)+dt*u(i,j,k,2)
       st(i,j,k,4)=st(i,j,k,4)+dt*u(i,j,k,2)**2.
       st(i,j,k,5)=st(i,j,k,5)+dt*u(i,j,k,3)
       st(i,j,k,6)=st(i,j,k,6)+dt*(u(i,j,k,3)*u(i,j,k,3))
       st(i,j,k,7)=st(i,j,k,7)+dt*(u(i,j,k,1)*u(i,j,k,2))
       st(i,j,k,8)=st(i,j,k,8)+dt*(u(i,j,k,1)*u(i,j,k,3))
       st(i,j,k,9)=st(i,j,k,9)+dt*(u(i,j,k,2)*u(i,j,k,3))
       st(i,j,k,10)=st(i,j,k,10)+dt*conc(i,j,k)
       st(i,j,k,11)=st(i,j,k,11)+dt*(temp(i,j,k)*u(i,j,k,1))
       st(i,j,k,12)=st(i,j,k,12)+dt*(temp(i,j,k)*u(i,j,k,2))
       st(i,j,k,13)=st(i,j,k,13)+dt*(temp(i,j,k)*u(i,j,k,3))
       st(i,j,k,14)=st(i,j,k,14)+dt*temp(i,j,k)**2.
       st(i,j,k,16)=st(i,j,k,16)+dt*pres(i,j,k)
       st(i,j,k,17)=st(i,j,k,17)+dt*pres(i,j,k)*pres(i,j,k)
      enddo
      enddo
      enddo

      RETURN
 110  FORMAT(f16.10,f16.10,f16.10,f16.10,f16.10,f16.10,f16.10,f16.10)
      END SUBROUTINE ESTADISTICAS

!________________________________________________________
      SUBROUTINE FINESTADISTICAS
!________________________________________________________

      use dimensiones
      use variables
      use velocidades
      use stat
      use tiempo
      IMPLICIT NONE
      integer i,j,k,l,m

       do k=1,nz
       do j=1,ny
       do i=1,nx
       do l=1,17
        st(i,j,k,l)=st(i,j,k,l)/(dto)
       enddo
       enddo
       enddo
       enddo
       OPEN(67,file='stat.data',form='unformatted')
       write(67)nx,ny,nz
       write(67)st
       CLOSE(67)

       RETURN
       END SUBROUTINE FINESTADISTICAS


