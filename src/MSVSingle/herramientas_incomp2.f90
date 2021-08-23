!2_______________________________________________________________________
      SUBROUTINE varc_to_varnc()
!_______________________________________________________________________
      use variables
      use velocidades
      use dimensiones
      use consadim
      use combustion
      use flow
      IMPLICIT NONE
      integer :: i,j,k,l,m
      real :: A

         

          DO k=1,nz
           DO j=1,ny
            DO i=1,nx
             conc(i,j,k)=um1(i,j,k,5)
             temp(i,j,k)=um1(i,j,k,5)
             pres(i,j,k)=um1(i,j,k,1)/(Mach*Mach)
             u(i,j,k,1)=um1(i,j,k,2)/um1(i,j,k,1)
             u(i,j,k,2)=um1(i,j,k,3)/um1(i,j,k,1)
             u(i,j,k,3)=um1(i,j,k,4)/um1(i,j,k,1)
            ENDDO
           ENDDO
          ENDDO

       RETURN
       END SUBROUTINE varc_to_varnc
!_______________________________________________________________________
      SUBROUTINE varnc_to_varc()
!_______________________________________________________________________
      use variables
      use velocidades
      use dimensiones
      use consadim
      use combustion
      use flow
      IMPLICIT NONE
      integer :: i,j,k,l,m
      real :: A

       DO k=1,nz
        DO j=1,ny
         DO i=1,nx
          um0(i,j,k,1)=pres(i,j,k)*Mach*Mach
          um0(i,j,k,5)=temp(i,j,k)
          um0(i,j,k,2)=u(i,j,k,1)*um0(i,j,k,1)
          um0(i,j,k,3)=u(i,j,k,2)*um0(i,j,k,1)
          um0(i,j,k,4)=u(i,j,k,3)*um0(i,j,k,1)
          pres0(i,j,k)=pres(i,j,k)
         ENDDO
        ENDDO
       ENDDO

       RETURN
       END SUBROUTINE varnc_to_varc



!-----------------------------------------------------------------------
      subroutine scatter(A,B,n)
!-----------------------------------------------------------------------

      use dimensiones

      implicit none

      integer :: n
      real, dimension(n) :: A,B

      B=A

      end subroutine scatter

!-----------------------------------------------------------------------
      subroutine scatter3d(A,B)
!-----------------------------------------------------------------------

      use dimensiones

      implicit none

      real, dimension(nx,ny,nz) :: A,B

      B=A

      end subroutine scatter3d

!_______________________________________________________
      SUBROUTINE MAXVAL()
!_______________________________________________________
   
      use dimensiones
      use velocidades
      use variables

      IMPLICIT NONE
      real umax,umin,zzz
      integer imin,jmin,kmin,imax,jmax,kmax
      integer i,j,k,l,m

      umin=1000000000.
      umax=-1000000000.
      do k=1,nz
      do j=1,ny
      do i=1,nx
       zzz=u(i,j,k,1)
       if(zzz.gt.umax)then
        umax=zzz
        imax=i
        jmax=j
        kmax=k
       endif
       if(zzz.lt.umin)then
        umin=zzz
        imin=i
        jmin=j
        kmin=k
       endif
      enddo
      enddo
      enddo
      write(6,100)'umin',' ' ,umin,' ',imin,' ',jmin,' ',kmin
      write(6,100)'umax',' ' ,umax,' ',imax,' ',jmax,' ',kmax
      umin=1000000000.
      umax=-1000000000.
     do k=1,nz
      do j=1,ny
      do i=1,nx
       zzz=u(i,j,k,2)
       if(zzz.gt.umax)then
        umax=zzz
        imax=i
        jmax=j
        kmax=k
       endif
       if(zzz.lt.umin)then
        umin=zzz
        imin=i
        jmin=j
        kmin=k
       endif
      enddo
      enddo
      enddo
      write(6,100)'vmin',' ' ,umin,' ',imin,' ',jmin,' ',kmin
      write(6,100)'vmax',' ' ,umax,' ',imax,' ',jmax,' ',kmax
      umin=1000000000.
      umax=-1000000000.
     do k=1,nz
      do j=1,ny
      do i=1,nx
       zzz=u(i,j,k,3)
       if(zzz.gt.umax)then
        umax=zzz
        imax=i
        jmax=j
        kmax=k
       endif
       if(zzz.lt.umin)then
        umin=zzz
        imin=i
        jmin=j
        kmin=k
       endif
      enddo
      enddo
      enddo
      write(6,100)'wmin',' ' ,umin,' ',imin,' ',jmin,' ',kmin
      write(6,100)'wmax',' ' ,umax,' ',imax,' ',jmax,' ',kmax
      umin=1000000.
      umax=0.0
     do k=1,nz
      do j=1,ny
      do i=1,nx
       zzz=pres(i,j,k)
       if(zzz.gt.umax)then
        umax=zzz
        imax=i
        jmax=j
        kmax=k
       endif
       if(zzz.lt.umin)then
        umin=zzz
        imin=i
        jmin=j
        kmin=k
       endif
      enddo
      enddo
      enddo
      write(6,100)'pmin',' ' ,umin,' ',imin,' ',jmin,' ',kmin
      write(6,100)'pmax',' ' ,umax,' ',imax,' ',jmax,' ',kmax
      umin=1000000.
      umax=0.0
     do k=1,nz
      do j=1,ny
      do i=1,nx
       zzz=temp(i,j,k)
       if(zzz.gt.umax)then
        umax=zzz
        imax=i
        jmax=j
        kmax=k
       endif
       if(zzz.lt.umin)then
        umin=zzz
        imin=i
        jmin=j
        kmin=k
       endif
      enddo
      enddo
      enddo
      write(6,100)'tmin',' ' ,umin,' ',imin,' ',jmin,' ',kmin
      write(6,100)'tmax',' ' ,umax,' ',imax,' ',jmax,' ',kmax
      umin=1000000.
      umax=0.00000001
 100  format(a4,a2,e12.4,a2,i3,a2,i3,a2,i3) 
      return
      end subroutine maxval
     

!-----------------------------------------------------------------------
      subroutine filtrado(itf,nmcit)
!-----------------------------------------------------------------------

      use dimensiones
      use variables
      use velocidades
      use consdernper_f
      use consderper_f
      use derivtools
      use jacobtools
      use consadim
      use flow

      implicit none
      integer :: i,j,k,l,m,itf,nmcit,nx1,nx2,mum,nfc
      real :: pr,tm

!     Se va a filtrar con el esquema con que se derivo el corrector, para cada direccion
!     no se filtra C

       DO m=1,nd


        IF((itf.eq.1).or.(itf.eq.2)) THEN


        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           du(i)=um1(i,j,k,m)*jbn(i,j,k,10)
          ENDDO
         call filtronper(nx,du,axnpf,bxnpf,cxnpf)
         DO i=1,nx
          um1(i,j,k,m)=du(i)*jbn(i,j,k,11)
         ENDDO
        END DO
       END DO


        ELSEIF((itf.eq.3).or.(itf.eq.4)) THEN


        DO k=1,nz
         DO i=1,nx
          DO j=1,ny
           dv(j)=um1(i,j,k,m)*jbn(i,j,k,10)
          ENDDO
          call filtronper(ny,dv,aynpf,bynpf,cynpf)
          DO j=1,ny
           um1(i,j,k,m)=dv(j)*jbn(i,j,k,11)
          ENDDO
         END DO
        END DO


        ELSEIF((itf.eq.5).or.(itf.eq.6)) THEN


        DO j=1,ny
         DO i=1,nx
          DO k=1,nz
           dw(k)=um1(i,j,k,m)*jbn(i,j,k,10)
          ENDDO
          call filtronper(nz,dw,azpf,bzpf,czpf)
          DO k=1,nz
           um1(i,j,k,m)=dw(k)*jbn(i,j,k,11)
          ENDDO
         END DO
        END DO

        ENDIF

        IF((itf.eq.3).or.(itf.eq.5)) THEN


        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           du(i)=um1(i,j,k,m)*jbn(i,j,k,10)
          ENDDO
         call filtronper(nx,du,axnpf,bxnpf,cxnpf)
         DO i=1,nx
          um1(i,j,k,m)=du(i)*jbn(i,j,k,11)
         ENDDO
        END DO
       END DO
 

        ELSEIF((itf.eq.1).or.(itf.eq.6)) THEN


        DO k=1,nz
         DO i=1,nx
          DO j=1,ny
           dv(j)=um1(i,j,k,m)*jbn(i,j,k,10)
          ENDDO
          call filtronper(ny,dv,aynpf,bynpf,cynpf)
          DO j=1,ny
           um1(i,j,k,m)=dv(j)*jbn(i,j,k,11)
          ENDDO
         END DO
        END DO

        ELSEIF((itf.eq.2).or.(itf.eq.4)) THEN


        DO j=1,ny
         DO i=1,nx
          DO k=1,nz
          dw(k)=um1(i,j,k,m)*jbn(i,j,k,10)
         ENDDO
          call filtronper(nz,dw,azpf,bzpf,czpf)
          DO k=1,nz
           um1(i,j,k,m)=dw(k)*jbn(i,j,k,11)
          ENDDO
         END DO
        END DO

        ENDIF

        IF((itf.eq.6).or.(itf.eq.4)) THEN


        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           du(i)=um1(i,j,k,m)*jbn(i,j,k,10)
          ENDDO
         call filtronper(nx,du,axnpf,bxnpf,cxnpf)
         DO i=1,nx
          um1(i,j,k,m)=du(i)*jbn(i,j,k,11)
         ENDDO
        END DO
       END DO

        ELSEIF((itf.eq.2).or.(itf.eq.5)) THEN

        DO k=1,nz
         DO i=1,nx
          DO j=1,ny
           dv(j)=um1(i,j,k,m)*jbn(i,j,k,10)
          ENDDO
          call filtronper(ny,dv,aynpf,bynpf,cynpf)
          DO j=1,ny
           um1(i,j,k,m)=dv(j)*jbn(i,j,k,11)
          ENDDO
         END DO
        END DO

        ELSEIF((itf.eq.1).or.(itf.eq.3)) THEN


        DO j=1,ny
         DO i=1,nx
          DO k=1,nz
           dw(k)=um1(i,j,k,m)*jbn(i,j,k,10)
          ENDDO
          call filtronper(nz,dw,azpf,bzpf,czpf)
          DO k=1,nz
           um1(i,j,k,m)=dw(k)*jbn(i,j,k,11)
          ENDDO
         END DO
        END DO

       ENDIF

       END DO



       DO j=1,ny
        DO i=1,nx
         DO k=1,nz
          u(i,j,k,1)  = um1(i,j,k,2) / um1(i,j,k,1)
          u(i,j,k,2)  = um1(i,j,k,3) / um1(i,j,k,1)
          u(i,j,k,3)  = um1(i,j,k,4) / um1(i,j,k,1)
          temp(i,j,k) = um1(i,j,k,5)
          pres(i,j,k) = um1(i,j,k,1) / (Mach*Mach)
         ENDDO
        ENDDO
       ENDDO

       RETURN
       END SUBROUTINE FILTRADO

