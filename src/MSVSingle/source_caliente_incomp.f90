
!____________________________________________________________________
      SUBROUTINE SOURCE (irk,neq)
!____________________________________________________________________

      use dimensiones
      use variables
      use fronterapl
      use tiempo
      use consadim
      use velocidades
      use mallagrid
!      use frontera_uval
      use ranaleo
      use flotacion
      use jacobtools
      use derivvel
      use acoustic
      use combustion
      use sgdmodel
      use source_calor
      use flow

      implicit none
      integer neq,irk,i,j,k,l,m,itc
      real :: pe,te,re,eps,bruit,ee,ce
      real :: ue,ve,we,ut,vt,wt,pt,ct
      real deltat,ds,zc,xc,xa,A,udif,vort,cmean,xm,vol
      real :: rm,dm,area,xnus1,xnus2,xnus3,xnus4,amean,theta
      real :: dtin,dtfin,dtotal,angulo,pm


       IF(irk.eq.1)THEN
        deltat=dt
       ELSEIF(irk.eq.2)THEN
        deltat=0.5*dt
       ENDIF

       if(neq.eq.3) then
        rm=0.0
        vol=0.0
        dm=0.0
        tm=0.0
       Do j=2,ny-1
       Do k=2,nz-1
       Do i=2,nx-1
        tm=tm+temp(i,j,k)*dxyzsgd(i,j,k)
        vol=vol+dxyzsgd(i,j,k)
       ENDDO
       ENDDO
       ENDDO
        tm=tm/vol


       Do k=1,nz
       Do j=2,ny-1
       Do i=2,nx-1
         A=(tm-temp(i,j,k))
         um1(i,j,k,neq)=um1(i,j,k,neq)-deltat*(A)*froude
       Enddo
       Enddo
       Enddo
       ENDIF

!       IF(irk.eq.2)THEN

!        vort=0.0
!        do k=1,nz
!        do j=1,ny
!        do i=1,nx
!         vort=vort+(((jbnm(i,j,k,5)*dcvel(i,j,k,6)-jbnm(i,j,k,9)*dcvel(i,j,k,8))**2+ &
!                   (jbnm(i,j,k,9)*dcvel(i,j,k,7)-jbnm(i,j,k,1)*dcvel(i,j,k,3))**2+ &
!                (jbnm(i,j,k,1)*dcvel(i,j,k,2)-jbnm(i,j,k,5)*dcvel(i,j,k,4))**2)**0.5)* &
!                dxyzsgd(i,j,k)
!         vol=vol+dxyzsgd(i,j,k)
!        enddo
!        enddo
!        enddo
!        write(76,100)dto,vort/vol,rm,tm,amean/vol
!
!        xnus1=0.
!        xnus2=0.
!        area=0.0
!        do k=1,nz
!        do j=1,ny
!         xnus1=xnus1+((conc(1,j,k)-conc(2,j,k))/x(2,j,k,1))*dxyzsgd(1,j,k)
!         xnus2=xnus2+((conc(nx,j,k)-conc(nx-1,j,k))/(x(nx,j,k,1)-x(nx-1,j,k,1)))  &
!               *dxyzsgd(nx,j,k)
!         area=area+dxyzsgd(nx,j,k)
!        enddo
!        enddo
!        xnus3=0.
!        xnus4=0.
!        do k=1,nz
!        do i=1,nx
!         xnus3=xnus3+((conc(i,1,k)-conc(i,2,k))/x(i,2,k,2))*dxyzsgd(i,1,k)
!         xnus4=xnus4+((conc(i,ny,k)-conc(i,ny-1,k))/(x(i,ny,k,2)-x(i,ny-1,k,2)))  &
!               *dxyzsgd(i,ny,k)
!        enddo
!        enddo
!        write(77,100)dto,xnus1/area,xnus2/area,xnus3/area,xnus4/area,alpha
 
         
!       ENDIF
       if(neq.eq.1) then

        tm=0.0
        vol=0.0
       Do j=2,ny-1
       Do k=2,nz-1
       Do i=2,nx-1
        tm=tm+um1(i,j,k,1)*dxyzsgd(i,j,k)
        vol=vol+dxyzsgd(i,j,k)
       ENDDO
       ENDDO
       ENDDO
        tm=tm/vol
        write(6,*)'Presion media=',tm
       Do i=1,nx
       Do k=1,nz
       Do j=1,ny
          um1(i,j,k,1)=um1(i,j,k,1)-deltat*(tm-1.)
       Enddo
       Enddo
       Enddo

        tm=0.0
        vol=0.0
       Do j=2,ny-1
       Do k=2,nz-1
       Do i=2,nx-1
        tm=tm+um1(i,j,k,1)*dxyzsgd(i,j,k)
        vol=vol+dxyzsgd(i,j,k)
       ENDDO
       ENDDO
       ENDDO
        tm=tm/vol
        write(6,*)'Presion media=',tm
       endif



 100   FORMAT(f16.8,f16.8,f16.8,f16.8,f16.8,f16.8)

       RETURN
       END SUBROUTINE SOURCE

