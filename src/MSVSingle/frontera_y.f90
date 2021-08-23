!________________________________________________
       SUBROUTINE frontera_y ()
!________________________________________________
      use dimensiones
      use velocidades
      use variables
      use deltas
      use consadim
      use flow
      use fronterapl
      use jacobtools
      use derivtools
      use consderper
      use consdernper
      use derivvel
      use afluid

      IMPLICIT NONE
      integer :: i,j,k,l,m
      real tii,pii,uii,vii,wii,rii,uci
      real tii_1,pii_1,uii_1,vii_1,wii_1
      real tii_2,pii_2,uii_2,vii_2,wii_2
      real tii_3,pii_3,uii_3,vii_3,wii_3
      real tii_4,pii_4,uii_4,vii_4,wii_4
      real tif,pif,uif,vif,wif,rif,ucf
      real tif_1,pif_1,uif_1,vif_1,wif_1
      real tif_2,pif_2,uif_2,vif_2,wif_2
      real tif_3,pif_3,uif_3,vif_3,wif_3
      real tif_4,pif_4,uif_4,vif_4,wif_4
      real dudyi,drdyi,dvdyi,dpdyi,dwdyi,dtdyi
      real dudyf,drdyf,dvdyf,dpdyf,dwdyf,dtdyf
      real pl_l1,pl_l2,pl_l3,pl_l4,pl_l5
      real pl_d1,pl_d2,pl_d3,pl_d4,pl_d5
      real p_inf,xmachi,xmachf

!______________________________________________
!    z=1
!______________________________________________

      DO i=1,nx
       DO k=1,nz
        uii=u(i,1,k,1)
        vii=u(i,1,k,2)
        wii=u(i,1,k,3)
        pii=pres(i,1,k)
        tii=temp(i,1,k)
        rii=1.0

        uii_1=u(i,2,k,1)
        vii_1=u(i,2,k,2)
        wii_1=u(i,2,k,3)
        pii_1=pres(i,2,k)
        tii_1=temp(i,2,k)

        uii_2=u(i,3,k,1)
        vii_2=u(i,3,k,2)
        wii_2=u(i,3,k,3)
        pii_2=pres(i,3,k)
        tii_2=temp(i,3,k)

        uii_3=u(i,4,k,1)
        vii_3=u(i,4,k,2)
        wii_3=u(i,4,k,3)
        pii_3=pres(i,4,k)
        tii_3=temp(i,4,k)

        uii_4=u(i,5,k,1)
        vii_4=u(i,5,k,2)
        wii_4=u(i,5,k,3)
        pii_4=pres(i,5,k)
        tii_4=temp(i,5,k)

!        dudyi=jbnp(i,1,k,5)*deltay*((-25./12.)*uii+4.*uii_1-   &
!              3.*uii_2+(4./3.)*uii_3-(1./4.)*uii_4)
!        dvdyi=jbnp(i,1,k,5)*deltay*((-25./12.)*vii+4.*vii_1-   &
!              3.*vii_2+(4./3.)*vii_3-(1./4.)*vii_4)
!        dwdyi=jbnp(i,1,k,5)*deltay*((-25./12.)*wii+4.*wii_1-   &
!              3.*wii_2+(4./3.)*wii_3-(1./4.)*wii_4)
!        dpdyi=jbnp(i,1,k,5)*deltay*((-25./12.)*pii+4.*pii_1-   &
!              3.*pii_2+(4./3.)*pii_3-(1./4.)*pii_4)
!        dtdyi=jbnp(i,1,k,5)*deltay*((-25./12.)*tii+4.*tii_1-   &
!              3.*tii_2+(4./3.)*tii_3-(1./4.)*tii_4)
        dudyi=jbnp(i,1,k,5)*deltay*(uii_1-uii)
        dvdyi=jbnp(i,1,k,5)*deltay*(vii_1-vii)
        dwdyi=jbnp(i,1,k,5)*deltay*(wii_1-wii)
        dpdyi=jbnp(i,1,k,5)*deltay*(pii_1-pii)
        dtdyi=jbnp(i,1,k,5)*deltay*(tii_1-tii)
        drdyi=0.0
!        write(6,*)'i',dudyi,dvdyi,dwdyi,dpdyi,drdyi

        uif=u(i,ny,k,1)
        vif=u(i,ny,k,2)
        wif=u(i,ny,k,3)
        pif=pres(i,ny,k)
        tif=temp(i,ny,k)
        rif=1.0

        uif_1=u(i,ny-1,k,1)
        vif_1=u(i,ny-1,k,2)
        wif_1=u(i,ny-1,k,3)
        pif_1=pres(i,ny-1,k)
        tif_1=temp(i,ny-1,k)

        uif_2=u(i,ny-2,k,1)
        vif_2=u(i,ny-2,k,2)
        wif_2=u(i,ny-2,k,3)
        pif_2=pres(i,ny-2,k)
        tif_2=temp(i,ny-2,k)

        uif_3=u(i,ny-3,k,1)
        vif_3=u(i,ny-3,k,2)
        wif_3=u(i,ny-3,k,3)
        pif_3=pres(i,ny-3,k)
        tif_3=temp(i,ny-3,k)

        uif_4=u(i,ny-4,k,1)
        vif_4=u(i,ny-4,k,2)
        wif_4=u(i,ny-4,k,3)
        pif_4=pres(i,ny-4,k)
        tif_4=temp(i,ny-4,k)

!        dudyf=jbn(i,ny,k,5)*deltay*((25./12.)*uif-4.*uif_1+   &
!              3.*uif_2-(4./3.)*uif_3+(1./4.)*uif_4)
!        dvdyf=jbn(i,ny,k,5)*deltay*((25./12.)*vif-4.*vif_1+   &
!              3.*vif_2-(4./3.)*vif_3+(1./4.)*vif_4)
!        dwdyf=jbn(i,ny,k,5)*deltay*((25./12.)*wif-4.*wif_1+   &
!              3.*wif_2-(4./3.)*wif_3+(1./4.)*wif_4)
!        dpdyf=jbn(i,ny,k,5)*deltay*((25./12.)*pif-4.*pif_1+   &
!              3.*pif_2-(4./3.)*pif_3+(1./4.)*pif_4)
!        dtdyf=jbn(i,ny,k,5)*deltay*((25./12.)*tif-4.*tif_1+   &
!              3.*tif_2-(4./3.)*tif_3+(1./4.)*tif_4)
        dudyf=jbnm(i,ny,k,5)*deltay*(uif-uif_1)
        dvdyf=jbnm(i,ny,k,5)*deltay*(vif-vif_1)
        dwdyf=jbnm(i,ny,k,5)*deltay*(wif-wif_1)
        dpdyf=jbnm(i,ny,k,5)*deltay*(pif-pif_1)
        dtdyf=jbnm(i,ny,k,5)*deltay*(tif-tif_1)
        drdyf=0.0

!        write(6,*)'f',dudyf,dvdyf,dwdyf,dpdyf,drdyf
        uci=SQRT(uii*uii+cs*cs)
        ucf=SQRT(uif*uif+cs*cs)
        xmachi=mach*vii/uci
        xmachf=mach*vif/ucf
!__________________________________________________
!       PARED CON DESLIZAMIENTO (w=0)
!__________________________________________________
         pl_l1=(vii-uci)*(c3*dpdyi-rii*uci*dvdyi)
         pl_l5=pl_l1
         pl_l2=0.0
         pl_l3=0.0
         pl_l4=0.0

         pl_d1=(pl_l2+0.5*(pl_l1+pl_l5))/(uci*uci)
         pl_d2=0.5*(pl_l1+pl_l2)
         pl_d3=(pl_l5-pl_l1)/(2.*rii*uci)
         pl_d4=0.
         pl_d5=0.

!         write(6,*)'i',pl_l1,pl_d1,pl_d2,pl_d3

         fronty(1,1,i,k)=pl_d1*cs*cs
         fronty(1,2,i,k)=uif*pl_d1
         fronty(1,3,i,k)=0.0
         fronty(1,4,i,k)=wif*pl_d1
         fronty(1,5,i,k)=0.5*(uii*uii+wii*wii)*pl_d1+pl_d2/c2  

         pl_l5=(vif+ucf)*(c3*dpdyf+rif*ucf*dvdyf)
         pl_l1=pl_l5
         pl_l2=0.0
         pl_l3=0.0
         pl_l4=0.0

         pl_d1=(pl_l2+0.5*(pl_l1+pl_l5))/(ucf*ucf)
         pl_d2=0.5*(pl_l1+pl_l2)
         pl_d3=(pl_l5-pl_l1)/(2.*rif*ucf)
         pl_d4=0.
         pl_d5=0.
!         write(6,*)'f',pl_l1,pl_d1,pl_d2,pl_d3

         fronty(2,1,i,k)=pl_d1*cs*cs
         fronty(2,2,i,k)=uif*pl_d1
         fronty(2,3,i,k)=0.0
         fronty(2,4,i,k)=wif*pl_d1
         fronty(2,5,i,k)=0.5*(uii*uii+wii*wii)*pl_d1+pl_d2/c2

       ENDDO
      ENDDO
      RETURN
      END SUBROUTINE FRONTERA_Y

!____________________________________________________________________
      SUBROUTINE FYSOURCE (irk,neq)
!____________________________________________________________________

      use dimensiones
      use variables
      use fronterapl
      use tiempo
      implicit none
      integer neq,irk,i,j,k,l,m
      real deltat

       IF(irk.eq.1)THEN
        deltat=dt
       ELSEIF(irk.eq.2)THEN
        deltat=0.5*dt
       ENDIF

       Do k=2,nz-1
        Do i=2,nx
         um1(i,1,k,neq)=um1(i,1,k,neq)-deltat*fronty(1,neq,i,k)
         um1(i,ny,k,neq)=um1(i,ny,k,neq)-deltat*fronty(2,neq,i,k)
        Enddo
       Enddo

       RETURN
       END SUBROUTINE FYSOURCE

!__________________________________________________
       SUBROUTINE frontera_y_val ()
!__________________________________________________

       use dimensiones
       use velocidades
       use variables
       use consadim

       IMPLICIT NONE
       integer :: i,j,k,l,m
       real :: ue,ve,we,re,pe

       DO k=1,nz
       DO i=1,nx
        um1(i,1,k,1)=um1(i,2,k,1)
        um1(i,ny,k,1)=um1(i,ny-1,k,1)
        um1(i,1,k,2)=0.0!um1(i,2,k,2)
        um1(i,ny,k,2)=0.0!um1(i,ny-1,k,2)
        um1(i,1,k,3)=0.0
        um1(i,ny,k,3)=0.0
        um1(i,1,k,4)=0.0!um1(i,2,k,4)
        um1(i,ny,k,4)=0.0!um1(i,ny-1,k,4)
        um1(i,1,k,5)=um1(i,2,k,5)
        um1(i,ny,k,5)=um1(i,ny-1,k,5)
       ENDDO
       ENDDO
       RETURN
       END SUBROUTINE frontera_y_val
