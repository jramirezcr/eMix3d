!__________________________________________________
       SUBROUTINE frontera_x ()
!__________________________________________________

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
      real dudxi,drdxi,dvdxi,dpdxi,dwdxi,dTdxi
      real dudxf,drdxf,dvdxf,dpdxf,dwdxf,dtdxf
      real p_inf,xmachi,xmachf
      real pl_l1,pl_l2,pl_l3,pl_l4,pl_l5
      real pl_d1,pl_d2,pl_d3,pl_d4,pl_d5

!______________________________________________
!    z=1
!______________________________________________

      DO k=1,nz
       DO j=1,ny

        uii=u(1,j,k,1)
        vii=u(1,j,k,2)
        wii=u(1,j,k,3)
        pii=pres(1,j,k)
        rii=1.0
        tii=temp(1,j,k)

        uii_1=u(2,j,k,1)
        vii_1=u(2,j,k,2)
        wii_1=u(2,j,k,3)
        pii_1=pres(2,j,k)
        tii_1=temp(2,j,k)

        uii_2=u(3,j,k,1)
        vii_2=u(3,j,k,2)
        wii_2=u(3,j,k,3)
        pii_2=pres(3,j,k)
        tii_2=temp(3,j,k)

        uii_3=u(4,j,k,1)
        vii_3=u(4,j,k,2)
        wii_3=u(4,j,k,3)
        pii_3=pres(4,j,k)
        tii_3=temp(4,j,k)

        uii_4=u(5,j,k,1)
        vii_4=u(5,j,k,2)
        wii_4=u(5,j,k,3)
        pii_4=pres(5,j,k)
        tii_4=temp(5,j,k)

!        dudxi=jbn(1,j,k,1)*deltax*((-25./12.)*uii+4.*uii_1-   &
!              3.*uii_2+(4./3.)*uii_3-(1./4.)*uii_4)
!        dvdxi=jbn(1,j,k,1)*deltax*((-25./12.)*vii+4.*vii_1-   &
!              3.*vii_2+(4./3.)*vii_3-(1./4.)*vii_4)
!        dwdxi=jbn(1,j,k,1)*deltax*((-25./12.)*wii+4.*wii_1-   &
!              3.*wii_2+(4./3.)*wii_3-(1./4.)*wii_4)
!        dpdxi=jbn(1,j,k,1)*deltax*((-25./12.)*pii+4.*pii_1-   &
!              3.*pii_2+(4./3.)*pii_3-(1./4.)*pii_4)
!        dtdxi=jbn(1,j,k,1)*deltax*((-25./12.)*tii+4.*tii_1-   &
!              3.*tii_2+(4./3.)*tii_3-(1./4.)*tii_4)
        drdxi=0.0
        dudxi=jbnp(1,j,k,1)*deltax*(uii_1-uii)
        dvdxi=jbnp(1,j,k,1)*deltax*(vii_1-vii)
        dwdxi=jbnp(1,j,k,1)*deltax*(wii_1-wii)
        dpdxi=jbnp(1,j,k,1)*deltax*(pii_1-pii)
        dtdxi=jbnp(1,j,k,1)*deltax*(tii_1-tii)


        uif=u(nx,j,k,1)
        vif=u(nx,j,k,2)
        wif=u(nx,j,k,3)
        pif=pres(nx,j,k)
        tif=temp(nx,j,k)
        rif=1.0

        uif_1=u(nx-1,j,k,1)
        vif_1=u(nx-1,j,k,2)
        wif_1=u(nx-1,j,k,3)
        pif_1=pres(nx-1,j,k)
        tif_1=temp(nx-1,j,k)

        uif_2=u(nx-2,j,k,1)
        vif_2=u(nx-2,j,k,2)
        wif_2=u(nx-2,j,k,3)
        pif_2=pres(nx-2,j,k)
        tif_2=temp(nx-2,j,k)

        uif_3=u(nx-3,j,k,1)
        vif_3=u(nx-3,j,k,2)
        wif_3=u(nx-3,j,k,3)
        pif_3=pres(nx-3,j,k)
        tif_3=temp(nx-3,j,k)

        uif_4=u(nx-4,j,k,1)
        vif_4=u(nx-4,j,k,2)
        wif_4=u(nx-4,j,k,3)
        pif_4=pres(nx-4,j,k)
        tif_4=temp(nx-4,j,k)

!        dudxf=jbn(nx,j,k,1)*deltax*((25./12.)*uif-4.*uif_1+   &
!              3.*uif_2-(4./3.)*uif_3+(1./4.)*uif_4)
!        dvdxf=jbn(nx,j,k,1)*deltax*((25./12.)*vif-4.*vif_1+   &
!              3.*vif_2-(4./3.)*vif_3+(1./4.)*vif_4)
!        dwdxf=jbn(nx,j,k,1)*deltax*((25./12.)*wif-4.*wif_1+   &
!              3.*wif_2-(4./3.)*wif_3+(1./4.)*wif_4)
!        dpdxf=jbn(nx,j,k,1)*deltax*((25./12.)*pif-4.*pif_1+   &
!              3.*pif_2-(4./3.)*pif_3+(1./4.)*pif_4)
!        dtdxf=jbn(nx,j,k,1)*deltax*((25./12.)*tif-4.*tif_1+   &
!              3.*tif_2-(4./3.)*tif_3+(1./4.)*tif_4)
        dudxf=jbnm(nx,j,k,1)*deltax*(uif-uif_1)
        dvdxf=jbnm(nx,j,k,1)*deltax*(vif-vif_1)
        dwdxf=jbnm(nx,j,k,1)*deltax*(wif-wif_1)
        dpdxf=jbnm(nx,j,k,1)*deltax*(pif-pif_1)
        dtdxf=jbnm(nx,j,k,1)*deltax*(tif-tif_1)
        drdxf=0.0

        uci=SQRT(uii*uii+cs*cs)
        ucf=SQRT(uif*uif+cs*cs)
        xmachi=mach*uii/uci
        xmachf=mach*uif/ucf

!__________________________________________________
!       ENTRADA 
!__________________________________________________
        pl_l1=(uii-uci)*(c3*dpdxi-rii*uci*dudxi)
        pl_l5=pl_l1
        pl_l2=0.5*c2*(pl_l1+pl_l5)
        pl_l3=0.0
        pl_l4=0.0

        pl_d1=(pl_l2+0.5*(pl_l1+pl_l5))/(uci*uci)
        pl_d2=0.0
        pl_d2=0.0
        pl_d3=0.0
        pl_d4=0.0

        frontx(1,1,j,k)=pl_d1*cs*cs
        frontx(1,2,j,k)=0.0         
        frontx(1,3,j,k)=0.0
        frontx(1,4,j,k)=0.0
        frontx(1,5,j,k)=0.0

!__________________________________________________
!       SALIDA 
!__________________________________________________

        p_inf=1.00
        pl_l5=(uif+ucf)*(c3*dpdxf+rif*ucf*dudxf)
        pl_l1=.1*(1.-xmachf*xmachf)*ucf*(pif-p_inf)*c3
        pl_l2=uif*(ucf*ucf*drdxf-c3*dpdxf)
        pl_l3=uif*dvdxf
        pl_l4=uif*dwdxf

        pl_d1=(pl_l2+0.5*(pl_l1+pl_l5))/(ucf*ucf)
        pl_d2=0.5*(pl_l1+pl_l5)
        pl_d3=(pl_l5-pl_l1)/(2.*rif*ucf)
        pl_d4=pl_l3
        pl_d5=pl_l4

        frontx(2,1,j,k)=pl_d1*cs*cs
        frontx(2,2,j,k)=uif*pl_d1+pl_d3
        frontx(2,3,j,k)=vif*pl_d1+pl_d4
        frontx(2,4,j,k)=wif*pl_d1+pl_d5
        frontx(2,5,j,k)=tif*pl_d1+uif*dtdxf

       ENDDO
      ENDDO
      RETURN
      END SUBROUTINE FRONTERA_X

!____________________________________________________________________
      SUBROUTINE FXSOURCE (irk,neq)
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

       Do k=1,nz
        Do j=1,ny
         um1(1,j,k,neq)=um1(1,j,k,neq)-deltat*frontx(1,neq,j,k)
         um1(nx,j,k,neq)=um1(nx,j,k,neq)-deltat*frontx(2,neq,j,k)
        Enddo
       Enddo

       RETURN
       END SUBROUTINE FXSOURCE

!__________________________________________________
       SUBROUTINE frontera_x_val (irk)
!__________________________________________________

       use dimensiones
       use velocidades
       use variables
       use consadim
       use fronterapl
       use frontera_uval
       use ranaleo
       use tiempo
       IMPLICIT NONE
       integer :: i,j,k,l,m,irk
       real :: pe,te,re,eps,bruit

       eps=0.2
       do k=1,nz
       do j=1,ny
           um1(1,j,k,2)=0.0
           um1(1,j,k,3)=0.0 !um1(2,j,k,3)
           um1(1,j,k,4)=0.0!um1(2,j,k,4)
           um1(1,j,k,5)=um1(2,j,k,5)
           um1(1,j,k,1)=um1(2,j,k,1)

           um1(nx,j,k,2)=0.0
           um1(nx,j,k,3)=0.0!um1(nx-1,j,k,3)
           um1(nx,j,k,4)=0.0!um1(nx-1,j,k,4)
           um1(nx,j,k,5)=um1(nx-1,j,k,5)
           um1(nx,j,k,1)=um1(nx-1,j,k,1)
        enddo
        enddo
        return
        end subroutine frontera_x_val

!__________________________________________________
       SUBROUTINE inifrontera_x_val ()
!__________________________________________________

       use dimensiones
       use velocidades
       use fronterapl
       IMPLICIT NONE
       integer :: i,j,k,l,m
       integer :: in_Nx,in_Ny,in_Nz,in_Nvar

      write(6,*)'LECTURA CAMPO DE ENTRADA X'
      open(11,file='frontx.in',form='unformatted')
      read(11,err=1002)in_Nvar
      if (Nd.ne.in_Nvar) goto 1002
      read(11,err=1003)in_Nx,in_Ny,in_Nz
      if (Nx.ne.in_Nx) goto 1003
      if (Ny.ne.in_Ny) goto 1003
      if (Nz.ne.in_Nz) goto 1003

      read(11,err=1004) frontin

      CLOSE (11)
      RETURN

 1002 Print*,in_Nvar,' leida ',Nd,' esperada'
      print*,'ERROR DETECTADO LEYENDO EN NUMERO DE VARIABLES'
      stop
 1003 print*,in_Nx,in_Ny,in_Nz,' leidas ',Nx,Ny,Nz,' eperadas'
      print*,'ERROR DETECTADO LEYENDO EN DIMENSIONES'
      stop
 1004 print*,'ERROR DETECTADO LEYENDO LA MALLA'
      stop

!       do k=1,nz
!       do j=1,ny
!        frontin(j,k,1)=temp(1,j,k)
!        frontin(j,k,2)=u(1,j,k,1)
!        frontin(j,k,3)=u(1,j,k,2)
!        frontin(j,k,4)=u(1,j,k,3)
!        frontin(j,k,5)=pres(1,j,k)
!       enddo
!       enddo
!       return
       end subroutine inifrontera_x_val

!__________________________________________________
       SUBROUTINE inialea ()
!__________________________________________________

       use ifport

       call seed(1995)
 
       return
       end subroutine inialea
