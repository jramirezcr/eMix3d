!____________________________________________________________________
      SUBROUTINE SOURCE (irk)
!____________________________________________________________________

      use dimensiones
      use variables
      use fronterapl
      use tiempo
      use consadim
      use velocidades
      use mallagrid
!      use fuente
      use flow
      use ranaleo
      use jacobtools
      use derivvel
!      use perdidas
!      use turbulencia
!      use mascara
      use mezclador

      implicit none
      integer neq,irk,i,j,itc,nmcit,k
      real :: ue,ve,rly
 
!tiempo: 
!      temp=3.2    ![s]

!radio del Aspa:
      rini=0.1 

!radio de la Manivela:
      RMAY=0.3

!velocidad angular de la Manivela:
      atheta=0.5  ![rad/s]

!velocidad Angular del Aspa:
      aphi=-1.5   ![rad/s]

!Angulo de la Manivela y Aspa(dado por el tiempo y la velocidad angular):
      thetas=atheta*dto
      phis=aphi*dto



!--------------------------------------------------------------------
!              Calculo del Radio del Aspa para los nodos                        
!--------------------------------------------------------------------   

!////Singularidades/////

      ly=1.0
      dym=ly/94.

      DO i=1,nx
       DO j=1,ny
        rtot(i,j)=0.0
        rtotx(i,j)=0.0
        rtoty(i,j)=0.0
        rsub(i,j)=0.0
        rmin(i,j)=0.0
       enddo
      enddo

      write(6,*)'aqui',1
      DO i=1,nx-1
       DO j=1,ny-1
          IF((abs(2*rini*sin(phis)).LE.dym))then
          med=INT(((ly/2)+RMAY*sin(thetas))/dym)+1
          rtot(i,med)=(x(i,med,1,1)-(RMAY*cos(thetas)))/(cos(phis))
          IF(abs(rtot(i,med)).GE.rini)rtot(i,med)=0.0
          ENDIF
          IF(abs(2*rini*sin(phis)).GE.dym) then
          rmin(i,j)=(x(i,j,1,1)-(RMAY*cos(thetas)))/(cos(phis))
          rsub(i,j)=(x(i,j,1,2)-(RMAY*sin(thetas)))/(sin(phis))
          ENDIF
       ENDDO
      ENDDO 

      write(6,*)'aqui',2
!///////Generales///////
      
      IF(abs(2*rini*sin(phis)).GE.dym) then
      DO i=1,nx
       DO j=1,ny
          IF(abs(rmin(i,j)).GE.rini)rmin(i,j)=0.0
          IF(abs(rsub(i,j)).GE.rini)rsub(i,j)=0.0
       ENDDO
      ENDDO
      
      DO i=1,nx
       DO j=1,ny
         IF((rmin(i,j).NE.0.0).AND.(rsub(i,j).NE.0.0))Erro(i,j)=&
         abs((rmin(i,j)-rsub(i,j))/(rmin(i,j)))*100
         IF((rsub(i,j).NE.0.0).AND.(rmin(i,j).NE.0.0))Erroy(i,j)=&
         abs((rmin(i,j)-rsub(i,j))/(rsub(i,j)))*100
        ENDDO
      ENDDO
      
      DO i=1,nx
        cont=0
       DO j=1,ny
         IF((rmin(i,j).NE.0.0).AND.(rsub(i,j).NE.0.0))then
         IF(cont.EQ.0)minerro=Erro(i,j)
         IF(Erro(i,j).LE.minerro)then
           minerro=Erro(i,j)
           eb=j 
           cont=1
         ENDIF        
        ENDIF
       ENDDO
       rtotx(i,eb)=rmin(i,eb)   
      ENDDO

      DO j=1,ny
          cont=0
        DO i=1,nx
         IF((rmin(i,j).NE.0.0).AND.(rsub(i,j).NE.0.0))then
         IF(cont.EQ.0)minerro=Erro(i,j)
         IF(Erroy(i,j).LE.minerro)then
           minerro=Erroy(i,j)
           ea=i
           cont=1
          ENDIF
         ENDIF
        ENDDO
        rtoty(ea,j)=rsub(ea,j)
      ENDDO

      
      DO i=1,nx
         cont=0
         DO j=1,ny
          IF((rtotx(i,j).NE.0.0).AND.(cont.EQ.0).AND.(abs(cos(phis)).GE.cos(3.1415926*0.25))) then
          phimod=phis
          hm=RMAY*sin(thetas)+rtotx(i,j)*sin(phis)
          rtot(i,j)=rtotx(i,j)
          IF(abs(rtot(i,j)).GE.rini)rtot(i,j)=0.0
          rtot(i,j+1)=rtotx(i,j)-((hm-x(i,j+1,1,2))/cos(phimod))
          IF(abs(rtot(i,j+1)).GE.rini)rtot(i,j+1)=0.0
          rtot(i,j-1)=rtotx(i,j)-((hm-x(i,j-1,1,2))/cos(phimod))
          IF(abs(rtot(i,j-1)).GE.rini)rtot(i,j-1)=0.0
          cont=1
          ENDIF
         ENDDO
      ENDDO


      DO j=1,ny
         cont=0
         DO i=1,nx
          IF((rtoty(i,j).NE.0.0).AND.(cont.EQ.0.0).AND.(abs(sin(phis)).GE.sin(3.1415926*0.25))) then
          phimod=phis
          hm=RMAY*cos(thetas)+rtoty(i,j)*cos(phis)
          rtot(i,j)=rtoty(i,j)
          IF(abs(rtot(i,j)).GE.rini)rtot(i,j)=0.0
          rtot(i+1,j)=rtoty(i,j)-(hm-x(i+1,j,1,1))/sin(phis)
          IF(abs(rtot(i+1,j)).GE.rini)rtot(i+1,j)=0.0
          rtot(i-1,j)=rtoty(i,j)-(hm-x(i-1,j,1,1))/sin(phis)
          IF(abs(rtot(i-1,j)).GE.rini)rtot(i-1,j)=0.0
          cont=1
         ENDIF
         ENDDO
      ENDDO



      ENDIF
!--------------------------------------------------------------------
!                    Calculo de las velocidades  
!--------------------------------------------------------------------

      DO i=1,nx
       DO j=1,ny
        DO k=1,nz-1
          IF(rtot(i,j).NE.0.0)then
            ue=-1*((atheta*RMAY*sin(thetas))+(rtot(i,j)*aphi&
            *sin(phis)))
            um1(i,j,k,2)=ue
            ve=((atheta*RMAY*cos(thetas))+(rtot(i,j)*aphi*&
            cos(phis)))
            um1(i,j,k,3)=ve
            um1(i,j,k,1)=0.0
            um1(i,j,k,4)=0.0
            um1(i,j,k,5)=1.0
!           write(6,*)i,j,ue,ve
          ENDIF
       ENDDO
       ENDDO
      ENDDO
      DO i=1,nx
       DO j=1,ny
          IF(rtot(i,j).NE.0.0)then
           write(6,*)i,j,ue,ve
          ENDIF
       ENDDO
       ENDDO


      RETURN
      END SUBROUTINE SOURCE
