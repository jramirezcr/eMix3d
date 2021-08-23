!____________________________________________________________________
      SUBROUTINE SOURCE(omegaCarrusel, omegaImpulsor, desfase1, desfase2)
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
      use mix_para

      implicit none
      integer neq,irk,i,j,itc,nmcit,k
      real :: ue,ve,rly
      real omegaCarrusel, omegaImpulsor, desfase1, desfase2
 
!tiempo: 
!      temp=3.2    ![s]

!radio del Aspa:
!      rini=0.1 

!radio de la Manivela:
!      RMAY=0.3

!velocidad angular de la Manivela:
!      atheta=0.5  ![rad/s]

!velocidad Angular del Aspa:
!      aphi=-1.5   ![rad/s]

!Angulo de la Manivela y Aspa(dado por el tiempo y la velocidad angular):
      thetas=omegaCarrusel*(dto) + desfase1
      phis=omegaImpulsor*(dto)  + desfase2



!--------------------------------------------------------------------
!              Calculo del Radio del Aspa para los nodos                        
!--------------------------------------------------------------------   

      write(6,*)'aqui',1,rini,RMAY,aphi,atheta
      DO j=1,ny
       DO i=1,nx
         xx(i,j)=x(i,j,1,1)
         yy(i,j)=x(i,j,1,2)           
       ENDDO
      ENDDO
      write(6,*)'aqui',2
          
      DO i=1,nx-1
       DO j=1,ny-1
          rtot(i,j)=0.0
          rtoty(i,j)=0.0
          rtotx(i,j)=0.0
          rsub(i,j)=0.0
          rmin(i,j)=0.0
          rmin(i,j)=(xx(i,j)-(RMAY*cos(thetas)))/(cos(phis))
          rsub(i,j)=(yy(i,j)-(RMAY*sin(thetas)))/(sin(phis))
       ENDDO
      ENDDO

      IF(abs(cos(phis)).GE.cos(3.1415926*0.25))then
       DO i=1,nx
        cont=0
        DO j=1,ny
          IF(abs(rmin(i,j)).GE.rini)rmin(i,j)=0.0
          hm=RMAY*sin(thetas)+rmin(i,j)*sin(phis)
          IF((cont.EQ.0).AND.(hm.LE.yy(i,j)))then
          delm=abs(yy(i,j)-yy(i,j-1))
          IF(abs(hm-yy(i,j)).LE.delm)eb=j
          IF(abs(hm-yy(i,j)).GE.delm)eb=j-1
          cont=1
          ENDIF
         ENDDO
        rtotx(i,eb)=rmin(i,eb)
       ENDDO
      ENDIF
      write(6,*)'aqui',3

      IF(abs(sin(phis)).GE.sin(3.1415926*0.25))then
        DO j=1,ny
          cont=0
           DO i=1,nx
            IF(abs(rsub(i,j)).GE.rini)rsub(i,j)=0.0
            hm=RMAY*cos(thetas)+rsub(i,j)*cos(phis)
            IF((cont.EQ.0).AND.(hm.LE.xx(i,j)))then
            delm=abs(xx(i,j)-xx(i-1,j))
            IF(abs(hm-xx(i,j)).LE.delm)ea=i
            IF(abs(hm-xx(i,j)).GE.delm)ea=i-1
             cont=1
            ENDIF
           ENDDO
          rtoty(ea,j)=rsub(ea,j)
        ENDDO
       ENDIF
      write(6,*)'aqui',3.5

        DO i=1,nx
         DO j=1,ny
           IF(sin(phis).EQ.0.0)rtot(i,j)=rtotx(i,j)
           IF(cos(phis).EQ.0.0)rtot(i,j)=rtoty(i,j)
!           IF(rtotx(i,j).NE.0.0)posdos(i,j)=1.0
!           IF(rtoty(i,j).NE.0.0)posdos(i,j)=1.0
         ENDDO
        ENDDO

      write(6,*)'aqui',4

      DO i=1,nx
         cont=0
         DO j=1,ny
          IF((rtotx(i,j).NE.0.0).AND.(cont.EQ.0).AND.(abs(cos(phis)).GE.cos(3.1415926*0.25))&
          .AND.(sin(phis).NE.0.0)) then
          phimod=phis
          hm=RMAY*sin(thetas)+rtotx(i,j)*sin(phis)
          rtot(i,j)=rtotx(i,j)
           IF(abs(rtot(i,j)).GE.rini)rtot(i,j)=0.0
           rtot(i,j+1)=rtotx(i,j)-((hm-yy(i,j+1))/cos(phimod))
           IF(abs(rtot(i,j+1)).GE.rini)rtot(i,j+1)=0.0
           rtot(i,j-1)=rtotx(i,j)-((hm-yy(i,j-1))/cos(phimod))
          IF(abs(rtot(i,j-1)).GE.rini)rtot(i,j-1)=0.0
          cont=1
          ENDIF
         ENDDO
      ENDDO


      write(6,*)'aqui',5
      DO j=1,ny
         cont=0
         DO i=1,nx
          IF((rtoty(i,j).NE.0.0).AND.(cont.EQ.0.0).AND.(abs(sin(phis)).GE.sin(3.1415926*0.25)&
          ).AND.(cos(phis).NE.0.0)) then
          phimod=phis
          hm=RMAY*cos(thetas)+rtoty(i,j)*cos(phis)
          rtot(i,j)=rtoty(i,j)
           IF(abs(rtot(i,j)).GE.rini)rtot(i,j)=0.0
           rtot(i+1,j)=rtoty(i,j)-(hm-xx(i+1,j))/sin(phis)
           IF(abs(rtot(i+1,j)).GE.rini)rtot(i+1,j)=0.0
           rtot(i-1,j)=rtoty(i,j)-(hm-xx(i-1,j))/sin(phis)
           IF(abs(rtot(i-1,j)).GE.rini)rtot(i-1,j)=0.0
          cont=1
         ENDIF
         ENDDO
        ENDDO


      write(6,*)'aqui',6

!--------------------------------------------------------------------
!                    Calculo de las velocidades  
!--------------------------------------------------------------------

      DO i=1,nx
       DO j=1,ny
        DO k=1,nz-1

            IF(&
            ((rtot(i,j).NE.0.0).AND.(x(i,j,k,3).GE.lo-lasp).AND.&
             (x(i,j,k,3).LE.(lo-lasp+gaa))) .OR. & 
            ( (abs(rtot(i,j)).LE.rini).AND.&
             (abs(rtot(i,j)).GE.(rini-gat)).AND.(x(i,j,k,3).GE.(lo-lasp+gaa)))&
             )then 



       ue= (-1.0*((omegaCarrusel*RMAY*sin(thetas))+(rtot(i,j)*omegaImpulsor &
            *sin(phis)))) / Rossby

!            um1(i,j,k,2)=ue
      ve= (((omegaCarrusel*RMAY*cos(thetas))+(rtot(i,j)*omegaImpulsor* &
            cos(phis)))) / Rossby  
!            um1(i,j,k,3)=ve

            um1(i,j,k,2)=ue
            um1(i,j,k,3)=ve

            um1(i,j,k,1)=0.0
            um1(i,j,k,4)=0.0
            um1(i,j,k,5)=1.0



          ENDIF
       ENDDO
       ENDDO
      ENDDO
      write(6,*)'aqui',7
      DO i=1,nx
       DO j=1,ny
          IF(rtot(i,j).NE.0.0)then
!           write(6,*)i,j,ue,ve
          ENDIF
       ENDDO
       ENDDO


      RETURN
      END SUBROUTINE SOURCE
