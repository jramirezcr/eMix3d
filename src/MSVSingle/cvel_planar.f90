!__________________________________________________________________
       SUBROUTINE cvel_imp (irk)
!__________________________________________________________________

      use dimensiones
      use variables
      use tiempo
      use sgdmodel
      use source_calor
      use mallagrid

      use mix_para
      use mezclador

      implicit none
      integer neq,irk,i,j,itc,nmcit,k,no_imp
      real :: ue,ve

!Angulo de la Manivela y Aspa(dado por el tiempo y la velocidad angular):
!      thetas=atheta*dto
!      phis=aphi*dto

      DO j=1,ny
       DO i=1,nx
         xx(i,j)=x(i,j,2,1)
         yy(i,j)=x(i,j,2,2)
         rtot(i,j) = sqrt(xx(i,j)*xx(i,j)+yy(i,j)*yy(i,j))
       ENDDO
      ENDDO


!--------------------------------------------------------------------
!              Calculo del Radio del Aspa para los nodos
!--------------------------------------------------------------------
    



!--------------------------------------------------------------------
!                    Calculo de las velocidades
!--------------------------------------------------------------------

      DO i=1,nx
       DO j=1,ny
        DO k=1,nz-1

            IF(xx(i,j).NE.0.0)then
              thetas = atan(abs(yy(i,j)/xx(i,j)))
            ELSEIF(xx(i,j).EQ.0.0)then
              thetas = 3.1415926/2.
            ENDIF

             IF((xx(i,j).LT.0.0).AND.(yy(i,j).GT.0.0))then
              thetas = 3.141592 - thetas
             ELSEIF((xx(i,j).LT.0.0).AND.(yy(i,j).LE.0.0))then
               thetas = thetas + 3.141592
             ELSEIF((xx(i,j).GE.0.0).AND.(yy(i,j).LT.0.0))then
               thetas = 2.0*3.141592 - thetas
             ENDIF

            ue =-1.0*rtot(i,j)*atheta*sin(thetas)
            ve =rtot(i,j)*atheta*cos(thetas)

          IF((rtot(i,j).LE.(RMAY/(2.0*RMAY))).AND.(x(i,j,k,3).LE.(lo/(2.0*RMAY))))then

!            Segunda Adimensionalizacion

!            um1(i,j,k,2) = 2.0*3.1415926*ue / atheta
!            um1(i,j,k,3) = 2.0*3.1415926*ve / atheta

!           primera adimensionalizacion

             um1(i,j,k,2) = 2.0*ue /( atheta)
             um1(i,j,k,3) = 2.0*ve /(atheta)

!           tercera  adimensionalizacion
!             um1(i,j,k,2) = 2.0*ue
!             um1(i,j,k,3) = 2.0*ve

!            um1(i,j,k,1)=0.0
            um1(i,j,k,4)=0.0
          ENDIF
       ENDDO
       ENDDO
      ENDDO

      RETURN
      END SUBROUTINE CVEL_IMP

