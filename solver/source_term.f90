
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
      use multiphase

      implicit none
      integer neq,irk,i,j,k,l,m,itc
      real :: pe,te,re,eps,bruit,ee,ce
      real :: ue,ve,we,ut,vt,wt,pt,ct
      real deltat,ds,zc,xc,xa,A,udif,vort,cmean,xm,vol
      real :: rm,dm,area,xnus1,xnus2,xnus3,xnus4,amean,theta
      real :: dtin,dtfin,dtotal,angulo,pm


       IF(irk.eq.1)THEN
         deltat = dt
       ELSEIF(irk.eq.2)THEN
         deltat = 0.25d0*dt
       ELSEIF(irk.eq.3)THEN
         deltat = (2.d0/3.d0)*dt
       ENDIF

       if(neq.eq.4) then

       write(*,*) 'EL FRUTEEEEEE',froude 

       Do k=2, nz-1
          Do j=2,ny-1
            Do i=2,nx-1
                  um1(i,j,k,neq)=um1(i,j,k,neq) - deltat*froude
            Enddo
           Enddo
       Enddo

       ENDIF



 100   FORMAT(f16.8,f16.8,f16.8,f16.8,f16.8,f16.8)

       RETURN
       END SUBROUTINE SOURCE

