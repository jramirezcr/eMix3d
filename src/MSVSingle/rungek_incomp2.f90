       SUBROUTINE RUNGEK()
       use dimensiones
       use consadim
       use velocidades
       use variables
       use tiempo
       use right
       use jacobtools
       use fieldrec
       use vieu
       use fronterapl
       use mac_it
       use combustion
       use afluid
       use flow
       use maskvar
       use mix_para
       IMPLICIT NONE
       integer i,j,k,l,m,agrava,avieu,n
       real xm,ue,dasum,dasum2,divmax

        dto=0.0
!        dt=0.0001
        it=1
        itf=1
        itf2=1
        igrava=0
        ivieu=1
        dgrava=dgrava1
        tvieu=tvieu1
        i_sample=0
        CALL MAXVAL()
        CALL GRAVA ()
        mncit=1
        DO WHILE (dto.le.tmax)
         CALL tstep()
         write(6,101)'it=',it,' ','dt=',dt,' ','nc=',nc,' ','dtotal=',dto
!______________________________________________________________________
!  1
         write(6,*)'1'
!          write(6,*)dt,cs,nc
         do m=1,nd
         do k=1,nz
         do j=1,ny
         do i=1,nx
          um(i,j,k,m)=um0(i,j,k,m)
         enddo
         enddo
         enddo
         enddo

!*
         DO n=1,nc

!          INTEGRATION IN THE STRETCHED PSEUDOTIME
!         LA PRESION (um(i,j,k,1)) AVANZA EN EL TIEMPO 
!         PERO EL RESTO DE LAS VARIABLES SOLO SE CORRIGEN (um(i,j,k,m) m=2,3,4,5)
          write(6,*)'n=',n
          IF(n.ge.2)THEN
           do m=2,nd
            do k=1,nz
             do j=1,ny
              do i=1,nx
               um(i,j,k,m)=um1(i,j,k,m)
              enddo
             enddo
            enddo
           enddo
           do m=1,1
            do k=1,nz
             do j=1,ny
              do i=1,nx
               um0(i,j,k,m)=um1(i,j,k,m)
               um(i,j,k,m)=um1(i,j,k,m)
              enddo
             enddo
            enddo
           enddo
          ENDIF




         CALL deriv_vel()
         CALL viscosidad_t()
         CALL viscosidad()
!         CALL esponja()
!         CALL frontera_x()
!         CALL frontera_z()
!         CALL frontera_y()
         DO m=1,nd
          CALL flujos(m)
          CALL divergencia(m,1,mncit)
          DO k=1,nz
           DO j=1,ny
            DO i=1,nx
             um1(i,j,k,m)=um0(i,j,k,m)+  &
                          dt*(rs(i,j,k)+rsv(i,j,k))
            ENDDO
           ENDDO
          ENDDO

!          CALL FXSOURCE (1,m)
!          CALL FZSOURCE (1,m)
!          CALL FYSOURCE (1,m)
         ENDDO
        
         CALL cvel_imp(1)
          !CALL SOURCE (atheta,aphi,0,0)
          !CALL SOURCE (atheta,(1.0*aphi),3.1415926, (0.5*3.1415926))


 

         CALL frontera_z_val (1)
         CALL frontera_y_val (1)
         CALL frontera_x_val (1)
          DO k=1,nz
           DO j=1,ny
            DO i=1,nx
             um1(i,j,k,1)=min(1.0,um1(i,j,k,1))
             um1(i,j,k,1)=max(-1.0,um1(i,j,k,1))
             um1(i,j,k,2)=min(1.0,um1(i,j,k,2))
             um1(i,j,k,2)=max(-1.0,um1(i,j,k,2))
             um1(i,j,k,3)=min(1.0,um1(i,j,k,3))
             um1(i,j,k,3)=max(-1.0,um1(i,j,k,3))
             um1(i,j,k,4)=min(1.0,um1(i,j,k,4))
             um1(i,j,k,4)=max(-1.0,um1(i,j,k,4))
             ENDDO
           ENDDO
          ENDDO


          DO k=1,nz
           DO j=1,ny
            DO i=1,nx
             if(mask(i,j,k).EQ.0.0)then
               um1(i,j,k,2)=0.0
               um1(i,j,k,3)=0.0
               um1(i,j,k,4)=0.0
             endif
             ENDDO
            ENDDO
           ENDDO



         CALL varc_to_varnc()

       
         do m=1,nd
         do k=1,nz
         do j=1,ny
         do i=1,nx
          um(i,j,k,m)=um1(i,j,k,m)
         enddo
         enddo
         enddo
         enddo


!______________________________________________________________________
!  2
         write(6,*)'2'
         CALL deriv_vel()
         CALL viscosidad_t()
         CALL viscosidad()
!         CALL esponja()
!         CALL frontera_x()
!         CALL frontera_z()
!         CALL frontera_y()
         DO m=1,nd
          CALL flujos(m)
          CALL divergencia(m,2,mncit)
          DO k=1,nz
           DO j=1,ny
            DO i=1,nx
             um1(i,j,k,m)=0.5*(um(i,j,k,m)+  &
                          um0(i,j,k,m))+ &
                          0.5*dt*(rs(i,j,k)+rsv(i,j,k))
            ENDDO
           ENDDO
          ENDDO

!          CALL FXSOURCE (2,m)
!          CALL FZSOURCE (2,m)
!          CALL FYSOURCE (2,m)
         ENDDO

         CALL cvel_imp(1)
         ! CALL SOURCE (atheta,aphi,0,0)
         ! CALL SOURCE (atheta,(1.0*aphi),3.1415926, (0.5*3.1415926))

         CALL frontera_z_val (2)
         CALL frontera_y_val (2)
         CALL frontera_x_val (2)
          DO k=1,nz
           DO j=1,ny
            DO i=1,nx
             um1(i,j,k,1)=min(1.0,um1(i,j,k,1))
             um1(i,j,k,1)=max(-1.0,um1(i,j,k,1))
             um1(i,j,k,2)=min(1.0,um1(i,j,k,2))
             um1(i,j,k,2)=max(-1.0,um1(i,j,k,2))
             um1(i,j,k,3)=min(1.0,um1(i,j,k,3))
             um1(i,j,k,3)=max(-1.0,um1(i,j,k,3))
             um1(i,j,k,4)=min(1.0,um1(i,j,k,4))
             um1(i,j,k,4)=max(-1.0,um1(i,j,k,4))
             ENDDO
           ENDDO
          ENDDO


          DO k=1,nz
           DO j=1,ny
            DO i=1,nx
             if(mask(i,j,k).EQ.0.0)then
               um1(i,j,k,2)=0.0
               um1(i,j,k,3)=0.0
               um1(i,j,k,4)=0.0
             endif
             ENDDO
            ENDDO
           ENDDO

        

!         if(it.eq.1)CALL FILTRADO(itf,mncit)
         IF(itf2.eq.1) THEN
         Write(6,*) 'FILTRADO' , itf
!         CALL FILTRADO(itf,mncit)
         itf=itf+1
         if(itf.eq.7) itf=1
         itf2=1
         ELSE
         itf2=itf2+1
         ENDIF

         CALL varc_to_varnc()

         ENDDO
!*
         do m=1,nd
         do k=1,nz
         do j=1,ny
         do i=1,nx
          um0(i,j,k,m)=um1(i,j,k,m)
         enddo
         enddo
         enddo
         enddo

         dasum=0.0
         dasum2=0.0
         divmax=0.0
         do k=1,nz
         do j=1,ny
         do i=1,nx
          dasum=dasum+((pres(i,j,k)-pres0(i,j,k))/pres(i,j,k))**2.
          dasum2=((pres(i,j,k)-pres0(i,j,k))/pres(i,j,k))**2.
          if(abs(dasum2).gt.divmax)divmax=dasum2
         enddo
         enddo
         enddo
          Write(6,*)'Diferencia presion Media  =',sqrt(dasum/float(nx*ny*nz))/(cs*cs*float(nc))
          Write(6,*)'Diferencia presion Maxima =',sqrt(divmax)/(cs*cs*float(nc))

         do k=1,nz
         do j=1,ny
         do i=1,nx
          pres0(i,j,k)=pres(i,j,k)
         enddo
         enddo
         enddo


         CALL ESTADISTICAS()

        mncit=mncit+1
        if(mncit.eq.9)mncit=1
        write(6,*)mncit

        avieu=int(dto/tvieu)
        IF(avieu.ge.1) THEN
         CALL MAXVAL()
         ivieu=ivieu+1
         tvieu=float(ivieu)*tvieu1
        endif

        agrava=int(dto/dgrava)
        IF(agrava.ge.1) THEN
        CALL GRAVA ()
        ENDIF
         
         dto=dto+sqrt(float(nc))*dt
         it=it+1
        ENDDO
        CALL GRAVA ()

        RETURN
 101    FORMAT(a3,i6,a1,a3,e12.4,a1,a3,i2,a1,a7,e12.4)
        END SUBROUTINE RUNGEK

