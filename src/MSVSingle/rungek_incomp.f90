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
       IMPLICIT NONE
       integer i,j,k,l,m,agrava,avieu
       real xm,ue

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
         write(6,101)'it=',it,' ','dt=',dt,' ','dtotal=',dto
!______________________________________________________________________
!  1
         write(6,*)'1'
          write(6,*)dt,cs,nc
         do m=1,nd
         do k=1,nz
         do j=1,ny
         do i=1,nx
          um(i,j,k,m)=um0(i,j,k,m)
         enddo
         enddo
         enddo
         enddo

         CALL deriv_vel()
         CALL viscosidad_t()
         CALL viscosidad()
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

          CALL SOURCE (1,m)
         ENDDO

         do k=1,nz
         do j=1,ny
         do i=1,nx
          if(um1(i,j,k,5).gt.1.+20./273.15)um1(i,j,k,5)=1.+20./273.15
          if(um1(i,j,k,5).lt.1.-20./273.15)um1(i,j,k,5)=1.-20./273.15
         enddo
         enddo
         enddo

 

!         CALL frontera_z_val (1)
         CALL frontera_y_val (1)
         CALL frontera_x_val (1)

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

          CALL SOURCE (2,m)
         ENDDO
         
         do k=1,nz
         do j=1,ny
         do i=1,nx
          if(um1(i,j,k,5).gt.1.+20./273.15)um1(i,j,k,5)=1.+20./273.15
          if(um1(i,j,k,5).lt.1.-20./273.15)um1(i,j,k,5)=1.-20./273.15
         enddo
         enddo
         enddo



!         CALL frontera_z_val (2)
         CALL frontera_y_val (2)
         CALL frontera_x_val (2)


        

!         if(it.eq.1)CALL FILTRADO(itf,mncit)
         IF(itf2.eq.25) THEN
         Write(6,*) 'FILTRADO' , itf
!         CALL FILTRADO(itf,mncit)
         itf=itf+1
         if(itf.eq.7) itf=1
         itf2=1
         ELSE
         itf2=itf2+1
         ENDIF

         CALL varc_to_varnc()

         do m=1,nd
         do k=1,nz
         do j=1,ny
         do i=1,nx
          um0(i,j,k,m)=um1(i,j,k,m)
         enddo
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
         
         dto=dto+dt
         it=it+1
        ENDDO
        CALL GRAVA ()

        RETURN
 101    FORMAT(a3,i6,a1,a3,e12.4,a1,a7,e12.4)
        END SUBROUTINE RUNGEK

