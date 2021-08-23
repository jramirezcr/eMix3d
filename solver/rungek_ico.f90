      SUBROUTINE RUNGEK()
      USE dimensiones
      USE consadim
      USE velocidades
      USE variables
      USE tiempo
      USE right
      USE jacobtools
      USE fieldrec
      USE vieu
      USE fronterapl
      USE mac_it
      USE combustion
      USE afluid
      USE flow
      USE maskvar
      USE mix_para
      USE multiphase
      USE sgdmodel
      USE mallagrid

      IMPLICIT NONE
      INTEGER :: i, j , k, l, m, agrava, avieu, n
      REAL    :: xm, ue, dasum, dasum2, divmax
      REAL    :: starttime, deltasmax, deltaprom, rad
      INTEGER :: sub_step


      itf    = 1
      itf2   = 1
      igrava = 0
      ivieu  = 1

      tvieu    = tvieu1
      i_sample = 0

      CALL MAXVAL()
      deltamin = minval(dmins)

      DO k=1,nz
         DO j=1,ny
            DO i=1,nx
               deltasmax = max(deltasmax,dxyzsgd(i,j,k))
            ENDDO
         ENDDO
      ENDDO
            
      deltaprom = (deltamin + deltasmax)*0.5d0
   
      call leer_time()
      starttime = dto

      dto = dto + dt
      it  = it + 1

      dgrava = float(int(dto/dgrava1))*dgrava1 + dgrava1

      CALL GRAVA ()
      mncit=1

      DO WHILE (dto.le.tmax)
         CALL tstep()
         write(6,101)'it=',it,' ','dt=',dt,' ','nc=',nc,' ','dtotal=',dto
!______________________________________________________________________
!  1

         write(6,*)'Runge Kutta substep 1'
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
            CALL divergencia(m)
            DO k=1,nz
               DO j=1,ny
                  DO i=1,nx
                     um1(i,j,k,m)=um0(i,j,k,m)+  &
                                  dt*(rs(i,j,k)+rsv(i,j,k))
                  ENDDO
               ENDDO
            ENDDO
            !CALL source(1,m)
         ENDDO
        

         CALL frontera_z_val (1)
         CALL frontera_y_val (1)
         CALL frontera_x_val (1)

         CALL rushton(1)

         CALL varc_to_varnc()

         DO k=1,nz
            DO j=1,ny
               DO i=1,nx
                  u(i,j,k,1)=min( 6.d0,u(i,j,k,1))
                  u(i,j,k,1)=max(-6.d0,u(i,j,k,1))
                  u(i,j,k,2)=min( 6.d0,u(i,j,k,2))
                  u(i,j,k,2)=max(-6.d0,u(i,j,k,2))
                  u(i,j,k,3)=min( 6.d0,u(i,j,k,3))
                  u(i,j,k,3)=max(-6.d0,u(i,j,k,3))

                  um1(i,j,k,2) = u(i,j,k,1)*um1(i,j,k,1)*mask(i,j,k)
                  um1(i,j,k,3) = u(i,j,k,2)*um1(i,j,k,1)*mask(i,j,k)
                  um1(i,j,k,4) = u(i,j,k,3)*um1(i,j,k,1)*mask(i,j,k)
               ENDDO
           ENDDO
         ENDDO


         CALL varc_to_varnc()
         CALL advectLevelSet(0)

         DO k=1,nz
            DO j=1,ny
               DO i=1,nx
                   um1(i,j,k,1) = pres(i,j,k)*Mach*Mach
                   um1(i,j,k,2) = u(i,j,k,1)*um1(i,j,k,1)  
                   um1(i,j,k,3) = u(i,j,k,2)*um1(i,j,k,1)
                   um1(i,j,k,4) = u(i,j,k,3)*um1(i,j,k,1)
             ENDDO
            ENDDO
          ENDDO
        
          um(:,:,:,:) = um1(:,:,:,:)


!______________________________________________________________________
!  2

         WRITE(6,*)'Runge Kutta substep 2'

         CALL deriv_vel()
         CALL viscosidad_t()
         CALL viscosidad()

         DO m=1,nd
            CALL flujos(m)
            CALL divergencia(m)
            um1(:,:,:,m)= 0.75d0*um0(:,:,:,m)  + 0.25d0*um(:,:,:,m) +  &
                          0.25d0*dt*(rs(:,:,:) + rsv(:,:,:))
            !CALL source(2,m)
         ENDDO


         CALL frontera_z_val (2)
         CALL frontera_y_val (2)
         CALL frontera_x_val (2)

         CALL rushton(2)

         CALL varc_to_varnc()

         DO k=1,nz
            DO j=1,ny
               DO i=1,nx
                  u(i,j,k,1)=min( 6.d0,u(i,j,k,1))
                  u(i,j,k,1)=max(-6.d0,u(i,j,k,1))
                  u(i,j,k,2)=min( 6.d0,u(i,j,k,2))
                  u(i,j,k,2)=max(-6.d0,u(i,j,k,2))
                  u(i,j,k,3)=min( 6.d0,u(i,j,k,3))
                  u(i,j,k,3)=max(-6.d0,u(i,j,k,3))

                  um1(i,j,k,2) = u(i,j,k,1)*um1(i,j,k,1)*mask(i,j,k)
                  um1(i,j,k,3) = u(i,j,k,2)*um1(i,j,k,1)*mask(i,j,k)
                  um1(i,j,k,4) = u(i,j,k,3)*um1(i,j,k,1)*mask(i,j,k)
               ENDDO
            ENDDO
         ENDDO

         CALL varc_to_varnc()
         CALL advectLevelSet(0)

         DO k=1,nz
            DO j=1,ny
               DO i=1,nx
                   um1(i,j,k,1) = pres(i,j,k)*Mach*Mach
                   um1(i,j,k,2) = u(i,j,k,1)*um1(i,j,k,1)  
                   um1(i,j,k,3) = u(i,j,k,2)*um1(i,j,k,1)
                   um1(i,j,k,4) = u(i,j,k,3)*um1(i,j,k,1)
               ENDDO
            ENDDO
         ENDDO

         um(:,:,:,:) = um1(:,:,:,:)


!____________________________________________________________________
!3

         WRITE(6,*)'Runge Kutta substep 3'

         CALL deriv_vel()
         CALL viscosidad_t()
         CALL viscosidad()

         DO m=1,nd
            CALL flujos(m)
            CALL divergencia(m)
            DO k=1,nz
               DO j=1,ny
                  DO i=1,nx

                     um1(i,j,k,m) = (1.d0/3.d0)*um0(i,j,k,m)  + &
                                    (2.d0/3.d0)*um(i,j,k,m)   + &  
                                    (2.d0/3.d0)*dt*(rs(i,j,k) + & 
                                    rsv(i,j,k))
                  ENDDO
               ENDDO
            ENDDO
         !   CALL source(3,m)
         ENDDO


         CALL frontera_z_val (3)
         CALL frontera_y_val (3)
         CALL frontera_x_val (3)

         CALL rushton(3)

         CALL varc_to_varnc()

         DO k=1,nz
            DO j=1,ny
               DO i=1,nx
                  u(i,j,k,1)=min( 6.d0,u(i,j,k,1))
                  u(i,j,k,1)=max(-6.d0,u(i,j,k,1))
                  u(i,j,k,2)=min( 6.d0,u(i,j,k,2))
                  u(i,j,k,2)=max(-6.d0,u(i,j,k,2))
                  u(i,j,k,3)=min( 6.d0,u(i,j,k,3))
                  u(i,j,k,3)=max(-6.d0,u(i,j,k,3))

                  um1(i,j,k,2) = u(i,j,k,1)*um1(i,j,k,1)*mask(i,j,k)
                  um1(i,j,k,3) = u(i,j,k,2)*um1(i,j,k,1)*mask(i,j,k)
                  um1(i,j,k,4) = u(i,j,k,3)*um1(i,j,k,1)*mask(i,j,k)
               ENDDO
            ENDDO
         ENDDO

         CALL varc_to_varnc()
         CALL advectLevelSet(0)

         DO k=1,nz
          DO j=1,ny
           DO i=1,nx
               um1(i,j,k,1) = pres(i,j,k)*Mach*Mach
               um1(i,j,k,2) = u(i,j,k,1)*um1(i,j,k,1)  
               um1(i,j,k,3) = u(i,j,k,2)*um1(i,j,k,1)
               um1(i,j,k,4) = u(i,j,k,3)*um1(i,j,k,1)
            ENDDO
           ENDDO
         ENDDO

         um(:,:,:,:)=um1(:,:,:,:)

         IF(itf.EQ.1) THEN
           Write(6,*) 'FILTRADO' , itf
           CALL filtrado( itf, mncit)
           itf = itf + 1
         ELSE
           itf = itf + 1
           if(itf.eq.10) itf = 1
         ENDIF

         um0(:,:,:,:)=um1(:,:,:,:)

!-------Runge Kutta ending

         DO k=1,nz
            DO j=1,ny
               DO i=1,nx
                  pres0(i,j,k)=pres(i,j,k)
               ENDDO
            ENDDO
         ENDDO


         CALL ESTADISTICAS()

         mncit =mncit + 1
         if(mncit.eq.9) mncit = 1
         write(6,*) mncit

         avieu = int((dto)/tvieu)

         IF(avieu.ge.1) THEN
           CALL MAXVAL()
           ivieu=ivieu+1
           tvieu=float(ivieu)*tvieu1
         ENDIF

         agrava=int((dto)/dgrava)
         IF(agrava.ge.1) THEN
         CALL GRAVA ()
         ENDIF
         
         dto = dto + sqrt(float(nc))*dt
         it = it + 1

       ENDDO
       !End Time Do-While    
       CALL GRAVA ()

       RETURN
 101   FORMAT(a3,i6,a1,a3,e12.4,a1,a3,i2,a1,a7,e12.4)
      END SUBROUTINE RUNGEK

