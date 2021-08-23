      SUBROUTINE leer_malla()
      use dimensiones
      use mallagrid
      use inputname
      IMPLICIT NONE
      integer :: i,j,k,l,m
      integer :: in_Nvar,in_Nx,in_Ny,in_Nz
      real, allocatable, dimension (:,:,:) :: buffer


      write(6,*)'LECTURA MALLA'
      open(11,file=inputgrid,form='unformatted')
      read(11,err=1002)in_Nvar
!      if (Nd.ne.in_Nvar) goto 1002
      read(11,err=1003)in_Nx,in_Ny,in_Nz
      if (Nx.ne.in_Nx) goto 1003
      if (Ny.ne.in_Ny) goto 1003
      if (Nz.ne.in_Nz) goto 1003

      allocate(buffer(nx,ny,nz))

      read(11,err=1004) buffer
      call scatter3d(buffer,x(:,:,:,1))
      read(11,err=1004) buffer
      call scatter3d(buffer,x(:,:,:,2))
      read(11,err=1004) buffer
      call scatter3d(buffer,x(:,:,:,3))

      deallocate(buffer)
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
      END SUBROUTINE LEER_MALLA



      SUBROUTINE leer_campos()
      use dimensiones
      use velocidades
      use variables
      use consadim
      use inputname
      use tiempo
      use multiphase
      IMPLICIT NONE
      integer :: i,j,k,l,m
      integer :: in_Nvar,in_Nx,in_Ny,in_Nz
      real, allocatable, dimension (:,:,:) :: buffer

      allocate(phi(nx,ny,nz))
      allocate(dmins(nx,ny,nz))
      allocate(hevi(nx,ny,nz))
      allocate(rhoTwoPhase(nx,ny,nz))
      allocate(muTwoPhase(nx,ny,nz))
      allocate(rhoint(nx,ny,nz))
      allocate(muint(nx,ny,nz))
      allocate(sobject(nx,ny,nz))

      write(6,*)'LECTURA CAMPOS'
      open(11,file=inputfield,form='unformatted')
      read(11,err=1002)in_Nvar
      if (Nd.ne.in_Nvar) goto 1002
      read(11,err=1003)in_Nx,in_Ny,in_Nz
      if (Nx.ne.in_Nx) goto 1003
      if (Ny.ne.in_Ny) goto 1003
      if (Nz.ne.in_Nz) goto 1003

      allocate(buffer(nx,ny,nz))

      read(11,err=1004) buffer
      call scatter3d(buffer,U(:,:,:,1))
      read(11,err=1004) buffer
      call scatter3d(buffer,U(:,:,:,2))
      read(11,err=1004) buffer
      call scatter3d(buffer,U(:,:,:,3))
      read(11,err=1004) buffer
      call scatter3d(buffer,temp(:,:,:))
      read(11,err=1004) buffer
      call scatter3d(buffer,pres(:,:,:))

      CALL varnc_to_varc ()
      deallocate(buffer)
      CLOSE (11)


      write(6,*)'LECTURA SOBJECT'
      open(11,file='sobject.field',form='unformatted')
      read(11,err=1003)in_Nx,in_Ny,in_Nz
      if (Nx.ne.in_Nx) goto 1003
      if (Ny.ne.in_Ny) goto 1003
      if (Nz.ne.in_Nz) goto 1003
      read(11,err=1004) sobject
      CLOSE (11)


      

      RETURN

 1002 Print*,in_Nvar,' leida ',Nd,' esperada'
      print*,'ERROR DETECTADO LEYENDO EN NUMERO DE VARIABLES,CAMPOS'
      stop
 1003 print*,in_Nx,in_Ny,in_Nz,' leidas ',Nx,Ny,Nz,' eperadas'
      print*,'ERROR DETECTADO LEYENDO EN DIMENSIONES, CAMPOS'
      stop
 1004 print*,'ERROR DETECTADO LEYENDO LOS CAMPOS'
      stop
      END SUBROUTINE LEER_CAMPOS

      SUBROUTINE leer_datos()
      use dimensiones
      use flow
      use inputname
      use tiempo
      use sgdmodel
      use fieldrec
      use vieu
      IMPLICIT NONE
      character (1) char1
      real cfilt


      open(10,file='data.in',form='formatted')

      read(10,'(a1)')char1
      read(10,101)Nx
      read(10,101)Ny
      read(10,101)Nz
      read(10,102)Reynolds
      read(10,102)Froude
      read(10,102)Prandtl
      read(10,107)sgs_model
      read(10,102)Prandtl_t
      read(10,103)inputgrid
      read(10,103)inputfield
      read(10,103)inputtime
      read(10,102)Mach
      read(10,102)gamma
      read(10,102)Tmax
      read(10,102)CFL
      read(10,102)Cvisc
      read(10,102)cfilt
      read(10,109)n_serie
      read(10,102)dgrava1
      read(10,102)tvieu1

        write(6,703)Reynolds,Prandtl,sgs_model,Prandtl_t &
                   ,inputgrid,inputfield,Mach,gamma &
                   ,Tmax,CFL,Cvisc,cfilt,n_serie,dgrava1, &
                    tvieu1 
 703  format('Reynolds number-------------: ',e10.4/ &
             'Prandtl number--------------: ',e10.4/ &
             'sgs model-------------------: ',a5/ &
             'Turbulent Prandtl number----: ',e10.4/ &
             'grid file-------------------: ',a32/ &
             'restart file----------------: ',a20/ &
             'Mach number-----------------: ',e10.4/ &
             'Gamma ----------------------: ',e10.4/ &
             'end time--------------------: ',e10.4/ &
             'CFL-------------------------: ',e10.4/ &
             'Cvisc-----------------------: ',e10.4/ &
             'Cfilt-----------------------: ',e10.4/ &
             'Serie-----------------------: ',i3/    &
             'T Gravado-------------------: ',e10.4/ &
             'T VIS. MIN/MAX--------------: ',e10.4)

 101  format(20x,i4)
 102  format(20x,e12.4)
 103  format(20x,a32)
 104  format(20x,i2)
 105  format(20x,i6)
 106  format(20x,a3)
 107  format(20x,a5)
 108  format(20x,a4)
 109  format(20x,i3)



      RETURN
      END


      SUBROUTINE leer_mix()

      use mix_para
      character char3

      OPEN(20,FILE='mixdat.in',FORM='formatted')
      READ(20,'(a1)') char3
      READ(20,110) rini
      READ(20,110) RMAY
      READ(20,110) aphi
      READ(20,110) atheta
      READ(20,110) lo
      READ(20,110) lasp
      READ(20,110) gat
      READ(20,110) gaa
      READ(20,110) rossby


110   FORMAT(20x,e12.4)
      CLOSE(20)
      RETURN

     END SUBROUTINE leer_mix

      SUBROUTINE leer_mascara()
      use maskvar

      IMPLICIT NONE
      integer mkx, mky, mkz, mkt 

      WRITE(6,*) 'LEYENDO MASCARA'
      OPEN(11,FILE = 'mez.mask',FORM ='unformatted')
      READ(11) mkt 
      WRITE(6,*) 'LEYENDO MASCARA'
      READ(11) mkx, mky ,mkz
      WRITE(6,*) 'LEYENDO MASCARA'
      READ(11) mask
      WRITE(6,*) 'LEYENDO MASCARA'
      CLOSE(11)

      ENDSUBROUTINE leer_mascara




      SUBROUTINE leer_time()
      use dimensiones
      use velocidades
      use variables
      use consadim
      use inputname
      use tiempo
      use multiphase
      IMPLICIT NONE

      write(6,*)'LECTURA TIME'
      open(11,file=inputtime,form='unformatted')
      read(11) dto, it
      CLOSE (11)
      write(6,*)'LECTURA TIME', dto, it

      ENDSUBROUTINE leer_time

      

