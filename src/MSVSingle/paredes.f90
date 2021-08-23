      SUBROUTINE leer_pared()
      use dimensiones
      use paredes
      use inputname
      IMPLICIT NONE
      integer :: i,j,k,l,m
      integer :: in_Nvar,in_Nx,in_Ny,in_Nz


      write(6,*)'LECTURA PARED'
      open(11,file='pared.in',form='unformatted')
      read(11,err=1002)in_Nvar
      if (Nd.ne.in_Nvar) goto 1002
      read(11,err=1003)in_Nx,in_Ny,in_Nz
      if (Nx.ne.in_Nx) goto 1003
      if (Ny.ne.in_Ny) goto 1003
      if (Nz.ne.in_Nz) goto 1003


      read(11,err=1004) xmask
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

      SUBROUTINE pared()
      use dimensiones
      use paredes
      use inputname
      IMPLICIT NONE
      integer :: i,j,k,l,m

      Do k=1,nz
       Do j=1,ny
        Do i=1,nx
         um1(i,j,k,2)=xmask(i,j,k)*um1(i,j,k,2)
         um1(i,j,k,3)=xmask(i,j,k)*um1(i,j,k,3)
         um1(i,j,k,4)=xmask(i,j,k)*um1(i,j,k,4)
        Enddo
       Enddo
      Enddo

      RETURN
      END SUBROUTINE PARED

      SUBROUTINE leer_extractores()
      use dimensiones
      use paredes
      use inputname
      IMPLICIT NONE
      integer :: i,j,k,l,m
      integer :: in_Nvar,in_Nx,in_Ny,in_Nz


      write(6,*)'LECTURA EXTRACTORES'
      open(11,file='pared.in',form='unformatted')
      read(11,err=1002)in_Nvar
      if (Nd.ne.in_Nvar) goto 1002
      read(11,err=1003)in_Nx,in_Ny,in_Nz
      if (Nx.ne.in_Nx) goto 1003
      if (Ny.ne.in_Ny) goto 1003
      if (Nz.ne.in_Nz) goto 1003


      read(11,err=1004) xextract
      read(11,err=1004) yextract
      read(11,err=1004) zextract
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
      END SUBROUTINE LEER_EXTRACTORES

      SUBROUTINE leer_inyectores()
      use dimensiones
      use paredes
      use inputname
      IMPLICIT NONE
      integer :: i,j,k,l,m
      integer :: in_Nvar,in_Nx,in_Ny,in_Nz

      
      write(6,*)'LECTURA EXTRACTORES'
      open(11,file='pared.in',form='unformatted')
      read(11,err=1002)in_Nvar
      if (Nd.ne.in_Nvar) goto 1002
      read(11,err=1003)in_Nx,in_Ny,in_Nz
      if (Nx.ne.in_Nx) goto 1003
      if (Ny.ne.in_Ny) goto 1003
      if (Nz.ne.in_Nz) goto 1003

      
      read(11,err=1004) xinyect
      read(11,err=1004) yinyect
      read(11,err=1004) zinyect
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
      END SUBROUTINE LEER_INYECTORES

      SUBROUTINE pared()
      use dimensiones
      use paredes
      use inputname
      IMPLICIT NONE
      integer :: i,j,k,l,m

      Do k=1,nz
       Do j=1,ny
        Do i=1,nx
         um1(i,j,k,2)=xmask(i,j,k)*um1(i,j,k,2)
         um1(i,j,k,3)=xmask(i,j,k)*um1(i,j,k,3)
         um1(i,j,k,4)=xmask(i,j,k)*um1(i,j,k,4)
        Enddo
       Enddo
      Enddo

      RETURN
      END SUBROUTINE PARED



