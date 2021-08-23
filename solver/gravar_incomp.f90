!______________________________________________________________
       SUBROUTINE GRAVA()
!______________________________________________________________

       use dimensiones
       use velocidades
       use fieldrec
       use multiphase
       use jacobtools
       use tiempo
       use maskvar
       IMPLICIT NONE
       character (20) outputfield
       real, allocatable, dimension (:,:,:) :: buffer
       integer i,j,k,l,m
       


       allocate(buffer(nx,ny,nz))
       WRITE(sample,'(i3.3)')i_sample
       WRITE(serie,'(i2.2)')n_serie
       outputfield='field_'//serie//'.'//sample
       write(6,*)'ESCRIBIENDO ARCHIVO','  ',outputfield
       OPEN(11,file=outputfield,form='unformatted')
       WRITE(11)nd
       WRITE(11)nx,ny,nz
       call scatter3d(u(:,:,:,1),buffer)
       WRITE(11)buffer
       call scatter3d(u(:,:,:,2),buffer)
       WRITE(11)buffer
       call scatter3d(u(:,:,:,3),buffer)
       WRITE(11)buffer
       WRITE(11)temp
       WRITE(11)pres
       CLOSE(11)

       outputfield='time_'//serie//'.'//sample
       write(6,*)'ESCRIBIENDO ARCHIVO','  ',outputfield
       OPEN(11,file=outputfield,form='unformatted')
       WRITE(11)dto, it
       CLOSE(11)

       outputfield='impellermask_'//serie//'.'//sample
       write(6,*)'ESCRIBIENDO ARCHIVO','  ',outputfield
       OPEN(11,file=outputfield,form='unformatted')
       WRITE(11)impellermask
       CLOSE(11)

       i_sample=i_sample+1
       igrava=igrava+1
       dgrava = float(int(dto/dgrava1))*dgrava1 + dgrava1
       RETURN
       END SUBROUTINE GRAVA
