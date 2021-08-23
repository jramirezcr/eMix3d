      SUBROUTINE INI_AFLUID()
      use acoustic
      use consadim
      use flow
      use tiempo
      use afluid
      IMPLICIT NONE
 
      dt=0.001
!      cs=cfl*dxyz_min/dt
      cs=sqrt(1./mach)     !3.0

      RETURN
      END SUBROUTINE INI_AFLUID






