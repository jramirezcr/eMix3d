      SUBROUTINE INI_CONSTANTES ()
      use dimensiones
      use consadim
      use flow
      use afluid
      IMPLICIT NONE

      c1=cs*cs
      c2=1.0
      c3=1.0
      c4=1./reynolds
      c5=1.0/(reynolds*eckert*prandtl)
      c6=1.0
      c7=1.0


      return
      end SUBROUTINE INI_CONSTANTES

      SUBROUTINE INI_MAC ()

      USE cons_mac

      cp1=-7./6.
      cp2=8./6.
      cp3=-1./6.

      cm1=7./6.
      cm2=-8./6.
      cm3=1./6.

      return
      end SUBROUTINE INI_MAC

