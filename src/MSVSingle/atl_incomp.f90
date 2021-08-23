!-----------------------------------------------------------------------
      program ehecatl
!-----------------------------------------------------------------------

      use dimensiones
      use tiempo
      use deltas
      use ranaleo
      use flotacion

      IMPLICIT NONE

      idum=565878
      nd=5
      ne=5
      itr=20
      itconc=0.
      CALL LEER_DATOS ()
      CALL LEER_MIX ()
      CALL FLOW_ALLOC ()
      CALL INI_MAC ()
      CALL LEER_MALLA ()
      CALL LEER_MASCARA()
      CALL LEER_CAMPOS ()
       write(6,*)'aqui'
      CALL INIESTADISTICAS()

      CALL INIFRONTERA_X_VAL ()

      deltax=float(nx-1)
      deltay=float(ny-1)
      deltaz=float(nz-1)

      CALL INIDERIV ()
      CALL INIFILTRO ()
      CALL INITSTEP ()
      CALL INISGDM()
!      CALL JACOBEANO_COMP ()
      CALL JACOBEANO ()
 
      CALL INI_AFLUID()
      CALL INI_ESPONJA()
      CALL INI_CONSTANTES ()
      CALL RUNGEK ()
      CALL FINESTADISTICAS()
!      CALL tp_maker ()
      write(6,*)'aqui'
      CALL FLOW_DEALLOC ()

      PRINT *,'FIN DEL CALCULO EN EL TIEMPO',' ',dto,' ','ITERACION',' ',it
      END

