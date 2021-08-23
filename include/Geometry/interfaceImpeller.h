!-----Interfaz to an external C function

      interface
         subroutine impellerRead(ntri) bind(C, name ='impellerRead')
          use iso_c_binding
          INTEGER(C_INT) :: ntri
         end subroutine 
      end interface

      interface
         subroutine impellerSource(neighbor, n, vrtx1, vrtx2, vrtx3, n_t) bind(C, name ='impellerSource')
         use iso_c_binding
         IMPLICIT NONE
          INTEGER(C_INT) :: neighbor(n_t)
          REAL(C_DOUBLE) :: n(n_t)
          REAL(C_DOUBLE) :: vrtx1(n_t)
          REAL(C_DOUBLE) :: vrtx2(n_t)
          REAL(C_DOUBLE) :: vrtx3(n_t)
          INTEGER(C_INT), value :: n_t
          
         end subroutine 
      end interface


!______________________________________________________________________

