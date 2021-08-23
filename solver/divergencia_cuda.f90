      SUBROUTINE divergencia(neq)
      use dimensiones
      use derivtools
      use mflujos
      use jacobtools
      use dmflujos
      use right
      use consderper
      use consdernper
      use deltas
      use bob
      use variables

      IMPLICIT NONE
      INTEGER :: neq
      integer i,j,k,l,m



!----------------------------------------------
      interface
        subroutine cuDivergence(rs, rsv,                &   
                                e, f, g, ev, fv, gv,    &   
                                jbn,                    &
                                deltax, deltay, deltaz, & 
                                Nx, Ny, Nz)bind(C, name='cuDivergence')
        use iso_c_binding
        IMPLICIT NONE
        real(c_double) :: rs(Nx,Ny,Nz)
        real(c_double) :: rsv(Nx,Ny,Nz)
        real(c_double) :: e(Nx,Ny,Nz)
        real(c_double) :: f(Nx,Ny,Nz)
        real(c_double) :: g(Nx,Ny,Nz)
        real(c_double) :: ev(Nx,Ny,Nz)
        real(c_double) :: fv(Nx,Ny,Nz)
        real(c_double) :: gv(Nx,Ny,Nz)
        real(c_double) :: jbn(Nx,Ny,Nz,11)
        real(c_double),value :: deltax, deltay, deltaz
        integer(c_int),value :: Nx, Ny, Nz
        end subroutine
 
      end interface

!----------------------------------------------


        call cuDivergence(rs, rsv,                             &   
                          e, f, g, ev, fv, gv,                 &   
                          jbn,                                 &
                          1.0/ deltax, 1.0/deltay, 1.0/deltaz, & 
                          Nx, Ny, Nz)


      end subroutine divergencia
