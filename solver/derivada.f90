
      SUBROUTINE inideriv()
      use consderper
      use consdernper
      use consrs
      use dimensiones
      IMPLICIT NONE
      integer i,j,k,l,m

      do i=1,nx
       bxnp(i)=1.0
      enddo

      do i=1,ny
       bynp(i)=1.0
      enddo

      do i=1,nz
       bznp(i)=1.0
      enddo

      axnp(1)=   0.0
      axnp(nx)=  3.0

      do i=2,nx-1
       axnp(i)= 1./4.
      enddo

      aynp(1)= 0.0
      aynp(ny)=3.0

      do i=2,ny-1
       aynp(i)=1./4.
      enddo

      aznp(1)= 0.0
      aznp(nz)=3.0

      do i=2,nz-1
       aznp(i)=1./4.
      enddo

      cxnp(1)= 3.0
      cxnp(nx)=0.0

      do i=2,nx-1
       cxnp(i)=1./4.
      enddo

      cynp(1)= 3.0
      cynp(ny)=0.0

      do i=2,ny-1
       cynp(i)=1./4.
      enddo

      cznp(1)= 3.0
      cznp(nz)=0.0

      do i=2,nz-1
       cznp(i)=1./4.
      enddo

!diagonales

      do i=1,nx
       bxp(i)=1.0
       axp(i)=1.0/4.0
       cxp(i)=1.0/4.0
      enddo
      do i=1,ny
       byp(i)=1.0
       ayp(i)=1.0/4.0
       cyp(i)=1.0/4.0
      enddo
      do i=1,nz
       bzp(i)=1.0
       azp(i)=1.0/4.0
       czp(i)=1.0/4.0
      enddo

!constantes nodos internos (4to orden)
     ars=3./4.
     brs=0.0

!constantes frontera
     ars1=-17./6.
     brs1=  3./2.
     crs1=  1.5


     return
     end subroutine inideriv

      SUBROUTINE derivnper(n,u,a,b,c,delta)
      use consrs
      IMPLICIT NONE
      integer i,j,k,l,m
      integer :: n
      real :: delta,bet
      real, dimension (n) ::  u,r,gam,a,b,c


        DO l=2,n-1
         r(l)=delta*(ars*(u(l+1)-u(l-1)))
        ENDDO

         r(1)=delta*(ars1*u(1)+brs1*u(2)+crs1*u(3)-1./6.*u(4))
         r(n)=delta*(-ars1*u(n)-brs1*u(n-1)-crs1*u(n-2)+1./6.*u(n-3))

        bet=b(1)
        u(1)=r(1)/bet

        do l=2,n
         gam(l)=c(l-1)/bet
         bet=b(l)-a(l)*gam(l)
         u(l)=(r(l)-a(l)*u(l-1))/bet
        enddo

        do l=n-1,1,-1
         u(l)=u(l)-gam(l+1)*u(l+1)
        enddo

        return
        end subroutine derivnper



