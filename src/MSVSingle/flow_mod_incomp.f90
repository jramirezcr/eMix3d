!-----------------------------------------------------------------------
      module dimensiones
!-----------------------------------------------------------------------

      integer nx,ny,nz,nd,ne

      end  module dimensiones

!-----------------------------------------------------------------------
      module deltas
!-----------------------------------------------------------------------

      real :: deltax,deltay,deltaz

      end  module deltas

!----------------------------------------------------------------------
      module cons_mac
!-----------------------------------------------------------------------

      real cp1,cp2,cp3,cm1,cm2,cm3

      end module cons_mac

!----------------------------------------------------------------------
      module mac_it
!-----------------------------------------------------------------------

      integer mncit

      end module mac_it

!----------------------------------------------------------------------
      module jacobtools
!-----------------------------------------------------------------------

      real, allocatable, dimension (:,:,:,:) :: jbnp,jbnm,jbn
      real t1,t2,t3,t4,t5,t6,t7,t8,t9

      end module jacobtools

!----------------------------------------------------------------------
      module mallagrid
!-----------------------------------------------------------------------

      real, allocatable, dimension (:,:,:,:) :: x

      end module mallagrid

!----------------------------------------------------------------------
      module derivvel
!-----------------------------------------------------------------------

      real, allocatable, dimension (:,:,:,:) :: ddvelp,ddvelm,dcvel
      real, allocatable, dimension (:,:,:,:) :: dconcp,dconcm
      real, allocatable, dimension (:,:,:,:) :: dconczp,dconczm
      real, allocatable, dimension (:,:,:,:) :: dtempp,dtempm

      end module derivvel

!----------------------------------------------------------------------
      module mflujos
!-----------------------------------------------------------------------

      real, allocatable, dimension (:,:,:) :: e,f,g
      real, allocatable, dimension (:,:,:) :: ex,fx,gx
      real, allocatable, dimension (:,:,:) :: evx,fvx,gvx
      real, allocatable, dimension (:,:,:) :: evp,fvp,gvp
      real, allocatable, dimension (:,:,:) :: evm,fvm,gvm


      end module mflujos

!----------------------------------------------------------------------
      module dmflujos
!-----------------------------------------------------------------------

      real, allocatable, dimension (:,:,:) :: dere,derf,derg
      real, allocatable, dimension (:,:,:) :: derev,derfv,dergv

      end module dmflujos
!----------------------------------------------------------------------
      module viscosidades
!-----------------------------------------------------------------------

      real, allocatable, dimension (:,:,:) :: vis,vis2,vis6,visl,vist,vis5,vis5z

      end module viscosidades

!----------------------------------------------------------------------
      module consadim
!-----------------------------------------------------------------------

      real c1,c2,c3,c4,c5,c6,c7,c8,c9

      end module consadim

!----------------------------------------------------------------------
      module velocidades
!-----------------------------------------------------------------------

      real, allocatable, dimension (:,:,:,:) :: u
      real, allocatable, dimension (:,:,:) :: conc,temp,pres,concz

      end module velocidades

!----------------------------------------------------------------------
      module variables
!-----------------------------------------------------------------------

      real, allocatable, dimension (:,:,:,:) :: um,um1,um0

      end module variables

!----------------------------------------------------------------------
      module right
!-----------------------------------------------------------------------

      real, allocatable, dimension (:,:,:) :: rs,rsv

      end module right

!----------------------------------------------------------------------
      module tiempo
!-----------------------------------------------------------------------

      real, allocatable, dimension (:,:,:) :: dx,dy,dz 
      real :: dt,dto,olddt,cflm,cfl,diffdt,up,down,tmax,cvisc
      real :: cfla
      integer :: it,itf,itr,itf2

      end module tiempo

!----------------------------------------------------------------------
      module inputname
!-----------------------------------------------------------------------

      character(32) :: inputgrid,inputfield

      end module inputname
!----------------------------------------------------------------------
      module flow
!-----------------------------------------------------------------------

      real :: reynolds,froude,prandtl,prandtl_t,gamma,mach

      end module flow

!----------------------------------------------------------------------
      module sgdmodel
!-----------------------------------------------------------------------

      character(32) :: sgs_model
      real :: amu0,cutoff
      real, allocatable, dimension (:,:,:) :: amut,dxyzsgd
      real, allocatable, dimension (:,:,:) :: dxmsgd,dxpsgd
      real, allocatable, dimension (:,:,:) :: dymsgd,dypsgd
      real, allocatable, dimension (:,:,:) :: dzmsgd,dzpsgd
      integer :: isgm

      end module sgdmodel
!----------------------------------------------------------------------
      module vorticidad
!-----------------------------------------------------------------------

      real, allocatable, dimension (:,:,:) :: wx,wy,wz

      end module vorticidad

!----------------------------------------------------------------------
      module consource
!-----------------------------------------------------------------------

      real :: debit,xydim,velm,deltau,derho,chight,cener

      end module consource

!----------------------------------------------------------------------
      module stat
!-----------------------------------------------------------------------

      real, allocatable, dimension (:,:,:,:) :: st
      real, allocatable, dimension (:,:,:) :: corr_esp_u,corr_esp_v
      real, allocatable, dimension (:,:,:) :: corr_esp_w,corr_esp_t

      end module stat

!----------------------------------------------------------------------
      module fieldrec
!-----------------------------------------------------------------------

      integer :: n_serie,i_sample,igrava
      real :: dgrava,dgrava1
      character(2) serie
      character(3) sample

      end module fieldrec

!----------------------------------------------------------------------
      module vieu
!-----------------------------------------------------------------------

      real :: tvieu,tvieu1
      integer :: ivieu


      end module vieu
!----------------------------------------------------------------------
      module fronterapl
!-----------------------------------------------------------------------

      real, allocatable, dimension (:,:,:,:) :: frontx,frontz,fronty
      real, allocatable, dimension (:,:,:) :: frontin
      integer, allocatable, dimension (:,:,:) :: imask
      integer, allocatable, dimension (:,:) :: nent

      end module fronterapl
!----------------------------------------------------------------------
      module frontera_uval
!-----------------------------------------------------------------------
      real :: ue,ve,we

     end module frontera_uval

!----------------------------------------------------------------------
      module bob
!-----------------------------------------------------------------------

      real, allocatable, dimension (:,:,:) :: esplayer
      real, allocatable, dimension (:) :: esource
      integer :: ief,iei,iefi,ieii


      end module bob

!----------------------------------------------------------------------
      module bob_yz 
!-----------------------------------------------------------------------
      
      real, allocatable, dimension (:,:,:) :: esplayer_yz,esplayer_zy
      real, allocatable, dimension (:) :: esource_yz,esource_zy
      integer :: iefyz,ieiyz,iiyz,ifyz


      end module bob_yz

!----------------------------------------------------------------------
      module ranaleo
!-----------------------------------------------------------------------

      integer :: idum

      end module ranaleo


!----------------------------------------------------------------------
      module consderper_f
!-----------------------------------------------------------------------

      real, allocatable, dimension (:) :: axpf,aypf,azpf
      real, allocatable, dimension (:) :: bxpf,bypf,bzpf
      real, allocatable, dimension (:) :: cxpf,cypf,czpf

      end  module consderper_f
!----------------------------------------------------------------------
      module consdernper_f
!-----------------------------------------------------------------------

      real, allocatable, dimension (:) :: axnpf,aynpf,aznpf
      real, allocatable, dimension (:) :: bxnpf,bynpf,bznpf
      real, allocatable, dimension (:) :: cxnpf,cynpf,cznpf
      
      end  module consdernper_f
      
!----------------------------------------------------------------------
      module consrs_f
!-----------------------------------------------------------------------

      real :: arsf,brsf,crsf,drsf,ersf,frsf
      real :: arsf1,brsf1,crsf1,drsf1,ersf1
      real :: arsf2,brsf2,crsf2,drsf2
      real :: arsf3,brsf3,crsf3,drsf3,ersf3,frsf3,grsf3
      real :: arsf4,brsf4,crsf4,drsf4,ersf4,frsf4,grsf4

      end module consrs_f

!----------------------------------------------------------------------
      module consfiltro
!-----------------------------------------------------------------------
      
      real :: cfilt
      
      end module consfiltro
      
!----------------------------------------------------------------------
      module derivtools
!-----------------------------------------------------------------------
      
      real, allocatable, dimension (:) :: du,dv,dw
      
      end module derivtools

!______________________________________________________________________
      module acoustic
!-----------------------------------------------------------------------

      integer :: nlm
      real :: alpha2,casr1,casr2,alpha

      end module acoustic
!----------------------------------------------------------------------
      module flotacion
!-----------------------------------------------------------------------

      real :: pmean,tmean,rmean,itconc,dtconc
      real, allocatable, dimension (:) :: cm


      end module flotacion

!----------------------------------------------------------------------
      module consderper
!-----------------------------------------------------------------------

      real, allocatable, dimension (:) :: axp,ayp,azp
      real, allocatable, dimension (:) :: bxp,byp,bzp
      real, allocatable, dimension (:) :: cxp,cyp,czp

      end  module consderper
!----------------------------------------------------------------------
      module consdernper
!-----------------------------------------------------------------------

      real, allocatable, dimension (:) :: axnp,aynp,aznp
      real, allocatable, dimension (:) :: bxnp,bynp,bznp
      real, allocatable, dimension (:) :: cxnp,cynp,cznp

      end  module consdernper

!----------------------------------------------------------------------
      module consrs
!-----------------------------------------------------------------------

      real ars,brs,ars1,brs1,crs1,ars2

      end module consrs

!----------------------------------------------------------------------
      module combustion
!-----------------------------------------------------------------------

      real :: zf,dyodz,dhcomb
      real, allocatable, dimension (:,:,:) :: hcomb,locflama
      real, allocatable, dimension (:,:,:,:) :: flama

      end module combustion

!-----------------------------------------------------------------------

      module source_calor
!----------------------------------------------------------------------

      real :: tm

      end module source_calor

!-----------------------------------------------------------------------

      module afluid
!----------------------------------------------------------------------

      real :: umax, dxyz_min,cs
      integer :: nc
      
      end module afluid
