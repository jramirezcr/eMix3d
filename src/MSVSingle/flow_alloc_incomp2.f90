!-----------------------------------------------------------------------
      subroutine flow_alloc()
!-----------------------------------------------------------------------

      use dimensiones
      use jacobtools
      use mallagrid
      use derivvel
      use mflujos
      use dmflujos
      use viscosidades
      use velocidades
      use variables
      use right
      use tiempo
      use sgdmodel
      use vorticidad
      use consderper_f
      use consdernper_f
      use stat
      use fronterapl
      use bob
      use bob_yz
      use derivtools
      use consderper
      use consdernper
      use flotacion
      use combustion
      use mezclador
      use maskvar
      use impellermotion
      IMPLICIT NONE
    


      allocate(axp(nx))
      allocate(bxp(nx))
      allocate(cxp(nx))
      allocate(ayp(ny))
      allocate(byp(ny))
      allocate(cyp(ny))
      allocate(azp(nz))
      allocate(bzp(nz))
      allocate(czp(nz))
      allocate(axnp(nx))
      allocate(bxnp(nx))
      allocate(cxnp(nx))
      allocate(aynp(ny))
      allocate(bynp(ny))
      allocate(cynp(ny))
      allocate(aznp(nz))
      allocate(bznp(nz))
      allocate(cznp(nz))

      allocate(du(nx))
      allocate(dv(ny))
      allocate(dw(nz))
      allocate(jbnp(nx,ny,nz,11))
      allocate(jbnm(nx,ny,nz,11))
      allocate(jbn(nx,ny,nz,11))
      allocate(x(nx,ny,nz,nd))
      allocate(ddvelp(nx,ny,nz,9))
      allocate(ddvelm(nx,ny,nz,9))
      allocate(dcvel(nx,ny,nz,9))
      allocate(dvel(nx,ny,nz,9))
      allocate(dconcp(nx,ny,nz,nd))
      allocate(dconcm(nx,ny,nz,nd))
      allocate(dconc(nx,ny,nz,nd))
      allocate(dtempp(nx,ny,nz,nd))
      allocate(dtempm(nx,ny,nz,nd))
      allocate(dtemp(nx,ny,nz,nd))
      allocate(e(nx,ny,nz))
      allocate(f(nx,ny,nz))
      allocate(g(nx,ny,nz))
      allocate(ex(nx,ny,nz))
      allocate(fx(nx,ny,nz))
      allocate(gx(nx,ny,nz))
      allocate(evx(nx,ny,nz))
      allocate(fvx(nx,ny,nz))
      allocate(gvx(nx,ny,nz))
      allocate(evp(nx,ny,nz))
      allocate(ev(nx,ny,nz))
      allocate(fv(nx,ny,nz))
      allocate(gv(nx,ny,nz))
      allocate(fvp(nx,ny,nz))
      allocate(gvp(nx,ny,nz))
      allocate(evm(nx,ny,nz))
      allocate(fvm(nx,ny,nz))
      allocate(gvm(nx,ny,nz))
      allocate(dere(nx,ny,nz))
      allocate(derf(nx,ny,nz))
      allocate(derg(nx,ny,nz))
      allocate(derev(nx,ny,nz))
      allocate(derfv(nx,ny,nz))
      allocate(dergv(nx,ny,nz))
      allocate(vis(nx,ny,nz))
      allocate(vis2(nx,ny,nz))
      allocate(vis6(nx,ny,nz))
      allocate(vis5(nx,ny,nz))
      allocate(u(nx,ny,nz,nd))
      allocate(conc(nx,ny,nz))
      allocate(temp(nx,ny,nz))
      allocate(pres(nx,ny,nz))
      allocate(pres0(nx,ny,nz))
      allocate(um(nx,ny,nz,nd))
      allocate(um0(nx,ny,nz,nd))
      allocate(um1(nx,ny,nz,nd))
      allocate(rs(nx,ny,nz))
      allocate(rsv(nx,ny,nz))
      allocate(dx(nx,ny,nz))
      allocate(dy(nx,ny,nz))
      allocate(dz(nx,ny,nz))
      allocate(amut(nx,ny,nz))
      allocate(dxmsgd(nx,ny,nz))
      allocate(dxpsgd(nx,ny,nz))
      allocate(dymsgd(nx,ny,nz))
      allocate(dypsgd(nx,ny,nz))
      allocate(dzmsgd(nx,ny,nz))
      allocate(dzpsgd(nx,ny,nz))
      allocate(dxyzsgd(nx,ny,nz))
      allocate(wx(nx,ny,nz))
      allocate(wy(nx,ny,nz))
      allocate(wz(nx,ny,nz))
      allocate(st(nx,ny,nz,17))
      allocate(frontx(2,nd,ny,nz))
      allocate(frontz(2,nd,nx,ny))
      allocate(fronty(2,nd,nx,nz))
      allocate(frontin(ny,nz,nd))
      allocate(imask(2,nx,ny))
      allocate(nent(ny,nz))
      allocate(esplayer(ny,nz,nd))
      allocate(esource(nx))
      allocate(esplayer_yz(nx,ny,nd))
      allocate(esplayer_zy(nx,nz,nd))
      allocate(esource_yz(ny))
      allocate(esource_zy(nz))
      allocate(axpf(nx))
      allocate(bxpf(nx))
      allocate(cxpf(nx))
      allocate(aypf(ny))
      allocate(bypf(ny))
      allocate(cypf(ny))
      allocate(azpf(nz))
      allocate(bzpf(nz))
      allocate(czpf(nz))
      allocate(axnpf(nx))
      allocate(bxnpf(nx))
      allocate(cxnpf(nx))
      allocate(aynpf(ny))
      allocate(bynpf(ny))
      allocate(cynpf(ny))
      allocate(aznpf(nz))
      allocate(bznpf(nz))
      allocate(cznpf(nz))
      allocate(cm(nx))

      allocate(flama(nx,ny,nz,3))
      allocate(locflama(nx,ny,nz))
      allocate(dconczp(nx,ny,nz,3))
      allocate(dconczm(nx,ny,nz,3))
      allocate(concz(nx,ny,nz))
      allocate(vis5z(nx,ny,nz))
      allocate(hcomb(nx,ny,nz))
      allocate(mask(nx,ny,nz))
      allocate(impellermask(nx,ny,nz))

      allocate(xx(nx,ny))
      allocate(yy(nx,ny))
      allocate(rmin(nx,ny))
      allocate(rsub(nx,ny))
!      allocate(erro(nx,ny))
      allocate(rtot(nx,ny))
      allocate(rtoty(nx,ny))
!      allocate(erroy(nx,ny))
      allocate(rtotx(nx,ny))
      allocate(triangleGeo(5000))


      
      return
      end subroutine flow_alloc
!-----------------------------------------------------------------------
      subroutine flow_dealloc()
!-----------------------------------------------------------------------
      use dimensiones
      use jacobtools
      use mallagrid
      use derivvel
      use mflujos
      use dmflujos
      use viscosidades
      use velocidades
      use variables
      use right
      use tiempo
      use sgdmodel
      use consderper_f
      use derivtools
      use consdernper_f
      use vorticidad
      use stat
      use fronterapl
      use bob
      use bob_yz
      use consderper
      use consdernper
      use flotacion
      use combustion
      use mezclador
      use maskvar
      use impellermotion
      IMPLICIT NONE

      deallocate(axp)
      deallocate(bxp)
      deallocate(cxp)
      deallocate(ayp)
      deallocate(byp)
      deallocate(cyp)
      deallocate(azp)
      deallocate(bzp)
      deallocate(czp)
      deallocate(axnp)
      deallocate(bxnp)
      deallocate(cxnp)
      deallocate(aynp)
      deallocate(bynp)
      deallocate(cynp)
      deallocate(aznp)
      deallocate(bznp)
      deallocate(cznp)

      deallocate(du)
      deallocate(dv)
      deallocate(dw)
      deallocate(jbnp)
      deallocate(jbnm)
      deallocate(jbn)
      deallocate(x)
      deallocate(ddvelp)
      deallocate(ddvelm)
      deallocate(dcvel)
      deallocate(dconcp)
      deallocate(dconcm)
      deallocate(dtempp)
      deallocate(dtempm)
      deallocate(e)
      deallocate(f)
      deallocate(g)
      deallocate(ex)
      deallocate(fx)
      deallocate(gx)
      deallocate(evx)
      deallocate(fvx)
      deallocate(gvx)
      deallocate(dere)
      deallocate(derf)
      deallocate(derg)
      deallocate(derev)
      deallocate(derfv)
      deallocate(dergv)
      deallocate(evp)
      deallocate(ev)
      deallocate(fv)
      deallocate(gv)
      deallocate(fvp)
      deallocate(gvp)
      deallocate(evm)
      deallocate(fvm)
      deallocate(gvm)
      deallocate(vis)
      deallocate(vis2)
      deallocate(vis6)
      deallocate(vis5)
      deallocate(u)
      deallocate(conc)
      deallocate(temp)
      deallocate(pres)
      deallocate(pres0)
      deallocate(um)
      deallocate(um0)
      deallocate(um1)
      deallocate(rs)
      deallocate(rsv)
      deallocate(dx)
      deallocate(dy)
      deallocate(dz)
      deallocate(amut)
      deallocate(dxmsgd)
      deallocate(dxpsgd)
      deallocate(dymsgd)
      deallocate(dypsgd)
      deallocate(dzmsgd)
      deallocate(dzpsgd)
      deallocate(dxyzsgd)
      deallocate(wx)
      deallocate(wy)
      deallocate(wz)
      deallocate(st)
      deallocate(frontx)
      deallocate(frontz)
      deallocate(fronty)
      deallocate(frontin)
      deallocate(imask)
      deallocate(nent)
      deallocate(esplayer)
      deallocate(esource)
      deallocate(esplayer_yz)
      deallocate(esplayer_zy)
      deallocate(esource_yz)
      deallocate(esource_zy)
      deallocate(axpf)
      deallocate(bxpf)
      deallocate(cxpf)
      deallocate(aypf)
      deallocate(bypf)
      deallocate(cypf)
      deallocate(azpf)
      deallocate(bzpf)
      deallocate(czpf)
      deallocate(axnpf)
      deallocate(bxnpf)
      deallocate(cxnpf)
      deallocate(aynpf)
      deallocate(bynpf)
      deallocate(cynpf)
      deallocate(aznpf)
      deallocate(bznpf)
      deallocate(cznpf)
      deallocate(cm)

      deallocate(flama)
      deallocate(locflama)
      deallocate(dconczp)
      deallocate(dconczm)
      deallocate(concz)
      deallocate(vis5z)
      deallocate(hcomb)
      deallocate(mask)
      deallocate(impellermask)

      deallocate(rmin)
      deallocate(rsub)
!      deallocate(erro)
      deallocate(rtot)
      deallocate(rtoty)
!      deallocate(erroy)
      deallocate(rtotx)
      deallocate(triangleGeo)
      return
      end subroutine flow_dealloc
