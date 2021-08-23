EXENAME  = ~/Procesamiento/ART01/RUSHTON/eMix.run
FC       =  ifort 
FFLAGS   =  -O3  -r8 -w

SRCPATH  =  ../src
COREPATH =  $(SRCPATH)/MSVSingle

INCLUDE = ../include/\
 -I$(SRCPATH)

#CUDA: COMPILER AND OPTIONS

CUDAC     =  nvcc
##LIBS      = -L/usr/local/cuda-8.0/lib64 -lcuda -lcudart
LIBS      = -lcuda -lcudart
##INCLUDECU = -I /usr/local/cuda-8.0/include 
CUOPTION  = -O3 --gpu-architecture=compute_75

OBJS=\
flow_mod_incomp2.o\
multiphase_mod.o\
aleatorio.o\
deriv_vel_cuda.o\
derivada.o\
divergencia_cuda.o\
atl_incomp.o\
esponja.o\
estadisticas_incomp.o\
flow_alloc_incomp2.o\
flujos_compact.o\
frontera_z.o\
frontera_x.o\
frontera_y.o\
gravar_incomp.o\
herramientas_incomp2.o\
iniconst_incomp.o\
jacob_malla.o\
jacob_malla_comp_incomp.o\
leer_incomp.o\
rungek_ico.o\
sgdm_incomp.o\
source_term.o\
source_incomp.o\
cvel_planar.o\
tstep_incomp.o\
viscosidad_incomp.o\
filtro_incomp.o\
rush_ruth_me.o\
rotating.o\
advectLS.o\
lsTransportEc.o\
lsTools.o\
extrapol.o\
interface.o\
cuDeriv.o\
derivCompact.o\
cuDivergence.o

$(EXENAME):  $(OBJS)
	$(FC)  $(OBJS) $(FFLAGS) -I$(INCLUDE) $(LIBS) -o $(EXENAME) -lstdc++
flow_mod_incomp2.o: $(COREPATH)/flow_mod_incomp2.f90
	$(FC) $(FFLAGS) -c  $<
multiphase_mod.o: multiphase_mod.f90
	$(FC) $(FFLAGS) -c  $<
aleatorio.o: $(COREPATH)/aleatorio.f90
	$(FC) $(FFLAGS) -c $<
deriv_vel_cuda.o:deriv_vel_cuda.f90
	$(FC) $(FFLAGS) -c $<
derivada.o:derivada.f90
	$(FC) $(FFLAGS) -c $<
divergencia_cuda.o:divergencia_cuda.f90
	$(FC) $(FFLAGS) -c $<
atl_incomp.o: $(COREPATH)/atl_incomp.f90
	$(FC) $(FFLAGS) -c $<
esponja.o: $(COREPATH)/esponja.f90
	$(FC) $(FFLAGS) -c $<
estadisticas_incomp.o:$(COREPATH)/estadisticas_incomp.f90
	$(FC) $(FFLAGS) -c $<
flow_alloc_incomp2.o: $(COREPATH)/flow_alloc_incomp2.f90
	$(FC) $(FFLAGS) -c $<
flujos_compact.o: flujos_compact.f90
	$(FC) $(FFLAGS) -c $<
frontera_z.o: $(COREPATH)/frontera_z.f90
	$(FC) $(FFLAGS) -c $<
frontera_x.o: $(COREPATH)/frontera_x.f90
	$(FC) $(FFLAGS) -c $<
frontera_y.o: $(COREPATH)/frontera_y.f90
	$(FC) $(FFLAGS) -c $<
gravar_incomp.o: gravar_incomp.f90
	$(FC) $(FFLAGS) -c $<
herramientas_incomp2.o: $(COREPATH)/herramientas_incomp2.f90
	$(FC) $(FFLAGS) -c $<
iniconst_incomp.o: $(COREPATH)/iniconst_incomp.f90
	$(FC) $(FFLAGS) -c $<
jacob_malla.o: jacob_malla.f90
	$(FC) $(FFLAGS) -c $<
jacob_malla_comp_incomp.o: $(COREPATH)/jacob_malla_comp_incomp.f90
	$(FC) $(FFLAGS) -c $<
leer_incomp.o: leer_incomp.f90
	$(FC) $(FFLAGS) -c $<
rungek_ico.o: rungek_ico.f90
	$(FC) $(FFLAGS) -c $<
sgdm_incomp.o: sgdm_incomp.f90
	$(FC) $(FFLAGS) -c $<
source_term.o: source_term.f90
	$(FC) $(FFLAGS) -c $<
source_incomp.o: $(COREPATH)/source_incomp.f90
	$(FC) $(FFLAGS) -c $<
cvel_planar.o: $(COREPATH)/cvel_planar.f90
	$(FC) $(FFLAGS) -c $<
tstep_incomp.o: $(COREPATH)/tstep_incomp.f90
	$(FC) $(FFLAGS) -c $<
viscosidad_incomp.o: $(COREPATH)/viscosidad_incomp.f90
	$(FC) $(FFLAGS) -c $<
filtro_incomp.o: $(COREPATH)/filtro_incomp.f90
	$(FC) $(FFLAGS) -c $<
rush_ruth_me.o: $(SRCPATH)/Geometry/rush_ruth_me.f90
	$(FC) $(FFLAGS) -c $<
rotating.o: $(SRCPATH)/Geometry/rotating.f90
	$(FC) $(FFLAGS) -c $<
advectLS.o: advectLS.f90
	$(FC) -I$(INCLUDE) $(FFLAGS) -c $<
lsTransportEc.o: $(SRCPATH)/LevelSet/lsTransportEc.cu
	$(CUDAC)  $(CUOPTION) -I$(INCLUDE) $(INCLUDECU) -c $<
lsTools.o: $(SRCPATH)/LevelSet/lsTools.cu
	$(CUDAC)  $(CUOPTION) -I$(INCLUDE) $(INCLUDECU) -c $<
extrapol.o: $(SRCPATH)/LevelSet/extrapol.cu
	$(CUDAC)  $(CUOPTION) -I$(INCLUDE) $(INCLUDECU) -c $<
interface.o: $(SRCPATH)/LevelSet/interface.cu
	$(CUDAC)  $(CUOPTION) -I$(INCLUDE) $(INCLUDECU) -c $<
cuDeriv.o:cuDeriv.cu $(SRCPATH)/Engine/cuArray.cu
	$(CUDAC)  $(CUOPTION) -I$(INCLUDE) $(INCLUDECU) -c $<
derivCompact.o:$(SRCPATH)/Schemes/derivCompact.cu $(SRCPATH)/Engine/cuArray.cu
	$(CUDAC)  $(CUOPTION) -I$(INCLUDE) $(INCLUDECU) -c $<
cuDivergence.o: cuDivergence.cu $(SRCPATH)/Engine/cuArray.cu
	$(CUDAC)  $(CUOPTION) -I$(INCLUDE) $(INCLUDECU) -c $<
CLEAN :
	/bin/rm *.o *.mod
clean :
	/bin/rm *.o *.mod
