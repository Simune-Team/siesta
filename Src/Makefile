
# Makefile for siesta
#
.SUFFIXES: .f .F .o .a  .f90 .F90
#
default: what siesta gen-basis  # denchar
#
include arch.make
#
# Uncomment the following line for debugging support
#
#FFLAGS=$(FFLAGS_DEBUG)
#
what:
	@echo
	@echo "Compilation architecture to be used: ${SIESTA_ARCH}"
	@echo "If this is not what you want, create the right"
	@echo "arch.make file using the models in Sys"
	@echo
	@echo "Hit ^C to abort..."
	@sleep 2
#
SYSOBJ=$(SYS).o
#
# Note that machine-specific files are now in top Src directory.
#
OBJS = arw.o  atomlwf.o atompop.o bands.o bessph.o cgwf.o chkdim.o chkgmx.o \
	chempot.o coceri.o conjgr.o constr.o coxmol.o cross.o \
	denmat.o detover.o dfscf.o dhscf.o diagon.o digcel.o fft3d.o \
	diagg.o diagk.o diagkp.o diag2g.o diag2k.o diagpol.o \
	diagsprl.o dipole.o dismin.o dnaefs.o dot.o dynamics.o \
	efield.o egandd.o ener3.o extrapol.o extrapolon.o \
	fixed.o forhar.o gradient.o grdsam.o \
	hsparse.o idiag.o  electrostatic.o initatom.o initdm.o inver.o \
	iodm.o iohs.o iolwf.o iorho.o ioxv.o ipack.o iopipes.o \
	kgrid.o kgridinit.o kinefsm.o ksv.o ksvinit.o linpack.o \
	madelung.o matel.o meshmatrix.o memory.o meshsubs.o \
        minvec.o mulliken.o naefs.o neighb.o nlefsm.o \
	ordern.o outcell.o outcoor.o overfsm.o \
	paste.o pdos.o pdosg.o pdosk.o phirphi.o pixmol.o plcharge.o \
	timestamp.o propor.o pulayx.o \
	ranger.o ran3.o recipes.o reclat.o redata.o redcel.o \
	reinit.o reord.o rhoofd.o rhoofdsp.o rhooda.o \
	savepsi.o shaper.o timer.o \
	vmb.o vmat.o vmatsp.o volcel.o xc.o xijorb.o \
	cellxc.o cdiag.o rdiag.o \
        cgvc.o iocg.o ioeig.o iofa.o iokp.o iomd.o repol.o typecell.o \
        ofc.o poison.o readsp.o radfft.o \
	siesta.o io.o \
        spin_init.o coor.o transfer.o \
        broadcast_basis.o sig.o eggbox.o dsyevds.o \
        optical.o phirphi_opt.o reoptical.o transition_rate.o
#--------------------------------------------------------------
#
WXML=libwxml.a
$(WXML): 
	(cd wxml ; $(MAKE) module)
#
# Routines using wxml calls.
#
pdos.o : $(WXML)
#--------------------------------------------------------------
FDF=libfdf.a
$(FDF): 
	(cd fdf ; $(MAKE) module)
#
# Routines using fdf calls.
#
dhscf.o initdm.o iodm.o iohs.o iolwf.o iorho.o grdsam.o : $(FDF)
recoor.o outcoor.o ioxv.o kgrid.o kgridinit.o ksv.o : $(FDF)
redata.o siesta.o diagon.o pdos.o plcharge.o: $(FDF)
coor.o cgvc.o iocg.o ioeig.o iofa.o iokp.o iomd.o repol.o ofc.o : $(FDF)
readsp.o eggbox.o iopipes.o : $(FDF)
#
atom.o basis_specs.o: $(FDF)
gen-basis.o: $(FDF)
#
#
# This is crude but will have to do for now.
# Note : precision must be the first module
#
MOD_OBJS=precision.o parallel.o parallelsubs.o atom.o atmparams.o atmfuncs.o listsc.o \
         memoryinfo.o numbvect.o sorting.o \
	 atomlist.o ionew.o atm_types.o old_atmfuncs.o radial.o parsing.o\
         alloc.o phonon.o spher_harm.o periodic_table.o version.o \
         basis_types.o xml.o pseudopotential.o basis_specs.o sys.o basis_io.o\
         chemical.o densematrix.o writewave.o on_subs.o fermid.o
COM_OBJS=$(OBJS) $(SYSOBJ)
ALL_OBJS=$(MOD_OBJS) $(COM_OBJS)
#
initatom.o: basis_specs.o basis_types.o basis_io.o
atom.o: basis_types.o
basis_specs.o: pseudopotential.o basis_types.o chemical.o
basis_types.o: pseudopotential.o atmparams.o
atom.o initatom.o transfer.o  plcharge.o   : old_atmfuncs.o
basis_types.o atom.o old_atmfuncs.o atmfuncs.o: atmparams.o
atomlwf.o dhscf.o dnaefs.o                 : atmfuncs.o
efield.o matel.o mulliken.o outcoor.o : atmfuncs.o
overfsm.o shaper.o                                : atmfuncs.o
atomlist.o                                        : atmfuncs.o 
kinefsm.o naefs.o nlefsm.o overfsm.o              : atmfuncs.o
siesta.o          : atmfuncs.o phonon.o ordern.o hsparse.o iopipes.o
redata.o          : phonon.o
dfscf.o hsparse.o : atmfuncs.o listsc.o
rhoofd.o rhoofdsp.o vmat.o vmatsp.o  : listsc.o
hsparse.o minvec.o : sorting.o
chempot.o         : numbvect.o
memory.o          : memoryinfo.o
dhscf.o dfscf.o forhar.o rhooda.o rhoofd.o rhoofdsp.o cellxc.o : meshsubs.o
vmat.o vmatsp.o : meshsubs.o
dfscf.o rhoofd.o rhoofdsp.o vmat.o vmatsp.o : meshmatrix.o
siesta.o: atomlist.o 
old_atmfuncs.o atomlist.o siesta.o atom.o arw.o pseudopotential.o: ionew.o
iodm.o pdos.o iopipes.o: ionew.o
coor.o: atomlist.o
pdos.o radial.o: xml.o
atm_types.o: radial.o
transfer.o: atm_types.o atmfuncs.o radial.o periodic_table.o
atmfuncs.o: atm_types.o radial.o
old_atmfuncs.o atmfuncs.o matel.o: spher_harm.o
#
diag2g.o diag2k.o diagg.o diagk.o diagkp.o diagsprl.o ionew.o: sys.o
atmfuncs.o atom.o basis_specs.o: sys.o
chkdim.o coor.o diagpol.o matel.o old_atmfuncs.o periodic_table.o: sys.o
pdos.o: atmfuncs.o atomlist.o
propor.o spher_harm.o pseudopotential.o: sys.o
redata.o dhscf.o: sys.o
diag2g.o diag2k.o diagg.o diagk.o diagkp.o diagsprl.o: fermid.o
optical.o transition_rate.o: fermid.o
#----------------------------------------------------
matel.o basis_io.o atom.o : radfft.o
#
basis_io.o: ionew.o atm_types.o atmfuncs.o radial.o sys.o chemical.o\
            basis_types.o xml.o
basis_io.o: $(FDF) $(NETCDF_INTERFACE)
siesta.o bands.o efield.o fixed.o grdsam.o writewave.o: parsing.o
initdm.o repol.o ksv.o: parsing.o
#
siesta.o atomlist.o coor.o dfscf.o matel.o: alloc.o
#
broadcast_basis.o: atm_types.o 
atom.o basis_specs.o: periodic_table.o
atompop.o:  dhscf.o
#
$(COM_OBJS): precision.o parallel.o $(MPI_INTERFACE)
#
# Interfaces to libraries
#
libmpi_f90.a: 
	@(cd MPI ; $(MAKE) )
libnetcdf_f90.a: 
	@(cd NetCDF ; $(MAKE) )
#
# Libraries that might need to be compiled
#
linalg.a: Libs/blas.f Libs/lapack.f
	@echo "==> Compiling linalg.a in Libs..."
	@(cd Libs ; $(MAKE) linalg.a)
dc_lapack.a:  Libs/dc_lapack.f
	@echo "==> Compiling dc_lapack.a in Libs..."
	@(cd Libs ; $(MAKE) dc_lapack.a)
#
version:
	@echo
	@echo "==> Incorporating information about present compilation (compiler and flags)"
	@sed  "s/SIESTA_ARCH/${SIESTA_ARCH}/g" version.F90 > temp0.F90
	@sed  "s/FFLAGS/${FC} ${FFLAGS}/g" temp0.F90 > compinfo.F90
	@($(MAKE) -s "FC=$(FC)" "FFLAGS=$(FFLAGS)" compinfo.o)
	@rm -f temp0.F90 compinfo.F90
	@mv compinfo.o version.o
	@echo

siesta: version $(MPI_INTERFACE) $(NETCDF_INTERFACE) $(FDF) $(WXML)\
                $(COMP_LIBS) $(ALL_OBJS)
	$(FC) -o siesta \
	       $(LDFLAGS) $(ALL_OBJS) $(FDF) $(WXML) $(MPI_INTERFACE)\
               $(NETCDF_INTERFACE) $(COMP_LIBS) $(LIBS) 
#
#------------------------------------------------------------------------
GEN-BASIS_OBJS=gen-basis.o basis_types.o precision.o parallel.o \
          basis_io.o  chemical.o transfer.o atm_types.o\
          atmparams.o old_atmfuncs.o radial.o  io.o \
	  paste.o ionew.o \
          basis_specs.o atom.o periodic_table.o\
          pseudopotential.o dot.o xc.o recipes.o arw.o  \
          sys.o timer.o  xml.o $(SYSOBJ)
# precision must be the first one here...
gen-basis.o: precision.o atom.o basis_types.o basis_specs.o basis_io.o $(FDF)
gen-basis: version $(NETCDF_INTERFACE) $(MPI_INTERFACE) $(FDF) $(WXML) \
                 $(GEN-BASIS_OBJS)
	$(FC) -o gen-basis \
	       $(LDFLAGS) $(GEN-BASIS_OBJS) $(MPI_INTERFACE) $(FDF) $(WXML) \
               $(NETCDF_INTERFACE) $(LIBS) 
#
#------------------------------------------------------------------------
# Denchar is now at the main source level...
#
DENCHAR_OBJS = m_denchar_init.o m_denchar_geom.o m_denchar_io.o \
               m_denchar_neighb.o m_denchar_work.o \
               denchar.o
#
m_denchar_io.o m_denchar_init.o denchar.o : $(FDF)
denchar.o: precision.o parallel.o basis_io.o listsc.o 
denchar.o: m_denchar_init.o m_denchar_io.o m_denchar_work.o m_denchar_geom.o
#
OTHER_DENCHAR_OBJS=  precision.o \
        f2kcli.o bessph.o chkdim.o dismin.o  dot.o iodm.o memory.o \
        paste.o  radfft.o recipes.o io.o \
        volcel.o parallel.o memoryinfo.o sys.o listsc.o  \
        atmparams.o atmfuncs.o ionew.o atm_types.o \
        radial.o spher_harm.o  basis_io.o basis_types.o \
        pseudopotential.o chemical.o xml.o $(SYSOBJ)

denchar: $(OTHER_DENCHAR_OBJS) $(DENCHAR_OBJS) 
	$(FC) -o denchar \
	   $(LDFLAGS) $(DENCHAR_OBJS) $(OTHER_DENCHAR_OBJS) \
           $(MPI_INTERFACE) $(FDF) $(NETCDF_INTERFACE)  $(LIBS) 
#----------------------------------------------------------------------------
#
clean: 
	@echo "==> Cleaning object, library, and executable files"
	rm -f siesta gen-basis denchar *.o  *.a *.pcl *.pc
	rm -f *.mod
	(cd fdf ; $(MAKE) clean)
	(cd wxml ; $(MAKE) clean)
	@if [ -d MPI ] ; then (cd MPI && $(MAKE) clean ) ; fi
	(cd NetCDF ; $(MAKE) clean )
	(cd Libs ; $(MAKE) clean )
#
pristine:  clean
	cp -fp Include/constr.f .
#



