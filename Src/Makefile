# Makefile for siesta
#
.SUFFIXES: .f .F .o .a  .f90 .F90
#
default: what siesta gen-basis
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
OBJS = arw.o  atomlwf.o bands.o bessph.o cgwf.o chkdim.o chkgmx.o \
	chempot.o coceri.o conjgr.o constr.o coxmol.o cross.o \
	denmat.o detover.o dfscf.o dhscf.o diagon.o digcel.o fft3d.o \
	diagg.o diagk.o diagkp.o diag2g.o diag2k.o diagpol.o \
	dipole.o dismin.o dnaefs.o dot.o dynamics.o \
	efield.o egandd.o ener3.o extrapol.o extrapolon.o \
	fermid.o fermispin.o fixed.o forhar.o gradient.o grdsam.o \
	hsparse.o idiag.o  initatom.o initdm.o inver.o  \
	iodm.o iohs.o iolwf.o iorho.o ioxv.o  \
	kgrid.o kgridinit.o kinefsm.o ksv.o ksvinit.o linpack.o \
	madelung.o matel.o meshmatrix.o memory.o meshsubs.o \
        minvec.o mulliken.o naefs.o neighb.o nlefsm.o \
	on_subs.o ordern.o outcell.o outcoor.o overfsm.o \
	paste.o pdos.o pdosk.o phirphi.o pixmol.o timestamp.o \
	propor.o pulayx.o  \
	ranger.o ran3.o recipes.o reclat.o redata.o redcel.o\
	reinit.o reord.o rhoofd.o rhooda.o  \
	savepsi.o shaper.o timer.o \
	vmb.o vmat.o volcel.o xc.o xijorb.o  \
	cellxc.o cdiag.o rdiag.o \
        cgvc.o iocg.o ioeig.o iofa.o iokp.o iomd.o repol.o typecell.o \
        ofc.o poison.o radfft.o  \
	siesta.o version.o io.o\
        spin_init.o coor.o transfer.o\
        broadcast_basis.o sig.o
#
FDF=libfdf.a
$(FDF): 
	(cd fdf ; $(MAKE) "FC=$(FC)" "FFLAGS=$(FFLAGS)" module)
#
# Routines using fdf calls.
#
dhscf.o initdm.o iodm.o iohs.o iolwf.o iorho.o grdsam.o : $(FDF)
recoor.o outcoor.o ioxv.o kgrid.o kgridinit.o ksv.o : $(FDF)
redata.o siesta.o diagon.o pdos.o: $(FDF)
coor.o cgvc.o iocg.o ioeig.o iofa.o iokp.o iomd.o repol.o ofc.o : $(FDF)
#
atom.o basis_specs.o: $(FDF)
gen-basis.o: $(FDF)
#
#
# This is crude but will have to do for now.
# Note : precision must be the first module
#
MOD_OBJS=precision.o atom.o atmparams.o atmfuncs.o listsc.o \
         memoryinfo.o numbvect.o  parallel.o sorting.o \
	 atomlist.o ionew.o atm_types.o old_atmfuncs.o radial.o parsing.o\
         interfaces.o  alloc.o phonon.o spher_harm.o periodic_table.o\
         basis_types.o pseudopotential.o basis_specs.o sys.o basis_io.o\
         chemical.o xml.o
COM_OBJS=$(OBJS) $(SYSOBJ)
ALL_OBJS=$(MOD_OBJS) $(COM_OBJS)
#
initatom.o: basis_specs.o basis_types.o basis_io.o
atom.o: basis_types.o
basis_specs.o: pseudopotential.o basis_types.o chemical.o
basis_types.o: pseudopotential.o atmparams.o
atom.o initatom.o transfer.o     : old_atmfuncs.o
basis_types.o atom.o old_atmfuncs.o atmfuncs.o: atmparams.o
atomlwf.o dhscf.o dnaefs.o                 : atmfuncs.o
efield.o matel.o mulliken.o outcoor.o : atmfuncs.o
overfsm.o shaper.o                                : atmfuncs.o
atomlist.o                                        : atmfuncs.o 
kinefsm.o naefs.o nlefsm.o overfsm.o              : atmfuncs.o
siesta.o          : atmfuncs.o phonon.o
redata.o          : phonon.o
dfscf.o hsparse.o : atmfuncs.o listsc.o
rhoofd.o vmat.o   : listsc.o
hsparse.o minvec.o : sorting.o
chempot.o         : numbvect.o
memory.o          : memoryinfo.o
dhscf.o dfscf.o forhar.o rhooda.o rhoofd.o vmat.o : meshsubs.o
dfscf.o rhoofd.o vmat.o : meshmatrix.o
siesta.o: atomlist.o 
old_atmfuncs.o atomlist.o siesta.o atom.o: ionew.o
coor.o: atomlist.o
radial.o: xml.o
atm_types.o: radial.o
transfer.o: atm_types.o atmfuncs.o radial.o periodic_table.o
atmfuncs.o: atm_types.o radial.o
old_atmfuncs.o atmfuncs.o matel.o: spher_harm.o
#
diag2g.o diag2k.o diagg.o diagk.o diagkp.o ionew.o: sys.o
atmfuncs.o atom.o basis_specs.o: sys.o
chkdim.o coor.o diagpol.o matel.o old_atmfuncs.o periodic_table.o: sys.o
pdos.o: xml.o atmfuncs.o atomlist.o
propor.o spher_harm.o pseudopotential.o: sys.o
redata.o dhscf.o: sys.o
#
#
basis_io.o: ionew.o atm_types.o atmfuncs.o radial.o sys.o chemical.o\
            basis_types.o xml.o
basis_io.o: $(FDF) $(NETCDF_INTERFACE)
siesta.o bands.o efield.o fixed.o grdsam.o: parsing.o
initdm.o repol.o ksv.o: parsing.o
#
siesta.o: interfaces.o
siesta.o atomlist.o coor.o dfscf.o matel.o: alloc.o
#
broadcast_basis.o: atm_types.o 
atom.o basis_specs.o: periodic_table.o
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
# Libraries that need to be compiled
#
linalg.a: Libs/blas.f Libs/lapack.f
	@(cd Libs ; $(MAKE) "FC=$(FC)" "FFLAGS=$(FFLAGS)")
#
version:
	sed  "s/SIESTA_ARCH/${SIESTA_ARCH}/g" version.F > temp0.F
	sed  "s/FFLAGS/${FC} ${FFLAGS}/g" temp0.F > temp.F
	$(FC) -c $(FFLAGS) $(DEFS) temp.F
	@rm -f temp0.F temp.F
	@mv temp.o version.o
	@echo

siesta: version $(MPI_INTERFACE) $(NETCDF_INTERFACE) $(FDF) \
                $(COMP_LIBS) $(ALL_OBJS)
	$(FC) -o siesta \
	       $(LDFLAGS) $(ALL_OBJS) $(FDF) $(MPI_INTERFACE)\
               $(NETCDF_INTERFACE) $(COMP_LIBS) $(LIBS) 
#
GEN-BASIS_OBJS=gen-basis.o basis_types.o precision.o\
          basis_io.o  chemical.o transfer.o atm_types.o\
          atmparams.o atmfuncs.o old_atmfuncs.o radial.o spher_harm.o io.o \
	  paste.o ionew.o \
          chkdim.o basis_specs.o atom.o periodic_table.o\
          pseudopotential.o dot.o xc.o recipes.o arw.o radfft.o\
          bessph.o sys.o timer.o  xml.o $(SYSOBJ)
# precision must be the first one here...
gen-basis.o: precision.o atom.o basis_types.o basis_specs.o basis_io.o $(FDF)
gen-basis: version $(NETCDF_INTERFACE) $(MPI_INTERFACE) $(FDF) \
                $(COMP_LIBS) $(GEN-BASIS_OBJS)
	$(FC) -o gen-basis \
	       $(LDFLAGS) $(GEN-BASIS_OBJS) $(MPI_INTERFACE) $(FDF) \
               $(NETCDF_INTERFACE) $(COMP_LIBS) $(LIBS) 
#
clean: 
	@echo "==> Cleaning object, library, and executable files"
	rm -f siesta gen-basis *.o  *.a
	rm -f *.mod
	(cd fdf ; $(MAKE) clean)
	@if [ -d MPI ] ; then (cd MPI && $(MAKE) clean ) ; fi
	(cd NetCDF ; $(MAKE) clean )
	(cd Libs ; $(MAKE) clean )
pristine:  clean
	cp -fp Include/constr.f .
#



