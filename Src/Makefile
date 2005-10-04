
# Makefile for siesta
#
# Remove implicit rules (including for Modula-2, which just confuses matters)
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#
default: what siesta gen-basis # denchar
#
configure: configure.ac aclocal.m4
	autoconf
M4_FILES = m4/ac_config_aux_dirs.m4  m4/TW_CHECK_FC_FPP.m4      m4/TW_FC_ID_FLAGS.m4  \
           m4/ACX_MPI.m4             m4/TW_CHECK_FC_TR15580.m4  m4/TW_FC_ID.m4        \
           m4/fortran.m4             m4/TW_CHECK_FC_TR15581.m4  m4/TW_FIND_FC_BLAS.m4 \
           m4/TW_CHECK_BLACS.m4      m4/TW_CHECK_SCALAPACK.m4   m4/TW_FIND_LAPACK.m4  \
           m4/TW_CHECK_FC_90.m4      m4/TW_FC_CHECK_ABORT.m4    m4/TW_PATH_NETCDF.m4  \
           m4/TW_CHECK_FC_95.m4      m4/TW_FC_CHECK_DCFUNS.m4   m4/TW_TRY_DC_LAPACK.m4\
           m4/TW_CHECK_FC_FPP_90.m4  m4/TW_FC_CHECK_FLUSH.m4    m4/TW_FC_KINDS.m4
aclocal.m4: $(M4_FILES)
	cat $(M4_FILES)  > aclocal.m4
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
SYSOBJ=$(SYS).o pxf.o
#
# Note that machine-specific files are now in top Src directory.
#
OBJS = arw.o  atomlwf.o bands.o bessph.o cgwf.o chkdim.o chkgmx.o \
	chempot.o coceri.o conjgr.o constr.o coxmol.o cross.o \
	denmat.o detover.o dfscf.o dhscf.o diagon.o digcel.o fft3d.o \
	diagg.o diagk.o diagkp.o diag2g.o diag2k.o diagpol.o \
	diagsprl.o dipole.o dismin.o dnaefs.o dot.o dynamics.o \
	efield.o egandd.o ener3.o ener3lomem.o extrapol.o extrapolon.o \
	fixed.o forhar.o gradient.o grdsam.o \
	hsparse.o idiag.o  initatom.o initdm.o inver.o \
	iodm.o iohs.o iolwf.o iorho.o ioxv.o ipack.o iopipes.o \
	kgrid.o kgridinit.o kinefsm.o ksv.o ksvinit.o linpack.o \
	madelung.o matel.o meshmatrix.o memory.o meshsubs.o \
        minvec.o mulliken.o naefs.o neighb.o nlefsm.o \
	ordern.o outcell.o outcoor.o overfsm.o \
	paste.o pdos.o pdosg.o pdosk.o pdoskp.o phirphi.o pixmol.o plcharge.o \
	propor.o pulayx.o \
	ranger.o ran3.o reclat.o redata.o redcel.o \
	reinit.o reord.o rhoofd.o rhoofdsp.o rhooda.o \
	savepsi.o shaper.o timer.o \
	vmb.o vmat.o vmatsp.o volcel.o xc.o xijorb.o \
	cellxc.o cdiag.o rdiag.o \
        cgvc.o iocg.o ioeig.o iofa.o iokp.o iomd.o repol.o typecell.o \
        ofc.o poison.o readsp.o radfft.o \
	siesta.o io.o \
        spin_init.o coor.o atm_transfer.o \
        broadcast_basis.o sig.o eggbox.o dsyevds.o zheevds.o \
        optical.o phirphi_opt.o reoptical.o transition_rate.o \
        initparallel.o setspatial.o setatomnodes.o uncell.o \
        cart2frac.o 
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
XMLPARSER=libxmlparser.a
$(XMLPARSER): 
	(cd xmlparser ; $(MAKE) module)
#
# Routines using xmlparser calls.
#
######  somefile.o : $(WXML)
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
MOD_OBJS=precision.o recipes.o \
         spatial.o parallel.o parallelsubs.o sys.o parsing.o xcmod.o atom.o atmparams.o \
         m_mpi_utils.o m_wrappers.o m_history.o \
         atmfuncs.o listsc.o memoryinfo.o m_memory.o numbvect.o sorting.o \
	 atomlist.o ionew.o atm_types.o old_atmfuncs.o radial.o m_smearing.o \
         alloc.o phonon.o spher_harm.o periodic_table.o version.o timestamp.o \
         basis_types.o xml.o pseudopotential.o basis_specs.o basis_io.o\
         onmod.o chemical.o densematrix.o writewave.o on_subs.o fermid.o m_broyddj.o \
	 electrostatic.o mneighb.o globalise.o siesta_cmlsubs.o siesta_cml.o units.o m_broyden_mixing.o m_wallclock.o
#
COM_OBJS=$(OBJS) $(SYSOBJ)
ALL_OBJS=$(MOD_OBJS) $(COM_OBJS)
#
initatom.o: basis_specs.o basis_types.o basis_io.o
atom.o: basis_types.o
basis_specs.o: pseudopotential.o basis_types.o chemical.o
basis_types.o: pseudopotential.o atmparams.o
atom.o initatom.o atm_transfer.o  plcharge.o   : old_atmfuncs.o
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
m_memory.o          : memoryinfo.o
dhscf.o dfscf.o forhar.o rhooda.o rhoofd.o rhoofdsp.o cellxc.o : meshsubs.o
vmat.o vmatsp.o : meshsubs.o
dfscf.o rhoofd.o rhoofdsp.o vmat.o vmatsp.o : meshmatrix.o
siesta.o: atomlist.o 
old_atmfuncs.o atomlist.o siesta.o atom.o arw.o pseudopotential.o: ionew.o
iodm.o pdos.o iopipes.o: ionew.o
coor.o: atomlist.o
pdos.o radial.o: xml.o
atm_types.o: radial.o
atm_transfer.o: atm_types.o atmfuncs.o radial.o periodic_table.o
atmfuncs.o: atm_types.o radial.o
old_atmfuncs.o atmfuncs.o matel.o: spher_harm.o
siesta.o parallel.o initparallel.o setspatial.o setatomnodes.o : spatial.o
ener3.o ener3lomem.o denmat.o gradient.o globalise : spatial.o
ordern.o atomlwf.o cgwf.o egandd.o ener3.o ener3lomem.o gradient.o denmat.o on_subs.o iolwf.o : onmod.o
atomlwf.o : mneighb.o
gradient.o ener3.o ener3lomem.o denmat.o : globalise.o
#
diag2g.o diag2k.o diagg.o diagk.o diagkp.o diagsprl.o ionew.o: sys.o
atmfuncs.o atom.o basis_specs.o: sys.o
chkdim.o coor.o diagpol.o matel.o old_atmfuncs.o periodic_table.o: sys.o
pdos.o: atmfuncs.o atomlist.o
propor.o spher_harm.o pseudopotential.o: sys.o
redata.o dhscf.o: sys.o
diag2g.o diag2k.o diagg.o diagk.o diagkp.o diagsprl.o: fermid.o
optical.o transition_rate.o: fermid.o
atom.o xc.o cellxc.o dhscf.o gen-basis.o siesta.o: xcmod.o

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
#
###$(COM_OBJS): precision.o parallel.o $(MPI_INTERFACE)
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
libblas.a:
	@echo "==> Compiling libblas.a in Libs..."
	@(cd Libs ; $(MAKE) libblas.a)
liblapack.a:
	@echo "==> Compiling liblapack.a in Libs..."
	@(cd Libs ; $(MAKE) liblapack.a)
linalg.a: Libs/blas.f Libs/lapack.f
	@echo "==> Compiling linalg.a in Libs..."
	@(cd Libs ; $(MAKE) linalg.a)
dc_lapack.a:  Libs/dc_lapack.f
	@echo "==> Compiling dc_lapack.a in Libs..."
	@(cd Libs ; $(MAKE) dc_lapack.a)
#

version: version.F90
	@echo
	@echo "==> Incorporating information about present compilation (compiler and flags)"
	@sed "s'SIESTA_ARCH'${SIESTA_ARCH}'g;s'FFLAGS'${FC} ${FFLAGS}'g;s'SIESTA_VERSION'$$(cat ../version.info)'g" version.F90 > compinfo.F90
	@($(MAKE) compinfo.o)
	@mv compinfo.o version.o
	#@rm -f compinfo.F90
	@echo

siesta: version $(MPI_INTERFACE) $(NETCDF_INTERFACE) $(FDF) $(WXML) $(XMLPARSER) \
                $(COMP_LIBS) $(ALL_OBJS) 
	$(FC) -o siesta \
	       $(LDFLAGS) $(ALL_OBJS) $(FDF) $(WXML) $(XMLPARSER) $(MPI_INTERFACE)\
               $(NETCDF_INTERFACE) $(COMP_LIBS) $(LIBS) 
#
#------------------------------------------------------------------------
GEN-BASIS_OBJS=gen-basis.o basis_types.o precision.o recipes.o parallel.o \
          parsing.o xcmod.o basis_io.o  chemical.o atm_transfer.o atm_types.o\
          atmparams.o old_atmfuncs.o radial.o io.o \
	  paste.o ionew.o \
          basis_specs.o atom.o periodic_table.o\
          pseudopotential.o dot.o xc.o  arw.o  \
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
OTHER_DENCHAR_OBJS=  precision.o recipes.o \
        f2kcli.o bessph.o chkdim.o dismin.o  dot.o iodm.o memory.o \
        paste.o  radfft.o recipes.o io.o  spatial.o\
        volcel.o parallel.o parallelsubs.o memoryinfo.o sys.o listsc.o  \
        atmparams.o atmfuncs.o ionew.o atm_types.o m_memory.o \
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
	rm -f aux_*.f aux_*.f90
	(cd fdf ; $(MAKE) clean)
	(cd wxml ; $(MAKE) clean)
	(cd xmlparser ; $(MAKE) clean)
	@if [ -d MPI ] ; then (cd MPI && $(MAKE) clean ) ; fi
	(cd NetCDF ; $(MAKE) clean )
	(cd Libs ; $(MAKE) clean )
#
pristine:  clean
	cp -fp Include/constr.f .
#



