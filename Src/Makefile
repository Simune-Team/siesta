# 
# This file is part of the SIESTA package.
#
# Copyright (c) Fundacion General Universidad Autonoma de Madrid:
# E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
# and J.M.Soler, 1996-2006.
# 
# Use of this software constitutes agreement with the full conditions
# given in the SIESTA license, as signed by all legitimate users.
#
#
# Makefile for siesta
#
# Remove implicit rules (including for Modula-2, which just confuses matters)
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#
default: what siesta # gen-basis denchar fsiesta.o
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
#
aclocal.m4: $(M4_FILES)
	cat $(M4_FILES)  > aclocal.m4
dep: 
	sfmakedepend --depend=obj --modext=o *.f *.f90 *.F *.F90
	@sed '/^siesta_cmlsubs.o:/s/version.o//'  Makefile > tmp.Makefile
	@mv tmp.Makefile Makefile
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
OBJS = automatic_cell.o \
        arw.o  atomlwf.o bands.o bessph.o bonds.o \
        born_charge.o cgwf.o chkdim.o chkgmx.o \
	chempot.o coceri.o constr.o coxmol.o cross.o \
	denmat.o denmatlomem.o detover.o dfscf.o dhscf.o diagon.o digcel.o fft3d.o \
	diagg.o diagk.o diagkp.o diag2g.o diag2k.o diagpol.o \
	diagsprl.o dipole.o dismin.o dnaefs.o dot.o  \
	efield.o egandd.o ener3.o ener3lomem.o extrapol.o extrapolon.o \
	fixed.o flib_spline.o gradient.o gradientlomem.o grdsam.o \
	hsparse.o idiag.o  initatom.o initdm.o inver.o \
	iodm.o iohs.o iolwf.o iozm.o ipack.o iopipes.o \
	kgrid.o kgridinit.o kinefsm.o ksv.o ksvinit.o linpack.o local_DOS.o\
	madelung.o matel.o meshmatrix.o memory.o meshsubs.o metaforce.o \
        m_check_supercell.o mulliken.o minvec.o naefs.o neighb.o \
	new_dm.o normalize_dm.o ordern.o outcell.o outcoor.o \
	paste.o pdos.o pdosg.o pdosk.o pdoskp.o phirphi.o pixmol.o plcharge.o \
	projected_DOS.o propor.o pulayx.o \
	ranger.o ran3.o reclat.o redcel.o \
	reinit.o reord.o rhoofd.o rhoofdsp.o rhooda.o \
	savepsi.o shaper.o timer.o \
	vmb.o vmat.o vmatsp.o volcel.o xc.o xijorb.o \
	cellxc.o cdiag.o rdiag.o \
        cgvc.o cgvc_zmatrix.o  \
        iocg.o ioeig.o iofa.o iokp.o iomd.o typecell.o \
        ofc.o poison.o readsp.o radfft.o \
	write_md_record.o kpoint_grid.o find_kgrid.o proximity_check.o\
	state_init.o siesta_move.o setup_hamiltonian.o compute_dm.o mixer.o\
	scfconvergence_test.o post_scf_work.o state_analysis.o write_subs.o \
	siesta_todo.o siesta_init.o struct_init.o siesta_options.o siesta_geom.o \
        siesta_analysis.o siesta_end.o siesta_forces.o \
        siesta.o io.o sparse_matrices.o \
        spin_init.o coor.o atm_transfer.o stress.o \
        broadcast_basis.o eggbox.o dsyevds.o zheevds.o \
        optical.o phirphi_opt.o reoptical.o transition_rate.o \
        initparallel.o setspatial.o setatomnodes.o uncell.o \
        cart2frac.o obc.o \
        precision.o sys.o m_cell.o recipes.o files.o \
        spatial.o parallel.o parallelsubs.o  parsing.o chemical.o xcmod.o atom.o atmparams.o \
        m_mpi_utils.o m_fdf_global.o m_history.o m_iorho.o  \
        atmfuncs.o listsc.o memoryinfo.o m_memory.o numbvect.o sorting.o \
        atomlist.o atm_types.o old_atmfuncs.o radial.o \
        alloc.o phonon.o spher_harm.o periodic_table.o version.o timestamp.o \
        basis_types.o xml.o pseudopotential.o basis_specs.o basis_io.o\
        onmod.o densematrix.o writewave.o on_subs.o fermid.o m_broyddj.o \
        electrostatic.o mneighb.o globalise.o siesta_cmlsubs.o siesta_cml.o \
        units.o zmatrix.o m_broyden_mixing.o forhar.o m_walltime.o m_wallclock.o m_iostruct.o \
        nlefsm.o overfsm.o overlap.o conjgr.o conjgr_old.o m_energies.o m_steps.o\
        m_broyddj_nocomm.o broyden_optim.o ioxv.o dynamics.o md_out.o molecularmechanics.o
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
dhscf.o initdm.o iodm.o iohs.o iolwf.o  grdsam.o : $(FDF)
recoor.o outcoor.o ioxv.o m_iostruct.o kgrid.o kgridinit.o ksv.o : $(FDF)
redata.o siesta.o diagon.o pdos.o plcharge.o: $(FDF)
coor.o cgvc.o iocg.o ioeig.o iofa.o iokp.o iomd.o  ofc.o : $(FDF)
readsp.o eggbox.o iopipes.o : $(FDF)
#
atom.o basis_io.o basis_specs.o: $(FDF)
gen-basis.o: $(FDF)
#
#
ALL_OBJS=$(OBJS) $(SYSOBJ)
#
basis_io.o md_out.o: $(NETCDF_INTERFACE)
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
	  paste.o  flib_spline.o \
          basis_specs.o atom.o periodic_table.o\
          pseudopotential.o dot.o xc.o  arw.o  \
          sys.o timer.o  xml.o m_walltime.o $(SYSOBJ)
# precision must be the first one here...
gen-basis.o: $(FDF)
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
#
OTHER_DENCHAR_OBJS=  precision.o recipes.o \
        f2kcli.o bessph.o chkdim.o dismin.o  dot.o iodm.o memory.o \
        paste.o  radfft.o io.o  spatial.o flib_spline.o\
        volcel.o parallel.o parallelsubs.o memoryinfo.o sys.o listsc.o  \
        atmparams.o atmfuncs.o atm_types.o m_memory.o \
        radial.o spher_harm.o  basis_io.o basis_types.o \
        pseudopotential.o chemical.o xml.o files.o $(SYSOBJ)

denchar: $(OTHER_DENCHAR_OBJS) $(DENCHAR_OBJS) $(COMP_LIBS)
	$(FC) -o denchar \
	   $(LDFLAGS) $(DENCHAR_OBJS) $(OTHER_DENCHAR_OBJS) \
           $(MPI_INTERFACE) $(FDF) $(NETCDF_INTERFACE) $(COMP_LIBS) $(LIBS) 
#----------------------------------------------------------------------------
MT_OBJS=precision.o m_mpi_utils.o sys.o parallel.o 
mpi_test: $(MPI_INTERFACE) $(MT_OBJS) mpi_test.o
	$(FC) -o mpi_test \
	       $(LDFLAGS) mpi_test.o $(MT_OBJS) $(MPI_INTERFACE) $(LIBS) 
#----------------------------------------------------------------------------
MIXPS_OBJS= f2kcli.o flib_spline.o mixps.o  pseudopotential.o sys.o \
            io.o pxf.o parallel.o periodic_table.o

mixps: $(MPI_INTERFACE)  $(MIXPS_OBJS)
	$(FC) -o mixps \
	       $(LDFLAGS) $(MIXPS_OBJS) $(MPI_INTERFACE) $(LIBS) 
#
#------------------------------------------------------------------------
#
clean: 
	@echo "==> Cleaning object, library, and executable files"
	rm -f siesta gen-basis denchar *.o  *.a *.pcl *.pc
	rm -f mixps mpi_test
	rm -f *.mod
	rm -f *.bck *~
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

# DO NOT DELETE THIS LINE - used by make depend
alloc.o: parallel.o precision.o sys.o
arw.o: parallel.o sys.o
atm_transfer.o: atm_types.o atmparams.o old_atmfuncs.o periodic_table.o
atm_transfer.o: radial.o sys.o
atm_types.o: precision.o radial.o
atmfuncs.o: atm_types.o precision.o radial.o spher_harm.o sys.o
atom.o: atmparams.o basis_specs.o basis_types.o old_atmfuncs.o periodic_table.o
atom.o: precision.o pseudopotential.o recipes.o sys.o xcmod.o
atomlist.o: alloc.o atm_types.o atmfuncs.o parallel.o precision.o siesta_geom.o
atomlwf.o: alloc.o atmfuncs.o mneighb.o onmod.o precision.o spatial.o sys.o
automatic_cell.o: atmfuncs.o parallel.o precision.o units.o
bands.o: alloc.o densematrix.o files.o parallel.o parallelsubs.o parsing.o
bands.o: precision.o sys.o
basis_io.o: atm_types.o atmparams.o basis_types.o chemical.o precision.o
basis_io.o: pseudopotential.o radial.o sys.o xml.o
basis_specs.o: basis_types.o chemical.o periodic_table.o precision.o
basis_specs.o: pseudopotential.o sys.o
basis_types.o: atmparams.o precision.o pseudopotential.o sys.o
bessph.o: parallel.o precision.o sys.o
bonds.o: atmfuncs.o precision.o recipes.o sorting.o units.o
born_charge.o: atomlist.o ksvinit.o parallel.o siesta_geom.o siesta_options.o
born_charge.o: siesta_todo.o sparse_matrices.o
broadcast_basis.o: atm_types.o parallel.o radial.o
broyden_optim.o: m_broyddj_nocomm.o m_fdf_global.o m_memory.o m_mpi_utils.o
broyden_optim.o: parallel.o precision.o units.o
cart2frac.o: sys.o
cdiag.o: alloc.o memoryinfo.o parallel.o precision.o sys.o
cellxc.o: meshsubs.o parallel.o parallelsubs.o precision.o sys.o xcmod.o
cgvc.o: conjgr_old.o m_fdf_global.o m_mpi_utils.o parallel.o precision.o
cgvc.o: units.o
cgvc_zmatrix.o: conjgr.o parallel.o precision.o zmatrix.o
cgwf.o: onmod.o precision.o sys.o
chemical.o: precision.o sys.o
chempot.o: numbvect.o parallel.o parallelsubs.o precision.o sys.o
chkdim.o: sys.o
chkgmx.o: minvec.o parallel.o precision.o sys.o
coceri.o: files.o periodic_table.o precision.o
compute_dm.o: atomlist.o kpoint_grid.o m_energies.o m_steps.o ordern.o
compute_dm.o: siesta_geom.o siesta_options.o siesta_todo.o sparse_matrices.o
compute_dm.o: sys.o
conjgr.o: precision.o
conjgr_old.o: precision.o
coor.o: alloc.o m_mpi_utils.o parallel.o precision.o siesta_geom.o sys.o
coor.o: zmatrix.o
coxmol.o: files.o periodic_table.o precision.o
denchar.o: basis_io.o f2kcli.o listsc.o m_denchar_geom.o m_denchar_init.o
denchar.o: m_denchar_io.o m_denchar_work.o parallel.o precision.o sys.o
denmat.o: globalise.o m_mpi_utils.o onmod.o precision.o
denmatlomem.o: globalise.o onmod.o precision.o
densematrix.o: precision.o
detover.o: parallel.o parallelsubs.o precision.o
dfscf.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o meshmatrix.o
dfscf.o: meshsubs.o meshsubs.o parallel.o precision.o sys.o
dhscf.o: atmfuncs.o forhar.o m_iorho.o meshsubs.o parallel.o parsing.o
dhscf.o: precision.o sys.o xcmod.o
diag2g.o: fermid.o parallel.o parallelsubs.o precision.o sys.o
diag2k.o: fermid.o parallel.o parallelsubs.o precision.o sys.o
diagg.o: fermid.o parallel.o parallelsubs.o precision.o sys.o writewave.o
diagk.o: fermid.o parallel.o parallelsubs.o precision.o sys.o writewave.o
diagkp.o: fermid.o parallel.o parallelsubs.o precision.o sys.o
diagon.o: alloc.o densematrix.o m_memory.o parallel.o parallelsubs.o
diagon.o: precision.o sys.o
diagpol.o: precision.o sys.o
diagsprl.o: atomlist.o fermid.o parallel.o parallelsubs.o precision.o
diagsprl.o: siesta_geom.o sys.o
dipole.o: parallel.o precision.o sys.o
dnaefs.o: atmfuncs.o precision.o
dynamics.o: atomlist.o files.o ioxv.o m_mpi_utils.o parallel.o precision.o
dynamics.o: sys.o units.o
efield.o: atmfuncs.o parallel.o parallelsubs.o parsing.o precision.o
efield.o: siesta_cml.o sys.o
egandd.o: alloc.o globalise.o m_mpi_utils.o onmod.o onmod.o precision.o sys.o
eggbox.o: parallel.o parsing.o precision.o
electrostatic.o: atm_types.o atmfuncs.o atmparams.o precision.o radfft.o
electrostatic.o: radial.o recipes.o sys.o
ener3.o: globalise.o m_mpi_utils.o onmod.o precision.o
ener3lomem.o: globalise.o m_mpi_utils.o onmod.o precision.o
extrapol.o: parallel.o precision.o sys.o
extrapolon.o: parallel.o precision.o sys.o
fermid.o: parallel.o precision.o recipes.o sys.o
fft3d.o: parallel.o parallelsubs.o precision.o sys.o
find_kgrid.o: alloc.o minvec.o parallel.o precision.o units.o
fixed.o: alloc.o parallel.o parsing.o precision.o
flib_spline.o: precision.o
forhar.o: meshsubs.o precision.o
gen-basis.o: atom.o basis_io.o basis_specs.o basis_types.o old_atmfuncs.o
gen-basis.o: parallel.o precision.o sys.o xcmod.o
globalise.o: alloc.o onmod.o onmod.o precision.o spatial.o
gradient.o: globalise.o m_mpi_utils.o onmod.o precision.o
gradientlomem.o: globalise.o onmod.o precision.o
grdsam.o: parallel.o parsing.o precision.o sys.o
hsparse.o: atmfuncs.o listsc.o parallel.o parallelsubs.o precision.o sorting.o
idiag.o: parallel.o sys.o
initatom.o: atom.o basis_io.o basis_specs.o basis_types.o electrostatic.o
initatom.o: old_atmfuncs.o precision.o
initdm.o: parallel.o parallelsubs.o parsing.o precision.o sys.o units.o
initparallel.o: alloc.o parallel.o parallelsubs.o precision.o spatial.o sys.o
io.o: sys.o
iocg.o: files.o precision.o sys.o
iodm.o: files.o parallel.o parallelsubs.o precision.o sys.o
ioeig.o: files.o precision.o siesta_cml.o units.o
iofa.o: files.o precision.o
iohs.o: files.o parallel.o parallelsubs.o precision.o sys.o
iokp.o: files.o precision.o
iolwf.o: alloc.o files.o onmod.o parallel.o parallelsubs.o precision.o sys.o
iomd.o: files.o precision.o
iopipes.o: parallel.o precision.o sys.o
ioxv.o: files.o parallel.o precision.o
iozm.o: files.o parallel.o precision.o siesta_geom.o zmatrix.o
ipack.o: sys.o
kgrid.o: parallel.o precision.o units.o
kgridinit.o: minvec.o parallel.o precision.o units.o
kinefsm.o: atmfuncs.o parallel.o parallelsubs.o precision.o
kpoint_grid.o: find_kgrid.o minvec.o parallel.o precision.o siesta_cml.o
kpoint_grid.o: siesta_options.o units.o
ksv.o: alloc.o atmfuncs.o densematrix.o ksvinit.o parallel.o parsing.o
ksv.o: precision.o sys.o
ksvinit.o: alloc.o parallel.o parsing.o precision.o
local_DOS.o: atomlist.o files.o kpoint_grid.o m_energies.o m_mpi_utils.o
local_DOS.o: parallel.o siesta_geom.o siesta_options.o siesta_todo.o
local_DOS.o: sparse_matrices.o sys.o
m_broyddj.o: alloc.o m_mpi_utils.o parallel.o precision.o sys.o
m_broyddj_nocomm.o: alloc.o parallel.o precision.o sys.o
m_broyden_mixing.o: m_broyddj.o m_memory.o m_mpi_utils.o parallel.o precision.o
m_cell.o: precision.o siesta_geom.o
m_check_supercell.o: minvec.o parallel.o precision.o
m_denchar_geom.o: precision.o sys.o
m_denchar_init.o: files.o
m_denchar_io.o: m_denchar_geom.o
m_denchar_neighb.o: precision.o sys.o
m_denchar_work.o: atmfuncs.o chemical.o listsc.o m_denchar_geom.o
m_denchar_work.o: m_denchar_io.o m_denchar_neighb.o parallel.o
m_energies.o: precision.o
m_fdf_global.o: m_mpi_utils.o parallel.o precision.o
m_iorho.o: parallel.o parallelsubs.o precision.o sys.o
m_iostruct.o: alloc.o files.o m_mpi_utils.o parallel.o precision.o
m_iostruct.o: siesta_geom.o sys.o units.o
m_memory.o: memoryinfo.o parallel.o sys.o
m_mpi_utils.o: precision.o sys.o
m_options.o: precision.o
m_wallclock.o: m_walltime.o
madelung.o: parallel.o precision.o sys.o
matel.o: alloc.o atmfuncs.o precision.o radfft.o recipes.o spher_harm.o
md_out.o: files.o precision.o sys.o units.o
md_utils.o: precision.o
memory.o: memoryinfo.o parallel.o precision.o
memoryinfo.o: precision.o
meshmatrix.o: atomlist.o meshsubs.o parallel.o parallelsubs.o precision.o
meshsubs.o: atm_types.o atmfuncs.o parallel.o parallelsubs.o precision.o
meshsubs.o: radial.o siesta_cml.o
metaforce.o: parallel.o parsing.o precision.o
minvec.o: precision.o sorting.o sys.o
mixer.o: atomlist.o m_broyden_mixing.o m_steps.o siesta_geom.o siesta_options.o
mixer.o: siesta_todo.o sparse_matrices.o
mixps.o: f2kcli.o flib_spline.o periodic_table.o precision.o pseudopotential.o
mixps.o: sys.o
mneighb.o: alloc.o precision.o sys.o
molecularmechanics.o: parallel.o parsing.o precision.o sys.o units.o
mpi_test.o: m_mpi_utils.o parallel.o precision.o
mulliken.o: atmfuncs.o parallel.o parallelsubs.o precision.o siesta_cml.o
naefs.o: atmfuncs.o precision.o
neighb.o: precision.o
new_dm.o: atomlist.o m_steps.o siesta_geom.o siesta_options.o siesta_todo.o
new_dm.o: sparse_matrices.o
nlefsm.o: atmfuncs.o parallel.o parallelsubs.o precision.o
normalize_dm.o: atomlist.o parallel.o siesta_todo.o sparse_matrices.o
ofc.o: files.o precision.o
old_atmfuncs.o: atmparams.o precision.o sys.o
on_subs.o: alloc.o onmod.o onmod.o spatial.o
onmod.o: precision.o
optical.o: alloc.o atmfuncs.o densematrix.o fermid.o files.o parallel.o
optical.o: parallelsubs.o parsing.o precision.o sys.o
ordern.o: alloc.o globalise.o m_mpi_utils.o onmod.o parallel.o precision.o
ordern.o: siesta_cml.o sys.o units.o
outcell.o: precision.o siesta_cml.o units.o
outcoor.o: atmfuncs.o precision.o siesta_geom.o sys.o units.o
overfsm.o: atmfuncs.o m_memory.o parallel.o parallelsubs.o precision.o
overlap.o: alloc.o atmfuncs.o parallel.o parallelsubs.o precision.o
parallelsubs.o: parallel.o spatial.o sys.o
pdos.o: alloc.o atmfuncs.o atomlist.o densematrix.o files.o parallel.o
pdos.o: parallelsubs.o precision.o siesta_geom.o units.o xml.o
pdosg.o: parallel.o parallelsubs.o precision.o sys.o
pdosk.o: parallel.o parallelsubs.o precision.o sys.o
pdoskp.o: parallel.o parallelsubs.o precision.o sys.o
periodic_table.o: precision.o sys.o
phirphi.o: atmfuncs.o parallel.o parallelsubs.o precision.o
phirphi_opt.o: alloc.o atm_types.o atmfuncs.o atmparams.o parallel.o
phirphi_opt.o: parallelsubs.o precision.o spher_harm.o sys.o
phonon.o: files.o parallel.o precision.o sys.o
pixmol.o: files.o periodic_table.o precision.o
plcharge.o: files.o old_atmfuncs.o parallel.o precision.o
poison.o: parallel.o precision.o sys.o
post_scf_work.o: atomlist.o m_energies.o setup_hamiltonian.o siesta_cml.o
post_scf_work.o: siesta_geom.o siesta_options.o siesta_todo.o sparse_matrices.o
projected_DOS.o: atomlist.o kpoint_grid.o m_mpi_utils.o parallel.o
projected_DOS.o: siesta_options.o siesta_todo.o sparse_matrices.o sys.o
propor.o: precision.o sys.o
proximity_check.o: atomlist.o parallel.o precision.o siesta_geom.o
proximity_check.o: siesta_options.o siesta_todo.o
pseudopotential.o: flib_spline.o precision.o sys.o
pulayx.o: alloc.o files.o parallel.o precision.o sys.o
radfft.o: precision.o recipes.o
radial.o: precision.o xml.o
ranger.o: precision.o sys.o
rdiag.o: alloc.o memoryinfo.o parallel.o precision.o sys.o
readsp.o: parallel.o precision.o sys.o
recipes.o: sys.o
redcel.o: precision.o sys.o
reinit.o: files.o parallel.o precision.o
reoptical.o: parallel.o parsing.o precision.o sys.o
reord.o: precision.o
rhooda.o: atmfuncs.o atomlist.o meshsubs.o meshsubs.o precision.o
rhoofd.o: atm_types.o atmfuncs.o atomlist.o listsc.o meshmatrix.o meshsubs.o
rhoofd.o: meshsubs.o parallel.o precision.o sys.o
rhoofdsp.o: atm_types.o atmfuncs.o atomlist.o listsc.o meshmatrix.o meshsubs.o
rhoofdsp.o: meshsubs.o precision.o sys.o
savepsi.o: parallel.o parallelsubs.o precision.o
scfconvergence_test.o: m_energies.o m_wallclock.o parallel.o siesta_cml.o
scfconvergence_test.o: siesta_options.o siesta_todo.o write_subs.o
setatomnodes.o: alloc.o parallel.o precision.o spatial.o
setspatial.o: alloc.o parallel.o precision.o spatial.o
setup_hamiltonian.o: atmfuncs.o atomlist.o m_energies.o m_steps.o metaforce.o
setup_hamiltonian.o: molecularmechanics.o nlefsm.o normalize_dm.o overfsm.o
setup_hamiltonian.o: siesta_cml.o siesta_geom.o siesta_options.o siesta_todo.o
setup_hamiltonian.o: sparse_matrices.o stress.o
setup_kscell.o: minvec.o parallel.o precision.o units.o
shaper.o: atmfuncs.o precision.o
siesta.o: m_steps.o siesta_analysis.o siesta_end.o siesta_forces.o
siesta.o: siesta_init.o siesta_move.o
siesta_analysis.o: atomlist.o bands.o files.o kpoint_grid.o ksvinit.o
siesta_analysis.o: local_DOS.o m_energies.o m_steps.o parallel.o
siesta_analysis.o: projected_DOS.o siesta_cml.o siesta_geom.o siesta_options.o
siesta_analysis.o: siesta_todo.o sparse_matrices.o write_subs.o writewave.o
siesta_analysis.o: zmatrix.o
siesta_cml.o: siesta_cmlsubs.o
siesta_cmlsubs.o: files.o parallel.o timestamp.o 
siesta_end.o: alloc.o m_wallclock.o parallel.o siesta_cml.o siesta_todo.o
siesta_end.o: timestamp.o
siesta_forces.o: compute_dm.o m_steps.o mixer.o parallel.o post_scf_work.o
siesta_forces.o: setup_hamiltonian.o siesta_options.o siesta_todo.o
siesta_forces.o: state_analysis.o state_init.o
siesta_geom.o: precision.o
siesta_init.o: alloc.o atomlist.o bands.o densematrix.o files.o kpoint_grid.o
siesta_init.o: ksvinit.o m_check_supercell.o m_energies.o m_fdf_global.o
siesta_init.o: m_iostruct.o m_mpi_utils.o m_steps.o m_wallclock.o metaforce.o
siesta_init.o: molecularmechanics.o parallel.o parallelsubs.o phonon.o
siesta_init.o: siesta_cmlsubs.o siesta_geom.o siesta_options.o siesta_todo.o
siesta_init.o: sparse_matrices.o sys.o timestamp.o writewave.o xcmod.o
siesta_move.o: atomlist.o broyden_optim.o dynamics.o iopipes.o ioxv.o
siesta_move.o: m_energies.o m_steps.o parallel.o phonon.o siesta_cml.o
siesta_move.o: siesta_geom.o siesta_options.o siesta_todo.o stress.o
siesta_move.o: write_subs.o zmatrix.o
siesta_options.o: files.o m_fdf_global.o m_mpi_utils.o memoryinfo.o parallel.o
siesta_options.o: precision.o siesta_cml.o sys.o units.o
siesta_todo.o: m_mpi_utils.o md_out.o precision.o
sparse_matrices.o: precision.o
spatial.o: precision.o
spher_harm.o: precision.o sys.o
spin_init.o: precision.o
state_analysis.o: atomlist.o born_charge.o m_energies.o m_wallclock.o
state_analysis.o: parallel.o siesta_cml.o siesta_geom.o siesta_options.o
state_analysis.o: siesta_todo.o sparse_matrices.o stress.o write_subs.o
state_analysis.o: zmatrix.o
state_init.o: alloc.o atomlist.o hsparse.o kpoint_grid.o m_check_supercell.o
state_init.o: m_energies.o m_steps.o new_dm.o overlap.o parallel.o siesta_cml.o
state_init.o: siesta_cmlsubs.o siesta_geom.o siesta_options.o siesta_todo.o
state_init.o: sparse_matrices.o write_subs.o zmatrix.o
stress.o: precision.o
struct_init.o: alloc.o atmfuncs.o atomlist.o files.o iopipes.o ioxv.o
struct_init.o: m_fdf_global.o m_iostruct.o parallel.o periodic_table.o
struct_init.o: siesta_cml.o siesta_geom.o siesta_options.o siesta_todo.o
struct_init.o: zmatrix.o
sys.o: parallel.o
timer.o: m_walltime.o parallel.o precision.o
transition_rate.o: fermid.o parallel.o parallelsubs.o precision.o sys.o
typecell.o: precision.o
uncell.o: precision.o
units.o: precision.o
vmat.o: atm_types.o atmfuncs.o atomlist.o listsc.o meshmatrix.o meshsubs.o
vmat.o: meshsubs.o parallel.o precision.o
vmatsp.o: atm_types.o atmfuncs.o atomlist.o listsc.o meshmatrix.o meshsubs.o
vmatsp.o: meshsubs.o precision.o
vmb.o: parallel.o precision.o sys.o
write_md_record.o: atomlist.o m_energies.o m_steps.o md_out.o parallel.o
write_md_record.o: siesta_geom.o siesta_options.o siesta_todo.o
write_subs.o: atomlist.o m_energies.o m_iostruct.o m_steps.o parallel.o
write_subs.o: siesta_cml.o siesta_geom.o siesta_options.o siesta_todo.o
write_subs.o: stress.o units.o
writewave.o: alloc.o atmfuncs.o atomlist.o densematrix.o files.o parallel.o
writewave.o: parallelsubs.o parsing.o precision.o siesta_geom.o sys.o
xc.o: precision.o sys.o xcmod.o
xcmod.o: parallel.o parsing.o precision.o sys.o
xijorb.o: parallel.o parallelsubs.o precision.o
xml.o: precision.o
zmatrix.o: alloc.o m_cell.o parallel.o parsing.o precision.o siesta_geom.o
zmatrix.o: sys.o units.o
band.o: bands.o
diagmemory.o: memoryinfo.o
listsc_module.o: listsc.o
m_born_charge.o: born_charge.o
m_broyden_optim.o: broyden_optim.o
m_compute_dm.o: compute_dm.o
m_conjgr.o: conjgr.o
m_conjgr_old.o: conjgr_old.o
m_dynamics.o: dynamics.o
m_fermid.o: fermid.o
m_find_kgrid.o: find_kgrid.o
m_forhar.o: forhar.o
m_hsparse.o: hsparse.o
m_ioxv.o: ioxv.o
m_ksvinit.o: ksvinit.o
m_local_dos.o: local_DOS.o
m_minvec.o: minvec.o
m_mixer.o: mixer.o
m_new_dm.o: new_dm.o
m_nlefsm.o: nlefsm.o
m_normalize_dm.o: normalize_dm.o
m_ordern.o: ordern.o
m_overfsm.o: overfsm.o
m_overlap.o: overlap.o
m_post_scf_work.o: post_scf_work.o
m_projected_dos.o: projected_DOS.o
m_radfft.o: radfft.o
m_recipes.o: recipes.o
m_scf_options.o: m_options.o
m_setup_hamiltonian.o: setup_hamiltonian.o
m_siesta_analysis.o: siesta_analysis.o
m_siesta_end.o: siesta_end.o
m_siesta_forces.o: siesta_forces.o
m_siesta_init.o: siesta_init.o
m_siesta_move.o: siesta_move.o
m_state_analysis.o: state_analysis.o
m_state_init.o: state_init.o
m_stress.o: stress.o
m_timestamp.o: timestamp.o
m_xyz_phiatm.o: matel.o
mesh.o: meshsubs.o
meshdscf.o: meshmatrix.o
meshphi.o: meshsubs.o
neighbour.o: mneighb.o
on_core.o: onmod.o
on_main.o: onmod.o
version_info.o: version.o
