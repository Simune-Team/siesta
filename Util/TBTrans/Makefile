#
# Makefile for tbtrans
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#
# Please see note about Linear-Algebra libraries below
#
default: what  tbtrans 
#
ABS_TOPSRC:=$(shell pwd)/../../Src
# path as seen from subdirectories
TOPSRC:=../../../Src
VPATH=$(ABS_TOPSRC)
#
OBJDIR=Obj
ARCH_MAKE=../../$(OBJDIR)/arch.make
include $(ARCH_MAKE)

# If compiling without system lapack, you might need something like this
# 
#TBTRANS_LINALG=linalg.a # Or maybe only: zgeev_lapack.a zgesv_lapack.a
LINALG:= $(TBTRANS_LINALG) $(LINALG)

SIESTA_SYS=${SIESTA_ARCH}

#
# Uncomment the following line for debugging support
#
#FFLAGS=$(FFLAGS_DEBUG) -g
#
what:
	@echo
	@echo "Compilation architecture to be used:  $(SIESTA_SYS)"
	@echo "If that is not what you want, give the correct"
	@echo "value to the variable SIESTA_SYS in your shell environment."
	@echo $(COMMENTS)
	@echo

SYSOBJ=$(SYS).o
# Tbtrans:
MOD_TOBJS=  m_tbt_options.o m_tbt_gf.o m_tbt_kpts.o hamil2.o tsread2.o

SIESTA_OBJS= precision.o sys.o m_io.o alloc.o  memoryinfo.o memory.o pxf.o \
             debugmpi.o \
             moreParallelSubs.o m_walltime.o m_timer.o timer.o m_mpi_utils.o \
             m_fdf_global.o paste.o \
             parallel.o  io.o files.o  m_walltime.o zheevds.o

TBTOBJS= $(SIESTA_OBJS) $(MOD_TOBJS) $(SYSOBJ) $(EXTRA) \
         tbtrans.o transmission.o   \
         calc_green.o cdiag_local.o getsfe.o \
         mkRealContour.o reinit_tb.o setupkham2.o TSiohs.o atompdos.o coop.o \
         csolve.o teigchana.o helpers.o \
         Lowdin.o netcdfwrite.o 

#
# Interfaces to libraries
#
FDF=libfdf.a
FDF_MAKEFILE=$(TOPSRC)/fdf/makefile
FDF_INCFLAGS:=-I $(TOPSRC)/fdf $(INCFLAGS)
$(FDF): 
	(cd fdf ; $(MAKE) -f $(FDF_MAKEFILE) "FC=$(FC)" "VPATH=$(TOPSRC)/fdf" \
                          "ARCH_MAKE=../$(ARCH_MAKE)" \
                          "INCFLAGS=$(FDF_INCFLAGS)" "FFLAGS=$(FFLAGS)" module)
#
INCFLAGS:=$(NETCDF_INCFLAGS) $(INCFLAGS)
#
MPI_MAKEFILE=$(TOPSRC)/MPI/Makefile
$(MPI_INTERFACE): 
	(cd MPI ; $(MAKE) -f $(MPI_MAKEFILE) "FC=$(FC)" "VPATH=$(TOPSRC)/MPI" \
                          "MAKEFILES=$(MPI_MAKEFILE)" "ARCH_MAKE=../$(ARCH_MAKE)" \
                          "FFLAGS=$(FFLAGS)" module_built)
#
LINALG_MAKEFILE=$(TOPSRC)/Libs/makefile
LINALG_INCFLAGS:=-I $(TOPSRC)/Libs $(INCFLAGS)
$(LINALG): 
	(cd Libs ; $(MAKE) -f $(LINALG_MAKEFILE) "FC=$(FC)" "VPATH=$(TOPSRC)/Libs" \
                          "ARCH_MAKE=../$(ARCH_MAKE)" \
                          "INCFLAGS=$(LINALG_INCFLAGS)" "FFLAGS=$(FFLAGS)" $(LINALG))
#
tbtrans: $(MPI_INTERFACE) $(TBTOBJS) $(FDF) $(SYSOBJ) $(NETCDF_INTERFACE) $(LINALG)
	$(FC) $(FFLAGS) $(LDFLAGS) -o tbtrans \
              $(TBTOBJS) $(MPI_INTERFACE) $(NETCDF_INTERFACE) \
              $(FDF) $(LIBS) $(LINALG)  || echo "***Please see top of Makefile"
#	@(cd ../../bin ; rm -f tsScan ; ln -s tbtrans tsScan)

#
clean: 
	@echo "==> Cleaning object, library, and executable files"
	rm -f tbtrans *.o *.a *.mod
	(cd fdf ; $(MAKE) -f $(FDF_MAKEFILE) "ARCH_MAKE=../$(ARCH_MAKE)" clean)
	(cd MPI ; $(MAKE) -f $(MPI_MAKEFILE) "ARCH_MAKE=../$(ARCH_MAKE)" clean)
	(cd Libs ; $(MAKE) -f $(MPI_MAKEFILE) "ARCH_MAKE=../$(ARCH_MAKE)" clean)

readTSHScdf: readTSHScdf.F $(NETCDF_INTERFACE)
	$(FC) -o readTSHScdf readTSHScdf.F $(NETCDF_INTERFACE) $(NETCDF_LIBS)

#
# Dependencies
# 
dep:
	sfmakedepend --depend=obj  --modext=o \
          $(ABS_TOPSRC)/*.f $(ABS_TOPSRC)/*.f90 \
          $(ABS_TOPSRC)/*.F $(ABS_TOPSRC)/*.F90 
#
$(TBTOBJS)      : $(FDF)
reinit.o        : $(MPI_INTERFACE)
sys.o           : $(MPI_INTERFACE) parallel.o
m_mpi_utils.o   : sys.o
m_tbt_gf.o      : m_tbt_kpts.o
m_tbt_kpts      : m_tbt_options.o $(MPI_INTERFACE)
TSprversion.o   : TSversion.h
setupkham2.o	: tsread2.o
tbtrans.o	: hamil2.o 
green4.o        : greensub.o
m_tbt_gf.o       : tsread2.o
netcdfwrite.o   : $(NETCDF_INTERFACE)
OutputRegionData.o: Lowdin.F $(NETCDF_INTERFACE)
# DO NOT DELETE THIS LINE - used by make depend
alloc.o: parallel.o precision.o m_io.o sys.o debugmpi.o
arw.o: alloc.o parallel.o precision.o sys.o
atm_transfer.o: atm_types.o atmparams.o old_atmfuncs.o periodic_table.o
atm_transfer.o: radial.o sys.o
atm_types.o: precision.o radial.o
atmfuncs.o: atm_types.o precision.o radial.o spher_harm.o sys.o
atom.o: atmparams.o atom_options.o basis_specs.o basis_types.o m_filter.o
atom.o: old_atmfuncs.o periodic_table.o precision.o pseudopotential.o recipes.o
atom.o: sys.o xcmod.o
atom_options.o: sys.o
atomlist.o: alloc.o atm_types.o atmfuncs.o parallel.o precision.o siesta_geom.o
atomlwf.o: alloc.o atmfuncs.o mneighb.o onmod.o precision.o spatial.o sys.o
automatic_cell.o: atmfuncs.o parallel.o precision.o units.o
bands.o: alloc.o densematrix.o files.o parallel.o parallelsubs.o precision.o
bands.o: sys.o
basis_enthalpy.o: atmfuncs.o atomlist.o precision.o siesta_geom.o units.o
basis_io.o: atm_types.o atmparams.o atom_options.o basis_types.o chemical.o
basis_io.o: precision.o pseudopotential.o radial.o sys.o xml.o
basis_specs.o: basis_types.o chemical.o periodic_table.o precision.o
basis_specs.o: pseudopotential.o sys.o
basis_types.o: alloc.o atmparams.o precision.o pseudopotential.o sys.o
bessph.o: parallel.o precision.o sys.o
bonds.o: alloc.o atmfuncs.o mneighb.o precision.o recipes.o sorting.o units.o
born_charge.o: atomlist.o ksv.o ksvinit.o m_mpi_utils.o m_spin.o parallel.o
born_charge.o: siesta_geom.o siesta_options.o sparse_matrices.o
broadcast_basis.o: atm_types.o parallel.o radial.o
broyden_optim.o: m_broyddj_nocomm.o m_memory.o parallel.o precision.o
broyden_optim.o: siesta_options.o sys.o units.o
cart2frac.o: sys.o
cdiag.o: alloc.o m_spin.o memoryinfo.o parallel.o precision.o sys.o
cell_broyden_optim.o: alloc.o m_broyddj_nocomm.o parallel.o precision.o sys.o
cell_broyden_optim.o: units.o zmatrix.o
cell_fire_optim.o: alloc.o m_fdf_global.o m_fire.o m_memory.o parallel.o
cell_fire_optim.o: precision.o siesta_options.o sys.o zmatrix.o
cellxc.o: alloc.o cellxc_mod.o meshsubs.o moremeshsubs.o parallel.o
cellxc.o: parallelsubs.o precision.o sys.o xcmod.o
cellxc_mod.o: sys.o xcmod.o
cg.o: alloc.o matrixop.o precision.o precond.o sys.o
cgvc.o: alloc.o conjgr_old.o m_mpi_utils.o parallel.o precision.o sys.o units.o
cgvc_zmatrix.o: alloc.o conjgr.o m_mpi_utils.o parallel.o precision.o sys.o
cgvc_zmatrix.o: units.o zmatrix.o
cgwf.o: onmod.o precision.o sys.o
chemical.o: parallel.o precision.o sys.o
chempot.o: alloc.o mneighb.o parallel.o parallelsubs.o precision.o sys.o
chkdim.o: sys.o
chkgmx.o: minvec.o parallel.o precision.o sys.o
coceri.o: files.o periodic_table.o precision.o
compute_dm.o: atomlist.o iodmhs_netcdf.o jacobi_davidson.o kpoint_grid.o
compute_dm.o: m_energies.o m_eo.o m_gamma.o m_mpi_utils.o m_rmaxh.o m_spin.o
compute_dm.o: m_steps.o m_ts_in_siesta.o m_ts_options.o ordern.o parallel.o
compute_dm.o: precision.o siesta_geom.o siesta_options.o sparse_matrices.o
compute_dm.o: sys.o
compute_norm.o: m_mpi_utils.o precision.o sparse_matrices.o
conjgr.o: precision.o
conjgr_old.o: precision.o
coor.o: alloc.o parallel.o precision.o siesta_geom.o sys.o zmatrix.o
coxmol.o: files.o periodic_table.o precision.o
debugmpi.o: parallel.o
denmat.o: alloc.o globalise.o m_mpi_utils.o onmod.o precision.o
denmatlomem.o: alloc.o globalise.o onmod.o precision.o
densematrix.o: alloc.o precision.o
detover.o: alloc.o parallel.o parallelsubs.o precision.o
dfscf.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o meshmatrix.o
dfscf.o: meshsubs.o meshsubs.o parallel.o parallelsubs.o precision.o sys.o
dhscf.o: alloc.o atmfuncs.o cellxc_mod.o doping.o files.o forhar.o
dhscf.o: iogrid_netcdf.o m_dhscf.o m_efield.o m_iorho.o m_ts_in_siesta.o
dhscf.o: m_ts_options.o meshsubs.o moremeshsubs.o parallel.o parsing.o
dhscf.o: precision.o siesta_options.o sys.o xcmod.o
diag2g.o: fermid.o m_spin.o parallel.o parallelsubs.o precision.o sys.o
diag2k.o: fermid.o m_spin.o parallel.o parallelsubs.o precision.o sys.o
diagg.o: alloc.o fermid.o parallel.o parallelsubs.o precision.o sys.o
diagg.o: writewave.o
diagk.o: compute_norm.o fermid.o parallel.o parallelsubs.o precision.o sys.o
diagk.o: writewave.o
diagk_file.o: fermid.o iowfs_netcdf.o parallel.o parallelsubs.o precision.o
diagk_file.o: sys.o writewave.o
diagkp.o: alloc.o fermid.o parallel.o parallelsubs.o precision.o sys.o
diagon.o: alloc.o densematrix.o m_diagon.o m_memory.o parallel.o parallelsubs.o
diagon.o: precision.o siesta_options.o sys.o
diagpol.o: precision.o sys.o
diagsprl.o: atomlist.o fermid.o parallel.o parallelsubs.o precision.o
diagsprl.o: siesta_geom.o sys.o
dipole.o: meshsubs.o parallel.o precision.o sys.o
distriphionmesh.o: alloc.o atmfuncs.o cellxc_mod.o meshsubs.o moremeshsubs.o
distriphionmesh.o: parallel.o precision.o
dnaefs.o: atmfuncs.o mneighb.o precision.o
domain_decom.o: alloc.o parallel.o precision.o printmatrix.o schecomm.o
domain_decom.o: sparse_matrices.o sys.o
doping.o: alloc.o meshsubs.o meshsubs.o parallel.o precision.o sys.o
dscfcomm.o: alloc.o parallel.o parallelsubs.o precision.o schecomm.o
dynamics.o: alloc.o atomlist.o files.o ioxv.o m_mpi_utils.o parallel.o
dynamics.o: precision.o sys.o units.o
egandd.o: alloc.o globalise.o onmod.o onmod.o precision.o sys.o
eggbox.o: parallel.o precision.o
electrostatic.o: atm_types.o atmfuncs.o atmparams.o precision.o radfft.o
electrostatic.o: radial.o recipes.o sys.o
ener3.o: alloc.o globalise.o m_mpi_utils.o onmod.o precision.o
ener3lomem.o: alloc.o globalise.o m_mpi_utils.o onmod.o precision.o
extrapolon.o: parallel.o precision.o sys.o
fermid.o: parallel.o precision.o recipes.o sys.o
fft3d.o: alloc.o parallel.o parallelsubs.o precision.o sys.o
final_H_f_stress.o: alloc.o atmfuncs.o atomlist.o dnaefs.o files.o kinefsm.o
final_H_f_stress.o: m_dipol.o m_energies.o m_forces.o m_gamma.o m_mpi_utils.o
final_H_f_stress.o: m_ntm.o m_spin.o m_steps.o m_stress.o m_ts_io.o
final_H_f_stress.o: m_ts_kpoints.o m_ts_options.o metaforce.o
final_H_f_stress.o: molecularmechanics.o naefs.o nlefsm.o normalize_dm.o
final_H_f_stress.o: overfsm.o siesta_geom.o siesta_options.o sparse_matrices.o
final_H_f_stress.o: sys.o
find_kgrid.o: alloc.o minvec.o parallel.o precision.o units.o
fire_optim.o: alloc.o m_fdf_global.o m_fire.o m_mpi_utils.o parallel.o
fire_optim.o: precision.o siesta_options.o units.o
fixed.o: alloc.o parallel.o precision.o
flib_spline.o: precision.o
forhar.o: alloc.o meshsubs.o moremeshsubs.o precision.o
get_target_stress.o: parallel.o precision.o sys.o units.o
globalise.o: alloc.o onmod.o onmod.o precision.o spatial.o
gmres.o: alloc.o matrixop.o precision.o precond.o sys.o
gradient.o: alloc.o globalise.o m_mpi_utils.o onmod.o precision.o
gradientlomem.o: alloc.o globalise.o onmod.o precision.o
grdsam.o: alloc.o m_dhscf.o parallel.o precision.o sys.o
hsparse.o: alloc.o atm_types.o atmfuncs.o atomlist.o domain_decom.o listsc.o
hsparse.o: mneighb.o parallel.o parallelsubs.o precision.o sorting.o
hsparse.o: sparse_matrices.o sys.o
idiag.o: parallel.o sys.o
initatom.o: atmparams.o atom.o atom_options.o basis_io.o basis_specs.o
initatom.o: basis_types.o electrostatic.o old_atmfuncs.o precision.o
initmesh.o: alloc.o atmfuncs.o meshsubs.o moremeshsubs.o parallel.o
initmesh.o: parallelsubs.o precision.o siesta_cml.o
initparallel.o: alloc.o domain_decom.o parallel.o parallelsubs.o precision.o
initparallel.o: siesta_options.o spatial.o sys.o
io.o: sys.o
iocg.o: files.o precision.o sys.o
iodm.o: alloc.o files.o parallel.o parallelsubs.o precision.o sys.o
iodm_netcdf.o: alloc.o parallel.o parallelsubs.o precision.o sys.o
iodmhs_netcdf.o: alloc.o parallel.o parallelsubs.o precision.o sys.o
ioeig.o: files.o precision.o siesta_cml.o units.o
iofa.o: files.o precision.o
iogrid_netcdf.o: alloc.o parallel.o precision.o sys.o
iohs.o: alloc.o atm_types.o atmfuncs.o atomlist.o files.o parallel.o
iohs.o: parallelsubs.o precision.o siesta_geom.o sys.o
iokp.o: files.o precision.o
iolwf.o: alloc.o files.o onmod.o parallel.o parallelsubs.o precision.o sys.o
iomd.o: files.o precision.o
iopipes.o: parallel.o precision.o sys.o
iowfs_netcdf.o: alloc.o parallel.o parallelsubs.o precision.o sys.o
ioxv.o: files.o parallel.o precision.o
iozm.o: files.o parallel.o precision.o siesta_geom.o zmatrix.o
ipack.o: sys.o
jacobi_davidson.o: alloc.o cg.o domain_decom.o fermid.o gmres.o matrixop.o
jacobi_davidson.o: parallel.o precision.o precond.o printmatrix.o sys.o
kgrid.o: parallel.o precision.o units.o
kgridinit.o: minvec.o parallel.o precision.o sys.o units.o
kinefsm.o: alloc.o atmfuncs.o mneighb.o parallel.o parallelsubs.o precision.o
kpoint_grid.o: find_kgrid.o minvec.o parallel.o precision.o siesta_cml.o
kpoint_grid.o: siesta_options.o sys.o units.o
kpoint_pdos.o: find_kgrid.o minvec.o parallel.o precision.o siesta_options.o
kpoint_pdos.o: sys.o units.o
ksv.o: alloc.o atmfuncs.o densematrix.o ksvinit.o parallel.o precision.o sys.o
ksvinit.o: alloc.o parallel.o precision.o sys.o
listsc.o: alloc.o
local_DOS.o: atomlist.o files.o kpoint_grid.o m_dhscf.o m_energies.o m_eo.o
local_DOS.o: m_forces.o m_gamma.o m_ntm.o m_spin.o parallel.o siesta_geom.o
local_DOS.o: siesta_options.o sparse_matrices.o sys.o
m_broyddj.o: alloc.o m_mpi_utils.o parallel.o precision.o sys.o
m_broyddj_nocomm.o: alloc.o m_fdf_global.o parallel.o precision.o sys.o
m_broyden_mixing.o: alloc.o m_broyddj.o m_mpi_utils.o parallel.o precision.o
m_cell.o: precision.o siesta_geom.o units.o
m_check_supercell.o: minvec.o parallel.o precision.o
m_dhscf.o: precision.o
m_diagon.o: precision.o
m_dipol.o: precision.o
m_efield.o: atmfuncs.o meshsubs.o parallel.o precision.o siesta_cml.o
m_efield.o: siesta_geom.o sys.o units.o
m_energies.o: precision.o
m_eo.o: precision.o
m_fdf_global.o: m_mpi_utils.o parallel.o precision.o
m_fire.o: parallel.o precision.o
m_fire_mixing.o: alloc.o m_fdf_global.o m_fire_para.o m_mpi_utils.o parallel.o
m_fire_mixing.o: precision.o
m_fire_para.o: m_mpi_utils.o parallel.o precision.o
m_forces.o: precision.o
m_getopts.o: f2kcli.o
m_history.o: alloc.o
m_hsx.o: atm_types.o atmfuncs.o atomlist.o files.o parallel.o parallelsubs.o
m_hsx.o: precision.o siesta_geom.o sys.o
m_iodm.o: alloc.o files.o m_mpi_utils.o parallel.o parallelsubs.o precision.o
m_iodm.o: sys.o
m_iorho.o: alloc.o parallel.o parallelsubs.o precision.o sys.o
m_iostruct.o: alloc.o files.o m_mpi_utils.o parallel.o precision.o
m_iostruct.o: siesta_geom.o sys.o units.o
m_kinetic.o: precision.o
m_memory.o: memoryinfo.o parallel.o sys.o
m_mpi_utils.o: precision.o sys.o
m_options.o: precision.o
m_pulay.o: alloc.o atomlist.o files.o m_mpi_utils.o m_spin.o parallel.o
m_pulay.o: precision.o siesta_options.o sparse_matrices.o sys.o
m_rmaxh.o: precision.o
m_sparse.o: alloc.o precision.o
m_spin.o: alloc.o parallel.o precision.o
m_stress.o: precision.o
m_ts_aux_rout.o: sys.o
m_ts_contour.o: files.o m_ts_aux_rout.o m_ts_options.o precision.o sys.o
m_ts_electrode.o: files.o m_ts_aux_rout.o m_ts_io.o m_ts_kpoints.o
m_ts_electrode.o: m_ts_options.o
m_ts_in_siesta.o: find_kgrid.o m_mpi_utils.o m_ts_contour.o m_ts_electrode.o
m_ts_in_siesta.o: m_ts_io.o m_ts_kpoints.o m_ts_options.o m_ts_scattering.o
m_ts_in_siesta.o: meshsubs.o parallel.o parallelsubs.o precision.o sys.o
m_ts_io.o: files.o m_spin.o m_ts_kpoints.o parallel.o parallelsubs.o
m_ts_io.o: precision.o sys.o
m_ts_kpoints.o: minvec.o parallel.o precision.o siesta_cml.o siesta_options.o
m_ts_kpoints.o: sys.o units.o
m_ts_options.o: m_fdf_global.o parallel.o precision.o siesta_options.o sys.o
m_ts_scattering.o: m_ts_aux_rout.o m_ts_options.o parallelsubs.o precision.o
m_ts_scattering.o: sys.o
m_wallclock.o: m_walltime.o
madelung.o: parallel.o precision.o sys.o
matel.o: alloc.o atmfuncs.o parallel.o precision.o radfft.o recipes.o
matel.o: spher_harm.o
matrixop.o: alloc.o domain_decom.o m_mpi_utils.o parallel.o precision.o
md_out.o: files.o m_energies.o precision.o sys.o units.o
md_utils.o: precision.o
memory.o: memoryinfo.o parallel.o precision.o
memoryinfo.o: precision.o
meshcomm.o: meshsubs.o precision.o
meshmatrix.o: alloc.o atomlist.o dscfcomm.o meshsubs.o parallel.o
meshmatrix.o: parallelsubs.o precision.o
meshsubs.o: alloc.o atm_types.o atmfuncs.o precision.o radial.o
metaforce.o: alloc.o parallel.o precision.o sys.o
minvec.o: precision.o sorting.o sys.o
mixer.o: atomlist.o iodm_netcdf.o m_broyden_mixing.o m_energies.o
mixer.o: m_fire_mixing.o m_iodm.o m_mpi_utils.o m_pulay.o m_spin.o m_steps.o
mixer.o: m_ts_io.o m_ts_options.o parallel.o siesta_geom.o siesta_options.o
mixer.o: sparse_matrices.o
mneighb.o: alloc.o precision.o sys.o
molecularmechanics.o: alloc.o parallel.o precision.o sys.o units.o
moremeshsubs.o: alloc.o m_dhscf.o meshsubs.o parallel.o precision.o schecomm.o
moremeshsubs.o: sys.o
mulliken.o: alloc.o atmfuncs.o parallel.o parallelsubs.o precision.o
mulliken.o: siesta_cml.o
naefs.o: atmfuncs.o mneighb.o precision.o
new_dm.o: alloc.o atomlist.o compute_norm.o m_energies.o m_iodm.o m_sparse.o
new_dm.o: m_spin.o m_steps.o m_ts_io.o m_ts_options.o parallel.o parallelsubs.o
new_dm.o: parsing.o precision.o siesta_geom.o siesta_options.o
new_dm.o: sparse_matrices.o sys.o units.o
nlefsm.o: alloc.o atmfuncs.o mneighb.o parallel.o parallelsubs.o precision.o
normalize_dm.o: atomlist.o m_mpi_utils.o m_spin.o parallel.o precision.o
normalize_dm.o: sparse_matrices.o
obc.o: alloc.o precision.o
ofc.o: alloc.o files.o precision.o units.o
old_atmfuncs.o: alloc.o atmparams.o precision.o sys.o
on_subs.o: alloc.o onmod.o onmod.o spatial.o
onmod.o: precision.o
optical.o: alloc.o atmfuncs.o densematrix.o fermid.o files.o parallel.o
optical.o: parallelsubs.o precision.o sys.o units.o
ordern.o: alloc.o chempot.o globalise.o m_mpi_utils.o onmod.o parallel.o
ordern.o: precision.o siesta_cml.o sys.o units.o
outcell.o: precision.o siesta_cml.o units.o
outcoor.o: alloc.o atmfuncs.o precision.o siesta_geom.o sys.o units.o
overfsm.o: alloc.o atmfuncs.o mneighb.o parallel.o parallelsubs.o precision.o
overlap.o: alloc.o atmfuncs.o mneighb.o parallel.o parallelsubs.o precision.o
parallelsubs.o: domain_decom.o parallel.o spatial.o sys.o
pdos.o: alloc.o atmfuncs.o atomlist.o densematrix.o files.o parallel.o
pdos.o: parallelsubs.o precision.o siesta_geom.o units.o xml.o
pdosg.o: alloc.o parallel.o parallelsubs.o precision.o sys.o units.o
pdosk.o: alloc.o parallel.o parallelsubs.o precision.o sys.o units.o
pdoskp.o: alloc.o parallel.o parallelsubs.o precision.o sys.o units.o
periodic_table.o: precision.o sys.o
phionmesh.o: alloc.o atmfuncs.o meshsubs.o meshsubs.o parallel.o precision.o
phirphi.o: alloc.o atmfuncs.o mneighb.o parallel.o parallelsubs.o precision.o
phirphi_opt.o: alloc.o atm_types.o atmfuncs.o atmparams.o mneighb.o parallel.o
phirphi_opt.o: parallelsubs.o precision.o spher_harm.o sys.o
phonon.o: alloc.o files.o parallel.o precision.o sys.o
pixmol.o: files.o periodic_table.o precision.o
plcharge.o: files.o old_atmfuncs.o parallel.o precision.o
poison.o: alloc.o parallel.o precision.o sys.o
post_scf_work.o: atomlist.o final_H_f_stress.o m_energies.o m_spin.o mneighb.o
post_scf_work.o: siesta_cml.o siesta_geom.o siesta_options.o sparse_matrices.o
precond.o: alloc.o matrixop.o precision.o
printmatrix.o: alloc.o
projected_DOS.o: alloc.o atomlist.o kpoint_grid.o kpoint_pdos.o m_eo.o
projected_DOS.o: m_gamma.o m_spin.o parallel.o precision.o siesta_geom.o
projected_DOS.o: siesta_options.o sparse_matrices.o sys.o
propor.o: precision.o sys.o
proximity_check.o: m_ts_options.o mneighb.o parallel.o precision.o
proximity_check.o: siesta_geom.o siesta_options.o
pseudopotential.o: alloc.o atom_options.o flib_spline.o precision.o sys.o
radfft.o: alloc.o precision.o recipes.o
radial.o: precision.o xml.o
rdiag.o: alloc.o m_diagon.o memoryinfo.o parallel.o precision.o sys.o
readsp.o: parallel.o precision.o sys.o
recipes.o: sys.o
redcel.o: parallel.o precision.o sys.o
reinit.o: files.o parallel.o
remove_intramol_pressure.o: precision.o zmatrix.o
reoptical.o: parallel.o precision.o sys.o
reord.o: alloc.o precision.o
rhooda.o: atmfuncs.o atomlist.o meshsubs.o meshsubs.o precision.o
rhoofd.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o meshmatrix.o
rhoofd.o: meshsubs.o meshsubs.o parallel.o parallelsubs.o precision.o sys.o
rhoofdsp.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o meshmatrix.o
rhoofdsp.o: meshsubs.o meshsubs.o precision.o sys.o
savepsi.o: alloc.o parallel.o parallelsubs.o precision.o
scfconvergence_test.o: m_convergence.o m_energies.o m_ts_options.o
scfconvergence_test.o: m_wallclock.o parallel.o siesta_cml.o siesta_options.o
scfconvergence_test.o: write_subs.o
schecomm.o: alloc.o
setatomnodes.o: alloc.o parallel.o precision.o spatial.o sys.o
setspatial.o: alloc.o parallel.o precision.o spatial.o
setup_H0.o: alloc.o atmfuncs.o atomlist.o dnaefs.o kinefsm.o m_dhscf.o
setup_H0.o: m_dipol.o m_energies.o m_forces.o m_gamma.o m_mpi_utils.o m_ntm.o
setup_H0.o: m_spin.o m_steps.o m_stress.o m_ts_options.o metaforce.o
setup_H0.o: molecularmechanics.o naefs.o nlefsm.o normalize_dm.o parallel.o
setup_H0.o: siesta_geom.o siesta_options.o sparse_matrices.o sys.o
setup_hamiltonian.o: alloc.o atmfuncs.o atomlist.o m_dhscf.o m_dipol.o
setup_hamiltonian.o: m_energies.o m_gamma.o m_hsx.o m_mpi_utils.o m_ntm.o
setup_hamiltonian.o: m_spin.o m_steps.o m_stress.o m_ts_options.o metaforce.o
setup_hamiltonian.o: molecularmechanics.o normalize_dm.o parallel.o
setup_hamiltonian.o: siesta_cml.o siesta_geom.o siesta_options.o
setup_hamiltonian.o: sparse_matrices.o sys.o
setup_kscell.o: minvec.o parallel.o precision.o units.o
shaper.o: atmfuncs.o mneighb.o precision.o
show_distribution.o: atomlist.o parallel.o parallelsubs.o siesta_geom.o sys.o
siesta.o: alloc.o m_steps.o siesta_analysis.o siesta_end.o siesta_forces.o
siesta.o: siesta_init.o siesta_move.o
siesta_analysis.o: alloc.o atomlist.o bands.o basis_enthalpy.o files.o
siesta_analysis.o: kpoint_grid.o ksv.o ksvinit.o local_DOS.o m_dhscf.o
siesta_analysis.o: m_dipol.o m_energies.o m_eo.o m_forces.o m_gamma.o
siesta_analysis.o: m_mpi_utils.o m_ntm.o m_spin.o m_steps.o parallel.o
siesta_analysis.o: projected_DOS.o siesta_cml.o siesta_geom.o siesta_options.o
siesta_analysis.o: sparse_matrices.o write_subs.o writewave.o zmatrix.o
siesta_cml.o: siesta_cmlsubs.o
siesta_cmlsubs.o: files.o parallel.o timestamp.o version.o
siesta_end.o: alloc.o bands.o densematrix.o domain_decom.o fixed.o m_pulay.o
siesta_end.o: m_wallclock.o meshsubs.o moremeshsubs.o parallel.o siesta_cml.o
siesta_end.o: sparse_matrices.o timestamp.o
siesta_forces.o: compute_dm.o m_convergence.o m_spin.o m_steps.o m_ts_options.o
siesta_forces.o: mixer.o parallel.o post_scf_work.o scfconvergence_test.o
siesta_forces.o: setup_H0.o setup_hamiltonian.o siesta_cml.o siesta_options.o
siesta_forces.o: sparse_matrices.o state_analysis.o state_init.o sys.o
siesta_geom.o: precision.o
siesta_init.o: alloc.o atomlist.o bands.o densematrix.o files.o ioxv.o
siesta_init.o: kpoint_grid.o kpoint_pdos.o ksvinit.o m_check_supercell.o
siesta_init.o: m_energies.o m_eo.o m_forces.o m_gamma.o m_iostruct.o
siesta_init.o: m_mpi_utils.o m_rmaxh.o m_spin.o m_steps.o m_ts_in_siesta.o
siesta_init.o: m_ts_options.o m_wallclock.o metaforce.o molecularmechanics.o
siesta_init.o: parallel.o parallelsubs.o phonon.o projected_DOS.o
siesta_init.o: siesta_cmlsubs.o siesta_geom.o siesta_options.o
siesta_init.o: sparse_matrices.o struct_init.o sys.o timestamp.o writewave.o
siesta_init.o: xcmod.o zmatrix.o
siesta_move.o: atomlist.o broyden_optim.o cell_broyden_optim.o
siesta_move.o: cell_fire_optim.o dynamics.o fire_optim.o iopipes.o ioxv.o
siesta_move.o: m_energies.o m_forces.o m_kinetic.o m_steps.o m_stress.o
siesta_move.o: parallel.o phonon.o siesta_cml.o siesta_geom.o siesta_options.o
siesta_move.o: sys.o write_subs.o zm_broyden_optim.o zm_fire_optim.o zmatrix.o
siesta_options.o: files.o memoryinfo.o parallel.o precision.o siesta_cml.o
siesta_options.o: sys.o units.o
sparse_matrices.o: alloc.o precision.o
spatial.o: precision.o
spher_harm.o: alloc.o precision.o sys.o
state_analysis.o: atomlist.o born_charge.o fixed.o m_energies.o m_forces.o
state_analysis.o: m_iostruct.o m_ntm.o m_spin.o m_stress.o m_wallclock.o
state_analysis.o: parallel.o remove_intramol_pressure.o siesta_cml.o
state_analysis.o: siesta_geom.o siesta_options.o sparse_matrices.o write_subs.o
state_analysis.o: zmatrix.o
state_init.o: alloc.o atomlist.o files.o hsparse.o iodm_netcdf.o
state_init.o: iodmhs_netcdf.o kpoint_grid.o m_check_supercell.o m_energies.o
state_init.o: m_eo.o m_gamma.o m_mpi_utils.o m_pulay.o m_rmaxh.o m_spin.o
state_init.o: m_steps.o m_ts_io.o m_ts_options.o new_dm.o overlap.o parallel.o
state_init.o: proximity_check.o siesta_cml.o siesta_cmlsubs.o siesta_geom.o
state_init.o: siesta_options.o sparse_matrices.o sys.o write_subs.o zmatrix.o
struct_init.o: alloc.o atmfuncs.o atomlist.o files.o iopipes.o ioxv.o
struct_init.o: m_iostruct.o parallel.o periodic_table.o siesta_cml.o
struct_init.o: siesta_geom.o siesta_options.o zmatrix.o
sys.o: parallel.o
timer.o: m_walltime.o parallel.o precision.o sys.o
transition_rate.o: alloc.o fermid.o parallel.o parallelsubs.o precision.o sys.o
typecell.o: precision.o
uncell.o: precision.o
units.o: precision.o
vmat.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o meshmatrix.o
vmat.o: meshsubs.o meshsubs.o parallel.o parallelsubs.o precision.o
vmatsp.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o meshmatrix.o
vmatsp.o: meshsubs.o meshsubs.o precision.o
vmb.o: fixed.o parallel.o precision.o sys.o
write_md_record.o: atomlist.o m_energies.o m_kinetic.o m_steps.o m_stress.o
write_md_record.o: md_out.o parallel.o siesta_geom.o siesta_options.o
write_raw_efs.o: atmfuncs.o atomlist.o precision.o siesta_geom.o
write_subs.o: alloc.o atomlist.o m_energies.o m_forces.o m_iostruct.o m_spin.o
write_subs.o: m_steps.o m_stress.o m_ts_options.o parallel.o precision.o
write_subs.o: siesta_cml.o siesta_geom.o siesta_options.o units.o zmatrix.o
writewave.o: alloc.o atmfuncs.o atomlist.o densematrix.o files.o kpoint_grid.o
writewave.o: parallel.o parallelsubs.o precision.o siesta_geom.o sys.o units.o
xc.o: alloc.o precision.o sys.o xcmod.o
xcmod.o: parallel.o precision.o sys.o
xml.o: precision.o
zm_broyden_optim.o: alloc.o m_broyddj_nocomm.o parallel.o precision.o sys.o
zm_broyden_optim.o: units.o zmatrix.o
zm_fire_optim.o: m_fdf_global.o m_fire.o parallel.o precision.o
zm_fire_optim.o: siesta_options.o units.o zmatrix.o
zmatrix.o: alloc.o m_cell.o parallel.o precision.o siesta_geom.o sys.o units.o
band.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/bands.o
diagmemory.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/memoryinfo.o
listsc_module.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/listsc.o
m_born_charge.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/born_charge.o
m_broyden_optim.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/broyden_optim.o
m_cell_broyden_optim.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/cell_broyden_optim.o
m_cell_fire_optim.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/cell_fire_optim.o
m_chempot.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/chempot.o
m_compute_dm.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/compute_dm.o
m_conjgr.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/conjgr.o
m_conjgr_old.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/conjgr_old.o
m_dnaefs.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/dnaefs.o
m_doping.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/doping.o
m_dscfcomm.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/dscfcomm.o
m_dynamics.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/dynamics.o
m_fermid.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/fermid.o
m_final_h_f_stress.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/final_H_f_stress.o
m_find_kgrid.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/find_kgrid.o
m_fire_optim.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/fire_optim.o
m_fixed.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/fixed.o
m_forhar.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/forhar.o
m_hsparse.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/hsparse.o
m_intramol_pressure.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/remove_intramol_pressure.o
m_ioxv.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/ioxv.o
m_kinefsm.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/kinefsm.o
m_ksv.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/ksv.o
m_ksvinit.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/ksvinit.o
m_local_dos.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/local_DOS.o
m_minvec.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/minvec.o
m_mixer.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/mixer.o
m_naefs.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/naefs.o
m_new_dm.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/new_dm.o
m_nlefsm.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/nlefsm.o
m_norm.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/compute_norm.o
m_normalize_dm.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/normalize_dm.o
m_ordern.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/ordern.o
m_overfsm.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/overfsm.o
m_overlap.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/overlap.o
m_post_scf_work.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/post_scf_work.o
m_projected_dos.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/projected_DOS.o
m_proximity_check.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/proximity_check.o
m_radfft.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/radfft.o
m_recipes.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/recipes.o
m_scf_options.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/m_options.o
m_scfconvergence_test.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/scfconvergence_test.o
m_setup_h0.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/setup_H0.o
m_setup_hamiltonian.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/setup_hamiltonian.o
m_siesta_analysis.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/siesta_analysis.o
m_siesta_end.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/siesta_end.o
m_siesta_forces.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/siesta_forces.o
m_siesta_init.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/siesta_init.o
m_siesta_move.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/siesta_move.o
m_state_analysis.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/state_analysis.o
m_state_init.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/state_init.o
m_struct_init.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/struct_init.o
m_timestamp.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/timestamp.o
m_xyz_phiatm.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/matel.o
m_zm_broyden_optim.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/zm_broyden_optim.o
m_zm_fire_optim.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/zm_fire_optim.o
mesh.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/meshsubs.o
meshdscf.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/meshmatrix.o
meshphi.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/meshsubs.o
neighbour.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/mneighb.o
on_core.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/onmod.o
on_main.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/onmod.o
printmat.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/printmatrix.o
tiempo.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/timer.o
version_info.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/version.o
write_subs_energies.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/write_subs.o
write_subs_positions.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/write_subs.o
write_subs_pressure.o: /Users/ag/code/SIESTA/BSC/bsc-trunk/Util/TBTrans/../../Src/write_subs.o
