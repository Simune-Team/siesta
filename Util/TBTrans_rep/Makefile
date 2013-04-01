#
# Makefile for tbtrans
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#
# Please see note about Linear-Algebra libraries below
#

default: what tbtrans 

VPATH:=$(shell pwd)/../../Src
OBJDIR=Obj
ARCH_MAKE=../../$(OBJDIR)/arch.make
include $(ARCH_MAKE)
# For removing any -ipo compilation in sub directories...
# Later compilers are also having this enabled, and it lets you
# specify what to not be used in libraries.
# This is so because inter-procedural compilation should NOT be performed in
# libraries (it removes names when a routine can be fully moved into others)
IPO_FLAG?= -ipo

# If compiling without system lapack, you might need something like this
# 
#TBTRANS_LINALG=linalg.a # Or maybe only: zgeev_lapack.a zgesv_lapack.a
LINALG:= $(TBTRANS_LINALG) $(LINALG)

SIESTA_SYS=$(SIESTA_ARCH)

#
# Uncomment the following line for debugging support
#
#FFLAGS=$(FFLAGS_DEBUG) -g

what:
	@echo
	@echo "Compilation architecture to be used:  $(SIESTA_SYS)"
	@echo "If that is not what you want, give the correct"
	@echo "value to the variable SIESTA_SYS in your shell environment."
	@echo $(COMMENTS)
	@echo

SYSOBJ=$(SYS).o

# Tbtrans:
MOD_TOBJS=  m_tbt_options.o m_tbt_kpoints.o m_tbt_out.o m_tbt_iotshs.o m_lowdin.o

SIESTA_OBJS= precision.o sys.o m_io.o alloc.o  memoryinfo.o memory.o pxf.o \
             moreParallelSubs.o m_walltime.o m_timer.o timer.o m_mpi_utils.o \
             m_fdf_global.o paste.o reclat.o idiag.o minvec.o \
             parallel.o io.o files.o cellsubs.o sorting.o m_hs_matrix.o \
	     kpoint_convert.o timestamp.o m_wallclock.o cdiag.o m_spin.o \
	     zheevds.o debugmpi.o \
	     m_ts_contour.o m_ts_gf.o m_ts_scattering.o \
	     m_ts_aux_rout.o m_ts_io.o m_ts_electrode.o m_ts_cctype.o

# Add specific class objects
SIESTA_OBJS += die.o m_uuid.o object_debug.o class_Sparsity.o \
        class_OrbitalDistribution.o \
	class_lsdczSpArr1D.o class_lsdczSpArr2D.o \
	class_lsdczArray1D.o class_lsdczArray2D.o


TBTOBJS= $(SIESTA_OBJS) $(MOD_TOBJS) $(SYSOBJ) $(EXTRA) \
         tbtrans.o transmission.o tbt_init.o tbt_reinit.o \
	 coop.o tt_eig.o atompdos.o


# Interfaces to libraries
FDF=libfdf.a
FDF_MAKEFILE=$(VPATH)/fdf/makefile
FDF_INCFLAGS:=-I$(VPATH)/fdf -I../MPI $(INCFLAGS)
$(FDF): 
	(mkdir -p fdf ; cd fdf ; $(MAKE) -f $(FDF_MAKEFILE) "FC=$(FC)" "VPATH=$(VPATH)/fdf" \
                          "ARCH_MAKE=../$(ARCH_MAKE)" \
                          "INCFLAGS=$(FDF_INCFLAGS)" "FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" module)

INCFLAGS:=$(INCFLAGS)

MPI_MAKEFILE=$(VPATH)/MPI/Makefile
$(MPI_INTERFACE): 
	(mkdir -p MPI; cd MPI ; $(MAKE) -f $(MPI_MAKEFILE) "FC=$(FC)" "VPATH=$(VPATH)/MPI" \
                          "MAKEFILES=$(MPI_MAKEFILE)" "ARCH_MAKE=../$(ARCH_MAKE)" \
                          "FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" module_built)

LINALG_MAKEFILE=$(VPATH)/Libs/makefile
LINALG_INCFLAGS:=-I$(VPATH)/Libs $(INCFLAGS)
$(LINALG): 
	(mkdir -p Libs; cd Libs ; $(MAKE) -f $(LINALG_MAKEFILE) "FC=$(FC)" "VPATH=$(VPATH)/Libs" \
                          "ARCH_MAKE=../$(ARCH_MAKE)" \
                          "INCFLAGS=$(LINALG_INCFLAGS)" "FFLAGS=$(FFLAGS:$(IPO_FLAG)=)" $(LINALG))

tbtrans: FPPFLAGS+=$(DEFS_PREFIX) -DTBTRANS
tbtrans: DEFS+=$(DEFS_PREFIX) -DTBTRANS
tbtrans: $(MPI_INTERFACE) $(FDF) $(TBTOBJS) $(SYSOBJ) $(LINALG)
	$(FC) $(FFLAGS) $(LDFLAGS) -o tbtrans \
              $(TBTOBJS) $(MPI_INTERFACE)  \
              $(FDF) $(LIBS) $(LINALG)  || echo "***Please see top of Makefile"

clean: 
	@echo "==> Cleaning object, library, and executable files"
	rm -f tbtrans *.o *.a *.mod
	rm -rf ./fdf
	rm -rf ./MPI
	rm -rf ./Libs

# Dependencies
dep:
	-sfmakedepend --depend=obj --modext=o \
		$(VPATH)/*.f $(VPATH)/*.f90 $(VPATH)/*.F $(VPATH)/*.F90 \
		*.f90 *.F90 *.f *.F
	@sed -i -e "/^m_ts_io.o:/s/m_ts_kpoints.o//" Makefile
	@sed -i -e "/^m_ts_io.o:/s/m_glob_sparse.o//" Makefile
	@sed -i -e "/^m_ts_electrode.o:/s/m_ts_kpoints.o//" Makefile

$(TBTOBJS): $(FDF)

# DO NOT DELETE THIS LINE - used by make depend
alloc.o: debugmpi.o m_io.o parallel.o precision.o sys.o
arw.o: alloc.o parallel.o precision.o sys.o
atm_transfer.o: atm_types.o atmparams.o old_atmfuncs.o periodic_table.o
atm_transfer.o: radial.o sys.o
atm_types.o: precision.o radial.o
atmfuncs.o: atm_types.o precision.o radial.o spher_harm.o sys.o
atom.o: atmparams.o atom_options.o basis_specs.o basis_types.o bsc_xcmod.o
atom.o: m_filter.o m_recipes.o old_atmfuncs.o periodic_table.o precision.o
atom.o: pseudopotential.o sys.o
atom_options.o: sys.o
atomlist.o: alloc.o atm_types.o atmfuncs.o parallel.o precision.o siesta_geom.o
atomlwf.o: alloc.o atmfuncs.o mneighb.o onmod.o precision.o spatial.o sys.o
automatic_cell.o: atm_types.o atmfuncs.o parallel.o precision.o units.o
bands.o: alloc.o atmfuncs.o atomlist.o densematrix.o files.o m_gamma.o
bands.o: parallel.o parallelsubs.o precision.o siesta_geom.o sys.o writewave.o
basis_enthalpy.o: atmfuncs.o atomlist.o precision.o siesta_geom.o units.o
basis_io.o: atm_types.o atmparams.o atom_options.o basis_types.o chemical.o
basis_io.o: precision.o pseudopotential.o radial.o sys.o xml.o
basis_specs.o: basis_types.o chemical.o periodic_table.o precision.o
basis_specs.o: pseudopotential.o sys.o
basis_types.o: alloc.o atmparams.o precision.o pseudopotential.o sys.o
bessph.o: precision.o sys.o
bonds.o: alloc.o atmfuncs.o m_recipes.o mneighb.o precision.o sorting.o units.o
born_charge.o: atomlist.o ksv.o ksvinit.o m_mpi_utils.o m_spin.o parallel.o
born_charge.o: siesta_geom.o siesta_options.o sparse_matrices.o
broadcast_basis.o: atm_types.o parallel.o radial.o
broyden_optim.o: m_broyddj_nocomm.o m_memory.o parallel.o precision.o
broyden_optim.o: siesta_options.o sys.o units.o
bsc_cellxc.o: alloc.o bsc_xcmod.o cellxc_mod.o mesh.o moremeshsubs.o parallel.o
bsc_cellxc.o: parallelsubs.o precision.o sys.o
bsc_xcmod.o: parallel.o precision.o sys.o
cart2frac.o: sys.o
cdiag.o: alloc.o m_spin.o memoryinfo.o parallel.o precision.o sys.o
cell_broyden_optim.o: alloc.o m_broyddj_nocomm.o parallel.o precision.o sys.o
cell_broyden_optim.o: units.o zmatrix.o
cell_fire_optim.o: alloc.o m_fdf_global.o m_fire.o m_memory.o parallel.o
cell_fire_optim.o: precision.o siesta_options.o sys.o zmatrix.o
cellsubs.o: precision.o
cellxc_mod.o: bsc_xcmod.o sys.o
cgvc.o: alloc.o conjgr_old.o m_mpi_utils.o parallel.o precision.o sys.o units.o
cgvc_zmatrix.o: alloc.o conjgr.o m_mpi_utils.o parallel.o precision.o sys.o
cgvc_zmatrix.o: units.o zmatrix.o
cgwf.o: egandd.o onmod.o parallel.o precision.o sys.o
check_scf_convergence.o: m_convergence.o m_energies.o m_wallclock.o parallel.o
check_scf_convergence.o: precision.o siesta_cml.o siesta_options.o units.o
check_scf_convergence.o: write_subs.o
chemical.o: parallel.o precision.o sys.o
chempot.o: alloc.o mneighb.o parallel.o parallelsubs.o precision.o sys.o
chkdim.o: sys.o
chkgmx.o: cellsubs.o minvec.o precision.o sys.o
class_Fstack_Pair_Geometry_SpMatrix.o: class_Pair_Geometry_SpMatrix.o
class_Geometry.o: alloc.o
class_Pair_Geometry_SpMatrix.o: class_Geometry.o class_lsdczSpArr2D.o
class_Sparsity.o: alloc.o
class_lsdczSpArr1D.o: class_lsdczArray1D.o class_lsdczArray1D.o
class_lsdczSpArr1D.o: class_lsdczArray1D.o class_lsdczArray1D.o
class_lsdczSpArr1D.o: class_lsdczArray1D.o class_lsdczArray1D.o
class_lsdczSpArr2D.o: class_lsdczArray2D.o class_lsdczArray2D.o
class_lsdczSpArr2D.o: class_lsdczArray2D.o class_lsdczArray2D.o
class_lsdczSpArr2D.o: class_lsdczArray2D.o class_lsdczArray2D.o
coceri.o: files.o periodic_table.o precision.o
compute_EKS_Harris.o: atomlist.o class_lsdczSpArr2D.o dhscf.o m_dipol.o
compute_EKS_Harris.o: m_energies.o m_mpi_utils.o m_ntm.o m_spin.o precision.o
compute_EKS_Harris.o: siesta_geom.o siesta_options.o sparse_matrices.o
compute_dm.o: atomlist.o iodmhs_netcdf.o kpoint_grid.o m_energies.o m_eo.o
compute_dm.o: m_gamma.o m_mpi_utils.o m_rmaxh.o m_spin.o m_steps.o
compute_dm.o: m_ts_global_vars.o m_ts_in_siesta.o m_ts_mem.o
compute_dm.o: m_ts_mem_sparsity.o m_ts_options.o ordern.o parallel.o
compute_dm.o: precision.o siesta_geom.o siesta_options.o sparse_matrices.o
compute_dm.o: sys.o units.o
compute_max_diff.o: atomlist.o m_mpi_utils.o m_spin.o precision.o
compute_max_diff.o: sparse_matrices.o
compute_norm.o: m_mpi_utils.o precision.o sparse_matrices.o
conjgr.o: precision.o
conjgr_old.o: precision.o
coor.o: alloc.o parallel.o precision.o siesta_geom.o sys.o zmatrix.o
coxmol.o: files.o periodic_table.o precision.o
create_Sparsity_SC.o: class_Sparsity.o geom_helper.o intrinsic_missing.o
debugmpi.o: parallel.o
denmat.o: alloc.o globalise.o m_mpi_utils.o onmod.o precision.o
denmatlomem.o: alloc.o globalise.o onmod.o precision.o
densematrix.o: alloc.o precision.o
detover.o: alloc.o parallel.o parallelsubs.o precision.o
dfscf.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o mesh.o meshdscf.o
dfscf.o: meshphi.o parallel.o parallelsubs.o precision.o sys.o
dhscf.o: alloc.o atmfuncs.o bsc_xcmod.o cellxc_mod.o doping_uniform.o files.o
dhscf.o: forhar.o iogrid_netcdf.o m_efield.o m_iorho.o m_partial_charges.o
dhscf.o: m_ts_global_vars.o m_ts_in_siesta.o m_ts_options.o m_ts_voltage.o
dhscf.o: mesh.o meshdscf.o meshsubs.o moremeshsubs.o parallel.o parsing.o
dhscf.o: precision.o siesta_options.o sys.o units.o vmat.o
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
die.o: parallel.o
dipole.o: mesh.o parallel.o precision.o sys.o
dnaefs.o: atmfuncs.o mneighb.o precision.o
domain_decom.o: alloc.o parallel.o precision.o printmatrix.o schecomm.o
domain_decom.o: sparse_matrices.o sys.o
doping_uniform.o: alloc.o m_ntm.o mesh.o parallel.o precision.o sys.o
dynamics.o: alloc.o atomlist.o files.o ioxv.o m_mpi_utils.o parallel.o
dynamics.o: precision.o sys.o units.o
egandd.o: alloc.o denmat.o ener3.o globalise.o gradient.o on_subs.o onmod.o
egandd.o: onmod.o precision.o sys.o
eggbox.o: parallel.o precision.o
electrostatic.o: atm_types.o atmfuncs.o atmparams.o bessph.o m_recipes.o
electrostatic.o: precision.o radfft.o radial.o sys.o
ener3.o: alloc.o globalise.o m_mpi_utils.o onmod.o precision.o
ener3lomem.o: alloc.o globalise.o m_mpi_utils.o onmod.o precision.o
extrapolateSpMatrix.o: class_OrbitalDistribution.o class_Sparsity.o
extrapolateSpMatrix.o: class_lsdczArray2D.o class_lsdczSpArr2D.o
extrapolon.o: parallel.o precision.o sys.o
fermid.o: m_recipes.o parallel.o precision.o sys.o
fft.o: alloc.o fft1d.o m_timer.o mesh.o parallel.o parallelsubs.o precision.o
fft.o: sys.o
fft1d.o: parallel.o precision.o sys.o
final_H_f_stress.o: alloc.o atmfuncs.o atomlist.o dnaefs.o files.o grdsam.o
final_H_f_stress.o: kinefsm.o m_dipol.o m_energies.o m_forces.o m_gamma.o
final_H_f_stress.o: m_hsx.o m_mpi_utils.o m_ntm.o m_spin.o m_steps.o m_stress.o
final_H_f_stress.o: m_ts_global_vars.o m_ts_io.o m_ts_kpoints.o m_ts_options.o
final_H_f_stress.o: metaforce.o molecularmechanics.o naefs.o nlefsm.o
final_H_f_stress.o: normalize_dm.o overfsm.o parallel.o siesta_geom.o
final_H_f_stress.o: siesta_options.o sparse_matrices.o sys.o
find_kgrid.o: alloc.o minvec.o parallel.o precision.o units.o
fire_optim.o: alloc.o m_fdf_global.o m_fire.o m_mpi_utils.o parallel.o
fire_optim.o: precision.o siesta_options.o units.o
fixed.o: alloc.o parallel.o precision.o
flib_spline.o: precision.o
forhar.o: alloc.o mesh.o moremeshsubs.o parallel.o precision.o
geom_helper.o: intrinsic_missing.o
get_target_stress.o: parallel.o precision.o sys.o units.o
globalise.o: alloc.o onmod.o onmod.o precision.o spatial.o
gradient.o: alloc.o globalise.o m_mpi_utils.o onmod.o precision.o
gradientlomem.o: alloc.o globalise.o onmod.o precision.o
grdsam.o: alloc.o dhscf.o m_partial_charges.o parallel.o precision.o
grdsam.o: siesta_geom.o siesta_options.o sys.o
hsparse.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o mneighb.o
hsparse.o: parallel.o parallelsubs.o precision.o siesta_options.o sorting.o
hsparse.o: sparse_matrices.o sys.o
idiag.o: parallel.o sys.o
initatom.o: atmparams.o atom.o atom_options.o basis_io.o basis_specs.o
initatom.o: basis_types.o electrostatic.o old_atmfuncs.o precision.o
initparallel.o: alloc.o class_OrbitalDistribution.o domain_decom.o
initparallel.o: m_fdf_global.o parallel.o parallelsubs.o precision.o
initparallel.o: siesta_options.o sparse_matrices.o spatial.o sys.o
io.o: sys.o
iocg.o: files.o precision.o sys.o
iodm.o: alloc.o files.o parallel.o parallelsubs.o precision.o sys.o
iodm_netcdf.o: alloc.o parallel.o parallelsubs.o precision.o sys.o
iodmhs_netcdf.o: alloc.o parallel.o parallelsubs.o precision.o sys.o
ioeig.o: files.o precision.o siesta_cml.o units.o
iofa.o: files.o precision.o
iogrid_netcdf.o: alloc.o parallel.o precision.o sys.o
iokp.o: files.o precision.o
iolwf.o: alloc.o files.o onmod.o parallel.o parallelsubs.o precision.o sys.o
iomd.o: files.o precision.o
iopipes.o: parallel.o precision.o sys.o
iowfs_netcdf.o: alloc.o parallel.o parallelsubs.o precision.o sys.o
ioxv.o: files.o parallel.o precision.o
iozm.o: files.o parallel.o precision.o siesta_geom.o zmatrix.o
ipack.o: sys.o
kgrid.o: parallel.o precision.o units.o
kgridinit.o: minvec.o parallel.o precision.o sys.o units.o
kinefsm.o: alloc.o atmfuncs.o mneighb.o parallel.o parallelsubs.o precision.o
kpoint_convert.o: precision.o sys.o units.o
kpoint_grid.o: find_kgrid.o minvec.o parallel.o precision.o siesta_cml.o
kpoint_grid.o: siesta_options.o sys.o units.o
kpoint_pdos.o: find_kgrid.o minvec.o parallel.o precision.o siesta_options.o
kpoint_pdos.o: sys.o units.o
ksv.o: alloc.o atmfuncs.o densematrix.o ksvinit.o parallel.o precision.o sys.o
ksvinit.o: alloc.o parallel.o precision.o sys.o
listsc.o: alloc.o
local_DOS.o: atomlist.o dhscf.o files.o kpoint_grid.o m_energies.o m_eo.o
local_DOS.o: m_forces.o m_gamma.o m_ntm.o m_spin.o parallel.o siesta_geom.o
local_DOS.o: siesta_options.o sparse_matrices.o sys.o
m_broyddj.o: alloc.o m_mpi_utils.o parallel.o precision.o sys.o
m_broyddj_nocomm.o: alloc.o m_fdf_global.o parallel.o precision.o sys.o
m_broyden_mixing.o: alloc.o m_broyddj.o m_mpi_utils.o parallel.o precision.o
m_cell.o: precision.o siesta_geom.o units.o
m_check_supercell.o: minvec.o parallel.o precision.o
m_diagon.o: precision.o
m_dipol.o: precision.o
m_dscfcomm.o: alloc.o parallel.o parallelsubs.o precision.o schecomm.o
m_efield.o: atmfuncs.o mesh.o parallel.o precision.o siesta_cml.o siesta_geom.o
m_efield.o: sys.o units.o
m_energies.o: precision.o
m_eo.o: precision.o
m_fdf_global.o: m_mpi_utils.o parallel.o precision.o
m_filter.o: bessph.o m_recipes.o radfft.o sorting.o
m_fire.o: parallel.o precision.o
m_fire_mixing.o: alloc.o m_fdf_global.o m_fire_para.o m_mpi_utils.o parallel.o
m_fire_mixing.o: precision.o
m_fire_para.o: m_mpi_utils.o parallel.o precision.o
m_forces.o: precision.o
m_getopts.o: f2kcli.o
m_glob_sparse.o: parallel.o parallelsubs.o precision.o
m_hs_matrix.o: alloc.o precision.o sys.o
m_hsx.o: atm_types.o atmfuncs.o atomlist.o files.o parallel.o parallelsubs.o
m_hsx.o: precision.o siesta_geom.o sys.o
m_io.o: sys.o
m_iodm.o: alloc.o files.o m_mpi_utils.o parallel.o parallelsubs.o precision.o
m_iodm.o: sys.o
m_iorho.o: alloc.o parallel.o parallelsubs.o precision.o sys.o
m_iostruct.o: alloc.o files.o m_mpi_utils.o parallel.o precision.o
m_iostruct.o: siesta_geom.o sys.o units.o
m_kinetic.o: precision.o
m_memory.o: memoryinfo.o parallel.o sys.o
m_mpi_utils.o: precision.o sys.o
m_new_dm.o: alloc.o atomlist.o class_Fstack_Pair_Geometry_SpMatrix.o
m_new_dm.o: class_Geometry.o class_OrbitalDistribution.o
m_new_dm.o: class_Pair_Geometry_SpMatrix.o class_Sparsity.o
m_new_dm.o: class_lsdczArray2D.o class_lsdczSpArr2D.o files.o m_energies.o
m_new_dm.o: m_mpi_utils.o m_spin.o m_steps.o m_svd.o m_ts_iodm.o m_ts_options.o
m_new_dm.o: parallel.o parallelsubs.o parsing.o precision.o readSpMatrix.o
m_new_dm.o: siesta_geom.o siesta_options.o sparse_matrices.o sys.o units.o
m_options.o: precision.o
m_partial_charges.o: atmfuncs.o atomlist.o mesh.o meshphi.o parallel.o
m_partial_charges.o: precision.o siesta_geom.o siesta_options.o sys.o
m_pulay.o: alloc.o atomlist.o m_spin.o m_svd.o parallel.o precision.o
m_pulay.o: siesta_options.o sparse_matrices.o sys.o
m_recipes.o: precision.o sys.o
m_rmaxh.o: precision.o
m_spin.o: alloc.o parallel.o precision.o
m_stress.o: precision.o
m_svd.o: precision.o
m_timer.o: m_io.o m_walltime.o moreParallelSubs.o parallel.o precision.o sys.o
m_ts_aux_rout.o: precision.o sys.o
m_ts_cctype.o: precision.o
m_ts_contour.o: m_ts_aux_rout.o m_ts_cctype.o parallel.o precision.o sys.o
m_ts_contour.o: units.o
m_ts_electrode.o: files.o m_hs_matrix.o m_ts_aux_rout.o m_ts_cctype.o m_ts_io.o
m_ts_electrode.o:  parallel.o precision.o sys.o units.o
m_ts_gf.o: m_ts_cctype.o m_ts_electrode.o parallel.o precision.o sys.o units.o
m_ts_global_vars.o: precision.o
m_ts_in_siesta.o: files.o m_glob_sparse.o m_hs_matrix.o m_ts_cctype.o
m_ts_in_siesta.o: m_ts_contour.o m_ts_gf.o m_ts_kpoints.o m_ts_options.o
m_ts_in_siesta.o: m_ts_scattering.o parallel.o parallelsubs.o precision.o sys.o
m_ts_in_siesta.o: units.o
m_ts_io.o:  parallel.o precision.o sys.o
m_ts_iodm.o: alloc.o files.o m_mpi_utils.o m_ts_global_vars.o parallel.o
m_ts_iodm.o: parallelsubs.o precision.o sys.o
m_ts_kpoints.o: files.o find_kgrid.o minvec.o parallel.o precision.o
m_ts_kpoints.o: siesta_cml.o siesta_options.o sys.o
m_ts_mem.o: class_OrbitalDistribution.o class_Sparsity.o class_lsdczSpArr1D.o
m_ts_mem.o: class_lsdczSpArr1D.o geom_helper.o intrinsic_missing.o
m_ts_mem.o: m_ts_cctype.o m_ts_contour.o m_ts_gf.o m_ts_kpoints.o
m_ts_mem.o: m_ts_mem_scat.o m_ts_mem_sparsity.o m_ts_options.o parallel.o
m_ts_mem.o: precision.o units.o
m_ts_mem_scat.o: class_Sparsity.o class_lsdczSpArr1D.o class_lsdczSpArr1D.o
m_ts_mem_scat.o: intrinsic_missing.o parallel.o precision.o units.o
m_ts_mem_sparsity.o: class_OrbitalDistribution.o class_Sparsity.o
m_ts_mem_sparsity.o: create_Sparsity_SC.o geom_helper.o m_glob_sparse.o
m_ts_mem_sparsity.o: m_ts_mem_scat.o m_ts_options.o parallel.o
m_ts_options.o: files.o m_ts_contour.o m_ts_global_vars.o m_ts_io.o parallel.o
m_ts_options.o: precision.o siesta_options.o sys.o units.o
m_ts_scattering.o: m_ts_aux_rout.o parallel.o precision.o sys.o units.o
m_ts_voltage.o: m_ts_options.o mesh.o parallel.o precision.o units.o
m_wallclock.o: m_walltime.o
madelung.o: parallel.o precision.o sys.o
matel.o: alloc.o atmfuncs.o m_recipes.o parallel.o precision.o radfft.o
matel.o: spher_harm.o
md_out.o: files.o m_energies.o precision.o sys.o units.o
md_utils.o: precision.o
memory.o: alloc.o memoryinfo.o parallel.o precision.o
memoryinfo.o: precision.o
mesh.o: precision.o
meshcomm.o: mesh.o precision.o
meshdscf.o: alloc.o atomlist.o m_dscfcomm.o meshphi.o parallel.o parallelsubs.o
meshdscf.o: precision.o
meshphi.o: alloc.o precision.o
meshsubs.o: alloc.o atm_types.o atmfuncs.o cellsubs.o cellxc_mod.o chkgmx.o
meshsubs.o: fft1d.o mesh.o meshphi.o moremeshsubs.o parallel.o parallelsubs.o
meshsubs.o: precision.o radial.o siesta_cml.o sys.o
metaforce.o: alloc.o parallel.o precision.o sys.o
minvec.o: cellsubs.o precision.o sorting.o sys.o
mixer.o: atomlist.o m_broyden_mixing.o m_energies.o m_fire_mixing.o m_pulay.o
mixer.o: m_spin.o m_ts_global_vars.o parallel.o precision.o
mixer.o: save_density_matrix.o siesta_geom.o siesta_options.o sparse_matrices.o
mneighb.o: alloc.o precision.o sys.o
molecularmechanics.o: alloc.o parallel.o precision.o sys.o units.o
moreParallelSubs.o: alloc.o m_io.o parallel.o precision.o sys.o
moremeshsubs.o: alloc.o mesh.o parallel.o precision.o schecomm.o sys.o
mulliken.o: alloc.o atmfuncs.o parallel.o parallelsubs.o precision.o
mulliken.o: siesta_cml.o
naefs.o: atmfuncs.o mneighb.o precision.o
nlefsm.o: alloc.o atmfuncs.o mneighb.o parallel.o parallelsubs.o precision.o
normalize_dm.o: atomlist.o m_mpi_utils.o m_spin.o parallel.o precision.o
normalize_dm.o: siesta_options.o sparse_matrices.o sys.o
obc.o: alloc.o precision.o
ofc.o: alloc.o files.o precision.o units.o
old_atmfuncs.o: alloc.o atmparams.o precision.o sys.o
on_subs.o: alloc.o onmod.o onmod.o
onmod.o: precision.o
optical.o: alloc.o atomlist.o densematrix.o fermid.o files.o parallel.o
optical.o: parallelsubs.o precision.o sys.o units.o
ordern.o: alloc.o cgwf.o chempot.o globalise.o m_mpi_utils.o onmod.o parallel.o
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
phirphi.o: alloc.o atmfuncs.o mneighb.o parallel.o parallelsubs.o precision.o
phirphi_opt.o: alloc.o atm_types.o atmfuncs.o atmparams.o mneighb.o parallel.o
phirphi_opt.o: parallelsubs.o precision.o spher_harm.o sys.o
phonon.o: alloc.o files.o parallel.o precision.o sys.o
pixmol.o: files.o periodic_table.o precision.o
plcharge.o: files.o old_atmfuncs.o parallel.o precision.o
poison.o: alloc.o cellsubs.o chkgmx.o fft.o parallel.o precision.o sys.o
post_scf_work.o: atomlist.o final_H_f_stress.o m_energies.o m_spin.o mneighb.o
post_scf_work.o: siesta_geom.o siesta_options.o sparse_matrices.o
printmatrix.o: alloc.o
projected_DOS.o: alloc.o atomlist.o kpoint_grid.o kpoint_pdos.o m_eo.o
projected_DOS.o: m_gamma.o m_spin.o parallel.o precision.o siesta_geom.o
projected_DOS.o: siesta_options.o sparse_matrices.o sys.o
propor.o: precision.o sys.o
proximity_check.o: m_ts_options.o mneighb.o parallel.o precision.o
proximity_check.o: siesta_geom.o siesta_options.o
pseudopotential.o: alloc.o atom_options.o flib_spline.o precision.o sys.o
radfft.o: alloc.o bessph.o m_recipes.o precision.o
radial.o: alloc.o precision.o xml.o
rdiag.o: alloc.o m_diagon.o memoryinfo.o parallel.o precision.o sys.o
readSpMatrix.o: class_OrbitalDistribution.o class_Sparsity.o
readSpMatrix.o: class_lsdczArray2D.o class_lsdczSpArr2D.o
read_xc_info.o: parallel.o precision.o sys.o
readsp.o: parallel.o precision.o sys.o
redcel.o: parallel.o precision.o sys.o
reinit.o: files.o parallel.o
remove_intramol_pressure.o: precision.o zmatrix.o
reoptical.o: parallel.o precision.o sys.o
reord.o: alloc.o precision.o
restructSpMatrix.o: class_OrbitalDistribution.o class_Sparsity.o
restructSpMatrix.o: class_lsdczArray2D.o class_lsdczSpArr2D.o
rhooda.o: atmfuncs.o atomlist.o mesh.o meshphi.o precision.o
rhoofd.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o mesh.o meshdscf.o
rhoofd.o: meshphi.o parallel.o parallelsubs.o precision.o sys.o
rhoofdsp.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o mesh.o
rhoofdsp.o: meshdscf.o meshphi.o precision.o sys.o
save_density_matrix.o: atomlist.o iodm_netcdf.o m_energies.o m_iodm.o m_spin.o
save_density_matrix.o: m_steps.o m_ts_global_vars.o m_ts_iodm.o precision.o
save_density_matrix.o: siesta_options.o sparse_matrices.o
savepsi.o: alloc.o parallel.o parallelsubs.o precision.o
schecomm.o: alloc.o
setatomnodes.o: alloc.o parallel.o precision.o spatial.o sys.o
setspatial.o: alloc.o parallel.o precision.o spatial.o
setup_H0.o: alloc.o atmfuncs.o atomlist.o class_lsdczSpArr2D.o dhscf.o dnaefs.o
setup_H0.o: kinefsm.o m_energies.o m_ntm.o m_spin.o metaforce.o
setup_H0.o: molecularmechanics.o naefs.o nlefsm.o siesta_geom.o
setup_H0.o: siesta_options.o sparse_matrices.o
setup_hamiltonian.o: alloc.o atmfuncs.o atomlist.o class_lsdczSpArr2D.o dhscf.o
setup_hamiltonian.o: m_dipol.o m_energies.o m_gamma.o m_hsx.o m_mpi_utils.o
setup_hamiltonian.o: m_ntm.o m_partial_charges.o m_spin.o m_steps.o m_stress.o
setup_hamiltonian.o: m_ts_global_vars.o metaforce.o molecularmechanics.o
setup_hamiltonian.o: normalize_dm.o parallel.o siesta_geom.o siesta_options.o
setup_hamiltonian.o: sparse_matrices.o sys.o
setup_kscell.o: minvec.o parallel.o precision.o units.o
setup_ordern_indexes.o: alloc.o domain_decom.o parallel.o spatial.o
shaper.o: atmfuncs.o mneighb.o precision.o
show_distribution.o: atomlist.o parallel.o parallelsubs.o siesta_geom.o sys.o
siesta.o: alloc.o m_steps.o siesta_analysis.o siesta_end.o siesta_forces.o
siesta.o: siesta_init.o siesta_move.o
siesta_analysis.o: alloc.o atomlist.o bands.o basis_enthalpy.o dhscf.o files.o
siesta_analysis.o: kpoint_grid.o ksv.o ksvinit.o local_DOS.o m_dipol.o
siesta_analysis.o: m_energies.o m_eo.o m_forces.o m_gamma.o m_mpi_utils.o
siesta_analysis.o: m_ntm.o m_partial_charges.o m_spin.o m_steps.o parallel.o
siesta_analysis.o: projected_DOS.o siesta_cml.o siesta_geom.o siesta_options.o
siesta_analysis.o: sparse_matrices.o write_subs.o writewave.o zmatrix.o
siesta_cmlsubs.o: files.o parallel.o siesta_cml.o timestamp.o version.o
siesta_end.o: alloc.o bands.o densematrix.o fixed.o m_dscfcomm.o m_pulay.o
siesta_end.o: m_wallclock.o meshdscf.o meshphi.o moremeshsubs.o parallel.o
siesta_end.o: siesta_cml.o siesta_cmlsubs.o sparse_matrices.o timestamp.o
siesta_forces.o: atomlist.o check_scf_convergence.o compute_EKS_Harris.o
siesta_forces.o: compute_dm.o compute_max_diff.o final_H_f_stress.o
siesta_forces.o: m_convergence.o m_forces.o m_pulay.o m_spin.o m_steps.o
siesta_forces.o: m_ts_global_vars.o m_ts_mem_sparsity.o m_ts_options.o mixer.o
siesta_forces.o: parallel.o post_scf_work.o precision.o save_density_matrix.o
siesta_forces.o: setup_H0.o setup_hamiltonian.o siesta_cml.o siesta_geom.o
siesta_forces.o: siesta_options.o sparse_matrices.o state_analysis.o
siesta_forces.o: state_init.o sys.o units.o write_subs.o
siesta_geom.o: precision.o
siesta_init.o: alloc.o atomlist.o bands.o bsc_xcmod.o
siesta_init.o: class_Fstack_Pair_Geometry_SpMatrix.o densematrix.o files.o
siesta_init.o: ioxv.o kpoint_grid.o kpoint_pdos.o ksvinit.o m_check_supercell.o
siesta_init.o: m_energies.o m_eo.o m_forces.o m_gamma.o m_iostruct.o
siesta_init.o: m_mpi_utils.o m_new_dm.o m_rmaxh.o m_spin.o m_steps.o m_timer.o
siesta_init.o: m_wallclock.o metaforce.o molecularmechanics.o object_debug.o
siesta_init.o: parallel.o parallelsubs.o phonon.o projected_DOS.o
siesta_init.o: siesta_cmlsubs.o siesta_geom.o siesta_options.o
siesta_init.o: sparse_matrices.o struct_init.o sys.o timestamp.o ts_init.o
siesta_init.o: writewave.o zmatrix.o
siesta_move.o: atomlist.o broyden_optim.o cell_broyden_optim.o
siesta_move.o: cell_fire_optim.o dynamics.o fire_optim.o iopipes.o ioxv.o
siesta_move.o: m_energies.o m_forces.o m_kinetic.o m_steps.o m_stress.o
siesta_move.o: parallel.o phonon.o siesta_cml.o siesta_geom.o siesta_options.o
siesta_move.o: sys.o write_subs.o zm_broyden_optim.o zm_fire_optim.o zmatrix.o
siesta_options.o: files.o memoryinfo.o parallel.o precision.o siesta_cml.o
siesta_options.o: sys.o units.o
sparse_matrices.o: alloc.o class_Fstack_Pair_Geometry_SpMatrix.o
sparse_matrices.o: class_OrbitalDistribution.o class_Sparsity.o
sparse_matrices.o: class_lsdczSpArr2D.o precision.o
spatial.o: precision.o
spher_harm.o: alloc.o precision.o sys.o
state_analysis.o: atomlist.o born_charge.o
state_analysis.o: class_Fstack_Pair_Geometry_SpMatrix.o class_Geometry.o
state_analysis.o: class_Pair_Geometry_SpMatrix.o fixed.o m_energies.o
state_analysis.o: m_forces.o m_iostruct.o m_ntm.o m_spin.o m_stress.o
state_analysis.o: m_wallclock.o parallel.o remove_intramol_pressure.o
state_analysis.o: siesta_cml.o siesta_geom.o siesta_options.o sparse_matrices.o
state_analysis.o: write_subs.o zmatrix.o
state_init.o: alloc.o atomlist.o class_Sparsity.o class_lsdczSpArr2D.o
state_init.o: domain_decom.o files.o hsparse.o iodm_netcdf.o iodmhs_netcdf.o
state_init.o: kpoint_grid.o m_check_supercell.o m_energies.o m_eo.o m_gamma.o
state_init.o: m_mpi_utils.o m_new_dm.o m_pulay.o m_rmaxh.o m_spin.o m_steps.o
state_init.o: m_ts_global_vars.o m_ts_io.o m_ts_options.o overlap.o parallel.o
state_init.o: proximity_check.o siesta_cml.o siesta_geom.o siesta_options.o
state_init.o: sparse_matrices.o sys.o write_subs.o zmatrix.o
struct_init.o: alloc.o atmfuncs.o atomlist.o files.o iopipes.o ioxv.o
struct_init.o: m_iostruct.o m_mpi_utils.o parallel.o periodic_table.o
struct_init.o: siesta_cml.o siesta_geom.o siesta_options.o zmatrix.o
sys.o: parallel.o
timer.o: m_timer.o sys.o
transition_rate.o: alloc.o fermid.o parallel.o parallelsubs.o precision.o sys.o
ts_init.o: files.o m_ts_contour.o m_ts_gf.o m_ts_kpoints.o m_ts_options.o
ts_init.o: parallel.o
ts_show_regions.o: parallel.o precision.o units.o
typecell.o: precision.o
uncell.o: precision.o
units.o: precision.o
vacuum_level.o: moreParallelSubs.o precision.o
vmat.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o mesh.o meshdscf.o
vmat.o: meshphi.o parallel.o parallelsubs.o precision.o
vmatsp.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o mesh.o meshdscf.o
vmatsp.o: meshphi.o precision.o
vmb.o: fixed.o parallel.o precision.o sys.o
write_md_record.o: atomlist.o m_energies.o m_kinetic.o m_steps.o m_stress.o
write_md_record.o: md_out.o parallel.o siesta_geom.o siesta_options.o
write_orb_indx.o: atmfuncs.o cellsubs.o files.o precision.o
write_raw_efs.o: atmfuncs.o atomlist.o precision.o siesta_geom.o
write_subs.o: alloc.o atomlist.o m_energies.o m_forces.o m_iostruct.o m_spin.o
write_subs.o: m_steps.o m_stress.o m_ts_global_vars.o parallel.o precision.o
write_subs.o: siesta_cml.o siesta_geom.o siesta_options.o units.o zmatrix.o
writewave.o: alloc.o atmfuncs.o atomlist.o densematrix.o kpoint_grid.o
writewave.o: parallel.o parallelsubs.o precision.o siesta_geom.o sys.o units.o
xc.o: alloc.o bsc_xcmod.o precision.o sys.o
xml.o: precision.o
zm_broyden_optim.o: alloc.o m_broyddj_nocomm.o parallel.o precision.o sys.o
zm_broyden_optim.o: units.o zmatrix.o
zm_fire_optim.o: m_fdf_global.o m_fire.o parallel.o precision.o
zm_fire_optim.o: siesta_options.o units.o zmatrix.o
zmatrix.o: alloc.o m_cell.o parallel.o precision.o siesta_geom.o sys.o units.o
atompdos.o: m_lowdin.o m_tbt_out.o parallel.o precision.o sys.o units.o
coop.o: m_tbt_out.o parallel.o precision.o sys.o units.o
m_lowdin.o: precision.o
m_tbt_iotshs.o: m_tbt_kpoints.o m_ts_io.o parallel.o precision.o sys.o
m_tbt_kpoints.o: alloc.o minvec.o parallel.o precision.o sys.o units.o
m_tbt_options.o: files.o m_ts_io.o memoryinfo.o parallel.o precision.o sys.o
m_tbt_options.o: units.o
m_tbt_out.o: parallel.o precision.o units.o
tbt_init.o: alloc.o files.o m_spin.o m_tbt_kpoints.o m_tbt_options.o m_timer.o
tbt_init.o: m_ts_contour.o m_wallclock.o parallel.o precision.o sys.o
tbt_init.o: timestamp.o
tbt_reinit.o: files.o parallel.o
tbtrans.o: alloc.o files.o m_hs_matrix.o m_tbt_iotshs.o m_tbt_kpoints.o
tbtrans.o: m_tbt_options.o m_tbt_out.o m_ts_contour.o m_ts_electrode.o
tbtrans.o: m_ts_gf.o m_ts_scattering.o m_wallclock.o parallel.o precision.o
tbtrans.o: sys.o timestamp.o units.o
transmission.o: m_ts_aux_rout.o precision.o sys.o
tt_eig.o: precision.o sys.o
band.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/bands.o
class_carray1d.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/class_lsdczArray1D.o
class_carray2d.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/class_lsdczArray2D.o
class_csparr1d.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/class_lsdczSpArr1D.o
class_csparr2d.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/class_lsdczSpArr2D.o
class_darray1d.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/class_lsdczArray1D.o
class_darray2d.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/class_lsdczArray2D.o
class_dsparr1d.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/class_lsdczSpArr1D.o
class_dsparr2d.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/class_lsdczSpArr2D.o
class_fstack_pair_geometry_dsparr2d.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/class_Fstack_Pair_Geometry_SpMatrix.o
class_iarray1d.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/class_lsdczArray1D.o
class_iarray2d.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/class_lsdczArray2D.o
class_isparr1d.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/class_lsdczSpArr1D.o
class_isparr2d.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/class_lsdczSpArr2D.o
class_larray1d.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/class_lsdczArray1D.o
class_larray2d.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/class_lsdczArray2D.o
class_lsparr1d.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/class_lsdczSpArr1D.o
class_lsparr2d.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/class_lsdczSpArr2D.o
class_pair_geometry_dsparr2d.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/class_Pair_Geometry_SpMatrix.o
class_sarray1d.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/class_lsdczArray1D.o
class_sarray2d.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/class_lsdczArray2D.o
class_ssparr1d.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/class_lsdczSpArr1D.o
class_ssparr2d.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/class_lsdczSpArr2D.o
class_zarray1d.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/class_lsdczArray1D.o
class_zarray2d.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/class_lsdczArray2D.o
class_zsparr1d.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/class_lsdczSpArr1D.o
class_zsparr2d.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/class_lsdczSpArr2D.o
diagmemory.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/memoryinfo.o
listsc_module.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/listsc.o
m_bessph.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/bessph.o
m_born_charge.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/born_charge.o
m_broyden_optim.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/broyden_optim.o
m_cell_broyden_optim.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/cell_broyden_optim.o
m_cell_fire_optim.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/cell_fire_optim.o
m_cgwf.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/cgwf.o
m_check_scf_convergence.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/check_scf_convergence.o
m_chempot.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/chempot.o
m_chkgmx.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/chkgmx.o
m_compute_dm.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/compute_dm.o
m_compute_max_diff.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/compute_max_diff.o
m_conjgr.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/conjgr.o
m_conjgr_old.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/conjgr_old.o
m_denmat.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/denmat.o
m_dhscf.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/dhscf.o
m_dnaefs.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/dnaefs.o
m_doping_uniform.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/doping_uniform.o
m_dynamics.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/dynamics.o
m_eandg.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/egandd.o
m_eks_harris.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/compute_EKS_Harris.o
m_ener3.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/ener3.o
m_fermid.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/fermid.o
m_fft.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/fft.o
m_final_h_f_stress.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/final_H_f_stress.o
m_find_kgrid.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/find_kgrid.o
m_fire_optim.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/fire_optim.o
m_fixed.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/fixed.o
m_forhar.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/forhar.o
m_gradient.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/gradient.o
m_grdsam.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/grdsam.o
m_hsparse.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/hsparse.o
m_intramol_pressure.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/remove_intramol_pressure.o
m_ioxv.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/ioxv.o
m_kinefsm.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/kinefsm.o
m_ksv.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/ksv.o
m_ksvinit.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/ksvinit.o
m_local_dos.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/local_DOS.o
m_minvec.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/minvec.o
m_mixer.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/mixer.o
m_naefs.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/naefs.o
m_nlefsm.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/nlefsm.o
m_norm.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/compute_norm.o
m_normalize_dm.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/normalize_dm.o
m_object_debug.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/object_debug.o
m_on_subs.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/on_subs.o
m_ordern.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/ordern.o
m_overfsm.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/overfsm.o
m_overlap.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/overlap.o
m_post_scf_work.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/post_scf_work.o
m_projected_dos.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/projected_DOS.o
m_proximity_check.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/proximity_check.o
m_radfft.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/radfft.o
m_readdsparr2d.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/readSpMatrix.o
m_save_density_matrix.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/save_density_matrix.o
m_scf_options.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/m_options.o
m_setup_h0.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/setup_H0.o
m_setup_hamiltonian.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/setup_hamiltonian.o
m_siesta_analysis.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/siesta_analysis.o
m_siesta_end.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/siesta_end.o
m_siesta_forces.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/siesta_forces.o
m_siesta_init.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/siesta_init.o
m_siesta_move.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/siesta_move.o
m_state_analysis.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/state_analysis.o
m_state_init.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/state_init.o
m_struct_init.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/struct_init.o
m_timestamp.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/timestamp.o
m_ts_init.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/ts_init.o
m_vmat.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/vmat.o
m_xyz_phiatm.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/matel.o
m_zm_broyden_optim.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/zm_broyden_optim.o
m_zm_fire_optim.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/zm_fire_optim.o
mtprng.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/m_uuid.o
neighbour.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/mneighb.o
on_core.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/onmod.o
on_main.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/onmod.o
printmat.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/printmatrix.o
version_info.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/version.o
write_subs_energies.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/write_subs.o
write_subs_positions.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/write_subs.o
write_subs_pressure.o: /home/nick/dtu/siesta-ts-scf/Util/TBTrans_rep/../../Src/write_subs.o
