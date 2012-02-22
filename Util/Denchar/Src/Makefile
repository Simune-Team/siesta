#
# Makefile for stand-alone Denchar
#
# The idea is to use the code in the top Src directory as much as possible.
# This is achieved by the VPATH directive below.
# Other points to note, until we switch to a better building system:
#
#  The arch.make file is supposed to be in $(OBJDIR). This is normally
#  the top Obj, but if you are using architecture-dependent build directories
#  you might want to change this. (If you do not understand this, you do not
#  need to change anything. Power users can do "make OBJDIR=Whatever".)
#
#  If your main Siesta build used an mpi compiler, you might need to
#  define an FC_SERIAL symbol in your top arch.make, to avoid linking
#  in the mpi libraries even if we explicitly undefine MPI below.
#  
#  The dependency list at the end is overly large, but harmless
#
OBJDIR=Obj
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#
TOP_LEVEL:=$(shell pwd)/../../../Src
VPATH:=.:$(TOP_LEVEL)
#
default: what denchar
#
dep:
	sfmakedepend --depend=obj  --modext=o \
          $(TOP_LEVEL)/*.f $(TOP_LEVEL)/*.f90 $(TOP_LEVEL)/*.F $(TOP_LEVEL)/*.F90 
	@sed '/^siesta_cmlsubs.o:/s/version.o//'  Makefile > tmp.Makefile
	@mv tmp.Makefile Makefile
#
#
ARCH_MAKE=../../../$(OBJDIR)/arch.make
include $(ARCH_MAKE)
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
DEFS:=$(DEFS) $(DEFS_PREFIX)-UMPI $(DEFS_PREFIX)-UCDF
FPPFLAGS:=$(FPPFLAGS) $(DEFS_PREFIX)-UMPI $(DEFS_PREFIX)-UCDF

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
OBJS_DENCHAR = atompla.o planed.o\
        bessph.o chkdim.o dismin.o \
	dot.o readpla.o readsts.o \
        io.o iodm.o itochar.o length.o \
        matvect.o memory.o colinear.o\
	neighb.o paste.o pxf.o \
	m_recipes.o radfft.o ranger.o \
        readwaves.o rhoofr.o wavofr.o \
        redata_denchar.o  local_reinit.o denchar.o \
	volcel.o wrout.o stsofr.o
#
# Use the makefile in Src/fdf and all the sources there.
#
FDF=libfdf.a
FDF_MAKEFILE=$(TOP_LEVEL)/fdf/makefile
FDF_INCFLAGS:=-I $(TOP_LEVEL)/fdf $(INCFLAGS)
$(FDF): 
	(mkdir -p fdf ; cd fdf ; $(MAKE) -f $(FDF_MAKEFILE) "FC=$(FC)" "VPATH=$(TOP_LEVEL)/fdf" \
                          "ARCH_MAKE=../$(ARCH_MAKE)" \
                          "DEFS=$(DEFS)" \
                          "FPPFLAGS=$(FPPFLAGS)" \
                          "MPI_INTERFACE= " \
                          "INCFLAGS=$(FDF_INCFLAGS)" "FFLAGS=$(FFLAGS)" module)
#
#
# This is crude but will have to do for now.
# Note : precision must be the first module
#
MOD_OBJS=precision.o parallel.o debugmpi.o sys.o m_io.o alloc.o memoryinfo.o listsc.o  \
         atmparams.o atom_options.o atmfuncs.o   atm_types.o \
         radial.o spher_harm.o  basis_io.o basis_types.o parallelsubs.o \
         schecomm.o sparse_matrices.o \
         files.o flib_spline.o pseudopotential.o chemical.o xml.o \
         domain_decom.o printmatrix.o qsort.o mmio.o pspltm1.o lenstr.o \
         spatial.o timer_local.o 
#
COM_OBJS_DENCHAR=$(OBJS_DENCHAR) $(SYSOBJ)
ALL_OBJS_DENCHAR=$(MOD_OBJS) $(COM_OBJS_DENCHAR)
#

denchar: $(FDF)  $(ALL_OBJS_DENCHAR)
	$(FC) -o denchar \
	       $(LDFLAGS) $(ALL_OBJS_DENCHAR) $(FDF) $(METIS_LIB)
#

clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f denchar *.o  *.a
	rm -f *.mod
	(cd fdf ; $(MAKE) -f $(FDF_MAKEFILE) "ARCH_MAKE=../$(ARCH_MAKE)" clean)

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
automatic_cell.o: atmfuncs.o parallel.o precision.o units.o
bands.o: alloc.o densematrix.o files.o parallel.o parallelsubs.o precision.o
bands.o: sys.o
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
chemical.o: parallel.o precision.o sys.o
chempot.o: alloc.o mneighb.o parallel.o parallelsubs.o precision.o sys.o
chkdim.o: sys.o
chkgmx.o: cellsubs.o minvec.o precision.o sys.o
coceri.o: files.o periodic_table.o precision.o
compute_dm.o: atomlist.o iodmhs_netcdf.o kpoint_grid.o m_energies.o m_eo.o
compute_dm.o: m_gamma.o m_mpi_utils.o m_rmaxh.o m_spin.o m_steps.o
compute_dm.o: m_ts_global_vars.o m_ts_in_siesta.o m_ts_options.o ordern.o
compute_dm.o: parallel.o precision.o siesta_geom.o siesta_options.o
compute_dm.o: sparse_matrices.o sys.o units.o
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
dfscf.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o mesh.o meshdscf.o
dfscf.o: meshphi.o parallel.o parallelsubs.o precision.o sys.o
dhscf.o: alloc.o atmfuncs.o bsc_xcmod.o cellxc_mod.o doping_uniform.o files.o
dhscf.o: forhar.o iogrid_netcdf.o m_efield.o m_iorho.o m_ts_global_vars.o
dhscf.o: m_ts_in_siesta.o m_ts_options.o mesh.o meshdscf.o meshsubs.o
dhscf.o: moremeshsubs.o parallel.o parsing.o precision.o siesta_options.o sys.o
dhscf.o: units.o vmat.o
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
get_target_stress.o: parallel.o precision.o sys.o units.o
globalise.o: alloc.o onmod.o onmod.o precision.o spatial.o
gradient.o: alloc.o globalise.o m_mpi_utils.o onmod.o precision.o
gradientlomem.o: alloc.o globalise.o onmod.o precision.o
grdsam.o: alloc.o dhscf.o parallel.o precision.o siesta_geom.o siesta_options.o
grdsam.o: sys.o
hsparse.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o mneighb.o
hsparse.o: parallel.o parallelsubs.o precision.o siesta_options.o sorting.o
hsparse.o: sparse_matrices.o sys.o
idiag.o: parallel.o sys.o
initatom.o: atmparams.o atom.o atom_options.o basis_io.o basis_specs.o
initatom.o: basis_types.o electrostatic.o old_atmfuncs.o precision.o
initparallel.o: alloc.o domain_decom.o m_fdf_global.o parallel.o parallelsubs.o
initparallel.o: precision.o siesta_options.o spatial.o sys.o
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
m_history.o: alloc.o
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
m_options.o: precision.o
m_pulay.o: alloc.o atomlist.o m_spin.o parallel.o precision.o siesta_options.o
m_pulay.o: sparse_matrices.o sys.o
m_recipes.o: precision.o sys.o
m_rmaxh.o: precision.o
m_sparse.o: alloc.o precision.o
m_spin.o: alloc.o parallel.o precision.o
m_stress.o: precision.o
m_timer.o: m_io.o m_walltime.o moreParallelSubs.o parallel.o precision.o sys.o
m_ts_aux_rout.o: precision.o sys.o
m_ts_contour.o: files.o m_ts_aux_rout.o m_ts_options.o precision.o sys.o
m_ts_electrode.o: files.o m_ts_aux_rout.o m_ts_io.o m_ts_kpoints.o
m_ts_electrode.o: m_ts_options.o precision.o
m_ts_global_vars.o: precision.o
m_ts_in_siesta.o: find_kgrid.o m_mpi_utils.o m_ts_contour.o m_ts_electrode.o
m_ts_in_siesta.o: m_ts_io.o m_ts_kpoints.o m_ts_options.o m_ts_scattering.o
m_ts_in_siesta.o: mesh.o parallel.o parallelsubs.o precision.o sys.o
m_ts_io.o: files.o m_spin.o m_ts_kpoints.o parallel.o parallelsubs.o
m_ts_io.o: precision.o sys.o
m_ts_iodm.o: alloc.o files.o m_mpi_utils.o m_ts_global_vars.o parallel.o
m_ts_iodm.o: parallelsubs.o precision.o sys.o
m_ts_kpoints.o: minvec.o parallel.o precision.o siesta_cml.o siesta_options.o
m_ts_kpoints.o: sys.o units.o
m_ts_options.o: m_fdf_global.o m_ts_global_vars.o parallel.o precision.o
m_ts_options.o: siesta_options.o sys.o units.o
m_ts_scattering.o: m_ts_aux_rout.o m_ts_options.o parallelsubs.o precision.o
m_ts_scattering.o: sys.o
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
mixer.o: atomlist.o iodm_netcdf.o m_broyden_mixing.o m_energies.o
mixer.o: m_fire_mixing.o m_iodm.o m_mpi_utils.o m_pulay.o m_spin.o m_steps.o
mixer.o: m_ts_global_vars.o m_ts_iodm.o m_ts_options.o parallel.o siesta_geom.o
mixer.o: siesta_options.o sparse_matrices.o
mneighb.o: alloc.o precision.o sys.o
molecularmechanics.o: alloc.o parallel.o precision.o sys.o units.o
moreParallelSubs.o: alloc.o m_io.o parallel.o precision.o sys.o
moremeshsubs.o: alloc.o mesh.o parallel.o precision.o schecomm.o sys.o
mulliken.o: alloc.o atmfuncs.o parallel.o parallelsubs.o precision.o
mulliken.o: siesta_cml.o
naefs.o: atmfuncs.o mneighb.o precision.o
new_dm.o: alloc.o atomlist.o m_energies.o m_fdf_global.o m_iodm.o m_sparse.o
new_dm.o: m_spin.o m_steps.o m_ts_global_vars.o m_ts_iodm.o m_ts_options.o
new_dm.o: parallel.o parallelsubs.o parsing.o precision.o siesta_geom.o
new_dm.o: siesta_options.o sparse_matrices.o sys.o units.o
nlefsm.o: alloc.o atmfuncs.o mneighb.o parallel.o parallelsubs.o precision.o
normalize_dm.o: atomlist.o m_mpi_utils.o m_spin.o parallel.o precision.o
normalize_dm.o: sparse_matrices.o
obc.o: alloc.o precision.o
ofc.o: alloc.o files.o precision.o units.o
old_atmfuncs.o: alloc.o atmparams.o precision.o sys.o
on_subs.o: alloc.o onmod.o onmod.o
onmod.o: precision.o
optical.o: alloc.o atmfuncs.o densematrix.o fermid.o files.o parallel.o
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
read_xc_info.o: parallel.o precision.o sys.o
readsp.o: parallel.o precision.o sys.o
redcel.o: parallel.o precision.o sys.o
reinit.o: files.o parallel.o
remove_intramol_pressure.o: precision.o zmatrix.o
reoptical.o: parallel.o precision.o sys.o
reord.o: alloc.o precision.o
rhooda.o: atmfuncs.o atomlist.o mesh.o meshphi.o precision.o
rhoofd.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o mesh.o meshdscf.o
rhoofd.o: meshphi.o parallel.o parallelsubs.o precision.o sys.o
rhoofdsp.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o mesh.o
rhoofdsp.o: meshdscf.o meshphi.o precision.o sys.o
savepsi.o: alloc.o parallel.o parallelsubs.o precision.o
scfconvergence_test.o: m_convergence.o m_energies.o m_pulay.o
scfconvergence_test.o: m_ts_global_vars.o m_ts_options.o m_wallclock.o
scfconvergence_test.o: parallel.o siesta_cml.o siesta_options.o write_subs.o
schecomm.o: alloc.o
setatomnodes.o: alloc.o parallel.o precision.o spatial.o sys.o
setspatial.o: alloc.o parallel.o precision.o spatial.o
setup_H0.o: alloc.o atmfuncs.o atomlist.o dhscf.o dnaefs.o kinefsm.o
setup_H0.o: m_energies.o m_ntm.o m_spin.o metaforce.o molecularmechanics.o
setup_H0.o: naefs.o nlefsm.o siesta_geom.o siesta_options.o sparse_matrices.o
setup_hamiltonian.o: alloc.o atmfuncs.o atomlist.o dhscf.o m_dipol.o
setup_hamiltonian.o: m_energies.o m_gamma.o m_hsx.o m_mpi_utils.o m_ntm.o
setup_hamiltonian.o: m_spin.o m_steps.o m_stress.o m_ts_global_vars.o
setup_hamiltonian.o: metaforce.o molecularmechanics.o normalize_dm.o parallel.o
setup_hamiltonian.o: siesta_geom.o siesta_options.o sparse_matrices.o sys.o
setup_kscell.o: minvec.o parallel.o precision.o units.o
setup_ordern_indexes.o: alloc.o domain_decom.o parallel.o spatial.o
shaper.o: atmfuncs.o mneighb.o precision.o
show_distribution.o: atomlist.o parallel.o parallelsubs.o siesta_geom.o sys.o
siesta.o: alloc.o m_steps.o siesta_analysis.o siesta_end.o siesta_forces.o
siesta.o: siesta_init.o siesta_move.o
siesta_analysis.o: alloc.o atomlist.o bands.o basis_enthalpy.o dhscf.o files.o
siesta_analysis.o: kpoint_grid.o ksv.o ksvinit.o local_DOS.o m_dipol.o
siesta_analysis.o: m_energies.o m_eo.o m_forces.o m_gamma.o m_mpi_utils.o
siesta_analysis.o: m_ntm.o m_spin.o m_steps.o parallel.o projected_DOS.o
siesta_analysis.o: siesta_cml.o siesta_geom.o siesta_options.o
siesta_analysis.o: sparse_matrices.o write_subs.o writewave.o zmatrix.o
siesta_cmlsubs.o: files.o parallel.o siesta_cml.o timestamp.o 
siesta_end.o: alloc.o bands.o densematrix.o fixed.o m_dscfcomm.o m_pulay.o
siesta_end.o: m_wallclock.o meshdscf.o meshphi.o moremeshsubs.o parallel.o
siesta_end.o: siesta_cml.o siesta_cmlsubs.o sparse_matrices.o timestamp.o
siesta_forces.o: compute_dm.o m_convergence.o m_steps.o m_ts_global_vars.o
siesta_forces.o: m_ts_options.o mixer.o parallel.o post_scf_work.o
siesta_forces.o: scfconvergence_test.o setup_H0.o setup_hamiltonian.o
siesta_forces.o: siesta_cml.o siesta_options.o state_analysis.o state_init.o
siesta_forces.o: sys.o
siesta_geom.o: precision.o
siesta_init.o: alloc.o atomlist.o bands.o bsc_xcmod.o densematrix.o files.o
siesta_init.o: ioxv.o kpoint_grid.o kpoint_pdos.o ksvinit.o m_check_supercell.o
siesta_init.o: m_energies.o m_eo.o m_forces.o m_gamma.o m_iostruct.o
siesta_init.o: m_mpi_utils.o m_rmaxh.o m_spin.o m_steps.o m_timer.o
siesta_init.o: m_ts_in_siesta.o m_ts_options.o m_wallclock.o metaforce.o
siesta_init.o: molecularmechanics.o parallel.o parallelsubs.o phonon.o
siesta_init.o: projected_DOS.o siesta_cmlsubs.o siesta_geom.o siesta_options.o
siesta_init.o: sparse_matrices.o struct_init.o sys.o timestamp.o writewave.o
siesta_init.o: zmatrix.o
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
state_init.o: alloc.o atomlist.o domain_decom.o files.o hsparse.o iodm_netcdf.o
state_init.o: iodmhs_netcdf.o kpoint_grid.o m_check_supercell.o m_energies.o
state_init.o: m_eo.o m_gamma.o m_mpi_utils.o m_pulay.o m_rmaxh.o m_spin.o
state_init.o: m_steps.o m_ts_global_vars.o m_ts_io.o m_ts_options.o new_dm.o
state_init.o: overlap.o parallel.o proximity_check.o siesta_cml.o siesta_geom.o
state_init.o: siesta_options.o sparse_matrices.o sys.o write_subs.o zmatrix.o
struct_init.o: alloc.o atmfuncs.o atomlist.o files.o iopipes.o ioxv.o
struct_init.o: m_iostruct.o parallel.o periodic_table.o siesta_cml.o
struct_init.o: siesta_geom.o siesta_options.o zmatrix.o
sys.o: parallel.o
timer.o: m_timer.o sys.o
transition_rate.o: alloc.o fermid.o parallel.o parallelsubs.o precision.o sys.o
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
writewave.o: alloc.o atmfuncs.o atomlist.o densematrix.o files.o kpoint_grid.o
writewave.o: parallel.o parallelsubs.o precision.o siesta_geom.o sys.o units.o
xc.o: alloc.o bsc_xcmod.o precision.o sys.o
xml.o: precision.o
zm_broyden_optim.o: alloc.o m_broyddj_nocomm.o parallel.o precision.o sys.o
zm_broyden_optim.o: units.o zmatrix.o
zm_fire_optim.o: m_fdf_global.o m_fire.o parallel.o precision.o
zm_fire_optim.o: siesta_options.o units.o zmatrix.o
zmatrix.o: alloc.o m_cell.o parallel.o precision.o siesta_geom.o sys.o units.o
band.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/bands.o
diagmemory.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/memoryinfo.o
listsc_module.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/listsc.o
m_bessph.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/bessph.o
m_born_charge.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/born_charge.o
m_broyden_optim.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/broyden_optim.o
m_cell_broyden_optim.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/cell_broyden_optim.o
m_cell_fire_optim.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/cell_fire_optim.o
m_cgwf.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/cgwf.o
m_chempot.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/chempot.o
m_chkgmx.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/chkgmx.o
m_compute_dm.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/compute_dm.o
m_conjgr.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/conjgr.o
m_conjgr_old.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/conjgr_old.o
m_denmat.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/denmat.o
m_dhscf.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/dhscf.o
m_dnaefs.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/dnaefs.o
m_doping_uniform.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/doping_uniform.o
m_dynamics.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/dynamics.o
m_eandg.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/egandd.o
m_ener3.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/ener3.o
m_fermid.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/fermid.o
m_fft.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/fft.o
m_final_h_f_stress.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/final_H_f_stress.o
m_find_kgrid.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/find_kgrid.o
m_fire_optim.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/fire_optim.o
m_fixed.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/fixed.o
m_forhar.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/forhar.o
m_gradient.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/gradient.o
m_grdsam.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/grdsam.o
m_hsparse.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/hsparse.o
m_intramol_pressure.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/remove_intramol_pressure.o
m_ioxv.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/ioxv.o
m_kinefsm.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/kinefsm.o
m_ksv.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/ksv.o
m_ksvinit.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/ksvinit.o
m_local_dos.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/local_DOS.o
m_minvec.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/minvec.o
m_mixer.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/mixer.o
m_naefs.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/naefs.o
m_new_dm.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/new_dm.o
m_nlefsm.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/nlefsm.o
m_norm.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/compute_norm.o
m_normalize_dm.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/normalize_dm.o
m_on_subs.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/on_subs.o
m_ordern.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/ordern.o
m_overfsm.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/overfsm.o
m_overlap.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/overlap.o
m_post_scf_work.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/post_scf_work.o
m_projected_dos.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/projected_DOS.o
m_proximity_check.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/proximity_check.o
m_radfft.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/radfft.o
m_scf_options.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/m_options.o
m_scfconvergence_test.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/scfconvergence_test.o
m_setup_h0.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/setup_H0.o
m_setup_hamiltonian.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/setup_hamiltonian.o
m_siesta_analysis.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/siesta_analysis.o
m_siesta_end.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/siesta_end.o
m_siesta_forces.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/siesta_forces.o
m_siesta_init.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/siesta_init.o
m_siesta_move.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/siesta_move.o
m_state_analysis.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/state_analysis.o
m_state_init.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/state_init.o
m_struct_init.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/struct_init.o
m_timestamp.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/timestamp.o
m_vmat.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/vmat.o
m_xyz_phiatm.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/matel.o
m_zm_broyden_optim.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/zm_broyden_optim.o
m_zm_fire_optim.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/zm_fire_optim.o
neighbour.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/mneighb.o
on_core.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/onmod.o
on_main.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/onmod.o
printmat.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/printmatrix.o
version_info.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/version.o
write_subs_energies.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/write_subs.o
write_subs_positions.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/write_subs.o
write_subs_pressure.o: /Users/ag/S/trunk/Util/Denchar/Src/../../../Src/write_subs.o
