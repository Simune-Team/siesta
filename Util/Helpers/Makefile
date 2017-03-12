# ---
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---
#
# Makefile for stand-alone get_chem_labels
#
OBJDIR=Obj
#
.SUFFIXES: .f .F .o .a  .f90 .F90
#
VPATH:=$(shell pwd)/../../Src
#
default: get_chem_labels
#
dep:
	sfmakedepend --depend=obj  --modext=o \
          $(VPATH)/*.f $(VPATH)/*.f90 $(VPATH)/*.F $(VPATH)/*.F90 
#
ARCH_MAKE=../../$(OBJDIR)/arch.make
include $(ARCH_MAKE)
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
DEFS:=$(DEFS_PREFIX)-UMPI $(DEFS) $(DEFS_PREFIX)-UMPI
FPPFLAGS:=$(DEFS_PREFIX)-UMPI $(FPPFLAGS) $(DEFS_PREFIX)-UMPI
INCFLAGS:=$(NETCDF_INCFLAGS) $(INCFLAGS)
#
# Uncomment the following line for debugging support
#
#FFLAGS=$(FFLAGS_DEBUG)
#
SYSOBJ=$(SYS).o
#
# Note that machine-specific files are now in top Src directory.
#
OBJS_GET_CHEM_LABELS=precision.o f2kcli.o parallel.o \
          chemical.o io.o m_io.o alloc.o \
          memoryinfo.o memory.o periodic_table.o\
          pxf.o atom_options.o  \
          local_sys.o get_chem_labels.o $(SYSOBJ)
#
#
# Use the makefile in Src/fdf and all the sources there.
#
FDF=libfdf.a
FDF_MAKEFILE=$(VPATH)/fdf/makefile
FDF_INCFLAGS:=-I $(VPATH)/fdf $(INCFLAGS)
$(FDF): 
	(mkdir -p fdf; cd fdf ; $(MAKE) -f $(FDF_MAKEFILE) "FC=$(FC)" "VPATH=$(VPATH)/fdf" \
                          "ARCH_MAKE=../$(ARCH_MAKE)" \
                          "DEFS=$(DEFS)" \
                          "FPPFLAGS=$(FPPFLAGS)" \
                          "MPI_INTERFACE= " \
                          "INCFLAGS=$(FDF_INCFLAGS)" "FFLAGS=$(FFLAGS)" module)
#
get_chem_labels: $(FDF)  $(OBJS_GET_CHEM_LABELS)
	$(FC) -o get_chem_labels \
	       $(LDFLAGS) $(OBJS_GET_CHEM_LABELS) $(FDF) $(LIBS)
#
clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f get_chem_labels ioncat *.o  *.a
	rm -f *.mod .tmp_fdflog
	rm -rf fdf

# DO NOT DELETE THIS LINE - used by make depend
alloc.o: parallel.o precision.o local_sys.o
arw.o: parallel.o local_sys.o
atm_transfer.o: atm_types.o atmparams.o old_atmfuncs.o periodic_table.o
atm_transfer.o: radial.o local_sys.o
atm_types.o: precision.o radial.o
atmfuncs.o: atm_types.o precision.o radial.o spher_harm.o local_sys.o
atom.o: atmparams.o atom_options.o basis_specs.o basis_types.o m_filter.o
atom.o: old_atmfuncs.o periodic_table.o precision.o pseudopotential.o recipes.o
atom.o: local_sys.o xcmod.o
atom_options.o: local_sys.o
atomlist.o: alloc.o atm_types.o atmfuncs.o parallel.o precision.o siesta_geom.o
atomlwf.o: alloc.o atmfuncs.o mneighb.o onmod.o precision.o spatial.o local_sys.o
automatic_cell.o: atmfuncs.o parallel.o precision.o units.o
bands.o: alloc.o densematrix.o files.o parallel.o parallelsubs.o parsing.o
bands.o: precision.o local_sys.o
basis_enthalpy.o: atmfuncs.o atomlist.o precision.o siesta_geom.o units.o
basis_io.o: atm_types.o atmparams.o atom_options.o basis_types.o chemical.o
basis_io.o: precision.o pseudopotential.o radial.o local_sys.o xml.o
basis_specs.o: basis_types.o chemical.o periodic_table.o precision.o
basis_specs.o: pseudopotential.o local_sys.o
basis_types.o: alloc.o atmparams.o precision.o pseudopotential.o local_sys.o
bessph.o: parallel.o precision.o local_sys.o
bonds.o: alloc.o atmfuncs.o mneighb.o precision.o recipes.o sorting.o units.o
born_charge.o: atomlist.o ksv.o ksvinit.o m_mpi_utils.o m_spin.o parallel.o
born_charge.o: siesta_geom.o siesta_options.o sparse_matrices.o
broadcast_basis.o: atm_types.o parallel.o radial.o
broyden_optim.o: m_broyddj_nocomm.o m_fdf_global.o m_memory.o m_mpi_utils.o
broyden_optim.o: parallel.o precision.o units.o
cart2frac.o: local_sys.o
cdiag.o: alloc.o m_spin.o memoryinfo.o parallel.o precision.o local_sys.o
cell_broyden_optim.o: m_broyddj_nocomm.o m_fdf_global.o m_memory.o
cell_broyden_optim.o: m_mpi_utils.o parallel.o precision.o local_sys.o units.o
cell_broyden_optim.o: zmatrix.o
cell_fire_optim.o: alloc.o m_fdf_global.o m_fire.o m_memory.o parallel.o
cell_fire_optim.o: precision.o siesta_options.o local_sys.o zmatrix.o
cellxc.o: alloc.o meshsubs.o parallel.o parallelsubs.o precision.o local_sys.o
cellxc.o: xcmod.o
cgvc.o: conjgr_old.o m_fdf_global.o m_mpi_utils.o parallel.o precision.o
cgvc.o: units.o
cgvc_zmatrix.o: conjgr.o m_fdf_global.o m_mpi_utils.o parallel.o precision.o
cgvc_zmatrix.o: units.o zmatrix.o
cgwf.o: onmod.o precision.o local_sys.o
chemical.o: precision.o local_sys.o
chempot.o: alloc.o mneighb.o parallel.o parallelsubs.o precision.o local_sys.o
chkdim.o: local_sys.o
chkgmx.o: minvec.o parallel.o precision.o local_sys.o
coceri.o: files.o periodic_table.o precision.o
compute_dm.o: atomlist.o iodmhs_netcdf.o kpoint_grid.o m_energies.o m_eo.o
compute_dm.o: m_gamma.o m_mpi_utils.o m_rmaxh.o m_spin.o m_steps.o
compute_dm.o: m_ts_in_siesta.o m_ts_options.o ordern.o parallel.o precision.o
compute_dm.o: siesta_geom.o siesta_options.o sparse_matrices.o local_sys.o
compute_norm.o: m_mpi_utils.o precision.o sparse_matrices.o
conjgr.o: precision.o
conjgr_old.o: precision.o
coor.o: alloc.o m_mpi_utils.o parallel.o precision.o siesta_geom.o local_sys.o
coor.o: zmatrix.o
coxmol.o: files.o periodic_table.o precision.o
denmat.o: alloc.o globalise.o m_mpi_utils.o onmod.o precision.o
denmatlomem.o: alloc.o globalise.o onmod.o precision.o
densematrix.o: precision.o
detover.o: alloc.o parallel.o parallelsubs.o precision.o
dfscf.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o meshmatrix.o
dfscf.o: meshsubs.o meshsubs.o parallel.o precision.o local_sys.o
dhscf.o: alloc.o atmfuncs.o files.o forhar.o iogrid_netcdf.o m_iorho.o
dhscf.o: m_ts_in_siesta.o m_ts_options.o meshsubs.o parallel.o parsing.o
dhscf.o: precision.o siesta_options.o local_sys.o xcmod.o
diag2g.o: fermid.o m_spin.o parallel.o parallelsubs.o precision.o local_sys.o
diag2k.o: fermid.o m_spin.o parallel.o parallelsubs.o precision.o local_sys.o
diagg.o: fermid.o parallel.o parallelsubs.o precision.o local_sys.o writewave.o
diagk.o: compute_norm.o fermid.o parallel.o parallelsubs.o precision.o local_sys.o
diagk.o: writewave.o
diagk_file.o: fermid.o iowfs_netcdf.o parallel.o parallelsubs.o precision.o
diagk_file.o: local_sys.o writewave.o
diagkp.o: fermid.o parallel.o parallelsubs.o precision.o local_sys.o
diagon.o: alloc.o densematrix.o m_memory.o parallel.o parallelsubs.o
diagon.o: precision.o siesta_options.o local_sys.o
diagpol.o: precision.o local_sys.o
diagsprl.o: atomlist.o fermid.o parallel.o parallelsubs.o precision.o
diagsprl.o: siesta_geom.o local_sys.o
dipole.o: parallel.o precision.o local_sys.o
dnaefs.o: atmfuncs.o mneighb.o precision.o
doping.o: meshsubs.o meshsubs.o parallel.o precision.o local_sys.o
dynamics.o: atomlist.o files.o ioxv.o m_mpi_utils.o parallel.o precision.o
dynamics.o: local_sys.o units.o
efield.o: atmfuncs.o parallel.o parallelsubs.o parsing.o precision.o
efield.o: siesta_cml.o local_sys.o
egandd.o: alloc.o globalise.o m_mpi_utils.o onmod.o onmod.o precision.o local_sys.o
eggbox.o: parallel.o parsing.o precision.o
electrostatic.o: atm_types.o atmfuncs.o atmparams.o precision.o radfft.o
electrostatic.o: radial.o recipes.o local_sys.o
ener3.o: globalise.o m_mpi_utils.o onmod.o precision.o
ener3lomem.o: globalise.o m_mpi_utils.o onmod.o precision.o
extrapolon.o: parallel.o precision.o local_sys.o
fermid.o: parallel.o precision.o recipes.o local_sys.o
fft3d.o: alloc.o parallel.o parallelsubs.o precision.o local_sys.o
find_kgrid.o: alloc.o minvec.o parallel.o precision.o units.o
fire_optim.o: alloc.o m_fdf_global.o m_fire.o m_mpi_utils.o parallel.o
fire_optim.o: precision.o siesta_options.o units.o
fixed.o: alloc.o parallel.o parsing.o precision.o
flib_spline.o: precision.o
forhar.o: alloc.o meshsubs.o precision.o
get_target_stress.o: m_mpi_utils.o parallel.o precision.o units.o
globalise.o: alloc.o onmod.o onmod.o precision.o spatial.o
gradient.o: globalise.o m_mpi_utils.o onmod.o precision.o
gradientlomem.o: globalise.o onmod.o precision.o
grdsam.o: parallel.o parsing.o precision.o local_sys.o
hsparse.o: alloc.o atm_types.o atmfuncs.o listsc.o mneighb.o parallel.o
hsparse.o: parallelsubs.o precision.o sorting.o local_sys.o
idiag.o: parallel.o local_sys.o
initatom.o: atom.o atom_options.o basis_io.o basis_specs.o basis_types.o
initatom.o: electrostatic.o old_atmfuncs.o precision.o
initparallel.o: alloc.o parallel.o parallelsubs.o precision.o spatial.o local_sys.o
io.o: local_sys.o
iocg.o: files.o precision.o local_sys.o
iodm.o: files.o parallel.o parallelsubs.o precision.o local_sys.o
iodm_netcdf.o: alloc.o parallel.o parallelsubs.o precision.o local_sys.o
iodmhs_netcdf.o: alloc.o parallel.o parallelsubs.o precision.o local_sys.o
ioeig.o: files.o precision.o siesta_cml.o units.o
iofa.o: files.o precision.o
iogrid_netcdf.o: alloc.o parallel.o precision.o local_sys.o
iohs.o: atm_types.o atmfuncs.o atomlist.o files.o parallel.o parallelsubs.o
iohs.o: precision.o siesta_geom.o local_sys.o
iokp.o: files.o precision.o
iolwf.o: alloc.o files.o onmod.o parallel.o parallelsubs.o precision.o local_sys.o
iomd.o: files.o precision.o
iopipes.o: parallel.o precision.o local_sys.o
iowfs_netcdf.o: alloc.o parallel.o parallelsubs.o precision.o local_sys.o
ioxv.o: files.o parallel.o precision.o
iozm.o: files.o parallel.o precision.o siesta_geom.o zmatrix.o
ipack.o: local_sys.o
kgrid.o: parallel.o precision.o units.o
kgridinit.o: minvec.o parallel.o precision.o units.o
kinefsm.o: alloc.o atmfuncs.o mneighb.o parallel.o parallelsubs.o precision.o
kpoint_grid.o: find_kgrid.o minvec.o parallel.o precision.o siesta_cml.o
kpoint_grid.o: siesta_options.o units.o
kpoint_pdos.o: find_kgrid.o minvec.o parallel.o precision.o siesta_options.o
kpoint_pdos.o: units.o
ksv.o: alloc.o atmfuncs.o densematrix.o ksvinit.o parallel.o parsing.o
ksv.o: precision.o local_sys.o
ksvinit.o: alloc.o parallel.o parsing.o precision.o
listsc.o: alloc.o
local_DOS.o: atomlist.o files.o kpoint_grid.o m_energies.o m_eo.o m_forces.o
local_DOS.o: m_gamma.o m_mpi_utils.o m_ntm.o m_spin.o parallel.o siesta_geom.o
local_DOS.o: siesta_options.o sparse_matrices.o local_sys.o
m_broyddj.o: alloc.o m_mpi_utils.o parallel.o precision.o local_sys.o
m_broyddj_nocomm.o: alloc.o m_fdf_global.o parallel.o precision.o local_sys.o
m_broyden_mixing.o: alloc.o m_broyddj.o m_mpi_utils.o parallel.o precision.o
m_cell.o: precision.o siesta_geom.o units.o
m_check_supercell.o: minvec.o parallel.o precision.o
m_dipol.o: precision.o
m_energies.o: precision.o
m_eo.o: precision.o
m_fdf_global.o: m_mpi_utils.o parallel.o precision.o
m_fire.o: parallel.o precision.o
m_fire_mixing.o: alloc.o m_fdf_global.o m_fire_para.o m_mpi_utils.o parallel.o
m_fire_mixing.o: precision.o
m_fire_para.o: m_mpi_utils.o parallel.o precision.o
m_forces.o: precision.o
m_getopts.o: f2kcli.o
m_hsx.o: atm_types.o atmfuncs.o atomlist.o files.o parallel.o parallelsubs.o
m_hsx.o: precision.o siesta_geom.o local_sys.o
m_iodm.o: alloc.o files.o m_mpi_utils.o parallel.o parallelsubs.o precision.o
m_iodm.o: local_sys.o
m_iorho.o: parallel.o parallelsubs.o precision.o local_sys.o
m_iostruct.o: alloc.o files.o m_mpi_utils.o parallel.o precision.o
m_iostruct.o: siesta_geom.o local_sys.o units.o
m_kinetic.o: precision.o
m_memory.o: memoryinfo.o parallel.o local_sys.o
m_mpi_utils.o: precision.o local_sys.o
m_options.o: precision.o
m_pulay.o: alloc.o atomlist.o files.o m_mpi_utils.o m_spin.o parallel.o
m_pulay.o: precision.o siesta_options.o sparse_matrices.o local_sys.o
m_rmaxh.o: precision.o
m_sparse.o: alloc.o precision.o
m_spin.o: alloc.o m_fdf_global.o parallel.o precision.o
m_stress.o: precision.o
m_ts_aux_rout.o: local_sys.o
m_ts_contour.o: files.o m_ts_aux_rout.o m_ts_options.o precision.o local_sys.o
m_ts_electrode.o: files.o m_ts_aux_rout.o m_ts_io.o m_ts_kpoints.o
m_ts_electrode.o: m_ts_options.o
m_ts_in_siesta.o: find_kgrid.o m_ts_contour.o m_ts_electrode.o m_ts_io.o
m_ts_in_siesta.o: m_ts_kpoints.o m_ts_options.o m_ts_scattering.o parallel.o
m_ts_in_siesta.o: parallelsubs.o precision.o local_sys.o
m_ts_io.o: files.o m_spin.o m_ts_kpoints.o parallel.o parallelsubs.o
m_ts_io.o: precision.o local_sys.o
m_ts_kpoints.o: minvec.o parallel.o precision.o siesta_cml.o siesta_options.o
m_ts_kpoints.o: local_sys.o units.o
m_ts_options.o: m_fdf_global.o parallel.o precision.o siesta_options.o local_sys.o
m_ts_scattering.o: m_ts_aux_rout.o m_ts_options.o parallelsubs.o precision.o
m_ts_scattering.o: local_sys.o
m_wallclock.o: m_walltime.o
madelung.o: parallel.o precision.o local_sys.o
matel.o: alloc.o atmfuncs.o precision.o radfft.o recipes.o spher_harm.o
md_out.o: files.o m_energies.o precision.o local_sys.o units.o
md_utils.o: precision.o
memory.o: memoryinfo.o parallel.o precision.o
memoryinfo.o: precision.o
meshmatrix.o: alloc.o atomlist.o meshsubs.o parallel.o parallelsubs.o
meshmatrix.o: precision.o
meshsubs.o: atm_types.o atmfuncs.o parallel.o parallelsubs.o precision.o
meshsubs.o: radial.o siesta_cml.o
metaforce.o: parallel.o parsing.o precision.o
minvec.o: precision.o sorting.o local_sys.o
mixer.o: atomlist.o iodm_netcdf.o m_broyden_mixing.o m_energies.o
mixer.o: m_fire_mixing.o m_iodm.o m_mpi_utils.o m_pulay.o m_spin.o m_steps.o
mixer.o: m_ts_io.o m_ts_options.o parallel.o siesta_geom.o siesta_options.o
mixer.o: sparse_matrices.o
mneighb.o: alloc.o precision.o local_sys.o
molecularmechanics.o: alloc.o parallel.o parsing.o precision.o local_sys.o units.o
mulliken.o: atmfuncs.o parallel.o parallelsubs.o precision.o siesta_cml.o
naefs.o: atmfuncs.o mneighb.o precision.o
new_dm.o: alloc.o atomlist.o compute_norm.o m_energies.o m_iodm.o m_sparse.o
new_dm.o: m_spin.o m_steps.o m_ts_io.o m_ts_options.o parallel.o parallelsubs.o
new_dm.o: parsing.o precision.o siesta_geom.o siesta_options.o
new_dm.o: sparse_matrices.o local_sys.o units.o
nlefsm.o: alloc.o atmfuncs.o mneighb.o parallel.o parallelsubs.o precision.o
normalize_dm.o: atomlist.o m_mpi_utils.o m_spin.o parallel.o precision.o
normalize_dm.o: sparse_matrices.o
obc.o: alloc.o precision.o
ofc.o: alloc.o files.o precision.o units.o
old_atmfuncs.o: alloc.o atmparams.o precision.o local_sys.o
on_subs.o: alloc.o onmod.o onmod.o spatial.o
onmod.o: precision.o
optical.o: alloc.o atmfuncs.o densematrix.o fermid.o files.o parallel.o
optical.o: parallelsubs.o parsing.o precision.o local_sys.o units.o
ordern.o: alloc.o chempot.o globalise.o m_mpi_utils.o onmod.o parallel.o
ordern.o: precision.o siesta_cml.o local_sys.o units.o
outcell.o: precision.o siesta_cml.o units.o
outcoor.o: alloc.o atmfuncs.o precision.o siesta_geom.o local_sys.o units.o
overfsm.o: alloc.o atmfuncs.o mneighb.o parallel.o parallelsubs.o precision.o
overlap.o: alloc.o atmfuncs.o mneighb.o parallel.o parallelsubs.o precision.o
parallelsubs.o: parallel.o spatial.o local_sys.o
pdos.o: alloc.o atmfuncs.o atomlist.o densematrix.o files.o parallel.o
pdos.o: parallelsubs.o precision.o siesta_geom.o units.o xml.o
pdosg.o: alloc.o parallel.o parallelsubs.o precision.o local_sys.o units.o
pdosk.o: alloc.o parallel.o parallelsubs.o precision.o local_sys.o units.o
pdoskp.o: alloc.o parallel.o parallelsubs.o precision.o local_sys.o units.o
periodic_table.o: precision.o local_sys.o
phirphi.o: alloc.o atmfuncs.o mneighb.o parallel.o parallelsubs.o precision.o
phirphi_opt.o: alloc.o atm_types.o atmfuncs.o atmparams.o mneighb.o parallel.o
phirphi_opt.o: parallelsubs.o precision.o spher_harm.o local_sys.o
phonon.o: files.o parallel.o precision.o local_sys.o
pixmol.o: files.o periodic_table.o precision.o
plcharge.o: files.o old_atmfuncs.o parallel.o precision.o
poison.o: alloc.o parallel.o precision.o local_sys.o
post_scf_work.o: atomlist.o m_energies.o m_spin.o mneighb.o setup_hamiltonian.o
post_scf_work.o: siesta_cml.o siesta_geom.o siesta_options.o sparse_matrices.o
projected_DOS.o: alloc.o atomlist.o kpoint_grid.o kpoint_pdos.o m_eo.o
projected_DOS.o: m_gamma.o m_mpi_utils.o m_spin.o parallel.o precision.o
projected_DOS.o: siesta_geom.o siesta_options.o sparse_matrices.o local_sys.o
propor.o: precision.o local_sys.o
proximity_check.o: m_ts_options.o mneighb.o parallel.o precision.o
proximity_check.o: siesta_geom.o siesta_options.o
pseudopotential.o: alloc.o atom_options.o flib_spline.o precision.o local_sys.o
radfft.o: alloc.o precision.o recipes.o
radial.o: precision.o xml.o
rdiag.o: alloc.o memoryinfo.o parallel.o precision.o local_sys.o
readsp.o: parallel.o precision.o local_sys.o
recipes.o: local_sys.o
redcel.o: precision.o local_sys.o
reinit.o: files.o parallel.o precision.o
remove_intramol_pressure.o: precision.o zmatrix.o
reoptical.o: parallel.o parsing.o precision.o local_sys.o
reord.o: alloc.o precision.o
rhooda.o: atmfuncs.o atomlist.o meshsubs.o meshsubs.o precision.o
rhoofd.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o meshmatrix.o
rhoofd.o: meshsubs.o meshsubs.o parallel.o precision.o local_sys.o
rhoofdsp.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o meshmatrix.o
rhoofdsp.o: meshsubs.o meshsubs.o precision.o local_sys.o
savepsi.o: alloc.o parallel.o parallelsubs.o precision.o
scfconvergence_test.o: m_convergence.o m_energies.o m_ts_options.o
scfconvergence_test.o: m_wallclock.o parallel.o siesta_cml.o siesta_options.o
scfconvergence_test.o: write_subs.o
setatomnodes.o: alloc.o parallel.o precision.o spatial.o local_sys.o
setspatial.o: alloc.o parallel.o precision.o spatial.o
setup_hamiltonian.o: alloc.o atmfuncs.o atomlist.o files.o m_dipol.o
setup_hamiltonian.o: m_energies.o m_forces.o m_gamma.o m_hsx.o m_mpi_utils.o
setup_hamiltonian.o: m_ntm.o m_spin.o m_steps.o m_stress.o m_ts_io.o
setup_hamiltonian.o: m_ts_kpoints.o m_ts_options.o metaforce.o
setup_hamiltonian.o: molecularmechanics.o nlefsm.o normalize_dm.o overfsm.o
setup_hamiltonian.o: siesta_cml.o siesta_geom.o siesta_options.o
setup_hamiltonian.o: sparse_matrices.o local_sys.o
setup_kscell.o: minvec.o parallel.o precision.o units.o
shaper.o: alloc.o atmfuncs.o mneighb.o precision.o
show_distribution.o: atomlist.o parallel.o parallelsubs.o siesta_geom.o local_sys.o
siesta.o: m_steps.o siesta_analysis.o siesta_end.o siesta_forces.o
siesta.o: siesta_init.o siesta_move.o
siesta_analysis.o: alloc.o atomlist.o bands.o basis_enthalpy.o files.o
siesta_analysis.o: kpoint_grid.o ksv.o ksvinit.o local_DOS.o m_dipol.o
siesta_analysis.o: m_energies.o m_eo.o m_forces.o m_gamma.o m_mpi_utils.o
siesta_analysis.o: m_ntm.o m_spin.o m_steps.o parallel.o projected_DOS.o
siesta_analysis.o: siesta_cml.o siesta_geom.o siesta_options.o
siesta_analysis.o: sparse_matrices.o write_subs.o writewave.o zmatrix.o
siesta_cml.o: siesta_cmlsubs.o
siesta_cmlsubs.o: files.o parallel.o timestamp.o version.o
siesta_end.o: alloc.o m_wallclock.o parallel.o siesta_cml.o timestamp.o
siesta_forces.o: compute_dm.o m_convergence.o m_spin.o m_steps.o m_ts_options.o
siesta_forces.o: mixer.o parallel.o post_scf_work.o scfconvergence_test.o
siesta_forces.o: setup_hamiltonian.o siesta_cml.o siesta_options.o
siesta_forces.o: sparse_matrices.o state_analysis.o state_init.o local_sys.o
siesta_geom.o: precision.o
siesta_init.o: alloc.o atomlist.o bands.o densematrix.o files.o ioxv.o
siesta_init.o: kpoint_grid.o kpoint_pdos.o ksvinit.o m_check_supercell.o
siesta_init.o: m_energies.o m_eo.o m_fdf_global.o m_forces.o m_gamma.o
siesta_init.o: m_iostruct.o m_mpi_utils.o m_rmaxh.o m_spin.o m_steps.o
siesta_init.o: m_ts_in_siesta.o m_ts_options.o m_wallclock.o metaforce.o
siesta_init.o: molecularmechanics.o parallel.o parallelsubs.o phonon.o
siesta_init.o: projected_DOS.o siesta_cmlsubs.o siesta_geom.o siesta_options.o
siesta_init.o: sparse_matrices.o struct_init.o local_sys.o timestamp.o writewave.o
siesta_init.o: xcmod.o zmatrix.o
siesta_move.o: atomlist.o broyden_optim.o cell_broyden_optim.o
siesta_move.o: cell_fire_optim.o dynamics.o fire_optim.o iopipes.o ioxv.o
siesta_move.o: m_energies.o m_forces.o m_kinetic.o m_steps.o m_stress.o
siesta_move.o: parallel.o phonon.o siesta_cml.o siesta_geom.o siesta_options.o
siesta_move.o: local_sys.o write_subs.o zm_broyden_optim.o zm_fire_optim.o zmatrix.o
siesta_options.o: files.o m_fdf_global.o m_mpi_utils.o memoryinfo.o parallel.o
siesta_options.o: precision.o siesta_cml.o local_sys.o units.o
sparse_matrices.o: precision.o
spatial.o: precision.o
spher_harm.o: alloc.o precision.o local_sys.o
state_analysis.o: atomlist.o born_charge.o m_energies.o m_forces.o m_iostruct.o
state_analysis.o: m_ntm.o m_spin.o m_stress.o m_wallclock.o parallel.o
state_analysis.o: remove_intramol_pressure.o siesta_cml.o siesta_geom.o
state_analysis.o: siesta_options.o sparse_matrices.o write_subs.o zmatrix.o
state_init.o: alloc.o atomlist.o files.o hsparse.o iodm_netcdf.o
state_init.o: iodmhs_netcdf.o kpoint_grid.o m_check_supercell.o m_energies.o
state_init.o: m_eo.o m_gamma.o m_mpi_utils.o m_pulay.o m_rmaxh.o m_spin.o
state_init.o: m_steps.o m_ts_io.o m_ts_options.o new_dm.o overlap.o parallel.o
state_init.o: proximity_check.o siesta_cml.o siesta_cmlsubs.o siesta_geom.o
state_init.o: siesta_options.o sparse_matrices.o local_sys.o write_subs.o zmatrix.o
struct_init.o: alloc.o atmfuncs.o atomlist.o files.o iopipes.o ioxv.o
struct_init.o: m_fdf_global.o m_iostruct.o parallel.o periodic_table.o
struct_init.o: siesta_cml.o siesta_geom.o siesta_options.o zmatrix.o
timer.o: m_walltime.o parallel.o precision.o local_sys.o
transition_rate.o: alloc.o fermid.o parallel.o parallelsubs.o precision.o local_sys.o
typecell.o: precision.o
uncell.o: precision.o
units.o: precision.o
vmat.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o meshmatrix.o
vmat.o: meshsubs.o meshsubs.o parallel.o precision.o
vmatsp.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o meshmatrix.o
vmatsp.o: meshsubs.o meshsubs.o precision.o
vmb.o: parallel.o precision.o local_sys.o
write_md_record.o: atomlist.o m_energies.o m_kinetic.o m_steps.o m_stress.o
write_md_record.o: md_out.o parallel.o siesta_geom.o siesta_options.o
write_raw_efs.o: atmfuncs.o atomlist.o precision.o siesta_geom.o
write_subs.o: atomlist.o m_energies.o m_forces.o m_iostruct.o m_spin.o
write_subs.o: m_steps.o m_stress.o m_ts_options.o parallel.o precision.o
write_subs.o: siesta_cml.o siesta_geom.o siesta_options.o units.o zmatrix.o
writewave.o: alloc.o atmfuncs.o atomlist.o densematrix.o files.o kpoint_grid.o
writewave.o: m_fdf_global.o parallel.o parallelsubs.o parsing.o precision.o
writewave.o: siesta_geom.o local_sys.o units.o
xc.o: alloc.o precision.o local_sys.o xcmod.o
xcmod.o: parallel.o parsing.o precision.o local_sys.o
xml.o: precision.o
zm_broyden_optim.o: m_broyddj_nocomm.o m_fdf_global.o m_mpi_utils.o parallel.o
zm_broyden_optim.o: precision.o units.o zmatrix.o
zm_fire_optim.o: m_fdf_global.o m_fire.o parallel.o precision.o
zm_fire_optim.o: siesta_options.o units.o zmatrix.o
zmatrix.o: alloc.o m_cell.o parallel.o parsing.o precision.o siesta_geom.o
zmatrix.o: local_sys.o units.o
