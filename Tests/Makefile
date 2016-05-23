#
# Makefile for SIESTA tests
#
# Usage:
#
#          make                      --- Run all tests
#          make  check               --- Run all tests and check with reference file
#          make  list                --- Show the list of tests
#
#          make tests=fe_broyden { check | }     --- Run a single test
#          (Not implemented) make batch                --- Dispatch to batch system
#                                        (site-dependent)
#                                        (needs extra configuration)
#
#          By default, the version of SIESTA used is the one in the
#          parent directory (Obj or custom OBJDIR), and is executed
#          in serial mode.
#
#          To use another one, or a parallel version
#
#          make SIESTA=/some/path/to/siesta-program ....
#
#           or
#
#          make SIESTA="mpiexec -np 4 /path/to/parallel-siesta-program" ....
#
#          It is also possible to have separate working directories,
#          by using the a "label". For example:
#
#          make label=finer 
#     
#          The default label is "work".
#
# See also test.mk for details, and the scripts in Scripts for 
# examples of other options for batch processing.

label=work

tests = h2o h2o_reparam h2o_basis h2o_dos h2o-pop \
	floating benzene bessel mgc-force mgco3 si2x1h force_2 born born_spin \
	var_cell constant_volume batio3 fe fe_broyden sih si64 \
	h2oZ sih_op_broyden h2o_op_broyden zmatrix md_pr md_npr \
	md_anneal md_verlet md_nose si_bandpoints sih_fire \
	graphite_c6 oxyn h2o_findp_bug h2o_radialgrid ag ch4 \
	si_coop pb_filter_cutoff pb_filter_tol fire_benzene \
	fe_clust_noncollinear fe_clust_noncollinear-gga \
	h_chain n_chain si64_coop si_pdos_kgrid \
	si2x1h-dipole ar2_vdw graphite_vdw_df sic-slab si-optical OMM_si \
	OMM_h2o \
	SOC_Pt2_xx SOC_Pt2_xz SOC_Pt2_zy SOC_Pt2_zz \
	wannier

# These tests are extremely time consuming and
# should only be runned sometimes
# Currently they may be runned individually.
#	SOC_FePt_xx SOC_FePt_xz SOC_FePt_zy SOC_FePt_zz


default: all

all:
	@for dir in $(tests); do\
	 $(MAKE) -C $$dir completed; \
	done

# Ensure that tests may be called individually
# (and that they are run, PHONY is required)
.PHONY: $(tests)
$(tests):
	$(MAKE) -C $@ completed

check:
	@for dir in $(tests); do\
	 $(MAKE) -C $$dir check; \
	done

show:
	@echo $(tests)

clean:
	@for i in $(tests) ; do \
		(cd $$i ; make clean) ; \
	done
