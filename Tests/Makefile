#
# Makefile for SIESTA tests
#
# Usage:
#
#          make                      --- Run all tests
#          make tests=fe_broyden     --- Run a single test
#          make batch                --- Dispatch to batch system
#                                        (site-dependent)
#                                        (needs extra configuration)
#
#          By default, the version of SIESTA used is the one in Src.
#          To use another one:
#
#          make SIESTA=/some/path/to/siesta-program ....
#
#
# See also test.mk for details, and the scripts in Scripts for 
# examples of other options for batch processing.
#
#
tests= h2o h2o_reparam h2o_basis h2o_dos h2o_orderN \
       floating bessel mgco3 si2x1h \
       var_cell constant_volume batio3 fe fe_broyden sih si64 \
       h2oZ sih_op_broyden h2o_op_broyden zmatrix md_pr md_npr \
       md_anneal md_verlet md_nose si_bandpoints sih_fire \
       graphite_c6 oxyn partial h2o_findp_bug
#
all: 
	@for i in ${tests} ; do \
		(cd $$i ; make completed) ; \
	done
#
xmlcheck: all
	@for i in ${tests} ; do \
		(cd $$i ; make xmlcheck) ; \
	done
#
batch:
	@for i in ${tests} ; do \
		(qsub $(HOME)/bin/sge_run.sh "make tests=$$i") ; \
	done
clean:
	@for i in ${tests} ; do \
		(cd $$i ; make clean) ; \
	done






