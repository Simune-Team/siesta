#
# Makefile for SIESTA tests
#
# See also test.mk for details
#
tests= h2o h2o_dos h2o_orderN floating bessel si mgco3 \
       var_cell constant_volume batio3 fe fe_broyden sih si64 \
       h2oZ sih_op_broyden h2o_op_broyden
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
clean:
	@for i in ${tests} ; do \
		(cd $$i ; make clean) ; \
	done






