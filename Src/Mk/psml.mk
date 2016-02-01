# PSML helper mk
#
PSML_ROOT=$(HOME)/lib/psml/0.9/gfortran-5.2.0
PSML_INCFLAGS=$(PSML_ROOT)/include
INCFLAGS:= $(INCFLAGS) $(PSML_INCFLAGS)
#
PSML_LIBS=$(PSML_ROOT)/lib/libpsml.a
#
