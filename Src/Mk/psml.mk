# PSML helper mk
#
# You need to define PSML_ROOT in the main arch.make file
# to point to your installation of libpsml.
#
PSML_INCFLAGS=-I $(PSML_ROOT)/include
PSML_LIBS=$(PSML_ROOT)/lib/libpsml.a
#
