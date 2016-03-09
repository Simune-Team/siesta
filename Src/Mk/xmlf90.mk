# XMLF90 helper mk
#
# You need to define XMLF90_ROOT in the main arch.make file
# to point to your installation of the xmlf90 library.
#
XMLF90_INCFLAGS= -I $(XMLF90_ROOT)/modules
XMLF90_LIBS=$(XMLF90_ROOT)/lib/libxmlf90.a
#
