#
# Makefile for PDOS XML processing
#
default: pdos
#
#---------------------------
MK=$(FLIB_ROOT)/fortran.mk
include $(MK)
#---------------------------
#
# Uncomment the following line for debugging support
#
FFLAGS=$(FFLAGS_DEBUG)
#
LIBS=$(LIB_PREFIX)$(LIB_STD) -lflib
#
OBJS=  m_orbital_chooser.o m_pdos.o f2kcli.o pdos.o
     
pdos: $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS)  $(LIBS)
#
clean: 
	rm -f *.o pdos *$(MOD_EXT)
#










