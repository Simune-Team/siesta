#
# Makefile for PDOS XML processing
#
.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
default: pdos
include ../../Src/arch.make
#
FFLAGS=$(FFLAGS_DEBUG) -I../../Src
#
LIBS= ../../Src/libxmlparser.a
#
OBJS=  m_orbital_chooser.o m_pdos.o f2kcli.o pdos.o
     
pdos: $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS)  $(LIBS)
#
clean: 
	rm -f *.o pdos *.o *.*d
#










