# Makefile for kbgen
#
.SUFFIXES: .f .F .o .a  .f90 .F90
#
.PHONY: ncps xmlparser clean
#
default: pseudooperator
#
OBJDIR=Obj
VPATH:=$(shell pwd)/../../Src
include ../../$(OBJDIR)/arch.make
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
INCFLAGS:=$(INCFLAGS) -I../../$(OBJDIR)/nc-pseudo # To pick up mod files
#
LIBS= ../../$(OBJDIR)/libncps.a ../../$(OBJDIR)/libxmlparser.a
#
# Note that machine-specific files are now in top Src directory.
#
OBJS = precision.o  \
	atmparams.o atom_options.o \
	pseudooperator.o module_pseudo.o schrodinger.o arw.o xc.o\
        local_die.o
#
#
MOD_OBJS=
COM_OBJS=$(OBJS) 
ALL_OBJS=$(MOD_OBJS) $(COM_OBJS)
#
$(COM_OBJS): 
#
what:
	echo "VPATH: " $(VPATH)
pseudooperator:  ncps xmlparser $(ALL_OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o pseudooperator $(ALL_OBJS) $(LIBS)
#
arw.o:  precision.o 
module_pseudo.o: atmparams.o precision.o schrodinger.o 
pseudooperator.o: precision.o atmparams.o module_pseudo.o
schrodinger.o: arw.o precision.o
#
clean: 
	@echo "==> Cleaning object, library, and executable files"
	rm -f pseudooperator *.o  *.a *.pcl *.pc *.mod
	rm -f *.mod
#
ncps:
	@echo "Making sure that the ncps library is compiled..."
	(cd ../../$(OBJDIR)/nc-pseudo ; make "VPATH=$(VPATH)/nc-pseudo")
#
xmlparser:
	@echo "Making sure that the xmlparser library is compiled..."
	(cd ../../$(OBJDIR)/xmlparser ; make "VPATH=$(VPATH)/xmlparser")
#


