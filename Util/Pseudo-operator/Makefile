# Makefile for macroave
#
.SUFFIXES: .f .F .o .a  .f90 .F90
#
default: pseudooperator
#
OBJDIR=Obj
VPATH:=$(shell pwd)/../../../Src
include ../../$(OBJDIR)/arch.make
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
SYSOBJ=pxf.o
#
#
# Note that machine-specific files are now in top Src directory.
#
OBJS = precision.o pseudopotential.o flib_spline.o \
	atmparams.o atom_options.o sys.o m_io.o io.o \
	pseudooperator.o module_pseudo.o schrodinger.o arw.o xc.o 
#
#
# This is crude but will have to do for now.
# Note : precision must be the first module
#
MOD_OBJS=
COM_OBJS=$(OBJS) 
ALL_OBJS=$(MOD_OBJS) $(COM_OBJS) $(SYSOBJ)
#
$(COM_OBJS): 
#
pseudooperator:  $(ALL_OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o pseudooperator $(ALL_OBJS) 
#
arw.o: sys.o precision.o 
atom_options.o: sys.o
flib_spline.o: precision.o
io.o: sys.o
m_io.o: sys.o
module_pseudo.o: atmparams.o precision.o schrodinger.o sys.o
pseudopotential.o: atom_options.o flib_spline.o precision.o sys.o m_io.o
pseudooperator.o: precision.o atmparams.o pseudopotential.o module_pseudo.o
schrodinger.o: arw.o precision.o
#
clean: 
	@echo "==> Cleaning object, library, and executable files"
	rm -f macroave permute *.o  *.a *.pcl *.pc *.mod
	rm -f *.mod
#



