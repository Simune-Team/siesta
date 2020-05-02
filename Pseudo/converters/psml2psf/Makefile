# Makefile for kbgen
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#
.PHONY: clean
#
default: psml2psf
#
OBJDIR=Obj
TOP_LEVEL=$(shell cd -P -- ../../../Src && pwd -P)
VPATH:=$(TOP_LEVEL)

MAIN_OBJDIR:=$(shell cd -P -- ../../../\$(OBJDIR) && pwd -P)
ARCH_MAKE=$(MAIN_OBJDIR)/arch.make
include $(ARCH_MAKE)
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
NCPS_INCFLAGS=-I$(MAIN_OBJDIR)/ncps/src

NCPS=$(MAIN_OBJDIR)/ncps/src/libncps.a

INCFLAGS:=$(INCFLAGS) $(PSML_INCFLAGS) $(NCPS_INCFLAGS)
INCFLAGS:=$(INCFLAGS) $(XMLF90_INCFLAGS) $(GRIDXC_INCFLAGS)
#
OBJS = 	psml2psf.o handlers.o
#
#
MOD_OBJS= m_getopts.o f2kcli.o
#
COM_OBJS=$(OBJS) 
ALL_OBJS=$(MOD_OBJS) $(COM_OBJS)
#
$(COM_OBJS): 
#
what:
	echo "VPATH: " $(VPATH)
psml2psf:   $(NCPS)  $(ALL_OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o psml2psf $(ALL_OBJS) \
        $(NCPS) $(PSML_LIBS) $(XMLF90_LIBS)
#
psml2psf.o: m_getopts.o 
m_getopts.o: f2kcli.o
#
clean: 
	@echo "==> Cleaning object, library, and executable files"
	rm -f psml2psf *.o  *.a *.pcl *.pc *.mod
	rm -f *.mod
#
$(NCPS):
	@echo "Making sure that the ncps library is compiled..."
	(cd $(MAIN_OBJDIR)/ncps/src ; make \
                                     "ARCH_MAKE=$(ARCH_MAKE)"\
                                     "PSML_INCFLAGS=-$(PSML_INCFLAGS)"\
                                     "VPATH=$(VPATH)/ncps/src")
#



