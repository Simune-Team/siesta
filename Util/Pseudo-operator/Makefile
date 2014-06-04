# Makefile for kbgen
#
.SUFFIXES: .f .F .o .a  .f90 .F90
#
.PHONY: ncps clean
#
default: pseudooperator
#
OBJDIR=Obj
VPATH:=$(shell pwd)/../../Src
MAIN_OBJDIR=../../$(OBJDIR)
include $(MAIN_OBJDIR)/arch.make
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
INCFLAGS:=$(INCFLAGS) -I$(MAIN_OBJDIR)/nc-pseudo # To pick up mod files
#
# If we have compiled libncps.a using XMLF90 instead of FoX,
# use the appropriate library
#
# First, check whether -DXMLF90 is already an option in the arch.make
# 
ifneq (,$(findstring DXMLF90,$(FPPFLAGS)))
 LIBS= $(MAIN_OBJDIR)/libncps.a $(MAIN_OBJDIR)/libxmlparser.a
else
#
# Check whether -DXMLF90 was added on-the-fly while compiling
# libncps.a. This happens if FoX is not available as an external
# library (which would include the sax parser), so we fall back to
# using XMLF90 in the Siesta distribution
#
ifeq ($(strip $(FOX_ROOT)),)
# FPPFLAGS:=$(FPPFLAGS) $(DEFS_PREFIX)-DXMLF90
# INCFLAGS:= $(INCFLAGS) -I$(MAIN_OBJDIR)/xmlparser
 LIBS= $(MAIN_OBJDIR)/libncps.a $(MAIN_OBJDIR)/libxmlparser.a
else
 include $(MAIN_OBJDIR)/Mk/FoX.mk
 LIBS= $(MAIN_OBJDIR)/libncps.a $(FOX_LIBS) 
endif
endif
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
pseudooperator:  ncps $(ALL_OBJS)
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



