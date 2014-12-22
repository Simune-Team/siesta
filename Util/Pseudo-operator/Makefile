# Makefile for kbgen
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#
.PHONY: ncps clean wxml psml
#
default: psop
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
INCFLAGS:=$(INCFLAGS) -I$(MAIN_OBJDIR)/psml # To pick up mod files
INCFLAGS:=$(INCFLAGS) -I$(MAIN_OBJDIR)/wxml # To pick up mod files
INCFLAGS:=$(INCFLAGS) -I$(MAIN_OBJDIR)/SiestaXC # To pick up mod files
#
# If we have compiled libncps.a using XMLF90 instead of FoX,
# use the appropriate library
#
FPPFLAGS:=$(FPPFLAGS) -DHAS_PSML_SUPPORT
#
ifneq (,$(findstring DPSML_USE_FOX,$(FPPFLAGS)))
 include $(MAIN_OBJDIR)/Mk/FoX.mk
 LIBS= $(MAIN_OBJDIR)/libpsml.a $(FOX_LIBS) 
else
#
 LIBS= $(MAIN_OBJDIR)/libpsml.a $(MAIN_OBJDIR)/libxmlparser.a
endif
LIBS:= $(LIBS) -L$(MAIN_OBJDIR)/SiestaXC -l SiestaXC
LIBS:= $(LIBS) -L$(MAIN_OBJDIR)/wxml -l wxml
#
# Note that machine-specific files are now in top Src directory.
#
OBJS = precision.o  \
	atmparams.o atom_options.o \
	psop.o module_pseudo.o schrodinger.o arw.o\
        local_die.o
#
#
MOD_OBJS=m_ncps.o m_ncps_froyen_ps_t.o	\
         m_ncps_froyen_reader.o	m_ncps_reader.o \
         m_ncps_translators.o m_ncps_writers.o \
         m_getopts.o f2kcli.o
#
COM_OBJS=$(OBJS) 
ALL_OBJS=$(MOD_OBJS) $(COM_OBJS)
#
$(COM_OBJS): 
#
what:
	echo "VPATH: " $(VPATH)
psop:  psml wxml $(ALL_OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o psop $(ALL_OBJS) $(LIBS)
TEST_OBJS= semicore_info.o test_psml.o local_die.o
#
test_psml:  psml $(TEST_OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o test_psml $(TEST_OBJS) $(LIBS)
#
arw.o:  precision.o 
module_pseudo.o: atmparams.o precision.o schrodinger.o 
psop.o: precision.o atmparams.o module_pseudo.o
schrodinger.o: arw.o precision.o
#
m_ncps.o: m_ncps_froyen_ps_t.o m_ncps_reader.o m_ncps_translators.o
m_ncps.o: m_ncps_writers.o
m_ncps_froyen_reader.o: flib_spline.o m_ncps_froyen_ps_t.o
m_ncps_reader.o: m_ncps_froyen_ps_t.o m_ncps_froyen_reader.o
m_ncps_reader.o: m_ncps_translators.o
m_ncps_translators.o: m_ncps_froyen_ps_t.o
m_ncps_writers.o: m_ncps_froyen_ps_t.o
psop.o: m_getopts.o
m_getopts.o: f2kcli.o
#
clean: 
	@echo "==> Cleaning object, library, and executable files"
	rm -f psop *.o  *.a *.pcl *.pc *.mod
	rm -f *.mod
#
psml:
	@echo "Making sure that the psml library is compiled..."
	(cd ../../$(OBJDIR)/psml ; make "VPATH=$(VPATH)/psml")
wxml:
	@echo "Making sure that the wxml library is compiled..."
	(cd ../../$(OBJDIR)/wxml ; make "VPATH=$(VPATH)/wxml")
#



