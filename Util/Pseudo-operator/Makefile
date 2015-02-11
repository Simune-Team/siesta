# Makefile for kbgen
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#
.PHONY: ncps clean wxml
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
# Note that PSML is now considered an external library
# Appropriate symbols should be put in arch.make
#
INCFLAGS:=$(INCFLAGS) $(PSML_INCFLAGS)
#
INCFLAGS:=$(INCFLAGS) -I$(MAIN_OBJDIR)/wxml # To pick up mod files
INCFLAGS:=$(INCFLAGS) -I$(MAIN_OBJDIR)/SiestaXC # To pick up mod files
#
LIBS_XC:= -L$(MAIN_OBJDIR)/SiestaXC -l SiestaXC
LIBS_WXML:= -L$(MAIN_OBJDIR)/wxml -l wxml
#
#
OBJS = precision.o  \
	psop_params.o psop_options.o\
	psop.o module_pseudo.o schrodinger.o arw.o\
        local_die.o psml_die.o semicore_info_froyen.o
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
psop:   wxml $(ALL_OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o psop $(ALL_OBJS) \
        $(PSML_LIBS) $(LIBS_XC) $(LIBS_WXML)

#
arw.o:  precision.o 
module_pseudo.o: psop_params.o precision.o schrodinger.o psop_options.o
psop.o: precision.o psop_params.o module_pseudo.o psop_options.o
schrodinger.o: arw.o precision.o psop_params.o
#
m_ncps.o: m_ncps_froyen_ps_t.o m_ncps_reader.o m_ncps_translators.o
m_ncps.o: m_ncps_writers.o
m_ncps_froyen_reader.o: flib_spline.o m_ncps_froyen_ps_t.o
m_ncps_reader.o: m_ncps_froyen_ps_t.o m_ncps_froyen_reader.o
m_ncps_reader.o: m_ncps_translators.o
m_ncps_translators.o: m_ncps_froyen_ps_t.o
m_ncps_writers.o: m_ncps_froyen_ps_t.o
psop.o: m_getopts.o semicore_info_froyen.o
m_getopts.o: f2kcli.o
#
clean: 
	@echo "==> Cleaning object, library, and executable files"
	rm -f psop *.o  *.a *.pcl *.pc *.mod
	rm -f *.mod
#
wxml:
	@echo "Making sure that the wxml library is compiled..."
	(cd ../../$(OBJDIR)/wxml ; make "VPATH=$(VPATH)/wxml")
#



