# Makefile for kbgen
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#
.PHONY: clean wxml
#
default: psop
#
OBJDIR=Obj
TOP_LEVEL=$(shell cd -P -- ../../Src && pwd -P)
VPATH:=$(TOP_LEVEL)

MAIN_OBJDIR:=$(shell cd -P -- ../../\$(OBJDIR) && pwd -P)
ARCH_MAKE=$(MAIN_OBJDIR)/arch.make
include $(ARCH_MAKE)
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
NCPS_INCFLAGS=-I$(MAIN_OBJDIR)/ncps/src
PSML_INCFLAGS=-I$(MAIN_OBJDIR)/psml/src
PSOP_INCFLAGS=-I$(MAIN_OBJDIR)/psoplib/src

PSML=$(MAIN_OBJDIR)/psml/src/libpsml.a
NCPS=$(MAIN_OBJDIR)/ncps/src/libncps.a
PSOP=$(MAIN_OBJDIR)/psoplib/src/libpsop.a

INCFLAGS:=$(INCFLAGS) $(PSML_INCFLAGS) $(NCPS_INCFLAGS) $(PSOP_INCFLAGS)
INCFLAGS:=$(INCFLAGS) -I$(MAIN_OBJDIR)/wxml -I./SiestaXC # To pick up mod files
#
LIBS_WXML:= -L$(MAIN_OBJDIR)/wxml -l wxml
#
#
OBJS = 	psop_options.o\
	psop.o local_xml.o write_proj_psml.o \
        local_die.o psml_die.o semicore_info_froyen.o\
        dpnint.o check_grid.o
#
#
MOD_OBJS= m_getopts.o f2kcli.o
#
COM_OBJS=$(OBJS) 
ALL_OBJS=$(MOD_OBJS) $(COM_OBJS)
#
$(COM_OBJS): 
#
XC=libSiestaXC.a
XC_MAKEFILE=$(TOP_LEVEL)/SiestaXC/makefile
XC_INCFLAGS:=-I $(TOP_LEVEL)/SiestaXC ##$(INCFLAGS)
$(XC): 
	(mkdir -p SiestaXC; cd SiestaXC ; $(MAKE) -f $(XC_MAKEFILE) \
                          "FC=$(FC)" "VPATH_ROOT=$(TOP_LEVEL)" \
                          "ARCH_MAKE=$(ARCH_MAKE)" \
                          "INCFLAGS=$(XC_INCFLAGS)" \
                          "DEFS=$(DEFS) $(DEFS_PREFIX)-UMPI" \
                          "FPPFLAGS=$(FPPFLAGS) $(DEFS_PREFIX)-UMPI" \
                          "MPI_INTERFACE= " \
                          "FFLAGS=$(FFLAGS)"  module )
#
what:
	echo "VPATH: " $(VPATH)
psop:   wxml $(NCPS) $(PSML) $(PSOP) $(XC) $(ALL_OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o psop $(ALL_OBJS) \
        $(NCPS) $(PSOP) $(PSML) $(XC) $(LIBS_WXML)
#
$(NCPS): $(PSML)
#
psop.o: psop_options.o local_xml.o write_proj_psml.o
#
psop.o: m_getopts.o semicore_info_froyen.o
write_proj_psml.o: local_xml.o
m_getopts.o: f2kcli.o
#
clean: 
	@echo "==> Cleaning object, library, and executable files"
	rm -f psop *.o  *.a *.pcl *.pc *.mod
	rm -f *.mod
#
wxml:
	@echo "Making sure that the wxml library is compiled..."
	(cd $(MAIN_OBJDIR)/wxml ; make \
                                     "ARCH_MAKE=$(ARCH_MAKE)"\
                                     "VPATH=$(VPATH)/wxml")
$(PSML):
	@echo "Making sure that the psml library is compiled..."
	(cd $(MAIN_OBJDIR)/psml/src ; make \
                                     "ARCH_MAKE=$(ARCH_MAKE)"\
                                     "VPATH=$(VPATH)/psml/src")
$(PSOP):
	@echo "Making sure that the psop library is compiled..."
	(cd $(MAIN_OBJDIR)/psoplib/src ; make \
                                     "ARCH_MAKE=$(ARCH_MAKE)"\
                                     "VPATH=$(VPATH)/psoplib/src")
$(NCPS):
	@echo "Making sure that the ncps library is compiled..."
	(cd $(MAIN_OBJDIR)/ncps/src ; make \
                                     "ARCH_MAKE=$(ARCH_MAKE)"\
                                     "PSML_INCFLAGS=-I$(MAIN_OBJDIR)/psml/src"\
                                     "VPATH=$(VPATH)/ncps/src")
#



