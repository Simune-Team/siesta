# Makefile for kbgen
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#
.PHONY: clean
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
include $(MAIN_OBJDIR)/Mk/psml.mk
include $(MAIN_OBJDIR)/Mk/xmlf90.mk
#
NCPS_INCFLAGS=-I$(MAIN_OBJDIR)/ncps/src
PSOP_INCFLAGS=-I$(MAIN_OBJDIR)/psoplib/src

NCPS=$(MAIN_OBJDIR)/ncps/src/libncps.a
PSOP=$(MAIN_OBJDIR)/psoplib/src/libpsop.a

INCFLAGS:=$(INCFLAGS) $(PSML_INCFLAGS) $(NCPS_INCFLAGS) $(PSOP_INCFLAGS)
INCFLAGS:=$(INCFLAGS) $(XMLF90_INCFLAGS) 
INCFLAGS:=$(INCFLAGS)  -I./SiestaXC # To pick up mod files
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
psop:   $(NCPS) $(PSML) $(PSOP) $(XC) $(ALL_OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o psop $(ALL_OBJS) \
        $(NCPS) $(PSOP) $(PSML_LIBS) $(XC) $(XMLF90_LIBS)
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
$(PSOP):
	@echo "Making sure that the psop library is compiled..."
	(cd $(MAIN_OBJDIR)/psoplib/src ; make \
                                     "ARCH_MAKE=$(ARCH_MAKE)"\
                                     "VPATH=$(VPATH)/psoplib/src")
$(NCPS):
	@echo "Making sure that the ncps library is compiled..."
	(cd $(MAIN_OBJDIR)/ncps/src ; make \
                                     "ARCH_MAKE=$(ARCH_MAKE)"\
                                     "PSML_INCFLAGS=-$(PSML_INCFLAGS)"\
                                     "VPATH=$(VPATH)/ncps/src")
#



