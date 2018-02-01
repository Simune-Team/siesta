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
NCPS_INCFLAGS=-I$(MAIN_OBJDIR)/ncps/src
PSOP_INCFLAGS=-I$(MAIN_OBJDIR)/psoplib/src

NCPS=$(MAIN_OBJDIR)/ncps/src/libncps.a
PSOP=$(MAIN_OBJDIR)/psoplib/src/libpsop.a

INCFLAGS:=$(INCFLAGS) $(PSML_INCFLAGS) $(NCPS_INCFLAGS) $(PSOP_INCFLAGS)
INCFLAGS:=$(INCFLAGS) $(XMLF90_INCFLAGS) $(GRIDXC_INCFLAGS)
#
OBJS = 	psop_options.o\
	psop.o m_kb.o store_proj_psml.o \
        local_die.o psml_die.o handlers.o semicore_info_froyen.o\
        check_grid.o m_uuid.o \
        m_libxc_list.o m_libxc_sxc_translation.o \
        m_siestaxc_list.o

#
#
MOD_OBJS= m_getopts.o f2kcli.o
#
COM_OBJS=$(OBJS) 
ALL_OBJS=$(MOD_OBJS) $(COM_OBJS)
#
$(COM_OBJS): 
#
XC=$(GRIDXC_LIBS)
#
what:
	echo "VPATH: " $(VPATH)
psop:   $(NCPS) $(PSOP) $(ALL_OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o psop $(ALL_OBJS) \
        $(NCPS) $(PSOP) $(PSML_LIBS) $(XC) $(XMLF90_LIBS)
#
psop.o: psop_options.o store_proj_psml.o m_uuid.o m_kb.o
psop.o: m_getopts.o semicore_info_froyen.o
psop.o: m_libxc_sxc_translation.o 
store_proj_psml.o: m_kb.o
m_libxc_sxc_translation.o: m_libxc_list.o  m_siestaxc_list.o
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



