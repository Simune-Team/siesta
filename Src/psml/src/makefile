.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
.PHONY: xmlf90-sax
#
default:  libpsml.a
module:  libpsml.a
#
# We allow for an external specification of the location 
# of the arch.make file. This must be an absolute path
#
#
ARCH_MAKE_DEFAULT=fortran.mk
ARCH_MAKE?=$(ARCH_MAKE_DEFAULT)
include $(ARCH_MAKE)
#
xmlf90-sax:
	@echo "Compiling the xmlf90-sax library..."
	(cd xmlf90-sax ; $(MAKE)\
            ARCH_MAKE=../$(ARCH_MAKE)\
            VPATH=$(VPATH)/xmlf90-sax)
#
OBJFILES= m_psml_core.o m_psml_api.o \
          m_psml_parsing_helpers.o m_psml_reader.o \
          assoc_list.o external_interfaces.o m_interp.o\
          m_psml.o
#
INCFLAGS= -I xmlf90-sax
#
libpsml.a: xmlf90-sax $(OBJFILES)
	$(AR) $(ARFLAGS_EXTRA) cru libpsml.a $(OBJFILES) xmlf90-sax/*.o
	-$(RANLIB) libpsml.a
#
#dep:
#	sfmakedepend --depend=obj --modext=o *.f *.f90 *.F *.F90 
clean:
	rm -f *.o *.*d *.a
	(cd xmlf90-sax ; $(MAKE) ARCH_MAKE=../$(ARCH_MAKE) clean)
#
# DO NOT DELETE THIS LINE - used by make depend
#m_interpolation.o: external_interfaces.o
m_psml_api.o: assoc_list.o external_interfaces.o m_psml_core.o 
m_psml_api.o: m_interp.o
m_psml_core.o: assoc_list.o  external_interfaces.o
m_psml_parsing_helpers.o: assoc_list.o external_interfaces.o 
m_psml_parsing_helpers.o: m_psml_core.o
m_psml_reader.o: m_psml_core.o m_psml_parsing_helpers.o m_interp.o
#
