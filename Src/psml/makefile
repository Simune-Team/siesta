.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
module:  libpsml.a
# This makefile can also be used "remotely", so we allow
# for an external specification of the (relative) location 
# of the arch.make file.
#
#
ARCH_MAKE_DEFAULT=../arch.make
ARCH_MAKE?=$(ARCH_MAKE_DEFAULT)
include $(ARCH_MAKE)
#
#
OBJFILES= m_psml_core.o \
          m_psml_parsing_helpers.o m_psml_reader.o \
          assoc_list.o external_interfaces.o m_interpolation.o\
          m_psml.o
#
# Choice of XML parsers.
#

ifneq (,$(findstring DPSML_USE_FOX,$(FPPFLAGS)))
# If the user has put -DPSML_USE_FOX in the arch.make...

ifeq ($(strip $(FOX_ROOT)),)
# If FoX is not available as an external library
# (which would include the sax parser), then
# fall back to using XMLF90 in the Siesta distribution
 FPPFLAGS:=$(FPPFLAGS) $(DEFS_PREFIX)-UPSML_USE_FOX
 INCFLAGS:= $(INCFLAGS) -I../xmlparser
else
 include ../Mk/FoX.mk
 OBJFILES:=$(OBJFILES) fox_extra.o
 m_psml_parsing_helpers.o: fox_extra.o
endif
else
# The user has not requested FoX
 INCFLAGS:= $(INCFLAGS) -I../xmlparser
endif
#
libpsml.a: $(OBJFILES)
	$(AR) $(ARFLAGS_EXTRA) cru libpsml.a $(OBJFILES)
	-$(RANLIB) libpsml.a
	cp libpsml.a ..
	@cp -p *.*d ..
#
dep:
	sfmakedepend --depend=obj --modext=o *.f *.f90 *.F *.F90 
#
clean:
	rm -f *.o *.*d *.a
# DO NOT DELETE THIS LINE - used by make depend
m_interpolation.o: external_interfaces.o
m_psml_core.o: assoc_list.o external_interfaces.o m_interpolation.o
m_psml_parsing_helpers.o: assoc_list.o external_interfaces.o 
m_psml_parsing_helpers.o: m_psml_core.o
m_psml_reader.o: m_psml_core.o m_psml_parsing_helpers.o
