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
# Choice of XML parsers.
#
# First, check whether -DXMLF90 is already an option in the arch.make
# 
ifneq (,$(findstring DXMLF90,$(FPPFLAGS)))
 INCFLAGS:= $(INCFLAGS) -I../xmlparser
else
#
# If FoX is not available as an external library
# (which would include the sax parser), then
# fall back to using XMLF90 in the Siesta distribution
# 
ifeq ($(strip $(FOX_ROOT)),)
 FPPFLAGS:=$(FPPFLAGS) $(DEFS_PREFIX)-DXMLF90
 INCFLAGS:= $(INCFLAGS) -I../xmlparser
else
 include ../Mk/FoX.mk
 # This adds the proper string to INCFLAGS
endif
endif
#
OBJFILES= m_ncps_xml_ps_t.o \
          m_ncps_parsing_helpers.o m_ncps_xmlreader.o \
          assoc_list.o\
          m_psml.o
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
m_ncps_parsing_helpers.o: m_ncps_xml_ps_t.o assoc_list.o
m_ncps_xmlreader.o: m_ncps_parsing_helpers.o m_ncps_xml_ps_t.o
m_ncps_xml_ps_t.o: assoc_list.o
m_psml.o: m_ncps_xml_ps_t.o m_ncps_xmlreader.o