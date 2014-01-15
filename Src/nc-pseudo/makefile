.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
module:  libncps.a
# This makefile can also be used "remotely", so we allow
# for an external specification of the (relative) location 
# of the arch.make file.
#
#
ARCH_MAKE_DEFAULT=../arch.make
ARCH_MAKE?=$(ARCH_MAKE_DEFAULT)
include $(ARCH_MAKE)
#
INCFLAGS:= $(INCFLAGS) -I../xmlparser
#
OBJFILES= m_ncps_froyen_ps_t.o m_ncps_xml_ps_t.o \
          m_ncps_parsing_helpers.o m_ncps_xmlreader.o \
          m_ncps_froyen_reader.o\
          m_ncps_writers.o m_ncps_reader.o\
          m_ncps_translators.o \
          m_ncps.o
#
libncps.a: $(OBJFILES)
	$(AR) $(ARFLAGS_EXTRA) cru libncps.a $(OBJFILES)
	-$(RANLIB) libncps.a
	cp libncps.a ..
	@cp -p *.*d ..
#
dep:
	sfmakedepend --depend=obj --modext=o *.f *.f90 *.F *.F90 
#
clean:
	rm -f *.o *.*d *.a
# DO NOT DELETE THIS LINE - used by make depend
m_ncps.o: m_ncps_froyen_ps_t.o m_ncps_reader.o m_ncps_translators.o
m_ncps.o: m_ncps_writers.o m_ncps_xml_ps_t.o
m_ncps_froyen_reader.o: m_ncps_froyen_ps_t.o
m_ncps_parsing_helpers.o: m_ncps_xml_ps_t.o
m_ncps_reader.o: m_ncps_froyen_ps_t.o m_ncps_froyen_reader.o
m_ncps_reader.o: m_ncps_translators.o m_ncps_xml_ps_t.o m_ncps_xmlreader.o
m_ncps_translators.o: m_ncps_froyen_ps_t.o m_ncps_xml_ps_t.o
m_ncps_writers.o: m_ncps_froyen_ps_t.o
m_ncps_xmlreader.o: m_ncps_parsing_helpers.o m_ncps_xml_ps_t.o
