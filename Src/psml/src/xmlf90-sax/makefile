.SUFFIXES: 
.SUFFIXES: .f .F .o .f90 .F90
#
PSML_MK_DEFAULT=../fortran.mk
PSML_MK?=$(PSML_MK_DEFAULT)
include $(PSML_MK)
#
OBJFILES= m_debug.o m_buffer.o m_elstack.o m_dictionary.o m_io.o \
          m_reader.o m_converters.o m_charset.o  m_entities.o \
          m_xml_error.o  m_fsm.o  m_xml_parser.o xmlf90_sax.o
MODFILES=$(OBJFILES:.o=)
#
default: objs_and_modules
#
objs_and_modules: $(OBJFILES)
#
#------------------------------
# Uncomment one of the following lines for debugging/profiling support
#
#FFLAGS=$(FFLAGS_DEBUG)
#FFLAGS=$(FFLAGS_PROFILE)
#FFLAGS=$(FFLAGS_CHECK)
#
m_reader.o: m_io.o m_debug.o
m_elstack.o: m_buffer.o
m_dictionary.o: m_buffer.o
m_entities.o: m_buffer.o m_xml_error.o
m_converters.o: m_debug.o
m_fsm.o: m_buffer.o m_dictionary.o m_charset.o m_entities.o m_elstack.o
m_xml_error.o: m_elstack.o
m_xml_parser.o: m_buffer.o m_elstack.o m_reader.o m_fsm.o m_dictionary.o \
                m_debug.o m_xml_error.o
#
xmlf90_sax.o: m_dictionary.o m_xml_parser.o m_converters.o m_xml_error.o
#
clean:
	rm -f *.o  *.mod

