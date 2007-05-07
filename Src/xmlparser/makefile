.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
module:  libxmlparser.a # Do not build by default: test-xml
include ../arch.make
#
OBJFILES= m_debug.o m_buffer.o m_elstack.o m_dictionary.o m_io.o \
          m_reader.o m_converters.o m_charset.o  m_entities.o \
          m_xml_error.o  m_fsm.o  m_xml_parser.o flib_sax.o

DOM_OBJ_FILES= flib_dom.o \
          m_dom_attribute.o m_dom_element.o m_dom_node.o m_dom_types.o \
          m_dom_debug.o m_dom_error.o m_dom_nodelist.o m_dom_utils.o \
          m_dom_document.o m_dom_namednodemap.o m_dom_parse.o \
          m_strings.o
#
libxmlparser.a: $(OBJFILES) 
	ar $(ARFLAGS_EXTRA) cru libxmlparser.a $(OBJFILES)
	-$(RANLIB) libxmlparser.a
	cp libxmlparser.a ..
	@cp -p *.*d ..
	@echo "xmlparser/work.pc" >> ../work.pcl
#
#         (Note trick for Intel Fortran Compiler)
#
INCFLAGS=-I../wxml

test.o: f2kcli.o
m_compare_tol.o: $(DOM_OBJ_FILES)

test-xml:test.o f2kcli.o $(DOM_OBJ_FILES) libxmlparser.a ../libwxml.a
	$(FC) $(FFLAGS) -o test-xml test.o f2kcli.o \
                           compare_m.o m_compare_tol.o \
                           $(DOM_OBJ_FILES) libxmlparser.a ../libwxml.a

list: list_nodes.o $(DOM_OBJ_FILES) libxmlparser.a ../libwxml.a
	$(FC) $(FFLAGS)  -o list list_nodes.o \
                            $(DOM_OBJ_FILES) libxmlparser.a ../libwxml.a

clean:
	rm -f *.o *.*d *.a test-xml list
#
# DO NOT DELETE THIS LINE - used by make depend
m_reader.o: m_io.o m_debug.o
m_elstack.o: m_buffer.o
m_dictionary.o: m_buffer.o
m_entities.o: m_buffer.o m_xml_error.o
m_converters.o: m_debug.o
m_fsm.o: m_buffer.o m_dictionary.o m_charset.o m_entities.o m_elstack.o
m_xml_error.o: m_elstack.o
m_xml_parser.o: m_buffer.o m_elstack.o m_reader.o m_fsm.o m_dictionary.o \
                m_debug.o m_xml_error.o
m_dom_document.o:  m_dom_types.o
m_dom_nodelist.o: m_dom_types.o
m_dom_namednodemap.o:  m_dom_types.o
m_dom_node.o:  m_dom_types.o m_dom_nodelist.o m_dom_namednodemap.o \
              m_dom_debug.o m_dom_error.o
m_dom_attribute.o: m_dom_types.o  m_dom_node.o
m_dom_element.o: m_dom_types.o m_dom_nodelist.o m_dom_namednodemap.o \
                 m_dom_attribute.o m_dom_document.o m_dom_debug.o m_dom_node.o
m_dom_parse.o: m_dom_types.o m_dom_element.o m_dom_document.o m_dom_node.o \
               m_dom_debug.o
m_dom_utils.o: m_dom_types.o m_dom_element.o m_dom_document.o m_dom_node.o \
                m_dom_namednodemap.o m_dom_debug.o
m_dom_types.o:m_strings.o
flib_dom.o:  m_dom_types.o m_dom_nodelist.o m_dom_namednodemap.o \
             m_dom_attribute.o m_dom_document.o m_dom_node.o \
             m_dom_node.o m_dom_element.o m_dom_parse.o m_dom_utils.o
compare_m.o: m_compare_tol.o
test.o: compare_m.o

