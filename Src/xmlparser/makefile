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
#
#
PARENT_VPATH=`grep '^VPATH' ../Makefile`
vpath_warning:
	@echo "** If you are using VPATH, compilations launched from this directory"
	@echo "** need to provide an explicit VPATH to the 'make' command."
	@echo "** Use the VPATH in Makefile in the top directory, followed by /xmlparser."
	@echo "** For example, to compile the XML tester program:"
	@echo
	@echo   make $(PARENT_VPATH)/xmlparser test-xml
	@echo

test.o: f2kcli.o m_getopts.o m_compare.o
m_compare.o: m_compare_tol.o corresponding_node.o string_utilities.o
corresponding_node.o: string_utilities.o
string_utilities.o: m_compare_tol.o	
m_getopts.o: f2kcli.o
m_compare_tol.o: $(DOM_OBJ_FILES)

OTHER_OBJS=../pxf.o
#
test-xml: vpath_warning  \
          test.o f2kcli.o m_getopts.o $(DOM_OBJ_FILES) corresponding_node.o string_utilities.o libxmlparser.a ../libwxml.a $(OTHER_OBJS)
	@echo
	@echo  "-- Building test-xml program. "
	@echo  "-- Some compilers might have trouble with the DOM subsystem. "
	@echo  "-- Make sure also that 'INCFLAGS' is part of the default compilation rules in arch.make"
	@echo
	$(FC) $(FFLAGS) -o test-xml test.o f2kcli.o m_getopts.o \
                           compare_m.o m_compare_tol.o corresponding_node.o string_utilities.o \
                           $(DOM_OBJ_FILES) libxmlparser.a ../libwxml.a $(OTHER_OBJS)

list: vpath_warning \
       list_nodes.o $(DOM_OBJ_FILES) libxmlparser.a ../libwxml.a $(OTHER_OBJS)
	$(FC) $(FFLAGS)  -o list list_nodes.o \
                            $(DOM_OBJ_FILES) libxmlparser.a ../libwxml.a $(OTHER_OBJS)

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

