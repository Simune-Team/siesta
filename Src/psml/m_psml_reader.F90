module m_psml_reader

  public :: psml_reader

  CONTAINS

  subroutine psml_reader(fname,ps)

  use m_psml_core,            only: ps_t, ps_destroy
  use m_psml_parsing_helpers, only: begin_element, end_element, pcdata_chunk
  use m_psml_parsing_helpers, only: pseudo

#ifdef PSML_USE_FOX
  use FoX_sax,           only: xml_t, open_xml_file, parse
#else
  use xmlf90_sax,        only: xml_t, open_xmlfile, xml_parse
#endif

  implicit none 

  character(len=*), intent(in) :: fname
  type(ps_t), pointer          :: ps

  type(xml_t)                     :: fxml
  integer :: iostat

  if (associated(ps)) then
     call ps_destroy(ps)
  endif

#ifdef PSML_USE_FOX
 call open_xml_file(fxml,fname,iostat)
 if (iostat /=0) call die("Cannot open XML file")
 call parse(fxml, startElement_handler=begin_element, &
                  endElement_handler=end_element,     &
                  characters_handler=pcdata_chunk) 
#else
 call open_xmlfile(fname,fxml,iostat)
 if (iostat /=0) call die("Cannot open XML file")
 call xml_parse(fxml, begin_element,end_element,pcdata_chunk,verbose=.false.)
#endif

 ps => pseudo

end subroutine psml_reader
end module m_psml_reader

