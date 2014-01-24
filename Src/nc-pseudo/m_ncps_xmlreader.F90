module m_ncps_xmlreader

  public :: ncps_xmlreader

  CONTAINS

  subroutine ncps_xmlreader(fname,psxml)

  use m_ncps_xml_ps_t,        only: xml_ps_t, dump_pseudo
  use m_ncps_parsing_helpers, only: begin_element, end_element, pcdata_chunk
  use m_ncps_parsing_helpers, only: pseudo

#ifdef XMLF90
  use flib_sax,        only: xml_t, open_xmlfile, xml_parse
#else
  use FoX_sax,         only: xml_t, open_xml_file, parse
#endif
  implicit none 

  character(len=*), intent(in) :: fname
  type(xml_ps_t), intent(out)  :: psxml

  type(xml_t)                     :: fxml
  integer :: iostat

#ifdef XMLF90
 call open_xmlfile(fname,fxml,iostat)
 if (iostat /=0) stop "Cannot open XML file"
 call xml_parse(fxml, begin_element,end_element,pcdata_chunk,verbose=.false.)
#else
 call open_xml_file(fxml,fname,iostat)
 if (iostat /=0) stop "Cannot open XML file"
 call parse(fxml, startElement_handler=begin_element, &
                  endElement_handler=end_element,     &
                  characters_handler=pcdata_chunk) 
#endif

 psxml = pseudo  ! should this be a pointer assignment?
 call dump_pseudo(pseudo,6)

end subroutine ncps_xmlreader
end module m_ncps_xmlreader

