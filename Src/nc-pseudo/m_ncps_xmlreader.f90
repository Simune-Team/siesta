module m_ncps_xmlreader

  public :: ncps_xmlreader

  CONTAINS

  subroutine ncps_xmlreader(fname,psxml)

  use m_ncps_xml_ps_t,        only: xml_ps_t
  use m_ncps_parsing_helpers, only: begin_element, end_element, pcdata_chunk
  use m_ncps_parsing_helpers, only: pseudo

  use flib_sax,        only: xml_t, open_xmlfile, xml_parse

  implicit none 

  character(len=*), intent(in) :: fname
  type(xml_ps_t), intent(out)  :: psxml

  type(xml_t)                     :: fxml
  integer :: iostat

 call open_xmlfile(fname,fxml,iostat)
 if (iostat /=0) stop "Cannot open XML file"
                                                                         
 call xml_parse(fxml, begin_element,end_element,pcdata_chunk,verbose=.false.)
 psxml = pseudo  ! should this be a pointer assignment?

end subroutine ncps_xmlreader
end module m_ncps_xmlreader

