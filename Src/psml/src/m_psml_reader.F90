module m_psml_reader

  public :: psml_reader

  CONTAINS

  subroutine psml_reader(fname,ps,debug)

  use m_psml_core,            only: ps_t, ps_destroy
  use m_psml_parsing_helpers, only: begin_element, end_element, pcdata_chunk
  use m_psml_parsing_helpers, only: pseudo, debug_parsing
  use m_interp,               only: set_default_interpolator
  use external_interfaces,    only: die => psml_die

#ifdef PSML_USE_FOX
  use FoX_sax,           only: xml_t, open_xml_file, close_xml_t, parse
#else
  use xmlf90_sax,        only: xml_t, open_xmlfile, xml_parse, close_xmlfile
#endif

  implicit none 

  character(len=*), intent(in)        :: fname
  type(ps_t), intent(out), target     :: ps
  logical, intent(in), optional       :: debug

  type(xml_t)                     :: fxml
  integer :: iostat

  ! Clean the object's internal data
  call ps_destroy(ps)

  ! Associate module pointer, so that the parsed data
  ! is written to ps
  pseudo => ps
  if (present(debug)) then
     debug_parsing = debug
  else
     debug_parsing = .false.
  endif

#ifdef PSML_USE_FOX
 call open_xml_file(fxml,fname,iostat)
 if (iostat /=0) call die("Cannot open XML file")
 call parse(fxml, startElement_handler=begin_element, &
                  endElement_handler=end_element,     &
                  characters_handler=pcdata_chunk) 
 call close_xml_t(fxml)
#else
 call open_xmlfile(fname,fxml,iostat)
 if (iostat /=0) call die("Cannot open XML file")
 call xml_parse(fxml, begin_element,end_element,pcdata_chunk,verbose=.false.)
 call close_xmlfile(fxml)
#endif

 ! Clean up association of module pointer
 pseudo => null()
 !
 ! Set default interpolator
 !
 call set_default_interpolator()

end subroutine psml_reader
end module m_psml_reader

