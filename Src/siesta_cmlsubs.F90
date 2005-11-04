Module siesta_cmlsubs

  Use flib_wxml, only: xmlf_t      ! help pgf95...
  Use flib_wxml
  Use flib_wcml

  Implicit None
  Private
  
  public :: siesta_cml_init, siesta_cml_exit

  Logical, public      :: cml_p = .False.
  Type(xmlf_t), public :: mainXML

  Contains

    Subroutine siesta_cml_init( )
      Use fdf, Only : fdf_boolean, fdf_string
      Use parallel, only : nodes, ionode
      Use version_info
      Use m_timestamp, only: datestring

      Character(len=25) :: sname
      Character(len=29) :: fname

      sname=''
      fname=''

      If (IOnode) Then
         cml_p = fdf_boolean( 'WriteXML', .True. )
      Else
         cml_p = .False.
      Endif !IOnode

      If (cml_p) Then
         sname = fdf_string('SystemLabel','siesta')
         Write(fname,'(a,a)') Trim(sname),'.xml'
         Call xml_OpenFile(trim(fname), mainXML, .True.)
         Call xml_AddXMLDeclaration(mainxml, 'UTF-8')
         Call xml_AddXMLStylesheet(mainXML, 'display.xsl', 'text/xsl')
         Call xml_NewElement(mainXML, 'cml')
         Call xml_AddAttribute(mainXML, 'xmlns', 'http://www.xml-cml.org/schema/CML2/Core')
         Call xml_AddAttribute(mainXML, name='xmlns:siesta', value='http://www.uam.es/siesta/namespace')
         Call cmlStartMetadataList(mainXML)
         Call cmlAddMetadata(mainXML, name='siesta:Program', content='Siesta')
         Call cmlAddMetadata(mainXML, name='siesta:Version', content=version_str)
         Call cmlAddMetadata(mainXML, name='siesta:Arch',    content=siesta_arch)
         Call cmlAddMetadata(mainXML, name='siesta:Flags',   content=fflags)
         Call cmlAddMetadata(mainXML, name='siesta:StartTime',content=datestring()) 
         If (nodes>1) Then
           Call cmlAddMetadata(mainXML, name='siesta:Mode', content='Parallel')
         Else
           Call cmlAddMetadata(mainXML, name='siesta:Mode', content='Serial')
         Endif
         Call cmlAddMetadata(mainXML, name='siesta:Nodes', content=nodes)
#ifdef CDF
         Call cmlAddMetadata(mainXML, name='siesta:NetCDF',  content='true')
#endif
         Call cmlEndMetadataList(mainXML)
      Endif !cml_p
      
    End Subroutine siesta_cml_init
       
    Subroutine siesta_cml_exit

      use m_timestamp, only : datestring


      If (cml_p) Then
        Call cmlAddMetadata(mainXML, name='siesta:EndTime',content=datestring())
        Call xml_EndElement(mainXML, 'cml')
        Call xml_Close(mainXML)
      Endif !cml_p

    End Subroutine siesta_cml_exit

End Module siesta_cmlsubs
