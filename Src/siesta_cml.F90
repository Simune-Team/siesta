Module siesta_cml

  Use flib_wxml
  Use flib_wcml

  Implicit None
  Public
  
  Logical  :: cml_p = .False.
  Type(xmlf_t) :: mainXML

  Contains

    Subroutine siesta_cml_init( )
      Use fdf, Only : fdf_boolean, fdf_string
      Use parallel, only : nodes, ionode
      Use version_info
      Use time

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
         Call xml_NewElement(mainXML, 'cml')
         Call cmlStartMetadataList(mainXML)
         Call cmlAddMetadata(mainXML, name='Program', content='Siesta')
         Call cmlAddMetadata(mainXML, name='Version', content=version_str)
         Call cmlAddMetadata(mainXML, name='Arch',    content=siesta_arch)
         Call cmlAddMetadata(mainXML, name='Flags',   content=fflags)
         Call cmlAddMetadata(mainXML, name='Initial Timestamp',content=datestring) 
         If (nodes>1) Then
           Call cmlAddMetadata(mainXML, name='Mode', content='Parallel')
         Else
           Call cmlAddMetadata(mainXML, name='Mode', content='Serial')
         Endif
         Call cmlAddMetadata(mainXML, name='Nodes', content=nodes)
#ifdef CDF
         Call cmlAddMetadata(mainXML, name='NetCDF',  content='true')
#endif
         Call cmlEndMetadataList(mainXML)
      Endif !cml_p
      
    End Subroutine siesta_cml_init
       
    Subroutine siesta_cml_exit

      use time, only : datestring


      If (cml_p) Then
        Call cmlAddMetadata(mainXML, name='Final Timestamp',content=datestring)
        Call xml_EndElement(mainXML, 'cml')
        Call xml_Close(mainXML)
      Endif !cml_p

    End Subroutine siesta_cml_exit

End Module siesta_cml
