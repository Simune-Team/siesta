! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
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
      Use fdf,   Only : fdf_boolean, fdf_string
      Use files, only : slabel, label_length
      Use parallel, only : nodes, ionode
      Use version_info
      Use m_timestamp, only: datestring

      Character(len=label_length+4) :: fname

      fname = ' '

      If (IOnode) Then
         cml_p = fdf_boolean( 'WriteXML', .True. )
      Else
         cml_p = .False.
      Endif !IOnode

      If (cml_p) Then
         Write(fname,'(a,a)') Trim(slabel),'.xml'
         Call xml_OpenFile(trim(fname), mainXML, .True.)
         Call xml_AddXMLDeclaration(mainxml, 'UTF-8')
         Call xml_AddXMLStylesheet(mainXML, 'http://www.eminerals.org/XSLT/display.xsl', 'text/xsl')
         Call xml_NewElement(mainXML, 'cml')
         Call xml_AddAttribute(mainXML, 'xmlns', 'http://www.xml-cml.org/schema')
         call xml_AddAttribute(mainXML, 'xmlns:xsd', 'http://www.w3.org/2001/XMLSchema')
         call xml_AddAttribute(mainXML, 'xmlns:dc', 'http://purl.org/dc/elements/1.1/title') 
         Call xml_AddAttribute(mainXML, 'xmlns:siesta', 'http://www.uam.es/siesta/namespace')
         Call xml_AddAttribute(mainXML, 'xmlns:siestaUnits', 'http://www.uam.es/siesta/namespace/units')
         Call xml_AddAttribute(mainXML, 'xmlns:eMinerals', 'http://www.eminerals.org/namespace') 
         Call cmlAddMetadata(mainXML, name='eMinerals:cmlSubsetVersion', content='0.9') 
         Call cmlAddMetadata(mainXML, name='dc:contributor', content='xmlf90 v1.99')
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
#else
         Call cmlAddMetadata(mainXML, name='siesta:NetCDF',  content='false')
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
