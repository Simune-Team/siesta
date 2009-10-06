! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
Module siesta_cmlsubs

  Use FoX_common, only: str
  Use FoX_wxml, only: xmlf_t      ! help pgf95...
  Use FoX_wcml

  Implicit None
  Private
  
  public :: siesta_cml_init, siesta_cml_exit

  Logical, public      :: cml_p = .False.
  Type(xmlf_t), public, save :: mainXML

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
         Write(fname,'(a)') Trim(slabel)//'.xml'
         Call cmlBeginFile(mainXML, trim(fname), unit=-1)
         Call cmlAddNamespace(mainXML, 'siesta', 'http://www.uam.es/siesta/namespace')
         Call cmlAddNamespace(mainXML, 'siestaUnits', 'http://www.uam.es/siesta/namespace/units')
         Call cmlStartCml(mainXML, convention="CMLComp")
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
         Call cmlAddMetadata(mainXML, name='siesta:Nodes', content=str(nodes))
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
        Call cmlEndCml(mainXML)
        Call cmlFinishFile(mainXML)
      Endif !cml_p

    End Subroutine siesta_cml_exit

End Module siesta_cmlsubs
