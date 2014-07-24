module m_psml

!
! In the current implementation, the PSML library code is
! scattered in three files belonging to the intermediate
! NCPS library. This module pulls in the relevant functionality.
!

! Notation: The 'psxml' prefix in the accessors should probably
! be changed to 'psml_', for consistency.


  use m_ncps_xmlreader, only: psml_reader => ncps_xmlreader

  use m_ncps_xml_ps_t

!!$  use m_ncps_xml_ps_t,  only: psml_t => ps_t
!!$  use m_ncps_xml_ps_t,  only: psml_destroy => ps_destroy
!!$
!!$  use m_ncps_xml_ps_t, only: ps_EvaluateValenceCharge
!!$  use m_ncps_xml_ps_t, only: ps_EvaluateCoreCharge
!!$  use m_ncps_xml_ps_t, only: ps_AtomicSymbol
!!$  use m_ncps_xml_ps_t, only: ps_AtomicNumber
!!$  use m_ncps_xml_ps_t, only: ps_Creator
!!$  use m_ncps_xml_ps_t, only: ps_Date
!!$  use m_ncps_xml_ps_t, only: ps_PseudoFlavor
!!$  use m_ncps_xml_ps_t, only: ps_PseudoZval
!!$  use m_ncps_xml_ps_t, only: ps_GenerationZval
!!$  use m_ncps_xml_ps_t, only: ps_XCFunctional
!!$  use m_ncps_xml_ps_t, only: ps_XCFunctionalType
!!$  use m_ncps_xml_ps_t, only: ps_IsRelativistic
!!$  use m_ncps_xml_ps_t, only: ps_IsSpinPolarized
!!$  use m_ncps_xml_ps_t, only: ps_HasCoreCorrections
!!$  use m_ncps_xml_ps_t, only: ps_HasGlobalLogGrid
!!$  use m_ncps_xml_ps_t, only: ps_GridNpoints
!!$  use m_ncps_xml_ps_t, only: ps_LogGridStep
!!$  use m_ncps_xml_ps_t, only: ps_LogGridScale
!!$  use m_ncps_xml_ps_t, only: ps_GridRmax
!!$  use m_ncps_xml_ps_t, only: ps_PotentialsUp
!!$  use m_ncps_xml_ps_t, only: ps_PotentialsDown
!!$  use m_ncps_xml_ps_t, only: ps_PotAngMomentum
!!$  use m_ncps_xml_ps_t, only: ps_Occupation
!!$  use m_ncps_xml_ps_t, only: ps_GenerationCutoff
!!$  use m_ncps_xml_ps_t, only: ps_PrincipalN
!!$  use m_ncps_xml_ps_t, only: ps_EvaluatePotential
!!$  use m_ncps_xml_ps_t, only: ps_EvaluatePseudoOrbital

  public

end module m_psml
