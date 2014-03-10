module m_psml

!
! In the current implementation, the PSML library code is
! scattered in three files belonging to the intermediate
! NCPS library. This module pulls in the relevant functionality.
!

! Notation: The 'psxml' prefix in the accessors should probably
! be changed to 'psml_', for consistency.

  use m_ncps_xml_ps_t,  only: psml_t => xml_ps_t
  use m_ncps_xml_ps_t,  only: psml_destroy => xml_ps_destroy

  use m_ncps_xmlreader, only: psml_reader => ncps_xmlreader

  use m_ncps_xml_ps_t, only: psxmlEvaluateValenceCharge
  use m_ncps_xml_ps_t, only: psxmlEvaluateCoreCharge
  use m_ncps_xml_ps_t, only: psxmlAtomicSymbol
  use m_ncps_xml_ps_t, only: psxmlAtomicNumber
  use m_ncps_xml_ps_t, only: psxmlCreator
  use m_ncps_xml_ps_t, only: psxmlDate
  use m_ncps_xml_ps_t, only: psxmlPseudoFlavor
  use m_ncps_xml_ps_t, only: psxmlPseudoZval
  use m_ncps_xml_ps_t, only: psxmlGenerationZval
  use m_ncps_xml_ps_t, only: psxmlXCFunctional
  use m_ncps_xml_ps_t, only: psxmlXCFunctionalType
  use m_ncps_xml_ps_t, only: psxmlIsRelativistic
  use m_ncps_xml_ps_t, only: psxmlIsSpinPolarized
  use m_ncps_xml_ps_t, only: psxmlHasCoreCorrections
  use m_ncps_xml_ps_t, only: psxmlHasGlobalLogGrid
  use m_ncps_xml_ps_t, only: psxmlGridNpoints
  use m_ncps_xml_ps_t, only: psxmlLogGridStep
  use m_ncps_xml_ps_t, only: psxmlLogGridScale
  use m_ncps_xml_ps_t, only: psxmlGridRmax
  use m_ncps_xml_ps_t, only: psxmlPotentialsUp
  use m_ncps_xml_ps_t, only: psxmlPotentialsDown
  use m_ncps_xml_ps_t, only: psxmlPotAngMomentum
  use m_ncps_xml_ps_t, only: psxmlOccupation
  use m_ncps_xml_ps_t, only: psxmlGenerationCutoff
  use m_ncps_xml_ps_t, only: psxmlPrincipalN
  use m_ncps_xml_ps_t, only: psxmlEvaluatePotential
  use m_ncps_xml_ps_t, only: psxmlEvaluatePseudoOrbital

  public

end module m_psml
