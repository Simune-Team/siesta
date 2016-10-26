! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
Module siesta_cml

  !
  ! Now using the CML module in xmlf90
  !
  Use xmlf90_cml, only: cmlStartModule, cmlEndModule
  Use xmlf90_cml, only: cmlNamespaceAttribute, cmlAddComment
  Use xmlf90_cml, only: cmlStartStep, cmlEndStep
  Use xmlf90_cml, only: cmlStartPropertyList, cmlEndPropertyList
  Use xmlf90_cml, only: cmlStartParameterList, cmlEndParameterList
  Use xmlf90_cml, only: cmlAddProperty, cmlAddLattice, cmlAddKPoint
  Use xmlf90_cml, only: cmlAddMolecule, cmlAddParameter, cmlAddCrystal

  Use xmlf90_wxml, only: xmlf_t  
  Use xmlf90_cml, only: cmlBeginFile, cmlStartCml
  Use xmlf90_cml, only: cmlStartMetadataList, cmlAddMetadata
  Use xmlf90_cml, only: cmlEndMetadataList, cmlEndCml, cmlFinishFile

  Implicit None
  Logical, public      :: cml_p = .False.
  Type(xmlf_t), public, save :: mainXML

  public :: cmlStartModule, cmlEndModule
  public :: cmlStartStep, cmlEndStep
  public :: cmlNamespaceAttribute, cmlAddComment
  public :: cmlStartPropertyList, cmlEndPropertyList
  public :: cmlStartParameterList, cmlEndParameterList
  public :: cmlAddProperty, cmlAddLattice, cmlAddKPoint
  public :: cmlAddMolecule, cmlAddParameter, cmlAddCrystal
  public :: cmlBeginFile, cmlStartCml
  public :: cmlStartMetadataList, cmlAddMetadata
  public :: cmlEndMetadataList, cmlEndCml, cmlFinishFile

  private
  
End Module siesta_cml
