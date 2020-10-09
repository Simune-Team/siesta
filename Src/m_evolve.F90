!
! Copyright (C) 1996-2016 The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module m_evolve

  implicit none
  private

  public :: evolve

contains

!*********************************************************************
  subroutine evolve( dt_tded )
! ********************************************************************
! Subroutine to time-evolve the eigenvectors, calculate the density
! and energy-density matrices, and occupation weights of each
! eigenvector, for given Hamiltonian and Overlap matrices (including
! spin polarization).
! Written by A. Tsolakidis, May 2000 after a suboutine
! by P. Ordejon and J. M. Soler.
! Modified by D. Sanchez-Portal, November 2002
! Modified by D. Sanchez-Portal, March 2008
! Modified by Rafi Ullah, October ,2015 making it parallel using
! Matrix Swtich.
! *******************************************************************

  use precision
  use sys,               only : die
  use m_spin,            only : nspin
  use cranknic_evolg,    only : cn_evolg
  use cranknic_evolk,    only : cn_evolk
  use atomlist,          only : no_s, no_u

  real(dp), intent(in)     ::   dt_tded

  logical :: not_using_auxcell

#ifdef DEBUG
  call write_debug( '    PRE evolve' )
#endif

  not_using_auxcell = (no_u == no_s)

  ! Call apropriate routine .............................................
  if ( nspin <= 2 ) then
    if ( not_using_auxcell ) then
      call cn_evolg( dt_tded )
    else
      call cn_evolk( dt_tded )
    end if
  else
    call die( 'evolve: TDDFT: non-collinear/spin-orbit spin not yet implemented' )
  endif

#ifdef DEBUG
  call write_debug( '    POS evolve' )
#endif

  END SUBROUTINE  evolve

END module m_evolve
