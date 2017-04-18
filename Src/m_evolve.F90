!
! Copyright (C) 1996-2016 The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
MODULE m_evolve
  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: evolve

  CONTAINS
!*********************************************************************    
      SUBROUTINE evolve ( dt_tded )         

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
  use parallel,          only : Node, Nodes, ParallelOverK
  use parallelsubs,      only : GlobalToLocalOrb, GetNodeOrbs
  use fdf
  use alloc
  use m_memory
  use sys,               only : die
  use MatrixSwitch
  use atomlist,          only : no_s, no_l, no_u
  use m_spin,            only : nspin
  use m_gamma,           only : gamma
  use sparse_matrices,   only : maxnh, Escf
  use kpoint_grid,       only : nkpnt, kpoint
#ifdef MPI
  use mpi_siesta,        only : MPI_Bcast, MPI_Comm_World,MPI_logical
#endif
  !
  implicit none
  !
  real(dp), intent(in)     ::   dt_tded
  !
  logical, save            ::  frstme  = .true.
  integer                  ::   nuo
  !
#ifdef MPI
  call GetNodeOrbs(no_u,Node,Nodes,nuo)
  if (frstme) then
    if (ParallelOverK) then
      call die('TDDFT: Not prepared for running parallel over Kpoints.')
    end if
  endif
#else
  Node = 0
  Nodes = 1
  nuo = no_u
#endif
  ! Start time counter ................................................
  call timer( 'evolve', 1 )
  !
  ! Call apropriate routine .............................................
  if (nspin.le.2 .and. gamma) then
    call evolg( nspin, nuo, no_s, no_u, maxnh,                   &
                Escf, no_u, dt_tded)
  elseif (nspin.le.2 .and. .not.gamma) then
    stop 'evolve: Error: kpoint for TDED is not implimented yet'
  !       call evolk( nspin, maxspn, nuo, no, maxo, maxnh, maxnd,           &
  !                   indxuo, nk, kpoint, wk, Dnew, Enew, Dk, Ek,           &
  !                   no_u, dt_tded,Haux,Saux,psi)

  elseif (nspin.eq.4 .and. gamma) then 
    stop 'evolve: ERROR: non-collinear spin not yet implemented'
  elseif (nspin.eq.4 .and. .not. gamma) then 
    stop 'evolve: ERROR: non-collinear spin not yet implemented'
  endif 
  ! Stop time counter ...................................................
  call timer( 'evolve', 2 )
  !
  END SUBROUTINE  evolve
  !
  END module m_evolve
