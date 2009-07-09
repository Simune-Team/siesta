! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2008.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
MODULE m_debug

! Initializes output files for debug info. J.M.Soler. Jan'2008

USE parallel,         only: node, Nodes  ! This node number, Number of nodes
USE moreParallelSubs, only: nodeString   ! Returns a string with a node number

PUBLIC &
  setDebugOutputUnit, &! Sets debug output unit and opens file for debug output
  udebug               ! Output file unit for debug info

PRIVATE ! Nothing is declared public beyond this point

  integer,save:: udebug=0  ! Output file unit for debug info

CONTAINS

subroutine setDebugOutputUnit()

  ! Sets debug output unit and opens file for debug output

  implicit none
  character(len=20):: fileName

  ! Set output file name, except node number, and find its name length
  fileName = 'debug.out'

  ! Append node number to file name
  if (Nodes>1) fileName = trim(fileName)//nodeString()

  ! Find an available output unit
  call io_assign( udebug )

  ! Open file and write header on it
  open( udebug, file=fileName )
  if (Nodes>1) &
    write(udebug,'(/,a)') 'Debug output for processor '//nodeString()

end subroutine setDebugOutputUnit

END MODULE m_debug
