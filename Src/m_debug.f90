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
USE moreParallelSubs, only: copyFile     ! Copies a file to another node

PUBLIC &
  setDebugOutputUnit,   &! Sets debug output unit and opens it for debug output
  closeDebugOutputFile, &! Close debug output file and copies it to node 0
  udebug                 ! Output file unit for debug info

PRIVATE ! Nothing is declared public beyond this point

  character(len=*),parameter:: filePrefix = 'debug.out' ! Prefix of file name
  character(len=32),    save:: fileName  ! Output file name
  integer,              save:: udebug=0  ! Output file unit for debug info

CONTAINS

!-----------------------------------------------------------------------------

subroutine setDebugOutputUnit()

  ! Sets debug output unit and opens file for debug output

  implicit none

  ! Set output file name, except node number, and find its name length
  fileName = filePrefix

  ! Append node number to file name
  if (Nodes>1) fileName = filePrefix//nodeString()

  ! Find an available output unit
  call io_assign( udebug )

  ! Open file and write header on it
  open( udebug, file=fileName )
  if (Nodes>1) &
    write(udebug,'(/,a)') 'Debug output for processor '//nodeString()

end subroutine setDebugOutputUnit

!-----------------------------------------------------------------------------

subroutine closeDebugOutputFile()

  ! Closes debug output file and copies it to node 0

  implicit none
  integer:: iNode
  character(len=32):: nodeFileName

  close( unit=udebug )

! Copy all debug.out* files to root node. Notice that in each call to
! copyFile, only node iNode will actually send its file to node=0
  do iNode = 1,Nodes-1
    nodeFileName = filePrefix//nodeString(iNode)
    call copyFile( srcNode=iNode, srcFile=nodeFileName, &
                   dstNode=0,     dstFile=nodeFileName, &
                   writeOption='overwrite' )
  end do

end subroutine closeDebugOutputFile

END MODULE m_debug
