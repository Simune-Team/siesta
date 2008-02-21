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

  USE parallel, only: node, Nodes  ! This node number, Number of nodes

PUBLIC

  integer,save:: udebug=0  ! Output file unit for debug info

CONTAINS

subroutine setDebugOutputUnit()

  ! Sets debug output unit and opens file for debug output

  implicit none
  character(len=20):: fileName, numName
  integer:: fileLen, maxLen, numLen

  ! Set output file name, except node number, and find its name length
  fileName = 'debug.out'
  fileLen = len_trim(fileName)

  ! Append node number to file name
  if (Nodes>1) then

    ! Find maximum name length of node numbers
    write(numName,*) Nodes-1   ! 0 <= node <= Nodes-1
    numName = adjustl(numName)
    maxLen = len_trim(numName)

    ! Find name of this node's number, and its name length
    write(numName,*) node
    numName = adjustl(numName)
    numLen = len_trim(numName)

    ! Set full file name, including node number
    ! e.g. debug.out00, debug.out01, ... debug.out63
    fileName = trim(fileName)//'000000000'
    fileName(fileLen+maxLen-numLen+1:fileLen+maxLen) = trim(numName)
    fileName(fileLen+maxLen+1:) = ' '

  end if ! (Nodes>1)

  ! Find an available output unit
  call io_assign( udebug )

  ! Open file and write header on it
  open( udebug, file=fileName )
  if (Nodes>1) &
    write(udebug,'(/,a)') 'Debug output for processor '//trim(numName)

end subroutine setDebugOutputUnit

END MODULE m_debug
