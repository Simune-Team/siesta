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
module parallelsubswan

  implicit none

  contains

  subroutine GetNodeProjs( numproj, Node, Nodes, NProjsNode)
  use parallel, only: thisNode => Node
!
!  Calculates the number of orbitals stored on the local Node.
!
!  Julian Gale, October 1998
!
!  Input :
!
!  integer numproj  = The total number of projections in the calculation
!  integer Node     = The local processor
!  integer Nodes    = The total number of processors
!
!  Output :
!
!  integer NProjsNode = The number of projections that will be computed 
!                       on this Node - if zero
!                       on input then calculated otherwise left unchanged
!
!

! Passed arguments
  integer numproj, NProjsNode, Node, Nodes

! Local variables
  integer Remainder, MinPerNode, iproj, globalindexproj

! Calculate the minimum number of projections per node
  MinPerNode = numproj / Nodes

! Find the remainder of unassigned projections
  Remainder = numproj - MinPerNode * Nodes 

! Workout the local number of orbitals
  NProjsNode = MinPerNode 
  if (Node.lt.Remainder) NProjsNode = NProjsNode + 1

  do iproj = 1, NProjsNode
    if ( Remainder .eq. 0 ) then
      globalindexproj = MinPerNode * Node + iproj
    else 
      if ( Node .lt. Remainder ) then 
         globalindexproj = (MinPerNode+1) * Node + iproj 
      else
         globalindexproj = ( MinPerNode + 1 ) * Remainder + &
 &                         ( Node - Remainder )+ 1
      endif
    endif
    write(6,*)' Node, iproj, globalindexproj = ',  &
                Node, iproj, globalindexproj 
  enddo

  end subroutine GetNodeProjs

end module parallelsubswan
