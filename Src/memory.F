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
      subroutine memory( Task, Type, NElements, CallingRoutine )
C 
C This subroutine keeps track of information relating to the use 
C of dynamic memory
C
C Input :
C
C character*1 Task  : job type = 'A' -> allocate
C                   :            'D' -> deallocate
C character*1 Type  : type of variable = 'I' = integer
C                   :                    'S' = single precision real
C                   :                    'D' = double precision real
C                   :                    'X' = grid precision real
C                   :                    'L' = logical
C                   :                    'C' = single precision complex
C                   :                    'Z' = double precision complex
C                   :                    'S' = character data (we assume takes one word)
C                   :                    'E' = double precision integer
C integer NElements : number of array elements being 
C                   : allocated/deallocated
C character         :
C   CallingRoutine  : string containing the name of the calling routine
C
C Created by J.D. Gale, October 1999
C
      use parallel, only: Node
      use alloc,    only: alloc_count

      implicit none

      integer, intent(in)                 :: NElements
      character(len=1), intent(in)        :: Task, Type
      character(len=*), intent(in)        :: CallingRoutine

C Local variables
      integer         :: allocSize
      character(len=1):: allocType

C Call alloc_count routine
      select case(Type)
      case('S')
        allocType = 'R'
        allocSize = NElements
      case('C')
        allocType = 'R'
        allocSize = NElements*2
      case('Z')
        allocType = 'D'
        allocSize = NElements*2
      case('X')
#ifdef GRID_DP
        allocType = 'D'
#else
        allocType = 'R'
#endif
        allocSize = NElements
      case default
        allocType = Type
        allocSize = NElements
      end select
      if (Task=='D') allocSize = -allocSize
      call alloc_count( allocSize, allocType, 
     .                  name=trim(CallingRoutine)//' unknown' )

      end subroutine memory
