! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module diagmemory
C
C  Stores the factor used to scale the default memory in rdiag/cdiag
C  By increasing this value it is possible to avoid failure to
C  converge eigenvalues.
C
C  real*8  MemoryFactor      : factor by which memory is scaled
C
      use precision, only : dp

      implicit none

      real(dp), save :: MemoryFactor

      end module diagmemory
