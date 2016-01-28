! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      MODULE m_diagon
      use precision
      implicit none
      logical           :: spiral, frstme = .true.,
     &                     AllInOne=.false., DivideConquer=.true.,
     &                     DoRead=.true., FirstCall=.true.,
     &                     NoExpert=.false., Serial=.true.,
     &                     PreRotate=.false., SaveEigenvectors=.false.,
     &                     Use2D=.true.
      real(dp)          :: qspiral(3)
      integer           :: ictxt, i2d_ctxt
      integer           :: maxclustersize = 12
#ifdef MPI
      logical           :: ParallelOverK
#endif
      END MODULE m_diagon
