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
