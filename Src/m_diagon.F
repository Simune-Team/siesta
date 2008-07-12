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
      real(dp), pointer :: h2d(:,:), s2d(:,:), z2d(:,:), Zsave(:,:)
#ifdef MPI
      logical           :: ParallelOverK
#endif

      CONTAINS

      subroutine resetDiagonPointers( )
      use alloc, only : de_alloc
      implicit none
      CHARACTER(LEN=*), PARAMETER :: MYNAME =  'resetDiagonPointers'
      if (associated(h2d))
     &  call de_alloc( h2d, name='h2d', routine=MYNAME )

      if (associated(s2d))
     &  call de_alloc( s2d, name='s2d', routine=MYNAME )

      if (associated(z2d))
     &  call de_alloc( z2d, name='z2d', routine=MYNAME )

      if (associated(Zsave))
     &  call de_alloc( Zsave, name='Zsave', routine=MYNAME )

      end subroutine resetDiagonPointers

      END MODULE m_diagon
