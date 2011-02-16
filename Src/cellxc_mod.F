      MODULE cellxc_mod
      PRIVATE
      PUBLIC          :: setGGA, GGA
      logical         :: GGA=.FALSE.
      CONTAINS
      subroutine setGGA( )
#ifndef BSC_CELLXC
      use siestaXC, only : getXC        ! Returns the XC functional used
#else /* BSC_CELLXC */
      use bsc_xcmod, only : nXCfunc, XCfunc
#endif /* BSC_CELLXC */
      use sys,   only : die
      implicit none
C     Local variables
#ifndef BSC_CELLXC
      integer         :: nf, nXCfunc
      character(len=20):: XCauth(10), XCfunc(10)
#else /* BSC_CELLXC */
      integer         :: nf
#endif /* BSC_CELLXC */

!---------------------------------------------------------------------- BEGIN
#ifndef BSC_CELLXC
      call getXC( nXCfunc, XCfunc, XCauth )
#endif /* ! BSC_CELLXC */
      do nf = 1,nXCfunc
        if ( XCfunc(nf).eq.'GGA' .or. XCfunc(nf).eq.'gga') then
          GGA = .true.
        else
          if ( XCfunc(nf).ne.'LDA' .and. XCfunc(nf).ne.'lda' .and.
     &         XCfunc(nf).ne.'LSD' .and. XCfunc(nf).ne.'lsd' ) then
            write(6,*) 'cellXC: Unknown functional ', XCfunc(nf)
            call die()
          endif
        endif
      enddo
!----------------------------------------------------------------------- END
      end subroutine setGGA

      END MODULE cellxc_mod
