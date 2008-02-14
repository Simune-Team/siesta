      MODULE cellxc_mod
      PRIVATE
      PUBLIC          :: setGGA, GGA
      logical         :: GGA=.FALSE.
      CONTAINS
      subroutine setGGA( )
      use xcmod, only : nXCfunc, XCfunc
      use sys,   only : die
      implicit none
C     Local variables
      integer         :: nf
      
!------------------------------------------------------------------------- BEGIN
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
!--------------------------------------------------------------------------- END
      end subroutine setGGA

      END MODULE cellxc_mod
