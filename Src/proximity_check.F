      subroutine proximity_check( )
      USE precision, only:dp
      USE SIESTA_TODO
      USE siesta_options
      use siesta_geom
      use atomlist, only: rco, iphkb, iphorb, indxua, no_s, rckb, rmaxo,
     .                    rmaxv, no_u, lastkb, lasto, rmaxkb
      use parallel, only: IOnode
      implicit none

      integer :: ii, jj, jamin, idxneighJ, nneigbI
      real(dp):: rmin, r2min 


      ! Check if any two atoms are unreasonably close
      do ii = 1,na_s
        r2min = huge(1._dp)
        jamin = 0
        nneigbI = maxna
        call neighb( scell, rmax, na_s, xa, ii, isel,
     .               nneigbI, jna, xij, r2ij )
        do jj = 1,nneigbI
          idxneighJ = jna(jj)
          if ( r2ij(jj).lt.r2min .and. idxneighJ.ge.ii ) then
!           Check that it is not the same atom
            if ( idxneighJ.ne.ii .or. r2ij(jj).gt.1.d-12 ) then
              r2min = r2ij(jj)
              jamin = idxneighJ
            endif
          endif
        enddo
        rmin = sqrt( r2min )
        if (IOnode) then
          if ( rmin .lt. rijmin ) write(6,'(a,2i6,a,f12.6,a)')
     .      'siesta: WARNING: Atoms', ii, jamin, ' too close: rij =',
     .       rmin/Ang, ' Ang'
        endif
      enddo

      end subroutine proximity_check
