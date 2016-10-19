! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine iokp( nk, points, weight )
c *******************************************************************
c Saves k-points (only writing) Bohr^-1
c Emilio Artacho, Feb. 1999
c ********** INPUT **************************************************
c integer nk           : Number k points
c real*8  points(3,nk) : k-point coordinates
c real*8  weight(3,nk) : k-point weight
c *******************************************************************

      use files,     only : slabel, label_length
      use precision, only : dp

      implicit          none

      integer                       :: nk
      real(dp)                      :: points(3,*), weight(*)
      external          io_assign, io_close

c Internal 
      character(len=label_length+3), save :: fname
      integer                             :: ik, iu, ix
      logical,                       save :: frstme = .true.
c -------------------------------------------------------------------

      if (frstme) then
        fname = trim(slabel)//'.KP'
        frstme = .false.
      endif

      call io_assign( iu )
      open( iu, file=fname, form='formatted', status='unknown' )      

      write(iu,'(i6)') nk
      write(iu,'(i6,3f12.6,3x,f12.6)')
     .     (ik, (points(ix,ik),ix=1,3), weight(ik), ik=1,nk)

      call io_close( iu )

      return
      end
