c $Id: iokp.f,v 1.1 1999/02/21 00:09:48 emilio Exp $

      subroutine iokp( nk, points, weight )

c *******************************************************************
c Saves k-points (only writing) Bohr^-1
c Emilio Artacho, Feb. 1999
c ********** INPUT **************************************************
c integer nk           : Number k points
c real*8  points(3,nk) : k-point coordinates
c real*8  weight(3,nk) : k-point weight
c *******************************************************************

      implicit          none
      character         paste*33
      integer           nk
      double precision  points(3,*), weight(*)
      external          io_assign, io_close, paste
      include          'fdf/fdfdefs.h'

c Internal 
      character  sname*30, fname*33
      integer    ik, iu, ix
      logical    frstme
      save       frstme, fname
      data frstme /.true./
c -------------------------------------------------------------------

      if (frstme) then
        sname = fdf_string( 'SystemLabel', 'siesta' )
        fname = paste( sname, '.KP' )
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
