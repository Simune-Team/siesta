c $Id: iofa.f,v 1.1 1999/02/21 00:09:49 emilio Exp $

      subroutine iofa( na, fa )

c *******************************************************************
c Writes forces in eV/Ang
c Emilio Artacho, Feb. 1999
c ********** INPUT **************************************************
c integer na           : Number atoms
c real*8  fa(3,na)     : Forces on the atoms
c *******************************************************************

      implicit          none
      character         paste*33
      integer           na
      double precision  fa(3,*)
      external          io_assign, io_close, paste
      include          'fdf/fdfdefs.h'

c Internal 
      character         sname*30, fname*33
      integer           ia, iu, ix
      logical           frstme
      double precision  Ang, eV
      save              frstme, fname, eV, Ang
      data frstme        /.true./
c -------------------------------------------------------------------

      if (frstme) then
        Ang    = 1.d0 / 0.529177d0
        eV     = 1.d0 / 13.60580d0
        sname  = fdf_string( 'SystemLabel', 'siesta' )
        fname  = paste( sname, '.FA' )
        frstme = .false.
      endif

      call io_assign( iu )
      open( iu, file=fname, form='formatted', status='unknown' )      

      write(iu,'(i6)') na
      write(iu,'(i6,3f12.6)') (ia, (fa(ix,ia)*Ang/eV,ix=1,3), ia=1,na)

      call io_close( iu )

      return
      end
