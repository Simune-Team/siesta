c $Id: pixmol.f,v 1.2 1999/04/13 11:20:49 emilio Exp $

      subroutine pixmol(iza, xa, na, slabel, last)

c *******************************************************************
c Writes and accumulates coordinates to be animated by Xmol
c Written by E. Artacho. February 1999.
c ********* INPUT ***************************************************
c integer   iza(na)   : Atomic numbers of different atoms
c double    xa(3,na)  : Atom coordinates
c integer   na        : Number of atoms
c character slabel*20 : Label for file naming
c logical   last      : true if last time step
c *******************************************************************

      implicit          none
      character         slabel*20, paste*24
      integer           na
      integer           iza(na)
      double precision  xa(3,na)
      logical           last
      external          io_assign, io_close, paste, symbol

c Internal variables and arrays
 
      character         fname*24, symbol*2
      logical           frstme
      integer           unit, i, ia
      double precision  Ang
      save              frstme, Ang, unit
      data              frstme /.true./
c -------------------------------------------------------------------

      if ( frstme ) then
        fname = paste(slabel,'.ANI')
        call io_assign(unit)
        open( unit, file=fname, form = 'formatted', status='unknown')
        call windf(unit)
        Ang  = 1.d0 / 0.529177d0
        frstme = .false.
      endif

      write(unit,'(i5)') na
      write(unit,*)
      write(unit,'(a2,2x,3f12.6)')
     .     ( symbol(iza(ia)), (xa(i,ia)/Ang,i=1,3), ia=1,na )

      if (last) call io_close(unit)
      
      return
      end

