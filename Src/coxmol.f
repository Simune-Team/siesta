c $Id: coxmol.f,v 1.3 1999/01/31 10:50:54 emilio Exp $

      subroutine coxmol(iza, xa, na, slabel)

c *******************************************************************
c Writes coordinates in format to be read by Xmol
c Written by E. Artacho. December 1997.
c ********* INPUT ***************************************************
c integer   iza(na)   : Atomic numbers of different atoms
c double    xa(3,na)  : Atom coordinates
c integer   na        : Number of atoms
c character slabel*20 : Label for file naming
c ******************************************************************

      implicit          none
      character         slabel*20, paste*24
      integer           na
      integer           iza(na)
      double precision  xa(3,na)
      external          io_assign, io_close, paste, symbol

c Internal variables and arrays
 
      character         fname*24, symbol*2
      integer           unit, i, ia
      double precision  Ang

      Ang  = 1.d0 / 0.529177d0

c Find file name

      fname = paste(slabel,'.xyz')

      write(6,'(/,2a)')'coxmol: Writing XMOL coordinates into file ',
     .                  fname

      call io_assign(unit)
      open( unit, file=fname, form = 'formatted', status='unknown')
      rewind(unit)

      write(unit,'(i5)') na
      write(unit,*)
      write(unit,'(a2,2x,3f12.6)')
     .     ( symbol(iza(ia)), (xa(i,ia)/Ang,i=1,3), ia=1,na )

      call io_close(unit)
      
      return
      end

