c $Id: prversion.f,v 1.3 1999/01/31 11:20:15 emilio Exp $

      subroutine prversion
c
c     Simple routine to print the version string. Could be extended to
c     provide more information, if needed.
c
c     Alberto Garcia, Feb. 23, 1998
c
      implicit none
      include 'version.h'

c     write(6,'(2a)') 'prversion:    ', version
      write(6,'(2a)') '              ', version

      end
