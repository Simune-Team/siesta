c $Id: ioeig.f,v 1.2 1999/11/26 18:28:16 wdpgaara Exp $

      subroutine ioeig(eo, ef, no, ns, nk, maxo, maxs, maxk)

c *******************************************************************
c Writes eigenvalues of Hamiltonian in k-points of sampling
c Emilio Artacho, Feb. 1999
c ********** INPUT **************************************************
c integer nk           : Number k points
c real*8  points(3,nk) : k-point coordinates
c real*8  weight(3,nk) : k-point weight
c *******************************************************************

      use fdf

      implicit          none
      character         paste*33
      integer           maxo, maxs, maxk, no, ns, nk
      double precision  eo(maxo, maxs, maxk), ef
      external          io_assign, io_close, paste

c Internal 
      character         sname*30, fname*33
      integer           ik, iu, io, is
      logical           frstme
      double precision  eV
      save              frstme, fname, eV
      data              frstme /.true./
c -------------------------------------------------------------------

      if (frstme) then
        sname = fdf_string( 'SystemLabel', 'siesta' )
        fname = paste( sname, '.EIG' )
        eV = 1.d0 / 13.60580d0
        frstme = .false.
      endif

      call io_assign( iu )
      open( iu, file=fname, form='formatted', status='unknown' )      

      write(iu,"(f14.4)") ef/eV
      write(iu,"(3i6)")   no, min(ns,2), nk
      do ik = 1,nk
        write(iu,"(i5,10f12.5,/,(5x,10f12.5))")
     .          ik, ((eo(io,is,ik)/eV,io=1,no),is=1,min(ns,2))
      enddo

      call io_close( iu )

      return
      end
