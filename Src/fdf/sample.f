      program sample
c
c     Shows fdf capabilities..
c
c $Id: sample.f,v 1.6 1997/04/02 17:31:57 wdpgaara Exp $
c-------
c $Log: sample.f,v $
c Revision 1.6  1997/04/02 17:31:57  wdpgaara
c Modifications to support change to fdf_block. Minor cosmetic
c improvements.
c
c Revision 1.5  1997/03/20 17:38:26  wdpgaara
c Support for declarations file.
c
c Revision 1.4  1997/03/03 08:33:23  wdpgaara
c Fixed a bug in fdf_physical (the default value was always used...)
c fdf_setdebug is now called from fdf_init: it looks for fdf-debug in the
c fdf file itsel.
c Added fdf_enabled and fdf_shutdown.
c Removed fdf_unit.
c Removed some unused variables.
c
c Revision 1.3  1997/02/03 15:40:30  wdpgaara
c Version 0.6 of fdf package.
c
c Revision 1.2  1997/01/03 19:44:10  wdpgaara
c Adapted to the new version of fdf.
c
c-------

      implicit none
      integer  maxa
      parameter ( maxa = 100 )

      character         fname*20, symbol(maxa)*2
      integer           i, ia, isa(maxa), na, na_default, iblk
      real              wmix
      double precision  factor, xa(3,maxa), cutoff, phonon_energy

      logical doit, debug

      include 'fdfdefs.h'
c
      call fdf_init('sample.fdf','sample.out')

      na_default = 0
      na = fdf_integer('NumberOfAtoms', na_default )
      write(6,*) 'examples: na =', na

      fname = fdf_string('NameOfFile','calimero')
      write(6,*) fname

      cutoff = fdf_physical('MeshCutoff',8.d0,'Ry')
      write(6,*) cutoff

      phonon_energy = fdf_physical('phonon-energy',0.01d0,'eV')
      write(6,*) phonon_energy

      i = fdf_integer('SomeInt',34)
      write(6,*) i

      wmix = fdf_single('WmixValue',0.55)
      write(6,*) wmix

      factor = fdf_double('FactorValue',1.d-10)
      write(6,*) factor

      debug = fdf_boolean('Debug',.true.)
      write(6,*) debug

      doit = fdf_boolean('DoIt',.false.)
      write(6,*) doit
      
      if (fdf_block('AtomicCoordinatesAndAtomicSpecies',iblk)) then
        do ia = 1,na
          read(iblk,*) (xa(i,ia),i=1,3), isa(ia)
        enddo
      endif
      
      if (fdf_block('AtomicSymbolsAndAtomicCoordinates',iblk)) then
         do ia = 1,na
          read(iblk,*) symbol(ia), (xa(i,ia),i=1,3)
        enddo
      endif

      do ia = 1,na
         write(6,*) (xa(i,ia),i=1,3)
      enddo

      if (fdf_block('AtomicInfo',iblk)) then
        do ia = 1,na
          read(iblk,*) (xa(i,ia),i=1,3)
        enddo
      endif

      do ia = 1,na
         write(6,*) (xa(i,ia),i=1,3)
      enddo

      end






