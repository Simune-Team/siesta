c Dimension parameters for siesta:
c  maxa   : Maximum number of atoms
c  maxbk  : Maximum number of band k points
c  maxk   : Maximum number of kpoints
c  maxkb  : Maximum total number of Kleinman-Bylander projectors
c  maxkba : Maximum number of KB projectors of one atom
c  maxl   : Maximum angular momentum
c  maxna  : Maximum number of neighbours of any atom
c  maxno  : Max basis orbitals connected to any given one,
c           either directly or through a KB projector
c  maxo   : Maximum total number of basis orbitals
c  maxos  : Maximum number of basis orbitals of any atom
c  maxs   : Maximum number of species
c  maxspn : Maximum number of spin polarizations (1 or 2)
c  maxuo  : Maximum number of basis orbitals in unit cell
c  maxzet : Maximum number of orbital zetas
c  dimaux : Dimension of Pulay auxiliary vectors

      integer maxa, maxbk, maxk, maxkb, maxkba, maxl
      integer maxna, maxno, maxo, maxos
      integer maxs, maxspn, maxuo, maxzet, dimaux
      parameter ( maxa   =  1000 )
      parameter ( maxbk  =   100 )
      parameter ( maxk   =   100 )
      parameter ( maxkb  =  9000 )
      parameter ( maxkba =     1 )
      parameter ( maxl   =     3 )
      parameter ( maxna  =   500 )
      parameter ( maxno  =     1 )
      parameter ( maxo   = 13000 )
      parameter ( maxos  =    60 )
      parameter ( maxs   =    10 )
      parameter ( maxspn =     2 )
      parameter ( maxuo  =   500 )
      parameter ( maxzet =     3 )
      parameter ( dimaux =     1 )

