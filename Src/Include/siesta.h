c $Id: siesta.h,v 1.7 1999/01/31 11:56:53 emilio Exp $

c Dimension parameters for siesta:
c  maxa   : Maximum total number of atoms
c  maxbk  : Maximum number of band k points
c  maxk   : Maximum number of kpoints
c  maxkb  : Maximum total number of Kleinman-Bylander projectors
c  maxna  : Maximum number of neighbours of any atom
c  maxno  : Max basis orbitals connected to any given one,
c           either directly or through a KB projector
c  maxo   : Maximum total number of basis orbitals
c  maxpul : Dimension of Pulay auxiliary vectors
c  maxspn : Maximum number of spin polarizations (1 or 2)
c  maxua  : Maximum number of atoms in unit cell
c  maxuo  : Maximum number of basis orbitals in unit cell
c  maxxij : Dimension of array xijo

      integer maxa, maxbk, maxk, maxkb
      integer maxna, maxno, maxo, maxpul
      integer maxspn, maxua, maxuo
      integer maxxij
      parameter ( maxa   =  1000 )
      parameter ( maxbk  =   100 )
      parameter ( maxk   =   100 )
      parameter ( maxkb  =  9000 )
      parameter ( maxna  =   500 )
      parameter ( maxno  =     1 )
      parameter ( maxo   = 13000 )
      parameter ( maxpul =     1 )
      parameter ( maxspn =     2 )
      parameter ( maxua  =  1000 )
      parameter ( maxuo  =   500 )
      parameter ( maxxij =     1 )

