C $Id: dhscf.h,v 1.3 1999/01/31 11:56:51 emilio Exp $

C DHSCF: Dimension parameters for DHSCF
C INTEGER MAXA   : MAXimum number of Atoms
C INTEGER MAXAUX : MAXimum size of AUXiliary array
C INTEGER MAXEP  : MAXimum number of Extended-mesh Points
C INTEGER MAXMP  : MAXimum number of Mesh Points in unit cell
C INTEGER MAXO   : MAXimum number of Orbitals
C INTEGER MAXOP  : MAXimum number of non-zero Orbital Points
C INTEGER MAXPCC : MAXimum Partial Core Corrections (0 or 1)
C INTEGER MAXSPN : MAXimum number of different spin polarizations
C INTEGER MAXTOP : MAXimum number of non-zero Total Orbital Points
C INTEGER MAXTP  : MAXimum number of Total mesh Points
C INTEGER MAXTY  : MAXimum number of orbital TYpes

      INTEGER MAXA, MAXAUX, MAXEP, MAXMP, MAXO, MAXOP
      INTEGER MAXPCC, MAXSPN, MAXTOP, MAXTP, MAXTY
      PARAMETER ( MAXAUX = 1 )
      PARAMETER ( MAXA   = 1 )
      PARAMETER ( MAXEP  = 1 )
      PARAMETER ( MAXMP  = 1 )
      PARAMETER ( MAXO   = 1 )
      PARAMETER ( MAXOP  = 1 )
      PARAMETER ( MAXPCC = 1 )
      PARAMETER ( MAXSPN = 1 )
      PARAMETER ( MAXTOP = 1 )
      PARAMETER ( MAXTP  = 1 )
      PARAMETER ( MAXTY  = 1 )


