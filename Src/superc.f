C $Id: superc.f,v 1.5 1999/11/26 18:28:35 wdpgaara Exp $

      subroutine superc( ucell, nua, nuo, nuokb,
     .                   maxa, maxo, maxkb,
     .                   xa, isa, iza, lasto, lastkb,
     .                   iaorb, iphorb, rco, iakb, iphkb, rckb,
     .                   scell, nsc, ncells,
     .                   na, no, nokb, indxua, indxuo )

C *********************************************************************
C Finds the supercell required to avoid multiple image overlaps,
C and expands arrays from unit cell to supercell
C Written by J.M.Soler. August 1998.
C *************** Input ***********************************************
C real*8  ucell(3,3) : Unit cell vectors
C integer nua        : Number of atoms in unit cell
C integer nuo        : Number of orbitals in unit cell
C integer nuokb      : Number of KB projectors in unit cell
C integer ncells     : Number of unit cells in supercell
C integer nsc(3)     : Diagonal elements of mscell
C integer maxa       : Dimension of several array arguments
C integer maxo       : Dimension of several array arguments
C integer maxkb      : Dimension of several array arguments
C ************** Input (unit cell) and Output (supercell) *************
C real*8  xa(2,maxa)     : Atomic coordinates
C integer isa(maxa)      : Species index of atoms
C integer iza(maxa)      : Atomic numbers
C integer lasto(0:maxa)  : Last orbital of each atom
C integer lastkb(0:maxa) : Last KB projector of each atom
C integer iaorb(maxo)    : Atom to which orbitals belong
C integer iphorb(maxo)   : Orbital index within atom
C real*8  rco(maxo)      : Orbital cutoff radiae
C integer iakb(maxo)     : Atom to which KB projectors belong
C integer iphKB(maxo)    : Index of KB projectors within atom
C real*8  rckb(maxo)     : Cutoff radiae of KB projectors
C *************** Output **********************************************
C real*8  scell(3,3)  : Supercell vectors
C integer na          : Number of atoms in supercell:    na=nua*ncells
C integer no          : Number of orbitals in supercell: no=nuo*ncells
C integer nokb        : Number of KB proj. in supercell: nokb=nuokb*ncells
C integer indxua(maxa): Index of equivalent atom in unit cell:
C                         indexua(ia) = mod(ia-1,nua)+1
C integer indxuo(maxo): Index of equivalent orbital in unit cell:
C                         indexuo(io) = mod(io-1,nuo)+1
C *********************************************************************

      implicit none
      integer
     .   maxa, maxkb, maxo, na, ncells, no, nokb, nua, nuo, nuokb
      integer
     .   iakb(maxkb), iaorb(maxo), indxua(maxa), indxuo(maxo),
     .   iphKB(maxkb), iphorb(maxo), isa(maxa), iza(maxa),
     .   lastkb(0:maxa), lasto(0:maxa), nsc(3)
      double precision
     .   rckb(maxkb), rco(maxo), ucell(3,3), scell(3,3),
     .   xa(3,maxa)
      external superx

C Internal variables
      integer           ia, io, iua, iuo, ja

C Find number of cells, atoms and orbitals in supercell
      ncells = nsc(1) * nsc(2) * nsc(3)
      na    = nua   * ncells
      no    = nuo   * ncells
      nokb  = nuokb * ncells

C Find supercell vectors and atomic coordinates in supercell 
      call superx( ucell, nsc, nua, maxa, xa, scell )

C Find indxua and expand isa, iza, lasto and lastkb to supercell 
      do ia = 1,na
        ja = mod(ia-1,nua) + 1
        indxua(ia) = ja
        isa(ia)    = isa(ja)
        iza(ia)    = iza(ja)
        lasto(ia)  = lasto(ia-1)  + lasto(ja)  - lasto(ja-1)
        lastkb(ia) = lastkb(ia-1) + lastkb(ja) - lastkb(ja-1)
      enddo

C Find indxuo and expand iaorb, iphorb, and rco 
      do io = 1,no
        indxuo(io) = mod(io-1,nuo) + 1
      enddo
      do ia = 1,na
        do io = lasto(ia-1)+1,lasto(ia)
          iuo = indxuo(io)
          iaorb(io)  = ia
          iphorb(io) = iphorb(iuo)
          rco(io)    = rco(iuo)
        enddo
      enddo

C Expand iakb and iphKB 
      do ia = 1,na
        iua = indxua(ia)
        iuo = lastkb(iua-1)
        do io = lastkb(ia-1)+1,lastkb(ia)
          iuo = iuo + 1
          iakb(io)  = ia
          iphKB(io) = iphKB(iuo)
          rckb(io)  = rckb(iuo)
        enddo
      enddo

      return
      end

