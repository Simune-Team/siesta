      subroutine dnaefs(na, ns, cell, xa, rmaxv,
     .                maxna, isa, izs, volume,
     .                jna, xij, r2ij,
     .                DEna, fa, stress)
C *********************************************************************
C Correction of Neutral Atom energies, forces and stress due to the
C overlap between ionic (bare pseudopotential) charges.
C (Energies in Ry.; distances in Bohr)
C
C Writen by J.Soler and P.Ordejon, March'97
C **************************** INPUT **********************************
C integer na               : Total number of atoms
C integer ns               : Total number of species
C real*8 cell(3,3)         : Unit cell vectors CELL(IXYZ,IVECT)
C real*8 xa(3,na)          : Atomic positions in cartesian coordinates
C real*8 rmaxv             : Maximum cutoff for NA potential
C integer maxna            : Maximum number of neighbours of any atom
C integer isa(na)          : Species index of each atom
C integer izs(ns)          : Atomic number of each species
C real*8 volume            : Cell volume 
C integer jna(maxna)       : Aux. space for overlaping neighbours (indexes)
C real*8 xij(3,maxna)      : Aux. space for overlaping neighbours (vectors)
C real*8 r2ij(maxna)       : Aux. space for overlaping neighbours (distances)
C **************************** OUTPUT *********************************
C real*8 DEna              : NA energy correction (Ry)
C real*8 fa(3,na)          : NA forces (Ry/Bohr) (added to input fa)
C real*8 stress(3,3)       : NA stress (Ry/Bohr) (added to input stress)
C *********************************************************************
      implicit none

      integer
     . maxna, na, ns

      integer
     . isa(na), izs(ns), jna(maxna)

      double precision
     . cell(3,3), DEna, fa(3,na), r2ij(maxna), rmaxv, 
     . stress(3,3), xa(3,na), xij(3,maxna), volume

C Internal variables ......................................................
      integer
     .  ia, is, isel, ix, ja, jn, js, jx, nnia

      double precision
     .  dvdr, fij(3), rij, r2min, vij

      parameter ( r2min = 1.d-15 )
C ......................


      DEna = 0.d0
      isel = 0
 
      do ia = 1,na
c       Find neighbour atoms
        nnia = maxna
        call neighb( cell, 2.d0*rmaxv, na, xa, ia, isel,
     .               nnia, jna, xij, r2ij )
        call chkdim( 'dnaefs', 'maxna', maxna, nnia, 1 )
        do jn = 1,nnia
          ja = jna(jn)
          is = isa(ia)
          js = isa(ja)
          if (r2ij(jn).gt.r2min .and.
     .        izs(is).gt.0 .and. izs(js).gt.0) then
            rij = sqrt( r2ij(jn) )
            call psover( is, js, rij, vij, dvdr )
            DEna = DEna + vij / 2.d0
            do ix = 1,3
              fij(ix) = dvdr * xij(ix,jn) / rij / 2.d0
              fa(ix,ia) = fa(ix,ia) + fij(ix)
              fa(ix,ja) = fa(ix,ja) - fij(ix)
              do jx = 1,3
                stress(jx,ix) = stress(jx,ix) +
     .                          xij(jx,jn) * fij(ix) / volume
              enddo
            enddo
          endif
        enddo
      enddo

      return
      end

