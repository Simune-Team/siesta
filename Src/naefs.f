      subroutine naefs(na, ns, cell, xa, rmaxv,
     .                maxna, isa, izs, volume,
     .                jna, xij, r2ij,
     .                Ena, fa, stress)
C *********************************************************************
C Routine to calculate Neutral Atom energies, forces and stress.
C This is the self energy of rho_na=-Laplacian(v_na(Ry))/(8*pi)
C (Energies in Ry.; distances in Bohr)
C
C Writen by J.Soler and P.Ordejon, August-October'96
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
C real*8 Ena               : NA energy (Ry)
C real*8 fa(3,na)          : NA forces (Ry/Bohr) (added to input fa)
C real*8 stress(3,3)      : NA stress (Ry/Bohr) (added to input stress)
C *********************************************************************
      implicit none

      integer
     . maxna, na, ns

      integer
     . isa(na), izs(ns), jna(maxna)

      double precision
     . cell(3,3), Ena, fa(3,na), r2ij(maxna), rmaxv, 
     . stress(3,3), xa(3,na), xij(3,maxna), volume

C Internal variables ......................................................
      integer
     .  ia, is, isel, ix, ja, jn, js, jx, nnia

      double precision
     .  fij(3), pi, vij
C ......................


      pi = 4.d0 * atan(1.d0)
      Ena = 0.d0
      isel = 0
 
      do ia = 1,na
c       Find neighbour atoms
        nnia = maxna
        call neighb( cell, 2.d0*rmaxv, na, xa, ia, isel,
     .               nnia, jna, xij, r2ij )
        call chkdim( 'naef', 'maxna', maxna, nnia, 1 )
        do jn = 1,nnia
          ja = jna(jn)
          is = isa(ia)
          js = isa(ja)
          if (izs(is).gt.0 .and. izs(js).gt.0) then
            call matel( 'T', is, js, 0, 0, xij(1,jn), vij, fij )
            Ena = Ena + vij / (16.d0*pi)
            do ix = 1,3
              fij(ix) = fij(ix) / (16.d0*pi)
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

