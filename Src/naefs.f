C $Id: naefs.f,v 1.5 1999/01/31 11:20:08 emilio Exp $

      subroutine naefs(nua, na, scell, xa, indxua, rmaxv,
     .                maxna, isa,
     .                jna, xij, r2ij,
     .                Ena, fa, stress)
C *********************************************************************
C Neutral Atom (NA) energy, forces and stress.
C This is the self energy of rho_na=-Laplacian(v_na(Ry))/(8*pi)
C Energies in Ry. Lengths in Bohr.
C Writen by J.Soler and P.Ordejon, August-October'96, June'98 
C Modified by DSP Aug., 1998
C **************************** INPUT **********************************
C integer nua              : Number of atoms in unit cell
C integer na               : Number of atoms in supercell
C real*8  scell(3,3)       : Supercell vectors SCELL(IXYZ,IVECT)
C real*8  xa(3,na)         : Atomic positions in cartesian coordinates
c integer indxua(na)       : Index of equivalent atom in unit cell
C real*8  rmaxv            : Maximum cutoff for NA potential
C integer maxna            : Maximum number of neighbours of any atom
C integer isa(na)          : Species index of each atom
C integer jna(maxna)       : Aux. space to find neighbours (indexes)
C real*8  xij(3,maxna)     : Aux. space to find neighbours (vectors)
C real*8  r2ij(maxna)      : Aux. space to find neighbours (distances)
C **************************** OUTPUT *********************************
C real*8 Ena               : NA energy
C ********************** INPUT and OUTPUT *****************************
C real*8 fa(3,na)          : Atomic forces (NA part added to input)
C real*8 stress(3,3)       : Stress tensor (NA part added to input)
C *********************************************************************
C The following function must exits:
C
C    INTEGER  FUNCTION  IZOFIS(IS) : Returns the atomic number
C Input:
C      INTEGER IS : Specie index
C
C *********************************************************************

      implicit none

      integer
     .  maxna, na, nua

      integer
     . indxua(na), isa(na), jna(maxna)

      double precision
     . scell(3,3), Ena, fa(3,nua), r2ij(maxna), rmaxv, 
     . stress(3,3), xa(3,na), xij(3,maxna)

C Internal variables ......................................................
      integer
     .  ia, is, ix, ja, jn, js, jx, jua, nnia, izofis

      double precision
     .  fij(3), pi, vij, volcel, volume 
      
      external 
     .  izofis
C ......................

C     Initialize neighb subroutine
      nnia = maxna
      call neighb( scell, 2.d0*rmaxv, na, xa, 0, 0,
     .             nnia, jna, xij, r2ij )

      pi = 4.d0 * atan(1.d0)
      volume = nua * volcel(scell) / na
      Ena = 0.d0
 
      do ia = 1,nua
c       Find neighbour atoms
        nnia = maxna
        call neighb( scell, 2.d0*rmaxv, na, xa, ia, 0,
     .               nnia, jna, xij, r2ij )
        call chkdim( 'naefs', 'maxna', maxna, nnia, 1 )
        do jn = 1,nnia
          ja = jna(jn)
          jua = indxua(ja)
          is = isa(ia)
          js = isa(ja)
          if (izofis(is).gt.0 .and. izofis(js).gt.0) then
            call matel( 'T', is, js, 0, 0, xij(1,jn), vij, fij )
            Ena = Ena + vij / (16.d0*pi)
            do ix = 1,3
              fij(ix) = fij(ix) / (16.d0*pi)
              fa(ix,ia)  = fa(ix,ia)  + fij(ix)
              fa(ix,jua) = fa(ix,jua) - fij(ix)
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

