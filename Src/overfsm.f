C $Id: overfsm.f,v 1.7 1999/05/05 17:25:35 emilio Exp $

      subroutine overfsm(nua, na, no, scell, xa, indxua, rmaxo, maxo,
     .                   maxna, maxnh, maxnd, lasto, iphorb, isa, 
     .                   numd, listd, numh, listh, nspin, Emat, 
     .                   jna, xij, r2ij,
     .                   fa, stress, S)
C *********************************************************************
C Overlap matrix and orthonormalization contribution to forces and stress.
C Energies in Ry. Lengths in Bohr.
C Writen by J.Soler and P.Ordejon, August-October'96, June'98
C **************************** INPUT **********************************
C integer nua              : Number of atoms in unit cell
C integer na               : Number of atoms in supercell
C integer no               : Number of orbitals in supercell
C real*8  scell(3,3)       : Supercell vectors SCELL(IXYZ,IVECT)
C real*8  xa(3,na)         : Atomic positions in cartesian coordinates
c integer indxua(na)       : Index of equivalent atom in unit cell
C real*8  rmaxo            : Maximum cutoff for atomic orbitals
C integer maxo             : Second dimension of Emat
C integer maxna            : Maximum number of neighbours of any atom
C integer maxnh            : First dimension of H and listh
C integer maxnd            : First dimension of Escf and listd
C integer lasto(0:na)      : Last orbital index of each atom
C integer iphorb(no)       : Orbital index of each orbital in its atom
C integer isa(na)          : Species index of each atom
C integer numd(no)         : Number of nonzero elements of each row
C                            of Emat
C integer listd(maxnd,no)  : Colum indexes of the nonzero elements  
C                            of each row of Emat
C integer numh(no)         : Number of nonzero elements of each row
C                            of the overlap matrix
C integer listh(maxnh,no)  : Colum indexes of the nonzero elements  
C                            of each row of the overlap matrix
C integer nspin            : Number of spin components of Emat
C integer Emat(maxnd,maxo,nspin) : Energy-Density matrix
C integer jna(maxna)       : Aux. space to find neighbours (indexes)
C real*8  xij(3,maxna)     : Aux. space to find neighbours (vectors)
C real*8  r2ij(maxna)      : Aux. space to find neighbours (distances)
C **************************** OUTPUT *********************************
C real*8  S(maxnh,no)      : Sparse overlap matrix
C ********************** INPUT and OUTPUT *****************************
C real*8  fa(3,nua)        : Atomic forces (Orthog. part added to input)
C real*8  stress(3,3)      : Stress tensor (Orthog. part added to input)
C *********************************************************************
      implicit none

      integer
     . maxna, maxnd, maxnh, maxo, na, no, nspin, nua

      integer
     . indxua(na), iphorb(no), isa(na), jna(maxna), lasto(0:na), 
     . listd(maxnd,no), numd(no), listh(maxnh,no), numh(no)

      double precision
     . scell(3,3), Emat(maxnd,maxo,nspin), 
     . fa(3,na), r2ij(maxna), rmaxo, 
     . stress(3,3), S(maxnh,no), xa(3,na), xij(3,maxna)

C Internal variables ......................................................
C nomax = maximum number of orbitals
      integer nomax
      parameter (nomax = 20000)
  
      integer
     .  ia, io, ioa, is, ispin, ix, 
     .  j, ja, jn, jo, joa, js, jua, jx, nnia

      double precision
     .  Di(nomax), fij(3), grSij(3) , rcut, rij, Si(nomax), Sij,
     .  volcel, volume

      external
     .  chkdim, matel, neighb, rcut
C ......................

      volume = nua * volcel(scell) / na

C     Initialize neighb subroutine 
      nnia = maxna
      call neighb( scell, 2.d0*rmaxo, na, xa, 0, 0,
     .             nnia, jna, xij, r2ij )

      call chkdim( 'overfsm', 'nomax', nomax, no, 1 )
      do jo = 1,no
        Si(jo) = 0.d0
        Di(jo) = 0.d0
      enddo

      do ia = 1,nua
        nnia = maxna
        call neighb( scell, 2.d0*rmaxo, na, xa, ia, 0,
     .               nnia, jna, xij, r2ij )
        call chkdim( 'overfsm', 'maxna', maxna, nnia, 1 )
        do io = lasto(ia-1)+1,lasto(ia)
          do j = 1,numd(io)
            jo = listd(j,io)
            do ispin = 1,nspin  
              Di(jo) = Di(jo) + Emat(j,io,ispin)
            enddo
          enddo
          do jn = 1,nnia
            ja = jna(jn)
            jua = indxua(ja)
            rij = sqrt( r2ij(jn) )
            do jo = lasto(ja-1)+1,lasto(ja)
              ioa = iphorb(io)
              joa = iphorb(jo)
              is = isa(ia)
              js = isa(ja)
              if (rcut(is,ioa)+rcut(js,joa) .gt. rij) then
                call matel( 'S', is, js, ioa, joa, xij(1,jn),
     .                      Sij, grSij )
                Si(jo) = Si(jo) + Sij
                do ix = 1,3
                  fij(ix) = (- Di(jo)) * grSij(ix)
                  fa(ix,ia)  = fa(ix,ia)  + fij(ix)
                  fa(ix,jua) = fa(ix,jua) - fij(ix)
                  do jx = 1,3
                    stress(jx,ix) = stress(jx,ix) +
     .                              xij(jx,jn) * fij(ix) / volume
                  enddo
                enddo
              endif
            enddo
          enddo
          do j = 1,numd(io)
            jo = listd(j,io)
            Di(jo) = 0.d0
          enddo
          do j = 1,numh(io)
            jo = listh(j,io)
            S(j,io) = Si(jo)
            Si(jo) = 0.d0
          enddo
        enddo
      enddo

      return
      end




