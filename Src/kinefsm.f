      subroutine kinefsm(na, no, cell, xa, rmaxo, maxo,
     .                  maxna, maxno, maxnd, lasto, iphorb, isa, volume,
     .                  numd, listd, numh, listh, nspin, Dscf, 
     .                  jna, xij, r2ij,
     .                  Ekin, fa, stress, H)
C *********************************************************************
C Routine to calculate Kinetic energies, forces, stress and matrix
C elements.
C (Energies in Ry.; distances in Bohr)
C
C Writen by J.Soler and P.Ordejon, August-October'96
C **************************** INPUT **********************************
C integer na               : Total number of atoms
C integer no               : Total number of orbitals
C real*8 cell(3,3)         : Unit cell vectors CELL(IXYZ,IVECT)
C real*8 xa(3,na)          : Atomic positions in cartesian coordinates
C real*8 rmaxo             : Maximum cutoff for atomic orbitals
C integer maxo             : Maximum number of atomic orbitals
C integer maxna            : Maximum number of neighbours of any atom
C integer maxno            : Maximum number of basis orbitals interacting
C                            either directly or thru a KB projector
C integer maxnd            : Maximum number of nonzero elements of each row
C                            of density matrix
C integer lasto(0:na)      : Position of last orbital of each atom
C integer iphorb(no)       : Orbital index of each orbital in its atom
C integer isa(na)          : Species index of each atom
C real*8 volume            : Cell volume 
C integer numd(no)         : Number of nonzero elements of each row of the
C                            density matrix
C integer listd(maxnd,no)  : Nonzero density-matrix element column 
C                            indexes for each matrix row
C integer numh(no)         : Number of nonzero elements of each row of the
C                            hamiltonian matrix
C integer listh(maxno,no)  : Nonzero hamiltonian-matrix element column 
C                            indexes for each matrix row
C integer nspin            : Spin polarization (1 or 2)
C integer Dscf(maxnd,maxo,nspin) : Density matrix
C integer jna(maxna)       : Aux. space for overlaping neighbours (indexes)
C real*8 xij(3,maxna)      : Aux. space for overlaping neighbours (vectors)
C real*8 r2ij(maxna)       : Aux. space for overlaping neighbours (distances)
C **************************** OUTPUT *********************************
C real*8 Ekin              : Kinetic energy (Ry)
C real*8 fa(3,na)          : NA forces (Ry/Bohr) (added to input fa)
C real*8 stress(3,3)       : NA stress (Ry/Bohr) (added to input stress)
C real*8 H(maxno,maxo,nspin) : Hamiltonian with Kinetic matrix elements added.
C *********************************************************************
      implicit none

      integer
     . maxna, maxnd, maxno, maxo, na, no, nspin

      integer
     . iphorb(no), isa(na), jna(maxna), lasto(0:na), 
     . listd(maxnd,no), listh(maxno,no), numd(no), numh(no)

      double precision
     . cell(3,3), Dscf(maxnd,maxo,nspin), Ekin, 
     . fa(3,na), H(maxno,maxo, nspin), r2ij(maxna), rmaxo, 
     . stress(3,3),  xa(3,na), xij(3,maxna), volume

C Internal variables ......................................................
C nomax = maximum number of orbitals
      integer nomax
      parameter (nomax = 10000)
  
      integer
     .  ia, io, ioa, is, isel, ispin, ix, 
     .  j, ja, jn, jo, joa, js, jx, nnia

      double precision
     .  Di(nomax), fij(3), grTij(3) , rcut, rij, Ti(nomax), Tij

      external rcut

      common/DV/Di,Ti
C ......................



      Ekin = 0.d0
      do jo = 1,no
        Di(jo) = 0.d0
        Ti(jo) = 0.d0
      enddo

      isel = 0

      do ia = 1,na
        nnia = maxna
        call neighb( cell, 2.d0*rmaxo, na, xa, ia, isel,
     .               nnia, jna, xij, r2ij )
        call chkdim( 'kb', 'maxna', maxna, nnia, 1 )
        do io = lasto(ia-1)+1,lasto(ia)
          do j = 1,numd(io)
            jo = listd(j,io)
            do ispin = 1,nspin
              Di(jo) = Di(jo) + Dscf(j,io,ispin)
            enddo
          enddo
          do jn = 1,nnia
            ja = jna(jn)
            rij = sqrt( r2ij(jn) )
            do jo = lasto(ja-1)+1,lasto(ja)
              ioa = iphorb(io)
              joa = iphorb(jo)
              is = isa(ia)
              js = isa(ja)
              if (rcut(is,ioa)+rcut(js,joa) .gt. rij) then
                call matel( 'T', is, js, ioa, joa, xij(1,jn),
     .                      Tij, grTij )
                Ti(jo) = Ti(jo) + Tij
                Ekin = Ekin + Di(jo) * Tij
                do ix = 1,3
                  fij(ix) = Di(jo) * grTij(ix)
                  fa(ix,ia) = fa(ix,ia) + fij(ix)
                  fa(ix,ja) = fa(ix,ja) - fij(ix)
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
            do ispin = 1,nspin
              H(j,io,ispin) = H(j,io,ispin) + Ti(jo)
            enddo
            Ti(jo) = 0.d0
          enddo
        enddo
      enddo

      return
      end
