      subroutine overfsm(nua, na, no, scell, xa, indxua, rmaxo, maxo,
     .                   maxna, maxnh, maxnd, lasto, iphorb, isa, 
     .                   numd, listdptr, listd, numh, listhptr, listh, 
     .                   nspin, Emat, jna, xij, r2ij, fa, stress, S)
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
C integer maxnh            : First dimension of S and listh
C integer maxnd            : First dimension of Escf and listd
C integer lasto(0:na)      : Last orbital index of each atom
C integer iphorb(no)       : Orbital index of each orbital in its atom
C integer isa(na)          : Species index of each atom
C integer numd(nuotot)     : Number of nonzero elements of each row
C                            of Emat
C integer listdptr(nuotot) : Pointer to start of rows (-1) of Emat
C integer listd(maxnh)     : Column indexes of the nonzero elements  
C                            of each row of Emat
C integer numh(nuotot)     : Number of nonzero elements of each row
C                            of the overlap matrix
C integer listhptr(nuotot) : Pointer to start of rows (-1) of overlap
C                            matrix
C integer listh(maxnh)     : Column indexes of the nonzero elements  
C                            of each row of the overlap matrix
C integer nspin            : Number of spin components of Emat
C integer Emat(maxnd,nspin): Energy-Density matrix
C integer jna(maxna)       : Aux. space to find neighbours (indexes)
C real*8  xij(3,maxna)     : Aux. space to find neighbours (vectors)
C real*8  r2ij(maxna)      : Aux. space to find neighbours (distances)
C **************************** OUTPUT *********************************
C real*8  S(maxnh)         : Sparse overlap matrix
C ********************** INPUT and OUTPUT *****************************
C real*8  fa(3,nua)        : Atomic forces (Orthog. part added to input)
C real*8  stress(3,3)      : Stress tensor (Orthog. part added to input)
C *********************************************************************
C
C  Modules
C
      use precision
      use parallel,      only : Node, Nodes
      use parallelsubs,  only : GlobalToLocalOrb
      use atmfuncs,      only : rcut

      implicit none

      integer
     . maxna, maxnd, maxnh, maxo, na, no, nspin, nua

      integer
     . indxua(na), iphorb(no), isa(na), jna(maxna), lasto(0:na), 
     . listd(*), numd(*), listh(maxnh), numh(*), listdptr(*),
     . listhptr(*)

      double precision
     . scell(3,3), Emat(maxnd,nspin), 
     . fa(3,nua), r2ij(maxna), rmaxo, 
     . stress(3,3), S(maxnh), xa(3,na), xij(3,maxna)

C Internal variables ......................................................
  
      integer
     .  ia, ind, io, ioa, is, ispin, ix, iio, 
     .  j, ja, jn, jo, joa, js, jua, jx, nnia

      double precision
     .  fij(3), grSij(3) , rij, Sij, volcel, volume

      double precision, dimension(:), allocatable, save ::
     .  Di, Si

      external
     .  neighb, timer, memory
C ......................

C Start timer
      call timer( 'overfsm', 1 )

      volume = nua * volcel(scell) / na

C Initialize neighb subroutine 
      nnia = maxna
      call neighb( scell, 2.d0*rmaxo, na, xa, 0, 0,
     .             nnia, jna, xij, r2ij )

C Allocate local memory
      allocate(Di(no))
      call memory('A','D',no,'overfsm')
      allocate(Si(no))
      call memory('A','D',no,'overfsm')

      Si(1:no) = 0.0d0
      Di(1:no) = 0.0d0

      do ia = 1,nua
        nnia = maxna
        call neighb( scell, 2.d0*rmaxo, na, xa, ia, 0,
     .               nnia, jna, xij, r2ij )
        do io = lasto(ia-1)+1,lasto(ia)

C Is this orbital on this Node?
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          if (iio.gt.0) then

C Valid orbital 
            do j = 1,numd(iio)
              ind = listdptr(iio)+j
              jo = listd(ind)
              do ispin = 1,nspin  
                Di(jo) = Di(jo) + Emat(ind,ispin)
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
            do j = 1,numd(iio)
              jo = listd(listdptr(iio)+j)
              Di(jo) = 0.d0
            enddo
            do j = 1,numh(iio)
              ind = listhptr(iio)+j
              jo = listh(ind)
              S(ind) = Si(jo)
              Si(jo) = 0.0d0
            enddo
          endif
        enddo
      enddo

C Deallocate local memory
      call memory('D','D',size(Di),'overfsm')
      deallocate(Di)
      call memory('D','D',size(Si),'overfsm')
      deallocate(Si)

C Finish timer
      call timer( 'overfsm', 2 )

      return
      end




