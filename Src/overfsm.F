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
C real*8  S(maxnh,nuo)     : Sparse overlap matrix
C ********************** INPUT and OUTPUT *****************************
C real*8  fa(3,nua)        : Atomic forces (Orthog. part added to input)
C real*8  stress(3,3)      : Stress tensor (Orthog. part added to input)
C *********************************************************************
C
C  Modules
C
      use precision
      use parallel
      use atmfuncs, only: rcut
      use matel_module, only: matel
#ifdef MPI
      use mpi
#endif

      implicit none

      integer
     . maxna, maxnd, maxnh, maxo, na, no, nspin, nua

      integer
     . indxua(na), iphorb(no), isa(na), jna(maxna), lasto(0:na), 
     . listd(maxnd,no), numd(no), listh(maxnh,no), numh(no)

      double precision
     . scell(3,3), Emat(maxnd,maxo,nspin), 
     . fa(3,nua), r2ij(maxna), rmaxo, 
     . stress(3,3), S(maxnh,*), xa(3,na), xij(3,maxna)

C Internal variables ......................................................
  
      integer
     .  ia, io, ioa, is, ispin, ix, iio, 
     .  j, ja, jn, jo, joa, js, jua, jx, nnia

      integer
     .  Node, Nodes

#ifdef MPI
      integer
     .  MPIerror
#endif

      double precision
     .  fij(3), grSij(3) , rij, Sij, volcel, volume

      double precision, dimension(:), allocatable ::
     .  Di, Si

      external
     .  chkdim, neighb, timer, memory
C ......................

C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
#else
      Node = 0
      Nodes = 1
#endif

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

C Is this orbital on this Node?
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          if (iio.gt.0) then

C Valid orbital 
            do j = 1,numd(io)
              jo = listd(j,io)
              do ispin = 1,nspin  
                Di(jo) = Di(jo) + Emat(j,iio,ispin)
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
              S(j,iio) = Si(jo)
              Si(jo) = 0.d0
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




