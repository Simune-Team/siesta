      subroutine diagg( nspin, nuo, no, maxuo, maxnh, maxnd, 
     .                  maxo, numh, listhptr, listh, numd, 
     .                  listdptr, listd, H, S,
     .                  getD, getPSI, fixspin, qtot, qs, temp, e1, e2,
     .                  eo, qo, Dnew, Enew, ef, efs, Entropy,
     .                  Haux, Saux, psi, aux, nuotot )
C *********************************************************************
C Subroutine to calculate the eigenvalues and eigenvectors, density
C and energy-density matrices, and occupation weights of each 
C eigenvector, for given Hamiltonian and Overlap matrices (including
C spin polarization). Gamma-point version.
C Writen by J.Soler, August 1998.
C **************************** INPUT **********************************
C integer nspin               : Number of spin components (1 or 2)
C integer nuo                 : Number of basis orbitals local to node
C integer no                  : Number of basis orbitals
C integer maxuo               : Last dimension of xij
C                               Must be at least max(indxuo)
C integer maxnh               : Maximum number of orbitals interacting  
C integer maxnd               : Maximum number of nonzero elements of 
C                               each row of density matrix
C integer maxo                : First dimension of eo and qo
C integer numh(nuo)           : Number of nonzero elements of each row 
C                               of hamiltonian matrix
C integer listhptr(nuo)       : Pointer to each row (-1) of the
C                               hamiltonian matrix
C integer listh(maxnh)        : Nonzero hamiltonian-matrix element  
C                               column indexes for each matrix row
C integer numd(nuo)           : Number of nonzero elements of each row 
C                               ofdensity matrix
C integer listdptr(nuo)       : Pointer to each row (-1) of the
C                               density matrix
C integer listd(maxnd)        : Nonzero density-matrix element column 
C                               indexes for each matrix row
C real*8  H(maxnh,nspin)      : Hamiltonian in sparse form
C real*8  S(maxnh)            : Overlap in sparse form
C logical getD                : Find occupations and density matrices?
C logical getPSI              : Find and print wavefunctions?
C logical fixspin             : Fix the spin of the system?
C real*8  qtot                : Number of electrons in unit cell
C real*8  qs(nspin)           : Number of electrons in unit cell for each
C                               spin component (if fixed spin option is used)
C real*8  temp                : Electronic temperature 
C real*8  e1, e2              : Energy range for density-matrix states
C                               (to find local density of states)
C                               Not used if e1 > e2
C integer nuotot              : total number of orbitals per unit cell
C                               over all processors
C *************************** OUTPUT **********************************
C real*8 eo(maxo,nspn)        : Eigenvalues
C ******************** OUTPUT (only if getD=.true.) *******************
C real*8 qo(maxo,nspn)        : Occupations of eigenstates
C real*8 Dnew(maxnd,nspin)    : Output Density Matrix
C real*8 Enew(maxnd,nspin)    : Output Energy-Density Matrix
C real*8 ef                   : Fermi energy
C real*8 efs(nspin)           : Fermi energy for each spin
C                               (for fixed spin calculations)
C real*8 Entropy              : Electronic entropy
C *************************** AUXILIARY *******************************
C real*8 Haux(nuotot,nuo)     : Auxiliary space for the hamiltonian matrix
C real*8 Saux(nuotot,nuo)     : Auxiliary space for the overlap matrix
C real*8 psi(nuotot,maxuo,nspin) : Auxiliary space for the eigenvectors
C real*8 aux(nuotot*5)        : Extra auxiliary space
C *************************** UNITS ***********************************
C eo, Enew and ef returned in the units of H.
C *************************** PARALLEL ********************************
C The auxiliary arrays are now no longer symmetry and so the order
C of referencing has been changed in several places to reflect this.
C *********************************************************************
C
C  Modules
C
      use precision
      use sys
      use parallel,      only : Node, Nodes, BlockSize
      use parallelsubs,  only : LocalToGlobalOrb
      use writewave
      use m_fermid,      only : fermid, fermispin, stepf

#ifdef MPI
      use mpi_siesta
#endif

      implicit none

#ifdef MPI
      integer 
     .  MPIerror
#endif

      integer
     .  maxnd, maxnh, maxuo, maxo, nuo, no, nspin, nuotot

      integer 
     .  listh(maxnh), numh(nuo), listhptr(nuo),
     .  listd(maxnd), numd(nuo), listdptr(nuo)

      real(dp)
     .  Dnew(maxnd,nspin), e1, e2, ef, Enew(maxnd,nspin), Entropy,
     .  eo(maxo,nspin), H(maxnh,nspin), qo(maxo,nspin), 
     .  qtot, qs(nspin), S(maxnh), temp, efs(nspin)
     
      real(dp)
     .  Haux(nuotot,nuo), Saux(nuotot,nuo), psi(nuotot,maxuo,nspin), 
     .  aux(5*nuotot)

      logical
     .  getD, getPSI, fixspin

      external rdiag

C  Internal variables .............................................
      integer           ie, io, iio, ispin, ix, j, jo, BNode, iie, ind,
     .                  BTest, ierror, nd
      integer, allocatable, save :: HighestOccupiedBand(:)
      logical           lfound
      real*8            ee, pipj, qe, t, k(3)
C  ....................

C Solve eigenvalue problem .........................................
      do ispin = 1,nspin
   10   do io = 1,nuo
          do jo = 1,nuotot
            Saux(jo,io) = 0.0d0
            Haux(jo,io) = 0.0d0
          enddo
        enddo
        do io = 1,nuo
          do j = 1,numh(io)
            ind = listhptr(io) + j
            jo = listh(ind)
            Saux(jo,io) = Saux(jo,io) + S(ind)
            Haux(jo,io) = Haux(jo,io) + H(ind,ispin)
          enddo
        enddo
        call rdiag( Haux, Saux, nuotot, nuo, nuotot, maxuo,
     .              eo(1,ispin), psi(1,1,ispin), nuotot, ierror )
        if (ierror.gt.0) then
          call die('Terminating due to failed diagonalisation')
        elseif (ierror.eq.-1) then
          goto 10
        endif

        if (getPSI) then
          do ix = 1,3
            k(ix) = 0.0d0
          enddo
          call writew(nuotot,nuo,1,k,ispin,
     .                eo(1,ispin),psi(1,1,ispin),.true.)
        endif

      enddo

C Check if we are done ................................................
      if (.not.getD) return

C Find new Fermi energy and occupation weights ........................
      if (fixspin) then
        call fermispin( nspin, nspin, 1, (/1.d0/), maxo, no, eo,
     .               temp, qs, qo, efs, Entropy )
      else
        call fermid( nspin, nspin, 1, (/1.d0/), maxo, no, eo, 
     .             temp, qtot, qo, ef, Entropy )
      endif

C Find weights for local density of states ............................
      if (e1 .lt. e2) then
*       e1 = e1 - ef
*       e2 = e2 - ef
        t = max( temp, 1.d-6 )
        do ispin = 1,nspin
          do io = 1,nuotot
            qo(io,ispin) = ( stepf((eo(io,ispin)-e2)/t) -
     .                       stepf((eo(io,ispin)-e1)/t)) * 2.0d0/nspin
          enddo
        enddo
      endif
      
C Create pointer to highest occupied band
      allocate(HighestOccupiedBand(nspin))
      do ispin = 1,nspin
        io = nuotot + 1
        lfound = .false.
        do while (io.gt.1.and..not.lfound)
          io = io - 1
          lfound = (abs(qo(io,ispin)).lt.1.0d-12)
        enddo
        HighestOccupiedBand(ispin) = io
      enddo

C New density and energy-density matrices of unit-cell orbitals .......
      nd = listdptr(nuo) + numd(nuo)
      Dnew(1:nd,1:nspin) = 0.0d0
      Enew(1:nd,1:nspin) = 0.0d0

C Global operation to form new density matrix
      do ispin = 1,nspin
        BNode = 0
        iie = 0
        do ie = 1,HighestOccupiedBand(ispin)
          if (Node.eq.BNode) then
            iie = iie + 1
            do j = 1,nuotot
              aux(j) = psi(j,iie,ispin)
            enddo
          endif
#ifdef MPI
#ifdef NODAT
          call MPI_Bcast(aux,nuotot,MPI_double_precision,BNode,
     .      MPI_Comm_World,MPIerror)
#else
          call MPI_Bcast(aux,nuotot,DAT_double,BNode,
     .      MPI_Comm_World,MPIerror)
#endif
#endif
          qe = qo(ie,ispin)
          ee = qo(ie,ispin) * eo(ie,ispin)
          do io = 1,nuo
            call LocalToGlobalOrb(io,Node,Nodes,iio)
            do j = 1,numd(io)
              ind = listdptr(io) + j
              jo = listd(ind)
              pipj = aux(iio) * aux(jo)
              Dnew(ind,ispin) = Dnew(ind,ispin) + qe * pipj
              Enew(ind,ispin) = Enew(ind,ispin) + ee * pipj
            enddo
          enddo
          BTest = ie/BlockSize
          if (BTest*BlockSize.eq.ie) then
            BNode = BNode + 1
            if (BNode .gt. Nodes-1) BNode = 0
          endif
        enddo
      enddo

C Free local memory
      deallocate(HighestOccupiedBand)

      end
