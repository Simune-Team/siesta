      subroutine diagg( nspin, nuo, no, maxuo, maxno, maxnd, 
     .                  maxo, numh, listh, numd, listd, H, S,
     .                  getD, fixspin, qtot, qs, temp, e1, e2,
     .                  eo, qo, Dnew, Enew, ef, efs,
     .                  Haux, Saux, psi, aux, nuotot,
     .                  Node, Nodes )
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
C integer maxno               : Maximum number of orbitals interacting  
C                               with any orbital
C integer maxnd               : Maximum number of nonzero elements of 
C                               each row of density matrix
C integer maxo                : First dimension of eo and qo
C integer numh(no)            : Number of nonzero elements of each row 
C                               of hamiltonian matrix
C integer listh(maxno,no)     : Nonzero hamiltonian-matrix element  
C                               column indexes for each matrix row
C integer numd(no)            : Number of nonzero elements of each row 
C                               ofdensity matrix
C integer listd(maxnd,no)     : Nonzero density-matrix element column 
C                               indexes for each matrix row
C real*8  H(maxno,maxuo,nspin): Hamiltonian in sparse form
C real*8  S(maxno,maxuo)      : Overlap in sparse form
C logical getD                : Find occupations and density matrices?
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
C integer Node                : local processor
C integer Nodes               : total number of processors
C *************************** OUTPUT **********************************
C real*8 eo(maxo,nspn)           : Eigenvalues
C ******************** OUTPUT (only if getD=.true.) *******************
C real*8 qo(maxo,nspn)           : Occupations of eigenstates
C real*8 Dnew(maxno,maxuo,nspin) : Output Density Matrix
C real*8 Enew(maxno,maxuo,nspin) : Output Energy-Density Matrix
C real*8 ef                      : Fermi energy
C real*8 efs(nspin)              : Fermi energy for each spin
C                                  (for fixed spin calculations)
C *************************** AUXILIARY *******************************
C real*8 Haux(nuotot,no)      : Auxiliary space for the hamiltonian matrix
C real*8 Saux(nuotot,no)      : Auxiliary space for the overlap matrix
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
      use parallel
#ifdef MPI
      use mpi
#endif

      implicit none

#ifdef MPI
      integer 
     .  MPIerror
#endif

      integer
     .  maxnd, maxno, maxuo, maxo, nuo, no, nspin, nuotot, 
     .  Node, Nodes

      integer 
     .  listh(maxno,no), numh(no),
     .  listd(maxnd,no), numd(no)

      double precision
     .  Dnew(maxnd,maxuo,nspin),
     .  e1, e2, ef, Enew(maxnd,maxuo,nspin), eo(maxo,nspin),
     .  H(maxno,maxuo,nspin), qo(maxo,nspin), qtot, qs(nspin),
     .  S(maxno,maxuo), stepf, temp, efs(nspin)
     
      double precision
     .  Haux(nuotot,nuo), Saux(nuotot,nuo), psi(nuotot,maxuo,nspin), 
     .  aux(5*nuotot)

      logical
     .  getD, fixspin

      external
     .  fermid, fermispin, rdiag, stepf

C  Internal variables .............................................
      integer           ie, io, iio, ispin, j, jo, BNode, iie, 
     .                  BTest, ierror
      double precision  ee, pipj, qe, t
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
          call LocalToGlobalOrb(io,Node,Nodes,iio)
          do j = 1,numh(iio)
            jo = listh(j,iio)
            Saux(jo,io) = Saux(jo,io) + S(j,io)
            Haux(jo,io) = Haux(jo,io) + H(j,io,ispin)
          enddo
        enddo
        call rdiag( Haux, Saux, no, nuo, nuotot, maxuo,
     .              eo(1,ispin), psi(1,1,ispin), aux, ierror )
        if (ierror.gt.0) then
          call die('Terminating due to failed diagonalisation')
        elseif (ierror.lt.0) then
          goto 10
        endif
      enddo
C ....................

C Check if we are done ................................................
      if (.not.getD) return
C ....................

C Find new Fermi energy and occupation weights ........................
      if (fixspin) then
        call fermispin( nspin, nspin, 1, 1.d0, maxuo, no, eo,
     .               temp, qs, qo, efs )
      else
        call fermid( nspin, nspin, 1, 1.d0, maxo, no, eo, 
     .             temp, qtot, qo, ef )
      endif
C ....................

C Find weights for local density of states ............................
      if (e1 .lt. e2) then
*       e1 = e1 - ef
*       e2 = e2 - ef
        t = max( temp, 1.d-6 )
        do ispin = 1,nspin
          do io = 1,nuotot
            qo(io,ispin) = ( stepf( (eo(io,ispin)-e2)/t ) -
     .                       stepf( (eo(io,ispin)-e1)/t ) ) / nspin
          enddo
        enddo
      endif
C ....................
      
C New density and energy-density matrices of unit-cell orbitals .......
      do ispin = 1,nspin
        do io = 1,nuo
          call LocalToGlobalOrb(io,Node,Nodes,iio)
          do j = 1,numd(iio)
            Dnew(j,io,ispin) = 0.d0
            Enew(j,io,ispin) = 0.d0
          enddo
        enddo
      enddo

C Global operation to form new density matrix
      do ispin = 1,nspin
        BNode = 0
        iie = 0
        do ie = 1,no
          if (Node.eq.BNode) then
            iie = iie + 1
            do j = 1,no
              aux(j) = psi(j,iie,ispin)
            enddo
          endif
#ifdef MPI
          call MPI_Bcast(aux,no,DAT_double,BNode,
     .      MPI_Comm_World,MPIerror)
#endif
          qe = qo(ie,ispin)
          ee = qo(ie,ispin) * eo(ie,ispin)
          do io = 1,nuo
            call LocalToGlobalOrb(io,Node,Nodes,iio)
            do j = 1,numd(iio)
              jo = listd(j,iio)
              pipj = aux(iio) * aux(jo)
              Dnew(j,io,ispin) = Dnew(j,io,ispin) + qe * pipj
              Enew(j,io,ispin) = Enew(j,io,ispin) + ee * pipj
            enddo
          enddo
          BTest = ie/BlockSize
          if (BTest*BlockSize.eq.ie) then
            BNode = BNode + 1
            if (BNode .gt. Nodes-1) BNode = 0
          endif
        enddo
      enddo
C ....................

      end
