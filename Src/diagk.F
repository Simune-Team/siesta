      subroutine diagk( nspin, nuo, no, maxspn, maxuo, maxnh, maxnd, 
     .                  maxo, numh, listhptr, listh, numd, listdptr, 
     .                  listd, H, S, getD, getPSI, fixspin, qtot, qs, 
     .                  temp, e1, e2, xij, indxuo, nk, kpoint, wk,
     .                  eo, qo, Dnew, Enew, ef, efs, Entropy,
     .                  Haux, Saux, psi, Dk, Ek, aux, nuotot)
C *********************************************************************
C Subroutine to calculate the eigenvalues and eigenvectors, density
C and energy-density matrices, and occupation weights of each 
C eigenvector, for given Hamiltonian and Overlap matrices (including
C spin polarization). K-sampling version.
C Writen by J.Soler, August 1998.
C **************************** INPUT **********************************
C integer nspin               : Number of spin components (1 or 2)
C integer nuo                 : Number of basis orbitals in unit cell
C                               local to this processor
C integer no                  : Number of basis orbitals in supercell
C integer maxspn              : Second dimension of eo and qo
C integer maxuo               : Last dimension of xij
C                               Must be at least max(indxuo)
C integer maxo                : First dimension of eo and qo
C integer maxnh               : Maximum number of orbitals interacting  
C integer maxnd               : First dimension of listd and DM
C integer numh(nuo)           : Number of nonzero elements of each row 
C                               of hamiltonian matrix
C integer listhptr(nuo)       : Pointer to each row (-1) of the
C                               hamiltonian matrix
C integer listh(maxnh)        : Nonzero hamiltonian-matrix element  
C                               column indexes for each matrix row
C integer numd(nuo)           : Number of nonzero elements of each row 
C                               of density matrix
C integer listdptr(nuo)       : Pointer to each row (-1) of the
C                               density matrix
C integer listd(maxnd)        : Nonzero density-matrix element column 
C                               indexes for each matrix row
C real*8  H(maxnh,nspin)      : Hamiltonian in sparse form
C real*8  S(maxnh)            : Overlap in sparse form
C logical getD                : Find occupations and density matrices?
C logical getPSI              : Find and print wave functions?
C real*8  qtot                : Number of electrons in unit cell
C real*8  temp                : Electronic temperature 
C real*8  e1, e2              : Energy range for density-matrix states
C                               (to find local density of states)
C                               Not used if e1 > e2
C real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
C                               (not used if only gamma point)
C integer indxuo(no)          : Index of equivalent orbital in unit cell
C                               Unit cell orbitals must be the first in
C                               orbital lists, i.e. indxuo.le.nuo, with
C                               nuo the number of orbitals in unit cell
C integer nk                  : Number of k points
C real*8  kpoint(3,nk)        : k point vectors
C real*8  wk(nk)              : k point weights (must sum one)
C integer nuotot              : total number of orbitals per unit cell
C                               over all processors
C *************************** OUTPUT **********************************
C real*8 eo(maxo,maxspn,nk)   : Eigenvalues
C ******************** OUTPUT (only if getD=.true.) *******************
C real*8 qo(maxo,maxspn,nk)   : Occupations of eigenstates
C real*8 Dnew(maxnd,nspin)    : Output Density Matrix
C real*8 Enew(maxnd,nspin)    : Output Energy-Density Matrix
C real*8 ef                   : Fermi energy
C real*8 Entropy              : Electronic entropy
C *************************** AUXILIARY *******************************
C real*8 Haux(2,nuotot,nuo) : Auxiliary space for the hamiltonian matrix
C real*8 Saux(2,nuotot,nuo) : Auxiliary space for the overlap matrix
C real*8 psi(2,nuotot,nuo)  : Auxiliary space for the eigenvectors
C real*8 aux(2,nuotot*5)    : Extra auxiliary space
C real*8 Dk(2,nuotot,nuo)   : Aux. space that may be the same as Haux
C real*8 Ek(2,nuotot,nuo)   : Aux. space that may be the same as Saux
C *************************** UNITS ***********************************
C xij and kpoint must be in reciprocal coordinates of each other.
C temp and H must be in the same energy units.
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

      implicit          none

#ifdef MPI
      integer 
     .  MPIerror
#endif

      integer           maxnd, maxnh, maxspn, maxuo, maxo, nk, no,
     .                  nspin, nuo, nuotot
      integer           indxuo(no), listh(maxnh), numh(nuo),
     .                  listd(maxnd), numd(nuo), listhptr(nuo),
     .                  listdptr(nuo)
      real(dp)          Dnew(maxnd,nspin),
     .                  e1, e2, ef, efs(nspin), Enew(maxnd,nspin),
     .                  Entropy, eo(maxo,maxspn,nk), H(maxnh,nspin),
     .                  kpoint(3,nk), qo(maxo,maxspn,nk), qtot, 
     .                  S(maxnh), temp, wk(nk),
     .                  xij(3,maxnh), qs(nspin)
      real(dp)          Dk(2,nuotot,nuo), Ek(2,nuotot,nuo),
     .                  Haux(2,nuotot,nuo), Saux(2,nuotot,nuo),
     .                  psi(2,nuotot,nuo), aux(2,nuotot*5)
      logical           getD, getPSI, fixspin
      external          cdiag

C  Internal variables .............................................
      integer
     .  BNode, BTest, ie, ierror, iie, ik, ind, io, iio,
     .  ispin, iuo, j, jo, juo, nd

      integer, allocatable, save :: HighestOccupiedBand(:,:)

      logical
     .  lfound

      real(dp)
     .  ckxij, ee, kxij, pipj1, pipj2, qe, skxij, t
C  ....................

C Skip eigenvalues if DM is not required and WF's are required ....
      if (getPSI .and. .not. getD) goto 10

C Find eigenvalues ................................................
      do ik = 1,nk
        do ispin = 1,nspin
          Saux = 0.0d0
          Haux = 0.0d0
          do iuo = 1,nuo
            do j = 1,numh(iuo)
              ind = listhptr(iuo) + j
              jo = listh(ind)
              juo = indxuo(jo)
              kxij = kpoint(1,ik) * xij(1,ind) +
     .               kpoint(2,ik) * xij(2,ind) +
     .               kpoint(3,ik) * xij(3,ind)
              ckxij = cos(kxij)
              skxij = sin(kxij)
C Note : sign of complex part changed to match change in order of iuo/juo
              Saux(1,juo,iuo)=Saux(1,juo,iuo)+S(ind)*ckxij
              Saux(2,juo,iuo)=Saux(2,juo,iuo)-S(ind)*skxij
              Haux(1,juo,iuo)=Haux(1,juo,iuo)+H(ind,ispin)*ckxij
              Haux(2,juo,iuo)=Haux(2,juo,iuo)-H(ind,ispin)*skxij
            enddo
          enddo
          call cdiag( Haux, nuotot, Saux, nuotot, nuo,
     .                eo(1,ispin,ik), psi, nuotot, 0, ierror )
          if (ierror.ne.0) then
            call die('Terminating due to failed diagonalisation')
          endif
        enddo
      enddo

C Check if we are done ................................................
      if (.not.getD .and. .not. getPSI) return

C Find new Fermi energy and occupation weights ........................
      if (fixspin) then
        call fermispin( nspin, nspin, nk, wk, maxo, nuotot, eo,
     .               temp, qs, qo, efs, Entropy )
      else
        call fermid( nspin, maxspn, nk, wk, maxo, nuotot, eo, 
     .             temp, qtot, qo, ef, Entropy )
      endif

C Find weights for local density of states ............................
      if (e1 .lt. e2) then
*       e1 = e1 - ef
*       e2 = e2 - ef
        t = max( temp, 1.d-6 )
        do ik = 1,nk
          do ispin = 1,nspin
            do io = 1,nuotot
              qo(io,ispin,ik) = wk(ik) * 
     .             ( stepf((eo(io,ispin,ik)-e2)/t) -
     .               stepf((eo(io,ispin,ik)-e1)/t)) * 2.0d0/nspin
            enddo
          enddo
        enddo
      endif

C Create pointer to highest occupied band
      if (getD) then
        allocate(HighestOccupiedBand(nspin,nk))
        do ik = 1,nk
          do ispin = 1,nspin
            io = nuotot + 1
            lfound = .false.
            do while (io.gt.1.and..not.lfound)
              io = io - 1
              lfound = (abs(qo(io,ispin,ik)).lt.1.0d-12)
            enddo
            HighestOccupiedBand(ispin,ik) = io
          enddo
        enddo
      endif
      
C New density and energy-density matrices of unit-cell orbitals .......
      nd = listdptr(nuo) + numd(nuo)
      Dnew(1:nd,1:nspin) = 0.0d0
      Enew(1:nd,1:nspin) = 0.0d0

  10  continue

      do ik = 1,nk
        do ispin = 1,nspin

C Find eigenvectors 
          Saux = 0.0d0
          Haux = 0.0d0
          do iuo = 1,nuo
            do j = 1,numh(iuo)
              ind = listhptr(iuo) + j
              jo = listh(ind)
              juo = indxuo(jo)
              kxij = kpoint(1,ik) * xij(1,ind) +
     .               kpoint(2,ik) * xij(2,ind) +
     .               kpoint(3,ik) * xij(3,ind)
              ckxij = cos(kxij)
              skxij = sin(kxij)
              Saux(1,juo,iuo)=Saux(1,juo,iuo)+S(ind)*ckxij
              Saux(2,juo,iuo)=Saux(2,juo,iuo)-S(ind)*skxij
              Haux(1,juo,iuo)=Haux(1,juo,iuo)+H(ind,ispin)*ckxij
              Haux(2,juo,iuo)=Haux(2,juo,iuo)-H(ind,ispin)*skxij
            enddo
          enddo
          call cdiag( Haux, nuotot, Saux, nuotot, nuo,
     .                eo(1,ispin,ik), psi, nuotot, nuo, ierror )
          if (ierror.ne.0) then
            call die('Terminating due to failed diagonalisation')
          endif

          if (getPSI) 
     .     call writew(nuotot,nuo,ik,kpoint(1,ik),ispin,
     .                 eo(1,ispin,ik),psi,.false.)
 
          if (getD) then

C Add contribution to density matrices of unit-cell orbitals
C WARNING: Dk and Ek may be EQUIVALENCE'd to Haux and Saux
            do iuo = 1,nuo
              do juo = 1,nuotot
                Dk(1,juo,iuo) = 0.0d0
                DK(2,juo,iuo) = 0.0d0
                Ek(1,juo,iuo) = 0.0d0
                Ek(2,juo,iuo) = 0.0d0
              enddo
            enddo

C Global operation to form new density matrix
            BNode = 0
            iie = 0
            do ie = 1,HighestOccupiedBand(ispin,ik)
              if (Node.eq.BNode) then
                iie = iie + 1
                do j = 1,nuotot
                  aux(1,j) = psi(1,j,iie)
                  aux(2,j) = psi(2,j,iie)
                enddo
              endif
#ifdef MPI
#ifdef NODAT
              call MPI_Bcast(aux(1,1),2*nuotot,MPI_double_precision,
     .          BNode,MPI_Comm_World,MPIerror)
#else
              call MPI_Bcast(aux(1,1),2*nuotot,DAT_double,
     .          BNode,MPI_Comm_World,MPIerror)
#endif
#endif
              qe = qo(ie,ispin,ik)
              ee = qo(ie,ispin,ik) * eo(ie,ispin,ik)
              do iuo = 1,nuo
                call LocalToGlobalOrb(iuo,Node,Nodes,iio)
                do juo = 1,nuotot
                  pipj1 = aux(1,iio) * aux(1,juo) +
     .                    aux(2,iio) * aux(2,juo)
                  pipj2 = aux(1,iio) * aux(2,juo) -
     .                    aux(2,iio) * aux(1,juo)
                  Dk(1,juo,iuo) = Dk(1,juo,iuo) + qe * pipj1
                  Dk(2,juo,iuo) = Dk(2,juo,iuo) + qe * pipj2
                  Ek(1,juo,iuo) = Ek(1,juo,iuo) + ee * pipj1
                  Ek(2,juo,iuo) = Ek(2,juo,iuo) + ee * pipj2
                enddo
              enddo
              BTest = ie/BlockSize
              if (BTest*BlockSize.eq.ie) then
                BNode = BNode + 1
                if (BNode .gt. Nodes-1) BNode = 0
              endif
            enddo
            do iuo = 1,nuo
              do j = 1,numd(iuo)
                ind = listdptr(iuo) + j
                jo = listd(ind)
                juo = indxuo(jo)
                kxij = kpoint(1,ik) * xij(1,ind) +
     .                 kpoint(2,ik) * xij(2,ind) +
     .                 kpoint(3,ik) * xij(3,ind)
                ckxij = cos(kxij)
                skxij = sin(kxij)
                Dnew(ind,ispin)=Dnew(ind,ispin)+ Dk(1,juo,iuo)*ckxij -
     .                                           Dk(2,juo,iuo)*skxij
                Enew(ind,ispin)=Enew(ind,ispin)+ Ek(1,juo,iuo)*ckxij -
     .                                           Ek(2,juo,iuo)*skxij
              enddo
            enddo

          endif

        enddo
      enddo

C Free local memory
      if (getD) deallocate(HighestOccupiedBand)

      end
