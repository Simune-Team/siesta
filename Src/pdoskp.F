      subroutine pdoskp(nspin, nuo, no, maxspn, maxuo, maxnh, 
     .                  maxo, numh, listhptr, listh, H, S,
     .                  E1, E2, nhist, sigma, 
     .                  xij, indxuo, nk, kpoint, wk, eo, 
     .                  Haux, Saux, psi, dtot, dpr, nuotot )

C **********************************************************************
C Find the density of states projected onto the atomic orbitals
C     D_mu(E) = Sum(n,k,nu) C(mu,n,k) C(nu,n,k) S(mu,nu,k) Delta(E-E(n,k))
C where n run over all the bands between two given energies
C Written by J. Junquera and E. Artacho. Nov' 99
C Modified version for parallel execution over K points by J.D. Gale
C March 2005
C ****  INPUT  *********************************************************
C integer nspin             : Number of spin components (1 or 2)
C integer nuo               : Number of atomic orbitals in the unit cell
C integer NO                : Number of atomic orbitals in the supercell
C integer maxspn            : Second dimension of eo and qo 
C                             (maximum number of differents spin polarizations)
C integer maxuo             : Maximum number of atomic orbitals in the unit cell
C integer maxnh             : Maximum number of orbitals interacting
C                             with any orbital
C integer maxo              : First dimension of eo
C integer numh(nuo)         : Number of nonzero elements of each row
C                             of hamiltonian matrix
C integer listhptr(nuo)     : Pointer to each row (-1) of the
C                             hamiltonian matrix
C integer listh(maxnh)      : Nonzero hamiltonian-matrix element
C                             column indexes for each matrix row
C real*8  H(maxnh,nspin)    : Hamiltonian in sparse format
C real*8  S(maxnh)          : Overlap in sparse format
C real*8  E1, E2            : Energy range for density-matrix states
C                             (to find local density of states)
C                             Not used if e1 > e2
C integer nhist             : Number of the subdivisions of the histogram
C real*8  sigma             : Width of the gaussian to expand the eigenvectors
C real*8  xij(3,maxnh)      : Vectors between orbital centers (sparse)
C                             (not used if only gamma point)
C integer indxuo(no)        : Index of equivalent orbital in unit cell
C integer nk                : Number of k points
C real*8  kpoint(3,nk)      : k point vectors
C real*8  wk(nk)            : Weights for k points
C real*8  eo(maxo,maxspn,nk): Eigenvalues
C integer nuotot            : Total number of orbitals per unit cell
C ****  AUXILIARY  *****************************************************
C real*8  Haux(2,nuo,nuo)   : Auxiliary space for the hamiltonian matrix
C real*8  Saux(2,nuo,nuo)   : Auxiliary space for the overlap matrix
C real*8  psi(2,nuo,nuo)    : Auxiliary space for the eigenvectors
C ****  OUTPUT  ********************************************************
C real*8  dtot(nhist,2)   : Total density of states
C real*8  dpr(nhist,nuo,2): Proyected density of states
C **********************************************************************

      use precision
      use parallel,     only : Node, Nodes
      use parallelsubs, only : GetNodeOrbs, LocalToGlobalOrb
      use parallelsubs, only : WhichNodeOrb, GlobalToLocalOrb
#ifdef MPI
      use mpi_siesta
#endif
      use sys,          only : die

      implicit none

      integer
     .  nspin, nuo, no, maxspn, maxuo, maxnh, NK, 
     .  maxo, nhist, nuotot

      integer
     .  numh(nuo), listhptr(nuo), listh(maxnh),
     .  indxuo(no)

      real(dp)
     .  H(maxnh,nspin), S(maxnh), E1, E2, sigma, 
     .  xij(3,maxnh), kpoint(3,nk), eo(maxo,maxspn,nk),
     .  Haux(2,nuotot,nuo), Saux(2,nuotot,nuo), psi(2,nuotot,nuo),
     .  dtot(nhist,2), dpr(nhist,nuotot,2), wk(nk)

C Internal variables ---------------------------------------------------
      integer
     .  ik, is, iio, io, iuo, juo, j, jo, ihist, iband, ind, ierror,
     .  maxnhg, nuog, BNode

      integer, dimension(:), allocatable, save ::
     .  numhg, listhptrg, listhg

      real(dp)
     .  kxij, Ckxij, Skxij, delta, ener, diff, pi, pipj1, pipj2, 
     .  pipjS1, pipjS2, gauss, norm, wksum

      real(dp), dimension(:), allocatable, save ::
     .  Snew, Dloc

      real(dp), dimension(:,:), allocatable, save ::
     .  Hnew, xijloc

#ifdef MPI
      integer ::
     .  MPIerror
      real(dp), dimension(:,:,:), allocatable, save ::
     .  Sloc
#endif

      external
     .  cdiag

C Initialize some variables
      delta = (E2 - E1)/nhist
      pi = 4.0d0 * atan(1.0D0)

C Globalise list arrays - assumes listh and listd are the same

C Allocate local memory for global list arrays
      allocate(numhg(nuotot))
      call memory('A','I',nuotot,'pdoskp')
      allocate(listhptrg(nuotot))
      call memory('A','I',nuotot,'pdoskp')

C Find maximum value in numh and create local storage
      allocate(Dloc(nuotot))
      call memory('A','D',nuotot,'pdoskp')

C Globalise numh
      do io = 1,nuotot
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          numhg(io) = numh(iio)
        endif
#ifdef MPI
        call MPI_Bcast(numhg(io),1,MPI_integer,BNode,
     .    MPI_Comm_World,MPIerror)
#endif
      enddo

C Build global listhptr
      listhptrg(1) = 0
      do io = 2,nuotot
        listhptrg(io) = listhptrg(io-1) + numhg(io-1)
      enddo

C Globalse listh
      maxnhg = listhptrg(nuotot) + numhg(nuotot)
      allocate(listhg(maxnhg))
      call memory('A','I',maxnhg,'pdoskp')
      do io = 1,nuotot
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          do jo = 1,numhg(io)
            listhg(listhptrg(io)+1:listhptrg(io)+numhg(io)) = 
     .        listh(listhptr(iio)+1:listhptr(iio)+numh(iio))
          enddo
        endif
#ifdef MPI
        call MPI_Bcast(listhg(listhptrg(io)+1),numhg(io),MPI_integer,
     .    BNode,MPI_Comm_World,MPIerror)
#endif
      enddo

C Create new distribution of H and S
      nuog = nuotot
      allocate(Snew(maxnhg))
      call memory('A','D',maxnhg,'diagkp')
      allocate(Hnew(maxnhg,nspin))
      call memory('A','D',maxnhg*nspin,'diagkp')
      allocate(xijloc(3,maxnhg))
      call memory('A','D',3*maxnhg,'diagkp')

      do io = 1,nuotot
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          do is = 1,nspin
            do jo = 1,numh(iio)
              Hnew(listhptrg(io)+jo,is) = H(listhptr(iio)+jo,is)
            enddo
          enddo
          do jo = 1,numh(iio)
            Snew(listhptrg(io)+jo) = S(listhptr(iio)+jo)
          enddo
          do jo = 1,numh(iio)
            xijloc(1:3,listhptrg(io)+jo) = xij(1:3,listhptr(iio)+jo)
          enddo
        endif
#ifdef MPI
        do is = 1,nspin
          call MPI_Bcast(Hnew(listhptrg(io)+1,is),numhg(io),
     .      MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
        enddo
        call MPI_Bcast(Snew(listhptrg(io)+1),numhg(io),
     .    MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
        call MPI_Bcast(xijloc(1,listhptrg(io)+1),3*numhg(io),
     .    MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
#endif
      enddo

C Solve eigenvalue problem for each k-point
      do is = 1,nspin

        do ik = 1+Node,nk,Nodes

C Initialize auxiliary variables 
          do iuo = 1,nuog
            do juo = 1,nuotot
              Saux(1,juo,iuo) = 0.0d0
              Saux(2,juo,iuo) = 0.0d0
              Haux(1,juo,iuo) = 0.0d0
              Haux(2,juo,iuo) = 0.0d0
            enddo
          enddo

          do io = 1,nuog
            do j = 1,numhg(io)
              ind = listhptrg(io) + j
              jo  = listhg(ind)
              iuo = indxuo(io)
              juo = indxuo(jo)
C Calculate the phases k*r_ij
              kxij = kpoint(1,ik) * xijloc(1,ind) +
     .               kpoint(2,ik) * xijloc(2,ind) +
     .               kpoint(3,ik) * xijloc(3,ind) 
              Ckxij = cos(kxij)
              Skxij = sin(kxij)
C Calculate the Hamiltonian and the overlap in k space
C H(k) = Sum(R) exp(i*k*R) * H(R)
              Saux(1,juo,iuo) = Saux(1,juo,iuo) + Snew(ind) * Ckxij
              Saux(2,juo,iuo) = Saux(2,juo,iuo) - Snew(ind) * Skxij
              Haux(1,juo,iuo) = Haux(1,juo,iuo) + Hnew(ind,is) * Ckxij
              Haux(2,juo,iuo) = Haux(2,juo,iuo) - Hnew(ind,is) * Skxij
            enddo
          enddo

C Diagonalize for each k point
          call cdiag( Haux, Saux, nuotot, nuog, nuotot,
     .                eo(1,is,ik), psi, nuotot, 1, ierror )

C Check error flag and take appropriate action
          if (ierror.gt.0) then
            call die('Terminating due to failed diagonalisation')
          elseif (ierror.lt.0) then
C Repeat diagonalisation with increased memory to handle clustering
            do iuo = 1,nuog
              do juo = 1,nuotot
                Saux(1,juo,iuo) = 0.0d0
                Saux(2,juo,iuo) = 0.0d0
                Haux(1,juo,iuo) = 0.0d0
                Haux(2,juo,iuo) = 0.0d0
              enddo
            enddo
            do io = 1, nuog
              do j = 1, numhg(io)
                ind = listhptrg(io) + j
                jo  = listhg(ind)
                iuo = indxuo(io)
                juo = indxuo(jo)
                kxij = kpoint(1,ik) * xijloc(1,ind) +
     .                 kpoint(2,ik) * xijloc(2,ind) +
     .                 kpoint(3,ik) * xijloc(3,ind) 
                Ckxij = cos(kxij)
                Skxij = sin(kxij)
                Saux(1,juo,iuo) = Saux(1,juo,iuo) + Snew(ind)*Ckxij
                Saux(2,juo,iuo) = Saux(2,juo,iuo) - Snew(ind)*Skxij
                Haux(1,juo,iuo) = Haux(1,juo,iuo) + Hnew(ind,is)*Ckxij
                Haux(2,juo,iuo) = Haux(2,juo,iuo) - Hnew(ind,is)*Skxij
              enddo
            enddo
            call cdiag( Haux, Saux, nuotot, nuog, nuotot,
     .                  eo(1,is,ik), psi, nuotot, 1, ierror )
          endif

C Recalculate again the overlap matrix in k-space
          do iuo = 1,nuog
            do juo = 1,nuotot
              Saux(1,juo,iuo) = 0.0d0
              Saux(2,juo,iuo) = 0.0d0
            enddo
          enddo

          do io = 1,nuog
            do j = 1,numhg(io)
              ind = listhptrg(io) + j
              jo  = listhg(ind)
              iuo = indxuo(io)
              juo = indxuo(jo)
C Calculates the phases k*r_ij
              kxij = kpoint(1,ik) * xijloc(1,ind) +
     .               kpoint(2,ik) * xijloc(2,ind) +
     .               kpoint(3,ik) * xijloc(3,ind) 
              ckxij = cos(kxij)
              skxij = sin(kxij)
C Calculates the hamiltonian and the overlap in k space
C H(k) = Sum(R) exp(i*k*R) * H(R)
              Saux(1,juo,iuo) = Saux(1,juo,iuo) + Snew(ind) * ckxij
              Saux(2,juo,iuo) = Saux(2,juo,iuo) - Snew(ind) * skxij
            enddo
          enddo

C Loop over all the energy range
          do ihist = 1,nhist
            ener = E1 + (ihist - 1) * delta
            do 170 iband = 1,nuog
              diff = (ener - eo(iband,is,ik))**2 / (sigma**2)
              if (diff .gt. 15.0d0) then
                goto 170
              else
                gauss = ( exp(-diff) )
                dtot(ihist,is) = dtot(ihist,is) + gauss*wk(ik)
                do iuo = 1, nuotot
C Solo para los Juo que satisfagan el criterio del record...
                  do juo = 1, nuotot
                    pipj1 = psi(1,iuo,iband) * psi(1,juo,iband) +
     .                      psi(2,iuo,iband) * psi(2,juo,iband)
                    pipj2 = psi(1,iuo,iband) * psi(2,juo,iband) -
     .                      psi(2,iuo,iband) * psi(1,juo,iband)
                    pipjS1= pipj1*Saux(1,iuo,juo)-pipj2*Saux(2,iuo,juo)
                    pipjS2= pipj1*Saux(2,iuo,juo)+pipj2*Saux(1,iuo,juo)
                    dpr(ihist,juo,is) = dpr(ihist,juo,is) + 
     .                                     pipjS1*gauss*wk(ik)
                  enddo
                enddo
              endif
 170        enddo

          enddo

        enddo

      enddo

C Free local memory from computation of dpr
      call memory('D','D',size(xijloc),'diagkp')
      deallocate(xijloc)
      call memory('D','D',size(Hnew),'diagkp')
      deallocate(Hnew)
      call memory('D','D',size(Snew),'diagkp')
      deallocate(Snew)
      call memory('D','D',size(Dloc),'diagkp')
      deallocate(Dloc)
      call memory('D','I',size(listhg),'diagkp')
      deallocate(listhg)
      call memory('D','I',size(listhptrg),'diagkp')
      deallocate(listhptrg)
      call memory('D','I',size(numhg),'diagkp')
      deallocate(numhg)

#ifdef MPI
C Allocate workspace array for global reduction
      allocate(Sloc(nhist,max(nuotot,nspin),nspin))
      call memory('A','D',nhist*nuotot*nspin,'pdoskp')

C Global reduction of dpr matrix
      Sloc(1:nhist,1:nuotot,1:nspin) = 0.0d0
      call MPI_AllReduce(dpr(1,1,1),Sloc(1,1,1),nhist*nuotot*nspin,
     .  MPI_double_precision,MPI_sum,MPI_Comm_World,MPIerror)
      dpr(1:nhist,1:nuotot,1:nspin) = Sloc(1:nhist,1:nuotot,1:nspin)

C Global reduction of dtot matrix
      Sloc(1:nhist,1:nspin,1) = 0.0d0
      call MPI_AllReduce(dtot(1,1),Sloc(1,1,1),nhist*nspin,
     .  MPI_double_precision,MPI_sum,MPI_Comm_World,MPIerror)
      dtot(1:nhist,1:nspin) = Sloc(1:nhist,1:nspin,1)

C Free workspace array for global reduction
      call memory('D','D',size(Sloc),'pdoskp')
      deallocate(Sloc)
#endif

      wksum = 0.0d0
      do ik = 1,nk
        wksum = wksum + wk(ik)
      enddo

      norm = sigma * sqrt(pi) * wksum

      do ihist = 1,nhist
        do is = 1,nspin
          dtot(ihist,is) = dtot(ihist,is) / norm
          do iuo = 1,nuotot
            dpr(ihist,iuo,is) = dpr(ihist,iuo,is) /norm
          enddo
        enddo
      enddo

      return
      end
