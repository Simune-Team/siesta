      subroutine ener3(c,grad,lam,eta,enum,h,s,
     .                 nbasis,nbands,ncmax,nctmax,
     .                 nfmax,nhmax,nhijmax,
     .                 numc,listc,numct,listct,cttoc,numf,listf,
     .                 numh,listh,numhij,listhij,
     .                 ener,nbasisloc)

C ************************************************************************
C Finds the energy at three points of the line passing thru C in the
C direction of GRAD. LAM is the distance (in units of GRAD) between 
C points.
C Uses the functional of Kim et al (PRB 52, 1640 (95))
C Works only with spin-unpolarized systems
C Written by P.Ordejon. October'96
C ****************************** INPUT ***********************************
C real*8 c(ncmax,nbasisloc)    : Current point (wave function coeff.
C                                  in sparse form)
C real*8 grad(ncmax,nbasisloc) : Direction of search (sparse)
C real*8 lam                   : Length of step
C real*8 eta                   : Fermi level parameter of Kim et al.
C real*8 enum                  : Total number of electrons
C real*8 h(nhmax,nbasisloc)    : Hamiltonian matrix (sparse)
C real*8 s(nhmax,nbasisloc)    : Overlap matrix (sparse)
C integer nbasis               : Number of basis orbitals
C integer nbands               : Number of LWF's
C integer ncmax                : Max num of <>0 elements of each row of C
C integer nctmax               : Max num of <>0 elements of each col of C
C integer nfmax                : Max num of <>0 elements of each row of 
C                                   F = Ct x H
C integer nhmax                : Max num of <>0 elements of each row of H
C integer nhijmax              : Max num of <>0 elements of each row of 
C                                   Hij=Ct x H x C
C integer numc(nbasis)         : Control vector of C matrix
C                                (number of <>0  elements of each row of C)
C integer listc(ncmax,nbasis)  : Control vector of C matrix 
C                               (list of <>0  elements of each row of C)
C integer numct(nbands)        : Control vector of C transpose matrix
C                               (number of <>0  elements of each col of C)
C integer listct(ncmax,nbands) : Control vector of C transpose matrix
C                               (list of <>0  elements of each col of C)
C integer cttoc(ncmax,nbands)  : Map from Ct to C indexing
C integer numf(nbands)         : Control vector of F matrix
C                                (number of <>0  elements of each row of F)
C integer listf(nfmax,nbands)  : Control vector of F matrix 
C                                (list of <>0  elements of each row of F)
C integer numh(nbasis)         : Control vector of H matrix
C                                (number of <>0  elements of each row of H)
C integer listh(nhmax,nbasis)  : Control vector of H matrix 
C                               (list of <>0  elements of each row of H)
C integer numhij(nbasis)       : Control vector of Hij matrix
C                                (number of <>0  elements of each row of Hij)
C integer listhij(nhijmax,nbands): Control vector of Hij matrix 
C                                (list of <>0  elements of each row of Hij)
C integer nbasisloc            : Number of local atomic orbitals
C ***************************** OUTPUT ***********************************
C real*8 ener(3)               : Energy at the three points:
C                                     C +     lam * GRAD
C                                     C + 2 * lam * GRAD
C                                     C + 3 * lam * GRAD
C ************************************************************************

C
C  Modules
C
      use precision
      use parallel
#ifdef MPI
      use mpi
#endif

      implicit none

      integer
     .  nbasis,nbands,ncmax,nctmax,nfmax,nhmax,nhijmax,nbasisloc

      integer
     .  cttoc(nctmax,*),listc(ncmax,nbasis),
     .  listct(nctmax,*),listf(nfmax,*),
     .  listh(nhmax,nbasis),listhij(nhijmax,*),
     .  numc(nbasis),numct(*),numf(*),
     .  numh(nbasis),numhij(*)

      double precision
     .  c(ncmax,nbasisloc),ener(3),eta,enum,grad(ncmax,nbasisloc),
     .  h(nhmax,nbasisloc),lam,s(nhmax,nbasisloc)

C Internal variables ......................................................

      integer
     .  i,in,j,jn,k,kl,kn,Nodes,nn,nbandsmin,nbandsmax,
     .  nbandspernode, iloc, nbandsloc
      integer, dimension(:), allocatable, save :: iptr
#ifdef MPI
      integer
     .  MPIerror, nbl, nremainder, ii, Node
      double precision 
     .  etmp1(6), etmp2(6)
      double precision, dimension(:,:), allocatable, save :: 
     .  cl,gl,hl,sl
#endif

      double precision
     .  a1,a2,a3,b1,b2,b3,c1,c2,c3,func1(3),func2(3),
     .  lam1,lam2,lam3,pp1,pp2,pp3

      double precision, dimension(:,:,:), allocatable, save :: 
     .  aux, bux1, bux2

C..................

C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
      call GetNodeOrbs(nbasis,0,Nodes,nbl)
      allocate(cl(ncmax,nbl))
      call memory('A','D',ncmax*nbl,'ener3')
      allocate(gl(ncmax,nbl))
      call memory('A','D',ncmax*nbl,'ener3')
      allocate(hl(nhmax,nbl))
      call memory('A','D',nhmax*nbl,'ener3')
      allocate(sl(nhmax,nbl))
      call memory('A','D',nhmax*nbl,'ener3')

C Work out range of bands for this Node
      nbandspernode = ((nbands-1)/Nodes) 
      nremainder = nbands - nbandspernode*Nodes
      nbandsmin = Node*nbandspernode + min(nremainder,Node) + 1
      nbandsmax = nbandsmin + nbandspernode - 1
      if (Node.lt.nremainder) nbandsmax = nbandsmax + 1
      nbandsloc = nbandsmax - nbandsmin + 1
#else
      Nodes = 1
      nbandsmin = 1
      nbandsmax = nbands
      nbandspernode = nbands
      nbandsloc = nbands
#endif

      call timer('ener3',1)

C Allocate dynamic memory
      allocate(aux(6,nbasis,nbandspernode+1))
      call memory('A','D',6*nbasis*(nbandspernode+1),'ener3')
      allocate(bux1(3,nbands,nbandspernode+1))
      call memory('A','D',3*nbands*(nbandspernode+1),'ener3')
      allocate(bux2(3,nbands,nbandspernode+1))
      call memory('A','D',3*nbands*(nbandspernode+1),'ener3')
      allocate(iptr(nbasis))
      call memory('A','I',nbasis,'ener3')

C Define points to compute energy ..........................................
      lam1 = lam
      lam2 = lam*2.0d0
      lam3 = lam*3.0d0
C..................

C Initialize output and auxiliary varialbles ...............................
      do i = 1,3
        ener(i) = 0.0d0
        func1(i) = 0.0d0
        func2(i) = 0.0d0
      enddo

      do nn = 1,nbandspernode+1
        do i = 1,nbasis
          do j = 1,6
            aux(j,i,nn) = 0.0d0
          enddo
        enddo
      enddo

      do nn = 1,nbandspernode+1
        do i = 1,nbands
          do j = 1,3
            bux1(j,i,nn) = 0.0d0
            bux2(j,i,nn) = 0.0d0
          enddo
        enddo
      enddo

C..................

C Loop over Nodes for broadcasting of c,g,h and s
      do nn = 1,Nodes
        
#ifdef MPI
C Get local size on broadcast node
        call GetNodeOrbs(nbasis,nn-1,Nodes,nbl)

C Place arrays in local storage on broadcast node
        if (Node.eq.nn-1) then
          do i = 1,nbasisloc
            do ii = 1,ncmax
              cl(ii,i) = c(ii,i)
              gl(ii,i) = grad(ii,i)
            enddo
            do ii = 1,nhmax
              hl(ii,i) = h(ii,i)
              sl(ii,i) = s(ii,i)
            enddo
          enddo
        endif

C Broadcast arrays
        call MPI_Bcast(cl(1,1),ncmax*nbl,DAT_double,nn-1,
     .        MPI_Comm_World,MPIerror)
        call MPI_Bcast(gl(1,1),ncmax*nbl,DAT_double,nn-1,
     .        MPI_Comm_World,MPIerror)
        call MPI_Bcast(hl(1,1),nhmax*nbl,DAT_double,nn-1,
     .        MPI_Comm_World,MPIerror)
        call MPI_Bcast(sl(1,1),nhmax*nbl,DAT_double,nn-1,
     .        MPI_Comm_World,MPIerror)
#endif

C Build pointer for condensing f and fs on local node
        do k = 1,nbasis
          call GlobalToLocalOrb(k,nn-1,Nodes,kl)
          iptr(k) = kl
        enddo

C Calculate Functional .....................................................
C F=CtH
C Fs=CtS

      do i = 1,nbandsloc
        do in = 1,numct(i)
          k = listct(in,i)
          kl = iptr(k)
          if (kl.gt.0) then
#ifdef MPI
            pp1 = cl(cttoc(in,i),kl) + lam1 * gl(cttoc(in,i),kl)
            pp2 = cl(cttoc(in,i),kl) + lam2 * gl(cttoc(in,i),kl)
            pp3 = cl(cttoc(in,i),kl) + lam3 * gl(cttoc(in,i),kl)

            do kn = 1,numh(k)
              aux(1,listh(kn,k),i) = aux(1,listh(kn,k),i) + 
     .          pp1 * hl(kn,kl)
              aux(2,listh(kn,k),i) = aux(2,listh(kn,k),i) + 
     .          pp2 * hl(kn,kl)
              aux(3,listh(kn,k),i) = aux(3,listh(kn,k),i) + 
     .          pp3 * hl(kn,kl)
              aux(4,listh(kn,k),i) = aux(4,listh(kn,k),i) + 
     .          pp1 * sl(kn,kl)
              aux(5,listh(kn,k),i) = aux(5,listh(kn,k),i) + 
     .          pp2 * sl(kn,kl)
              aux(6,listh(kn,k),i) = aux(6,listh(kn,k),i) + 
     .          pp3 * sl(kn,kl)
            enddo
#else
            pp1 = c(cttoc(in,i),kl) + lam1 * grad(cttoc(in,i),kl)
            pp2 = c(cttoc(in,i),kl) + lam2 * grad(cttoc(in,i),kl)
            pp3 = c(cttoc(in,i),kl) + lam3 * grad(cttoc(in,i),kl)

            do kn = 1,numh(k)
              aux(1,listh(kn,k),i) = aux(1,listh(kn,k),i) + 
     .          pp1 * h(kn,kl)
              aux(2,listh(kn,k),i) = aux(2,listh(kn,k),i) + 
     .          pp2 * h(kn,kl)
              aux(3,listh(kn,k),i) = aux(3,listh(kn,k),i) + 
     .          pp3 * h(kn,kl)
              aux(4,listh(kn,k),i) = aux(4,listh(kn,k),i) + 
     .          pp1 * s(kn,kl)
              aux(5,listh(kn,k),i) = aux(5,listh(kn,k),i) + 
     .          pp2 * s(kn,kl)
              aux(6,listh(kn,k),i) = aux(6,listh(kn,k),i) + 
     .          pp3 * s(kn,kl)
            enddo
#endif
          endif
        enddo

      enddo

C End of loop over nodes
      enddo

C Loop over Nodes for broadcasting of c,g,h and s
      do nn = 1,Nodes

#ifdef MPI
C Get local size on broadcast node
        call GetNodeOrbs(nbasis,nn-1,Nodes,nbl)

C Place arrays in local storage on broadcast node
        if (Node.eq.nn-1) then
          do i = 1,nbasisloc
            do ii = 1,ncmax
              cl(ii,i) = c(ii,i)
              gl(ii,i) = grad(ii,i)
            enddo
          enddo
        endif

C Broadcast arrays
        call MPI_Bcast(cl(1,1),ncmax*nbl,DAT_double,nn-1,
     .        MPI_Comm_World,MPIerror)
        call MPI_Bcast(gl(1,1),ncmax*nbl,DAT_double,nn-1,
     .        MPI_Comm_World,MPIerror)
#endif

C Build pointer for condensing f and fs on local node
        do k = 1,nbasis
          call GlobalToLocalOrb(k,nn-1,Nodes,kl)
          iptr(k) = kl
        enddo

      do i = 1,nbandsloc

        do in = 1,numf(i)
          k = listf(in,i)
          kl = iptr(k)
          if (kl.gt.0) then
            a1 = aux(1,k,i)
            a2 = aux(2,k,i)
            a3 = aux(3,k,i)
            b1 = aux(4,k,i)
            b2 = aux(5,k,i)
            b3 = aux(6,k,i)

C Hij=CtHC
C Sij=CtSC
C multiply FxC and FsxC
            do kn = 1,numc(k)
#ifdef MPI
              c1 = cl(kn,kl) + lam1 * gl(kn,kl)
              c2 = cl(kn,kl) + lam2 * gl(kn,kl)
              c3 = cl(kn,kl) + lam3 * gl(kn,kl)
#else
              c1 = c(kn,kl) + lam1 * grad(kn,kl)
              c2 = c(kn,kl) + lam2 * grad(kn,kl)
              c3 = c(kn,kl) + lam3 * grad(kn,kl)
#endif
              bux1(1,listc(kn,k),i) = bux1(1,listc(kn,k),i) + 
     .          a1 * c1
              bux1(2,listc(kn,k),i) = bux1(2,listc(kn,k),i) + 
     .          a2 * c2
              bux1(3,listc(kn,k),i) = bux1(3,listc(kn,k),i) + 
     .          a3 * c3
              bux2(1,listc(kn,k),i) = bux2(1,listc(kn,k),i) + 
     .          b1 * c1
              bux2(2,listc(kn,k),i) = bux2(2,listc(kn,k),i) + 
     .          b2 * c2
              bux2(3,listc(kn,k),i) = bux2(3,listc(kn,k),i) + 
     .          b3 * c3
            enddo
          endif
        enddo
      enddo

C End of loop over nodes
      enddo

      do i = nbandsmin,nbandsmax
        iloc = i - nbandsmin + 1
        func1(1) = func1(1) + bux1(1,i,iloc) - eta * bux2(1,i,iloc)
        func1(2) = func1(2) + bux1(2,i,iloc) - eta * bux2(2,i,iloc)
        func1(3) = func1(3) + bux1(3,i,iloc) - eta * bux2(3,i,iloc)
        do jn = 1,numhij(iloc)
          func2(1) = func2(1) + (bux1(1,listhij(jn,iloc),iloc) -
     .             eta * bux2(1,listhij(jn,iloc),iloc)) * 
     .             bux2(1,listhij(jn,iloc),iloc)
          func2(2) = func2(2) + (bux1(2,listhij(jn,iloc),iloc) -
     .             eta * bux2(2,listhij(jn,iloc),iloc)) * 
     .             bux2(2,listhij(jn,iloc),iloc)
          func2(3) = func2(3) + (bux1(3,listhij(jn,iloc),iloc) -
     .             eta * bux2(3,listhij(jn,iloc),iloc)) * 
     .             bux2(3,listhij(jn,iloc),iloc)
        enddo
      enddo


#ifdef MPI
C Globalisation of func1 and func2
      do i=1,3
        etmp1(i) = func1(i)
        etmp1(3+i) = func2(i)
      enddo
      call MPI_AllReduce(etmp1,etmp2,6,DAT_double,MPI_sum,
     .  MPI_Comm_World,MPIerror)
      do i=1,3
        func1(i) = etmp2(i)
        func2(i) = etmp2(3+i)
      enddo
#endif

C This is valid for an spin-unpolarized sytem
      do i=1,3
        ener(i) = 2.0d0 * func1(i) - func2(i) + eta * enum / 2.0d0
      enddo

C ...................

C Free dynamic memory
#ifdef MPI
      call memory('D','D',size(cl),'ener3')
      deallocate(cl)
      call memory('D','D',size(gl),'ener3')
      deallocate(gl)
      call memory('D','D',size(hl),'ener3')
      deallocate(hl)
      call memory('D','D',size(sl),'ener3')
      deallocate(sl)
#endif
      call memory('D','D',size(aux),'ener3')
      deallocate(aux)
      call memory('D','D',size(bux1),'ener3')
      deallocate(bux1)
      call memory('D','D',size(bux2),'ener3')
      deallocate(bux2)
      call memory('D','I',size(iptr),'ener3')
      deallocate(iptr)

      call timer('ener3',2)
      return
      end
