      subroutine eandg(iopt,eta,enum,lam,maxnh,numh,listhptr,listh,
     .                 ncmax,numc,listc,h,s,c,nbasis,nbands,e3,e,
     .                 grad,dm,edm,nbasisloc,maxnct)
C **************************************************************************
C Subroutine to link the CG algorithms to the calculation of
C the functional energy and gradients, and to set up control
C vectors for auxiliary sparse matrices.
C It also computes the density matrix.
C This routine works with the funcional of Kim et al. (PRB 52, 1640 (95)).
C Written by P.Ordejon. October'96
C ****************************** INPUT *************************************
C integer iopt             : Input option parameter
C                              iopt = 0  => Set up control vectors
C                              iopt = 1  => Call energy routine for line. min.
C                              iopt = 2  => Call gradient routine
C                              iopt = 3  => Call density matrix routine
C real*8 eta               : Fermi level parameter of Kim et al
C real*8 enum              : Total number of electrons
C real*8 lam               : Length of step for line minimization
C integer maxnh            : First dimension of listh and H, and maximum
C                            number of nonzero elements of each row of H
C integer numh(nbasisloc)  : Control vector of H matrix
C                            (number of <>0 element of each row)
C integer listhptr(nbasisloc) : Control vector of H matrix
C                            (pointer to start of each row (-1) of H matrix)
C integer listh(maxnh)     : Control vector of H matrix
C                            (list of <>0 element of each row)
C integer ncmax            : First dimension of listc and C, and maximum
C                            number of nonzero elements of each row of C
C integer numc(nbasis)     : Control vector of C matrix
C                            (number of <>0 element of each row)
C integer listc(ncmax,nbasis): Control vector of C matrix
C                            (list of <>0 element of each row)
C real*8 h(maxnh)          : Hamiltonian in sparse form
C real*8 s(maxnh)          : Overlap in sparse form
C real*8 c(ncmax,nbasis)   : Current point (wave func. coeffs. in sparse)
C integer nbasis           : Number of atomic orbitals
C integer nbands           : Number of Localized Wave Functions
C integer nbasisloc        : Number of local atomic orbitals
C integer maxnct           :) needed
C ******** INPUT OR OUTPUT (DEPENDING ON ARGUMENT IOPT) ********************
C real*8 grad(ncmax,nbasis)  : Gradient of the functional
C                            (input if iopt = 1)
C                            (output if iopt = 2)
C ***************************** OUTPUT *************************************
C real*8 e(3)              : Value of the energy in three points C+LAM_i*GRAD
C real*8 e                 : Value of the energy at point C
C real*8 dm(maxnh)         : Density matrix in sparse form
C real*8 edm(maxnh)        : Energy density matrix in sparse form
C **************************** BEHAVIOUR ***********************************
C The overlap matrix 'o' must be in the same sparse format as the 
C Hamiltonian matrix 'h', even if the overlap is more sparse than h
C (as due to the KB projectors, for instance). It will, in general,
C contain some zeros, therefore.
C **************************************************************************

C
C  Modules
C
      use precision
#ifdef MPI
      use mpi
#endif

      implicit none

      integer
     .  iopt,nbands,nbasis,ncmax,maxnh,nbasisloc,maxnct

      integer
     .  listc(ncmax,nbasis),listh(maxnh),listhptr(nbasisloc),
     .  numc(nbasis),numh(nbasisloc)

      double precision
     .  c(ncmax,nbasisloc),e,e3(3),dm(maxnh),edm(maxnh),enum,eta,
     .  h(maxnh),grad(ncmax,nbasisloc),lam,s(maxnh)

      external
     .  axb_build1,axb_build2,axb_build1s,axb_build2s,
     .  ctrans1,ctrans2,ctrans2s,ener3,gradient,memory

C Internal variables .....................................................
      integer 
     .  Node, Nodes, nbandsmin, nbandsmax, nbandspernode, nbandsloc,
     .  nremainder, maxnf, maxnft, maxnhij

#ifdef MPI
      integer MPIerror
#endif

      integer, dimension(:), allocatable, save ::
     .  numct, numf, numft, numhij

      integer, dimension(:,:), allocatable, save ::
     .  cttoc, fttof, indgf, listct, listf, listft, listhij

      integer, save ::
     .  maxo, maxoloc, maxlwf, maxnfl, maxnftl, maxnctl,
     .  maxncl, maxnhijl, maxlwfl

      data maxo     / 0 /
      data maxoloc  / 0 /
      data maxlwf  / 0 /
      data maxlwfl / 0 /
      data maxnfl  / 0 /
      data maxnftl / 0 /
      data maxncl  / 0 /
      data maxnctl / 0 /
      data maxnhijl / 0 /
C .....................

C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
#else
      Node = 0
      Nodes = 1
#endif

*     call timer('eandg',1)

C Work out range of bands for this Node
      nbandspernode = ((nbands-1)/Nodes)
      nremainder = nbands - nbandspernode*Nodes
      nbandsmin = Node*nbandspernode + min(nremainder,Node) + 1
      nbandsmax = nbandsmin + nbandspernode - 1
      if (Node.lt.nremainder) nbandsmax = nbandsmax + 1
      nbandsloc = nbandsmax - nbandsmin + 1

C Check array dimensions/handle dynamic memory
      if (nbasis .gt. maxo) then
        if (allocated(numft)) then
          call memory('D','I',size(numft),'eandg')
          deallocate(numft)
        endif
        allocate(numft(nbasis))
        call memory('A','I',nbasis,'eandg')
      endif
      if (nbasisloc .gt. maxoloc .or. ncmax .gt. maxncl) then
        if (allocated(indgf)) then
          call memory('D','I',size(indgf),'eandg')
          deallocate(indgf)
        endif
        allocate(indgf(ncmax,nbasisloc))
        call memory('A','I',ncmax*nbasisloc,'eandg')
      endif
      if (nbandsloc .gt. maxlwfl) then
        if (allocated(numct)) then
          call memory('D','I',size(numct),'eandg')
          deallocate(numct)
        endif
        if (allocated(numf)) then
          call memory('D','I',size(numf),'eandg')
          deallocate(numf)
        endif
        if (allocated(numhij)) then
          call memory('D','I',size(numhij),'eandg')
          deallocate(numhij)
        endif
        allocate(numct(nbandsloc))
        call memory('A','I',nbandsloc,'eandg')
        allocate(numf(nbandsloc))
        call memory('A','I',nbandsloc,'eandg')
        allocate(numhij(nbandsloc))
        call memory('A','I',nbandsloc,'eandg')
      endif
      if (nbandsloc .gt. maxlwfl .or. maxnct .gt. maxnctl) then
        if (allocated(cttoc)) then
          call memory('D','I',size(cttoc),'eandg')
          deallocate(cttoc)
        endif
        if (allocated(listct)) then
          call memory('D','I',size(listct),'eandg')
          deallocate(listct)
        endif
        allocate(cttoc(maxnct,nbandsloc))
        call memory('A','I',maxnct*nbandsloc,'eandg')
        allocate(listct(maxnct,nbandsloc))
        call memory('A','I',maxnct*nbandsloc,'eandg')
      endif
      maxo = max(maxo,nbasis)
      maxlwf = max(maxlwf,nbands)
      maxncl = max(maxncl,ncmax)
      maxnctl = max(maxnctl,maxnct)

C .....................
C Set up index lists for sparse matrices ..................................
      if (iopt .eq. 0) then
C GET Ct LISTS
        call ctrans1(nbasis,nbandsmin,nbandsmax,ncmax,maxnct,
     .               numc,listc,numct,listct,cttoc)
C GET F list sizes
        call axb_build1s(nbandsmin,nbandsmax,nbasis,maxnct,numct,
     .                   listct,nbasis,nbasis,maxnh,numh,listhptr,
     .                   listh,maxnf,numf,Node,Nodes)
      else
        maxnf = maxnfl
      endif

C Adjust sizes based on maxnf
      if (nbandsloc .gt. maxlwfl .or. maxnf .gt. maxnfl) then
        if (allocated(listf)) then
          call memory('D','I',size(listf),'eandg')
          deallocate(listf)
        endif
        allocate(listf(maxnf,nbandsloc))
        call memory('A','I',maxnf*nbandsloc,'eandg')
      endif
      maxnfl = max(maxnfl,maxnf)

      if (iopt .eq. 0) then
C GET F LISTS
        call axb_build1(nbandsmin,nbandsmax,nbasis,maxnct,numct,listct,
     .                 nbasis,nbasis,maxnh,numh,listhptr,listh,maxnf,
     .                 numf,listf,Node,Nodes)
C GET indgf map
        call ind_gf(nbasis,nbasisloc,nbandsmin,nbandsmax,nbands,
     .              ncmax,maxnf,numc,listc,numf,listf,indgf,Node,
     .              Nodes)
C GET Ft list sizes
        call ctrans2s(nbands,nbasis,maxnf,maxnft,numf,listf,
     .                numft,Node,Nodes)
      else
        maxnft = maxnftl
      endif

C Adjust sizes based on information from ctrans2s
      if (nbasisloc .gt. maxoloc .or. maxnft .gt. maxnftl) then
        if (allocated(fttof)) then
          call memory('D','I',size(fttof),'eandg')
          deallocate(fttof)
        endif
        if (allocated(listft)) then
          call memory('D','I',size(listft),'eandg')
          deallocate(listft)
        endif
        allocate(fttof(maxnft,nbasisloc))
        call memory('A','I',maxnft*nbasisloc,'eandg')
        allocate(listft(maxnft,nbasisloc))
        call memory('A','I',maxnft*nbasisloc,'eandg')
      endif
      maxnftl = max(maxnftl,maxnft)
      maxoloc = max(maxoloc,nbasisloc)

      if (iopt .eq. 0) then
C GET Ft LISTS
        call ctrans2(nbands,nbasis,maxnf,maxnft,numf,listf,
     .               numft,listft,fttof,Node,Nodes)
C GET Hij list sizes
        call axb_build2s(nbandsmin,nbandsmax,nbasis,maxnf,numf,
     .                   listf,nbasis,nbands,ncmax,numc,listc,
     .                   maxnhij,numhij)
      else
        maxnhij = maxnhijl
      endif

C Adjust sizes based on information from axb_build2s
      if (nbandsloc .gt. maxlwfl .or. maxnhij .gt. maxnhijl) then
        if (allocated(listhij)) then
          call memory('D','I',size(listhij),'eandg')
          deallocate(listhij)
        endif
        allocate(listhij(maxnhij,nbandsloc))
        call memory('A','I',maxnhij*nbandsloc,'eandg')
      endif
      maxnhijl = max(maxnhijl,maxnhij)
      maxlwfl = max(maxlwfl,nbandsloc)

      if (iopt .eq. 0) then
C GET Hij LISTS
        call axb_build2(nbandsmin,nbandsmax,nbasis,maxnf,numf,
     .                  listf,nbasis,nbands,ncmax,numc,listc,
     .                  maxnhij,numhij,listhij)
        goto 999
      endif
C.........................

C Calculate the energy at three points of the line, for the
C CG line minimization .....................................................
      if (iopt .eq. 1) then
        call ener3(c,grad,lam,eta,enum,h,s,nbasis,nbands,ncmax,maxnct,
     .             maxnf,maxnh,maxnhij,numc,listc,numct,listct,cttoc,
     .             numf,listf,numh,listhptr,listh,numhij,listhij,
     .             e3,nbasisloc)
        goto 999
      endif
C.........................

C Calculate the energy and Gradient at current point .......................
      if (iopt .eq. 2) then
        call gradient(c,eta,enum,h,s,nbasis,nbands,ncmax,maxnct,maxnf,
     .                maxnft,maxnh,maxnhij,numc,listc,numct,listct,
     .                cttoc,numf,listf,numft,listft,fttof,numh,listhptr,
     .                listh,numhij,listhij,indgf,grad,e,nbasisloc)
        goto 999
      endif
C.........................

C Calculate density matrix .................................................
C-JMS Modified denmat argument list
      if (iopt .eq. 3) then
        call denmat(c,eta,h,s,enum,nbasis,nbands,ncmax,maxnct,maxnf,
     .              maxnft,maxnh,numc,listc,numct,listct,cttoc,numf,
     .              listf,numft,listft,fttof,numh,listhptr,listh,dm,
     .              edm,nbasisloc)
        goto 999
      endif

C Trap iopt out of range
      if (Node.eq.0) then
        stop 'Error in eandg: incorrect iopt'
      else
        stop
      endif
      
  999 continue
*     call timer('eandg',2)
      return
      end
