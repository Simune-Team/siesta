      subroutine eandg(iopt,eta,enum,lam,
     .                 nhmax,numh,listh,ncmax,numc,listc,h,s,c,
     .                 nbasis,nbands,
     .                 e3,e,grad,dm,edm,nbasisloc,
     .                 maxnf,maxnft,maxnhij,maxnct)
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
C integer nhmax            : First dimension of listh and H, and maximum
C                            number of nonzero elements of each row of H
C integer numh(nbasis)       : Control vector of H matrix
C                            (number of <>0 element of each row)
C integer listh(nhmax,nbasis): Control vector of H matrix
C                            (list of <>0 element of each row)
C integer ncmax            : First dimension of listc and C, and maximum
C                            number of nonzero elements of each row of C
C integer numc(nbasis)       : Control vector of C matrix
C                            (number of <>0 element of each row)
C integer listc(ncmax,nbasis): Control vector of C matrix
C                            (list of <>0 element of each row)
C real*8 h(nhmax,nbasis)     : Hamiltonian in sparse form
C real*8 s(nhmax,nbasis)     : Overlap in sparse form
C real*8 c(ncmax,nbasis)     : Current point (wave func. coeffs. in sparse)
C integer nbasis           : Number of atomic orbitals
C integer nbands           : Number of Localized Wave Functions
C integer nbasisloc        : Number of local atomic orbitals
C integer maxnf            :)
C integer maxnft           :) Array dimensions to check against local values
C integer maxnhij          :) for saved arrays to see if re-dimensioning is
C integer maxnct           :) needed
C ******** INPUT OR OUTPUT (DEPENDING ON ARGUMENT IOPT) ********************
C real*8 grad(ncmax,nbasis)  : Gradient of the functional
C                            (input if iopt = 1)
C                            (output if iopt = 2)
C ***************************** OUTPUT *************************************
C real*8 e(3)              : Value of the energy in three points C+LAM_i*GRAD
C real*8 e                 : Value of the energy at point C
C real*8 dm(nhmax,nbasis)  : Density matrix in sparse form
C real*8 edm(nhmax,nbasis) : Energy density matrix in sparse form
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
     .  iopt,nbands,nbasis,ncmax,nhmax,nbasisloc,maxnf,maxnft,maxnhij,
     .  maxnct

      integer
     .  listc(ncmax,nbasis),listh(nhmax,nbasis),
     .  numc(nbasis),numh(nbasis)

      double precision
     .  c(ncmax,nbasisloc),e,e3(3),
     .  dm(nhmax,nbasisloc),edm(nhmax,nbasisloc),enum,eta,
     .  h(nhmax,nbasisloc),grad(ncmax,nbasisloc),lam,
     .  s(nhmax,nbasisloc)

      external
     .  axb_build,ctrans1,ctrans2,ener3,gradient,memory

C Internal variables .....................................................
      integer 
     .  Node, Nodes, nbandsmin, nbandsmax, nbandspernode,
     .  nremainder

#ifdef MPI
      integer MPIerror
#endif

      integer, dimension(:), allocatable ::
     .  nindv, ind

      integer, dimension(:), allocatable, save ::
     .  numct, numf, numft, numhij

      integer, dimension(:,:), allocatable, save ::
     .  cttoc, fttof, listct, listf, listft, listhij

      logical, save ::
     .  frstme

      integer, save ::
     .  maxo, maxoloc, maxlwf, maxnfl, maxnftl, maxnctl,
     .  maxnhijl

      data frstme   /.true./
      data maxo     / 0 /
      data maxoloc  / 0 /
      data maxlwf  / 0 /
      data maxnfl  / 0 /
      data maxnftl / 0 /
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

C Check array dimensions/handle dynamic memory
      if (nbasis .gt. maxo) then
        if (allocated(numft)) then
          call memory('D','I',size(numft),'eandg')
          deallocate(numft)
        endif
        allocate(numft(nbasis))
        call memory('A','I',nbasis,'eandg')
      endif
      if (nbands .gt. maxlwf) then
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
        allocate(numct(nbands))
        call memory('A','I',nbands,'eandg')
        allocate(numf(nbands))
        call memory('A','I',nbands,'eandg')
        allocate(numhij(nbands))
        call memory('A','I',nbands,'eandg')
      endif
      if (nbasisloc .gt. maxoloc .or. maxnft .gt. maxnftl) then
        if (allocated(fttof)) then
          call memory('D','I',size(fttof),'eandg')
          deallocate(fttof)
        endif
        if (allocated(listft)) then
          call memory('D','I',size(listft),'eandg')
          deallocate(listft)
        endif
        allocate(fttof(maxnf,nbasisloc))
        call memory('A','I',maxnf*nbasisloc,'eandg')
        allocate(listft(maxnf,nbasisloc))
        call memory('A','I',maxnf*nbasisloc,'eandg')
      endif
      if (nbands .gt. maxlwf .or. maxnct .gt. maxnctl) then
        if (allocated(cttoc)) then
          call memory('D','I',size(cttoc),'eandg')
          deallocate(cttoc)
        endif
        if (allocated(listct)) then
          call memory('D','I',size(listct),'eandg')
          deallocate(listct)
        endif
        allocate(cttoc(maxnct,nbands))
        call memory('A','I',maxnct*nbands,'eandg')
        allocate(listct(maxnct,nbands))
        call memory('A','I',maxnct*nbands,'eandg')
      endif
      if (nbands .gt. maxlwf .or. maxnf .gt. maxnfl) then
        if (allocated(listf)) then
          call memory('D','I',size(listf),'eandg')
          deallocate(listf)
        endif
        allocate(listf(maxnf,nbands))
        call memory('A','I',maxnf*nbands,'eandg')
      endif
      if (nbands .gt. maxlwf .or. maxnhij .gt. maxnhijl) then
        if (allocated(listhij)) then
          call memory('D','I',size(listhij),'eandg')
          deallocate(listhij)
        endif
        allocate(listhij(maxnhij,nbands))
        call memory('A','I',maxnhij*nbands,'eandg')
      endif
      maxo = max(maxo,nbasis)
      maxoloc = max(maxoloc,nbasisloc)
      maxlwf = max(maxlwf,nbands)
      maxnftl = max(maxnftl,maxnft)
      maxnfl = max(maxnfl,maxnf)
      maxnhijl = max(maxnhijl,maxnhij)
      maxnctl = max(maxnctl,maxnct)

C Work out range of bands for this Node
      nbandspernode = ((nbands-1)/Nodes)
      nremainder = nbands - nbandspernode*Nodes
      nbandsmin = Node*nbandspernode + min(nremainder,Node) + 1
      nbandsmax = nbandsmin + nbandspernode - 1
      if (Node.lt.nremainder) nbandsmax = nbandsmax + 1

      if (frstme) then
c        if (Node.eq.0) then
c          call prmem( 0, 'eandg', 'cttoc',   'i', maxnct*maxlwf  )
c          call prmem( 0, 'eandg', 'fttof',   'i', maxnft*maxoloc )
c          call prmem( 0, 'eandg', 'listct',  'i', maxnct*maxlwf  )
c          call prmem( 0, 'eandg', 'listf',   'i', maxnf*maxlwf   )
c          call prmem( 0, 'eandg', 'listft',  'i', maxnft*maxoloc )
c          call prmem( 0, 'eandg', 'listhij', 'i', maxnhij*maxlwf )
c          call prmem( 0, 'eandg', ' ',       ' ', 0              )
c        endif
        frstme = .false.
      endif

C .....................
C Set up index lists for sparse matrices ..................................
      if (iopt .eq. 0) then
        allocate(ind(nbasis))
        call memory('A','I',nbasis,'eandg')
        allocate(nindv(maxnf))
        call memory('A','I',maxnf,'eandg')
C GET Ct LISTS
        call ctrans1(nbasis,nbandsmin,nbandsmax,nbands,ncmax,maxnct,
     .               numc,listc,numct,listct,cttoc)

C GET F LISTS
        call axb_build(nbandsmin,nbandsmax,nbasis,maxnct,numct,listct,
     .                 nbasis,nbasis,nhmax,numh,listh,
     .                 ind,nindv,maxnf,numf,listf)
C GET Ft LISTS
        call ctrans2(nbands,nbasis,maxnf,maxnft,numf,listf,
     .              numft,listft,fttof,Node,Nodes)
C GET Hij LISTS
        call axb_build(nbandsmin,nbandsmax,nbasis,maxnf,numf,listf,
     .                 nbasis,nbands,ncmax,numc,listc,
     .                 ind,nindv,
     .                 maxnhij,numhij,listhij)
        call memory('D','I',size(nindv),'eandg')
        deallocate(nindv)
        call memory('D','I',size(ind),'eandg')
        deallocate(ind)
        goto 999
      endif
C.........................

C Calculate the energy at three points of the line, for the
C CG line minimization .....................................................
      if (iopt .eq. 1) then
        call ener3(c,grad,lam,eta,enum,h,s,
     .             nbasis,nbands,ncmax,maxnct,
     .             maxnf,nhmax,maxnhij,
     .             numc,listc,numct,listct,cttoc,numf,listf,
     .             numh,listh,numhij,listhij,
     .             e3,nbasisloc)
*       return
        goto 999
      endif
C.........................

C Calculate the energy and Gradient at current point .......................
      if (iopt .eq. 2) then
        call gradient(c,eta,enum,h,s,
     .                nbasis,nbands,ncmax,maxnct,
     .                maxnf,maxnft,nhmax,maxnhij,
     .                numc,listc,numct,listct,cttoc,numf,listf,
     .                numft,listft,
     .                numh,listh,numhij,listhij,
     .                grad,e,nbasisloc)
*       return
        goto 999
      endif
C.........................

C Calculate density matrix .................................................
C-JMS Modified denmat argument list
      if (iopt .eq. 3) then
        call denmat(c,eta,h,s,enum,
     .                nbasis,nbands,ncmax,maxnct,
     .                maxnf,maxnft,nhmax,
     .                numc,listc,numct,listct,cttoc,
     .                numf,listf,numft,listft,fttof,
     .                numh,listh,dm,edm,nbasisloc)
*       return
        goto 999
      endif
C.........................
      if (Node.eq.0) then
        stop 'Error in eandg: incorrect iopt'
      else
        stop
      endif
      
  999 continue
*     call timer('eandg',2)
      return
      end
