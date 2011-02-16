! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      subroutine eandg(iopt,eta,qs,lam,nhmax,numh,listh,listhptr,
     .                 ncmax,numc,listc,h,s,c,nbasis,nbasisloc,nbands,
     .                 e3,e,grad,dm,edm,nbasisCloc,nspin,Node,Nodes)
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
C real*8 qs(nspin)         : Total number of electrons
C real*8 lam               : Length of step for line minimization
C integer nhmax            : First dimension of listh and H, and maximum
C                            number of nonzero elements of each row of H
C integer numh(nbasis)     : Control vector of H matrix
C                            (number of <>0 element of each row)
C integer listh(nhmax)     : Control vector of H matrix
C                            (list of <>0 element of each row)
C integer listhptr(nbasis) : Control vector of H matrix
C                            (pointer to start of row in listh/H)
C integer ncmax            : First dimension of listc and C, and maximum
C                            number of nonzero elements of each row of C
C integer numc(nbasis)     : Control vector of C matrix
C                            (number of <>0 element of each row)
C integer listc(ncmax,nbasis): Control vector of C matrix
C                            (list of <>0 element of each row)
C real*8 h(nhmax)          : Hamiltonian in sparse form
C real*8 s(nhmax)          : Overlap in sparse form
C real*8 c(ncmax,nbasis)   : Current point (wave func. coeffs. in sparse)
C integer nbasis           : Global number of atomic orbitals
C integer nbasisloc        : Local number of atomic orbitals
C integer nbands           : Number of Localized Wave Functions
C ******** INPUT OR OUTPUT (DEPENDING ON ARGUMENT IOPT) ********************
C real*8 grad(ncmax,nbasis)  : Gradient of the functional
C                            (input if iopt = 1)
C                            (output if iopt = 2)
C ***************************** OUTPUT *************************************
C real*8 e(3)              : Value of the energy in three points C+LAM_i*GRAD
C real*8 e                 : Value of the energy at point C
C real*8 dm(nhmax)         : Density matrix in sparse form
C real*8 edm(nhmax)        : Energy density matrix in sparse form
C **************************** BEHAVIOUR ***********************************
C The overlap matrix 'o' must be in the same sparse format as the 
C Hamiltonian matrix 'h', even if the overlap is more sparse than h
C (as due to the KB projectors, for instance). It will, in general,
C contain some zeros, therefore.
C **************************************************************************

      use precision
      use alloc
      use fdf
      use on_core
      use on_main,   only : nbL2G, nbG2L, nbandsloc, ncP2T
      use sys,       only : die
#ifdef MPI
      use globalise, only : setglobaliseB, setglobalise, setglobaliseF 
      use globalise, only : globalinitB
#endif

      implicit none

      integer
     .  iopt,nbands,nbasis,nbasisloc,nbasisCloc,ncmax,nhmax,nspin,
     .  Node,Nodes

      integer
     .  listc(ncmax,nbasisCloc),listh(nhmax),listhptr(nbasis),
     .  numc(nbasis),numh(nbasis)

      real(dp) ::
     .  c(ncmax,nbasisCloc,nspin),e(nspin),e3(3,nspin),dm(nhmax,nspin),
     .  edm(nhmax,nspin),qs(nspin),eta(nspin),h(nhmax,nspin),
     .  grad(ncmax,nbasisCloc,nspin),lam,s(nhmax)

      external
     .  axb_build1,axb_build2,ctrans1,ctrans2,ener3,gradient

C Internal variables .....................................................
      integer, save :: maxnft  = 1
      integer, save :: maxnct  = 1
      integer, save :: maxnf   = 1
      integer, save :: maxnhij = 1

      integer       :: i
      integer       :: m, mt, n

      logical, pointer, save :: lNeeded(:)
      logical,          save :: frstme = .true.
      logical,          save :: frstme2 = .true.
      logical,          save :: LoMem = .false.

C .....................

*     call timer('eandg',1)

      if (frstme) then
C Nullify pointers
        nullify(cttoc)
        nullify(fttof)
        nullify(listct)
        nullify(listf)
        nullify(listft)
        nullify(listhij)
        nullify(numct)
        nullify(numf)
        nullify(numft)
        nullify(numhij)
        nullify(indon)
        nullify(nindv)
        nullify(nbL2G)
        nullify(nbG2L)
        nullify(f)
        nullify(fs)

C Set flag for memory algorithm choice
        LoMem = fdf_boolean('ON.LowerMemory',.false.)

        frstme = .false.

      endif

      if (iopt.eq.0) then
C Find nbandsloc and pointers
        call re_alloc( nbL2G, 1, nbands, 'nbL2G', 'eandg' )
        call re_alloc( nbG2L, 1, nbands, 'nbG2L', 'eandg' )
        call re_alloc( lNeeded, 1, nbands, 'lNeeded', 'eandg' )
        do i = 1,nbands
          lNeeded(i) = .false.
        enddo
        do m = 1,nbasisloc
          mt = ncP2T(m)
          do n = 1,numc(mt)
            lNeeded(listc(n,mt)) = .true.
          enddo
        enddo
        nbandsloc = 0
        nbG2L(1:nbands) = 0
        do i = 1,nbands
          if (lNeeded(i)) then
            nbandsloc = nbandsloc + 1
            nbL2G(nbandsloc) = i
            nbG2L(i) = nbandsloc
          endif
        enddo
        deallocate(lNeeded)

C Allocation of memory
        call re_alloc(numct,1,nbandsloc,'numct','eandg',.true.,.false.)
        call re_alloc(numf,1,nbandsloc,'numf','eandg',.true.,.false.)
        call re_alloc(numft,1,nbasisCloc,'numft','eandg',.true.,
     &               .false.)
        call re_alloc(numhij,1,nbandsloc,'numhij','eandg',.true.,
     &               .false.)
        call re_alloc(indon,1,nbasis,'indon','eandg',.true.,.false.)
        call re_alloc(nindv,1,nbasis,'nindv','eandg',.true.,.false.)
        call re_alloc(cttoc,1,maxnct,1,nbandsloc,'cttoc','eandg',
     &               .true.,.false.)
        call re_alloc(fttof,1,maxnft,1,nbasisCloc,'fttof','eandg',
     &               .true.,.false.)
        call re_alloc(listct,1,maxnct,1,nbandsloc,'listct','eandg',
     &               .true.,.false.)
        call re_alloc(listf,1,maxnf,1,nbandsloc,'listf','eandg',
     &               .true.,.false.)
        call re_alloc(listft,1,maxnft,1,nbasisCloc,'listft','eandg',
     &               .true.,.false.)
        call re_alloc(listhij,1,maxnhij,1,nbandsloc,'listhij','eandg',
     &               .true.,.false.)
        call re_alloc(f,1,maxnf,1,nbandsloc,'f','eandg',.true.,.false.)
        call re_alloc(fs,1,maxnf,1,nbandsloc,'fs','eandg',.true.,
     &               .false.)

C Set up index lists for sparse matrices ..................................

C GET Ct LISTS
        call ctrans1(nbasisCloc,maxnct)

C GET F LISTS
        call axb_build1(nbasis,maxnct,numct,listct,nhmax,numh,listh,
     .                  listhptr,indon,maxnf)

C Allocate arrays that depend on variables set above
        call re_alloc( listf, 1, maxnf, 1, nbandsloc, 'listf', 'eandg' )

C GET Ft LISTS
        call ctrans2(nbasisCloc,maxnf,maxnft,numf,listf)

C GET Hij LISTS
        call axb_build2(nbands,maxnf,numf,listf,
     .                  nbasisCloc,ncmax,numc,listc,indon,maxnhij)

C Allocate arrays that depend on variables set above
        call re_alloc( f, 1, maxnf, 1, nbandsloc, 'f','eandg' )
        call re_alloc( fs, 1, maxnf, 1, nbandsloc, 'fs','eandg' )

#ifdef MPI
C Set up communication information
        call setglobalise(nbasis,nbasisloc,nbasisCloc,nhmax,numh,listh,
     .                    listhptr,Node,Nodes)
        call setglobaliseB(nbands,Node)
        call setglobaliseF(nbasis,nbands,nbasisCloc,maxnft,Node)
#endif
        goto 999
      endif
C.........................

      if (frstme2) then
#ifdef MPI
C Initialise communication arrays - call here to avoid increasing size of arrays later
        if (LoMem) then
          call globalinitB(1)
        else
          call globalinitB(3)
        endif
#endif
        frstme2 = .false.
      endif

C Calculate the energy at three points of the line, for the
C CG line minimization .....................................................
      if (iopt .eq. 1) then
        if (LoMem) then
          call ener3lomem(c,grad,lam,eta,qs,h,s,nbasis,nbasisloc,nbands,
     .                    ncmax,maxnct,maxnf,nhmax,maxnhij,numc,listc,
     .                    numct,listct,cttoc,numf,listf,numh,listh,
     .                    listhptr,numhij,listhij,e3,nbasisCloc,nspin,
     .                    Node)
        else
          call ener3(c,grad,lam,eta,qs,h,s,nbasis,nbasisloc,nbands,
     .               ncmax,maxnct,maxnf,nhmax,maxnhij,numc,listc,numct,
     .               listct,cttoc,numf,listf,numh,listh,listhptr,
     .               numhij,listhij,e3,nbasisCloc,nspin,Node)
        endif
        goto 999
      endif
C.........................

C Calculate the energy and Gradient at current point .......................
      if (iopt .eq. 2) then
        if (LoMem) then
          call gradientlomem(c,eta,qs,h,s,nbasis,nbasisloc,nbands,ncmax,
     .                  maxnct,maxnf,maxnft,nhmax,maxnhij,numc,
     .                  listc,numct,listct,cttoc,numf,listf,numft,
     .                  listft,fttof,numh,listh,listhptr,numhij,
     .                  listhij,f,fs,grad,e,nbasisCloc,nspin,Node)
        else
          call gradient(c,eta,qs,h,s,nbasis,nbasisloc,nbands,ncmax,
     .                  maxnct,maxnf,maxnft,nhmax,maxnhij,numc,
     .                  listc,numct,listct,cttoc,numf,listf,numft,
     .                  listft,fttof,numh,listh,listhptr,numhij,
     .                  listhij,f,fs,grad,e,nbasisCloc,nspin,Node)
        endif
        goto 999
      endif
C.........................

C Calculate density matrix .................................................
C-JMS Modified denmat argument list
      if (iopt .eq. 3) then
        if (LoMem) then
          call denmatlomem(c,eta,h,s,qs,nbasis,nbasisloc,nbands,ncmax,
     .                     maxnct,maxnf,maxnft,nhmax,maxnhij,numc,listc,
     .                     numct,listct,cttoc,numf,listf,numft,listft,
     .                     numh,listh,listhptr,numhij,listhij,dm,edm,
     .                     nbasisCloc,nspin,Node)
        else
          call denmat(c,eta,h,s,qs,nbasis,nbasisloc,nbands,ncmax,
     .                maxnct,maxnf,maxnft,nhmax,maxnhij,numc,listc,
     .                numct,listct,cttoc,numf,listf,numft,listft,
     .                numh,listh,listhptr,numhij,listhij,dm,edm,
     .                nbasisCloc,nspin,Node)
        endif
        goto 999
      endif
C.........................
      call die('Error in eandg: incorrect iopt')
      
  999 continue
*     call timer('eandg',2)
      end
