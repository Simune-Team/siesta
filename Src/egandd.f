      subroutine eandg(iopt,eta,enum,lam,
     .                 nhmax,numh,listh,ncmax,numc,listc,h,s,c,
     .                 nbasis,nbands,
     .                 e3,e,grad,dm,edm)
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

      implicit none

      include 'ordern.h'

      integer
     .  iopt,nbands,nbasis,ncmax,nhmax

      integer
     .  listc(ncmax,nbasis),listh(nhmax,nbasis),
     .  numc(nbasis),numh(nbasis)

      double precision
     .  c(ncmax,nbasis),e,e3(3),
     .  dm(nhmax,nbasis),edm(nhmax,nbasis),enum,eta,
     .  h(nhmax,nbasis),grad(ncmax,nbasis),lam,s(nhmax,nbasis)

      external
     .  axb_build,ctrans,ener3,gradient,ind_gf

C Internal variables .....................................................
      integer i

      integer
     .  cttoc(maxnct,maxlwf),fttof(maxnft,maxo),
     .  indgf(maxnc,maxo),listct(maxnct,maxlwf),
     .  listf(maxnf,maxlwf),listft(maxnft,maxo),
     .  listhij(maxnhij,maxlwf),
     .  numct(maxlwf),numf(maxlwf),
     .  numft(maxo),numhij(maxlwf),
     .  ind(maxo),nindv(maxnhf)

      double precision
     .  f(maxnf*maxlwf),fs(maxnf*maxlwf)
     
      logical
     .  frstme

      save
     .  cttoc,frstme,fttof,indgf,listct,listf,listft,listhij,
     .  numct,numf,numft,numhij
     
      data frstme /.true./
C .....................

*     call timer('eandg',1)

      if (frstme) then
        call prmem( 0, 'eandg', 'cttoc',   'i', maxnct*maxlwf  )
        call prmem( 0, 'eandg', 'f',       'd', maxnf*maxlwf   )
        call prmem( 0, 'eandg', 'fs',      'd', maxnf*maxlwf   )
        call prmem( 0, 'eandg', 'fttof',   'i', maxnft*maxo    )
        call prmem( 0, 'eandg', 'indgf',   'i', maxnc*maxo     )
        call prmem( 0, 'eandg', 'listct',  'i', maxnct*maxlwf  )
        call prmem( 0, 'eandg', 'listf',   'i', maxnf*maxlwf   )
        call prmem( 0, 'eandg', 'listft',  'i', maxnft*maxo    )
        call prmem( 0, 'eandg', 'listhij', 'i', maxnhij*maxlwf )
        call prmem( 0, 'eandg', ' ',       ' ', 0              )
        frstme = .false.
      endif


C Check matrix dimensions .................................................
      call chkdim('eandg','ncmax',ncmax,maxnc,0)
C .....................
C Set up index lists for sparse matrices ..................................
      if (iopt .eq. 0) then
C GET Ct LISTS
        call ctrans(nbasis,nbands,maxnc,maxnct,numc,listc,
     .              numct,listct,cttoc)

        do i = 1,nbands
        enddo
C GET F LISTS
        call axb_build(nbands,nbasis,maxnct,numct,listct,
     .                 nbasis,nbasis,nhmax,numh,listh,
     .                 ind,nindv,
     .                 maxnf,numf,listf)
        do i = 1,nbands
        enddo
C GET indgf map
        call ind_gf(nbasis,nbands,maxnc,maxnf,numc,listc,numf,listf,
     .              indgf)
C GET Ft LISTS
        call ctrans(nbands,nbasis,maxnf,maxnft,numf,listf,
     .              numft,listft,fttof)
        do i = 1,nbasis
        enddo
C GET Hij LISTS
        call axb_build(nbands,nbasis,maxnf,numf,listf,
     .                 nbasis,nbands,maxnc,numc,listc,
     .                 ind,nindv,
     .                 maxnhij,numhij,listhij)
        do i = 1,nbands
        enddo
*       return
        goto 999
      endif
C.........................

C Calculate the energy at three points of the line, for the
C CG line minimization .....................................................
      if (iopt .eq. 1) then
        call ener3(c,grad,lam,eta,enum,h,s,
     .             nbasis,nbands,maxnc,maxnct,
     .             maxnf,nhmax,maxnhij,
     .             numc,listc,numct,listct,cttoc,numf,listf,
     .             numh,listh,numhij,listhij,
     .             e3)
*       return
        goto 999
      endif
C.........................

C Calculate the energy and Gradient at current point .......................
      if (iopt .eq. 2) then
        call gradient(c,eta,enum,h,s,
     .                nbasis,nbands,maxnc,maxnct,
     .                maxnf,maxnft,nhmax,maxnhij,
     .                numc,listc,numct,listct,cttoc,numf,listf,
     .                numft,listft,fttof,
     .                numh,listh,numhij,listhij,indgf,f,fs,
     .                grad,e)
*       return
        goto 999
      endif
C.........................

C Calculate density matrix .................................................
C-JMS Modified denmat argument list
      if (iopt .eq. 3) then
        call denmat(c,eta,h,s,enum,
     .                nbasis,nbands,maxnc,maxnct,
     .                maxnf,maxnft,nhmax,maxnhij,
     .                numc,listc,numct,listct,cttoc,
     .                numf,listf,numft,listft,fttof,
     .                numh,listh,numhij,listhij,f,fs,
     .                dm,edm)
*       return
        goto 999
      endif
C.........................
      stop 'Error in eandg: incorrect iopt'
      
  999 continue
*     call timer('eandg',2)
      end
