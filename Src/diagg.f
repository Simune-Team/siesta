C $Id: diagg.f,v 1.3 1999/05/18 17:17:32 ordejon Exp $

      subroutine diagg( nspin, no, maxuo, maxno, maxnd, 
     .                  numh, listh, numd, listd, H, S,
     .                  getD, fixspin, qtot, qs, temp, e1, e2,
     .                  eo, qo, Dnew, Enew, ef, efs, 
     .                  Haux, Saux, psi, aux )
C *********************************************************************
C Subroutine to calculate the eigenvalues and eigenvectors, density
C and energy-density matrices, and occupation weights of each 
C eigenvector, for given Hamiltonian and Overlap matrices (including
C spin polarization). Gamma-point version.
C Writen by J.Soler, August 1998.
C **************************** INPUT **********************************
C integer nspin               : Number of spin components (1 or 2)
C integer no                  : Number of basis orbitals
C integer maxuo               : First dimension of eo, qo, last of xij
C                               Must be at least max(indxuo)
C integer maxno               : Maximum number of orbitals interacting  
C                               with any orbital
C integer maxnd               : Maximum number of nonzero elements of 
C                               each row of density matrix
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
C *************************** OUTPUT **********************************
C real*8 eo(maxuo,nspn)          : Eigenvalues
C ******************** OUTPUT (only if getD=.true.) *******************
C real*8 qo(maxuo,nspn)          : Occupations of eigenstates
C real*8 Dnew(maxno,maxuo,nspin) : Output Density Matrix
C real*8 Enew(maxno,maxuo,nspin) : Output Energy-Density Matrix
C real*8 ef                      : Fermi energy
C real*8 efs(nspin)              : Fermi energy for each spin
C                                  (for fixed spin calculations)
C *************************** AUXILIARY *******************************
C real*8 Haux(no,no)      : Auxiliary space for the hamiltonian matrix
C real*8 Saux(no,no)      : Auxiliary space for the overlap matrix
C real*8 psi(no,no,nspin) : Auxiliary space for the eigenvectors
C real*8 aux(no*5)        : Extra auxiliary space
C *************************** UNITS ***********************************
C eo, Enew and ef returned in the units of H.
C *********************************************************************
      implicit none

      integer
     .  maxnd, maxno, maxuo, no, nspin

      integer 
     .  listh(maxno,no), numh(no),
     .  listd(maxnd,no), numd(no)

      double precision
     .  Dnew(maxnd,maxuo,nspin),
     .  e1, e2, ef, efs(nspin), Enew(maxnd,maxuo,nspin), 
     .  eo(maxuo,nspin),
     .  H(maxno,maxuo,nspin), qo(maxuo,nspin), qs(nspin), qtot, 
     .  S(maxno,maxuo), stepf, temp
     
      double precision
     .  Haux(no,no), Saux(no,no), psi(no,no,nspin), aux(no*5)

      logical
     .  fixspin, getD

      external
     .  fermid, fermispin, rdiag, stepf

C  Internal variables .............................................
      integer           ie, io, ispin, j, jo
      double precision  ee, pipj, qe, t
C  ....................

C Solve eigenvalue problem .........................................
      do ispin = 1,nspin
        do io = 1,no
          do jo = 1,no
            Saux(jo,io) = 0.d0
            Haux(jo,io) = 0.d0
          enddo
        enddo
        do io = 1,no
          do j = 1,numh(io)
            jo = listh(j,io)
            Saux(io,jo) = Saux(io,jo) + S(j,io)
            Haux(io,jo) = Haux(io,jo) + H(j,io,ispin)
          enddo
        enddo
        do io = 1,no
          do jo = 1,io-1
            Saux(jo,io) = 0.5d0 * ( Saux(jo,io) + Saux(io,jo) )
            Haux(jo,io) = 0.5d0 * ( Haux(jo,io) + Haux(io,jo) )
            Saux(io,jo) = Saux(jo,io)
            Haux(io,jo) = Haux(jo,io)
          enddo
        enddo
        call rdiag( Haux, Saux, no, no,
     .              eo(1,ispin), psi(1,1,ispin), aux )
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
        call fermid( nspin, nspin, 1, 1.d0, maxuo, no, eo, 
     .               temp, qtot, qo, ef )
      endif
C ....................

C Find weights for local density of states ............................
      if (e1 .lt. e2) then
*       e1 = e1 - ef
*       e2 = e2 - ef
        t = max( temp, 1.d-6 )
        do ispin = 1,nspin
          do io = 1,no
            qo(io,ispin) = ( stepf( (eo(io,ispin)-e2)/t ) -
     .                       stepf( (eo(io,ispin)-e1)/t ) ) / nspin
          enddo
        enddo
      endif
C ....................
      
C New density and energy-density matrices of unit-cell orbitals .......
      do ispin = 1,nspin
        do io = 1,no
          do j = 1,numd(io)
            Dnew(j,io,ispin) = 0.d0
            Enew(j,io,ispin) = 0.d0
          enddo
        enddo
      enddo
      do ispin = 1,nspin
        do ie = 1,no
          qe = qo(ie,ispin)
          ee = qo(ie,ispin) * eo(ie,ispin)
          do io = 1,no
            do j = 1,numd(io)
              jo = listd(j,io)
              pipj = psi(io,ie,ispin) * psi(jo,ie,ispin)
              Dnew(j,io,ispin) = Dnew(j,io,ispin) + qe * pipj
              Enew(j,io,ispin) = Enew(j,io,ispin) + ee * pipj
            enddo
          enddo
        enddo
      enddo
C ....................

      end


