C $Id: diagon.f,v 1.15 1999/05/18 17:17:32 ordejon Exp $

      subroutine diagon(no, nspin, maxspn, maxuo, maxno, maxnd, 
     .                  numh, listh, numd, listd, H, S,
     .                  qtot, fixspin, qs, temp, e1, e2,
     .                  gamma, xij, indxuo, nk, kpoint, wk,
     .                  eo, qo, Dnew, Enew, ef, efs)
C *********************************************************************
C Subroutine to calculate the eigenvalues and eigenvectors, density
C and energy-density matrices, and occupation weights of each 
C eigenvector, for given Hamiltonian and Overlap matrices (including
C spin polarization).
C Writen by J.Soler and P.Ordejon, August 1997.
C **************************** INPUT **********************************
C integer no                  : Number of basis orbitals
C integer nspin               : Spin polarization (1 or 2)
C integer maxspn              : Second dimension of eo and qo
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
C real*8  qtot                : Number of electrons in unit cell
C logical fixspin             : Fix the spin of the system?
C real*8  qs(nspin)           : Number of electrons in unit cell for each
C                               spin component (if fixed spin option is used)
C real*8  temp                : Electronic temperature 
C real*8  e1, e2              : Energy range for density-matrix states
C                               (to find local density of states)
C                               Not used if e1 > e2
C logical gamma               : Only gamma point?
C real*8  xij(3,maxno,maxuo)  : Vectors between orbital centers (sparse)
C                               (not used if only gamma point)
C integer indxuo(no)          : Index of equivalent orbital in unit cell
C                               Unit cell orbitals must be the first in
C                               orbital lists, i.e. indxuo.le.nuo, with
C                               nuo the number of orbitals in unit cell
C integer nk                  : Number of k points
C real*8  kpoint(3,nk)        : k point vectors
C real*8  wk(nk)              : k point weights (must sum one)
C *************************** OUTPUT **********************************
C real*8 eo(maxuo,maxspn,nk)   : Eigenvalues
C real*8 qo(maxuo,maxspn,nk)   : Occupations of eigenstates
C real*8 Dnew(maxno,maxuo,nspin) : Output Density Matrix
C real*8 Enew(maxno,maxuo,nspin) : Output Energy-Density Matrix
C real*8 ef                    : Fermi energy
C real*8 efs(nspin)            : Fermi energy for each spin
C                                (for fixed spin calculations)
C *************************** UNITS ***********************************
C xij and kpoint must be in reciprocal coordinates of each other.
C temp and H must be in the same energy units.
C eo, Enew and ef returned in the units of H.
C *********************************************************************
      implicit none

      integer
     .  maxnd, maxno, maxspn, maxuo, nk, no, nspin

      integer 
     .  indxuo(no), listh(maxno,no), numh(no),
     .  listd(maxnd,no), numd(no)

      double precision
     .  Dnew(maxnd,maxuo,nspin),
     .  e1, e2, ef, efs(nspin), 
     .  Enew(maxnd,maxuo,nspin), eo(maxuo,maxspn,nk),
     .  H(maxno,maxuo,nspin), kpoint(3,nk), qo(maxuo,maxspn,nk), 
     .  qs(nspin),
     .  qtot, S(maxno,maxuo), temp, wk(nk), xij(3,maxno,maxuo)

      logical
     .  fixspin, gamma

      external
     .  diagg, diagk, diag2g, diag2k, io_assign, io_close, prmem

C  Internal variables .............................................
      include 'diagon.h'

      double precision tiny
      parameter (tiny = 1.d-15)

      logical
     .  frstme, getD, overflow

      integer
     .  io, iu, iuo, muo(maxorb), naux, nhs, npsi, nuo

C     Common auxiliary arrays shared with routine bands
      double precision
     .    Dk(maxhs),   Ek(maxhs), 
     .  Haux(maxhs), Saux(maxhs), psi(maxpsi), aux(maxaux)
      equivalence (Dk,Haux), (Ek,Saux)
      common /diacom/ Haux, Saux, psi, aux

      save frstme
      data frstme /.true./
C  ....................

C Start time counter ................................................
      call timer( 'diagon', 1 )
C ......................

C Print array sizes .................................................
      if (frstme) then
        frstme = .false.
        call prmem( 0, 'diagon', 'aux',  'd', maxaux )
        call prmem( 0, 'diagon', 'Haux', 'd', maxhs  )
        call prmem( 0, 'diagon', 'psi',  'd', maxpsi )
        call prmem( 0, 'diagon', 'Saux', 'd', maxhs  )
        call prmem( 0, 'diagon', ' ',    ' ', 0      )     
      endif
C  ....................

C Find number of orbitals per unit cell and check argument sizes .....
      nuo = 0
      do io = 1,no
        nuo = max( nuo, indxuo(io) )
      enddo
      call chkdim( 'diagon', 'maxuo',  maxuo,  nuo,   1 )
      call chkdim( 'diagon', 'maxspn', maxspn, nspin, 1 )
C ....................

C Check internal dimensions ..........................................
      if (nspin.le.2 .and. gamma) then
        nhs  = nuo * nuo
        npsi = nuo * nuo * nspin
        naux = nuo * 5
      elseif (nspin.le.2 .and. .not.gamma) then
        nhs  = 2 * nuo * nuo
        npsi = 2 * nuo * nuo
        naux = nuo * 5
      elseif (nspin.eq.4) then
        nhs  = 2 * (2*nuo) * (2*nuo)
        npsi = 2 * (2*nuo) * (2*nuo)
        naux = (2*nuo) * 5
      else
        stop 'diagon: ERROR: incorrect value of nspin'
      endif
      overflow = .false.
      if (naux  .gt. maxaux)   overflow = .true.
      if (nhs   .gt. maxhs )   overflow = .true.
C  ....................

C Print good dimensions ..............................................
      if (overflow) then
        call io_assign(iu)
        open( iu, file='diagon.h', status='unknown' )
        write(iu,'(A)') 'C rdiagon: Dimension parameters for diagon'
        write(iu,'(6X,A)') 'integer maxaux,maxhs,maxorb,maxpsi'
        write(iu,'(6X,A,I10,A)') 'parameter ( maxaux  = ',  naux, ' )'
        write(iu,'(6X,A,I10,A)') 'parameter ( maxhs   = ',   nhs, ' )'
        write(iu,'(6X,A,I10,A)') 'parameter ( maxorb  = ',   nuo, ' )'
        write(iu,'(6X,A,I10,A)') 'parameter ( maxpsi  = ',  npsi, ' )'
        call io_close(iu)
        write(6,'(/,a,/)') 'diagon: BAD DIMENSIONS. RECOMPILE.'
        stop 'diagon: BAD DIMENSIONS. RECOMPILE.'
      endif
C .....................

C Check indxuo .......................................................
      do iuo = 1,nuo
        muo(iuo) = 0
      enddo
      do io = 1,no
        iuo = indxuo(io)
        if (indxuo(io).le.0 .or. indxuo(io).gt.no) then
          write(6,*) 'diagon: ERROR: invalid index: io, indxuo =',
     .      io, indxuo(io)
          stop 'diagon: ERROR: invalid indxuo'
        endif
        muo(iuo) = muo(iuo) + 1
      enddo
      do iuo = 1,nuo
        if (muo(iuo) .ne. muo(1)) then
          write(6,'(/,2a,3i6)') 'diagon: ERROR: inconsistent indxuo.',
     .     ' iuo, muo(iuo), muo(1) =', iuo, muo(iuo), muo(1)
          stop 'diagon: ERROR: inconsistent indxuo.'
        endif
      enddo
C ......................

C Call apropriate routine .............................................
      getD = .true.
      if (nspin.le.2 .and. gamma) then
        call diagg( nspin, nuo, maxuo, maxno, maxnd, 
     .              numh, listh, numd, listd, H, S,
     .              getD, fixspin, qtot, qs, temp, e1, e2,
     .              eo, qo, Dnew, Enew, ef, efs,
     .              Haux, Saux, psi, aux )
      elseif (nspin.le.2 .and. .not.gamma) then
        call diagk( nspin, nuo, no, maxspn, maxuo, maxno, maxnd, 
     .              numh, listh, numd, listd, H, S,
     .              getD, fixspin, qtot, qs, temp, e1, e2,
     .              xij, indxuo, nk, kpoint, wk,
     .              eo, qo, Dnew, Enew, ef, efs,
     .              Haux, Saux, psi, Dk, Ek, aux )
      elseif (nspin.eq.4 .and. gamma) then
        call diag2g( no, maxuo, maxno, maxnd, 
     .               numh, listh, numd, listd, H, S,
     .               getD, qtot, temp, e1, e2,
     .               eo, qo, Dnew, Enew, ef,
     .               Haux, Saux, psi, aux )
      elseif (nspin.eq.4 .and. .not.gamma) then
        call diag2k( nuo, no, maxuo, maxno, maxnd, 
     .               numh, listh, numd, listd, H, S,
     .               getD, qtot, temp, e1, e2,
     .               xij, indxuo, nk, kpoint, wk,
     .               eo, qo, Dnew, Enew, ef,
     .               Haux, Saux, psi, Dk, Ek, aux )
      endif
C ....................

C Stop time counter ...................................................
      call timer( 'diagon', 2 )
C .......................

      end


