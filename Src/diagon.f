      subroutine diagon(no, nspin, maxspn, maxo, maxno, maxnd, 
     .                  numh, listh, numd, listd, H, S,
     .                  qtot, temp, e1, e2,
     .                  xij, maxuo, indxuo, nk, kpoint, wk,
     .                  eo, qo, Dnew, Enew, ef)
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
C integer maxo                : Maximum number of basis  orbitals
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
C real*8  H(maxno,maxo,nspin) : Hamiltonian in sparse form
C real*8  S(maxno,maxo)       : Overlap in sparse form
C real*8  qtot                : Total number of electrons
C real*8  temp                : Electronic temperature 
C real*8  e1, e2              : Energy range for density-matrix states
C                               (to find local density of states)
C                               Not used if e1 > e2
C real*8  xij(3,maxno,maxuo)  : Vectors between orbital centers (sparse)
C                               (not used if only gamma point)
C integer maxuo               : First dimension of eo, qo, last of xij
C                               Must be at least max(indxuo)
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
C real*8 Dnew(maxo,maxo,nspin) : Output Density Matrix
C real*8 Enew(maxo,maxo,nspin) : Output Energy-Density Matrix
C real*8 ef                    : Fermi energy
C *************************** UNITS ***********************************
C xij and kpoint must be in reciprocal coordinates of each other.
C temp and H must be in the same energy units.
C eo, Enew and ef returned in the units of H.
C *********************************************************************
      implicit none

      integer
     .  maxo, maxnd, maxno, maxspn, maxuo, nk, no, nspin

      integer 
     .  indxuo(no), listh(maxno,no), numh(no),
     .  listd(maxnd,no), numd(no)

      double precision
     .  Dnew(maxnd,maxo,nspin), dot,
     .  e1, e2, ef, Enew(maxnd,maxo,nspin), eo(maxuo,maxspn,nk),
     .  H(maxno,maxo,nspin), kpoint(3,nk), qo(maxuo,maxspn,nk), qtot,
     .  S(maxno,maxo), stepf, temp, wk(nk), xij(3,maxno,maxuo)
     
      external
     .  cdiag, dot, fermid, io_assign, io_close, iodm, prmem, rdiag,
     .  stepf

C  Internal variables .............................................

      double precision tiny
      parameter (tiny = 1.d-15)

      include 'diagon.h'

      integer
     .  i1, i2, ie, ik, io, ispin, iu, iuo, j, jo, juo,
     .  muo(maxorb), naux, ncells, nuo, n1or2

      double precision
     .  aux(maxaux), deltaD, ee, Haux(max1or2,maxorb,maxorb), kxij,
     .  pipj, pipj1, pipj2, psi(max1or2,maxorb,maxorb,maxspin), 
     .  qe, Saux(max1or2,maxorb,maxorb), t

      logical
     .  frstme, gamma, overflow
C    .  found

      save frstme
      data frstme /.true./
C  ....................

C Start time counter ................................................
      call timer( 'diagon', 1 )
C ......................

C Print array sizes .................................................
      if (frstme) then
        frstme = .false.
        call prmem( 0, 'diagon', 'aux',  'd', maxaux                )
        call prmem( 0, 'diagon', 'Haux', 'd', max1or2*maxorb*maxorb )
        call prmem( 0, 'diagon', 'psi',  'd', max1or2*maxorb*
     .                                               maxorb*maxspin )
        call prmem( 0, 'diagon', 'Saux', 'd', max1or2*maxorb*maxorb )
        call prmem( 0, 'diagon', ' ',    ' ', 0                     )     
      endif
C  ....................

C Find if only gamma point is used .................................
      if (nk.eq.1 .and. abs(kpoint(1,1)).lt.tiny .and.
     .                  abs(kpoint(2,1)).lt.tiny .and.
     .                  abs(kpoint(3,1)).lt.tiny) then
        gamma = .true.
      else
        gamma = .false.
      endif
C ....................

C Find number of orbitals per unit cell and check argument sizes .....
      nuo = 0
      do io = 1,no
        nuo = max( nuo, indxuo(io) )
      enddo
      call chkdim( 'diagon', 'maxuo',  maxuo,  nuo,   1 )
      call chkdim( 'diagon', 'maxo',   maxo,   no,    1 )
      call chkdim( 'diagon', 'maxspn', maxspn, nspin, 1 )
C ....................

C Check internal dimensions ..........................................
      if (gamma) then
        n1or2 = 1
      else
        n1or2 = 2
      endif
      naux  = nuo*5
      overflow = .false.
      if (naux  .gt. maxaux)   overflow = .true.
      if (nuo   .gt. maxorb)   overflow = .true.
      if (nspin .gt. maxspin)  overflow = .true.
      if (n1or2 .ne. max1or2)  overflow = .true.
C  ....................

C Print good dimensions ..............................................
      if (overflow) then
        call io_assign(iu)
        open( iu, file='diagon.h', status='unknown' )
        write(iu,'(A)') 'C rdiagon: Dimension parameters for diagon'
        write(iu,'(6X,A)') 'integer maxaux,maxorb,maxspin,max1or2'
        write(iu,'(6X,A,I10,A)') 'parameter ( maxaux  = ',  naux, ' )'
        write(iu,'(6X,A,I10,A)') 'parameter ( maxorb  = ',   nuo, ' )'
        write(iu,'(6X,A,I10,A)') 'parameter ( maxspin = ', nspin, ' )'
        write(iu,'(6X,A,I10,A)') 'parameter ( max1or2 = ', n1or2, ' )'
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
          write(6,*) 'diagon: invalid index: io, indxuo =',
     .      io, indxuo(io)
          stop 'diagon: invalid indxuo'
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

C Find number of unit cells ........................................
      ncells = no / nuo
C ....................

C Impose the translational symmetry of the supercell ...............
C (possibly broken by the mesh)
      if (ncells .lt. 1) then
        do ispin = 1,nspin
          do io = nuo+1,no
            iuo = indxuo(io)
            do j = 1,numh(iuo)
              H(j,iuo,ispin) = H(j,iuo,ispin) + H(j,io,ispin)
            enddo
          enddo
          do io = 1,nuo
            do j = 1,numh(iuo)
              H(j,io,ispin) = H(j,io,ispin) / ncells
            enddo
          enddo
          do io = nuo+1,no
            iuo = indxuo(io)
            do j = 1,numh(iuo)
              H(j,io,ispin) = H(j,iuo,ispin)
            enddo
          enddo
        enddo
      endif
C ....................

C Solve eigenvalue problem .........................................
      if (gamma) then
        do ispin = 1,nspin
          do iuo = 1,nuo
            do juo = 1,nuo
              Saux(1,juo,iuo) = 0.d0
              Haux(1,juo,iuo) = 0.d0
            enddo
          enddo
          do io = 1,nuo
            do j = 1,numh(io)
              jo = listh(j,io)
              iuo = indxuo(io)
              juo = indxuo(jo)
              Saux(1,iuo,juo) = Saux(1,iuo,juo) + S(j,io)
              Haux(1,iuo,juo) = Haux(1,iuo,juo) + H(j,io,ispin)
            enddo
          enddo
          do iuo = 1,nuo
            do juo = 1,iuo-1
              Saux(1,juo,iuo) = 0.5d0 * ( Saux(1,juo,iuo) +
     .                                    Saux(1,iuo,juo) )
              Saux(1,iuo,juo) = Saux(1,juo,iuo)
              Haux(1,juo,iuo) = 0.5d0 * ( Haux(1,juo,iuo) +
     .                                    Haux(1,iuo,juo) )
              Haux(1,iuo,juo) = Haux(1,juo,iuo)
            enddo
          enddo
          call rdiag( Haux, Saux, nuo, maxorb,
     .                eo(1,ispin,1), psi(1,1,1,ispin), aux )
        enddo
      else
C       We use i1,i2 instead of 1,2 to avoid compilation warnings
        i1 = 1
        i2 = 2
        do ik = 1,nk
          do ispin = 1,nspin
            do iuo = 1,nuo
              do juo = 1,nuo
                Saux(i1,juo,iuo) = 0.d0
                Saux(i2,juo,iuo) = 0.d0
                Haux(i1,juo,iuo) = 0.d0
                Haux(i2,juo,iuo) = 0.d0
              enddo
            enddo
            do io = 1,nuo
              do j = 1,numh(io)
                jo = listh(j,io)
                iuo = indxuo(io)
                juo = indxuo(jo)
                kxij = dot( kpoint(1,ik), xij(1,j,io), 3 )
                Saux(i1,iuo,juo) = Saux(i1,iuo,juo) +
     .                             S(j,io) * cos(kxij)
                Saux(i2,iuo,juo) = Saux(i2,iuo,juo) +
     .                             S(j,io) * sin(kxij)
                Haux(i1,iuo,juo) = Haux(i1,iuo,juo) +
     .                             H(j,io,ispin) * cos(kxij)
                Haux(i2,iuo,juo) = Haux(i2,iuo,juo) +
     .                             H(j,io,ispin) * sin(kxij)
              enddo
            enddo
            do iuo = 1,nuo
              do juo = 1,iuo-1
                Saux(i1,juo,iuo) = 0.5d0 * ( Saux(i1,juo,iuo) +
     .                                       Saux(i1,iuo,juo) )
                Saux(i1,iuo,juo) = Saux(i1,juo,iuo)
                Saux(i2,juo,iuo) = 0.5d0 * ( Saux(i2,juo,iuo) -
     .                                       Saux(i2,iuo,juo) )
                Saux(i2,iuo,juo) = -Saux(i2,juo,iuo)
                Haux(i1,juo,iuo) = 0.5d0 * ( Haux(i1,juo,iuo) +
     .                                       Haux(i1,iuo,juo) )
                Haux(i1,iuo,juo) = Haux(i1,juo,iuo)
                Haux(i2,juo,iuo) = 0.5d0 * ( Haux(i2,juo,iuo) -
     .                                       Haux(i2,iuo,juo) )
                Haux(i2,iuo,juo) = -Haux(i2,juo,iuo)
              enddo
              Saux(i2,iuo,iuo) = 0.d0
              Haux(i2,iuo,iuo) = 0.d0
            enddo
            call cdiag( Haux, maxorb, Saux, maxorb, nuo,
     .                  eo(1,ispin,ik), psi(1,1,1,ispin),
     .                  maxorb, aux )
          enddo
        enddo
      endif
C ....................

C Find new Fermi energy and occupation weights ........................
      call fermid( nspin, maxspn, nk, wk, maxuo, nuo, eo, 
     .             temp, qtot/ncells, qo, ef )
C ....................

C Find weights for local density of states ............................
      if (e1 .lt. e2) then
*       e1 = e1 - ef
*       e2 = e2 - ef
        t = max( temp, 1.d-6 )
        do ik = 1,nk
          do ispin = 1,nspin
            do io = 1,nuo
              qo(io,ispin,ik) = wk(ik) * 
     .             ( stepf( (eo(io,ispin,ik)-e2)/t ) -
     .               stepf( (eo(io,ispin,ik)-e1)/t ) ) / nspin
            enddo
          enddo
        enddo
      endif
C ....................
      

c New density and energy-density matrices of unit-cell orbitals .......
      do ispin = 1,nspin
        do io = 1,nuo
          do j = 1,numd(io)
            Dnew(j,io,ispin) = 0.d0
            Enew(j,io,ispin) = 0.d0
          enddo
        enddo
      enddo

      if (gamma) then

        do ispin = 1,nspin
          do ie = 1,nuo
            qe = qo(ie,ispin,1)
            ee = qo(ie,ispin,1) * eo(ie,ispin,1)
            do io = 1,nuo
              do j = 1,numd(io)
                jo = listd(j,io)
                iuo = indxuo(io)
                juo = indxuo(jo)
                pipj = psi(1,iuo,ie,ispin) * psi(1,juo,ie,ispin)
                Dnew(j,io,ispin) = Dnew(j,io,ispin) + qe * pipj
                Enew(j,io,ispin) = Enew(j,io,ispin) + ee * pipj
              enddo
            enddo
          enddo
        enddo

      else
        i1 = 1
        i2 = 2

        do ik = 1,nk
          do ispin = 1,nspin

C           Find eigenvectors again (were stored only for one k point)
            if (nk .gt. 1) then
              do iuo = 1,nuo
                do juo = 1,nuo
                  Saux(i1,juo,iuo) = 0.d0
                  Saux(i2,juo,iuo) = 0.d0
                  Haux(i1,juo,iuo) = 0.d0
                  Haux(i2,juo,iuo) = 0.d0
                enddo
              enddo
              do io = 1,nuo
                do j = 1,numh(io)
                  jo = listh(j,io)
                  iuo = indxuo(io)
                  juo = indxuo(jo)
                  kxij = dot( kpoint(1,ik), xij(1,j,io), 3 )
                  Saux(i1,iuo,juo) = Saux(i1,iuo,juo) +
     .                               S(j,io) * cos(kxij)
                  Saux(i2,iuo,juo) = Saux(i2,iuo,juo) +
     .                               S(j,io) * sin(kxij)
                  Haux(i1,iuo,juo) = Haux(i1,iuo,juo) + 
     .                               H(j,io,ispin) * cos(kxij)
                  Haux(i2,iuo,juo) = Haux(i2,iuo,juo) + 
     .                               H(j,io,ispin) * sin(kxij)
                enddo
              enddo
              do iuo = 1,nuo
                do juo = 1,iuo-1
                  Saux(i1,juo,iuo) = 0.5d0 * ( Saux(i1,juo,iuo) +
     .                                         Saux(i1,iuo,juo) )
                  Saux(i1,iuo,juo) = Saux(i1,juo,iuo)
                  Saux(i2,juo,iuo) = 0.5d0 * ( Saux(i2,juo,iuo) -
     .                                         Saux(i2,iuo,juo) )
                  Saux(i2,iuo,juo) = -Saux(i2,juo,iuo)
                  Haux(i1,juo,iuo) = 0.5d0 * ( Haux(i1,juo,iuo) +
     .                                         Haux(i1,iuo,juo) )
                  Haux(i1,iuo,juo) = Haux(i1,juo,iuo)
                  Haux(i2,juo,iuo) = 0.5d0 * ( Haux(i2,juo,iuo) -
     .                                         Haux(i2,iuo,juo) )
                  Haux(i2,iuo,juo) = -Haux(i2,juo,iuo)
                enddo
                Saux(i2,iuo,iuo) = 0.d0
                Haux(i2,iuo,iuo) = 0.d0
              enddo
              call cdiag( Haux, maxorb, Saux, maxorb, nuo,
     .                    eo(1,ispin,ik), psi(1,1,1,ispin),
     .                    maxorb, aux )
            endif

C           Add contribution to density matrices of unit-cell orbitals
            call timer( 'diagon3', 1 )
            do io = 1,nuo
              do j = 1,numd(io)
                jo = listd(j,io)
                iuo = indxuo(io)
                juo = indxuo(jo)
                kxij = dot( kpoint(1,ik), xij(1,j,io), 3 )
                do ie = 1,nuo
                  pipj1 = psi(i1,iuo,ie,ispin) * psi(i1,juo,ie,ispin) +
     .                    psi(i2,iuo,ie,ispin) * psi(i2,juo,ie,ispin)
                  pipj2 = psi(i1,iuo,ie,ispin) * psi(i2,juo,ie,ispin) -
     .                    psi(i2,iuo,ie,ispin) * psi(i1,juo,ie,ispin)
                  deltaD = qo(ie,ispin,ik) * ( pipj1 * cos(kxij) -
     .                                         pipj2 * sin(kxij) )
                  Dnew(j,io,ispin) = Dnew(j,io,ispin) + deltaD
                  Enew(j,io,ispin) = Enew(j,io,ispin) + 
     .                                deltaD * eo(ie,ispin,ik)
                enddo
              enddo
            enddo
            call timer( 'diagon3', 2 )

          enddo
        enddo

      endif
C ....................

C Copy density matrices to the remaining orbitals .....................
      do ispin = 1,nspin
        do io = nuo+1,no
          iuo = indxuo(io)
          do j = 1,numd(iuo)
            Dnew(j,io,ispin) = Dnew(j,iuo,ispin)
            Enew(j,io,ispin) = Enew(j,iuo,ispin)
          enddo
        enddo
      enddo
C ....................

C Save DM to disk .....................................................
*     call iodm( 'write', maxnd, maxo, no, nspin,
*    .           numd, listd, Dnew, found )
C .......................

C Stop time counter ...................................................
      call timer( 'diagon', 2 )
C .......................

      end


