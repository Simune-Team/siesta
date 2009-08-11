! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      module writewave

      USE precision, only : dp
      use alloc,     only : re_alloc, de_alloc
      implicit none

      private

      integer, save                   :: nwf
      integer, save, pointer          :: nwflist(:)
      integer, save, pointer          :: iwf(:,:)

      integer, public, save           :: nwk
      real(dp), public, pointer, save :: wfk(:,:)
      logical, public, save           :: gamma_wavefunctions

      logical,    save :: wwf = .false.

      public :: setup_wf_kpoints, wwave, writew

      CONTAINS

      subroutine setup_wf_kpoints()
      USE alloc, only: re_alloc
      USE sys, only: die
      USE atomlist, only: no_u

      implicit none

      ! Find number of k-points for wavefunction printout

      nullify(wfk)

      nwk = 0
      call initwave( no_u, nwk, wfk )

      if (nwk .eq. 0) then
         gamma_wavefunctions = .true.
      else if (nwk .eq. 1) then
         if (dot_product(wfk(:,1),wfk(:,1)) < 1.0e-20_dp) then
            gamma_wavefunctions = .true.
         else
            gamma_wavefunctions = .false.
         endif
      else  
         gamma_wavefunctions = .false.
      endif

      end subroutine setup_wf_kpoints

!-----------------------------------------------------------------------

      subroutine initwave( norb, nk, kpoint )
C *********************************************************************
C Finds k-points for wavefunction printout
C Based on initband routine by J.Soler
C Written by P. Ordejon, June 2003
C **************************** INPUT **********************************
C integer norb           : Number of orbitals
C *************************** OUTPUT **********************************
C integer nk             : Number k points to compute wavefunctions
C real*8  kpoint(3,maxk) : k point vectors
C *************************** UNITS ***********************************
C Lengths in atomic units (Bohr).
C k vectors in reciprocal atomic units.
C ***************** BEHAVIOUR *****************************************
C - If nk=0 on input, k-points are read from labels WaveFuncKPoints and 
C   WaveFuncKPointsScale from the input fdf data file. If these labels 
C   are not present, it returns with nk=0.
C - Allowed values for WaveFuncKPointsScale are ReciprocalLatticeVectors 
C   and pi/a (default). If another value is given, it returns with nk=0
C   after printing a warning.
C - If nk>maxk, k-points and wavefunctions are not calculated and no 
C   warning is printed before return
C ***************** USAGE *********************************************
C Example of fdf wavefunction k-points specification for an FCC lattice.
C 
C     WaveFuncKPointsScale  pi/a
C     %block WaveFuncKPoints              # These are comments
C     0.000  0.000  0.000  from 1 to 10   # eigenstates 1-10 of Gamma
C     2.000  0.000  0.000  1 3 5          # eigenstates 1,3,5 of X
C     1.500  1.500  1.500                 # all eigenstates of K
C     %endblock WaveFuncKPoints
C
C If only given points (not lines) are desired, simply specify 1 as 
C the number of points along the line.
C *********************************************************************
C
C  Modules
C
      use precision
      use parallel,     only : Node
      use fdf
      use sys,          only : die
      use alloc,        only : re_alloc
      implicit          none

      integer           :: norb, nk
      real(dp), pointer :: kpoint(:,:)
      character         :: paste*30
      external          :: paste, redcel, memory
C *********************************************************************

      character 
     .  name1*10, name2*10, name3*10,
     .  scale*30, msg*120

      logical       :: outlng
      logical       :: WaveFuncPresent

      integer
     .  i, ix, iw, iw1, iw2, iw3,
     .  mscell(3,3), ni, nn, nr, nv

      real(dp)
     .  alat, caux(3,3), pi,
     .  rcell(3,3), ucell(3,3)

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline
 
C Start time counter 
      call timer( 'writewave', 1 )

C Find if there are k-points data
      outlng = fdf_boolean('LongOutput', .false.)
      wwf = fdf_boolean('WriteWaveFunctions',outlng)

C Find if there are k-points data
      WaveFuncPresent = fdf_block('WaveFuncKPoints',bfdf)
      if ( WaveFuncPresent ) then

C Find lattice constant
        alat = fdf_physical( 'LatticeConstant', 0.d0, 'Bohr' )
        if (alat .eq. 0.d0) then
          call die('writewave: ERROR: Lattice constant required')
        endif

C Find scale used in k point data
        scale = fdf_string( 'WaveFuncKPointsScale', 'pi/a' )
        if (scale .eq. 'pi/a') then
          pi = 4.d0 * atan(1.d0)
        elseif (scale .eq. 'ReciprocalLatticeVectors') then
          call redcel( alat, caux, ucell, mscell )
          call reclat( ucell, rcell, 1 )
        elseif (Node.eq.0) then
          write(6,'(a,/,2a,/,a)')
     .    'writewave: WARNING: Invalid value for WaveFuncKPointsScale',
     .    'writewave: Allowed values are pi/a and',
     .              ' ReciprocalLatticeVectors',
     .    'writewave: No band calculation performed'
        endif

C Find max number of k points and max number of bands per k-point
C Loop on data lines
        nk = 0
        nwf = 0
        do while(fdf_bline(bfdf,pline))
C Read and parse data line
          nn = fdf_bnnames(pline)
          nv = fdf_bnvalues(pline)
          ni = fdf_bnintegers(pline)
          nr = fdf_bnreals(pline)

C Check if syntax line is correct
          if (nv .ge. 3) then

C Check syntax
            if (nr .ne. 3) then
              write(msg,'(a,/,a)')
     .          'writewave: syntax ERROR in %block WaveFuncKPoints:',
     .          trim(fdf_getline(bfdf%mark))
              call die(msg)
            endif


C Add this point to total number of k points
            nk = nk + 1

C Find which eigenvectors should be printed
C Store indexes of wave functions to printout

C #X.XXX #Y.YYY #Z.ZZZ from #start to #end [step] #step
            if (nn .ge. 1) then

C Check that line contains 'from', 'to' and maybe 'step'
              if ((nn .ne. 2) .and. (nn .ne. 3)) then
                write(msg,'(a,/,a)')
     .            'writewave: syntax ERROR in %block WaveFuncKPoints:',
     .            trim(fdf_getline(bfdf%mark))
                call die(msg)
              endif
              name1 = fdf_bnames(pline,1)
              name2 = fdf_bnames(pline,2)
              if (nn .eq. 3) name3 = fdf_bnames(pline,3)
              if (name1 .ne. 'from' .or. name2 .ne. 'to') then
                write(msg,'(a,/,a)')
     .            'writewave: syntax ERROR in %block WaveFuncKPoints:',
     .            trim(fdf_getline(bfdf%mark))
                call die(msg)
              endif
              if (nn .eq. 3) then
                if (name3 .ne. 'step') then
                  write(msg,'(a,/,a)')
     .             'writewave: syntax ERROR in %block WaveFuncKPoints:',
     .             trim(fdf_getline(bfdf%mark))
                  call die(msg)
                endif
              endif

              iw1 = fdf_bintegers(pline,1)
              iw2 = fdf_bintegers(pline,2)
              if (iw1 .lt. 0) iw1 = norb + iw1 + 1
              if (iw2 .lt. 0) iw2 = norb + iw2 + 1
              if (nn .eq. 3) then
                iw3 = abs(fdf_bintegers(pline,3))
              else
                iw3 = 1
              endif
              ni = 0
              do iw = min(iw1,iw2), max(iw1,iw2), iw3
                ni = ni + 1
              enddo
              nwf = max(nwf,ni)

C #X.XXX #Y.YYY #Z.ZZZ #eigen1 [#eigen2 ...]
            elseif (ni .ne. 0) then
              nwf = max(nwf,ni)

C #X.XXX #Y.YYY #Z.ZZZ
            elseif (ni .eq. 0) then
              nwf = max(nwf,norb)
            endif
          else
C Bad syntax
            write(msg,'(a,/,a)')
     .        'writewave: syntax ERROR in %block WaveFuncKPoints:',
     .        trim(fdf_getline(bfdf%mark))
            call die(msg)
          endif
        enddo


C Allocate nwflist, iwf and kpoints structures according to nk and nwf
        call re_alloc( nwflist, 1, nk, 'nwflist', 'initwave' )
        call re_alloc( iwf, 1, nk, 1, nwf, 'iwf', 'initwave' )

        call re_alloc( kpoint, 1, 3, 1, nk, 'kpoint', 'initwave' )

C Fill nwflist, iwf and kpoints structures
C Loop on data lines
        nk = 0
        call fdf_brewind(bfdf)
        do while(fdf_bline(bfdf,pline))
C Read and parse data line
          nn = fdf_bnnames(pline)
          nv = fdf_bnvalues(pline)
          ni = fdf_bnintegers(pline)
          nr = fdf_bnreals(pline)

C Add this point to total number of k points
          nk = nk + 1

C Find coordinates of k point
          if (scale .eq. 'pi/a') then
            kpoint(1,nk) = fdf_bvalues(pline,1) * pi / alat
            kpoint(2,nk) = fdf_bvalues(pline,2) * pi / alat
            kpoint(3,nk) = fdf_bvalues(pline,3) * pi / alat
          elseif (scale .eq. 'ReciprocalLatticeVectors') then
            do ix = 1,3
              kpoint(ix,nk) = rcell(ix,1) * fdf_bvalues(pline,1) +
     .                        rcell(ix,2) * fdf_bvalues(pline,2) +
     .                        rcell(ix,3) * fdf_bvalues(pline,3)
            enddo
          endif


C Find which eigenvectors should be printed
C Store indexes of wave functions to printout

C #X.XXX #Y.YYY #Z.ZZZ from #start to #end [step] #step
          if (nn .ge. 2) then

            name1 = fdf_bnames(pline,1)
            name2 = fdf_bnames(pline,2)
            if (nn .eq. 3) name3 = fdf_bnames(pline,3)

            iw1 = fdf_bintegers(pline,1)
            iw2 = fdf_bintegers(pline,2)
            if (iw1 .lt. 0) iw1 = norb + iw1 + 1
            if (iw2 .lt. 0) iw2 = norb + iw2 + 1
            if (nn .eq. 3) then
              iw3 = abs(fdf_bintegers(pline,3))
            else
              iw3 = 1
            endif
            ni = 0
            do iw = min(iw1,iw2), max(iw1,iw2), iw3
              ni = ni + 1
              if (iw .lt. 0) then
                iwf(nk,ni) = norb + iw + 1
              else
                iwf(nk,ni) = iw
              endif
            enddo

C #X.XXX #Y.YYY #Z.ZZZ #eigen1 [#eigen2 ...]
          elseif (ni .ne. 0) then
            do i= 1, ni
              iw = fdf_bintegers(pline,i)
              if (iw .lt. 0) then
                iwf(nk,i) = norb + iw + 1
              else
                iwf(nk,i) = iw
              endif
            enddo

C #X.XXX #Y.YYY #Z.ZZZ
          elseif (ni .eq. 0) then
            ni = norb
            do i= 1, ni
              iwf(nk,i) = i
            enddo
          endif

          nwflist(nk) = ni
        enddo

      endif

      call timer( 'writewave', 2 )

      end subroutine initwave


      subroutine wwave( no, nspin, maxspn, maxo, maxuo, maxnh, maxk,
     .                  numh, listhptr, listh, H, S, ef, xij, indxuo,
     .                  nk, kpoint, nuotot, gamma, occtol)
C *********************************************************************
C Finds wavefunctions at selected k-points.
C Written by P. Ordejon, June 2003
C from routine 'bands' written by J.M.Soler
C **************************** INPUT **********************************
C integer no                  : Number of basis orbitals
C integer nspin               : Number of spin components
C integer maxspn              : Second dimension of ek
C integer maxo                : First dimension of ek
C integer maxuo               : Second dimension of H and S
C integer maxnh               : Maximum number of orbitals interacting  
C                               with any orbital
C integer maxk                : Last dimension of kpoint and ek
C integer numh(nuo)           : Number of nonzero elements of each row 
C                               of hamiltonian matrix
C integer listhptr(nuo)       : Pointer to start of each row of the
C                               hamiltonian matrix
C integer listh(maxlh)        : Nonzero hamiltonian-matrix element  
C                               column indexes for each matrix row
C real*8  H(maxnh,nspin)      : Hamiltonian in sparse form
C real*8  S(maxnh)            : Overlap in sparse form
C real*8  ef                  : Fermi energy
C real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
C                               (not used if only gamma point)
C integer indxuo(no)          : Index of equivalent orbital in unit cell
C                               Unit cell orbitals must be the first in
C                               orbital lists, i.e. indxuo.le.nuo, with
C                               nuo the number of orbitals in unit cell
C integer nk                  : Number of band k points
C real*8  kpoint(3,maxk)      : k point vectors
C integer nuotot              : Total number of orbitals in unit cell
C logical gamma               : Indicates if this is a Gamma point run
C real*8  occtol              : Occupancy threshold for DM build
C *************************** OUTPUT **********************************
C None; output is dumped to wave functions file SystemLabel.WFS
C *************************** UNITS ***********************************
C Lengths in atomic units (Bohr).
C k vectors in reciprocal atomic units.
C Energies in Rydbergs.
C
C  Modules
C
      use precision
      use parallel,    only : Node, Nodes
      use fdf
      use files,       only : slabel, label_length
      use densematrix, only : Haux, Saux, psi
      use alloc
      use sys,         only : die
#ifdef MPI
      use parallelsubs, only : GetNodeOrbs
#endif

      implicit          none

      integer           maxk, maxnh, maxspn, maxo, maxuo, nk, no, nspin,
     .                  nuotot, indxuo(no), listh(maxnh), numh(*), 
     .                  listhptr(*)
      logical           gamma
      real(dp)          ef, H(maxnh,nspin), kpoint(3,maxk), 
     .                  S(maxnh), xij(3,maxnh), occtol
      character(len=label_length+4) :: paste
      external          io_assign, io_close, memory
C *********************************************************************

      logical
     .  fixspin

      integer
     .  io, iu, iuo, naux, nuo

      real(dp)
     .  Dnew, qs(2), e1, e2, efs(2), Enew, eV, qk, qtot, 
     .  temp, wk, Entropy

C Dynamic arrays
      logical,           save :: getD = .false.
      logical,           save :: getPSI = .true.
      integer,        pointer :: muo(:)
      real(dp),       pointer :: aux(:)
      real(dp),       pointer :: ek(:,:,:)

      character(len=label_length+4) :: fname
      character(len=11)             :: fform

      parameter ( eV = 1.d0 / 13.60580d0 )
      save Dnew, Enew, e1, e2, qk, qtot, temp, wk
      data Dnew, Enew, e1, e2, qk, qtot, temp, wk /8*0.d0/

C Get local number of orbitals
#ifdef MPI
      call GetNodeOrbs(nuotot,Node,Nodes,nuo)
#else
      nuo = nuotot
#endif

C Start time counter 
      call timer( 'writewave', 1 )

C Check parameter maxk 
      if (nk .gt. maxk) then
        if (Node.eq.0) then
          write(6,'(/,a,/,a)')
     .       'writewave: WARNING: parameter maxk too small',
     .       'writewave: No wavefunction calculation performed'
        endif
        goto 999
      endif

C Check spin
      if (nspin .gt. 2) then
        if (Node.eq.0) then
          write(6,'(/,a,/,a)')
     .       'writewave: WARNING: wavefunctions printed is only '
          write(6,'(/,a,/,a)')
     .       '           implemented for nspin = 1 or 2'
        endif
        goto 999
      endif

C Allocate local arrays - only aux is relevant here
      naux = 10*nuotot
      call re_alloc( Haux, 1, 2*nuotot*nuo, 'Haux', 'densematrix' )
      call re_alloc( Saux, 1, 2*nuotot*nuo, 'Saux', 'densematrix' )
      call re_alloc( psi,  1, 2*nuotot*nuo, 'psi',  'densematrix' )
      nullify( ek, aux, muo )
      call re_alloc( ek,  1, nspin, 1, nuo, 1, nuotot,
     &              'ek', 'writewave')
      call re_alloc( aux, 1, naux, 'aux', 'writewave' )
      call re_alloc( muo, 1, nuotot, 'muo', 'writewave' )

C Check indxuo 
      do iuo = 1,nuotot
        muo(iuo) = 0
      enddo
      do io = 1,no
        iuo = indxuo(io)
        if (indxuo(io).le.0 .or. indxuo(io).gt.nuotot) then
          if (Node.eq.0) then
            write(6,*) 'writewave: invalid index: io, indxuo =',
     .        io, indxuo(io)
            call die('writewave: invalid indxuo')
          else
            call die()
          endif
        endif
        muo(iuo) = muo(iuo) + 1
      enddo
      do iuo = 1,nuotot
        if (muo(iuo) .ne. muo(1)) then
          if (Node.eq.0) then
          write(6,'(/,2a,3i6)') 'writewave: ERROR: inconsistent indxuo.'
     .     ,' iuo, muo(iuo), muo(1) =', iuo, muo(iuo), muo(1)
            call die('writewave: ERROR: inconsistent indxuo.')
          else
            call die()
          endif
        endif
      enddo


C Find output file name, and open for Node 0

      if (Node.eq.0) then

        fform = 'unformatted'

        fname = paste(slabel,'.WFS')

        call io_assign( iu )
        open(iu, file=fname, form=fform, status='unknown' )

        rewind (iu)

        if (fform .eq. 'formatted') then
          write(iu,*)
          write(iu,'(a22,2x,i6)') 'Nr of k-points = ',nk
          write(iu,'(a22,2x,i6)') 'Nr of Spins = ',nspin
          write(iu,'(a22,2x,i6)') 'Nr of basis orbs = ',nuotot
          write(iu,*)
        endif
        if (wwf) then
          write(6,*)
          write(6,'(a)') 'writewave: Wave Functions Coefficients'
          write(6,*)
          write(6,'(a26,2x,i6)') 'Number of k-points = ',nk
          write(6,'(a26,2x,i6)') 'Number of Spins = ',nspin
          write(6,'(a26,2x,i6)') 'Number of basis orbs = ',nuotot
          write(6,*)
        endif
        if (fform .eq. 'unformatted') then
          write(iu) nk
          write(iu) nspin
          write(iu) nuotot
        endif

        endfile (iu)
        backspace (iu)
        close (iu)
        call io_close(iu)

      endif

C Find the eigenvectors
c fixspin and qs are not used in diagk, since getD=.false. ...
      fixspin = .false.
      qs(1) = 0.0d0
      qs(2) = 0.0d0

C Call appropriate diagonalization routine

      if (gamma) then
         call diagg( nspin, nuo, maxuo, maxnh, maxnh, 
     .               maxo, numh, listhptr, listh, numh, listhptr, 
     .               listh, H, S, getD, getPSI,
     .               fixspin, qtot, qs, temp,
     .               e1, e2, ek, qk, Dnew, Enew, ef, efs, Entropy,
     .               Haux, Saux, psi, nuotot, occtol, 1, nuotot )
      else
         call diagk( nspin, nuo, no, maxspn, maxnh, maxnh, 
     .               maxo, numh, listhptr, listh, numh, listhptr, 
     .               listh, H, S, getD, getPSI,
     .               fixspin, qtot, qs, temp,
     .               e1, e2, xij, indxuo, nk, kpoint, wk,
     .               ek, qk, Dnew, Enew, ef, efs, Entropy,
     .               Haux, Saux, psi, Haux, Saux, aux, nuotot,
     .               occtol, 1 )
      endif

C Free local arrays
      call de_alloc( muo, 'muo', 'writewave' )
      call de_alloc( aux, 'aux', 'writewave' )
      call de_alloc( ek, 'ek', 'writewave' )

C This is the only exit point 
  999 continue
      call timer( 'writewave', 2 )

      end subroutine wwave


      subroutine writew(nuotot,nuo,ik,k,ispin,eo,psi,gamma)

      use precision
      use sys,          only : die
      use parallel,     only : Node, Nodes, BlockSize
      use parallelsubs, only : GlobalToLocalOrb, WhichNodeOrb
      use fdf
      use files,        only : slabel, label_length
      use atmfuncs,     only : symfio, cnfigfio, labelfis, nofis
      use siesta_geom,     only : isa
      use atomlist,     only : iaorb, iphorb

#ifdef MPI
      use mpi_siesta
#endif

      implicit          none

#ifdef MPI
      integer
     .  MPIerror
#endif

      integer nuotot, nuo, ispin, ik

      real(dp) eo(*), psi(*), k(3)

      logical gamma

C  Internal variables .............................................
      integer
     .  BNode, ntot, iie, iw, indwf, j, ind, iu

      logical, save :: frstime = .true.

      character(len=label_length+4), save :: fname
      character(len=label_length+4)       :: paste
      character(len=11),             save :: fform
      real(dp),                   pointer :: aux(:,:)

      external io_assign, io_close

C ...................

C Fix whether formatted or unformatted files wil be used
      fform = 'unformatted'

C Allocate auxiliary arrays
      call re_alloc( aux, 1, 2, 1, nuotot, 'aux', 'writewave' )

C Find file name, and open for Node 0
      if (frstime) then

        if (Node.eq.0) then
          fname = paste(slabel,'.WFS')
        endif

        frstime = .false.
      endif

      if (Node .eq. 0) then
        call io_assign( iu )
        open( iu, file=fname, form=fform, position='append',
     .        status='old' )
      endif

C Check that the total number of orbitals is correct

#ifdef MPI
      if (Nodes.gt.1) then
        call MPI_AllReduce(nuo,ntot,1,MPI_integer,MPI_sum,
     .    MPI_Comm_World,MPIerror)
      else
        ntot = nuo
      endif
#else
      ntot = nuo
#endif

      if (ntot .ne. nuotot) call die('Inconsistent number of orbitals')

C First print the index and value of k-point

      if (Node .eq. 0) then
        if (fform .eq. 'formatted') then
          write(iu,*)
          write(iu,'(a72)')    ' ***************************************
     .********************************'
          write(iu,'(a22,2x,i6,2x,3f10.6)') 'k-point = ',ik,
     .                                       k(1),k(2),k(3)
          write(iu,'(a22,2x,i6)') 'Spin component = ',ispin
          write(iu,'(a22,2x,i6)') 'Num. wavefunctions = ',nwflist(ik)
        endif
        if (wwf) then
          write(6,*)
          write(6,'(a72)')     ' ***************************************
     .********************************'
          write(6,'(a22,2x,i6,2x,3f10.6)') 'k-point = ',ik,
     .                                       k(1),k(2),k(3)
          write(6,'(a22,2x,i6)') 'Spin component = ',ispin
          write(6,'(a22,2x,i6)') 'Num. wavefunctions = ',nwflist(ik)
        endif
        if (fform .eq. 'unformatted') then
          write(iu) ik,k(1),k(2),k(3)
          write(iu) ispin
          write(iu) nwflist(ik)
        endif
      endif

C Loop over wavefunctions that should be printed

      do iw = 1,nwflist(ik)
        indwf = iwf(ik,iw)

C Determine which node handles this wavefunction
        call WhichNodeOrb(indwf,Nodes,BNode)

        if (Node.eq.BNode) then

C Determine the index of the orbital in the local node

          call GlobalToLocalOrb( indwf, BNode, Nodes, iie)

C Save wavefunction in aux array

C psi has different structure in diagk and diagg, so the indexing
C must be handled differently

          if (gamma) then
            do j = 1,ntot
              ind = j + (iie-1)*nuotot 
              aux(1,j) = psi(ind)
              aux(2,j) = 0.0d0
            enddo
          else
            do j = 1,ntot
              ind = 1+(j-1)*2+(iie-1)*2*nuotot
              aux(1,j) = psi(ind)
              ind = 2+(j-1)*2+(iie-1)*2*nuotot
              aux(2,j) = psi(ind)
            enddo
          endif

        endif

C Pass the wf to the other processors
#ifdef MPI
        if (BNode.ne.0) then
          call MPI_Bcast(aux(1,1),2*nuotot,MPI_double_precision,
     .          BNode,MPI_Comm_World,MPIerror)
        endif
#endif


C eigenvector is now stored in aux in all processors, and can be printed

        if (Node .eq. 0) then
          if (fform .eq. 'formatted') then
            write(iu,*)
            write(iu,'(a22,2x,i6)') 'Wavefunction = ', indwf
            write(iu,'(a22,2x,f10.6)') 'Energy (eV) = ', 
     .                                  eo(indwf)*13.60580d0
            write(iu,'(a72)')  ' ---------------------------------------
     .--------------------------------'
            write(iu,'(a72)')  '  Atom  Species Orb-global  Orb-in-atom
     . Orb-type     Re(psi)   Im(psi)'
            do j = 1,ntot
              write(iu,'(i6,5x,a3,1x,i10,8x,i3,7x,i1,a7,1x,2(f10.6))') 
     .                 iaorb(j),labelfis(isa(iaorb(j))),j,
     .                 iphorb(j), cnfigfio(isa(iaorb(j)),iphorb(j)),
     .                 symfio(isa(iaorb(j)),iphorb(j)),
     .                 aux(1,j), aux(2,j)
            enddo
            write(iu,'(a72)')  ' ---------------------------------------
     .--------------------------------'
          endif
          if (wwf) then
            write(6,*)
            write(6,'(a22,2x,i6)') 'Wavefunction = ', indwf
            write(6,'(a22,2x,f10.6)') 'Energy (eV) = ', 
     .                                  eo(indwf)*13.60580d0
            write(6,'(a72)')   ' ---------------------------------------
     .--------------------------------'
            write(6,'(a72)')  '  Atom  Species Orb-global  Orb-in-atom
     . Orb-type     Re(psi)   Im(psi)'
            do j = 1,ntot
              write(6,'(i6,5x,a3,1x,i10,8x,i3,7x,i1,a7,1x,2(f10.6))') 
     .                 iaorb(j),labelfis(isa(iaorb(j))),j,
     .                 iphorb(j), cnfigfio(isa(iaorb(j)),iphorb(j)),
     .                 symfio(isa(iaorb(j)),iphorb(j)),
     .                 aux(1,j), aux(2,j)
            enddo
            write(6,'(a72)')   ' ---------------------------------------
     .--------------------------------'
          endif
          if (fform .eq. 'unformatted') then
            write(iu) indwf
            write(iu) eo(indwf)*13.60580d0
            do j = 1,ntot
              write(iu) 
     .                 iaorb(j),labelfis(isa(iaorb(j))),j,
     .                 iphorb(j), cnfigfio(isa(iaorb(j)),iphorb(j)),
     .                 symfio(isa(iaorb(j)),iphorb(j)),
     .                 aux(1,j), aux(2,j)
            enddo
          endif
        endif

      enddo

C Close output file

      if (Node .eq. 0) then
        close (iu)
        call io_close(iu)
      endif

C Deallocate auxiliary arrays
      call de_alloc( aux, 'aux', 'writewave' )
        
      end subroutine writew

      end module writewave

