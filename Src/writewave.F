	      module writewave

      implicit none

      integer :: maxlin
      parameter (maxlin = 1000)

      character, save :: label(maxlin)*8
      integer, save   :: lastk(maxlin)

      integer, save   :: nwf
      integer, save, dimension(:), allocatable :: nwflist
      integer, save, dimension(:,:), allocatable :: iwf

      logical, save   :: wwf

      data nwf /1/
      data wwf /.false./
      

      contains

      subroutine initwave( maxk, norb, nk, kpoint, overflow)
C *********************************************************************
C Finds k-points for wavefunction printout
C Based on initband routine by J.Soler
C Written by P. Ordejon, June 2003
C **************************** INPUT **********************************
C integer maxk           : Last dimension of kpoint
C integer norb           : Number of orbitals
C *************************** OUTPUT **********************************
C integer nk             : Number k points to compute wavefunctions
C real*8  kpoint(3,maxk) : k point vectors
C logical overflow       : true if the internal memory was not sufficient
C                          in which case it should be called again
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
      use fdf
      use parsing
#ifdef MPI
      use mpi_siesta
      use parallel
#endif
      implicit          none
      integer           maxk, norb, nk
      double precision  kpoint(3,maxk)
      character         paste*30
      logical           overflow
      external          paste, redcel, memory
C *********************************************************************

      character 
     .  line*130, names*80, name1*10, name2*10, name3*10,
     .  scale*30

      logical
     .  WaveFuncPresent, Frstime, outlng

      integer
     .  i, ik, il, integs(4), iu, ix, iw, iw1, iw2, iw3,
     .  lastc, lc(0:3), mscell(3,3),
     .  ni, nkl, nn, nr, nv, Node, nmk

#ifdef MPI
      integer
     .  MPIerror
#endif

      double precision
     .  alat, caux(3,3), pi,
     .  rcell(3,3), reals(4), ucell(3,3), values(4)

      save
     .  Frstime

      data
     .  Frstime /.true./

C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
#else
      Node = 0
#endif

C Start time counter 
      call timer( 'writewave', 1 )

      overflow = .false.

      if (Frstime) then
         allocate(nwflist(maxk))
         call memory('A','I',maxk,'writewave')
         allocate(iwf(maxk,1))
         call memory('A','I',maxk,'writewave')
      endif


C Find k points if they are not given in argument 
      if (nk .le. 0) then

C       Find if there are k-points data
        if (Node.eq.0) then
          outlng = fdf_boolean('LongOutput', .false.)
          wwf = fdf_boolean('WriteWaveFunctions',outlng)
        endif

C       Find if there are k-points data
        if (Node.eq.0) then
          WaveFuncPresent = fdf_defined('WaveFuncKPoints')
        endif
#ifdef MPI
        call MPI_Bcast(WaveFuncPresent,1,MPI_logical,0,MPI_Comm_World,
     .    MPIerror)
#endif
        if ( WaveFuncPresent ) then

C         Find lattice constant
          if (Node.eq.0) then
            alat = fdf_physical( 'LatticeConstant', 0.d0, 'Bohr' )
          endif
#ifdef MPI
#ifdef NODAT
          call MPI_Bcast(alat,1,MPI_double_precision,0,
     .      MPI_Comm_World,MPIerror)
#else
          call MPI_Bcast(alat,1,DAT_double,0,
     .      MPI_Comm_World,MPIerror)
#endif
#endif
          if (alat .eq. 0.d0) then
            if (Node.eq.0) then
              write(6,'(a)') 
     .          'writewave: ERROR: Lattice constant required'
            endif
            goto 999
          endif

C         Find scale used in k point data
          if (Node.eq.0) then
            scale = fdf_string( 'WaveFuncKPointsScale', 'pi/a' )
          endif
#ifdef MPI
          call MPI_Bcast(scale,30,MPI_character,0,
     .      MPI_Comm_World,MPIerror)
#endif

          if (scale .eq. 'pi/a') then
            pi = 4.d0 * atan(1.d0)
          elseif (scale .eq. 'ReciprocalLatticeVectors') then
            if (Node.eq.0) then
              call redcel( alat, caux, ucell, mscell )
            endif
#ifdef MPI
#ifdef NODAT
            call MPI_Bcast(alat,1,MPI_double_precision,0,
     .        MPI_Comm_World,MPIerror)
            call MPI_Bcast(ucell(1,1),9,MPI_double_precision,0,
     .        MPI_Comm_World,MPIerror)
#else
            call MPI_Bcast(alat,1,DAT_double,0,
     .        MPI_Comm_World,MPIerror)
            call MPI_Bcast(ucell(1,1),9,DAT_double,0,
     .        MPI_Comm_World,MPIerror)
#endif
#endif
            call reclat( ucell, rcell, 1 )
          elseif (Node.eq.0) then
            write(6,'(a,/,2a,/,a)')
     .     'writewave: WARNING: Invalid value for WaveFuncKPointsScale',
     .     'writewave: Allowed values are pi/a and',
     .                ' ReciprocalLatticeVectors',
     .     'writewave: No band calculation performed'
          endif

C         Loop on data lines
          if (Node.eq.0) then
            nk = 0
            WaveFuncPresent = fdf_block('WaveFuncKPoints',iu)
            do il = 1,maxlin

C           Read and parse data line
              read(iu,'(a)',end=50) line
              lastc = index(line,'#') - 1
              if (lastc .le. 0) lastc = len(line)
              call parse( line(1:lastc), nn, lc, names, nv, values,
     .                  ni, integs, nr, reals )

C           Check if data are already finished
              if (nv .ge. 3) then

c             Check syntax
                if (nr.ne.0 .and. nr.ne.3) then
                  write(6,'(a,/,a)')
     .             'writewave: syntax ERROR in %block WaveFuncKPoints:',
     .              line(1:lastc)
                  goto 10
                endif


C             Add this point to total number of k points
                nk = nk + 1

C             If there is room to store k points
                if (nk .le. maxk) then

C               Find coordinates of k point
                  if (scale .eq. 'pi/a') then
                    kpoint(1,nk) = values(1) * pi / alat
                    kpoint(2,nk) = values(2) * pi / alat
                    kpoint(3,nk) = values(3) * pi / alat
                  elseif (scale .eq. 'ReciprocalLatticeVectors') then
                    do ix = 1,3
                      kpoint(ix,nk) = rcell(ix,1) * values(1) +
     .                                rcell(ix,2) * values(2) +
     .                                rcell(ix,3) * values(3)
                    enddo
                  endif
                else
                  overflow = .true.
                  if (.not. Frstime) then
                    stop 'writewave: invalid initialization'
                  endif
                endif


C             Do the next even if there is no room to store k points
C             so that the information on number of vectors to print
C             is available on next call

C             Find which eigenvectors should be printed

                if (nn .ge. 1) then

C               Check that line contains 'from', 'to' and maybe 'step'
                  if (nn .ne. 2 .and. nn .ne. 3) then
                    write(6,'(a,/,a)')
     .             'writewave: syntax ERROR in %block WaveFuncKPoints:',
     .              line(1:lastc)
                    goto 10
                  endif
                  name1 = names(lc(0)+1:lc(1))
                  name2 = names(lc(1)+1:lc(2))
                  if (nn .eq. 3) name3=names(lc(2)+1:lc(3))
                  if (name1 .ne. 'from' .or. name2 .ne. 'to') then
                    write(6,'(a,/,a)')
     .             'writewave: syntax ERROR in %block WaveFuncKPoints:',
     .              line(1:lastc)
                    goto 10
                  endif
                  if (nn .eq. 3) then
                    if (name3 .ne. 'step') then
                      write(6,'(a,/,a)')
     .             'writewave: syntax ERROR in %block WaveFuncKPoints:',
     .                line(1:lastc)
                      goto 10
                    endif
                  endif

                  iw1=integs(1)
                  iw2=integs(2)
                  if (iw1 .lt. 0) iw1 = norb + iw1 +1
                  if (iw2 .lt. 0) iw2 = norb + iw2 +1
                  if (nn .eq. 3) then
                    iw3 = abs(integs(3))
                  else
                    iw3 = 1
                  endif
                  ni = 0
                  do iw = min(iw1,iw2),max(iw1,iw2),iw3
                    ni = ni + 1
                    integs(ni) = iw
                  enddo

                endif

C               Store indexes of wave functions to printout

                if (ni .eq. 0) then
                  ni = norb
                  iw1=1
                  iw2=norb
                  do i=1,ni
                    integs(i)=i
                  enddo
                endif

                if (ni .le. nwf) then

                  if (nk .le. maxk) then
                    nwflist(nk) = ni
                    do i = 1,ni
                      iw = integs(i)
                      if (iw .lt. 0) iw = norb +iw +1
                      iwf(nk,i) = iw
                    enddo

                  else 

                    overflow = .true.

                    call memory('D','I',size(nwflist),'writewave')
                    deallocate(nwflist)
                    call memory('D','I',size(iwf),'writewave')
                    deallocate(iwf)

                    nmk=nk
                    allocate(nwflist(nmk))
                    call memory('A','I',nmk,'writewave')
                    allocate(iwf(nmk,nwf))
                    call memory('A','I',nwf*nmk,'writewave')

                  endif

                else 

                  overflow = .true.
                  if (.not. Frstime) then
                    stop 'writewave: invalid initialization'
                  endif

                  call memory('D','I',size(nwflist),'writewave')
                  deallocate(nwflist)
                  call memory('D','I',size(iwf),'writewave')
                  deallocate(iwf)

                  nwf = ni
                  nmk = max(nk,maxk)
                  allocate(nwflist(nmk))
                  call memory('A','I',nmk,'writewave')
                  allocate(iwf(nmk,nwf))
                  call memory('A','I',nwf*nmk,'writewave')

                endif
              else
C             No more lines to read => Exit do loop
                goto 50
              endif
   10         continue
            enddo
            write(6,'(a)') 'writewave: ERROR Parameter maxlin too small'
   50       continue
          endif
        else
C         No k-point data available => go to exit point
          goto 999
        endif
      endif

C Global broadcast of values
#ifdef MPI
      call MPI_Bcast(nk,1,MPI_integer,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(nwf,1,MPI_integer,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(nmk,1,MPI_integer,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(overflow,1,MPI_logical,0,MPI_Comm_World,MPIerror)

      if (.not. overflow) then
#ifdef NODAT
         call MPI_Bcast(kpoint(1,1),3*nk,MPI_double_precision,0,
     .     MPI_Comm_World,MPIerror)
#else
         call MPI_Bcast(kpoint(1,1),3*nk,DAT_double,0,
     .     MPI_Comm_World,MPIerror)
#endif
      endif

      if (overflow) then
         if (Node .ne. 0) then
           call memory('D','I',size(nwflist),'writewave')
           deallocate(nwflist)
           call memory('D','I',size(iwf),'writewave')
           deallocate(iwf)

           allocate(nwflist(nmk))
           call memory('A','I',nmk,'writewave')
           allocate(iwf(nmk,nwf))
           call memory('A','I',nwf*nmk,'writewave')
         endif
       endif

      if (.not. overflow) then
        call MPI_Bcast(nwflist(1),maxk,MPI_integer,0,
     .     MPI_Comm_World,MPIerror)
        call MPI_Bcast(iwf(1,1),maxk*nwf,MPI_integer,0,
     .     MPI_Comm_World,MPIerror)
      endif
 
#endif


C This is the only exit point 
  999 continue

      Frstime = .false.
      call timer( 'writewave', 2 )

      end subroutine initwave


      subroutine wwave( no, nspin, maxspn, maxo, maxuo, maxnh, maxk,
     .                  numh, listhptr, listh, H, S, ef, xij, indxuo,
     .                  nk, kpoint, nuotot, gamma)
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
      use fdf
#ifdef MPI
      use mpi_siesta
      use parallel
#endif
      implicit          none
      integer           maxk, maxnh, maxspn, maxo, maxuo, nk, no, nspin,
     .                  nuotot, indxuo(no), listh(maxnh), numh(*), 
     .                  listhptr(*)
      logical           gamma
      double precision  dot, ef, 
     .                  H(maxnh,nspin), kpoint(3,maxk), 
     .                  S(maxnh), xij(3,maxnh)
      character         paste*30
      external          dot, io_assign, io_close, paste, memory
C *********************************************************************

      logical
     .  getD, getPSI, fixspin

      integer
     .  ik, il, io, ispin, iu, iuo, 
     .  naux, nuo, Node, Nodes

#ifdef MPI
      integer
     .  MPIerror
#endif

      double precision
     .  Dnew, qs(2),
     .  e1, e2, efs(2), emax, emin, Enew, eV, qk, qtot, path, pi,
     .  temp, wk, Entropy

C     Dynamic arrays
      double precision, dimension(:), allocatable, save :: Haux, Saux
      double precision, dimension(:), allocatable, save :: aux
      double precision, dimension(:,:,:), allocatable :: psi
      double precision, dimension(:,:,:), allocatable, save :: ek
      integer         , dimension(:), allocatable, save :: muo

      character sname*30, fname*33, fform*11


      parameter ( eV = 1.d0 / 13.60580d0 )
      save getD, getPSI, Dnew, Enew, e1, e2, qk, qtot, temp, wk
      data getD /.false./
      data getPSI /.true./
      data Dnew, Enew, e1, e2, qk, qtot, temp, wk /8*0.d0/

C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
      call GetNodeOrbs(nuotot,Node,Nodes,nuo)
#else
      Node = 0
      Nodes = 1
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
      naux  = 2*nuotot*5
      allocate(Haux(2*nuo*nuotot))
      call memory('A','D',2*nuo*nuotot,'writewave')
      allocate(Saux(2*nuo*nuotot))
      call memory('A','D',2*nuo*nuotot,'writewave')
      allocate(psi(2,nuotot,nuo))
      call memory('A','D',2*nuo*nuotot,'writewave')
      allocate(ek(nspin,nuo,nuotot))
      call memory('A','D',nspin*nuo*nuotot,'writewave')
      allocate(aux(naux))
      call memory('A','D',naux,'writewave')
      allocate(muo(nuotot))
      call memory('A','I',nuotot,'writewave')

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
            stop 'writewave: invalid indxuo'
          else
            stop
          endif
        endif
        muo(iuo) = muo(iuo) + 1
      enddo
      do iuo = 1,nuotot
        if (muo(iuo) .ne. muo(1)) then
          if (Node.eq.0) then
          write(6,'(/,2a,3i6)') 'writewave: ERROR: inconsistent indxuo.'
     .     ,' iuo, muo(iuo), muo(1) =', iuo, muo(iuo), muo(1)
            stop 'writewave: ERROR: inconsistent indxuo.'
          else
            stop 
          endif
        endif
      enddo


C Find output file name, and open for Node 0

      if (Node.eq.0) then

        fform = 'unformatted'

        sname = fdf_string('SystemLabel','siesta')
        fname = paste(sname,'.WFS')

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
      qs(1)=0.0d0
      qs(2)=0.0d0
c ...
c

C Call appropriate diagonalization routine

      if (gamma) then
         call diagg( nspin, nuo, no, maxuo, maxnh, maxnh, 
     .               maxo, numh, listhptr, listh, numh, listhptr, 
     .               listh, H, S, getD, getPSI,
     .               fixspin, qtot, qs, temp,
     .               e1, e2, ek, qk, Dnew, Enew, ef, efs, Entropy,
     .               Haux, Saux, psi, aux,
     .               nuotot, Node, Nodes)
      else
         call diagk( nspin, nuo, no, maxspn, maxuo, maxnh, maxnh, 
     .               maxo, numh, listhptr, listh, numh, listhptr, 
     .               listh, H, S, getD, getPSI,
     .               fixspin, qtot, qs, temp,
     .               e1, e2, xij, indxuo, nk, kpoint, wk,
     .               ek, qk, Dnew, Enew, ef, efs, Entropy,
     .               Haux, Saux, psi, Haux, Saux, aux,
     .               nuotot, Node, Nodes)
      endif
 

C Free local arrays 
      call memory('D','D',size(Haux),'writewave')
      deallocate(Haux)
      call memory('D','D',size(Saux),'writewave')
      deallocate(Saux)
      call memory('D','D',size(psi),'writewave')
      deallocate(psi)
      call memory('D','I',size(ek),'writewave')
      deallocate(ek)
      call memory('D','D',size(aux),'writewave')
      deallocate(aux)
      call memory('D','I',size(muo),'writewave')
      deallocate(muo)

C This is the only exit point 
  999 continue
      call timer( 'writewave', 2 )

      end subroutine wwave


      subroutine writew(nuotot,nuo,ik,k,ispin,eo,psi,gamma)

      use precision
      use sys
      use parallel
      use fdf
      use atmfuncs, only: symfio, cnfigfio, labelfis, nofis
      use atomlist, only: isa, iaorb, iphorb

#ifdef MPI
      use mpi_siesta
#endif

      implicit          none

#ifdef MPI
      integer
     .  MPIerror
#endif

      integer nuotot, nuo, ispin, ik

      double precision eo(*), psi(*), k(3)
C      double precision eo(*), psi(2,nuotot,nuo), k(3)

      logical gamma

C  Internal variables .............................................
      integer
     .  BNode, BTest, Node, Nodes, ntot, iie, ie, iw, indwf, j, ind,
     .  iu

      logical frstime, writeie

      character sname*30, fname*33, paste*33, fform*11

      double precision, dimension(:,:), allocatable :: aux

      save frstime, sname, fform, fname
      data frstime /.true./

      external paste, io_assign, io_close

C ...................

c Fix whether formatted or unformatted files wil be used
      fform = 'unformatted'

C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
#else
      Node = 0
      Nodes = 1
#endif

C Allocate auxiliary arrays

      allocate(aux(2,nuotot))
      call memory('A','D',2*nuotot,'writewave')

C Find file name, and open for Node 0
      if (frstime) then

        if (Node.eq.0) then
          sname = fdf_string('SystemLabel','siesta')
          fname = paste(sname,'.WFS')
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

c Determine which node handles this wavefunction

        BNode = indwf/BlockSize
        if (BNode*BlockSize.eq.indwf) then
          BNode = BNode + 1
11        continue
          if (BNode .gt. Nodes-1) then
            BNode =BNode - Nodes
            goto 11
          endif
        endif

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
#ifdef NODAT
        call MPI_Bcast(aux(1,1),2*nuotot,MPI_double_precision,
     .        BNode,MPI_Comm_World,MPIerror)
#else
        call MPI_Bcast(aux(1,1),2*nuotot,DAT_double,
     .        BNode,MPI_Comm_World,MPIerror)
#endif
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
     . Orb-type      Re(psi)   Im(psi)'
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
     . Orb-type      Re(psi)   Im(psi)'
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

      call memory('D','D',size(aux),'writewave')
      deallocate(aux)

        
      end subroutine writew

      end module writewave
