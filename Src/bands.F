      module band

      implicit none

      integer :: maxlin
      parameter (maxlin = 1000)

      character, save :: label(maxlin)*8
      integer, save   :: lastk(maxlin), nlines

      contains

      subroutine initbands( maxk, nk, kpoint )
C *********************************************************************
C Finds band k-points
C Based on initialisation part of subroutine bands from original code.
C Written by J.Soler, August 1997 and August 1998.
C **************************** INPUT **********************************
C integer maxk           : Last dimension of kpoint
C *************************** OUTPUT **********************************
C integer nk             : Number of band k points
C real*8  kpoint(3,maxk) : k point vectors
C *************************** UNITS ***********************************
C Lengths in atomic units (Bohr).
C k vectors in reciprocal atomic units.
C ***************** BEHAVIOUR *****************************************
C - If nk=0 on input, k-points are read from labels BandLines and 
C   BandLinesScale of the input fdf data file. If these labels are 
C   not present, it returns with nk=0.
C - Allowed values for BandLinesScale are ReciprocalLatticeVectors and
C   pi/a (default). If another value is given, it returns with nk=0
C   after printing a warning.
C - If nk>maxk, k points and bands are not calculated and no warning
C   is printed before return
C ***************** USAGE *********************************************
C Example of fdf band lines specification for an FCC lattice.
C Last column is an optional LaTex label (for plot)
C     BandLinesScale  pi/a
C     %block BandLines                  # These are comments
C      1  0.000  0.000  0.000  \Gamma   # Begin at Gamma
C     25  2.000  0.000  0.000     X     # 25 points from Gamma to X
C     10  2.000  1.000  0.000     W     # 10 points from X to W
C     15  1.000  1.000  1.000     L     # 15 points from W to L
C     20  0.000  0.000  0.000  \Gamma   # 20 points from L to Gamma
C     25  1.500  1.500  1.500     K     # 25 points from Gamma to K
C     %endblock BandLines
C
C Example for BCC:
C     BandLinesScale  pi/a
C     %block BandLines
C      1  0.000  0.000  0.000  \Gamma
C     20  2.000  0.000  0.000     H
C     15  1.000  1.000  0.000     N
C     15  0.000  0.000  0.000  \Gamma
C     20  1.000  1.000  1.000     P
C     10  1.000  1.000  0.000     N
C     10  1.000  1.000  1.000     P
C     20  2.000  2.000  2.000     H
C     %endblock BandLines
C
C Example for HCP (an angle of 120 deg is assumed between reciprocal
C lattice vectors, what implies an angle of 60 deg between the first 
C two vectors of cell argument):
C     BandLinesScale  ReciprocalLatticeVectors
C     %block BandLines
C      1  0.000000000  0.000000000  0.000000000  \Gamma
C     20  0.666666667  0.333333333  0.000000000     K 
C     10  0.500000000  0.000000000  0.000000000     M
C     20  0.000000000  0.000000000  0.000000000  \Gamma
C     15  0.000000000  0.000000000  0.500000000     A
C     20  0.666666667  0.333333333  0.500000000     H
C     10  0.500000000  0.000000000  0.500000000     L
C     %endblock BandLines
C
C If only given points (not lines) are desired, simply specify 1 as 
C the number of points along the line.
C *********************************************************************
C
C  Modules
C
      use precision
      use fdf
#ifdef MPI
      use mpi
      use parallel
#endif
      implicit          none
      integer           maxk, nk
      double precision  kpoint(3,maxk)
      character         paste*30
      external          parse, paste, redcel, memory
C *********************************************************************

      character 
     .  fname*30, line*130, names*80,
     .  scale*30, sname*25, string*10

      logical
     .  BandLinesPresent

      integer
     .  ik, il, integs(4), iu, ix, 
     .  lastc, lc(0:3), 
     .  naux, ni, nkl, nn, nr, nv,
     .  Node

#ifdef MPI
      integer
     .  MPIerror
#endif

      double precision
     .  alat, caux(3,3), pi,
     .  rcell(3,3), reals(4), ucell(3,3), values(4)

C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
#else
      Node = 0
#endif

C Start time counter 
      call timer( 'bands', 1 )

C Initialise the number of band lines
      nlines = 0

C Find k points if they are not given in argument 
      if (nk .le. 0) then

C       Find if there are band-lines data
        if (Node.eq.0) then
          BandLinesPresent = fdf_defined('BandLines')
        endif
#ifdef MPI
        call MPI_Bcast(BandLinesPresent,1,MPI_logical,0,MPI_Comm_World,
     .    MPIerror)
#endif
        if ( BandLinesPresent ) then

C         Find lattice constant
          if (Node.eq.0) then
            alat = fdf_physical( 'LatticeConstant', 0.d0, 'Bohr' )
          endif
#ifdef MPI
          call MPI_Bcast(alat,1,DAT_double,0,
     .      MPI_Comm_World,MPIerror)
#endif
          if (alat .eq. 0.d0) then
            if (Node.eq.0) then
              write(6,'(a)') 'bands: ERROR: Lattice constant required'
            endif
            goto 999
          endif

C         Find scale used in k point data
          if (Node.eq.0) then
            scale = fdf_string( 'BandLinesScale', 'pi/a' )
          endif
#ifdef MPI
          call MPI_Bcast(scale,1,DAT_double,0,
     .      MPI_Comm_World,MPIerror)
#endif
          if (scale .eq. 'pi/a') then
            pi = 4.d0 * atan(1.d0)
          elseif (scale .eq. 'ReciprocalLatticeVectors') then
            if (Node.eq.0) then
              call redcel( alat, caux, ucell )
            endif
#ifdef MPI
            call MPI_Bcast(alat,1,DAT_double,0,
     .        MPI_Comm_World,MPIerror)
            call MPI_Bcast(ucell(1,1),9,DAT_double,0,
     .        MPI_Comm_World,MPIerror)
#endif
            call reclat( ucell, rcell, 1 )
          elseif (Node.eq.0) then
            write(6,'(a,/,2a,/,a)')
     .        'bands: WARNING: Invalid value for BandLinesScale',
     .        'bands: Allowed values are pi/a and',
     .              ' ReciprocalLatticeVectors',
     .        'bands: No band calculation performed'
          endif

C         Loop on data lines
          if (Node.eq.0) then
            nk = 0
            BandLinesPresent = fdf_block('BandLines',iu)
            do il = 1,maxlin

C           Read and parse data line
              read(iu,'(a)',end=50) line
              lastc = index(line,'#') - 1
              if (lastc .le. 0) lastc = len(line)
              call parse( line(1:lastc), nn, lc, names, nv, values,
     .                  ni, integs, nr, reals )

C           Check if data are already finished
              if (nv .ge. 3) then

C             Add to total number of k points
                nkl = integs(1)
                nk = nk + nkl

C             If there is room to store k points
                if (nk .le. maxk) then

C               Find last point in line
                  if (scale .eq. 'pi/a') then
                    kpoint(1,nk) = values(2) * pi / alat
                    kpoint(2,nk) = values(3) * pi / alat
                    kpoint(3,nk) = values(4) * pi / alat
                  elseif (scale .eq. 'ReciprocalLatticeVectors') then
                    do ix = 1,3
                      kpoint(ix,nk) = rcell(ix,1) * values(2) +
     .                                rcell(ix,2) * values(3) +
     .                                rcell(ix,3) * values(4)
                    enddo
                  endif

C               Find points along the line
                  do ik = 1,nkl-1
                    do ix = 1,3
                      kpoint(ix,nk-nkl+ik) =
     .                  kpoint(ix,nk-nkl) * dble(nkl-ik) / dble(nkl) + 
     .                  kpoint(ix,nk)     * dble(ik)     / dble(nkl)
                    enddo
                  enddo

C               Find point label
                  if (nn .gt. 0) then
                    label(il) = names(1:lc(1))
                  else
                    label(il) = ' '
                  endif
                  lastk(il) = nk

                endif
              else
C             No more lines to read => Exit do loop
                goto 50
              endif
            enddo
            write(6,'(a)') 'bands: ERROR. Parameter maxlin too small'
   50       continue
            nlines = il - 1
          endif
        else
C         No k-point data available => go to exit point
          goto 999
        endif
      endif

C Global broadcast of values
#ifdef MPI
      call MPI_Bcast(nk,1,MPI_integer,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(kpoint(1,1),3*nk,DAT_double,0,
     .  MPI_Comm_World,MPIerror)
#endif

C This is the only exit point 
  999 continue
      call timer( 'bands', 2 )

      end subroutine initbands


      subroutine bands( no, nspin, maxspn, maxo, maxuo, maxnh, maxk,
     .                  numh, listhptr, listh, H, S, ef, xij, indxuo,
     .                  writeb, nk, kpoint, ek, nuotot )
C *********************************************************************
C Finds band energies at selected k-points.
C Written by J.Soler, August 1997 and August 1998.
C Initialisation moved into a separate routine, JDG Jan 2000.
C **************************** INPUT **********************************
C integer no                  : Number of basis orbitals
C integer nspin               : Number of spin components
C integer maxspn              : Second dimension of ek
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
C integer maxo                : First dimension of ek
C integer maxuo               : Second dimension of H and S
C integer indxuo(no)          : Index of equivalent orbital in unit cell
C                               Unit cell orbitals must be the first in
C                               orbital lists, i.e. indxuo.le.nuo, with
C                               nuo the number of orbitals in unit cell
C real*8  ef                  : Fermi energy
C logical writeb              : This routine must write bands?
C integer nuotot              : Total number of orbitals in unit cell
C integer nk                  : Number of band k points
C real*8  kpoint(3,maxk)      : k point vectors
C *************************** OUTPUT **********************************
C real*8  ek(maxo,maxspn,maxk) : Eigenvalues
C *************************** UNITS ***********************************
C Lengths in atomic units (Bohr).
C k vectors in reciprocal atomic units.
C Energies in Rydbergs.
C ***************** BEHAVIOUR *****************************************
C - When writeb=true, bands are saved in file sys_name.bands, where
C   sys_name is the value of fdf label SystemLabel, or 'siesta'
C   by default.
C *********************************************************************
C
C  Modules
C
      use precision
      use fdf
#ifdef MPI
      use mpi
      use parallel
#endif
      implicit          none
      integer           maxk, maxnh, maxspn, maxo, maxuo, nk, no, nspin,
     .                  nuotot, indxuo(no), listh(maxnh), numh(*), 
     .                  listhptr(*)
      logical           writeb
      double precision  dot, ef, ek(maxo,maxspn,maxk),
     .                  H(maxnh,nspin), kpoint(3,maxk), 
     .                  S(maxnh), xij(3,maxnh)
      character         paste*30
      external          dot, io_assign, io_close, paste, memory
C *********************************************************************

      character 
     .  fname*30, sname*25, string*10

      logical
     .  getD, fixspin

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
     .  temp, wk

C     Dynamic arrays
      double precision, dimension(:), allocatable, save :: Haux, Saux
      double precision, dimension(:), allocatable, save :: psi, aux
      integer         , dimension(:), allocatable, save :: muo

      parameter ( eV = 1.d0 / 13.60580d0 )
      save getD, Dnew, Enew, e1, e2, qk, qtot, temp, wk
      data getD /.false./
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
      call timer( 'bands', 1 )

C Check parameter maxk 
      if (nk .gt. maxk) then
        if (Node.eq.0) then
          write(6,'(/,a,/,a)')
     .       'bands: WARNING: parameter maxk too small',
     .       'bands: No bands calculation performed'
        endif
        goto 999
      endif

C Allocate local arrays - only aux is relevant here
      naux  = nuotot*5
      allocate(Haux(2*nuo*nuotot))
      call memory('A','D',2*nuo*nuotot,'bands')
      allocate(Saux(2*nuo*nuotot))
      call memory('A','D',2*nuo*nuotot,'bands')
      allocate(psi(2*nuo*nuotot))
      call memory('A','D',2*nuo*nuotot,'bands')
      allocate(aux(naux))
      call memory('A','D',naux,'bands')
      allocate(muo(nuotot))
      call memory('A','I',nuotot,'bands')

C Check indxuo 
      do iuo = 1,nuotot
        muo(iuo) = 0
      enddo
      do io = 1,no
        iuo = indxuo(io)
        if (indxuo(io).le.0 .or. indxuo(io).gt.nuotot) then
          if (Node.eq.0) then
            write(6,*) 'bands: invalid index: io, indxuo =',
     .        io, indxuo(io)
            stop 'bands: invalid indxuo'
          else
            stop
          endif
        endif
        muo(iuo) = muo(iuo) + 1
      enddo
      do iuo = 1,nuotot
        if (muo(iuo) .ne. muo(1)) then
          if (Node.eq.0) then
          write(6,'(/,2a,3i6)') 'bands: ERROR: inconsistent indxuo.',
     .     ' iuo, muo(iuo), muo(1) =', iuo, muo(iuo), muo(1)
            stop 'bands: ERROR: inconsistent indxuo.'
          else
            stop 
          endif
        endif
      enddo

C Find the band energies 
      if (nspin.le.2) then
c fixspin and qs are not used in diagk, since getD=.false. ...
        fixspin = .false.
        qs(1)=0.0d0
        qs(2)=0.0d0
c ...
        call diagk( nspin, nuo, no, maxspn, maxuo, maxnh, maxnh, 
     .              maxo, numh, listhptr, listh, numh, listhptr, 
     .              listh, H, S, getD, fixspin, qtot, qs, temp, 
     .              e1, e2, xij, indxuo, nk, kpoint, wk,
     .              ek, qk, Dnew, Enew, ef, efs,
     .              Haux, Saux, psi, Haux, Saux, aux,
     .              nuotot, Node, Nodes)
      elseif (nspin.eq.4) then
        call diag2k(nuo, no, maxuo, maxnh, maxnh, maxo,
     .              numh, listhptr, listh, numh, listhptr, 
     .              listh, H, S, getD, qtot, temp, e1, e2,
     .              xij, indxuo, nk, kpoint, wk,
     .              ek, qk, Dnew, Enew, ef,
     .              Haux, Saux, psi, Haux, Saux, aux,
     .              nuotot, Node, Nodes)
      else
        if (Node.eq.0) then
          stop 'bands: ERROR: incorrect value of nspin'
        else
          stop
        endif
      endif

C Write bands 
      if (writeb.and.Node.eq.0) then

C       Find name of output file and open it
        sname = fdf_string('SystemLabel','siesta')
        fname = paste(sname,'.bands')
        call io_assign(iu)
        open( iu, file=fname, status='unknown')

C       Write Fermi energy
        write(iu,*) ef/eV

C       Find and write the ranges of k and ek
        path = 0.d0
        emax = ek(1,1,1)
        emin = ek(1,1,1)
        do ik = 1,nk
          if (ik .gt. 1)
     .      path = path + sqrt( (kpoint(1,ik)-kpoint(1,ik-1))**2 +
     .                          (kpoint(2,ik)-kpoint(2,ik-1))**2 +
     .                          (kpoint(3,ik)-kpoint(3,ik-1))**2 )
          do ispin = 1,min(nspin,2)
            do io = 1, nuo
              emax = max( emax, ek(io,ispin,ik) )
              emin = min( emin, ek(io,ispin,ik) )
            enddo
          enddo
        enddo
        write(iu,*) 0.d0, path
        write(iu,*) emin/eV, emax/eV

C       Write eigenvalues
        write(iu,*) nuo, min(nspin,2), nk
        path = 0.d0
        do ik = 1,nk
          if (ik .gt. 1)
     .      path = path + sqrt( (kpoint(1,ik)-kpoint(1,ik-1))**2 +
     .                          (kpoint(2,ik)-kpoint(2,ik-1))**2 +
     .                          (kpoint(3,ik)-kpoint(3,ik-1))**2 )
          write(iu,'(f10.6,10f12.4,/,(10x,10f12.4))')
     .      path, ((ek(io,ispin,ik)/eV,io=1,nuo),ispin=1,min(nspin,2))
        enddo

C       Write abscisas of line ends and their labels
        write(iu,*) nlines
        il = 1
        path = 0.d0
        do ik = 1,nk
          if (ik .gt. 1)
     .      path = path + sqrt( (kpoint(1,ik)-kpoint(1,ik-1))**2 +
     .                          (kpoint(2,ik)-kpoint(2,ik-1))**2 +
     .                          (kpoint(3,ik)-kpoint(3,ik-1))**2 )
          if (ik .eq. lastk(il)) then
C           Put label between quotes
            if (label(il) .eq. ' ') then
              string = ''' '''
            else
              string = paste( ''''//label(il),'''' )
            endif
            write(iu,'(f12.6,3x,a)') path, string
            il = il + 1
          endif
        enddo

C       Close output file
        call io_close(iu)
      endif

C Free local arrays 
      call memory('D','D',size(Haux),'bands')
      deallocate(Haux)
      call memory('D','D',size(Saux),'bands')
      deallocate(Saux)
      call memory('D','D',size(psi),'bands')
      deallocate(psi)
      call memory('D','D',size(aux),'bands')
      deallocate(aux)
      call memory('D','I',size(muo),'bands')
      deallocate(muo)

C This is the only exit point 
  999 continue
      call timer( 'bands', 2 )

      end subroutine bands

      end module band
