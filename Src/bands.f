C $Id: bands.f,v 1.10 1999/05/25 13:56:03 ordejon Exp $

      subroutine bands( no, nspin, maxspn, maxuo, maxno, maxk,
     .                  numh, listh, H, S, ef, xij, indxuo,
     .                  writeb, nk, kpoint, ek )
C *********************************************************************
C Finds band energies at selected k-points.
C Written by J.Soler, August 1997 and August 1998.
C **************************** INPUT **********************************
C integer no                  : Number of basis orbitals
C integer nspin               : Number of spin components
C integer maxspn              : Second dimension of ek
C integer maxno               : Maximum number of orbitals interacting  
C                               with any orbital
C integer maxk                : Last dimension of kpoint and ek
C integer numh(no)            : Number of nonzero elements of each row 
C                               of hamiltonian matrix
C integer listh(maxno,no)     : Nonzero hamiltonian-matrix element  
C                               column indexes for each matrix row
C real*8  H(maxno,maxuo,nspin): Hamiltonian in sparse form
C real*8  S(maxno,maxuo)      : Overlap in sparse form
C real*8  ef                  : Fermi energy
C real*8  xij(3,maxno,maxuo)  : Vectors between orbital centers (sparse)
C                               (not used if only gamma point)
C integer maxuo               : First dimension of ek
C integer indxuo(no)          : Index of equivalent orbital in unit cell
C                               Unit cell orbitals must be the first in
C                               orbital lists, i.e. indxuo.le.nuo, with
C                               nuo the number of orbitals in unit cell
C real*8                      : Fermi energy
C logical writeb              : This routine must write bands?
C *************************** INPUT OR OUTPUT *************************
C integer nk             : Number of band k points
C real*8  kpoint(3,maxk) : k point vectors
C *************************** OUTPUT **********************************
C real*8  ek(maxuo,maxspn,maxk) : Eigenvalues
C *************************** UNITS ***********************************
C Lengths in atomic units (Bohr).
C k vectors in reciprocal atomic units.
C Energies in Rydbergs.
C ***************** BEHAVIOUR *****************************************
C - When writeb=true, bands are saved in file sys_name.bands, where
C   sys_name is the value of fdf label SystemLabel, or 'siesta'
C   by default.
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
      implicit          none
      integer           maxk, maxno, maxspn, maxuo, nk, no, nspin
      integer           indxuo(no), listh(maxno,no), numh(no)
      logical           writeb
      double precision  dot, ef, ek(maxuo,maxspn,maxk),
     .                  H(maxno,maxuo,nspin), kpoint(3,maxk), 
     .                  S(maxno,maxuo), xij(3,maxno,maxuo)
      character         paste*30
      external          cdiag, dot, io_assign, io_close,
     .                  parse, paste, prmem, rdiag, redcel
      include          'fdf/fdfdefs.h'
C *********************************************************************

C Internal variables 
      include 'diagon.h'

      integer maxlin
      parameter (maxlin = 1000)

      character 
     .  fname*30, label(maxlin)*8, line*130, names*80,
     .  scale*30, sname*25, string*10

      logical
     .  fixspin, found, frstme, getD

      integer
     .  ik, il, integs(4), io, ispin, iu, iuo, ix, 
     .  lastc, lastk(maxlin), lc(0:3), 
     .  mscell(3,3), muo(maxorb),
     .  naux, ni, nkl, nlines, nn, nr, nuo, nv

      double precision
     .  alat, Dnew, e1, e2, efs(2), emax, emin, Enew, eV,
     .  qk, qs(2), qtot, path, pi, rcell(3,3), reals(4),
     .  scell(3,3), temp, ucell(3,3), values(4), wk

C     Common auxiliary arrays shared with routine diagon
      double precision
     .  Haux(maxhs), Saux(maxhs), psi(maxpsi), aux(maxaux)
      common /diacom/ Haux, Saux, psi, aux

      parameter ( eV = 1.d0 / 13.60580d0 )
      save frstme, getD, Dnew, Enew, e1, e2, qk, qtot, temp, wk
      data frstme /.true./
      data getD /.false./
      data Dnew, Enew, e1, e2, qk, qtot, temp, wk /8*0.d0/

C Start time counter 
      call timer( 'bands', 1 )
      pi = 4.d0 * atan(1.d0)

C Find k points if they are not given in argument 
      if (nk .le. 0) then

C       Find if there are band-lines data
        if ( fdf_block('BandLines',iu) ) then

C         Find lattice constant
          alat = fdf_physical( 'LatticeConstant', 0.d0, 'Bohr' )
          if (alat .eq. 0.d0) then
            write(6,'(a)') 'bands: ERROR: Lattice constant required'
            goto 999
          endif

C         Find scale used in k point data
          scale = fdf_string( 'BandLinesScale', 'pi/a' )
          if (scale .eq. 'pi/a') then
            pi = 4.d0 * atan(1.d0)
          elseif (scale .eq. 'ReciprocalLatticeVectors') then
            call redcel( alat, ucell, scell, mscell )
            call reclat( ucell, rcell, 1 )
          else
            write(6,'(a,/,2a,/,a)')
     .        'bands: WARNING: Invalid value for BandLinesScale',
     .        'bands: Allowed values are pi/a and',
     .              ' ReciprocalLatticeVectors',
     .        'bands: No band calculation performed'
          endif

C         Loop on data lines
          nk = 0
          do il = 1,maxlin

C           Read and parse data line
            if (il .eq. 1) found = fdf_block('BandLines',iu)
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
     .                              rcell(ix,2) * values(3) +
     .                              rcell(ix,3) * values(4)
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
   50     continue
          nlines = il - 1
        else
C         No k-point data available => go to exit point
          goto 999
        endif
      endif

C Check parameter maxk 
      if (nk .gt. maxk) then
*       write(6,'(/,a,/,a)')
*    .    'bands: WARNING: parameter maxk too small',
*    .    'bands: No bands calculation performed'
        goto 999
      endif

C Find number of orbitals per unit cell and check argument sizes 
      nuo = 0
      do io = 1,no
        nuo = max( nuo, indxuo(io) )
      enddo
      call chkdim( 'bands', 'maxuo', maxuo, nuo, 1 )

C Check internal dimensions 
      naux  = nuo*5
      call chkdim( 'bands', 'maxaux', maxaux, naux, 1 )
      call chkdim( 'bands', 'maxorb', maxorb, nuo,  1 )

C Check indxuo 
      do iuo = 1,nuo
        muo(iuo) = 0
      enddo
      do io = 1,no
        iuo = indxuo(io)
        if (indxuo(io).le.0 .or. indxuo(io).gt.no) then
          write(6,*) 'bands: invalid index: io, indxuo =',
     .      io, indxuo(io)
          stop 'bands: invalid indxuo'
        endif
        muo(iuo) = muo(iuo) + 1
      enddo
      do iuo = 1,nuo
        if (muo(iuo) .ne. muo(1)) then
          write(6,'(/,2a,3i6)') 'bands: ERROR: inconsistent indxuo.',
     .     ' iuo, muo(iuo), muo(1) =', iuo, muo(iuo), muo(1)
          stop 'bands: ERROR: inconsistent indxuo.'
        endif
      enddo

C Find the band energies 
      if (nspin.le.2) then
c fixspin and qs are not used in diagk, since getD=.false. ...
        fixspin = .false.
        qs(1)=0.0d0
        qs(2)=0.0d0
c ...
        call diagk( nspin, nuo, no, maxspn, maxuo, maxno, maxno, 
     .              numh, listh, numh, listh, H, S,
     .              getD, fixspin, qtot, qs, temp, e1, e2, 
     .              xij, indxuo, nk, kpoint, wk,
     .              ek, qk, Dnew, Enew, ef, efs,
     .              Haux, Saux, psi, Haux, Saux, aux )
      elseif (nspin.eq.4) then
        call diag2k( nuo, no, maxuo, maxno, maxno, 
     .               numh, listh, numh, listh, H, S,
     .               getD, qtot, temp, e1, e2,
     .               xij, indxuo, nk, kpoint, wk,
     .               ek, qk, Dnew, Enew, ef,
     .               Haux, Saux, psi, Haux, Saux, aux )
      else
        stop 'bands: ERROR: incorrect value of nspin'
      endif

C Write bands 
      if (writeb) then

C       Find name of output file and open it
        sname = fdf_string('SystemLabel','siesta')
        fname = paste(sname,'.bands')
        call io_assign(iu)
*       write(6,*) 'bands: iu,fname=', iu, fname
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
          write(iu,'(f10.6,10f12.5,/,(10x,10f12.5))')
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

C This is the only exit point 
  999 continue
      frstme = .false.
      call timer( 'bands', 2 )
      end


