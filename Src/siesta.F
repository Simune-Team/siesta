      Program SIESTA

C ***
C SIESTA Density Functional LCAO program package.
C Copyright by E.Artacho, P.Ordejon, D.Sanchez-Portal and 
C   J.M.Soler, 1996, 1997.
C The use of this program is allowed for non-for-profit research only.
C Copy or disemination of all or part of this package is not
C permitted without prior and explicit authorization by the authors.
C Send comments/suggestions/bug-reports to:
C   Jose M. Soler: jose.soler@uam.es
C   Pablo Ordejon: ordejon@icmab.es
C ***
C This program uses the flexible data format (FDF) package,
C copyrighted by A.Garcia and J.M.Soler, 1996, 1997
C ***
C This program uses exchange-correlation routines copyrighted
C by C.Balbas, J.L.Martins and J.M.Soler, 1996, 1997.
C ***
C Routine idiag copyrighted by J.Moreno and J.M.Soler, 1990.
C ***
C Routine pulayx copyrighted by In-Ho Lee and P.Ordejon, 1997.
C ***
C This program uses routines from 'Numerical Recipes', copyrighted
C by W.H.Press, A.A.Teukolsky, W.T.Veterlig and B.P.Flannery.
C Some of these routines have been modified or their name changed.
C ***
C The atomic part uses old routines from unknown authors. If you wrote
C any of these routines, please let us know, to include the apropriate
C acknowledgements.
C ***
C This module was written by P.Ordejon and J.M.Soler, 1996, 1997, 1998.
C ***
C Parallelisation changes introduced JDG, 1998, 1999.
C ***
C General architectural issues: Alberto Garcia, 1997, 1998, 1999, 2000
C ***
C
C  Modules
C
      use precision
      use sys
      use interfaces   ! 
      use parallel
      use atmfuncs, only: rcut, uion
      use atomlist
      use band
      use fdf
      use parse
      use ionew, only: io_setup, IOnode
#ifdef MPI
      use mpi
#endif
      use alloc
      use phonon

      implicit none

      integer
     .  fincoor,
     .  i, ia, ia1, ia2, iadispl, ianneal,
     .  idyn, ifa, ifinal, ihmat, iiscf, ik, ikb, in, iio,
     .  ind, inicoor, io, ioa, ioptlwf, iord,
     .  iquench, is, isel, iscf, 
     .  isolve, ispin, istp, istpsave, istart, istep, istr,
     .  iu, iunit,  iv, ix, ixdispl, 
     .  j, ja, jamin, jo, jx, kscell(3,3),
     .  lastc, lc(0:1), level, maxbk, maxk, maxnh, maxnhold,
     .  maxna, maxsav, mscell(3,3), mullipop, ncells,
     .  nkick, nnamax, nauxpul, nbcell, nbk, ncgmax, nh,
     .  ni, nkpnt, nmove, nn, nnia, nr,
     .  ns, nsc(3), nscf, nspin, ntm(3), no_l, nv, 
     .  nxij, pmax, nkpol, nhist

      integer
     .  Node, Nodes
#ifdef MPI
      integer
     .  MPIerror, ntmp, procYval, procYdefault
#endif

      integer, dimension(:), allocatable, save ::
     .  ibuffer1, 
     .  jna, listh, listhold, listhptr, listhptrold,
     .  numh, numhold

      real(dp)
     .  amu, Ang, bcell(3,3), beta, bulkm, 
     .  charnet, cfmax, cftem, const, cstress(3,3),
     .  Debye, dipol(3), dDmax, dDtol, DEna,
     .  dt, DUext, DUscf, Dxc, dx, dxmax, e1, e2, sigma,
     .  Ecorrec, ef, Eharrs, Eions, Ekin, Ekinion, Emad, Ena, 
     .  Enaatm, Enascf, Enl, Entrop, 
     .  eta, etol, Etot, eV, Exc, E0,
     .  factor,
     .  fmax, fmean, FreeE, fres, ftem, ftol, ftot(3),
     .  g2cut, g2max, getot,
     .  kBar, kcutof, kdispl(3), Kelvin,
     .  kn, kpr, 
     .  mn, mpr, Pint, Pmol, Psol
      real(dp)
     .  qaux, qspin(4), qsol,
     .  rc, rcoor, rcoorcp, rijmin, 
     .  rmax, rmaxh, rmin, r2min,
     .  scell(3,3), stot, stress(3,3), strtol, svec(3),
     .  taurelax, temp, tempinit, tempion, tiny, tp, ts, tt,
     .  Uatm, ucell(3,3), Uscf,
     .  vcell(3,3), virial, vn,
     .  volcel, volume, vpr, we, wmix, wmixkick, wo,
     .  xmax, xmin,
     .  stressl(3,3), veclen

      real(dp), dimension(:), allocatable, save ::
     .  efs, kweight, qs, r2ij, 
     .  wgthpol, dpbuffer1, H0, S

      real(dp), dimension(:,:), allocatable, save :: 
     .  auxpul, bk, cfa, fa, fal, kpoint, kpol,
     .  polR, polxyz, H, va, xij,
     .  xijo, dpbuffer2, Dscf, Dscfsave, Dold,
     .  Eold, Escf

#ifdef MPI
      real(dp)
     .  buffer1, stresstmp(3,3)

      real(dp), dimension(:,:), allocatable, save :: 
     .  fatmp
#endif

      real(dp), dimension(:,:,:), allocatable, save :: 
     .  ebk, eo, qo

      logical
     .  chebef, dminit, default, dumpcharge, 
     .  first, fixspin, found, foundxv, 
     .  gamma, inspn, itest, last, lastst, mix, mmix, negl, noeta, 
     .  outlng, overflow, pulfile, relaxd,
     .  same, savehs, savevh, savevt, savdrh, savrho,
     .  usesavecg, usesavelwf, usesavedm, usesavexv, 
     .  writbk, writmd, writpx, writb, writec, writef, 
     .  writek, writic, varcel, genlogic, do_pdos, writedm, atmonly,
     $  harrisfun

#ifdef MPI
      logical
     .  writei
#endif

      character
     .  filevh*25, filevt*25, fildrh*25, filrho*25,
     .  line*150, paste*25,
     .  slabel*20, sname*150, shape*10, message*79

      external
     .  anneal, chkdim, cgvc, fixed,
     .  dhscf, diagon, dnaefs, extrapol, initatom,
     .  iodm, ioxv, 
     .  kgrid, kgridinit, kinefsm, mulliken, naefs, neighb, nlefsm, 
     .  nose, npr, overfsm, paste, pulayx, 
     .  redata, reinit, shaper, spnvec, 
     .  timer, verlet2, volcel, xijorb, memory,
     .  ioeig, iofa, iokp, iomd, prversion

      type(parsed_line), pointer  :: p

      data
     .  e1, e2 / 1._dp, -1._dp /
     .  kcutof, kscell, kdispl / 0._dp, 9*0, 3*0._dp /
     .  nkpnt / 1 /
     .  relaxd /.false./
     .  tiny /1.e-15_dp/
     .  maxnh / 1 /
     .  maxna / 200 /
     .  tempion / 0.0_dp /
     .  no_l, nnamax /2*1/
     .  nauxpul, nbk, ns, nspin, nxij /5*1/
     .  ncells, nsc, mscell / 1,   1,1,1,   1,0,0, 0,1,0, 0,0,1 /
c---------------------------------------------------------------------

C Initialise MPI and set processor number
#ifdef MPI
      call MPI_Init( MPIerror )
      call MPI_Comm_Rank( MPI_Comm_World, Node, MPIerror )
      call MPI_Comm_Size( MPI_Comm_World, Nodes, MPIerror )
#else
      Node =  0
      Nodes = 1
#endif

      call io_setup     ! Initialize IOnode

C Print version information ...........................................
      if (IOnode) then
        call prversion
#ifdef MPI
         write(6,'(/,a,i4,a)')
     $        '* Running on ', Nodes, ' nodes in parallel'
#else
         write(6,'(/,a,i4,a)')
     $        '* Running in serial mode'
#endif
      endif
C ..................

C Start time counter ..................................................
      call timer( 'siesta', 0 )
      call timer( 'siesta', 1 )
      call timer( 'Setup', 1 )
C ..................

      if (IOnode)  call timestamp('Start of run')

C Factors of conversion to internal units (Bohr,fs,Ry) ................
      Ang    = 1._dp / 0.529177_dp
      eV     = 1._dp / 13.60580_dp
      kBar   = 1._dp / 1.47108e5_dp
      Kelvin = eV / 11604.45_dp
      Debye  = 0.393430_dp
      amu    = 2.133107_dp
C ..................

C Initialize some variables
      DUext = 0.0_dp
      Eharrs = 0.0_dp
      Eions = 0.0_dp
      Emad =0.0_dp
      FreeE = 0.0_dp

C Initialise read .....................................................
      call reinit(slabel, sname)

C Set allocation report level .........................................
      level = fdf_integer( 'alloc_report_level', 0 )
      call alloc_report( level=level, file=trim(slabel)//'.alloc' )
C ..................

! Initialize pseudopotentials and atomic orbitals

      if (IOnode) call initatom
      call broadcast_basis

      if (IOnode) atmonly = fdf_defined('Atom-Setup-Only')
#ifdef MPI
      call MPI_Bcast(atmonly,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
      if (atmonly) call bye("End of atom setup")

C Read simulation sizes ...............................................

!! Read number of atoms and coordinates, and unit cell
      if (IOnode) call coor(na_u,ucell)  ! Sets na_u, xa, and isa

#ifdef MPI
      call MPI_Bcast(na_u,1,MPI_integer,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(ucell(1,1),9,DAT_double,0,MPI_Comm_World,
     .  MPIerror)

      if (Node.ne.0) then    ! Allocate isa and xa in rest of nodes
         nullify(isa,xa)
         call realloc(isa,1,na_u,name='isa',routine='siesta')
         call realloc(xa,1,3,1,na_u,name='xa',routine='siesta')
      endif
      call MPI_Bcast(isa,na_u,MPI_integer,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(xa(1,1),3*na_u,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif

C Allocate arrays based on read sizes ................
!
!     Those living in module atomlists
!
      nullify(indxua,iza,lastkb,lasto,qa,amass,xalast)
      call realloc(indxua,1,na_u,name='indxua',routine='siesta')
      call realloc(iza,1,na_u,name='iza',routine='siesta')
      call realloc(lastkb,0,na_u,name='lastkb',routine='siesta')
      call realloc(lasto,0,na_u,name='lasto',routine='siesta')
      call realloc(qa,1,na_u,name='qa',routine='siesta')
      call realloc(xalast,1,3,1,na_u,name='xalast',routine='siesta')
      call realloc(amass,1,na_u,name='amass',routine='siesta')

!
!     Others
!
      allocate(cfa(3,na_u))
      call memory('A','D',3*na_u,'siesta')
      allocate(fa(3,na_u))
      call memory('A','D',3*na_u,'siesta')
      allocate(fal(3,na_u))
      call memory('A','D',3*na_u,'siesta')
      allocate(va(3,na_u))
      call memory('A','D',3*na_u,'siesta')

C Initialise those arrays that must be pre-initialised
      fal(1:3,1:na_u) = 0.0_dp
      stress(1:3,1:3) = 0.0_dp

      if (IOnode) call spin_init(nspin)
#ifdef MPI
      call MPI_Bcast(nspin,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif
      allocate(efs(nspin))
      call memory('A','D',nspin,'siesta')
      allocate(qs(nspin))
      call memory('A','D',nspin,'siesta')
      allocate(polR(3,nspin))
      call memory('A','D',3*nspin,'siesta')
      allocate(polxyz(3,nspin))
      call memory('A','D',3*nspin,'siesta')

C Read simulation data ................................................
      call redata(na_u, ns, nspin, outlng, 
     .            g2cut, charnet, negl, nscf, dDtol, mix, wmix, 
     .            isolve, temp, fixspin, ts, ncgmax, ftol, strtol, 
     .            eta, etol, rcoor, 
     .            ioptlwf, chebef, noeta, rcoorcp, beta, pmax,
     .            idyn, istart, ifinal, nmove, ianneal, iquench,
     .            dt, ia1, ia2, dx, dxmax, tt, tp, mn, mpr, 
     .            bulkm, taurelax,
     .            usesavelwf, usesavedm, usesavecg,
     .            mullipop, inspn, maxsav, nkick, wmixkick, 
     .            pulfile, tempinit, dumpcharge, varcel, harrisfun )


C Set ProcessorY and BlockSize in module parallel
#ifdef MPI
      if (IOnode) then
        call set_processorYdefault(Nodes,procYdefault)
        procYval = fdf_integer('processorY',procYdefault)
        BlockSize = fdf_integer('blocksize',8)
      endif
      call MPI_Bcast(procYval,1,MPI_integer,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(BlockSize,1,MPI_integer,0,MPI_Comm_World,MPIerror)
      call set_processorY(procYval)
      if (IOnode) then
         write(6,'(a,2i4,/)')
     $        '* ProcessorY, Blocksize: ', procYval, Blocksize
      endif
#else
      call set_processorY(1)
      BlockSize = 8
#endif

C Find some switches ..................................................
      if (IOnode) then
        writek    = fdf_boolean('WriteKpoints'    , outlng )
        writef    = fdf_boolean('WriteForces'     , outlng )
        writedm   = fdf_boolean('WriteDM'         , .true. )
#ifdef MPI
        writei    = fdf_boolean('WriteEigenvalues', outlng )
#endif
        writb     = fdf_boolean('WriteBands'      , outlng )
        writbk    = fdf_boolean('WriteKbands'     , outlng )
        writec    = fdf_boolean('WriteCoorStep'   , outlng )
        writic    = fdf_boolean('WriteCoorInitial', .true. )
        writmd    = fdf_boolean('WriteMDhistory'  , .false.)
        writpx    = fdf_boolean('WriteMDXmol'     , .not. writec)
        default   = fdf_boolean('UseSaveData'     , .false.)
        usesavexv = fdf_boolean('MD.UseSaveXV'    , default)
        savehs    = fdf_boolean('SaveHS'          , .false.)

        rijmin    = fdf_physical( 'WarningMinimumAtomicDistance',
     .                                              1._dp, 'Bohr' )

      endif
#ifdef MPI
      call MPI_Bcast(writek,1,MPI_logical,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(writef,1,MPI_logical,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(writedm,1,MPI_logical,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(writei,1,MPI_logical,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(writb,1,MPI_logical,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(writbk,1,MPI_logical,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(writec,1,MPI_logical,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(writic,1,MPI_logical,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(writmd,1,MPI_logical,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(writpx,1,MPI_logical,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(default,1,MPI_logical,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(usesavexv,1,MPI_logical,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(savehs,1,MPI_logical,0,MPI_Comm_World,MPIerror)

      call MPI_Bcast(rijmin,1,DAT_double,0,MPI_Comm_World,MPIerror)
#endif

C .....................

C Read cell shape and atomic positions from a former run ..............
      foundxv = .false.
      if (usesavexv) then
        call ioxv('read', ucell, vcell, na_u, isa, iza, xa, va, foundxv)
        if (IOnode.and..not.foundxv) write(6,'(/,a)')
     .    'siesta: WARNING: XV file not found'
      endif
C ..................

C Dump initial coordinates to output ..................................
      if ( writic.and.(IOnode) ) then
        write(6,'(/a)') 'siesta: Atomic coordinates (Bohr) and species'
        write(6,"('siesta: ',2x,3f10.5,i3,3x,i4)")
     .           ( (xa(ix,ia), ix=1,3), isa(ia), ia, ia=1, na_u)
      endif
C ..................


C Initialize atom lists 

      call initatomlists
      qtot = qtot - charnet
!

C Find number of locally stored orbitals and allocated related arrays
      call GetNodeOrbs(no_u,Node,Nodes,no_l)

      allocate(listhptr(no_l))
      call memory('A','I',no_l,'siesta')
      listhptr(:) = 0
      allocate(listhptrold(no_l))
      call memory('A','I',no_l,'siesta')
      listhptrold(:) = 0
      allocate(numh(no_l))
      call memory('A','I',no_l,'siesta')
      numh(:) = 0
      allocate(numhold(no_l))
      call memory('A','I',no_l,'siesta')
      numhold(:) = 0

!

C  calculate spin populations for fixed spin case...
      if (fixspin) then
        if (nspin .ne. 2) then
           write(message,'(2a)')
     $       'siesta: ERROR: You can only fix the spin of the system',
     $       'for collinear spin polarized calculations.'
             call die(message)
        endif
        do i=1,2
          qs(i) = (qtot + (3-2*i)*ts) / 2.0_dp
        enddo
      endif
C ..................


C Find maximum interaction range ......................................
      if (negl) then
        rmaxh = 2.0_dp*rmaxo
      else
        rmaxh = 2.0_dp*rmaxo + 2.0_dp*rmaxkb
      endif
C ......................

C Automatic cell generation ...........................................
      if (volcel(ucell) .lt. 1.0d-8) then
        ucell(1:3,1:3) = 0.0_dp
        scell(1:3,1:3) = 0.0_dp
        do ix = 1,3
          xmin =  huge(1._dp)
          xmax = -xmin
          do ia = 1,na_u
            is = isa(ia)
            rc = rcut(is,0)
            xmin = min( xmin, xa(ix,ia)-rc )
            xmax = max( xmax, xa(ix,ia)+rc )
          enddo
C         Use a 10% margin for atomic movements
          ucell(ix,ix) = 1.10_dp * (xmax - xmin)
          scell(ix,ix) = ucell(ix,ix)
        enddo
C build cubic cell if system is charged
        if (charnet .ne. 0.0_dp) then
          xmax = -huge(1._dp)
          do ix = 1,3
            if (ucell(ix,ix) .gt. xmax) xmax = ucell(ix,ix)
          enddo
          do ix = 1,3
            ucell(ix,ix) = xmax
            scell(ix,ix) = xmax
          enddo
        endif
C
        volume = volcel( ucell )
        if (IOnode) then
          write(6,'(/,a,3(/,a,3f12.6))')
     .      'siesta: Automatic unit cell vectors (Ang):',
     .      ('siesta:', (ucell(ix,iv)/Ang,ix=1,3), iv =1,3)
        endif
      endif
C ......................

C Find system shape ...................................................
      call shaper( ucell, na_u, isa, xa, shape, nbcell, bcell )
      if (IOnode) then
        write(6,'(/,2a)') 'siesta: System type = ', shape
      endif
C ......................

C Madelung correction for charged systems .............................
      if (charnet .ne. 0.0_dp) then
        call madelung(ucell, shape, charnet, Emad)
      endif
C ......................

C Find k-grid for Brillouin zone integration ..........................
C Call kgrid once to get the dimensions and a second time to actually
C initialise the data - could be made more efficient in future!
      nkpnt = 1
      call kgridinit( ucell, kscell, kdispl, kcutof, nkpnt )
      maxk=nkpnt
      allocate(kpoint(3,maxk))
      call memory('A','D',3*maxk,'siesta')
      allocate(kweight(maxk))
      call memory('A','D',maxk,'siesta')
      kpoint(1:3,1)=0.0_dp
      kweight(1)=1.0_dp
      allocate(eo(no_u,nspin,maxk))
      call memory('A','D',no_u*nspin*maxk,'siesta')
      allocate(qo(no_u,nspin,maxk))
      call memory('A','D',no_u*nspin*maxk,'siesta')
      call kgrid( ucell, kscell, kdispl,
     .            nkpnt, kpoint, kweight )
C ......................

C Find number of band k-points ........................................
      nbk = 0
      maxbk = 1
      allocate(bk(3,maxbk))
      call memory('A','D',3*maxbk,'siesta')
C
      call initbands( maxbk, nbk, bk )
C
      if (nbk .gt. maxbk) then
C If there wasn't enough space to store bands on first call correct
C the dimensions and repeat the initialisation
        maxbk=max(nbk,1)
        call memory('D','D',size(bk),'siesta')
        deallocate(bk)
        allocate(bk(3,maxbk))
        call memory('A','D',3*maxbk,'siesta')
        nbk = 0
        call initbands( maxbk, nbk, bk )
      endif
      allocate(ebk(no_u,nspin,maxbk))
      call memory('A','D',no_u*nspin*maxbk,'siesta')

C ......................

C Find the grid for the calculation of the polarization..............
      nkpol=1
      allocate(kpol(3,1))
      call memory('A','D',3,'siesta')
      allocate(wgthpol(1))
      call memory('A','D',1,'siesta')
      call KSV_init(ucell, 0, nkpol, kpol, wgthpol)
      call memory('D','D',size(kpol),'siesta')
      deallocate(kpol)
      call memory('D','D',size(wgthpol),'siesta')
      deallocate(wgthpol)
      if (nkpol.gt.0) then
        allocate(kpol(3,nkpol))
        call memory('A','D',3*nkpol,'siesta')
        allocate(wgthpol(nkpol))
        call memory('A','D',nkpol,'siesta')
      endif

C Find if only gamma point is used ....................................
      if (nkpnt.eq.1 .and. abs(kpoint(1,1)).lt.tiny .and.
     .                     abs(kpoint(2,1)).lt.tiny .and.
     .                     abs(kpoint(3,1)).lt.tiny) then
        gamma = .true.
      else
        gamma = .false.
      endif
      if (nbk .gt. 0) gamma = .false.
      if (nkpol.gt.0) gamma = .false.
C ....................

C Print k-points ......................................................
      if (.not.gamma .and. IOnode) then
        if ( writek ) then
          write(6,'(/,a)')
     .     'siesta: k-point coordinates (Bohr**-1) and weights:'
          write(6,'(a,i4,3f12.6,3x,f12.6)')
     .      ('siesta: ', ik, (kpoint(ix,ik),ix=1,3), kweight(ik),
     .      ik=1,nkpnt)
        else
          call iokp( nkpnt, kpoint, kweight )
        endif
        write(6,'(/a,i6)')
     .    'siesta: k-grid: Number of k-points =', nkpnt
        write(6,'(a,f10.3,a)')
     .    'siesta: k-grid: Cutoff             =', kcutof/Ang, ' Ang'
        write(6,'(a)')
     .    'siesta: k-grid: Supercell and displacements'
        write(6,'(a,3i4,3x,f8.3)') 'siesta: k-grid: ',
     .                             (kscell(i,1),i=1,3), kdispl(1)
        write(6,'(a,3i4,3x,f8.3)') 'siesta: k-grid: ',
     .                             (kscell(i,2),i=1,3), kdispl(2)
        write(6,'(a,3i4,3x,f8.3)') 'siesta: k-grid: ',
     .                             (kscell(i,3),i=1,3), kdispl(3)
      endif
C ....................

C Find required supercell
      do i = 1,3
        mscell(1:3,i) = 0.0_dp
        if (gamma) then
          nsc(i) = 1
          mscell(i,i) = 1
        else
          veclen = sqrt( ucell(1,i)**2 + ucell(2,i)**2 + ucell(3,i)**2 )
          nsc(i) = 1 + 2 * rmaxh / veclen
          mscell(i,i) = nsc(i)
        endif
      enddo

C Find auxiliary supercell (required only for k sampling) ............

      call superc( ucell, scell, nsc)

C Initialize atomic velocities to zero ................................
      if (.not. foundxv) then
        va(1:3,1:na_u) = 0.0_dp
        vcell(1:3,1:3)=0.0_dp
      endif
C ..................

C Begin of coordinate relaxation iteration ============================
C Notice that this loop is not indented
      if (idyn .eq. 0) then
        inicoor = 0
        fincoor = nmove
      else if (idyn .ge. 1 .and. idyn .le. 5) then
        inicoor = istart
        fincoor = ifinal
      else if (idyn .eq. 6) then
        inicoor = 0
        fincoor = (ia2-ia1+1)*3*2
      else if (idyn .eq. 7) then
        call phonon_setup
        inicoor = 1
        fincoor = phonon_num_disps
      else
         call die('siesta: wrong idyn')
      endif

C Build initial velocities according to Maxwell-Bolzmann distribution....
      if (idyn .ne. 0 .and. idyn .ne. 6 .and. (.not. foundxv)) 
     .    call vmb(na_u,tempinit,amass,va)
C ..................

      istp = 0
      call timer( 'Setup', 2 )

C Output memory use before main loop
      call printmemory( 6, 0 )

C Start loop over coordinate changes 

      do istep = inicoor,fincoor
      call timer( 'IterMD', 1 )
      istp = istp + 1
      if (IOnode) then
        write(6,'(/,a)') 'siesta:    ==============================='
        if (idyn .eq. 0) 
     .   write(6,'(a,i6)') 'siesta:        Begin CG move = ',istep
        if (idyn .gt. 0 .and. idyn .lt. 6) 
     .   write(6,'(a,i6)') 'siesta:        Begin MD step = ',istep
        if (idyn .eq. 6)  then
          write(6,'(a,i6)') 'siesta:        Begin FC step = ',istep
          if (istep .eq. 0) then
            write(6,'(a)') 'siesta:        Undisplaced coordinates'
          else
            iadispl = (istep-mod(istep-1,6))/6+ia1
            write(6,'(a,i6)') 'siesta:        displace atom = ',
     .        iadispl
            ix = mod(istep-1,6)+1
            ixdispl = (ix - mod(ix-1,2) +1)/2
            write(6,'(a,i6)') 'siesta:        in direction  = ',
     .        ixdispl
            dx=-dx
            write(6,'(a,f8.4,a)') 'siesta:        by            = ',
     .                      dx, ' Bohr'
c displace atom by dx...
            xa(ixdispl,iadispl)=xa(ixdispl,iadispl)+dx
          endif
        endif
        if (idyn .eq. 7)  then
           call phonon_set_coords(istep,xa,ucell)
        endif
        write(6,'(a)')   'siesta:    ==============================='
      else
        if (idyn .eq. 6)  then
          if (istep .ne. 0) then
            iadispl = (istep-mod(istep-1,6))/6+ia1
            ix = mod(istep-1,6)+1
            ixdispl = (ix - mod(ix-1,2) +1)/2
            dx=-dx
c displace atom by dx...
            xa(ixdispl,iadispl)=xa(ixdispl,iadispl)+dx
          endif
        else if (idyn .eq. 7)  then
           call phonon_set_coords(istep,xa,ucell)
        endif
      endif

      if (IOnode) then
C Print atomic coordinates ............................................
        call outcoor( ucell, xa, isa, na_u, ' ', writec )
      endif
C ...................


C Actualize things if variable cell ...................................

      if ( varcel .and. (istep.ne.inicoor) .and. (.not.gamma) ) then

C k-grid 
        nkpnt = maxk
        call kgridinit( ucell, kscell, kdispl, kcutof, nkpnt )

C If number of k points has changed - re-size arrays
        if (nkpnt .gt. maxk) then
          call memory('D','D',size(kpoint),'siesta')
          deallocate(kpoint)
          call memory('D','D',size(kweight),'siesta')
          deallocate(kweight)
          allocate(kpoint(3,nkpnt))
          call memory('A','D',3*nkpnt,'siesta')
          allocate(kweight(nkpnt))
          call memory('A','D',nkpnt,'siesta')
          kpoint(1:3,1)=0.0_dp
          kweight(1)=1.0_dp
          call kgrid( ucell, kscell, kdispl,
     .      nkpnt, kpoint, kweight )
          call memory('D','D',size(eo),'siesta')
          deallocate(eo)
          call memory('D','D',size(qo),'siesta')
          deallocate(qo)
          allocate(eo(no_u,nspin,nkpnt))
          call memory('A','D',no_u*nspin*nkpnt,'siesta')
          allocate(qo(no_u,nspin,nkpnt))
          call memory('A','D',no_u*nspin*nkpnt,'siesta')
          maxk=nkpnt
       else
          call kgrid( ucell, kscell, kdispl,
     $                nkpnt, kpoint, kweight )

       endif
 
C Print k-points 
        if (IOnode) then
          if ( writek ) then
            write(6,'(/,a)')
     .       'siesta: k-point coordinates (Bohr**-1) and weights:'
            write(6,'(a,i4,3f12.6,3x,f12.6)')
     .        ('siesta: ', ik, (kpoint(ix,ik),ix=1,3), kweight(ik),
     .        ik=1,nkpnt)
          else
            call iokp( nkpnt, kpoint, kweight )
          endif
          write(6,'(/a,i6)')
     .      'siesta: k-grid: Number of k-points =', nkpnt
          write(6,'(a,f10.3,a)')
     .      'siesta: k-grid: Cutoff             =', kcutof/Ang, ' Ang'
          write(6,'(a)')
     .      'siesta: k-grid: Supercell and displacements'
          write(6,'(a,3i4,3x,f8.3)') 'siesta: k-grid: ',
     .                               (kscell(i,1),i=1,3), kdispl(1)
          write(6,'(a,3i4,3x,f8.3)') 'siesta: k-grid: ',
     .                               (kscell(i,2),i=1,3), kdispl(2)
          write(6,'(a,3i4,3x,f8.3)') 'siesta: k-grid: ',
     .                               (kscell(i,3),i=1,3), kdispl(3)
        endif
 
C Find required supercell
        do i = 1,3
          mscell(1:3,i) = 0.0_dp
          if (gamma) then
            nsc(i) = 1
            mscell(i,i) = 1
          else
            veclen = sqrt(ucell(1,i)**2+ucell(2,i)**2+ucell(3,i)**2)
            nsc(i) = 1 + 2 * rmaxh / veclen
            mscell(i,i) = nsc(i)
          endif
        enddo

C Madelung correction for charged systems .............................
        if (charnet .ne. 0.0_dp) then
          call madelung(ucell, shape, charnet, Emad)
        endif

      endif
C End variable cell actualization

C auxiliary supercell

        call superc(ucell,  scell, nsc)


C Print unit cell for variable cell ...................................
      if ( varcel.and.(IOnode) ) call outcell(ucell)
C ...................

C Find unit cell volume ...............................................
      volume = volcel( ucell )
C ...................

C Initialize neighb subroutine ........................................
  144 ia = 0
      isel = 0
      rmax = max( 2._dp*rmaxv, 2._dp*rmaxo, rmaxo+rmaxkb )
      nnia = maxna
      if (allocated(jna)) then
        call memory('D','I',size(jna),'siesta')
        deallocate(jna)
      endif
      if (allocated(r2ij)) then
        call memory('D','D',size(r2ij),'siesta')
        deallocate(r2ij)
      endif
      if (allocated(xij)) then
        call memory('D','D',size(xij),'siesta')
        deallocate(xij)
      endif
      allocate(jna(maxna))
      call memory('A','I',maxna,'siesta')
      allocate(r2ij(maxna))
      call memory('A','D',maxna,'siesta')
      allocate(xij(3,maxna))
      call memory('A','D',3*maxna,'siesta')
      call neighb( scell, rmax, na_s, xa, ia, isel,
     .             nnia, jna, xij, r2ij )
      nnamax = 0
      do ia = 1,na_s
        nnia = 0
        call neighb( scell, rmax, na_s, xa, ia, isel,
     .               nnia, jna, xij, r2ij )
        nnamax = max( nnamax, nnia )
      enddo
      if (nnamax .gt. maxna) then
C       Increase maxna with safety margin when atoms move
        maxna = nnamax + 0.10 * nnamax + 10
        overflow = .true.
      else
        overflow = .false.
      endif
      if (overflow) goto 144
C ..................

C Check if any two atoms are unreasonably close .......................
      do ia = 1,na_s
        r2min = huge(1._dp)
        jamin = 0
        nnia = maxna
        call neighb( scell, rmax, na_s, xa, ia, isel,
     .               nnia, jna, xij, r2ij )
        do j = 1,nnia
          ja = jna(j)
          if ( r2ij(j).lt.r2min .and. ja.ge.ia ) then
C           Check that it is not the same atom
            if ( ja.ne.ia .or. r2ij(j).gt.1.d-12 ) then
              r2min = r2ij(j)
              jamin = ja
            endif
          endif
        enddo
        rmin = sqrt( r2min )
        if (IOnode) then
          if ( rmin .lt. rijmin ) write(6,'(a,2i6,a,f12.6,a)')
     .      'siesta: WARNING: Atoms', ia, jamin, ' too close: rij =',
     .       rmin/Ang, ' Ang'
        endif
      enddo
C ..................

C List of nonzero hamiltonian matrix elements .........................
      if (.not.allocated(listh)) then
        allocate(listh(maxnh))
        call memory('A','I',maxnh,'siesta')
      endif
   44 nh = maxnh
      maxnhold = maxnh
      call hsparse( negl, scell, nsc, na_s, isa, xa, lasto, lastkb, 
     .              iphorb, iphKB, nh, numh, listhptr, listh, 
     .              Node, Nodes )
      overflow = .false.
      if (nh .gt. maxnh) then
C Increase maxnh with safety margin when atoms move
        maxnh = 1.05 * nh + 40
        overflow = .true.
      endif
      if (overflow) then
C For the first step there is no need to preserve the contents
C but for subsequent steps listhold, Dscf and Dscfsave must be conserved
        if (istp.eq.1) then
          if (allocated(listhold)) then
            call memory('D','I',size(listhold),'siesta')
            deallocate(listhold)
          endif
          if (allocated(Dscfsave)) then
            call memory('D','D',size(Dscfsave),'siesta')
            deallocate(Dscfsave)
          endif
          if (allocated(Dscf)) then
            call memory('D','D',size(Dscf),'siesta')
            deallocate(Dscf)
          endif
          allocate(listhold(maxnh))
          listhold(1:maxnh) = 0
          call memory('A','I',maxnh,'siesta')
          allocate(Dscfsave(maxnh,nspin))
          call memory('A','D',maxnh*nspin,'siesta')
          allocate(Dscf(maxnh,nspin))
          call memory('A','D',maxnh*nspin,'siesta')
C Initialise Dscfsave to avoid problems in extrapol
          Dscfsave(1:maxnh,1:nspin) = 0.0_dp
        else
C
          allocate(ibuffer1(maxnhold))
          call memory('A','I',maxnhold,'siesta')
          ibuffer1(1:maxnhold)=listhold(1:maxnhold)
          if (allocated(listhold)) then
            call memory('D','I',size(listhold),'siesta')
            deallocate(listhold)
          endif
          allocate(listhold(maxnh))
          listhold(1:maxnh) = 0
          call memory('A','I',maxnh,'siesta')
          listhold(1:maxnhold)=ibuffer1(1:maxnhold)
          call memory('D','I',size(ibuffer1),'siesta')
          deallocate(ibuffer1)
C
          allocate(dpbuffer2(maxnhold,nspin))
          call memory('A','D',maxnhold*nspin,'siesta')
          dpbuffer2(1:maxnhold,1:nspin)=Dscfsave(1:maxnhold,1:nspin)
          if (allocated(Dscfsave)) then
            call memory('D','D',size(Dscfsave),'siesta')
            deallocate(Dscfsave)
          endif
          allocate(Dscfsave(maxnh,nspin))
          call memory('A','D',maxnh*nspin,'siesta')
          Dscfsave(1:maxnhold,1:nspin)=dpbuffer2(1:maxnhold,1:nspin)
          dpbuffer2(1:maxnhold,1:nspin)=Dscf(1:maxnhold,1:nspin)
          if (allocated(Dscf)) then
            call memory('D','D',size(Dscf),'siesta')
            deallocate(Dscf)
          endif
          allocate(Dscf(maxnh,nspin))
          call memory('A','D',maxnh*nspin,'siesta')
          Dscf(1:maxnhold,1:nspin)=dpbuffer2(1:maxnhold,1:nspin)
          call memory('D','D',size(dpbuffer2),'siesta')
          deallocate(dpbuffer2)
C
        endif
        if (allocated(Dold)) then
          call memory('D','D',size(Dold),'siesta')
          deallocate(Dold)
        endif
        if (allocated(Eold)) then
          call memory('D','D',size(Eold),'siesta')
          deallocate(Eold)
        endif
        if (allocated(Escf)) then
          call memory('D','D',size(Escf),'siesta')
          deallocate(Escf)
        endif
        if (allocated(listh)) then
          call memory('D','I',size(listh),'siesta')
          deallocate(listh)
        endif
        allocate(Dold(maxnh,nspin))
        call memory('A','D',maxnh*nspin,'siesta')
        allocate(Eold(maxnh,nspin))
        call memory('A','D',maxnh*nspin,'siesta')
        allocate(Escf(maxnh,nspin))
        call memory('A','D',maxnh*nspin,'siesta')
        allocate(listh(maxnh))
        call memory('A','I',maxnh,'siesta')
        goto 44
      endif

C Allocate/reallocate storage associated with Hamiltonian/Overlap matrix
      if (allocated(S).and.size(S).lt.maxnh) then
        call memory('D','D',size(H),'siesta')
        deallocate(H)
        call memory('D','D',size(H0),'siesta')
        deallocate(H0)
        call memory('D','D',size(S),'siesta')
        deallocate(S)
      endif
      if (.not.allocated(S)) then
        allocate(H(maxnh,nspin))
        call memory('A','D',maxnh*nspin,'siesta')
        allocate(H0(maxnh))
        call memory('A','D',maxnh,'siesta')
        allocate(S(maxnh))
        call memory('A','D',maxnh,'siesta')
      endif

C ..................

C Some printout for debugging ........................................
*     if (IOnode) then
*       write(6,'(/,a)') 'siesta: connected orbitals'
*       do io = 1,no_u
*         call GlobalToLocalOrb(io,Node,Nodes,iio)
*         if (iio.gt.0) then
*           write(6,'(i6,4x,15i4)') 
*    .        io, (listh(listhptr(iio)+j),j=1,numh(iio))
*         endif
*#ifdef MPI
*         call MPI_Barrier(MPI_Comm_World,MPIerror)
*#endif
*       enddo
*       write(6,*) ' '
*     endif
C ..................

C Find vectors between orbital centers ................................
      if (allocated(xijo)) then
        call memory('D','D',size(xijo),'siesta')
        deallocate(xijo)
      endif
      if (.not.gamma) then
        nxij = maxnh
        allocate(xijo(3,nxij))
        call memory('A','D',3*nxij,'siesta')
        call xijorb( negl, scell, na_u, na_s, xa,
     .               lasto, lastkb, rco, rckb,
     .               maxnh, numh, listhptr,
     .               listh, xijo, Node, Nodes )
      else
        nxij = 1
        allocate(xijo(3,1))
        call memory('A','D',3,'siesta')
      endif
C ..................

C Initialize density matrix ...........................................
C set density matrix for first step
      found = .false.
      dminit = .false.
      if (istp .eq. 1) dminit = .true.
      if (istp .ne. 1 .and. harrisfun) dminit = .true.
      if (istp .ne. 1 .and. (idyn .eq. 6)
     $         .and. usesavedm .and. writedm)  dminit = .true.
      if (istp .ne. 1 .and. (idyn .eq. 7)
     $         .and. usesavedm)  dminit = .true.

      if (dminit)
     .   call initdm(Datm, Dscf, Dold, lasto, na_s,
     .               maxnh, no_s, no_l, nspin, na_u, no_l, nspin,
     .               numh, numhold, listhptr, listhptrold,
     .               listh, listhold, iaorb, found, inspn, 
     .               usesavedm, no_u)


C Initialize energy-density matrix to zero for first call to overfsm
      Escf = 0.0_dp

C Extrapolate density matrix between steps
      itest = .false.
      istpsave = 0
      iord = 1
      if (idyn .eq. 0) iord = 0
      if (idyn .eq. 6) iord = 0
      if (idyn .eq. 7) iord = 0
C  If DM has just been read from disk, 
C  call extrapol with istep = 2 and iord = 0
C  to make it update the structure of DM, if needed
      if (found .and. ((istp .eq. 1) .or. (idyn .eq. 6)
     $                               .or. (idyn .eq. 7))) then
        istpsave = istp
        istp = 2
        iord = 0
        itest = .true.
      endif
      if( .not. harrisfun)
     . call extrapol(istp, iord, nspin, no_s, no_u, no_l, maxnh, 
     .              numh, listhptr, listh, numhold, listhptrold,
     .              listhold, Dscfsave, Dscf)
C  If DM have just been read, restore istp
      if (itest) istp = istpsave
      itest = .false.
C ..................

C Check for Pulay auxiliary matrices sizes ...................................
      if (pulfile .or. maxsav .le. 0) then
        nauxpul = 1
        if (.not.allocated(auxpul)) then
          allocate(auxpul(nauxpul,2))
          call memory('A','D',2*nauxpul,'siesta')
        endif
      else
        nauxpul = 0
        do io = 1,no_l
          nauxpul = nauxpul + numh(io)
        enddo
        nauxpul = nauxpul * nspin * maxsav
#ifdef MPI
        call MPI_AllReduce(nauxpul,ntmp,1,MPI_integer,MPI_max,
     .    MPI_Comm_World,MPIerror)
        nauxpul = ntmp
#endif
C       Increase nauxpul with safety margin when atoms move
        nauxpul = 1.1 * nauxpul + 10
        if (allocated(auxpul)) then
          if (size(auxpul,1).ne.nauxpul) then
            call memory('D','D',size(auxpul),'siesta')
            deallocate(auxpul)
            allocate(auxpul(nauxpul,2))
            call memory('A','D',2*nauxpul,'siesta')
          endif
        else
          allocate(auxpul(nauxpul,2))
          call memory('A','D',2*nauxpul,'siesta')
        endif
      endif
C ....................

C Find overlap matrix ...............................................
      call overfsm(na_u, na_s, no_s, scell, xa, indxua, rmaxo, no_l,
     .             maxna, maxnh, maxnh, lasto, iphorb, isa, 
     .             numh, listhptr, listh, numh, listhptr, listh, 
     .             min(nspin,2), Escf, jna, xij, r2ij,
     .             fal, stress, S, Node, Nodes)
C ..................

C Start of SCF iteration _____________________________________________
      first = .true.
      last  = .false.
      if (wmix .le. 0._dp) then
        if (IOnode) then
          write(6,'(/,a,f15.8)')
     .     'siesta: WARNING: Mixing weight for SCF loop =', wmix
        endif
        last = .true.
      endif

      do iscf = 1, nscf
        if (iscf .eq. nscf) last = .true.
        call timer( 'IterSCF', 1 )

C Normalize density matrix to exact charge ...........................
        qsol = 0.0_dp
        do ispin = 1,min(nspin,2)
          do io = 1,nh
            qsol = qsol + Dscf(io,ispin) * s(io)
          enddo
        enddo
#ifdef MPI
        call MPI_AllReduce(qsol,buffer1,1,DAT_double,
     .    MPI_sum,MPI_Comm_World,MPIerror)
        qsol=buffer1
#endif
        if (IOnode) then
          if (.not.first .and.
     .       abs(qsol/qtot-1._dp).gt.1.d-2) write(6,'(a,2f15.6)')
     .      'siesta: WARNING: Qtot, Tr[D*S] =', qtot, qsol
        endif
        do ispin = 1,nspin
          do io = 1,nh
            Dscf(io,ispin) = Dscf(io,ispin) * qtot/qsol
            Escf(io,ispin) = Escf(io,ispin) * qtot/qsol
          enddo
        enddo
C ..................

C Initialize Hamiltonian ........................................
        H = 0.0_dp

C Initialize forces and stress ...................
        if (first.or.last) then
          fa(1:3,1:na_u) = 0.0_dp
          fal(1:3,1:na_u) = 0.0_dp
          stress(1:3,1:3) = 0.0_dp
          stressl(1:3,1:3) = 0.0_dp
        endif
C ..................

C Self-energy of isolated ions ........................................
        if (first) then
          Eions = 0.0_dp
          do ia = 1,na_u
            is = isa(ia)
            Eions = Eions + uion(is)
          enddo
        endif
C ..................

C Neutral-atom: energy, forces and stress ............................
C First time for energy, last time for forces
        if (first.or.last) then
          call naefs(na_u, na_s, scell, xa, indxua, rmaxv,
     .               maxna, isa, jna, xij, r2ij,
     .               Ena, fa, stress)
          call dnaefs(na_u, na_s, scell, xa, indxua, rmaxv,
     .               maxna, isa, jna, xij, r2ij,
     .               DEna, fa, stress) 
          Ena = Ena + DEna
        endif
C ..................

C Kinetic: energy, forces, stress and matrix elements .................
        if (first.or.last) then
          call kinefsm(na_u, na_s, no_s, scell, xa, indxua, rmaxo, no_l,
     .                 maxna, maxnh, maxnh, lasto, iphorb, isa, 
     .                 numh, listhptr, listh, numh, listhptr, listh, 
     .                 min(nspin,2), Dscf, jna, xij, r2ij,
     .                 Ekin, fal, stressl, H, Node, Nodes) 
#ifdef MPI
C Global reduction of energy terms
          call MPI_AllReduce(Ekin,buffer1,1,DAT_double,
     .      MPI_sum,MPI_Comm_World,MPIerror)
          Ekin = buffer1
#endif
        endif
C ..................

C Non-local-pseudop: energy, forces, stress and matrix elements .......
        if (first.or.last) then
          call nlefsm(scell, na_u, na_s, isa, xa, indxua, maxnh,
     .                maxnh, lasto, lastkb, iphorb, iphKB, numh,
     .                listhptr, listh, numh, listhptr, listh, 
     .                min(nspin,2), Dscf, Enl, fal, stressl, H,
     .                Node, Nodes) 
#ifdef MPI
C Global reduction of energy terms
          call MPI_AllReduce(Enl,buffer1,1,DAT_double,MPI_sum,
     .      MPI_Comm_World,MPIerror)
          Enl = buffer1
#endif
        endif
C ..................

C Save or get partial Hamiltonian (non-SCF part) ......................
        if (first.or.last) then
          H0(1:nh) = H(1:nh,1)
        else
          do ispin = 1,nspin
            if (ispin .le. 2) then
              H(1:nh,ispin) = H0(1:nh)
            else
              H(1:nh,ispin) = 0.0_dp
            endif
          enddo          
        endif
C ..................

C Non-SCF part of total energy .......................................
        if (first.or.last) then
          E0 = -Eions + Ena + Ekin + Enl
        else
          E0 = 0.0_dp
          do ispin = 1,min(nspin,2)
            do io = 1,nh
              E0 = E0 + H0(io) * Dscf(io,ispin)
            enddo
          enddo
#ifdef MPI
C Global reduction of E0
          call MPI_AllReduce(E0,buffer1,1,DAT_double,
     .      MPI_sum,MPI_Comm_World,MPIerror)
          E0 = buffer1
#endif
          E0 = E0 - Eions + Ena
        endif
C ..................

C Non-local-pseudop: energy, forces, stress and matrix elements .......
C Add SCF contribution to energy and matrix elements ..................
        g2max = g2cut
        if (last) then
c         Last call to dhscf and grid-cell sampling if requested
          ifa  = 1
          istr = 1
          call grdsam( nspin, no_s, iaorb, iphorb, indxuo,
     .                 no_l, no_u, na_u, na_s, isa, xa, indxua,
     .                 ucell, mscell, g2max, ntm, ifa, istr, maxnh,
     .                 maxnh, numh, listhptr, listh, Dscf, Datm, H,
     .                 Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     .                 Exc, Dxc, dipol, fa, stress, fal, stressl)
        else
          ifa  = 0
          istr = 0
          ihmat = 1
          call dhscf( nspin, no_s, iaorb, iphorb, indxuo, no_l,
     .                no_u, na_u, na_s, isa, xa, indxua, 
     .                ucell, mscell, g2max, ntm,
     .                ifa, istr, ihmat, ' ', ' ', ' ', ' ',
     .                maxnh, numh, listhptr, listh, Dscf, Datm,
     .                maxnh, numh, listhptr, listh, H,
     .                Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     .                Exc, Dxc, dipol, fa, stress, fal, stressl)
        endif
            
C Output memory use after first call to dhscf
        if (istp.eq.1 .and. iscf.eq.1) call printmemory( 6, 0 )

*       if (istp.eq.1 .and. iscf.eq.1) write(6,'(/,a,f10.3,a)')
*    .    'siesta: dhscf mesh cutoff =', g2max, ' Ry'

C ..................

C Orthonormalization forces ...........................................
        if (last) then
           call overfsm(na_u, na_s, no_s, scell, xa, indxua,rmaxo,no_l,
     .                 maxna, maxnh, maxnh, lasto, iphorb, isa, 
     .                 numh, listhptr, listh, numh, listhptr, listh, 
     .                 min(nspin,2), Escf, jna, xij, r2ij,
     .                 fal, stressl, S, Node, Nodes) 
        endif
C ..................

C Find entropy ........................................................
        Entrop = 0.0_dp
        if (isolve .eq. 0) then
          if (istp.gt.1 .or. iscf.gt.1) then
            do ik = 1,nkpnt
              do io = 1,no_u
                if (nspin .eq. 1) then
                  wo = qo(io,1,ik) / kweight(ik) * 0.5_dp
                  we = 1._dp - wo
                  wo = max( wo, tiny )
                  we = max( we, tiny )
                  Entrop = Entrop - 2._dp * kweight(ik) *
     .                             ( wo*log(wo) + we*log(we) )
                else
                  do ispin = 1,2
                    wo = qo(io,ispin,ik) / kweight(ik)
                    we = 1._dp - wo
                    wo = max( wo, tiny )
                    we = max( we, tiny )
                    Entrop = Entrop - kweight(ik) *
     .                               ( wo*log(wo) + we*log(we) )
                  enddo
                endif
              enddo
            enddo
          endif
        endif

C Save present density matrix ........................................
        Dold(1:nh,1:nspin) = Dscf(1:nh,1:nspin)
        Eold(1:nh,1:nspin) = Escf(1:nh,1:nspin)

C Save hamiltonian and overlap matrices ............................
        if (savehs) then
          call iohs( 'write', gamma, no_u, nspin, maxnh,
     .               numh, listhptr, listh, H, S, qtot, temp, xijo )
        endif
C ..................

C Solve eigenvalue problem .........................................
        if( .not. last) then
          if (isolve .eq. 0) then
            call diagon(no_s, nspin, nspin, no_l, maxnh, maxnh, no_u,
     .                  numh, listhptr, listh, numh, listhptr, listh, 
     .                  H, S, qtot, fixspin, qs, temp, e1, e2,
     .                  gamma, xijo, indxuo, nkpnt, kpoint, kweight,
     .                  eo, qo, Dscf, Escf, ef, efs, no_u)
            Ecorrec = 0._dp
          elseif (isolve .eq. 1) then
            if (.not. gamma) call die("Cannot do O(N) with k-points.")
            call ordern(usesavelwf,ioptlwf,na_s,no_s,lasto,isa,qa,rcoor,
     .                  rmaxh,ucell,xa,iscf,istp,ncgmax,etol,eta,
     .                  qtot,no_l,maxnh,numh,listhptr,listh,H,S,
     .                  chebef,noeta,rcoorcp,beta,pmax,Dscf,Escf,
     .                  Ecorrec)
          else
            call die('siesta: ERROR: wrong solution method')
          endif
C ..................

C Harris-functional energy ............................................
          Eharrs = 0.0_dp
          do ispin = 1,nspin
C const factor takes into account that there are two nondiagonal
C elements in non-collinear spin density matrix, stored as
C ispin=1 => D11; ispin=2 => D22, ispin=3 => Real(D12);
C ispin=4 => Imag(D12)
            const = 1._dp
            if (ispin .gt. 2) const = 2._dp
            do io = 1,nh
              Eharrs = Eharrs + H(io,ispin) * const * 
     .                     ( Dscf(io,ispin) - Dold(io,ispin) )
            enddo
          enddo
#ifdef MPI
C Global reduction of Eharrs
          call MPI_AllReduce(Eharrs,buffer1,1,DAT_double,
     .      MPI_sum,MPI_Comm_World,MPIerror)
          Eharrs = buffer1
#endif
C ..................

C Mix input and output energy-density and density matrices ............
C Following line for using and saving the density matrix without mix ..
          if (wmix.ne.0._dp) then
C Pulay mixing
            mmix  = mix
            iiscf = iscf
            if (maxsav .le. 0) then
              iiscf = 1
              if (iscf .ne. 1) mmix = .true.
            endif
            call pulayx( pulfile, iiscf, mmix, no_l, no_s, maxnh,
     .                   numh, listhptr, nspin, maxsav, wmix, nkick, 
     .                   wmixkick, auxpul(1,1), auxpul(1,2), nauxpul,
     .                   Dscf, Dold, dDmax)
          endif

C Ensure that dDmax is the same on all nodes for convergence test/output
#ifdef MPI
          call MPI_AllReduce(dDmax,buffer1,1,DAT_double,
     .      MPI_max,MPI_Comm_World,MPIerror)
          dDmax = buffer1
#endif
C ...................

C Save density matrix on disk, after mixing ...........................
          if (writedm) then
            if ((idyn .eq. 6) .or. (idyn .eq. 7)) then
              if (istp .eq.1)
     .        call iodm( 'write', maxnh, no_l, nspin,
     .                   numh, listhptr, listh, Dscf, found )
            else
              call iodm( 'write', maxnh, no_l, nspin,
     .                   numh, listhptr, listh, Dscf, found )
            endif
          endif
C ...................

C Print energies ......................................................
          DEna = Enascf - Enaatm
          Etot = E0 + DEna + DUscf + DUext + Exc + Ecorrec + Emad
          Eharrs = Etot + Eharrs
          FreeE  = Etot - Temp * Entrop
          if (IOnode) then
            if (istp.eq.1 .and. first) write(6,'(/,a,/,(a,f15.6))')
     .        'siesta: Program''s energy decomposition (eV):',
     .        'siesta: Eions   =', Eions/eV,
     .        'siesta: Ena     =', Ena/eV,
     .        'siesta: Ekin    =', Ekin/eV,
     .        'siesta: Enl     =', Enl/eV,
     .        'siesta: DEna    =', DEna/eV,
     .        'siesta: DUscf   =', DUscf/eV,
     .        'siesta: DUext   =', DUext/eV,
     .        'siesta: Exc     =', Exc/eV,
     .        'siesta: eta*DQ  =', Ecorrec/eV,
     .        'siesta: Emadel  =', Emad/eV,
     .        'siesta: Eharris =', Eharrs/eV,
     .        'siesta: Etot    =', Etot/eV,
     .        'siesta: FreeEng =', FreeE/eV
C ...................

C Print total energy and density matrix error .........................
            if (isolve .eq. 0) then
              if (fixspin) then
                if (iscf .eq. 1)  write(6,'(/,a12,3a14,a8,a7,a11)')
     .            'siesta: iscf', '   Eharris(eV)', 
     .            '      E_KS(eV)', '   FreeEng(eV)', 
     .            '   dDmax', '  Ef_up', '  Ef_dn(eV)'
                write(6,'(a8,i4,3f14.4,f8.4,2f9.4)')
     .            'siesta: ',iscf, Eharrs/eV, Etot/eV, FreeE/eV, dDmax, 
     .                   (Efs(i)/eV,i=1,2)

              else
                if (iscf .eq. 1) write(6,'(/,a12,3a14,2a8)')
     .            'siesta: iscf', '   Eharris(eV)', 
     .            '      E_KS(eV)', '   FreeEng(eV)', 
     .            '   dDmax', '  Ef(eV)'
                write(6,'(a8,i4,3f14.4,2f8.4)')
     .            'siesta: ',iscf, Eharrs/eV, Etot/eV, FreeE/eV, 
     .            dDmax, Ef/eV
              endif
            endif
            if (isolve .eq. 1) then
              write(6,'(/,a15,i4)') 'siesta: iscf = ',iscf
              write(6,'(a14,f15.4,a13,f15.4,a10,f7.4/)') 
     .            'Eharris(eV) = ',Eharrs/eV,
     .            '  E_KS(eV) = ',Etot/eV,'  dDmax = ',dDmax
            endif
          endif
        endif
C ...................

C If last iteration, exit SCF loop ....................................
        if (last) then
          Dscf(1:nh,1:nspin) = Dold(1:nh,1:nspin)
          Escf(1:nh,1:nspin) = Eold(1:nh,1:nspin)
          if (dumpcharge) then
             if (IOnode) write(6,*) 'Charge info in .RHO file'
          endif
          call timer( 'IterSCF', 2 )
          goto 50
        endif
C ...................

C If converged, make last iteration to find forces ....................
        if (dDmax.lt.dDtol) last = .true.
C ...................

        call timer( 'IterSCF', 2 )
        if (istep.eq.inicoor .and. first) call timer( 'IterSCF', 3 )
        first = .false.
      enddo
   50 continue
C End of SCF iteration_________________________________________________

C Write final Kohn-Sham Energy ........................................
      if (IOnode) then
        write(6,"(/a,f14.4)") 'siesta: E_KS(eV) = ', Etot/eV
      endif

#ifdef MPI
C Global reduction of forces and stresses
      allocate(fatmp(3,na_u))
      call memory('A','D',3*na_u,'siesta')
      call MPI_AllReduce(stressl(1,1),stresstmp(1,1),9,
     .  DAT_double,MPI_sum,MPI_Comm_World,MPIerror)
      call MPI_AllReduce(fal(1,1),fatmp(1,1),3*na_u,
     .  DAT_double,MPI_sum,MPI_Comm_World,MPIerror)
      stress(1:3,1:3) = stress(1:3,1:3) + stresstmp(1:3,1:3)
      fa(1:3,1:na_u) = fa(1:3,1:na_u) + fatmp(1:3,1:na_u)
      call memory('D','D',size(fatmp),'siesta')
      deallocate(fatmp)
#else
      stress(1:3,1:3) = stress(1:3,1:3) + stressl(1:3,1:3)
      fa(1:3,1:na_u) = fa(1:3,1:na_u) + fal(1:3,1:na_u)
#endif

C Impose constraints to atomic movements by changing forces ...........
      call fixed(ucell,stress,na_u,isa, amass, xa, fa, cstress, cfa )
C ...................

C Write atomic forces .................................................
      fmax = 0.0_dp
      cfmax = 0.0_dp
      fres = 0.0_dp
      do ix = 1,3
        ftot(ix) = 0.0_dp
        do ia = 1,na_u
          ftem = fa(ix,ia)
          cftem = cfa(ix,ia)
          ftot(ix) = ftot(ix) + ftem
          fres = fres + ftem*ftem
          fmax = max( fmax, dabs(ftem) )
          cfmax = max( cfmax, dabs(cftem) )
        enddo
      enddo
      fres = dsqrt( fres / (3.0_dp*na_u) )
      if (IOnode) then
        write(6,'(/,a)') 'siesta: Atomic forces (eV/Ang):'
        if ( writef ) then
        write(6,'(i4,3f12.6)') (ia,(fa(ix,ia)*Ang/eV,ix=1,3),ia=1,na_u)
        else
          call iofa( na_u, fa )
        endif
        write(6,'(40(1h-),/,a4,3f12.6)') 'Tot',(ftot(ix)*Ang/eV,ix=1,3)
        write(6,'(40(1h-),/,a4, f12.6)') 'Max',fmax*Ang/eV
        write(6,'(a4,f12.6,a)')'Res',fres*Ang/eV,
     .                   '    sqrt( Sum f_i^2 / 3N )'
        write(6,'(40(1h-),/,a4, f12.6,a)') 'Max',cfmax*Ang/eV, 
     .                   '    constrained'

C Write Force Constant matrix if FC calculation ...

        if (idyn .eq. 6) then
          call ofc(fa,dx,na_u)
        endif
        if (idyn .eq. 7) then
          call phonon_write_forces(fa,na_u,ns,ucell,istep)
        endif

C Write stress tensor for any variable cell dynamics ..................
        if ( varcel ) then
           write(6,'(/,a,3(/,a,3f12.6))')
     .       'siesta: Stress tensor (eV/Ang**3):',
     .       ('     ',(stress(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)

C       Pressure (only for the solid)
          Psol = - ((stress(1,1) + stress(2,2) + stress(3,3))/3.0_dp)
          write(6,'(/,a,f20.8,a)')
     .      'siesta: Pressure:', Psol/kBar, '  kBar'

        endif
      endif
C ...................

C Mulliken population analysis .......................................
      call mulliken( mullipop, nspin, na_u, no_l, maxnh,
     .               numh, listhptr, listh, S, Dscf, isa, 
     .               lasto, iaorb, iphorb )

C ...................

C Move atoms ..........................................................
      if (idyn .eq. 0) then 
        xalast(1:3,1:na_s)=xa(1:3,1:na_s)
        if (nmove .ne. 0) then
          call cgvc( na_u, xa, cfa, ucell, cstress, volume, dxmax,
     .               tp, ftol, strtol, varcel, relaxd, usesavecg )
!         Propagate the new structure to the virtual supercell
          call superx( ucell, nsc, na_u, na_s, xa, scell )
          if (relaxd) then
            xalast(1:3,1:na_s)=xa(1:3,1:na_s)
C Exit coordinate relaxation loop
            call timer( 'IterMD', 2 )
            goto 60
          endif
        endif
      endif

! Save the last coordinates for which the density matrix has been calculated
      if (istep.eq.fincoor) then
        xalast(1:3,1:na_s)=xa(1:3,1:na_s)
      endif

      Ekinion  = 0.0_dp
      vn       = 0.0_dp
      vpr      = 0.0_dp
      kn       = 0.0_dp
      kpr      = 0.0_dp

      iunit = 2

      if (idyn .eq. 1) then
        call verlet2(istp, iunit, iquench, na_u, cfa, dt, 
     .               amass, va, xa, Ekinion, tempion)
!         Propagate the new structure to the virtual supercell
          call superx( ucell, nsc, na_u, na_s, xa, scell )
      
      elseif (idyn .eq. 2) then
        call nose(istp, iunit, na_u, cfa, tt, dt, amass, mn, 
     .            va, xa, Ekinion, kn, vn, tempion)
!         Propagate the new structure to the virtual supercell
          call superx( ucell, nsc, na_u, na_s, xa, scell )

      elseif (idyn .eq. 3) then
        call pr(istp, iunit, iquench, na_u, cfa, cstress, tp, dt, amass,
     .          mpr, va, xa, vcell, ucell, Ekinion, kpr, vpr, tempion,
     .          Pint)
!         Propagate the new structure to the virtual supercell
          call superx( ucell, nsc, na_u, na_s, xa, scell )
        if (IOnode) then
          write(6,'(/,a,f12.3,a)')
     .      'siesta: E_kin PR =', kpr/Kelvin, ' K'
        endif
      
      elseif (idyn .eq. 4) then
        call npr(istp, iunit, na_u, cfa, cstress, tp, tt, dt, amass, mn,
     .           mpr, va, xa, vcell, ucell, Ekinion, kn, kpr, vn, vpr,
     .           tempion, Pint)
!         Propagate the new structure to the virtual supercell
          call superx( ucell, nsc, na_u, na_s, xa, scell )

      elseif (idyn .eq. 5) then
        call anneal(istp, iunit, ianneal, taurelax, bulkm,
     .              na_u, cfa, cstress, tp, tt, dt, amass,
     .              va, xa, ucell, Ekinion, tempion, Pint)
!         Propagate the new structure to the virtual supercell
          call superx( ucell, nsc, na_u, na_s, xa, scell )
      endif

      if (IOnode) then
        if (idyn .gt. 0 .and. idyn .lt. 6) then
          write(6,'(/,a,f12.3,a)')
     .      'siesta: Temp_ion =', tempion, ' K'
        endif
      endif
C ...................

C Save last atomic positions and velocities ..........................
      call ioxv( 'write', ucell, vcell, na_u, isa, iza, xa, va, foundxv)
C ...................

      getot = Etot + Ekinion + kn + kpr + vn + vpr

C Add kinetic term to stress tensor ..................................
*     do ia = 1,na_u
*       do jx = 1,3
*         do ix = 1,3
*           stress(ix,jx) = stress(ix,jx) -
*    .            amu * amass(ia) * va(ix,ia) * va(jx,ia) / volume
*         enddo
*       enddo
*     enddo
C ...................

c restore original coordinates after FC displacements ...
      if (idyn .eq. 6 .and. istep .ne. 0) then
        xa(ixdispl,iadispl)=xa(ixdispl,iadispl)-dx
      endif
      if (idyn .eq. 7) then
           call phonon_restore_coords(istep,xa,ucell)
      endif
c ...

C Save atomic positions and velocities accumulatively ................
      if (writmd.and.IOnode) call iomd( na_u, isa, iza, 
     .         xa, va, ucell, vcell, varcel, istep, inicoor, 
     .         fincoor, tempion, Etot, getot)

C Accumulate coor in Xmol file for animation .........................
      lastst = fincoor .le. istep
      if (writpx.and.IOnode) 
     .  call pixmol(iza, xa, na_u, slabel, lastst)
C ...................

C Output memory use at the end of this geometry step
      call printmemory( 6, 0 )

      call timer( 'IterMD', 2 )
      enddo
   60 continue
C End of coordinate-relaxation loop ==================================

      if (IOnode) then
C Print atomic coordinates (and also unit cell for ParrRah.)
        if (nmove .ne. 0) then
          if (relaxd) 
     .      call outcoor(ucell, xa, isa, na_u,'Relaxed', .true. )
          if (.not.relaxd) 
     .      call outcoor(ucell, xa, isa, na_u,
     .                 'Final (unrelaxed)', .true. )
          if ( varcel ) call outcell(ucell)
        endif


C Print coordinates in xmol format in a separate file

        if (fdf_boolean('WriteCoorXmol',.false.)) 
     .     call coxmol(iza, xa, na_u, slabel)

C Print coordinates in cerius format in a separate file

        if (fdf_boolean('WriteCoorCerius',.false.))
     .     call coceri(iza, xa, ucell, na_u, sname, slabel)
      endif


C Find and print bands
      if (nbk.gt.0) then
        call bands( no_s, nspin, nspin, no_u, no_l, maxnh, maxbk,
     .              numh, listhptr, listh, H, S, Ef, xijo, indxuo,
     .              .true., nbk, bk, ebk, no_u )
        if (IOnode) then
          if ( writbk ) then
            write(6,'(/,a,/,a4,a12)')
     .       'siesta: Band k vectors (Bohr**-1):', 'ik', 'k'
            do ik = 1,nbk
              write(6,'(i4,3f12.6)') ik, (bk(ix,ik),ix=1,3)
            enddo
          endif
        
          if ( writb ) then
            write(6,'(/,a,/,a4,a3,a7)')
     .       'siesta: Band energies (eV):', 'ik', 'is', 'eps'
            do ispin = 1,min(nspin,2)
              do ik = 1,nbk
                write(6,'(i4,i3,10f7.2)')
     .            ik, ispin, (ebk(io,ispin,ik)/eV,io=1,min(10,no_u))
                if (no_u.gt.10) write(6,'(7x,10f7.2)')
     .              (ebk(io,ispin,ik)/eV,io=11,no_u)
              enddo
            enddo
          endif
        endif
      endif

C Print eigenvalues
      if (IOnode) then
        if (isolve.eq.0 .and. no_l.lt.1000) then
          if (nspin .le. 2) then
            write(6,'(/,a,/,a4,a3,a7)')
     .       'siesta: Eigenvalues (eV):', 'ik', 'is', 'eps'
            do ik = 1,nkpnt
              do ispin = 1,nspin
                write(6,'(i4,i3,10f7.2)')
     .            ik, ispin, (eo(io,ispin,ik)/eV,io=1,min(10,no_u))
                if (no_u.gt.10) write(6,'(7x,10f7.2)')
     .            (eo(io,ispin,ik)/eV,io=11,no_u)
              enddo
            enddo
          else
            write(6,'(/,a)') 'siesta: Eigenvalues (eV):'
            do ik = 1,nkpnt
              write(6,'(a,i6)') 'ik =', ik
              write(6,'(10f7.2)')
     .          ((eo(io,ispin,ik)/eV,io=1,no_u),ispin=1,2)
            enddo
          endif
          write(6,'(a,f15.6,a)') 'siesta: Fermi energy =', ef/eV, ' eV'
        endif
      endif

      if (isolve.eq.0.and.IOnode) 
     .     call ioeig(eo, ef, no_l, nspin, nkpnt, no_l, nspin, maxk)
!--
C Compute the projected density of states
!
      if (IOnode) then
         do_pdos = fdf_block('ProjectedDensityOfStates',iu)
         if (isolve .ne. 0) then
            write(6,*)
     .           'siesta: ERROR: PDOS implemented only with diagon'
            do_pdos = .false.
         endif
      endif
#ifdef MPI
      call MPI_Bcast(do_pdos,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif

      if (do_pdos) then

c     Find the desired energy range
         if (IOnode) then
            read(iu,'(a)') line
            p=>digest(line)
            if (nvalues(p).lt.3 .or. nnames(p).ne.1)
     $            call die("Wrong format in PDOS block")
            factor = fdf_convfac( names(p,1), 'Ry' )
            e1 = values(p,1) * factor
            e2 = values(p,2) * factor
            sigma = values(p,3) * factor
            nhist = integers(p,1)
            write(6,*) '>>'
            write(6,*) 'PDOS info: '
            write(6,*) 'e1, e2, sigma, nhist: ', e1, e2, sigma, nhist
            write(6,*) '>>'
!
         endif  ! IOnode
#ifdef MPI
        call MPI_Bcast(e1,1,DAT_double,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(e2,1,DAT_double,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(sigma,1,DAT_double,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(nhist,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif
        if (IOnode) write(6,*) 'PDOS not implemented yet'

      endif                     ! PDOS calc (do_pdos)

C Print program's energy decomposition
      if (IOnode) then
        write(6,'(/,a,/,(a,f15.6))')
     .     'siesta: Program''s energy decomposition (eV):',
     .     'siesta:-Eions   =', (-Eions)/eV,
     .     'siesta: Ena     =', Ena/eV,
     .     'siesta: Ekin    =', Ekin/eV,
     .     'siesta: Enl     =', Enl/eV,
     .     'siesta: DEna    =', DEna/eV,
     .     'siesta: DUscf   =', DUscf/eV,
     .     'siesta: DUext   =', DUext/eV,
     .     'siesta: Exc     =', Exc/eV,
     .     'siesta: eta*DQ  =', Ecorrec/eV,
     .     'siesta: Emadel  =', Emad/eV,
     .     'siesta: Ekinion =', Ekinion/eV,
     .     'siesta: Eharris =', (Eharrs+Ekinion)/eV,
     .     'siesta: Etot    =', (Etot+Ekinion)/eV,
     .     'siesta: FreeEng =', (FreeE+Ekinion)/eV

C Print standard energy decomposition
        write(6,'(/,a)') 'siesta: Final energy (eV):'
        write(6,'(a,a15,f15.6)')
     .    'siesta: ',      'Kinetic =', Ekin/eV,
     .    'siesta: ',      'Hartree =', Uscf/eV,
     .    'siesta: ',   'Ext. field =', DUext/eV,
     .    'siesta: ',  'Exch.-corr. =', Exc/eV,
     .    'siesta: ', 'Ion-electron =', (Enascf+Enl+DUscf-Uscf-Uatm)/eV,
     .    'siesta: ',      'Ion-ion =', (Ena+Uatm-Enaatm-Eions)/eV,
     .    'siesta: ',      'Ekinion =', Ekinion/eV,
     .    'siesta: ',        'Total =', (Etot+Ekinion)/eV
      endif

C Print atomic forces
      fmax = 0.0_dp
      do ix = 1,3
        ftot(ix) = 0.0_dp
        do ia = 1,na_u
          fmax = max( fmax, dabs(fa(ix,ia)) )
          ftot(ix) = ftot(ix) + fa(ix,ia)
        enddo
      enddo
      if (IOnode) then
        if (fmax .gt. ftol) then
          write(6,'(/,a)') 'siesta: Atomic forces (eV/Ang):'
          write(6,'(a,i4,3f12.6)')
     .      ('siesta: ', ia,(fa(ix,ia)*Ang/eV,ix=1,3),ia=1,na_u)
          write(6,'(a,40(1h-),/,a,a4,3f12.6)')
     .      'siesta: ','siesta: ','Tot',(ftot(ix)*Ang/eV,ix=1,3)
        endif
      endif

C Print constrained atomic forces
      same = .true.
      do ia = 1,na_u
        do ix = 1,3
          if (cfa(ix,ia) .ne. fa(ix,ia)) same = .false.
        enddo
      enddo
      if (.not.same) then
        fmax = 0.0_dp
        do ix = 1,3
          ftot(ix) = 0.0_dp
          do ia = 1,na_u
            fmax = max( fmax, dabs(cfa(ix,ia)) )
            ftot(ix) = ftot(ix) + cfa(ix,ia)
          enddo
        enddo
        if (IOnode) then
          if (fmax .gt. ftol) then
            write(6,'(/,a)') 'siesta: Constrained forces (eV/Ang):'
            write(6,'(a,i4,3f12.6)')
     .        ('siesta: ', ia,(cfa(ix,ia)*Ang/eV,ix=1,3),ia=1,na_u)
            write(6,'(a,40(1h-),/,a,a4,3f12.6)')
     .        'siesta: ','siesta: ','Tot',(ftot(ix)*Ang/eV,ix=1,3)
          endif
        endif
      endif

C Print stress tensor
      if (IOnode) then
        write(6,'(/,a,3(/,a,3f12.6))')
     .    'siesta: Stress tensor (eV/Ang**3):',
     .    ('siesta: ',(stress(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)
      endif

C Print constrained stress tensor
      same = .true.
      do ix = 1,3
        do jx = 1,3
          if (cstress(jx,ix) .ne. stress(jx,ix)) same = .false.
        enddo
      enddo
      if (IOnode) then
        if (.not.same)
     .    write(6,'(/,a,3(/,a,3f12.6))')
     .      'siesta: Constrained stress tensor (eV/Ang**3):',
     .      ('siesta: ',(cstress(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)
      endif

C Find pressure
      virial = 0.0_dp
      do ix = 1,3
        fmean = 0.0_dp
        do ia = 1,na_u
          fmean = fmean + fa(ix,ia) / na_u
        enddo
        do ia = 1,na_u
          virial = virial + xalast(ix,ia) * (fa(ix,ia) - fmean)
        enddo
      enddo
      Psol = - (( stress(1,1) + stress(2,2) + stress(3,3) )/3.0_dp)
      Pmol = Psol - virial / volume / 3.0_dp
      if (IOnode) then
        write(6,'(/,a,f18.6,a)')
     .    'siesta: Cell volume =', volume/Ang**3, ' Ang**3'
        write(6,'(/,a,/,a,2a20,a,3(/,a,2f20.8,a))')
     .    'siesta: Pressure:',
     .    'siesta: ','Solid',        'Molecule',      '  Units',
     .    'siesta: ', Psol,           Pmol,           '  Ry/Bohr**3',
     .    'siesta: ', Psol*Ang**3/eV, Pmol*Ang**3/eV, '  eV/Ang**3',
     .    'siesta: ', Psol/kBar,      Pmol/kBar,      '  kBar'
      endif

c Print spin polarization
      if (nspin .ge. 2) then
        do ispin = 1,nspin
          qspin(ispin) = 0.0_dp
          do io = 1,no_l
            do j = 1,numh(io)
              ind = listhptr(io)+j
              jo = listh(ind)
              qspin(ispin) = qspin(ispin) + Dscf(ind,ispin) * S(ind)
            enddo
          enddo
        enddo
        if (nspin .eq. 2) then
          if (IOnode) then
            write(6,'(/,a,f12.6)')
     .       'siesta: Total spin polarization (Qup-Qdown) =', 
     .       qspin(1) - qspin(2)
          endif
        elseif (nspin .eq. 4) then
          call spnvec( nspin, qspin, qaux, stot, svec )
          if (IOnode) then
            write(6,'(/,a,f12.6)')
     .       'siesta: Total spin polarization (Qup-Qdown) =', stot
            write(6,'(a,3f12.6)') 'siesta: Spin vector =', svec
          endif
        endif
      endif

c Print electric dipole
      if (IOnode) then
        if (shape .ne. 'bulk') then
          write(6,'(/,a,3f12.6)')
     .      'siesta: Electric dipole (a.u.)  =', dipol
          write(6,'(a,3f12.6)')
     .      'siesta: Electric dipole (Debye) =', 
     .      (dipol(ix)/Debye,ix=1,3)
        endif
      endif

C Calculation of the bulk polarization using the Berry phase
C formulas by King-Smith and Vanderbilt
C Attention H0 is used as an auxiliar array (It does not contain
C any useful information any more ) !!!!
      if (nkpol.gt.0) then
        call KSV_pol(na_u, na_s, xalast, rmaxo, scell, ucell,
     .               no_u, no_l, no_s, nspin, qspin, maxna, 
     .               maxnh, nkpol, numh, listhptr, listh, 
     .               H, S, H0, xijo, indxuo, isa, iphorb, 
     .               iaorb, lasto, jna, xij, r2ij,shape,
     .               nkpol,kpol,wgthpol, polR, polxyz) 
      endif
      if (nkpol.gt.0.and.IOnode) then
        do ispin=1,nspin
          if (nspin.gt.1) then
            if (ispin.eq.1) write(6,'(/,a)')
     .       'siesta: Macroscopic polarization for spin Up:'
            if (ispin.eq.2) write(6,'(/,a)')
     .       'siesta: Macroscopic polarization for spin Down:'
          endif
          write(6,'(/,a)')
     .     'siesta: Macroscopic polarization per unit cell (a.u.):'
          write(6,'(a,3f12.6)')
     .     'siesta: Along the lattice vectors  ',
     .       (polR(ix,ispin),ix=1,3)
          write(6,'(a,3f12.6)')
     .     'siesta: Along cartesian directions ',
     .      (polxyz(ix,ispin),ix=1,3)
          write(6,'(/,a)')
     .     'siesta: Macroscopic polarization per unit cell (Debye):'
          write(6,'(a,3f12.6)')
     .     'siesta: Along the lattice vectors  ',
     .       (polR(ix,ispin)/Debye,ix=1,3)
          write(6,'(a,3f12.6)')
     .     'siesta: Along cartesian directions ',
     .      (polxyz(ix,ispin)/Debye,ix=1,3) 
        enddo 
        if (nspin.gt.1) then 
C Modified so that compiler is happy when nspin = 1 and bounds checking
C is turned on. JDG
          write(6,'(/,a,/a,3f12.6)')
     .      'siesta: Sum along cartesian directions (a.u.): ',
     .      'siesta: ',(polxyz(ix,1)+polxyz(ix,min(nspin,2)),ix=1,3)
          write(6,'(/,a,/a,3f12.6)')
     .      'siesta: Sum along cartesian directions (Debye): ',
     .      'siesta: ',((polxyz(ix,1)+polxyz(ix,min(nspin,2)))/Debye,
     .      ix=1,3)
        endif

      endif
c...................................

c Save electron density and potential
      if (IOnode) then
        savrho = fdf_boolean( 'SaveRho',    dumpcharge .or. .false.)
        savdrh = fdf_boolean( 'SaveDeltaRho',               .false.)
        savevh = fdf_boolean( 'SaveElectrostaticPotential', .false.)
        savevt = fdf_boolean( 'SaveTotalPotential',         .false.)
      endif
#ifdef MPI
      call MPI_Bcast(savrho,1,MPI_logical,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(savdrh,1,MPI_logical,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(savevh,1,MPI_logical,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(savevt,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
      if (savrho .or. savdrh .or. savevh .or. savevt) then
        filrho = ' '
        fildrh = ' '
        filevh = ' '
        filevt = ' '
        if (savrho) filrho = paste( slabel, '.RHO' )
        if (savdrh) fildrh = paste( slabel, '.DRHO' )
        if (savevh) filevh = paste( slabel, '.VH'  )
        if (savevt) filevt = paste( slabel, '.VT'  )
        g2max = g2cut
        call dhscf( nspin, no_s, iaorb, iphorb, indxuo, no_l,
     .              no_u, na_u, na_s, isa, xa, indxua, 
     .              ucell, mscell, g2max, ntm,
     .              0, 0, 0, filrho, fildrh, filevh, filevt,
     .              maxnh, numh, listhptr, listh, Dscf, Datm,
     .              maxnh, numh, listhptr, listh, H,
     .              Enaatm, Enascf, Uatm, Uscf, DUscf, DUext, Exc, Dxc,
     .              dipol, fa, stress, fal, stressl )
      endif

c Find local density of states
      if (IOnode) then
        genlogic = fdf_block('LocalDensityOfStates',iu)
      endif
#ifdef MPI
      call MPI_Bcast(genlogic,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
      if ( genlogic ) then

c       Find the desired energy range
        if (IOnode) then
          read(iu,'(a)') line
          p=>digest(line)
          if (.not. match(p,"vvn"))
     $       call die("Wrong format in LocalDensityofStates")
          factor = fdf_convfac( names(p,1), 'Ry' )
          e1 = values(p,1) * factor
          e2 = values(p,2) * factor
          call destroy(p)
        endif
#ifdef MPI
        call MPI_Bcast(e1,1,DAT_double,0,
     .    MPI_Comm_World,MPIerror)
        call MPI_Bcast(e2,1,DAT_double,0,
     .    MPI_Comm_World,MPIerror)
#endif

c       Find the density matrix for states between e1 and e2
        if (isolve .eq. 0) then
          call diagon(no_s, nspin, nspin, no_l, maxnh, maxnh, no_u,
     .                numh, listhptr, listh, numh, listhptr, listh, 
     .                H, S, qtot, fixspin, qs, temp, e1, e2,
     .                gamma, xijo, indxuo, nkpnt, kpoint, kweight,
     .                eo, qo, Dscf, Escf, ef, efs, no_u)
        else
          if (IOnode) then
            write(6,*)
     .       'siesta: ERROR: LDOS implemented only with diagon'
          endif
          goto 70
        endif
     
c       Find the LDOS in the real space mesh
        filrho = paste( slabel, '.LDOS' )
        g2max = g2cut
        call dhscf( nspin, no_s, iaorb, iphorb, indxuo, no_l,
     .              no_u, na_u, na_s, isa, xa, indxua, 
     .              ucell, mscell, g2max, ntm,
     .              0, 0, 0, filrho, ' ', ' ', ' ',
     .              maxnh, numh, listhptr, listh, Dscf, Datm,
     .              maxnh, numh, listhptr, listh, H,
     .              Enaatm, Enascf, Uatm, Uscf, DUscf, DUext, Exc, Dxc,
     .              dipol, fa, stress, fal, stressl )
      endif
   70 continue

C Output memory use up to the end of the program
#ifdef MPI
      call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif
      call printmemory( 6, 1 )
#ifdef MPI
      call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif

c Print allocation report
      call alloc_report( printNow=.true. )

c Stop time counter
      call timer( 'siesta', 2 )
      call timer( 'all', 3 )

C Print final date and time
      if (IOnode) then
        call timestamp('End of run')
      endif

C Finalize MPI
#ifdef MPI
      call MPI_Finalize( MPIerror )
#endif

      end









