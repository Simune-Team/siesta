      Program SIESTA

C ***
C SIESTA Density Functional LCAO program package.
C Copyright by Fundacion General Universidad Autonoma de Madrid:
C E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal 
C and J.M.Soler, 1996-2003
C ***
C Copy or disemination of all or part of this package is not permitted 
C without prior and explicit authorization by the authors.
C Send comments/suggestions/bug-reports to siesta@uam.es
C ***

C  Modules
      use precision, only: dp
      use parallel, only: Node, Nodes, IOnode
      use atmfuncs, only: rcut, uion, izofis
      use atomlist, only: xa, xalast, indxuo, na_u
      use atomlist, only: isa, elem, indxua, iza, lastkb, lasto
      use atomlist, only: qa, amass, qtot, no_u, rmaxo, rmaxkb
      use atomlist, only: rmaxv, na_s, iphorb, iphkb, rco, rckb
      use atomlist, only: datm, no_s, iaorb
      use atomlist, only: superx, superc, initatomlists
!
      use band, only: initbands, bands
      use fdf
      use sys, only: die, bye
      use parse
      use ionew, only: io_setup
      use periodic_table, only : symbol
      use xcmod, only: setXC
#ifdef MPI
      use mpi_siesta, only: mpi_comm_world, mpi_sum, mpi_max
      use mpi_siesta, only: mpi_integer, mpi_logical
#ifdef NODAT
      use mpi_siesta, only: mpi_double_precision
#else
      use mpi_siesta, only: DAT_double
#endif
#endif
      use alloc, only: re_alloc, alloc_report

      use phonon, only: phonon_num_disps, phonon_setup
      use phonon, only: phonon_write_forces, phonon_restore_coords
      use phonon, only: phonon_set_coords

      use densematrix

      use m_ordern,     only: ordern
      use m_hsparse,    only: hsparse

      use parallelsubs, only: getnodeorbs
      use writewave,    only: initwave, wwave
      use iopipes,      only: forcestopipe, coordsfrompipe

      use siesta_cml

      use time, only : timestamp
      use units

      implicit none

      integer
     .  fincoor,
     .  i, ia, ia1, ia2, iadispl, ianneal,
     .  idyn, ifa, ifinal, ihmat, ihuge, iiscf, ik, 
     .  ind, inicoor, io, ioptlwf, iord,
     .  iquench, is, isel, iscf, 
     .  isolve, ispin, istp, istpsave, istart, istep, istr,
     .  iu, iunit,  iv, ix, ixdispl, 
     .  j, ja, jamin, jo, jx, kscell(3,3),
     .  level, maxbk, maxk, maxnh, 
     .  maxna, maxsav, maxwk, mscell(3,3), mullipop, 
     .  nkick, nnamax, nauxpul, nbcell, nbk, ncgmax, nh,
     .  nkpnt, nmove, nnia, ns, nsc(3), nscold(3), nscf, nspin, 
     .  ntm(3), ntcon, no_l, nxij, nwk, pmax, nkpol, nhist, io_sig,
     .  neigwanted, neigmin

#ifdef MPI
      integer
     .  MPIerror, ntmp
#endif

      integer, dimension(:), allocatable ::
     .  jna

      integer, pointer, save ::
     .  listh(:), listhold(:), listhptr(:), listhptrold(:),
     .  numh(:), numhold(:)

      real(dp)
     .  bcell(3,3), beta, bulkm, 
     .  charnet, cfmax, cftem, const, cstress(3,3),
     .  dipol(3), dDmax, dDtol, DEharr, DEna,
     .  dt, DUext, DUscf, Dxc, dx, dxmax, e1, e2, sigma,
     .  Ecorrec, ef, Eharrs, Eharrs1, Eions, Ekin, Ekinion, Emad, Ena, 
     .  Enaatm, Enascf, Enl, Entrop, Entropy,
     .  eta(2), etol, Etot, Exc, E0,
     .  factor,
     .  fmax, fmean, FreeE, fres, ftem, ftol, ftot(3),
     .  g2cut, g2max, getot,
     .  kcutof, kdispl(3),
     .  kn, kpr, 
     .  mn, mpr, occtol, Pint, Pmol, Psol
      real(dp)
     .  qaux, qspin(4), qsol,
     .  rc, rcoor, rcoorcp, rijmin, rmax, rmaxh, rmin, r2min,
     .  scell(3,3), stot, stress(3,3), strtol, svec(3),
     .  taurelax, temp, tempinit, tempion, tiny, tp, ts, tt,
     .  tstress(3,3), Uatm, ucell(3,3), Uscf,
     .  vcell(3,3), virial, vn,
     .  volcel, volume, vpr, wmix, wmixkick, 
     .  xmax, xmin, stressl(3,3), veclen

      real(dp), dimension(:), allocatable ::
     .  efs, kweight, qs, r2ij

      real(dp), pointer, save :: 
     .  H0(:), S(:), wgthpol(:)

      real(dp), dimension(:,:), allocatable :: 
     .  auxpul, cfa, fa, fal, kpoint, 
     .  polR, polxyz, va, xij, xijo, wfk

      real(dp), pointer, save ::
     .  Dold(:,:), Dscf(:,:), Dscfsave(:,:), Eold(:,:), Escf(:,:),
     .  bk(:,:), H(:,:), kpol(:,:)

#ifdef MPI
      real(dp)
     .  buffer1, stresstmp(3,3), qtmp(4)

      real(dp), dimension(:,:), allocatable :: 
     .  fatmp
#endif

      real(dp), dimension(:,:,:), allocatable :: 
     .  ebk, eo, qo

      logical
     .  auxchanged, chebef, dminit, default, dumpcharge, 
     .  first, fixauxcell, fixspin, found, foundxv, gamma,
     .  initdmaux, inspn, itest, last, lastst, mix, mmix, negl, noeta, 
     .  outlng, overflow, overflowed, pulfile, relaxd,
     .  savehs, savevh, savevt, savdrh, savrho,
     .  savepsch, savetoch,
     .  usesavecg, usesavelwf, usesavedm, usesavedmloc, usesavexv, 
     .  writbk, writmd, writpx, writb, writec, writef, 
     .  writek, writic, varcel, genlogic, do_pdos, writedm, atmonly,
     .  harrisfun, signat, muldeb, eggbox_block, final

      logical  :: printeigen = .false.       ! AG: Initialize

#ifdef MPI
      logical
     .  writei
#endif

      character
     .  filevh*25, filevt*25, fildrh*25, filrho*25,
     .  filepsch*25, filetoch*25,
     .  line*150, paste*25,
     .  slabel*20, sname*150, shape*10, message*79

      external
     .  anneal, cgvc, fixed,
     .  dhscf, diagon, dnaefs, extrapol, initatom,
     .  iodm, ioxv, 
     .  kgrid, kgridinit, kinefsm, mulliken, naefs, neighb, nlefsm, 
     .  nose, npr, overfsm, paste, pulayx, 
     .  redata, reinit, shaper, spnvec, 
     .  timer, verlet2, volcel, xijorb, memory,
     .  ioeig, iofa, iokp, iomd, prversion, eggbox

      type(parsed_line), pointer  :: p

      data
     .  e1, e2 / 1._dp, -1._dp /
     .  kcutof, kscell, kdispl / 0._dp, 9*0, 3*0._dp /
     .  nkpnt / 1 /
     .  eggbox_block /.true./
     .  relaxd /.false./
     .  final /.false./
     .  tiny /1.e-15_dp/
     .  ihuge /1073741823/
      ! maxnh is (initially) the minimum size of all the sparse arrays.
      ! It must be 1 (not 0) since we frequently use the F77 idiom of
      ! passing the first element of the array instead of all of it.
     .  maxnh / 1 /
     .  maxna / 200 /
     .  tempion / 0.0_dp /
     .  no_l, nnamax /2*1/
     .  nauxpul, nbk, ns, nspin, nxij /5*1/
     .  nscold / 3*0 /
     .  nsc, mscell / 1,1,1,   1,0,0, 0,1,0, 0,0,1 /
c---------------------------------------------------------------------

C Initialise MPI and set processor number
#ifdef MPI
      call MPI_Init( MPIerror )
      call MPI_Comm_Rank( MPI_Comm_World, Node, MPIerror )
      call MPI_Comm_Size( MPI_Comm_World, Nodes, MPIerror )
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
         call timestamp('Start of run')
      endif
C ..................

C Start time counter ..................................................
      call timer( 'siesta', 0 )
      call timer( 'siesta', 1 )
      call timer( 'Setup', 1 )

C Nullify arrays
      nullify(Haux,Saux,psi)
      nullify(bk,H,H0,S,kpol,wgthpol)
      nullify(Dold,Dscf,Dscfsave,Eold,Escf)
      nullify(listh,listhold,listhptr,listhptrold)
      nullify(numh,numhold)

C Initialize some variables
      DUext = 0.0_dp
      Eharrs = 0.0_dp
      Eharrs1 = 0.0_dp
      Eions = 0.0_dp
      Ekinion = 0.0_dp
      Emad = 0.0_dp
      Entrop = 0.0_dp
      Entropy = 0.0_dp
      FreeE = 0.0_dp

C Initialise read .....................................................
      call reinit(slabel, sname)

      call siesta_cml_init() ! Initialize CML (relies on reinit)

C Set allocation report level .........................................
      if (IOnode) level = fdf_integer( 'alloc_report_level', 0 )
#ifdef MPI
      call MPI_Bcast(level,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif
      call alloc_report( level=level, file=trim(slabel)//'.alloc' )
C ..................

C Initialise exchange-correlation functional information
      call setXC

C Initialize pseudopotentials and atomic orbitals
      if (IOnode) call initatom
      call broadcast_basis

      if (IOnode) atmonly = fdf_defined('Atom-Setup-Only')
#ifdef MPI
      call MPI_Bcast(atmonly,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
      if (atmonly) call bye("End of atom setup")

      if (Node.eq.0) then
        write(6,'(/,a,20("*"),a,28("*"))')
     .    'siesta: ', ' Simulation parameters '
        write(6,'(a)')  'siesta:'
        write(6,'(a)')  'siesta: The following are some of the '//
     .                           'parameters of the simulation.'
        write(6,'(a)')  'siesta: A complete list of the parameters '//
     .                           'used, including defect values,'
        write(6,'(a,a)')'siesta: can be found in file out.fdf'
        write(6,'(a)')  'siesta:'
      endif

C Read simulation sizes ...............................................

!! Read number of atoms and coordinates, and unit cell
      if (IOnode) call coor(na_u,ucell)  ! Sets na_u, xa, and isa
      if (cml_p) then
        ! We need the names of the elements on node 0
        allocate(elem(na_u))
        call memory('A','S',2*na_u,'siesta')
        do i = 1, na_u
          elem(i) = symbol(izofis(isa(i)))
        enddo
      endif

#ifdef MPI
      call MPI_Bcast(na_u,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#ifdef NODAT
      call MPI_Bcast(ucell(1,1),9,MPI_double_precision,0,MPI_Comm_World,
     .  MPIerror)
#else
      call MPI_Bcast(ucell(1,1),9,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif

      if (Node.ne.0) then    ! Allocate isa and xa in rest of nodes
         nullify(isa,xa)
         call re_alloc(isa,1,na_u,name='isa',routine='siesta')
         call re_alloc(xa,1,3,1,na_u,name='xa',routine='siesta')
      endif
      call MPI_Bcast(isa,na_u,MPI_integer,0,MPI_Comm_World,MPIerror)
#ifdef NODAT
      call MPI_Bcast(xa(1,1),3*na_u,MPI_double_precision,0,
     .  MPI_Comm_World,MPIerror)
#else
      call MPI_Bcast(xa(1,1),3*na_u,DAT_double,0,MPI_Comm_World,
     .  MPIerror)
#endif
#endif

C Allocate arrays based on read sizes ................
!
!     Those living in module atomlists
!
      nullify(indxua,iza,lastkb,lasto,qa,amass,xalast)
      call re_alloc(indxua,1,na_u,name='indxua',routine='siesta')
      call re_alloc(iza,1,na_u,name='iza',routine='siesta')
      call re_alloc(lastkb,0,na_u,name='lastkb',routine='siesta')
      call re_alloc(lasto,0,na_u,name='lasto',routine='siesta')
      call re_alloc(qa,1,na_u,name='qa',routine='siesta')
      call re_alloc(xalast,1,3,1,na_u,name='xalast',routine='siesta')
      call re_alloc(amass,1,na_u,name='amass',routine='siesta')

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
     .            pulfile, tempinit, dumpcharge, varcel, harrisfun,
     .            occtol)

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
        fixauxcell= fdf_boolean('FixAuxillaryCell', .false.)
        initdmaux = fdf_boolean('ReInitialiseDM'  , .false.)
        signat    = fdf_boolean('SignatureRecords', .false.)
        muldeb    = fdf_boolean('MullikenInSCF'   , .false.)

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
      call MPI_Bcast(fixauxcell,1,MPI_logical,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(initdmaux,1,MPI_logical,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(muldeb,1,MPI_logical,0,MPI_Comm_World,MPIerror)
      if (.Not.IONode) signat=.false.
#ifdef NODAT
      call MPI_Bcast(rijmin,1,MPI_double_precision,0,MPI_Comm_World,
     .  MPIerror)
#else
      call MPI_Bcast(rijmin,1,DAT_double,0,MPI_Comm_World,MPIerror)
#endif
#endif

C Initialise file SIG for Signature records (whatever they are...) ...
      if (signat) then
        call io_assign(io_sig)
        call sig_setup(io_sig)
      endif
C .....................

C Read cell shape and atomic positions from a former run ..............
      foundxv = .false.
      if (usesavexv) then
        call ioxv('read', ucell, vcell, na_u, isa, iza, xa, va, foundxv)
        if (IOnode.and..not.foundxv) write(6,'(/,a)')
     .    'siesta: WARNING: XV file not found'
      endif
C ..................

C Read cell shape and atomic positions from driver program through pipe
      if (idyn.eq.8) then
        call coordsFromPipe( na_u, xa, ucell )
      end if
C .....................

C Dump initial coordinates to output ..................................
      if ( writic.and.(IOnode) ) then
        write(6,'(/a)') 'siesta: Atomic coordinates (Bohr) and species'
        write(6,"('siesta: ',2x,3f10.5,i3,3x,i6)")
     .           ( (xa(ix,ia), ix=1,3), isa(ia), ia, ia=1, na_u)
      endif
C ..................


C Initialize atom lists 
      call initatomlists
      qtot = qtot - charnet

C Calculate spin populations for fixed spin case...
      if (fixspin) then
        if (nspin .ne. 2) then
          if (Node.eq.0) then
            write(message,'(2a)')
     .       'siesta: ERROR: You can only fix the spin of the system',
     .       'for collinear spin polarized calculations.'
             call die(message)
          endif
        endif
        do i = 1,2
          qs(i) = (qtot + (3-2*i)*ts) / 2.0_dp
        enddo
      else
        qs(1) = qtot
        qs(2:nspin) = 0.0_dp
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

C Find system shape ...................................................
      call shaper( ucell, na_u, isa, xa, shape, nbcell, bcell )
      if (IOnode) then
        write(6,'(/,2a)') 'siesta: System type = ', shape
      endif

C Output of initial system details:
      if (cml_p) then
        call cmlStartModule(xf=mainXML, title='Initial System')
        call cmlAddMolecule(xf=mainXML, natoms=na_u, 
     .       coords=xa/Ang, elements=elem, refs=isa,
     .       style='x3', fmt='(f12.6)')
        call cmlAddLattice(xf=mainXML, cell=ucell, units='Ang', 
     .       dictref='siesta:ucell')
        call cmlAddProperty(xf=mainXML, property=trim(shape), 
     .       dictref='siesta:shape')
        call cmlEndModule(xf=mainXML)
      endif

C Madelung correction for charged systems .............................
      if (charnet .ne. 0.0_dp) then
        call madelung(ucell, shape, charnet, Emad)
      endif

C Parallel initialisation
      call initparallel(no_u,na_u,lasto,xa,ucell,rmaxh,rcoor,isolve)

C Find number of locally stored orbitals and allocated related arrays
      call GetNodeOrbs(no_u,Node,Nodes,no_l)

C Initialise arrays
      call re_alloc(listhptr,1,no_l,name='listhptr',routine='siesta',
     .              copy=.false.)
      call re_alloc(listhptrold,1,no_l,name='listhptrold',
     .              routine='siesta',copy=.false.)
      call re_alloc(numh,1,no_l,name='numh',routine='siesta',
     .              copy=.false.)
      call re_alloc(numhold,1,no_l,name='numhold',routine='siesta',
     .              copy=.false.)
      listhptr(1:no_l) = 0
      listhptrold(:) = 0
      numh(:) = 0
      numhold(:) = 0

C Get number of eigenstates that need to be calculated
      if (IOnode) then
        neigwanted = fdf_integer('NumberOfEigenStates',no_u)
      endif
#ifdef MPI
      call MPI_Bcast(neigwanted,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif

C Check number of eigenstates - cannot be larger than number of
C basis functions or smaller than number of occupied states + 1
C so that the Fermi level can be estimated
      do is = 1,nspin
        neigmin = nint(qs(is)/real(3 - min(nspin,2), kind=dp)) + 1
        neigwanted = max(neigwanted,neigmin)
      enddo
      neigwanted = min(neigwanted,no_u)

C Find k-grid for Brillouin zone integration ..........................
C Call kgrid once to get the dimensions and a second time to actually
C initialise the data - could be made more efficient in future!
      nkpnt = 1
      call kgridinit( ucell, kscell, kdispl, kcutof, nkpnt )
      maxk = nkpnt
      allocate(kpoint(3,maxk))
      call memory('A','D',3*maxk,'siesta')
      allocate(kweight(maxk))
      call memory('A','D',maxk,'siesta')
      kpoint(1:3,1) = 0.0_dp
      kweight(1) = 1.0_dp
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
      call re_alloc(bk,1,3,1,maxbk,name='bk',routine='siesta',
     .              copy=.false.)
C
      call initbands( maxbk, nbk, bk )
C
      if (nbk .gt. maxbk) then
C If there wasn't enough space to store bands on first call correct
C the dimensions and repeat the initialisation
        maxbk = max(nbk,1)
        call re_alloc(bk,1,3,1,maxbk,name='bk',routine='siesta',
     .                copy=.false.)
        nbk = 0
        call initbands( maxbk, nbk, bk )
      endif
      allocate(ebk(no_u,nspin,maxbk))
      call memory('A','D',no_u*nspin*maxbk,'siesta')

C ......................

C Find number of k-points for wavefunction printout ....................
      nwk = 0
      maxwk = 1
      allocate(wfk(3,maxwk))
      call memory('A','D',3*maxwk,'siesta')

      call initwave( maxwk, no_u, nwk, wfk, overflow )

      if (overflow) then
        if (nwk .gt. maxwk) then
C If there wasn't enough space to store bands on first call correct
C the dimensions and repeat the initialisation
          maxwk = max(nwk,1)
          call memory('D','D',size(wfk),'siesta')
          deallocate(wfk)
          allocate(wfk(3,maxwk))
          call memory('A','D',3*maxwk,'siesta')
        endif
        nwk = 0
        call initwave( maxwk, no_u, nwk, wfk, overflow )
        if (overflow) then
          write(message,'(2a)')
     .      'siesta: ERROR: Unsuccessful initialization of',
     .      '        list of wavefunctions to print'
           call die(message)
         endif

      endif
C ......................

C Find the grid for the calculation of the polarization..............
      nkpol = 1
      call re_alloc(kpol,1,3,1,nkpol,name='kpol',routine='siesta',
     .              copy=.false.)
      call re_alloc(wgthpol,1,nkpol,name='wgthpol',routine='siesta',
     .              copy=.false.)

      call KSV_init(ucell, 0, nkpol, kpol, wgthpol)

      call re_alloc(kpol,1,3,1,nkpol,name='kpol',routine='siesta',
     .              shrink=.false.,copy=.false.)
      call re_alloc(wgthpol,1,nkpol,name='wgthpol',routine='siesta',
     .              shrink=.false.,copy=.false.)

C Find if only gamma point is used ....................................
      if (nkpnt.eq.1 .and. abs(kpoint(1,1)).lt.tiny .and.
     .                     abs(kpoint(2,1)).lt.tiny .and.
     .                     abs(kpoint(3,1)).lt.tiny) then
        gamma = .true.
      else
        gamma = .false.
      endif
      if (nbk .gt. 0) gamma = .false.
      if (nwk .gt. 1) gamma = .false.
      if (nwk .eq. 1) then
        if (abs(wfk(1,1)).gt.tiny .and.
     .      abs(wfk(2,1)).gt.tiny .and.
     .      abs(wfk(3,1)).gt.tiny) then
          gamma = .false.
        endif
      endif
      if (nkpol.gt.0) gamma = .false.
C ....................

C Print k-points ......................................................
      if (.not.gamma .and. IOnode) call siesta_write_k_points()

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
        nscold(i) = nsc(i)
      enddo

C Find auxiliary supercell (required only for k sampling) ............
      call superc( ucell, scell, nsc)

C Initialize atomic velocities to zero ................................
      if (.not. foundxv) then
        va(1:3,1:na_u) = 0.0_dp
        vcell(1:3,1:3) = 0.0_dp
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
      else if (idyn .eq. 8) then
        inicoor = 0
        fincoor = ihuge
      else
         call die('siesta: wrong idyn')
      endif

C Build initial velocities according to Maxwell-Bolzmann distribution....
      if (idyn .ne. 0 .and. idyn .ne. 6 .and. (.not. foundxv)) 
     .    call vmb(na_u,tempinit,amass,xa,isa,va)
C ..................

      istp = 0
      call timer( 'Setup', 2 )

C Output memory use before main loop
      call printmemory( 6, 0 )

! Initialization now complete. Flush stdout.
      if (ionode) call pxfflush(6)

C Start loop over coordinate changes 

      do istep = inicoor,fincoor
      call timer( 'IterMD', 1 )
      istp = istp + 1
      if (IOnode) then
        write(6,'(/2a)') 'siesta:                 ',
     .                    '=============================='
        select case (idyn)
        case (0)
          write(6,'(28(" "),a,i6)') 'Begin CG move = ',istep
          if (cml_p) call cmlStartStep(mainXML, type='CG')
        case (1:5)
          write(6,'(28(" "),a,i6)') 'Begin MD step = ',istep
          if (cml_p) call cmlStartStep(mainXML, type='MD')
        case (6)
          write(6,'(28(" "),a,i6)') 'Begin FC step = ',istep
          if (cml_p) call cmlStartStep(mainXML, type='FC')
          if (istep .eq. 0) then
            write(6,'(28(" "),a)') 'Undisplaced coordinates'
          else
            iadispl = (istep-mod(istep-1,6))/6+ia1
            write(6,'(28(" "),a,i6)') 'displace atom   ',
     .        iadispl
            ix = mod(istep-1,6)+1
            ixdispl = (ix - mod(ix-1,2) +1)/2
            write(6,'(28(" "),a,i6)') 'in direction    ',
     .        ixdispl
            dx=-dx
            write(6,'(28(" "),a,f8.4,a)') 'by       ',
     .                      dx, ' Bohr'
C Displace atom by dx...
            xa(ixdispl,iadispl) = xa(ixdispl,iadispl) + dx
          endif
        case (7)
          call phonon_set_coords(istep,xa,ucell)
        end select

        write(6,'(2a)') '                        ',
     .                    '=============================='
      else ! not IOnode
        select case (idyn)
        ! We don't need to do anything for 0<idyn<6
        case(6)
          if (istep .ne. 0) then
            iadispl = (istep-mod(istep-1,6))/6 + ia1
            ix = mod(istep-1,6) + 1
            ixdispl = (ix - mod(ix-1,2) +1)/2
            dx = - dx
C Displace atom by dx...
            xa(ixdispl,iadispl) = xa(ixdispl,iadispl) + dx
          endif
        case(7)
          call phonon_set_coords(istep,xa,ucell)
        end select
      endif

C Get coordinates from driver program trough pipe
      if (idyn.eq.8 .and. istep.ne.inicoor) then
        call coordsFromPipe( na_u, xa, ucell )
      end if

      if (IOnode) then
C Print atomic coordinates ............................................
        call outcoor( ucell, xa, isa, na_u, ' ', writec )
      endif
C ...................


C Actualize things if variable cell ...................................
      auxchanged = .false.
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
          maxk = nkpnt
       else
          call kgrid( ucell, kscell, kdispl,
     .                nkpnt, kpoint, kweight )

       endif
 
C Print k-points 
        if (IOnode) call siesta_write_k_points()
 
C Find required supercell
        auxchanged = .false.
        do i = 1,3
          mscell(1:3,i) = 0.0_dp
          if (gamma) then
            nsc(i) = 1
            mscell(i,i) = 1
          else
            if (fixauxcell) then
              nsc(i) = nscold(i)
              mscell(i,i) = nsc(i)
            else
              veclen = sqrt(ucell(1,i)**2+ucell(2,i)**2+ucell(3,i)**2)
              nsc(i) = 1 + 2 * rmaxh / veclen
              mscell(i,i) = nsc(i)
            endif
          endif
          if (nsc(i).ne.nscold(i)) auxchanged = .true.
          nscold(i) = nsc(i)
        enddo

C Madelung correction for charged systems .............................
        if (charnet .ne. 0.0_dp) then
          call madelung(ucell, shape, charnet, Emad)
        endif

      endif
C End variable cell actualization

C Auxiliary supercell
      call superc(ucell, scell, nsc)

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
C Increase maxna with safety margin when atoms move
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

C List of nonzero Hamiltonian matrix elements .........................
      overflow=.true.
      overflowed=.false.
      do while (overflow)
         nh = maxnh
         call re_alloc(listh,1,maxnh,name='listh',routine='siesta',
     .                 shrink=.false.)
         call hsparse( negl, scell, nsc, na_s, isa, xa, lasto, lastkb, 
     .                 iphorb, iphKB, nh, numh, listhptr, listh )
         if (nh .gt. maxnh) then
      ! Increase maxnh with safety margin for when atoms move
            maxnh = 1.05 * nh + 40
            overflowed=.true.
         else
            overflow=.false.
         endif
      enddo
      ! In first step, allocate anyway (to catch corner case
      ! where one node has nh=0, and doesn't overflow)
      if (istp==1) then
        call re_alloc(Dscf,1,maxnh,1,nspin,name='Dscf',
     .                routine='siesta',copy=.false.)
        call re_alloc(Dscfsave,1,maxnh,1,nspin,name='Dscfsave',
     .                routine='siesta',copy=.false.)
        call re_alloc(listhold,1,maxnh,name='listhold',
     .                routine='siesta',copy=.false.)
        ! Initialise Dscfsave to avoid problems in extrapol
        Dscfsave(1:maxnh,1:nspin) = 0.0_dp
        listhold(1:maxnh) = 0
      elseif (overflowed) then
        ! We need to preserve the contents of these arrays
        ! when reallocating.
        call re_alloc(Dscf,1,maxnh,1,nspin,name='Dscf',
     .                routine='siesta',copy=.true.)
        call re_alloc(Dscfsave,1,maxnh,1,nspin,name='Dscfsave',
     .                routine='siesta',copy=.true.)
        call re_alloc(listhold,1,maxnh,name='listhold',
     .                routine='siesta',copy=.true.)
      endif
      if (istp==1.or.overflowed) then
        call re_alloc(Dold,1,maxnh,1,nspin,name='Dold',
     .                routine='siesta',copy=.false.)
        call re_alloc(Eold,1,maxnh,1,nspin,name='Eold',
     .                routine='siesta',copy=.false.)
        call re_alloc(Escf,1,maxnh,1,nspin,name='Escf',
     .                routine='siesta',copy=.false.)
      endif

C Allocate/reallocate storage associated with Hamiltonian/Overlap matrix
      call re_alloc(H,1,maxnh,1,nspin,name='H',
     .              routine='siesta',shrink=.false.,copy=.false.)
      call re_alloc(H0,1,maxnh,name='H0',routine='siesta',
     .              shrink=.false.,copy=.false.)
      call re_alloc(S,1,maxnh,name='S',routine='siesta',
     .              shrink=.false.,copy=.false.)

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
     .               maxnh, numh, listhptr, listh, xijo )
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

C If auxiliary cell has changed, optionally reset density matrix
C and set usesavedata to false to avoid reading back saved copy
      if (initdmaux.and.auxchanged) then
        dminit = .true.
        usesavedmloc = .false.
      else
        usesavedmloc = usesavedm
      endif

      if (dminit)
     .   call initdm(Datm, Dscf, Dold, lasto, na_s,
     .               maxnh, no_s, no_l, nspin, na_u, no_l, nspin,
     .               numh, numhold, listhptr, listhptrold,
     .               listh, listhold, iaorb, found, inspn, 
     .               usesavedmloc, no_u)


C Initialize energy-density matrix to zero for first call to overfsm
      Escf(1:maxnh,1:nspin) = 0.0_dp

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
     .                               .or. (idyn .eq. 7))) then
        istpsave = istp
        istp = 2
        iord = 0
        itest = .true.
      endif
      if (.not.harrisfun)
     . call extrapol(istp, iord, nspin, no_s, no_l, maxnh, 
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
C Increase nauxpul with safety margin when atoms move
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
     .             fal, stress, S )
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

        if (cml_p) call cmlStartStep(xf=mainXML, type='SCF')

C Normalize density matrix to exact charge ...........................
        qsol = 0.0_dp
        do ispin = 1,min(nspin,2)
          do io = 1,nh
            qsol = qsol + Dscf(io,ispin) * s(io)
          enddo
        enddo
#ifdef MPI
#ifdef NODAT
        call MPI_AllReduce(qsol,buffer1,1,MPI_double_precision,
     .    MPI_sum,MPI_Comm_World,MPIerror)
#else
        call MPI_AllReduce(qsol,buffer1,1,DAT_double,
     .    MPI_sum,MPI_Comm_World,MPIerror)
#endif
        qsol = buffer1
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
     .                 Ekin, fal, stressl, H ) 
#ifdef MPI
C Global reduction of energy terms
#ifdef NODAT
          call MPI_AllReduce(Ekin,buffer1,1,MPI_double_precision,
     .      MPI_sum,MPI_Comm_World,MPIerror)
#else
          call MPI_AllReduce(Ekin,buffer1,1,DAT_double,
     .      MPI_sum,MPI_Comm_World,MPIerror)
#endif
          Ekin = buffer1
#endif
        endif
C ..................

C Non-local-pseudop: energy, forces, stress and matrix elements .......
        if (first.or.last) then
          call nlefsm(scell, na_u, na_s, isa, xa, indxua, maxna,
     .                maxnh, maxnh, lasto, lastkb, iphorb, iphKB, 
     .                numh, listhptr, listh, numh, listhptr, listh, 
     .                min(nspin,2), Dscf, Enl, fal, stressl, H)
        
C Check whether maxna has been increased during nlefsm & resize if needed
          if (maxna.gt.size(jna)) then
            call memory('D','I',size(jna),'siesta')
            deallocate(jna)
            call memory('D','D',size(r2ij),'siesta')
            deallocate(r2ij)
            call memory('D','D',size(xij),'siesta')
            deallocate(xij)
            allocate(jna(maxna))
            call memory('A','I',maxna,'siesta')
            allocate(r2ij(maxna))
            call memory('A','D',maxna,'siesta')
            allocate(xij(3,maxna))
            call memory('A','D',3*maxna,'siesta')
          endif

#ifdef MPI
C Global reduction of energy terms
#ifdef NODAT
          call MPI_AllReduce(Enl,buffer1,1,MPI_double_precision,
     .      MPI_sum,MPI_Comm_World,MPIerror)
#else
          call MPI_AllReduce(Enl,buffer1,1,DAT_double,MPI_sum,
     .      MPI_Comm_World,MPIerror)
#endif
          Enl = buffer1
#endif
        endif
C ..................

C Save or get partial Hamiltonian (non-SCF part) ......................
        if (first.or.last) then
          do io = 1,nh
            H0(io) = H(io,1)
          enddo
        else
          do ispin = 1,nspin
            if (ispin .le. 2) then
              do io = 1,nh
                H(io,ispin) = H0(io)
              enddo
            else
              do io = 1,nh
                H(io,ispin) = 0.0_dp
              enddo
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
#ifdef NODAT
          call MPI_AllReduce(E0,buffer1,1,MPI_double_precision,
     .      MPI_sum,MPI_Comm_World,MPIerror)
#else
          call MPI_AllReduce(E0,buffer1,1,DAT_double,
     .      MPI_sum,MPI_Comm_World,MPIerror)
#endif
          E0 = buffer1
#endif
          E0 = E0 - Eions + Ena
        endif
C ..................

C Non-local-pseudop: energy, forces, stress and matrix elements .......
C Add SCF contribution to energy and matrix elements ..................
        g2max = g2cut
        if (last) then
C Last call to dhscf and grid-cell sampling if requested
          ifa  = 1
          istr = 1
          call grdsam( nspin, no_s, iaorb, iphorb, 
     .                 no_l, no_u, na_u, na_s, isa, xa, indxua,
     .                 ucell, mscell, g2max, ntm, ifa, istr, maxnh,
     .                 maxnh, numh, listhptr, listh, Dscf, Datm, H,
     .                 Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     .                 Exc, Dxc, dipol, fa, stress, fal, stressl)
        else
          ifa  = 0
          istr = 0
          ihmat = 1
          call dhscf( nspin, no_s, iaorb, iphorb, no_l,
     .                no_u, na_u, na_s, isa, xa, indxua, 
     .                ucell, mscell, g2max, ntm,
     .                ifa, istr, ihmat, ' ', ' ', ' ', ' ', ' ', ' ',
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
     .                 fal, stressl, S ) 
        endif
C ..................

C Find entropy ........................................................
C Entropy is returned from the call to diagon. To add to the energy
C the entropy computed from the input charge, here it is assigned to the one
C of the former SCF step

        Entropy = 0.0_dp
        if (isolve .eq. 0) then
          if (istp.gt.1 .or. iscf.gt.1) then
            Entropy = Entrop
          endif
        endif

C Save present density matrix ........................................
        do is = 1,nspin
          do io = 1,nh
            Dold(io,is) = Dscf(io,is)
            Eold(io,is) = Escf(io,is)
          enddo
        enddo

C Save Hamiltonian and overlap matrices ............................
        if (savehs) then
          call iohs( 'write', gamma, no_u, nspin, maxnh,
     .               numh, listhptr, listh, H, S, qtot, temp, xijo )
        endif

C Solve eigenvalue problem .........................................
        if (.not.last) then
          if (isolve .eq. 0) then
            call diagon(no_s, nspin, nspin, no_l, maxnh, maxnh, no_u,
     .                  numh, listhptr, listh, numh, listhptr, listh, 
     .                  H, S, qtot, fixspin, qs, temp, e1, e2,
     .                  gamma, xijo, indxuo, nkpnt, kpoint, kweight,
     .                  eo, qo, Dscf, Escf, ef, efs, Entrop, no_u,
     .                  occtol, iscf, neigwanted)
            Ecorrec = 0.0_dp
            printeigen = .true.
          elseif (isolve .eq. 1) then
            if (.not. gamma) call die("Cannot do O(N) with k-points.")
            call ordern(usesavelwf,ioptlwf,na_s,no_s,no_l,lasto,isa,qa,
     .                  rcoor,rmaxh,ucell,xa,iscf,istp,ncgmax,etol,eta,
     .                  qtot,maxnh,numh,listhptr,listh,H,S,
     .                  chebef,noeta,rcoorcp,beta,pmax,Dscf,Escf,
     .                  Ecorrec,nspin,qs)
            Entrop = 0.0_dp
          else
            call die('siesta: ERROR: wrong solution method')
          endif

C Harris-functional energy ............................................
          DEharr = 0.0_dp
          do ispin = 1,nspin
C const factor takes into account that there are two nondiagonal
C elements in non-collinear spin density matrix, stored as
C ispin=1 => D11; ispin=2 => D22, ispin=3 => Real(D12);
C ispin=4 => Imag(D12)
            const = 1._dp
            if (ispin .gt. 2) const = 2._dp
            do io = 1,nh
              DEharr = DEharr + H(io,ispin) * const * 
     .                     ( Dscf(io,ispin) - Dold(io,ispin) )
            enddo
          enddo
#ifdef MPI
C Global reduction of DEharr
#ifdef NODAT
          call MPI_AllReduce(DEharr,buffer1,1,MPI_double_precision,
     .      MPI_sum,MPI_Comm_World,MPIerror)
#else
          call MPI_AllReduce(DEharr,buffer1,1,DAT_double,
     .      MPI_sum,MPI_Comm_World,MPIerror)
#endif
          DEharr = buffer1
#endif
C ..................

C Print populations at each SCF step if requested before mixing ......

          if (muldeb) then
             write (6,"(/a)")
     .             'siesta: Mulliken populations before mixing'
             call mulliken( mullipop, nspin, na_u, no_u, maxnh,
     .                      numh, listhptr, listh, S, Dscf, isa,
     .                      lasto, iaorb, iphorb )
          endif
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
#ifdef NODAT
          call MPI_AllReduce(dDmax,buffer1,1,MPI_double_precision,
     .      MPI_max,MPI_Comm_World,MPIerror)
#else
          call MPI_AllReduce(dDmax,buffer1,1,DAT_double,
     .      MPI_max,MPI_Comm_World,MPIerror)
#endif
          dDmax = buffer1
#endif
C ...................

C Print populations at each SCF step, if requested, after mixing ......

          if (muldeb) then 
             write (6,"(/a)")
     .             'siesta: Mulliken populations after mixing'
             call mulliken( mullipop, nspin, na_u, no_u, maxnh,
     .                      numh, listhptr, listh, S, Dscf, isa,
     .                      lasto, iaorb, iphorb )
          endif
C ..................

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
        endif !last
C ...................
        DEna = Enascf - Enaatm
        Etot = E0 + DEna + DUscf + DUext + Exc + Ecorrec + Emad
        Eharrs = Etot + DEharr
        FreeE  = Etot - Temp * Entropy
C       recalculating the energy in the last iter (for gridcellsampling)
C       but preserving the value of Eharrs1
        if (.not.last) Eharrs1 = Eharrs

        if (IOnode.and..not.last) then
          call siesta_write_energies()

          if (harrisfun) then 
            write(6,"(/a,f14.6,/)") 'siesta: Eharris(eV) = ', Eharrs/eV
            if (cml_p) then
              call cmlStartPropertyList(mainXML, title='SCF Cycle')
              call cmlAddProperty(xf=mainXML, property=Eharrs/eV, 
     .             units="eV", dictRef="siesta:Eharrs", fmt="(f14.7)")
              call cmlEndPropertyList(mainXML)
            endif
          endif
        endif
C ...................

! End of one SCF step - flush stdout
        if (ionode) call pxfflush(6)

C If last iteration, exit SCF loop ....................................
        if (last) then
          do ispin = 1,nspin
            do io = 1,nh
              Dscf(io,ispin) = Dold(io,ispin)
              Escf(io,ispin) = Eold(io,ispin)
            enddo
          enddo
          if (dumpcharge) then
             call plcharge( no_s, na_s, no_u, maxnh, maxna, nspin,
     .                      isa, iphorb, indxuo, lasto,
     .                      scell, nsc, xa, rmaxo, datm )
          endif
          call timer( 'IterSCF', 2 )
          if (cml_p) call cmlEndStep(mainXML)
          goto 50
        endif
C ...................

C If converged, make last iteration to find forces ....................
        if (dDmax.lt.dDtol) last = .true.
C ...................

        call timer( 'IterSCF', 2 )
        if (istep.eq.inicoor .and. first) call timer( 'IterSCF', 3 )
        first = .false.
        if (cml_p) call cmlEndStep(mainXML)
      enddo
   50 continue
C End of SCF iteration_________________________________________________

C Write final Kohn-Sham Energy ........................................
      if (cml_p) call cmlStartPropertyList(mainXML, 
     .                                   title='Final KS Energy')
      if (IOnode) then
        if ( .not. harrisfun) 
     .    write(6,"(/a,f14.4)")  'siesta: E_KS(eV) =        ', Etot/eV
        if (signat) call sig(io_sig,'Final KS energy',Etot/eV)
        if (cml_p) call cmlAddProperty(xf=mainXML, property=Etot/eV,
     .       dictref='siesta:E_KS', units='eV', fmt='(f14.6)')
      endif

C Substract egg box effect form energy ................................
      if (eggbox_block) then
        call eggbox('energy',ucell,na_u,isa,ntm,xa,fa,Etot,
     .               eggbox_block)
        if (IOnode) 
     .    write(6,"(/a,f14.4)") 'siesta: E_KS - E_eggbox = ',Etot/eV
        if (cml_p) call cmlAddProperty(xf=mainXML, property=Etot/eV, 
     .         dictref='siesta:E_KS_egg', units='eV', fmt='(f14.6)')
      endif
      if (cml_p) call cmlEndPropertyList(mainXML)

#ifdef MPI
C Global reduction of forces and stresses
      allocate(fatmp(3,na_u))
      call memory('A','D',3*na_u,'siesta')
#ifdef NODAT
      call MPI_AllReduce(stressl(1,1),stresstmp(1,1),9,
     .  MPI_double_precision,MPI_sum,MPI_Comm_World,MPIerror)
      call MPI_AllReduce(fal(1,1),fatmp(1,1),3*na_u,
     .  MPI_double_precision,MPI_sum,MPI_Comm_World,MPIerror)
#else
      call MPI_AllReduce(stressl(1,1),stresstmp(1,1),9,
     .  DAT_double,MPI_sum,MPI_Comm_World,MPIerror)
      call MPI_AllReduce(fal(1,1),fatmp(1,1),3*na_u,
     .  DAT_double,MPI_sum,MPI_Comm_World,MPIerror)
#endif
      stress(1:3,1:3) = stress(1:3,1:3) + stresstmp(1:3,1:3)
      fa(1:3,1:na_u) = fa(1:3,1:na_u) + fatmp(1:3,1:na_u)
      call memory('D','D',size(fatmp),'siesta')
      deallocate(fatmp)
#else
      stress(1:3,1:3) = stress(1:3,1:3) + stressl(1:3,1:3)
      fa(1:3,1:na_u) = fa(1:3,1:na_u) + fal(1:3,1:na_u)
#endif

C Substract egg box effect from the forces ............................
      if (eggbox_block) then
        call eggbox('forces',ucell,na_u,isa,ntm,xa,fa,Etot,eggbox_block)
      endif
C ...................

C Impose constraints to atomic movements by changing forces ...........
      call fixed(ucell,stress,na_u,isa, amass, xa, fa, cstress, cfa, 
     .           ntcon )
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

C Add kinetic term to stress tensor ..................................
      do ia = 1,na_u
        do jx = 1,3
          do ix = 1,3
            tstress(ix,jx) = stress(ix,jx) -
     .             amu * amass(ia) * va(ix,ia) * va(jx,ia) / volume
          enddo
        enddo
      enddo
C ...................

      if (IOnode) then
        call siesta_write_forces()
        call siesta_write_stress_pressure()
      endif
C ...................

C Mulliken population analysis .......................................
      call mulliken( mullipop, nspin, na_u, no_u, maxnh,
     .               numh, listhptr, listh, S, Dscf, isa, 
     .               lasto, iaorb, iphorb )

C ...................

C Move atoms ..........................................................

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

      select case(idyn)
      case(0)
        if (nmove .ne. 0) then
          call cgvc( na_u, xa, cfa, ucell, cstress, volume, dxmax,
     .         tp, ftol, strtol, varcel, relaxd, usesavecg )
          ! Propagate the new structure to the virtual supercell
          call superx( ucell, nsc, na_u, na_s, xa, scell )
          if (relaxd) goto 60
            ! Exit coordinate relaxation loop
        endif
!----------
      case(1)
        call verlet2(istp, iunit, iquench, na_u, cfa, dt, 
     .       amass, ntcon, va, xa, Ekinion, tempion)
        ! Propagate the new structure to the virtual supercell
        call superx( ucell, nsc, na_u, na_s, xa, scell )
!-----------
      case (2)
        call nose(istp, iunit, na_u, cfa, tt, dt, amass, mn, 
     .       ntcon, va, xa, Ekinion, kn, vn, tempion)
        ! Propagate the new structure to the virtual supercell
        call superx( ucell, nsc, na_u, na_s, xa, scell )
!-----------
      case (3)
        call pr(istp, iunit, iquench, na_u, cfa, cstress, tp, dt, amass,
     .       mpr, ntcon, va, xa, vcell, ucell, Ekinion, kpr, vpr, 
     .       tempion, Pint)
        ! Propagate the new structure to the virtual supercell
        call superx( ucell, nsc, na_u, na_s, xa, scell )
        if (IOnode) write(6,'(/,a,f12.3,a)')
     .       'siesta: E_kin PR =', kpr/Kelvin, ' K'
!-----------
      case (4)
        call npr(istp, iunit, na_u, cfa, cstress, tp, tt, dt, amass, mn,
     .       mpr, ntcon, va, xa, vcell, ucell, Ekinion, kn, kpr, vn,
     .       vpr, tempion, Pint)
        ! Propagate the new structure to the virtual supercell
        call superx( ucell, nsc, na_u, na_s, xa, scell )
!-----------
      case (5)
        call anneal(istp, iunit, ianneal, taurelax, bulkm,
     .       na_u, cfa, cstress, tp, tt, dt, amass, ntcon,
     .       va, xa, ucell, Ekinion, tempion, Pint)
        ! Propagate the new structure to the virtual supercell
        call superx( ucell, nsc, na_u, na_s, xa, scell )
!-----------
      case (6:7)
        continue !We can't go until after ioxv - see below
!-----------
      case (8)
        call forcesToPipe( na_u, Etot, cfa, cstress )
      end select

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

c restore original coordinates after FC displacements ...
      if (idyn .eq. 6 .and. istep .ne. 0) then
        xa(ixdispl,iadispl)=xa(ixdispl,iadispl)-dx
      endif
      if (idyn .eq. 7) then
           call phonon_restore_coords(istep,xa,ucell)
      endif
c ...

C Save atomic positions and velocities accumulatively ................
      if (writmd.and.IOnode) then
         if ( .not. harrisfun) then
            getot = Etot + Ekinion + kn + kpr + vn + vpr
            call iomd( na_u, isa, iza, 
     .        xa, va, ucell, vcell, varcel, istep, inicoor, 
     .        fincoor, tempion, Etot, getot, volume/Ang**3, Psol/kbar)
         else
            getot = Eharrs1 + Ekinion + kn + kpr + vn + vpr
            call iomd( na_u, isa, iza, 
     .        xa, va, ucell, vcell, varcel, istep, inicoor, 
     .        fincoor, tempion, Eharrs1, getot, volume/Ang**3,Psol/kbar) 
         endif
      endif

C Accumulate coor in Xmol file for animation .........................
      lastst = fincoor .le. istep
      if (writpx.and.IOnode) 
     .  call pixmol(iza, xa, na_u, slabel, lastst)
C ...................

   60 continue
C Output memory use at the end of this geometry step
      if (cml_p) call cmlEndStep(mainXML)
      call printmemory( 6, 0 )
      call timer( 'IterMD', 2 )

! End of one MD step - flush stdout
      if (ionode) call pxfflush(6)

      if (relaxd) exit

      enddo
C End of coordinate-relaxation loop ==================================
      final=.true.

      ! We want xalast to equal xa for coordinate relaxation only.
      if (idyn==0) xalast(1:3,1:na_s)=xa(1:3,1:na_s)

      if (cml_p) then
        call cmlStartModule(xf=mainXML, title='Finalization')
      endif

      if (IOnode) then
C Print atomic coordinates (and also unit cell for ParrRah.)
        if (nmove .ne. 0) then
          if (relaxd) 
     .      call outcoor(ucell, xa, isa, na_u, 'Relaxed', .true. )
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

C Find and print wavefunctions at selected k-points
      if (nwk.gt.0) then
        call wwave( no_s, nspin, nspin, no_u, no_l, maxnh, maxwk,
     .              numh, listhptr, listh, H, S, Ef, xijo, indxuo,
     .              nwk, wfk, no_u, gamma, occtol )
      endif

C Find and print bands
      if (nbk.gt.0) then
        call bands( no_s, nspin, nspin, no_u, no_l, maxnh, maxbk,
     .              numh, listhptr, listh, H, S, Ef, xijo, indxuo,
     .              .true., nbk, bk, ebk, no_u, occtol )
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
      if (IOnode .and. printeigen) then
        if (isolve.eq.0 .and. no_l.lt.1000) then
          if (nspin .le. 2) then
            write(6,'(/,a,/,a4,a3,a7)')
     .       'siesta: Eigenvalues (eV):', 'ik', 'is', 'eps'
            do ik = 1,nkpnt
              do ispin = 1,nspin
                write(6,'(i4,i3,10f7.2)')
     .            ik,ispin,(eo(io,ispin,ik)/eV,io=1,min(10,neigwanted))
                if (no_u.gt.10) write(6,'(7x,10f7.2)')
     .            (eo(io,ispin,ik)/eV,io=11,neigwanted)
              enddo
            enddo
          else
            write(6,'(/,a)') 'siesta: Eigenvalues (eV):'
            do ik = 1,nkpnt
              write(6,'(a,i6)') 'ik =', ik
              write(6,'(10f7.2)')
     .          ((eo(io,ispin,ik)/eV,io=1,neigwanted),ispin=1,2)
            enddo
          endif
          write(6,'(a,f15.6,a)') 'siesta: Fermi energy =', ef/eV, ' eV'
        endif
      endif

      if (isolve.eq.0.and.IOnode.and.printeigen) 
     .     call ioeig(eo,ef,neigwanted,nspin,nkpnt,no_u,nspin,maxk,
     .                kpoint, kweight)

C Compute the projected density of states
      if (IOnode) then
        do_pdos = fdf_block('ProjectedDensityOfStates',iu)
        if (isolve.ne.0.and.do_pdos) then
          write(6,*)
     .         'siesta: ERROR: PDOS implemented only with diagon'
          do_pdos = .false.
        endif
      endif
#ifdef MPI
      call MPI_Bcast(do_pdos,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif

      if (do_pdos) then
C Find the desired energy range
        if (IOnode) then
          read(iu,'(a)') line
          p=>digest(line)
          if (nvalues(p).lt.3 .or. nnames(p).ne.1)
     $          call die("Wrong format in PDOS block")
          factor = fdf_convfac( names(p,1), 'Ry' )
          e1 = values(p,1) * factor
          e2 = values(p,2) * factor
          sigma = values(p,3) * factor
          nhist = integers(p,1)
          write(6,'(a)') 'siesta: PDOS info: '
          write(6,'(a,3(f8.2,a),2x,i5)')
     $           'siesta: e1, e2, sigma, nhist: ',
     $           e1/eV,' eV',e2/eV,' eV',sigma/eV,' eV', nhist
        endif

#ifdef MPI
C Broadcast values to all nodes
#ifdef NODAT
        call MPI_Bcast(e1,1,MPI_double_precision,0,MPI_Comm_World,
     .                 MPIerror)
        call MPI_Bcast(e2,1,MPI_double_precision,0,MPI_Comm_World,
     .                 MPIerror)
        call MPI_Bcast(sigma,1,MPI_double_precision,0,MPI_Comm_World,
     .                 MPIerror)
#else
        call MPI_Bcast(e1,1,DAT_double,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(e2,1,DAT_double,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(sigma,1,DAT_double,0,MPI_Comm_World,MPIerror)
#endif
        call MPI_Bcast(nhist,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif

        call pdos( no_s, nspin, nspin, no_l, maxnh,
     .             no_u, numh, listhptr, listh, H, S,
     .             e1, e2, sigma, nhist,
     .             gamma, xijo, indxuo, nkpnt, kpoint, kweight, eo,
     .             no_u)

      endif                     ! PDOS calc (do_pdos)

C Print program's energy decomposition and final forces
      if (IOnode) then
        call siesta_write_energies()
        call siesta_write_forces()
        call siesta_write_stress_pressure()
      endif

C Print spin polarization
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
#ifdef MPI
C Global reduction of spin components
#ifdef NODAT
        call MPI_AllReduce(qspin,qtmp,nspin,
     .    MPI_double_precision,MPI_sum,MPI_Comm_World,MPIerror)
#else
        call MPI_AllReduce(qspin,qtmp,nspin,
     .    DAT_double,MPI_sum,MPI_Comm_World,MPIerror)
#endif
        qspin(1:nspin) = qtmp(1:nspin)
#endif
        if (nspin .eq. 2) then
          if (IOnode) then
            write(6,'(/,a,f12.6)')
     .       'siesta: Total spin polarization (Qup-Qdown) =', 
     .       qspin(1) - qspin(2)
          endif
          if (cml_p) call cmlAddProperty(xf=mainXML, 
     .         property=qspin(1)-qspin(2), dictref='siesta:qspin')
        elseif (nspin .eq. 4) then
          call spnvec( nspin, qspin, qaux, stot, svec )
          if (IOnode) then
            write(6,'(/,a,f12.6)')
     .       'siesta: Total spin polarization (Qup-Qdown) =', stot
            if (signat) call sig(io_sig,'Tot spin pol',stot)
            write(6,'(a,3f12.6)') 'siesta: Spin vector =', svec
            if (cml_p) then
              call cmlAddProperty(xf=mainXML, property=stot, 
     .             dictref='siesta:stot')
              call cmlAddProperty(xf=mainXML, property=svec,
     .             dictref='siesta:svec')
            endif !cml_p
          endif
        endif
      endif

C Print electric dipole
      if (shape .ne. 'bulk') then
        if (IOnode) then
          write(6,'(/,a,3f12.6)')
     .      'siesta: Electric dipole (a.u.)  =', dipol
          write(6,'(a,3f12.6)')
     .      'siesta: Electric dipole (Debye) =', 
     .      (dipol(ix)/Debye,ix=1,3)
        endif
        if (cml_p) then
          call cmlAddProperty(xf=mainXML, property=dipol,
     .         title='Electric dipole', dictref='siesta:dipol', 
     .         units='a.u')
        endif !cml_p
      endif

C Calculation of the bulk polarization using the Berry phase
C formulas by King-Smith and Vanderbilt
C Attention H0 is used as an auxiliary array
      if (nkpol.gt.0) then
        call KSV_pol(na_u, na_s, xalast, rmaxo, scell, ucell,
     .               no_u, no_l, no_s, nspin, qspin, maxna, 
     .               maxnh, nkpol, numh, listhptr, listh, 
     .               H, S, H0, xijo, indxuo, isa, iphorb, 
     .               iaorb, lasto, jna, xij, r2ij,shape,
     .               nkpol,kpol,wgthpol, polR, polxyz ) 
      endif

C Calculation of the optical conductivity
C Attention H0, Eold, Dold are used as auxiliary arrays
      call optical(na_u, na_s, xa, scell, ucell,
     .             no_u, no_l, no_s, nspin, qspin,
     .             maxna, maxnh, numh, listhptr, listh, H, S, H0,
     .             Eold(1,1), Dold(1,1),
     .             xijo, indxuo, indxua, ebk, ef, temp,
     .             isa, iphorb, iphKB, iaorb, lasto, lastkb,
     .             jna, xij, r2ij, shape, slabel)

c...................................

C Save electron density and potential
      if (IOnode) then
        savrho = fdf_boolean( 'SaveRho',    dumpcharge .or. .false.)
        savdrh = fdf_boolean( 'SaveDeltaRho',               .false.)
        savevh = fdf_boolean( 'SaveElectrostaticPotential', .false.)
        savevt = fdf_boolean( 'SaveTotalPotential',         .false.)
        savepsch = fdf_boolean( 'SaveIonicCharge',          .false.)
        savetoch = fdf_boolean( 'SaveTotalCharge',          .false.)
      endif
#ifdef MPI
      call MPI_Bcast(savrho,1,MPI_logical,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(savdrh,1,MPI_logical,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(savevh,1,MPI_logical,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(savevt,1,MPI_logical,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(savepsch,1,MPI_logical,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(savetoch,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
      if (savrho .or. savdrh .or. savevh .or. savevt .or.
     .    savepsch .or. savetoch ) then
        filrho = ' '
        fildrh = ' '
        filevh = ' '
        filevt = ' '
        filepsch = ' '
        filetoch = ' '
        if (savrho) filrho = paste( slabel, '.RHO' )
        if (savdrh) fildrh = paste( slabel, '.DRHO' )
        if (savevh) filevh = paste( slabel, '.VH'  )
        if (savevt) filevt = paste( slabel, '.VT'  )
        if (savepsch) filepsch = paste( slabel, '.IOCH'  )
        if (savetoch) filetoch = paste( slabel, '.TOCH'  )
        g2max = g2cut
        call dhscf( nspin, no_s, iaorb, iphorb, no_l,
     .              no_u, na_u, na_s, isa, xa, indxua, 
     .              ucell, mscell, g2max, ntm,
     .              0, 0, 0, filrho, fildrh, filevh, filevt,
     .              filepsch, filetoch, 
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

C Find the desired energy range
        if (IOnode) then
          read(iu,'(a)') line
          p=>digest(line)
          if (.not. match(p,"vvn"))
     .       call die("Wrong format in LocalDensityofStates")
          factor = fdf_convfac( names(p,1), 'Ry' )
          e1 = values(p,1)*factor
          e2 = values(p,2)*factor
          call destroy(p)
        endif
#ifdef MPI
#ifdef NODAT
        call MPI_Bcast(e1,1,MPI_double_precision,0,
     .    MPI_Comm_World,MPIerror)
        call MPI_Bcast(e2,1,MPI_double_precision,0,
     .    MPI_Comm_World,MPIerror)
#else
        call MPI_Bcast(e1,1,DAT_double,0,
     .    MPI_Comm_World,MPIerror)
        call MPI_Bcast(e2,1,DAT_double,0,
     .    MPI_Comm_World,MPIerror)
#endif
#endif

C Find the density matrix for states between e1 and e2
        if (isolve .eq. 0) then
          call diagon(no_s, nspin, nspin, no_l, maxnh, maxnh, no_u,
     .                numh, listhptr, listh, numh, listhptr, listh, 
     .                H, S, qtot, fixspin, qs, temp, e1, e2,
     .                gamma, xijo, indxuo, nkpnt, kpoint, kweight,
     .                eo, qo, Dscf, Escf, ef, efs, Entrop, no_u,
     .                occtol, iscf, neigwanted)
        else
          if (IOnode) then
            write(6,*)
     .       'siesta: ERROR: LDOS implemented only with diagon'
          endif
          goto 70
        endif
     
C Find the LDOS in the real space mesh
        filrho = paste( slabel, '.LDOS' )
        g2max = g2cut
        call dhscf( nspin, no_s, iaorb, iphorb, no_l,
     .              no_u, na_u, na_s, isa, xa, indxua, 
     .              ucell, mscell, g2max, ntm,
     .              0, 0, 0, filrho, ' ', ' ', ' ', ' ', ' ',
     .              maxnh, numh, listhptr, listh, Dscf, Datm,
     .              maxnh, numh, listhptr, listh, H,
     .              Enaatm, Enascf, Uatm, Uscf, DUscf, DUext, Exc, Dxc,
     .              dipol, fa, stress, fal, stressl )
      endif
   70 continue

C Output memory use up to the end of the program
      call printmemory( 6, 1 )

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

      if (cml_p) then
        call cmlEndModule(mainXML)
        call siesta_cml_exit()
      endif

! End of program
! Internal subroutines follow.

      contains
      
      subroutine siesta_write_k_points()
      if ( writek ) then
        write(6,'(/,a)')
     .       'siesta: k-point coordinates (Bohr**-1) and weights:'
        write(6,'(a,i4,3f12.6,3x,f12.6)')
     .       ('siesta: ', ik, (kpoint(ix,ik),ix=1,3), kweight(ik),
     .       ik=1,nkpnt)
        if (cml_p) then
          call cmlAddProperty(xf=mainXML, property=kpoint, 
     .         dictref='siesta:kpoint')
          call cmlAddProperty(xf=mainXML, property=kweight,
     .         dictref='siesta:kweight')
        endif
      else
        call iokp( nkpnt, kpoint, kweight )
      endif
      write(6,'(/a,i6)')
     .     'siesta: k-grid: Number of k-points =', nkpnt
      write(6,'(a,f10.3,a)')
     .     'siesta: k-grid: Cutoff             =', kcutof/Ang, ' Ang'
      write(6,'(a)')
     .     'siesta: k-grid: Supercell and displacements'
      write(6,'(a,3i4,3x,f8.3)') 'siesta: k-grid: ',
     .     (kscell(i,1),i=1,3), kdispl(1)
      write(6,'(a,3i4,3x,f8.3)') 'siesta: k-grid: ',
     .     (kscell(i,2),i=1,3), kdispl(2)
      write(6,'(a,3i4,3x,f8.3)') 'siesta: k-grid: ',
     .     (kscell(i,3),i=1,3), kdispl(3)
      if (cml_p) then
        call cmlStartPropertyList(mainXML, title='k-points')
        call cmlAddProperty(xf=mainXML, property=nkpnt, 
     .       dictref='siesta:nkpnt')
        call cmlAddProperty(xf=mainXML, property=kcutof/Ang, 
     .       dictref='siesta:kcutof', units='Ang')
        call cmlEndPropertyList(mainXML)
        call cmlAddProperty(xf=mainXML, property=kscell,
     .       dictref='siesta:kscell')
        call cmlAddProperty(xf=mainXML, property=kdispl,
     .       dictref='siesta:kdispl')
      endif
      end subroutine siesta_write_k_points


      subroutine siesta_write_forces()

      ! Almost the same forces output whether during simulation
      ! or at the end. Unfortunately not quite, therefore slightly
      ! tortuous logic below. If we are content to change format
      ! of output file slightly, this can be simplified.
      if (.not.final) then
        ! print forces to xml every step. 
        ! output forces to stdout depending on writef
        if (cml_p) then
          call cmlStartPropertyList(mainXML, title='Forces')
          call cmlAddProperty(xf=mainXML, property=fa*Ang/eV,
     .         dictref='siesta:forces', units='evpa')
          call cmlAddProperty(xf=mainXML, property=ftot,
     .         dictref='siesta:ftot')
          call cmlAddProperty(xf=mainXML, property=fmax, 
     .         dictref='siesta:fmax')
          call cmlAddProperty(xf=mainXML, property=fres,
     .         dictref='siesta:fres')
          call cmlAddProperty(xf=mainXML, property=cfmax, 
     .         dictref='siesta:cfmax')
          call cmlEndPropertyList(mainXML)
        endif
        write(6,'(/,a)') 'siesta: Atomic forces (eV/Ang):'
        if (writef) then
          write(6,'(i6,3f12.6)')(ia,(fa(ix,ia)*Ang/eV,ix=1,3),ia=1,na_u)
        else
          call iofa( na_u, fa )
        endif
        write(6,'(40("-"),/,a4,3f12.6)') 'Tot',(ftot(ix)*Ang/eV,ix=1,3)
        write(6,'(40("-"),/,a4, f12.6)') 'Max',fmax*Ang/eV
        if (signat) call sig(io_sig,'Maximum force',fmax*Ang/eV)
        write(6,'(a4,f12.6,a)')'Res',fres*Ang/eV,
     .       '    sqrt( Sum f_i^2 / 3N )'
        write(6,'(40("-"),/,a4, f12.6,a)') 'Max',cfmax*Ang/eV, 
     .       '    constrained'
      else !not final
C In finalization, only print forces if sufficiently large.
        fmax = maxval(abs(fa))
        ftot = sum(fa, dim=2)
        if (fmax .gt. ftol) then
          write(6,'(/,a)') 'siesta: Atomic forces (eV/Ang):'
          write(6,'(a,i6,3f12.6)')
     .         ('siesta: ', ia,(fa(ix,ia)*Ang/eV,ix=1,3),ia=1,na_u)
          write(6,'(a,40("-"),/,a,a4,3f12.6)')
     .         'siesta: ','siesta: ','Tot',(ftot(ix)*Ang/eV,ix=1,3)
          if (cml_p) then
            call cmlStartPropertyList(mainXML, title='Force Summary')
            call cmlAddProperty(xf=mainXML, property=fa*Ang/eV,
     .           dictref='siesta:forces', units='evpa')
            call cmlAddProperty(xf=mainXML, property=ftot*Ang/eV, 
     .           dictref='siesta:ftot', units='evpa')
            call cmlEndPropertyList(mainXML)
          endif !cml_p
        endif
        if (Any(cfa /= fa)) then
          fmax = maxval(abs(cfa))
          ftot = sum(cfa, dim=2)
          if (fmax .gt. ftol) then
            write(6,'(/,a)') 'siesta: Constrained forces (eV/Ang):'
            write(6,'(a,i6,3f12.6)')
     .           ('siesta: ',ia,(cfa(ix,ia)*Ang/eV,ix=1,3),ia=1,na_u)
            write(6,'(a,40("-"),/,a,a4,3f12.6)')
     .           'siesta: ','siesta: ','Tot',(ftot(ix)*Ang/eV,ix=1,3)
            if (cml_p) then
              call cmlStartPropertyList(mainXML,
     .             title='Constrained Force Summary')
              call cmlAddProperty(xf=mainXML, property=cfa*Ang/eV, 
     .             dictref='siesta:cforces', units='evpa')
              call cmlAddProperty(xf=mainXML, property=ftot*Ang/eV, 
     .             dictref='siesta:cftot', units='evpa')
              call cmlEndPropertyList(mainXML)
            endif !cml_p
          endif
        endif
      endif !final for forces

      end subroutine siesta_write_forces


      subroutine siesta_write_stress_pressure()
! Stress tensor and pressure:
      
      if (.not.final) then
      ! Output depends on dynamics option
        select case (idyn)
        case(0:5)
          if (idyn==0 .and..not.varcel) then
            continue
          else
            write(6,'(/,a,3(/,a,3f12.6))')
     .           'siesta: Stress tensor (static) (eV/Ang**3):',
     .           ('     ',(stress(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)
            Psol = - ((stress(1,1) + stress(2,2) + stress(3,3))/3.0_dp)
            write(6,'(/,a,f20.8,a)')
     .           'siesta: Pressure (static):', Psol/kBar, '  kBar'
            if (cml_p) then
              call cmlAddProperty(xf=mainXML, property=stress*Ang**3, 
     .             dictref='siesta:stress')
              call cmlAddProperty(xf=mainXML, property=Psol, 
     .             dictref='siesta:psol', title='Pressure (Static)')
            endif !cml_p
            if (signat) call sig(io_sig,'Pressure kbar',Psol/kBar)
            write(6,'(/,a,3(/,a,3f12.6))')
     .           'siesta: Stress tensor (total) (eV/Ang**3):',
     .           ('     ',(tstress(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)
            Psol = - ((tstress(1,1)+tstress(2,2) +tstress(3,3))/3.0_dp)
            write(6,'(/,a,f20.8,a)')
     .           'siesta: Pressure (total):', Psol/kBar, '  kBar'
            if (cml_p) then
              call cmlAddProperty(xf=mainXML, property=tstress*Ang**3, 
     .             dictref='siesta:tstress')
              call cmlAddProperty(xf=mainXML, property=Psol,
     .             dictref='siesta:tpsol', title='Pressure (Total)')
            endif !cml_p
          endif !varcel
        ! Write Force Constant matrix if FC calculation ...
        case(6)
          call ofc(fa,dx,na_u)
        case(7)
          call phonon_write_forces(fa,na_u,ns,ucell,istep)
        end select !idyn

      else !final

C Print stress tensor unconditionally
        write(6,'(/,a,3(/,a,3f12.6))')
     .       'siesta: Stress tensor (static) (eV/Ang**3):',
     .       ('siesta: ',(stress(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)
        if (cml_p) then
          call cmlAddProperty(xf=mainXML, property=stress*Ang**3/eV, 
     .         dictref='siesta:stress', units='eV/Ang^3')
        endif !cml_p

C Print constrained stress tensor if different from unconstrained
        if (Any(cstress /= stress )) then
             write(6,'(/,a,3(/,a,3f12.6))')
     .       'siesta: Constrained stress tensor (static) (eV/Ang**3):',
     .       ('siesta: ',(cstress(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)
          if (cml_p) then
            call cmlAddProperty(xf=mainXML, property=cstress*Ang**3/eV, 
     .           dictref='siesta:cstress', units='eV/Ang^3')
          endif !cml_p
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
        write(6,'(/,a,f18.6,a)')
     .       'siesta: Cell volume =', volume/Ang**3, ' Ang**3'
        write(6,'(/,a,/,a,2a20,a,3(/,a,2f20.8,a))')
     .       'siesta: Pressure (static):',
     .       'siesta: ','Solid',        'Molecule',      '  Units',
     .       'siesta: ', Psol,           Pmol,           '  Ry/Bohr**3',
     .       'siesta: ', Psol*Ang**3/eV, Pmol*Ang**3/eV, '  eV/Ang**3',
     .       'siesta: ', Psol/kBar,      Pmol/kBar,      '  kBar'
        if (signat) call sig(io_sig,'Pressure kbar',Psol/kBar)
        if (cml_p) then
          call cmlStartPropertyList(mainXML, title='Final Pressure')
          call cmlAddProperty(xf=mainXML, property=volume/Ang**3, 
     .         title='cell volume', dictref='siesta:cellvol', 
     .         units='Ang^3')
          call cmlAddProperty(xf=mainXML, property=Psol/kBar, 
     .         title='Pressure of Solid', dictref='siesta:pressSol', 
     .         units='kbar')
          call cmlAddProperty(xf=mainXML, property=Pmol/kBar,       
     .         title='Pressure of Molecule', dictref='siesta:pressMol', 
     .         units='kbar')
          call cmlEndPropertyList(mainXML)
        endif !cml_p
        
      endif !final for stress & pressure

      end subroutine siesta_write_stress_pressure
      
      
      subroutine siesta_write_energies()
      ! Only print out full decomposition at very beginning and end.
      if ((istp==1.and.first).or.final) then
        write(6,'(/,a,/,(a,f17.6))')
     .     'siesta: Program''s energy decomposition (eV):',
     .     'siesta: Eions   =', Eions/eV,
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
     .     'siesta: Eharris =', (Eharrs1+Ekinion)/eV,
     .     'siesta: Etot    =', (Etot+Ekinion)/eV,
     .     'siesta: FreeEng =', (FreeE+Ekinion)/eV
        if (cml_p) then
          call cmlStartPropertyList(mainXML,
     .         title='Energy Decomposition')
          call cmlAddProperty(xf=mainXML, property=Eions/eV, units='eV', 
     .         dictref='siesta:Eions', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=Ena/eV, units='eV', 
     .         dictref='siesta:Ena', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=Ekin/eV, units='eV', 
     .         dictref='siesta:Ekin', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=Enl/eV, units='eV', 
     .         dictref='siesta:Enl', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=DEna/eV, units='eV', 
     .         dictref='siesta:DEna', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=DUscf/eV, units='eV',
     .         dictref='siesta:DUscf', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=DUext/eV, units='eV', 
     .         dictref='siesta:DUext', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=Exc/eV, units='eV', 
     .         dictref='siesta:Exc', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML,property=Ecorrec/eV,units='eV', 
     .         dictref='siesta:Ecorrec', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=Emad/eV, units='eV', 
     .         dictref='siesta:Emad', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML,property=Ekinion/eV,units='eV', 
     .         dictref='siesta:Ekinion', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=(Eharrs1+Ekinion)/eV,
     .         units='eV', dictref='siesta:EharrsK', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=(Etot+Ekinion)/eV, 
     .         units='eV', dictref='siesta:EtotK', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=(FreeE+Ekinion)/eV, 
     .         units='eV', dictref='siesta:FreeEK', fmt='(f17.6)')
          call cmlEndPropertyList(mainXML)
        endif
      endif
      ! On all SCF steps, print out the current energy (format depending on type of run)
      if (.not.final) then
        ! Print total energy and density matrix error .........................
        if (cml_p) then
          call cmlStartPropertyList(mainXML, title='SCF Cycle')
          ! Eharrs is always output
          call cmlAddProperty(xf=mainXML, property=Eharrs/eV, 
     .         units="eV", dictRef="siesta:Eharrs", fmt="(f14.7)")
        endif
        ! This chain of if statements determines which properties are output.
        if (harrisfun) then 
          write(6,"(/a,f14.6,/)") 'siesta: Eharris(eV) = ', Eharrs/eV
            ! No need for further cml output
        elseif (isolve==0) then
          if (cml_p) 
     .         call cmlAddProperty(xf=mainXML, property=FreeE/eV, 
     .         units="eV", dictRef="siesta:FreeE",  fmt="(f14.7)")
          if (fixspin) then
            if (cml_p) then
              call cmlAddProperty(xf=mainXML, property=Etot/eV, 
     .             units="eV", dictRef="siesta:Etot",   fmt="(f14.7)")
              call cmlAddProperty(xf=mainXML, property=FreeE/eV, 
     .             units="eV", dictRef="siesta:FreeE", fmt="(f14.7)")
              call cmlAddProperty(xf=mainXML, property=dDmax/eV, 
     .             units="eV", dictRef="siesta:dDmax",  fmt="(f14.7)")
            endif
            if ((iscf .eq. 1).or.muldeb)
     .           write(6,'(/,a12,3a14,a8,a7,a11)')
     .           'siesta: iscf', '   Eharris(eV)', 
     .           '      E_KS(eV)', '   FreeEng(eV)', 
     .           '   dDmax', '  Ef_up', '  Ef_dn(eV)'
            write(6,'(a8,i4,3f14.4,f8.4,2f9.4)')
     .           'siesta: ',iscf, Eharrs/eV, Etot/eV, FreeE/eV, dDmax, 
     .           (Efs(i)/eV,i=1,2)
            if (cml_p) then
              call cmlAddProperty(xf=mainXML, property=Efs(1)/eV, 
     .             units="eV", dictRef="siesta:Efs",    fmt="(f14.7)")
              call cmlAddProperty(xf=mainXML, property=Efs(2)/eV, 
     .             units="eV", dictRef="siesta:Efs",    fmt="(f14.7)")
            endif
          else !fixspin
            if ((iscf .eq. 1).or.muldeb)
     .           write(6,'(/,a12,3a14,2a8)')
     .           'siesta: iscf', '   Eharris(eV)', 
     .           '      E_KS(eV)', '   FreeEng(eV)', 
     .           '   dDmax', '  Ef(eV)'
            write(6,'(a8,i4,3f14.4,2f8.4)')
     .           'siesta: ',iscf, Eharrs/eV, Etot/eV, FreeE/eV, 
     .           dDmax, Ef/eV
            if (cml_p) then
              call cmlAddProperty(xf=mainXML,property=Ef/eV, 
     .             units="eV", dictRef="siesta:Ef", fmt="(f14.7)")
            endif !cml_p
          endif !fixspin
        elseif (isolve==1) then
          write(6,'(/,a15,i4)') 'siesta: iscf = ',iscf
          write(6,'(a14,f15.4,a13,f15.4,a10,f7.4/)') 
     .         'Eharris(eV) = ',Eharrs/eV,
     .         '  E_KS(eV) = ',Etot/eV,'  dDmax = ',dDmax
          if (cml_p) then
            call cmlAddProperty(xf=mainXML, property=Etot/eV, 
     .           units="eV", dictRef="siesta:Etot",   fmt="(f14.7)")
            call cmlAddProperty(xf=mainXML, property=dDmax/eV, 
     .           units="eV", dictRef="siesta:dDmax",  fmt="(f14.7)")
          endif
        endif !harrisfun/isolve
        
        if (cml_p) then
          call cmlEndPropertyList(mainXML)
        endif
        
        if (signat) call sig(io_sig,'Harris in scf',Eharrs/eV)
        if (signat) call sig(io_sig,'Etot in scf',Etot/eV)
        if (signat) call sig(io_sig,'dDmax in scf',dDmax)

      else !final
      ! Print out additional information in finalization.
        ! Signature records
        if (signat) then
          call sig(io_sig,'siesta: Eions   =', Eions/eV)
          call sig(io_sig,'siesta: Ena     =', Ena/eV)
          call sig(io_sig,'siesta: Ekin    =', Ekin/eV)
          call sig(io_sig,'siesta: Enl     =', Enl/eV)
          call sig(io_sig,'siesta: DEna    =', DEna/eV)
          call sig(io_sig,'siesta: DUscf   =', DUscf/eV)
          call sig(io_sig,'siesta: DUext   =', DUext/eV)
          call sig(io_sig,'siesta: Exc     =', Exc/eV)
          call sig(io_sig,'siesta: eta*DQ  =', Ecorrec/eV)
          call sig(io_sig,'siesta: Emadel  =', Emad/eV)
          call sig(io_sig,'siesta: Ekinion =', Ekinion/eV)
          call sig(io_sig,'siesta: Eharris =', (Eharrs1+Ekinion)/eV)
          call sig(io_sig,'siesta: Etot    =', (Etot+Ekinion)/eV)
          call sig(io_sig,'siesta: FreeEng =', (FreeE+Ekinion)/eV)
        endif

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
        if (cml_p) then
          call cmlStartPropertyList(xf=mainXML, title='Final Energy')
          call cmlAddProperty(xf=mainXML, property=Ekin/eV,  units='eV', 
     .         dictref='siesta:Ekin',  fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=Uscf/eV,  units='eV', 
     .         dictref='siesta:Uscf',  fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=DUext/eV, units='eV', 
     .         dictref='siesta:DUext', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=Exc/eV,   units='eV', 
     .         dictref='siesta:Exc',   fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, 
     .         property=(Enascf+Enl+DUscf-Uscf-Uatm)/eV, units='eV', 
     .         dictref='siesta:I-e', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, 
     .         property=(Ena+Uatm-Enaatm-Eions)/eV,      units='eV', 
     .         dictref='siesta:I-I', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML,property=Ekinion/eV,units='eV', 
     .         dictref='siesta:Ekinion', fmt='(f17.6)')
          call cmlAddProperty(xf=mainXML, property=(Etot+Ekinion)/eV,
     .         units='eV', dictref='siesta:Etot', fmt='(f17.6)')
          call cmlEndPropertyList(mainXML)
        endif !cml_p
      endif !final

      end subroutine siesta_write_energies
      
      end program siesta
