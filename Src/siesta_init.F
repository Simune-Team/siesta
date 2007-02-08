      MODULE m_siesta_init
      private
      public :: siesta_init

      CONTAINS

      subroutine siesta_init()
      USE SIESTA_TODO
      USE Kpoint_grid, only: setup_Kpoint_grid
      USE siesta_options
      implicit none


      external :: automatic_cell,
     .  bonds, cgvc, cgvc_zmatrix, fixed,
     .  dhscf, diagon, dnaefs, extrapol, initatom,
     .  iodm, iozm,
     .  kinefsm, mulliken, naefs, neighb, !kgrid, kgridinit,
     .  pulayx, 
     .  reinit, shaper, spnvec, 
     .  timer, xijorb, memory,
     .  ioeig, iofa, iokp, iomd, prversion, eggbox


!------------------------------------------------------------------------- BEGIN
! Initialise MPI and set processor number
#ifdef MPI
      call MPI_Init( MPIerror )
      call MPI_Comm_Rank( MPI_Comm_World, Node, MPIerror )
      call MPI_Comm_Size( MPI_Comm_World, Nodes, MPIerror )
#endif

      IOnode = (Node .eq. 0)

! Print version information ...........................................
      if (IOnode) then
        call prversion()
#ifdef MPI
        if (Nodes.gt.1) then
          write(6,'(/,a,i4,a)')
     .        '* Running on ', Nodes, ' nodes in parallel'
        else
          write(6,'(/,a,i4,a)')
     .        '* Running in serial mode with MPI'
        endif
#else
        write(6,'(/,a,i4,a)')
     .        '* Running in serial mode'
#endif
         call timestamp('Start of run')
         call wallclock('Start of run')
      endif
! ..................

! Start time counter ..................................................
      call timer( 'siesta', 0 )
      call timer( 'siesta', 1 )
      call timer( 'Setup', 1 )

! Nullify arrays
      nullify(Haux,Saux,psi)
      nullify(kpol,wgthpol)
      nullify(listhptr,listhptrold)
      nullify(numh,numhold)

! Initialize some variables
      DUext = 0.0_dp
      Eharrs = 0.0_dp
      Eharrs1 = 0.0_dp
      Eions = 0.0_dp
      Ekinion = 0.0_dp
      Elast = 0.0_dp
      Emad = 0.0_dp
      Emeta = 0.0_dp
      Emm = 0.0_dp
      Entrop = 0.0_dp
      Entropy = 0.0_dp
      FreeE = 0.0_dp

! Initialise read .....................................................
      call reinit(sname)

      call siesta_cml_init() ! Initialize CML (relies on reinit)

! Set allocation report level .........................................
      call fdf_global_get(level, 'alloc_report_level', 0 )
      call alloc_report( level=level, file=trim(slabel)//'.alloc' )
! ..................

! Initialise exchange-correlation functional information
      call setXC()

! Initialise force field component
      call inittwobody()

! Initialize pseudopotentials and atomic orbitals
      if (IOnode) call initatom(ns)
      call broadcast_basis()

      call fdf_global_get(atmonly,'Atom-Setup-Only',.false.)
      if (atmonly) call bye("End of atom setup")

! Read geometry
      call struct_init()      ! Sets na_u, isa, ucell

! Initialize atom lists

      call initatomlists()    ! Sets iza

! early exit if only checking the structure

      call fdf_global_get(struct_only,'Output-Structure-Only',.false.)
      if (struct_only) then
         if (IONode) call write_struct( ucell, na_u, isa, iza, xa )
         call bye("End of structure processing")
      endif
      
!  --   End of Initial Structure Processing --

      if (Node.eq.0) then
        write(6,'(/,a,20("*"),a,28("*"))')
     .    'siesta: ', ' Simulation parameters '
        write(6,'(a)')  'siesta:'
        write(6,'(a)')  'siesta: The following are some of the '//
     .                           'parameters of the simulation.'
        write(6,'(a)')  'siesta: A complete list of the parameters '//
     .                           'used, including default values,'
        write(6,'(a,a)')'siesta: can be found in file out.fdf'
        write(6,'(a)')  'siesta:'
      endif

! Allocate other arrays based on read sizes ................
!
      allocate(cfa(3,na_u))
      call memory('A','D',3*na_u,'siesta')
      allocate(fa(3,na_u))
      call memory('A','D',3*na_u,'siesta')
      allocate(fal(3,na_u))
      call memory('A','D',3*na_u,'siesta')

! Initialise those arrays that must be pre-initialised
      fal(1:3,1:na_u) = 0.0_dp
      stress(1:3,1:3) = 0.0_dp

      if (IOnode) call spin_init(nspin)
      call broadcast(nspin)

      allocate(efs(nspin))
      call memory('A','D',nspin,'siesta')
      allocate(qs(nspin))
      call memory('A','D',nspin,'siesta')
      allocate(polR(3,nspin))
      call memory('A','D',3*nspin,'siesta')
      allocate(polxyz(3,nspin))
      call memory('A','D',3*nspin,'siesta')

! Read simulation data ................................................
      call read_options( na_u, ns, nspin )


      qtot = qtot - charnet     ! qtot set in initatomlists
                                ! charnet set in redata
      write(6,fmt="(a,f12.6)") "Total number of electrons: ", qtot
      write(6,fmt="(a,f12.6)") "Total ionic charge: ", zvaltot


! Calculate spin populations for fixed spin case...
      if (fixspin) then
        if (nspin .ne. 2)
     $   call die('siesta: ERROR: ' //
     $        'You can only fix the spin of the system' //
     $        ' for collinear spin polarized calculations.')
        do i = 1,2
          qs(i) = (qtot + (3-2*i)*ts) / 2.0_dp
        enddo
      else
        qs(1:nspin) = 0.0_dp
        if (nspin .le. 2) then
          do ispin = 1,nspin
            qs(ispin) = qtot/nspin
          enddo
        endif
      endif
! ..................

! Find maximum interaction range ......................................
      if (negl) then
        rmaxh = 2.0_dp*rmaxo
      else
        rmaxh = 2.0_dp*rmaxo + 2.0_dp*rmaxkb
      endif
! ......................

! Madelung correction for charged systems .............................
      if (charnet .ne. 0.0_dp) then
        call madelung(ucell, shape, charnet, Emad)
      endif

! Parallel initialisation
      call initparallel(no_u,na_u,lasto,xa,ucell,rmaxh,rcoor,isolve)

! Find number of locally stored orbitals and allocated related arrays
      call GetNodeOrbs(no_u,Node,Nodes,no_l)

! Initialise arrays
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

! Get number of eigenstates that need to be calculated

      call fdf_global_get(neigwanted,'NumberOfEigenStates',no_u)


! Check number of eigenstates - cannot be larger than number of
! basis functions or smaller than number of occupied states + 1
! so that the Fermi level can be estimated
      do is = 1,nspin
        neigmin = nint(qs(is)/real(3 - min(nspin,2), kind=dp)) + 1
        neigwanted = max(neigwanted,neigmin)
      enddo
      neigwanted = min(neigwanted,no_u)

      ! Find k-grid for Brillouin zone integration 
      ! Call kgridinit to get the dimensions and then kgrid to actually
      ! initialise the data 
      call setup_Kpoint_grid( ucell )

      allocate(eo(no_u,nspin,maxk))
      call memory('A','D',no_u*nspin*maxk,'siesta')
      allocate(qo(no_u,nspin,maxk))
      call memory('A','D',no_u*nspin*maxk,'siesta')



! Find number of k-points for wavefunction printout ....................
      nwk = 0
      maxwk = 1
      allocate(wfk(3,maxwk))
      call memory('A','D',3*maxwk,'siesta')

      call initwave( maxwk, no_u, nwk, wfk, overflow )

      if (overflow) then
        if (nwk .gt. maxwk) then
! If there wasn't enough space to store bands on first call correct
! the dimensions and repeat the initialisation
          maxwk = max(nwk,1)
          call memory('D','D',size(wfk),'siesta')
          deallocate(wfk)
          allocate(wfk(3,maxwk))
          call memory('A','D',3*maxwk,'siesta')
        endif
        nwk = 0
        call initwave( maxwk, no_u, nwk, wfk, overflow )
        if (overflow)
     $   call die('siesta: ERROR: Unsuccessful initialization of' //
     $            ' list of wavefunctions to print')
      endif
! ......................

! Find the grid for the calculation of the polarization..............
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

! Find if only gamma point is used ....................................
      if (nkpnt.eq.1 .and. abs(kpoint(1,1)).lt.tiny .and.
     .                     abs(kpoint(2,1)).lt.tiny .and.
     .                     abs(kpoint(3,1)).lt.tiny) then
        gamma = .true.
      else
        gamma = .false.
      endif
      ! eliminado por conflicto con la diagonalizacion
      ! codigo de bandas movido a siesta_analysis
!      if (nbk .gt. 0) gamma = .false.

      if (nwk .gt. 1) gamma = .false.
      if (nwk .eq. 1) then
        if (abs(wfk(1,1)).gt.tiny .and.
     .      abs(wfk(2,1)).gt.tiny .and.
     .      abs(wfk(3,1)).gt.tiny) then
          gamma = .false.
        endif
      endif
      if (nkpol.gt.0) gamma = .false.
! ....................

! Print k-point information (even if gamma point)
      if (IOnode) call siesta_write_k_points()

! Find required supercell
! 2*rmaxh is used to guarantee that two given orbitals in the
! supercell can only overlap once

      if (gamma) then
         nsc(1:3) = 1
      else
         do i=1,3
            veclen = sqrt(ucell(1,i)**2+ucell(2,i)**2+ucell(3,i)**2)
            nsc(i) = ceiling( 2 * rmaxh / veclen )
         enddo
         if (.not. naiveauxcell)
     $        call check_sc_factors(ucell,nsc,2*rmaxh)
      endif

      mscell = 0.0_dp
      do i = 1, 3
         mscell(i,i) = nsc(i)
         nscold(i) = nsc(i)
      enddo


! Find auxiliary supercell (required only for k sampling) ............
      call superc( ucell, scell, nsc)

! Initialise metadynamic forces if required
      call initmeta()

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
        call phonon_setup()
        inicoor = 1
        fincoor = phonon_num_disps
      else if (idyn .eq. 8) then
        inicoor = 0
        fincoor = ihuge
      else
         call die('siesta: wrong idyn')
      endif

! Build initial velocities according to Maxwell-Bolzmann distribution....
      if (idyn .ne. 0 .and. idyn .ne. 6 .and. (.not. foundxv)) 
     .    call vmb(na_u,tempinit,amass,xa,isa,va)
! ..................

      istp = 0
      call timer( 'Setup', 2 )

! Output memory use before main loop
      call printmemory( 6, 0 )

! Initialization now complete. Flush stdout.
      if (ionode) call pxfflush(6)
!--------------------------------------------------------------------------- END

      END subroutine siesta_init

      End MODULE m_siesta_init
