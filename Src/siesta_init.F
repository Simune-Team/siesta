! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      MODULE m_siesta_init
      private
      public :: siesta_init

      CONTAINS

      subroutine siesta_init()
      use Kpoint_grid,        only: setup_Kpoint_grid, gamma_scf, maxk
      USE Kpoint_pdos,        only: gamma_pdos
      use Band,               only: gamma_bands, setup_bands
      use m_ksvinit,          only: gamma_polarization,
     &                              estimate_pol_kpoints
      use m_struct_init,      only: struct_init
      use siesta_options
      use timer_options,      only: use_tree_timer
      use sparse_matrices
      use siesta_geom
      use atomlist,           only: no_u, rmaxkb, amass, lasto, qtot,
     &                              iza, rmaxo, zvaltot, superc,
     &                              initatomlists, no_l
      use fdf
      use sys,                only: die, bye
#ifdef BSC_CELLXC
      use bsc_xcmod,          only: setXC
#endif /* BSC_CELLXC */
      use molecularmechanics, only: inittwobody
      use metaforce,          only: initmeta
      use m_mpi_utils,        only: broadcast
      use alloc,              only: re_alloc, alloc_report
      use parallelsubs,       only: getnodeorbs
      use m_iostruct,         only: write_struct, read_struct
      use zmatrix,            only: lUseZmatrix
      use zmatrix,            only: write_canonical_ucell_and_Zmatrix
      use m_check_supercell,  only: check_sc_factors
      use files,              only: slabel
      use siesta_cmlsubs,     only: siesta_cml_init
      use m_timestamp,        only: timestamp
      use m_wallclock,        only: wallclock
      use parallel,           only: Node, Nodes, IOnode
      use densematrix,        only: Haux, Saux, Psi
      use m_energies
      use m_steps
      use m_spin
      use m_rmaxh
      use m_forces 
      use m_eo
      use m_gamma
      use m_ioxv,             only: xv_file_read
      use m_projected_DOS,    only: init_projected_DOS
      use writewave,          only: gamma_wavefunctions,
     &                              setup_wf_kpoints
#ifdef DEBUG_XC
      USE siestaXC, only: setDebugOutputUnit
#endif /* DEBUG_XC */
      USE m_timer, only: timer_report ! Sets options for report of CPU times
      use m_check_walltime
#ifdef MPI
      use mpi_siesta
#endif

#ifdef TRANSIESTA
      use m_ts_init, only : ts_init
#endif

      implicit none
#ifdef MPI
      logical :: initialized
      integer :: MPIerror  ! Return error code in MPI routines
#endif
      real(dp):: veclen       ! Length of a unit-cell vector
      integer :: i, is, ispin
      integer :: neigmin      ! Min. number of eigenstates (per k point)
      integer :: ns           ! Number of species
      logical :: user_wants_supercell
      real(dp) :: walltime_m, walltime_s

      external :: read_options
      external :: reset_messages_file
      
!------------------------------------------------------------------------- BEGIN
! Initialise MPI unless siesta is running as a subroutine 
! of a driver program that may have initialized MPI already
#ifdef MPI
      call MPI_Initialized( initialized, MPIerror )
      if (.not.initialized) then
        call MPI_Init( MPIerror )
#ifdef _TRACE_
        call MPItrace_shutdown( )
#endif
      endif ! (.not.initialized)
      call MPI_Comm_Rank( MPI_Comm_World, Node, MPIerror )
      call MPI_Comm_Size( MPI_Comm_World, Nodes, MPIerror )

      call reset_messages_file()
      
      if (.not. fdf_parallel()) then
        call die('siesta_init: ERROR: FDF module ' //
     &           'doesn''t have parallel support')
      endif
#else
      call reset_messages_file()
#endif
      if (Node==0) then
         call system("rm -f 0_NORMAL_EXIT")
      endif
#ifdef DEBUG
!     Generates a debug file for every process to track the execution.
!     The file is called debug.$(PID) Where PID is the process number.
!     It also works in secuencial mode.
      call debugMpiOn( )
      call write_debug( 'PRE siesta_init' )
#endif
      IOnode = (Node .eq. 0)

!     Print version information
      if (IOnode) then
        call prversion()
#ifdef MPI
        if (Nodes.gt.1) then
          write(6,'(/,a,i4,a)')
     &        '* Running on ', Nodes, ' nodes in parallel'
        else
          write(6,'(/,a,i4,a)')
     &        '* Running in serial mode with MPI'
        endif
#else
        write(6,'(/,a,i4,a)')
     &        '* Running in serial mode'
#endif
         call timestamp('Start of run')
         call wallclock('Start of run')
      endif


!     Nullify arrays
      nullify(Haux,Saux,psi)
      nullify(Dscf,Dold,Eold,Escf)
      nullify(listh)
      nullify(H,S)
      nullify(xijo)

!     Initialize some variables
      DUext   = 0.0_dp
      Eharrs  = 0.0_dp
      Eions   = 0.0_dp
      Ekinion = 0.0_dp
      Emad    = 0.0_dp
      Emeta   = 0.0_dp
      Emm     = 0.0_dp
      Entropy = 0.0_dp
      FreeE   = 0.0_dp

#ifdef DEBUG_XC
! Set output file unit for debug info
      call setDebugOutputUnit()
#endif /* DEBUG_XC */

!     Initialise read .................................................
      call reinit( sname )

!     Initialize CML (relies on reinit)
      call siesta_cml_init( )

! Set timer report file and threshold .................................
      threshold = fdf_get('timer_report_threshold', 0._dp)
      call timer_report( file=trim(slabel)//'.times', 
     .                   threshold=threshold )
      ! This variable lives in siesta_options
      use_tree_timer = fdf_get('UseTreeTimer', .false.)

!     Start time counter
!     Note new placement of this first use, so that
!     initialization is done after the setup of relevant variables

      call timer( 'siesta', 0 )
      call timer( 'siesta', 1 )
      call timer( 'Setup', 1 )

! Set allocation report level .........................................
! variables level and threshold imported from module siesta_options
      level = fdf_get('alloc_report_level', 0)
      threshold = fdf_get('alloc_report_threshold', 0._dp)
      call alloc_report( level=level, file=trim(slabel)//'.alloc',
     .                   threshold=threshold,
     &                   printNow=.false. )

!     Initialise exchange-correlation functional information
#ifndef BSC_CELLXC
      call read_xc_info()
#else /* BSC_CELLXC */
      call setXC( )
#endif /* BSC_CELLXC */

!     Initialise force field component
      call inittwobody( )

!     Initialize pseudopotentials and atomic orbitals
      if (IOnode) call initatom( ns )
      call broadcast( ns )

      call broadcast_basis( )

      call register_rfs()
!     call test_register()

      atmonly = fdf_get( 'Atom-Setup-Only', .false. )
      if (atmonly) call bye( 'End of atom setup' )

!     Read geometry
      call struct_init( )      ! Sets na_u, isa, ucell

!     Initialize atom lists
      call initatomlists( )    ! Sets iza

!     early exit if only checking the structure
      struct_only = fdf_get( 'Output-Structure-Only', .false. )
      if (IONode) then
         call write_struct( ucell, na_u, isa, iza, xa )
         if (fdf_boolean('WriteCoorXmol',.false.)) then
            call coxmol(iza, xa, na_u)
         endif
         if (lUseZmatrix) then
            call write_canonical_ucell_and_Zmatrix(
     &                        filename="OUT.UCELL.ZMATRIX.INITIAL")
         endif
      endif
      if (struct_only) then
         call bye('End of structure processing')
      endif

!     Walltime control (numbers in seconds)
      if ( fdf_isphysical("MaxWalltime") ) then
         walltime_m = fdf_get("MaxWalltime", huge(1._dp), 's')
      else
         walltime_m = fdf_get("MaxWalltime", huge(1._dp))
      end if
      ! Period for clean-up operations
      if ( fdf_isphysical("MaxWalltime.Slack") ) then
         walltime_s = fdf_get("MaxWalltime.Slack", 5._dp, 's')
      else
         walltime_s = fdf_get("MaxWalltime.Slack", 5._dp)
      end if
      ! Note that the slack detracts from net available time
      walltime_warning = walltime_m - walltime_s
      
!     End of Initial Structure Processing
      if (Node.eq.0) then
        write(6,'(/,a,20("*"),a,28("*"))')
     &    'siesta: ', ' Simulation parameters '
        write(6,'(a)')  'siesta:'
        write(6,'(a)')  'siesta: The following are some of the '//
     &                           'parameters of the simulation.'
        write(6,'(a)')  'siesta: A complete list of the parameters '//
     &                           'used, including default values,'
        write(6,'(a,a)')'siesta: can be found in file out.fdf'
        write(6,'(a)')  'siesta:'
      endif

!     Allocate other arrays based on read sizes
      nullify(fa,cfa)
      call re_alloc( fa, 1, 3, 1, na_u, 'fa', 'siesta_init' )
      call re_alloc( cfa, 1, 3, 1, na_u, 'cfa', 'siesta_init' )

      call init_spin( )

!     Read simulation data
      call read_options( na_u, ns, nspin )

      qtot = qtot - charnet     ! qtot set in initatomlists
                                ! charnet set in redata
      if (IOnode) then
        write(6,fmt="(a,f12.6)") 'Total number of electrons: ', qtot
        write(6,fmt="(a,f12.6)") 'Total ionic charge: ', zvaltot
      endif
!
!     Warn the user: if not doing a direct optimization, the Zmatrix
!     coordinates are no longer updated. Only coordinates are treated.
!
      if (lUseZmatrix) then
         if (idyn .ne. 0) then
            write(6,"(a)")
     &      " *** WARNING: Zmatrix form will be used only for input !!"
            write(0,"(a)")
     &      " *** WARNING: Zmatrix form will be used only for input !!"
         endif
      endif

! Calculate spin populations for fixed spin case...
      if (fixspin) then
        if (nspin .ne. 2)
     &    call die( 'siesta: ERROR: ' //
     &              'You can only fix the spin of the system' //
     &              ' for collinear spin polarized calculations.' )
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

!     Find maximum interaction range
      if (negl) then
        rmaxh = 2.0_dp*rmaxo
      else
        rmaxh = 2.0_dp*rmaxo + 2.0_dp*rmaxkb
      endif

!     Madelung correction for charged systems
      if (charnet .ne. 0.0_dp) then
        call madelung(ucell, shape, charnet, Emad)
      endif

!     Parallel initialisation
      call initparallel( no_u, na_u, lasto, xa, ucell, rmaxh )
!     if (IOnode) call show_distribution( )

!     Find number of locally stored orbitals and allocated related arrays
      call GetNodeOrbs(no_u,Node,Nodes,no_l)

!     Initialise arrays
      nullify( listhptr, numh )
      call re_alloc( listhptr, 1, max(1,no_l), 'listhptr', 'sparseMat' )
      call re_alloc( numh,     1, max(1,no_l), 'numh',     'sparseMat' )
      listhptr    = 0
      numh        = 0

!     Get number of eigenstates that need to be calculated
      neigwanted = fdf_get('NumberOfEigenStates',no_u)

!     Check number of eigenstates - cannot be larger than number of
!     basis functions or smaller than number of occupied states + 1
!     so that the Fermi level can be estimated
      do is = 1,nspin
        neigmin = nint(qs(is)/real(3 - min(nspin,2), kind=dp)) + 1
        neigwanted = max(neigwanted,neigmin)
      enddo
      neigwanted = min(neigwanted,no_u)

!     Find k-grid for Brillouin zone integration 
!     NOTE: We need to know whether gamma is .true. or
!     not early, in order to decide whether to use an 
!     auxiliary supercell for the calculation of matrix elements.
      call setup_Kpoint_grid( ucell )
      gamma = gamma_scf

!     Call initialisation of PDOS here since we need to check if 
!     the auxiliary supercell is needed for a non-gamma calculation
      call init_projected_DOS( )
      if (do_pdos) then
        gamma = gamma .and. gamma_pdos
      endif

#ifdef TRANSIESTA
      call ts_init(nspin,ucell,na_u,xa)
#endif

      nullify(eo,qo)
      call re_alloc(eo, 1, no_u, 1, nspin, 1, maxk, 'eo', 'siesta_init')
      call re_alloc(qo, 1, no_u, 1, nspin, 1, maxk, 'qo', 'siesta_init')

      call setup_bands( )
      gamma = gamma .and. gamma_bands

      call setup_wf_kpoints( )
      gamma = gamma .and. gamma_wavefunctions

      call estimate_pol_kpoints( ucell )
      gamma = gamma .and. gamma_polarization
!      print *, "gamma, gamma_pol: ", gamma, gamma_polarization
!
!     User can request that the calculation is done with an explicit
!     auxiliary supercell and hermitian version, even if using only
!     the gamma point
!
      user_wants_supercell = fdf_get( 'ForceAuxCell', .false. )
      if (user_wants_supercell) gamma = .false.

!     Find required supercell
!     2*rmaxh is used to guarantee that two given orbitals in the
!     supercell can only overlap once
      if (gamma) then
        nsc(1:3) = 1
      else
        do i=1,3
          veclen = sqrt(ucell(1,i)**2+ucell(2,i)**2+ucell(3,i)**2)
          nsc(i) = ceiling( 2 * rmaxh / veclen )
        enddo
        if (.not. naiveauxcell)
     &    call check_sc_factors( ucell, nsc, 2*rmaxh )
      endif

      mscell = 0.0_dp
      do i = 1, 3
        mscell(i,i) = nsc(i)
        nscold(i)   = nsc(i)
      enddo

!     Find auxiliary supercell (required only for k sampling)
      call superc( ucell, scell, nsc )

!     Initialise metadynamic forces if required
      call initmeta( )

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
        call die( "'PHONON' support is deprecated" )
      else if (idyn .eq. 8) then
        inicoor = 0
        fincoor = huge(1)
      else
        call die( 'siesta: wrong idyn' )
      endif

!     Build initial velocities according to Maxwell-Bolzmann distribution....
      if (idyn .ne. 0 .and. idyn .ne. 6 .and. (.not. xv_file_read)) 
     &    call vmb( na_u, tempinit, amass, xa, isa, va )

      istp = 0

!     Output memory use before main loop
!!      call printmemory( 6, 0 )

!     Initialization now complete. Flush stdout.
      if (ionode) call pxfflush( 6 )

      call timer( 'Setup', 2 )
#ifdef DEBUG
      call write_debug( 'POS siesta_init' )
#endif
!--------------------------------------------------------------------------- END
      END subroutine siesta_init

      END MODULE m_siesta_init
