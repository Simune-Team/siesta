      MODULE m_state_init

      private
      public :: state_init

      CONTAINS

      subroutine state_init( istep )
      use m_new_dm,          only: new_dm
      use Kpoint_grid,       only: setup_Kpoint_grid, maxk
      use m_proximity_check, only: proximity_check
      use siesta_options
      use sparse_matrices
      use siesta_geom
      use atomlist,          only: iphorb, iphkb, rco, rckb, indxua,
     &                             no_s, rmaxkb, rmaxo, no_u, lastkb,
     &                             lasto, superc, rmaxv, indxuo
      use atomlist,          only: no_l
      use alloc,             only: re_alloc, de_alloc, alloc_report
      use m_hsparse,         only: hsparse
      use m_overlap,         only: overlap
      use m_check_supercell, only: check_sc_factors
      use siesta_cml
      use siesta_cmlsubs,    only: siesta_cml_init
      use zmatrix,           only: lUseZmatrix, write_zmatrix
      use m_energies,        only: Emad
      use write_subs
      use m_steps
      use parallel,          only: IOnode, node, nodes
      use m_spin,            only: nspin
      use m_rmaxh
      use m_pulay,           only: init_pulay_arrays
      use m_eo
      use m_gamma
      use files,             only: slabel
#ifdef MPI
      use m_mpi_utils,       only: globalize_max
#endif
#ifdef CDF
      use iodm_netcdf, only: setup_dm_netcdf_file
      use iodmhs_netcdf, only: setup_dmhs_netcdf_file
#endif

      implicit none

      integer            :: istep
      real(dp)           :: veclen      ! Length of a unit-cell vector
      real(dp)           :: rmax
      logical            :: cell_can_change
      integer            :: i, ix, iadispl, ixdispl
      logical            :: auxchanged   ! Auxiliary supercell changed?
      external           ::  madelung, timer
      real(dp), external :: volcel
!------------------------------------------------------------------------- BEGIN
#if defined(MPI_TRACE) && defined(MPI)
      call mpitrace_eventandcounters( 1000, 11 )
#endif
      call timer( 'IterMD', 1 )

      istp = istp + 1

      if (IOnode) then
        write(6,'(/2a)') 'siesta:                 ',
     &                    '=============================='
        select case (idyn)
        case (0)
          write(6,'(28(" "),a,i6)') 'Begin CG move = ',istep
          if (cml_p) call cmlStartStep(mainXML, type='CG', index=istp)
          ! Print Z-matrix coordinates
          if (lUseZmatrix) then
             call write_Zmatrix()
          endif
        case (1:5)
          write(6,'(28(" "),a,i6)') 'Begin MD step = ',istep
          if (cml_p) call cmlStartStep(mainXML, type='MD', index=istep)
        case (6)
          write(6,'(28(" "),a,i6)') 'Begin FC step = ',istep
          if (cml_p) call cmlStartStep(mainXML, type='FC', index=istep)

          if (istep .eq. 0) then
            write(6,'(28(" "),a)') 'Undisplaced coordinates'
          else
            iadispl = (istep-mod(istep-1,6))/6+ia1
            write(6,'(28(" "),a,i6)') 'displace atom   ', iadispl
            ix = mod(istep-1,6)+1
            ixdispl = (ix - mod(ix-1,2) +1)/2
            write(6,'(28(" "),a,i6)') 'in direction    ', ixdispl
            write(6,'(28(" "),a,f8.4,a)') 'by       ', dx, ' Bohr'
          endif

        case (8)
          write(6,'(28(" "),a,i6)') 'Begin Server step = ',istep
          if (cml_p) call cmlStartStep(mainXML, type='FS', index=istep)
          ! Find interatomic distances (output in file BONDS)
          call bonds( ucell, na_u, isa, xa,
     $                rmax_bonds, trim(slabel)// ".BONDS" )

        end select
        write(6,'(2a)') '                        ',
     &                    '=============================='


C       Print atomic coordinates
        call outcoor( ucell, xa, na_u, ' ', writec )
        call siesta_write_positions()
      endif

C     Actualize things if variable cell 
      auxchanged = .false.
      cell_can_change = ( varcel .or.
     &                    (idyn .eq. 8)  ! Force/stress evaluation
     &                  )
      if (change_kgrid_in_md) then
         cell_can_change = cell_can_change .or.
     &                     (idyn .eq. 3)   .or. ! Parrinello-Rahman
     &                     (idyn .eq. 4)   .or. ! Nose-Parrinello-Rahman
     &                     (idyn .eq. 5)        ! Anneal
      endif

      if ( cell_can_change .and.
     &     (istep.ne.inicoor) .and. (.not.gamma) ) then

C       Will print k-points also
        call setup_Kpoint_grid( ucell )

        call re_alloc( eo, 1, no_u, 1, nspin, 1, maxk, 'eo',
     &                 'state_init')
        call re_alloc( qo, 1, no_u, 1, nspin, 1, maxk, 'qo',
     &                 'state_init' )

C       Find required supercell
        if (gamma) then
           nsc(1:3) = 1
        else if (fixauxcell) then
           nsc(1:3) = nscold(1:3)
        else
           do i=1,3
              veclen = sqrt(ucell(1,i)**2+ucell(2,i)**2+ucell(3,i)**2)
              nsc(i) = ceiling( 2 * rmaxh / veclen )
           enddo
           if (.not. naiveauxcell)
     &         call check_sc_factors(ucell,nsc,2*rmaxh)
        endif

        mscell = 0.0_dp
        do i = 1, 3
           mscell(i,i) = nsc(i)
           if (nsc(i).ne.nscold(i)) auxchanged = .true.
           nscold(i) = nsc(i)
        enddo

C       Madelung correction for charged systems 
        if (charnet .ne. 0.0_dp) then
          call madelung(ucell, shape, charnet, Emad)
        endif

      endif
C     End variable cell actualization

C     Auxiliary supercell
C     Do not move from here, as the coordinates might have changed
C     even if not the unit cell
      call superc(ucell, scell, nsc)

C     Print unit cell and compute cell volume
C     Possible BUG: 
C     Note that this volume is later used in write_subs and the md output
C     routines, even if the cell later changes.
      if (IOnode) call outcell(ucell)
      volume_of_some_cell = volcel(ucell)

C     Use largest possible range in program, except hsparse...
C     2 * rmaxv: Vna overlap
C     2 * rmaxo: orbital overlap
C     rmaxo + rmaxkb: Non-local KB action
C     2.0_dp * (rmaxo+rmaxkb) : Orbital interaction through KB projectors
      rmax = max( 2._dp*rmaxv, 2._dp*rmaxo, rmaxo+rmaxkb)

C     AG: Why not use the largest range and be done with it?
      if (.not. negl) then
        rmax = max(rmax, 2.0_dp * (rmaxo+rmaxkb) )
      endif

C     Check if any two atoms are unreasonably close
      call proximity_check(rmax)

      if (istp==1) then 
        nullify( Dscf, Dscfsave, Dold, Eold, Escf )
        nullify( listhold, listh )
        nullify( H, H0, S )
        nullify( xijo )
        nullify( listhptrold, numhold )
      endif

C     List of nonzero Hamiltonian matrix elements
C     and, if applicable,  vectors between orbital centers
C     Listh (and xijo if .not. gamma) are allocated inside hsparse
C     nh is set by the routine also
   
      nh = maxnh
#ifdef DEBUG
      call write_debug( '    PRE hsparse' )
#endif
      call hsparse( negl, scell, nsc, na_s, isa, xa, lasto,
     &              lastkb, iphorb, iphKB, nh, gamma )
#ifdef DEBUG
      call write_debug( '    POS hsparse' )
#endif
      maxnh = nh

      call re_alloc( listhptrold, 1, no_l, 'listhptrold', 'sparseMat')
      call re_alloc( numhold, 1, no_l, 'numhold', 'sparseMat' )
!!    Do Not initialize here. These are managed in extrapol
!!      listhptrold = 0
!!      numhold     = 0

C     OBSOLETE: 
C       safety margin for when atoms move 
C       if (nh > maxnh)  maxnh = 1.05 * nh + 40 

C     Use re_alloc regardless of whether nh is greater than 
C     the old maxnh or not, BUT DO NOT SHRINK, as some of
C     these arrays are later used in extrapol.
      call re_alloc( Dscf, 1, maxnh, 1, nspin, 'Dscf',
     &               'sparseMat', .true., .false.)
      call re_alloc( Dscfsave, 1, maxnh, 1, nspin, 'Dscfsave',
     &               'sparseMat', .true., .false. )
      call re_alloc( listhold, 1, maxnh, 'listhold',
     &               'sparseMat', .true., .false. )

      if (istp==1) then
C       Initialise Dscfsave to avoid problems in extrapol
C       This might be superfluous, as alloc fills array with zeroes
        Dscfsave(1:maxnh,1:nspin) = 0.0_dp
        listhold(1:maxnh)         = 0
      endif

      call re_alloc( Dold, 1, maxnh, 1, nspin, 'Dold',
     &               'sparseMat', copy=.false. )
      call re_alloc( Eold, 1, maxnh, 1, nspin, 'Eold',
     &               'sparseMat', copy=.false. )
      call re_alloc( Escf, 1, maxnh, 1, nspin, 'Escf',
     &               'sparseMat', copy=.false. )

C     Allocate/reallocate storage associated with Hamiltonian/Overlap matrix
      call re_alloc( H, 1, maxnh, 1, nspin, 'H', 'sparseMat',
     &               .false., .false.)
      call re_alloc( H0, 1, maxnh, 'H0', 'sparseMat',
     &               .false., .false. )
      call re_alloc( S, 1, maxnh, 'S', 'sparseMat',
     &               .false., .false. )

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

C     Initialize density matrix
      call new_dm( auxchanged )

C     Check for size of Pulay auxiliary matrices
      call init_pulay_arrays()

C     Find overlap matrix 
#ifdef DEBUG
      call write_debug( '    PRE overlap' )
#endif
      call overlap(na_u, na_s, no_s, scell, xa, indxua, rmaxo,
     &             maxnh, lasto, iphorb, isa, 
     &             numh, listhptr, listh, S )
#ifdef DEBUG
      call write_debug( '    POS overlap' )
#endif

#if defined(MPI_TRACE) && defined(MPI)
      call mpitrace_eventandcounters( 1000, 0 )
#endif

#ifdef CDF
      if (writedm_cdf .or. writedm_cdf_history) then
         call setup_dm_netcdf_file( maxnh, no_l, nspin,
     $                              no_s, indxuo,
     $                              numh, listhptr, listh)
      endif
      if (writedmhs_cdf .or. writedmhs_cdf_history) then
         call setup_dmhs_netcdf_file( maxnh, no_l, nspin,
     $                              no_s, indxuo,
     $                              numh, listhptr, listh,
     $                              s )
      endif
#endif
!--------------------------------------------------------------------------- END
      END subroutine state_init
      END module m_state_init
