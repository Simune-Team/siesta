      MODULE m_state_init

      private
      public :: state_init

      CONTAINS

      subroutine state_init( istep )
      use Kpoint_grid,       only: setup_Kpoint_grid, maxk
      use m_io_s,            only: file_exist
      use m_new_dm,          only: new_dm
      use m_proximity_check, only: proximity_check
      use siesta_options
      use units,             only: Ang

      use sparse_matrices, only: maxnh, numh, listh, listhptr
      use sparse_matrices, only: Dscf, Dold, Escf, Hold
      use sparse_matrices, only: xijo, H, S_1D, S, H_kin_1D
      use sparse_matrices, only: H_2D, xij_2D
      use sparse_matrices, only: sparse_pattern, block_dist
      use sparse_matrices, only: H_vkb_1D, DM_2D, EDM_2D, DM_history

      use siesta_geom
      use atomlist,          only: iphorb, iphkb, indxua,
     &                             no_s, rmaxkb, rmaxo, no_u, lastkb,
     &                             lasto, superc, rmaxv, indxuo,
     &                             no_l
      use alloc,             only: re_alloc, de_alloc, alloc_report
      use m_hsparse,         only: hsparse
      use m_overlap,         only: overlap
      use m_supercell,       only: exact_sc_ag
      use siesta_cml,        only: cml_p, cmlStartStep, mainXML
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
      use m_mpi_utils,       only: globalize_or
      use m_mpi_utils,       only: globalize_max
      use domain_decom,      only: domainDecom, use_dd, use_dd_perm

      use m_sparse, only : xij_offset

#ifdef TRANSIESTA
      use m_ts_charge, only : TS_RHOCORR_METHOD, TS_RHOCORR_FERMI
      use m_ts_options, only : opt_TriMat_method
      use m_ts_options, only : onlyS
      use m_ts_options, only : TS_Analyze
      use m_ts_options, only : N_Elec, Elecs, IsVolt
      use m_ts_electype
      use m_ts_global_vars, only: TSrun, TSmode
      use sys, only : bye
      use m_ts_io, only : fname_TSHS, ts_write_tshs
      use m_ts_sparse, only : ts_sparse_init, ts_sp_uc
      use m_ts_tri_init, only : ts_tri_init, ts_tri_analyze
      use files, only: slabel, label_length
#endif
#ifdef CDF
      use iodm_netcdf, only: setup_dm_netcdf_file
      use iodmhs_netcdf, only: setup_dmhs_netcdf_file
#endif
      use class_Sparsity
      use class_dSpData1D
      use class_dSpData2D
      use class_dData2D
#ifdef TEST_IO
      use m_test_io
#endif

      implicit none

      integer            :: istep
      real(dp)           :: veclen      ! Length of a unit-cell vector
      real(dp)           :: rmax
      logical            :: cell_can_change
      integer            :: i, ix, iadispl, ixdispl
      logical            :: auxchanged   ! Auxiliary supercell changed?
      logical            :: folding, folding1
      external           ::  madelung, timer
      real(dp), external :: volcel
#ifdef TRANSIESTA
      integer                       :: ts_kscell_file(3,3) = 0
      real(dp)                      :: ts_kdispl_file(3) = 0.0
      logical                       :: ts_Gamma_file = .true.
      character(len=label_length+6) :: fname
      integer                       :: fnlength
      real(dp)                      :: dummyef=0.0, dummyqtot=0.0
#endif
      
      character(len=256)            :: oname

      type(dData2D) :: tmp_2D


!------------------------------------------------------------------- BEGIN
      call timer( 'IterMD', 1 )
#ifdef DEBUG
      call write_debug( '  PRE state_init' )
#endif
      call timer( 'STinit', 1 )

      istp = istp + 1

      if (IOnode) then
        write(6,'(/2a)') 'siesta:                 ',
     &                    '=============================='
        select case (idyn)
        case (0)
          write(6,'(28(" "),a,i6)') 'Begin CG move = ',istep
          if (cml_p) call cmlStartStep(mainXML, type='CG', index=istp)
!         Print Z-matrix coordinates
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
            write(6,'(28(" "),a,f8.4,a)') 'by       ', dx/Ang, ' Ang'
          endif

        case (8)
          write(6,'(28(" "),a,i6)') 'Begin Server step = ',istep
          if (cml_p) call cmlStartStep(mainXML, type='FS', index=istep)
!         Find interatomic distances (output in file BONDS)
          call bonds( ucell, na_u, isa, xa,
     &                rmax_bonds, trim(slabel)// ".BONDS" )

       case (9)
          if ( istep == 0 ) then
             write(6,'(t24,a,i7)')'Explicit coord. initialization'
          else
             write(6,'(t24,a,i7)')'Explicit coord. step =',istep
          end if
          if (cml_p) call cmlStartStep(mainXML, type='ECS', index=istep)

        end select
        write(6,'(2a)') '                        ',
     &                    '=============================='

!       Print atomic coordinates
        call outcoor( ucell, xa, na_u, ' ', writec )
        call siesta_write_positions()
      endif

!     Actualize things if variable cell 
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

!       Will print k-points also
        call setup_Kpoint_grid( ucell )

        call re_alloc( eo, 1, no_u, 1, nspin, 1, maxk, 'eo',
     &                 'state_init')
        call re_alloc( qo, 1, no_u, 1, nspin, 1, maxk, 'qo',
     &                 'state_init' )

!       Find required supercell
        if (gamma) then
           nsc(1:3) = 1
        else if (fixauxcell) then
           nsc(1:3) = nscold(1:3)
        else
           do i = 1 , 3
              veclen = sqrt(ucell(1,i)**2+ucell(2,i)**2+ucell(3,i)**2)
              nsc(i) = 1 + 2 * ceiling( rmaxh / veclen )
           end do
           if ( .not. naiveauxcell )
     &          call exact_sc_ag(negl,ucell,na_u,isa,xa,nsc)
        endif

        mscell = 0.0_dp
        do i = 1, 3
           mscell(i,i) = nsc(i)
           if (nsc(i) /= nscold(i)) auxchanged = .true.
           nscold(i) = nsc(i)
        enddo

!       Madelung correction for charged systems 
        if (charnet .ne. 0.0_dp) then
          call madelung(ucell, shape, charnet, Emad)
        endif

      endif
!     End variable cell actualization

!     Auxiliary supercell
!     Do not move from here, as the coordinates might have changed
!     even if not the unit cell
      call superc(ucell, scell, nsc)

!     Print unit cell and compute cell volume
!     Possible BUG: 
!     Note that this volume is later used in write_subs and the md output
!     routines, even if the cell later changes.
      if (IOnode) call outcell(ucell)
      volume_of_some_cell = volcel(ucell)

!     Use largest possible range in program, except hsparse...
!     2 * rmaxv: Vna overlap
!     2 * rmaxo: orbital overlap
!     rmaxo + rmaxkb: Non-local KB action
!     2.0_dp * (rmaxo+rmaxkb) : Orbital interaction through KB projectors
      rmax = max( 2._dp*rmaxv, 2._dp*rmaxo, rmaxo+rmaxkb)

      if (.not. negl) then
        rmax = max(rmax, 2.0_dp * (rmaxo+rmaxkb) )
      endif

!     Check if any two atoms are unreasonably close
      call proximity_check(rmax)

      ! Ensure sparsity pattern is empty
      call delete(sparse_pattern)
      ! sadly deleting the sparse pattern does not necessarily
      ! mean that the arrays are de-associated.
      ! Remember that the reference counter could (in MD)
      ! be higher than 1, hence we need to create "fake"
      ! containers and let the new<class> delete the old
      ! sparsity pattern
      nullify(numh,listhptr,listh)
      allocate(numh(no_l),listhptr(no_l))
      ! We do not need to allocate listh
      ! that will be allocated in hsparse

!     List of nonzero Hamiltonian matrix elements
!     and, if applicable,  vectors between orbital centers

!     Listh and xijo are allocated inside hsparse
!     Note: We always generate xijo now, for COOP and other
!           analyses.
      call delete(xij_2D) ! as xijo will be reallocated
      nullify(xijo)
      call hsparse( negl, scell, nsc, na_s, isa, xa, lasto,
     &              lastkb, iphorb, iphKB, maxnh, gamma,
     $                 set_xijo=.true., folding=folding1)
!
      call globalize_or(folding1,folding)
      if (folding) then
         if (IOnode) then
            print *, "Folding of H and S is implicitly performed"
         endif
      endif
      !
      ! If using domain decomposition, redistribute orbitals
      ! for this geometry, based on the hsparse info. 
      ! The first time round, the initial distribution is a
      ! simple block one (given by preSetOrbitLimits).
      !
      ! Any DM, etc, read from file will be redistributed according
      ! to the new pattern. 
      ! Inherited DMs from a previous geometry cannot be used if the
      ! orbital distribution changes. For now, we avoid changing the
      ! distribution (the variable use_dd_perm is .true. if domain
      ! decomposition is in effect). Names should be changed...

      if (use_dd .and. (.not. use_dd_perm)) then
         call domainDecom( no_u, no_l, maxnh )  ! maxnh intent(in) here
         maxnh = sum(numh(1:no_l))
         ! We still need to re-create Julian Gale's
         ! indexing for O(N) in parallel.
         print "(a5,i3,a20,3i8)",
     $         "Node: ", Node, "no_u, no_l, maxnh: ", no_u, no_l, maxnh
         call setup_ordern_indexes(no_l, no_u, Nodes)
      endif

      ! I would like to skip this alloc/move/dealloc/attach
      ! by allowing sparsity to have pointer targets.
      ! However, this poses a problem with intel compilers,
      ! as it apparently errors out when de-allocating a target pointer
      write(oname,"(a,i0)") "sparsity for geom step ", istep
      call newSparsity(sparse_pattern,no_l,no_u,maxnh,
     &     numh,listhptr,listh, name = oname)
      deallocate(numh,listhptr,listh)
      call attach(sparse_pattern, 
     &     n_col = numh, list_ptr = listhptr, list_col = listh )
      
      ! Copy over xijo array (we can first do it here... :( )
      call newdData2D(tmp_2D,xijo,'xijo')
      deallocate(xijo)
      write(oname,"(a,i0)") "xijo in state_init at geom step ", istep
      call newdSpData2D(sparse_pattern,tmp_2D,block_dist,xij_2D, 
     &     name=oname)
      call delete(tmp_2D) ! decrement container...
      xijo => val(xij_2D)

      write(oname,"(a,i0)") "H_vkb in state_init at geom step ", istep
      call newdSpData1D(sparse_pattern,block_dist,H_vkb_1D,name=oname)
      !if (ionode) call print_type(H_vkb_1D)

      write(oname,"(a,i0)") "EDM in state_init at geom step ", istep
      call newdSpData2D(sparse_pattern,nspin,block_dist,EDM_2D,
     &     name=oname)
      !if (ionode) call print_type(EDM_2D)
      Escf => val(EDM_2D)

      call re_alloc(Dold,1,maxnh,1,nspin,name='Dold',
     .     routine='state_init',copy=.false.,shrink=.true.)
      call re_alloc(Hold,1,maxnh,1,nspin,name='Hold',
     .     routine='state_init',copy=.false.,shrink=.true.)

!     Allocate/reallocate storage associated with Hamiltonian/Overlap matrix
      write(oname,"(a,i0)") "H in state_init at geom step ", istep
      call newdSpData2D(sparse_pattern,nspin,block_dist,H_2D,name=oname)
      !if (ionode) call print_type(H_2D)
      H => val(H_2D)

      write(oname,"(a,i0)") "H_kin in state_init at geom step ", istep
      call newdSpData1D(sparse_pattern,block_dist,H_kin_1D,name=oname)
      !if (ionode) call print_type(H_kin_1D)

      write(oname,"(a,i0)") "S in state_init at geom step ", istep
      call newdSpData1D(sparse_pattern,block_dist,S_1D,name=oname)
      if (ionode) call print_type(S_1D)
      S => val(S_1D)

!     Find overlap matrix 
      call overlap( na_u, na_s, no_s, scell, xa, indxua, rmaxo, maxnh,
     &              lasto, iphorb, isa, numh, listhptr, listh, S )

      if ( Gamma ) then
         ! Here we create the super-cell offsets
         call re_alloc(isc_off,1,3,1,1)
         isc_off(:,:) = 0
      else
         call xij_offset(ucell,nsc, na_u,xa,lasto, 
     &        xij_2D, isc_off,
     &        Bcast=.true.)
      end if

      ! Initialize density matrix
      call new_dm( auxchanged,DM_history,DM_2D,EDM_2D )
      Dscf => val(DM_2D)
      Escf => val(EDM_2D)

      ! Initialize energy-density matrix to zero for first call to overfsm
#ifdef TRANSIESTA
      ! Only part of Escf is updated in TS, so if it is put as zero here
      ! a continuation run gives bad forces.
      if(.not.TSrun) then
!$OMP parallel workshare default(shared)
         Escf(:,:) = 0.0_dp
!$OMP end parallel workshare
      end if
#else
!$OMP parallel workshare default(shared)
      Escf(:,:) = 0.0_dp
!$OMP end parallel workshare
#endif /* TRANSIESTA */

!     Check for size of Pulay auxiliary matrices
      call init_pulay_arrays()

#ifdef TEST_IO
      ! We test the io-performance here
      call time_io(nspin,H_2D)
#endif

#ifdef TRANSIESTA
!     If onlyS, Save overlap matrix and exit                       ! MPN
      if (onlyS) then           ! MPN 
         i = 0
         ix = 0
         fname = fname_TSHS(slabel,i,ix,.true.)
         ! We include H as S, well-knowing that we only write one of
         ! them, there is no need to allocate space for no reason!
         call ts_write_tshs(fname, 
     &        .true., Gamma, ts_Gamma_file,
     &        ucell, nsc, isc_off, na_u, no_s, nspin,
     &        ts_kscell_file, ts_kdispl_file,
     &        xa, lasto, 
     &        H_2D, S_1D, indxuo, 
     &        dummyEf, dummyQtot, Temp,0,0)
        call bye( 'Save overlap matrix and exit' ) ! Exit siesta   ! MPN 
      endif                                                        ! MPN

      ! When the user requests to only do an analyzation, we can call
      ! appropriate routines and quit
      if ( TS_Analyze ) then

         ! Force the creation of the full sparsity pattern
         call ts_sparse_init(slabel,IsVolt, N_Elec, Elecs, 
     &        ucell, nsc, na_u, xa, lasto, block_dist, sparse_pattern, 
     &        Gamma, isc_off)

         ! create the tri-diagonal matrix
         call ts_tri_analyze( block_dist, sparse_pattern , N_Elec,
     &        Elecs, ucell, na_u, lasto, nsc, isc_off,
     &        opt_TriMat_method )

         ! Bye also waits for all processors
         call bye('transiesta analyzation performed')
      end if

      ! In case the user is requesting a Fermi-correction
      ! we need to delete the TS_FERMI file after each iteration
      if ( TSmode .and. TS_RHOCORR_METHOD == TS_RHOCORR_FERMI
     &     .and. IONode ) then
          ! Delete the TS_FERMI file (enables
          ! reading it in and improve on the convergence)
          if ( file_exist('TS_FERMI') ) then
             i = 23455 ! this should just not be used any were...
             ! Delete the file...
             open(unit=i,file='TS_FERMI')
             close(i,status='delete')
          end if
       end if
      
#endif /* TRANSIESTA */

#ifdef CDF
      if (writedm_cdf) then
         call setup_dm_netcdf_file( maxnh, no_l, nspin,
     &                              no_s, indxuo,
     &                              numh, listhptr, listh)
      endif
      if (writedm_cdf_history) then
         call setup_dm_netcdf_file( maxnh, no_l, nspin,
     &                              no_s, indxuo,
     &                              numh, listhptr, listh,
     &                              istep)
      endif
      if (writedmhs_cdf) then
         call setup_dmhs_netcdf_file( maxnh, no_l, nspin,
     &                              no_s, indxuo,
     &                              numh, listhptr, listh,
     &                              s)
      endif
      if(writedmhs_cdf_history) then
         call setup_dmhs_netcdf_file( maxnh, no_l, nspin,
     &                              no_s, indxuo,
     &                              numh, listhptr, listh,
     &                              s,
     &                              istep)
      endif
#endif
      call timer( 'STinit', 2 )
#ifdef DEBUG
      call write_debug( '  POS state_init' )
#endif
!--------------------------------------------------------------------------- END
      END subroutine state_init
      END module m_state_init
