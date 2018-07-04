! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      MODULE m_state_init

      private
      public :: state_init

      CONTAINS

      subroutine state_init( istep )
      use Kpoint_grid,       only: setup_Kpoint_grid, maxk
      use m_new_dm,          only: new_dm
      use m_proximity_check, only: proximity_check
      use siesta_options
      use units, only: Ang
      use sparse_matrices, only: maxnh, numh, listh, listhptr
      use sparse_matrices, only: Dscf, Dold, Escf, Eold, Hold
      use sparse_matrices, only: xijo, H, S, H_kin, H_vkb
      use siesta_geom
      use atomlist,          only: iphorb, iphkb, rco, rckb, indxua,
     &                             no_s, rmaxkb, rmaxo, no_u, lastkb,
     &                             lasto, superc, rmaxv, indxuo
      use atomlist,          only: no_l, iza
      use alloc,             only: re_alloc, de_alloc, alloc_report
      use m_hsparse,         only: hsparse
      use m_overlap,         only: overlap
      use m_check_supercell, only: check_sc_factors
      use siesta_cml,   only: cml_p, cmlStartStep, mainXML
      use siesta_cml,   only: cmlStartPropertyList, cmlEndPropertyList
      use siesta_cml,   only: cmlAddProperty
      use zmatrix,           only: lUseZmatrix, write_zmatrix
      use m_energies,        only: Emad
      use write_subs
      use m_ioxv,            only: ioxv
      use m_steps
      use parallel,          only: IOnode, node, nodes
      use m_spin,            only: nspin
      use m_rmaxh
      use m_pulay,           only: init_pulay_arrays
      use m_eo
      use m_gamma
      use files,             only: slabel
      use m_mpi_utils,       only: globalize_or
      use m_mpi_utils,       only: globalize_sum
      use domain_decom,      only: domainDecom, use_dd, use_dd_perm
      use fdf,               only: fdf_get
      use sys,               only: message, die
#ifdef TRANSIESTA
      use m_ts_options, only : onlyS
      use sys, only : bye
      use m_ts_io, only : ts_iohs
      use files, only: slabel, label_length
#endif
#ifdef CDF
      use iodm_netcdf, only: setup_dm_netcdf_file
      use iodmhs_netcdf, only: setup_dmhs_netcdf_file
#endif

      implicit none

      integer            :: istep, nnz
      real(dp)           :: veclen      ! Length of a unit-cell vector
      real(dp)           :: rmax
      logical            :: cell_can_change
      integer            :: i, ix, iadispl, ixdispl
      logical            :: auxchanged   ! Auxiliary supercell changed?
      logical            :: folding, folding1
      logical            :: diag_folding, diag_folding1
      logical            :: foundxv  ! dummy for call to ioxv
      real(dp)           :: dummy_qspin(8)
      external           ::  madelung, timer
      real(dp), external :: volcel

#ifdef TRANSIESTA
      integer                       :: ts_kscell_file(3,3) = 0
      real(dp)                      :: ts_kdispl_file(3) = 0.0
      logical                       :: ts_gamma_scf_file = .true.
      character(len=label_length+6) :: fname
      integer                       :: fnlength
      real(dp)                      :: dummyef=0.0, dummyqtot=0.0
      real(dp),             pointer :: dummyH(:,:)
#endif
!------------------------------------------------------------------------- BEGIN
      call timer( 'IterGeom', 1 )

      call timer( 'state_init', 1 )

      istp = istp + 1

      if (IOnode) then
        write(6,'(/,t22,a)') repeat('=',36)
        select case (idyn)
        case (0)
           if (nmove == 0) then
              write(6,'(t25,a)') 'Single-point calculation'
              if (cml_p) call cmlStartStep(mainXML, type='Single-Point',
     $                        index=istp)
           else
              if (broyden_optim) then
                 write(6,'(t25,a,i6)') 'Begin Broyden opt. move = ',
     $                istep
              else if (fire_optim) then
                 write(6,'(t25,a,i6)') 'Begin FIRE opt. move = ',
     $                istep
              else
                 write(6,'(t25,a,i6)') 'Begin CG opt. move = ',
     $                istep
              end if
              if (cml_p) call cmlStartStep(mainXML, type='Geom. Optim',
     $             index=istp)
           endif

!         Print Z-matrix coordinates
          if (lUseZmatrix) then
             call write_Zmatrix()
          endif
        case (1, 3)
           if (iquench > 0 ) then
              write(6,'(t25,a,i6)') 'Begin MD quenched step = ',
     $             istep
              if (cml_p) call cmlStartStep(mainXML, type='MD-quenched',
     $             index=istep)
           else
              write(6,'(t25,a,i6)') 'Begin MD step = ',
     $             istep
              if (cml_p) call cmlStartStep(mainXML, type='MD',
     $             index=istep)
           endif
        case (2,4,5)
           write(6,'(t25,a,i6)') 'Begin MD step = ', istep
           if (cml_p) call cmlStartStep(mainXML, type='MD', index=istep)
        case (6)
          write(6,'(t25,a,i6)') 'Begin FC step = ',istep
          if (cml_p) call cmlStartStep(mainXML, type='FC', index=istep)

          if (istep .eq. 0) then
            write(6,'(t25,a)') 'Undisplaced coordinates'
          else
            iadispl = (istep-mod(istep-1,6))/6+ia1
            ix = mod(istep-1,6)+1
            ixdispl = (ix - mod(ix-1,2) +1)/2
            write(6,'(t26,a,i0,/,t26,a,i1,a,f10.6,a)') 'displace atom ',
     &           iadispl,'in direction ',ixdispl,' by', dx/Ang,' Ang'
          endif

        case (8)
          write(6,'(t25,a,i6)') 'Begin Server step = ',istep
          if (cml_p) call cmlStartStep(mainXML, type='FS', index=istep)

       end select
       write(6,'(t22,a)') repeat('=',36)

!      Print atomic coordinates
       call outcoor( ucell, xa, na_u, ' ', writec )

!       Save structural information in crystallographic format
!       (in file SystemLabel.STRUCT_OUT),
!       canonical Zmatrix (if applicable), and CML record

       call siesta_write_positions(moved=.false.)
       
      endif

      ! Write the XV file for single-point calculations, so that
      ! it is there at the end for those users who rely on it
      call ioxv( 'write', ucell, vcell, na_u, isa, iza, xa, va, 
     &           foundxv)

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

      if ( cell_can_change .and. istep /= inicoor ) then

!       Madelung correction for charged systems 
        if (charnet .ne. 0.0_dp) then
          call madelung(ucell, shape, charnet, Emad)
        end if

        if ( .not. gamma ) then

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

        end if
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

!     List of nonzero Hamiltonian matrix elements
!     and, if applicable,  vectors between orbital centers

!     Listh and xijo are allocated inside hsparse
!     Note: We always generate xijo now, for COOP and other
!           analyses.
      call hsparse( negl, scell, nsc, na_s, isa, xa, lasto,
     &              lastkb, iphorb, iphKB, maxnh, gamma,
     $              set_xijo=.true., folding=folding1,
     $              diagonal_folding=diag_folding1,
     $              debug_folding=fdf_get('debug-folding',.false.))
!
      call globalize_or(diag_folding1,diag_folding)
      call globalize_or(folding1,folding)
      if (diag_folding .and. gamma) then
         call message("WARNING","Gamma-point calculation " //
     $         "with interaction between periodic images")
         call message("WARNING",
     $           "Some features might not work optimally:")
         call message("WARNING",
     $           "e.g. DM initialization from atomic data")
         if (harrisfun) call die("Harris functional run needs " //
     $                              "'force-aux-cell T'")

      else if (folding) then
         if (gamma) then
               call message("INFO","Gamma-point calculation " //
     $                      "with multiply-connected orbital pairs")
               call message("INFO",
     $              "Folding of H and S implicitly performed")
               call check_cohp()
         else
            write(6,"(a,/,a)") "Non Gamma-point calculation " //
     $           "with multiply-connected orbital pairs " //
     $           "in auxiliary supercell.",
     $           "Possible internal error. " //
     $           "Use 'debug-folding T' to debug."
            call die("Inadequate auxiliary supercell")
         endif
      endif
!
      call globalize_sum(maxnh,nnz)
      if (cml_p) then
         call cmlStartPropertyList(mainXML,title='Orbital info')
         call cmlAddProperty(xf=mainXML, value=no_u,
     $        title='Number of orbitals in unit cell',
     $        dictref='siesta:no_u', units="cmlUnits:countable")
         call cmlAddProperty(xf=mainXML, value=nnz,
     $        title='Number of non-zeros',
     $        dictref='siesta:nnz', units="cmlUnits:countable")
         call cmlEndPropertyList(mainXML)
      endif
      !
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

      ! Note that the resizing of Dscf is done inside new_dm below
      call re_alloc(Dold,1,maxnh,1,nspin,name='Dold',
     .     routine='sparseMat',copy=.false.,shrink=.true.)
      call re_alloc( Eold, 1, maxnh, 1, nspin, 'Eold',
     .     routine='sparseMat',copy=.false.,shrink=.true.)
      call re_alloc( Escf, 1, maxnh, 1, nspin, 'Escf',
     .     routine='sparseMat',copy=.false.,shrink=.true.)

      if(mixH) then
         call re_alloc(Hold,1,maxnh,1,nspin,name='Hold',
     .     routine='sparseMat',copy=.false.,shrink=.true.)
      end if

      
!     Allocate/reallocate storage associated with Hamiltonian/Overlap matrix
      call re_alloc( H, 1, maxnh, 1, nspin, 'H',
     .              routine='sparseMat',shrink=.true.,copy=.false.)
      call re_alloc(H_kin,1,maxnh,name='H_kin',routine='sparseMat',
     &              shrink=.true., copy=.false. )
      call re_alloc(H_vkb,1,maxnh,name='H_vkb',routine='sparseMat',
     &              shrink=.true., copy=.false. )
      call re_alloc(S,1,maxnh,name='S',routine='sparseMat',
     .              copy=.false.,shrink=.true.)

!     Find overlap matrix 
      call overlap( na_u, na_s, no_s, scell, xa, indxua, rmaxo, maxnh,
     &              lasto, iphorb, isa, numh, listhptr, listh, S )

!
!     Here we could also read a Hamiltonian, either to proceed to
!     the analysis section (with nscf=0) or to start a mix-H scf cycle.
!
!     Initialize density matrix
      ! The resizing of Dscf is done inside new_dm
      call new_dm( auxchanged, numh, listhptr, listh, Dscf)
      if (nspin > 1) call print_spin(dummy_qspin)

!     Check for size of Pulay auxiliary matrices
      call init_pulay_arrays()

#ifdef TRANSIESTA
!     If onlyS, Save overlap matrix and exit                       ! MPN
      if (onlyS) then                                              ! MPN    
        fname = trim(slabel) // '.onlyS'
        fnlength = len_trim(fname)
        nullify(dummyH)
        call re_alloc( dummyH, 1, maxnh, 1, 1, 'dummyH', 'state_init' )

        call ts_iohs('write   ', gamma, .true., no_u, no_s, 
     &       1, indxuo, maxnh, numh, listhptr, listh, dummyH, 
     &       S, dummyqtot, temp, xijo, fnlength, fname, na_u, 
     &       lasto, isa, dummyef, 
     &       ucell, ts_kscell_file, ts_kdispl_file,            
     &       ts_gamma_scf_file, xa, 0, 0 )
        call de_alloc( dummyH, 'dummyH', 'state_init' )
        call bye( 'Save overlap matrix and exit' ) ! Exit siesta   ! MPN 
      endif                                                        ! MPN    
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
      call timer( 'state_init', 2 )

!--------------------------------------------------------------------------- END
      END subroutine state_init

      subroutine check_cohp()
      use siesta_options, only: write_coop
      use sys,            only: message
      
      if (write_coop) then
         call message("WARNING","There are multiply-connected "//
     $                          "orbitals.")
         call message("WARNING","Your COOP/COHP analysis might " //
     $                          "be affected by folding.")
         call message("WARNING",'Use "force-aux-cell T "' //
     $                          'or k-point sampling')
      endif
      end subroutine check_cohp

      END module m_state_init
