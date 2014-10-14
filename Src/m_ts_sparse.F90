!
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
!
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2013, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

! A module that supplements the reduced memory TranSIESTA version.
! It greatly reduces the memory requirement of transiesta, as well
! as making the code more clearer in intent.
! This module is probably the most *hard* to understand part.

module m_ts_sparse

  use class_Sparsity
  use class_iSpData1D

  use precision, only : dp

  implicit none

  ! The full Transiesta sparsity region in the unit-cell equivalent
  ! region, this is used to accomodate the Hamiltonian and overlap
  ! matrices. In principle this could be made obsolete. However,
  ! that task seems more cumbersome than worthy of notice (for
  ! the moment).
  type(Sparsity), save :: ts_sp_uc ! TS-GLOBAL (UC)

  ! We will save the "update region sparsity"
  ! Note that this is NOT the same as the sparsity pattern
  ! provided by SIESTA!
  ! After using this sparsity pattern, we do not need the listud[g] arrays
  ! as this sparsity pattern is a reflection of the mask used by listud[g]
  ! Lastly the usage of a MASKED sparsity pattern, will reduce the clutter
  ! between MPI and non-MPI codes as it will be the same array in both 
  ! circumstances.
  type(Sparsity), save :: tsup_sp_uc ! TS-update-GLOBAL (UC)

  ! We will save the local "update region sparsity"
  ! Note that this is NOT the same as the sparsity pattern
  ! provided by SIESTA!
  ! After using this sparsity pattern, we do not need the listud[g] arrays
  ! as this sparsity pattern is a reflection of the mask used by listud[g]
  ! This reflects the local sparsity pattern updated elements.
  type(Sparsity), save :: ltsup_sp_sc ! TS-update-local (SC)

  ! This is an index array which points from ltsup_sp_sc to the local siesta
  ! sparsity pattern.
  ! TODO: check how much speed we gain from not searching in the sparsity pattern
  type(iSpData1D), save :: ltsup_sc_pnt

  ! The offsets for the supercells
  real(dp), pointer, save :: sc_off(:,:) => null()

#ifdef MPI
  ! The reduction of the calculated sparse patterns of the Green's function
  ! at non-zero bias.
  ! It can be read as this:
  !  d_dist(Node,1) is the starting orbital that is to be transferred
  !  d_dist(Node,2) is the ending orbital that is to be transferred
  !  d_dist(Node,3) is the number of elements to be transferred
  !  d_dist(Node,4) is the stride in the elements
!  integer, save, allocatable :: d_dist(:,:)
#endif

  integer, parameter :: CUTHILL_MCKEE = 0
  integer, parameter :: CUTHILL_MCKEE_Z_PRIORITY = 1
  integer, parameter :: PAPIOR = 2

contains

! This routine setups what-ever is needed to do the
! memory reduced TranSIESTA code.
! This means collecting information about which region needs
! update, etc.
  subroutine ts_sparse_init(slabel, &
       IsVolt, N_Elec, Elecs, &
       ucell, nsc, na_u,xa,lasto, block_dist,sparse_pattern, Gamma, &
       isc_off)

    use class_OrbitalDistribution

    use alloc
#ifdef MPI
    use mpi_siesta, only : MPI_Comm_Self
#endif 
    use parallel, only: IONode

    use m_ts_electype
    use m_ts_method
#ifdef TRANSIESTA_DEBUG
    use m_ts_debug
    use parallel, only: Node
#endif
    use m_region
    use m_sparsity_handling

!    use m_monitor
! **********************
! * INPUT variables    *
! **********************
    character(len=*), intent(in) :: slabel
    logical, intent(in) :: IsVolt ! bias calculation
    integer, intent(in) :: N_Elec
    type(Elec), intent(inout) :: Elecs(N_Elec)
    ! Unit cell
    real(dp), intent(in) :: ucell(3,3)
    ! Number of super-cells in each direction
    integer, intent(in) :: nsc(3)
    ! Number of atoms in the unit-cell
    integer, intent(in) :: na_u
    ! Atomic coordinates
    real(dp), intent(in) :: xa(3,na_u)
    ! Last orbital of the equivalent unit-cell atom
    integer, intent(in) :: lasto(0:na_u)
    ! The distribution for the sparsity-pattern
    type(OrbitalDistribution), intent(inout) :: block_dist
    ! SIESTA local sparse pattern (not changed)
    type(Sparsity), intent(inout) :: sparse_pattern
    ! Whether we have xij or not
    logical, intent(in) :: Gamma
    ! supercell indices
    integer, intent(in) :: isc_off(3,product(nsc))

! **********************
! * LOCAL variables    *
! **********************
    type(OrbitalDistribution) :: dit
    ! We create a temporary sparsity pattern which removes
    ! all cross-connections across the electrode transport direction.
    type(Sparsity) :: tmp_sp
    type(tRegion) :: r_oE(N_Elec), r_tmp1, r_tmp2
    ! Temporary arrays for knowing the electrode size
    logical :: bool
    integer :: no_u_TS, iEl, i

    ! Number of orbitals in TranSIESTA
    no_u_TS = nrows_g(sparse_pattern) - no_Buf

    ! Do a crude check of the sizes
    if ( no_u_TS <= sum(TotUsedOrbs(Elecs)) ) then
       call die("The contact region size is &
            &smaller than the electrode size. Please correct.")
    end if

    ! Create the ts-offsets
    if ( Gamma ) then
       ! Initialize the sc_off array
       call re_alloc(sc_off,1,3,1,1)
       sc_off(:,:) = 0._dp
    else
       call re_alloc(sc_off,1,3,1,product(nsc))
       sc_off(:,:) = matmul(ucell,isc_off)
    end if

    ! Some of the routines does the exact same thing, however,
    ! here we do it to ensure that any terms crossing the unit-cell
    ! boundary are removed.
    if ( no_Buf > 0 ) then
       ! Remove buffer atoms...
       call Sp_remove_region(block_dist,sparse_pattern,r_oBuf,tmp_sp)
    else
       tmp_sp = sparse_pattern
    end if

#ifdef TRANSIESTA_DEBUG
    if(IONode)write(*,*)'As is sparsity pattern (20)'
    call sp_to_file(20+Node,sparse_pattern)
#endif

    ! Remove all electrode to other side connections
    ! this only has effect when we cross 
    do iEl = 1 , N_Elec

       ! Create electrode region
       call region_range(r_oE(iEl), Elecs(iEl)%idx_o, &
            Elecs(iEl)%idx_o+TotUsedOrbs(Elecs(iEl))-1)

       ! Remove the connections that cross the boundary
       ! starting from this electrode
       call region_connect(r_oE(iEl), block_dist, tmp_sp, r_tmp1)
       call region_union(r_oE(iEl), r_tmp1, r_tmp2)

       ! Remove connections from this electrode across the boundary...
       call Sp_remove_crossterms(block_dist,tmp_sp,product(nsc),isc_off, &
            Elecs(iEl)%t_dir, &
            tmp_sp, r = r_tmp2)

#ifdef TRANSIESTA_DEBUG
       if(IONode)write(*,*)'Created TS-NO Cross (50)'
       call sp_to_file(40+Node+10*iEl,tmp_sp)
#endif

    end do

    if ( IONode .and. Gamma ) then
       write(*,'(a,a)')'transiesta: ',&
            'We cannot assure cross-boundary connections in Gamma calculations.'
       write(*,'(a,a)')'transiesta: ',&
            'Ensure that the electrode says: Principal cell is perfect'
    end if

    do iEl = 1 , N_Elec - 1
       call region_delete(r_tmp1)
       call region_delete(r_tmp2)
       do i = iEl + 1 , N_Elec
          call region_copy(r_tmp2,r_tmp1)
          call region_union(r_oE(i),r_tmp1,r_tmp2)
       end do
       call Sp_remove_region2region(block_dist,tmp_sp,r_oE(iEl),r_tmp2,tmp_sp)

    end do
    call region_delete(r_tmp1)
    call region_delete(r_tmp2)
    do iEl = 1 , N_Elec
       call region_delete(r_oE(iEl))
    end do

#ifdef TRANSIESTA_DEBUG
    if(IONode)write(*,*)'Created TS-NO ELEC Cross (90)'
    call sp_to_file(90+Node,tmp_sp)
#endif

    if ( IsVolt ) then 

       ! Create the update region (a direct subset of the local sparsity pattern)
       ! Hence it is still a local sparsity pattern.
       call ts_Sparsity_Update(block_dist,tmp_sp, N_Elec, Elecs, &
            ltsup_sp_sc)
       
       ! assign distribution array
       !call ts_init_distribution(block_dist,sparse_pattern)
       
       if ( IONode ) then
          write(*,'(/,a)') 'Created the TranSIESTA local update sparsity pattern:'
          call print_type(ltsup_sp_sc)
       end if

#ifdef TRANSIESTA_DEBUG
       if(IONode)write(*,*)'Created TS-local UP (300)'
       call sp_to_file(300+Node,ltsup_sp_sc)
#endif

       ! Create the pointer from the local transiesta update sparsity 
       ! to the local siesta sparsity
       call ts_Sparsity_Subset_pointer(block_dist,sparse_pattern,ltsup_sp_sc, &
            ltsup_sc_pnt)

    end if

    ! Create the global transiesta H(k), S(k) sparsity pattern
    call ts_Sparsity_Global(block_dist,tmp_sp, &
         N_Elec, Elecs, &
         ts_sp_uc)

    ! The update sparsity pattern can be simplied to the H,S sparsity
    ! pattern, if all electrodes have certain options to be the same.
    bool = all(Elecs(:)%Bulk) .and. all(Elecs(:)%DM_update == 1)
    bool = bool .or. all( .not. Elecs(:)%Bulk )
    bool = bool .or. all( Elecs(:)%DM_update == 2 )

    if ( IONode .and. .not. bool ) then
       write(*,'(/,a)') 'Created the TranSIESTA H,S sparsity pattern:'
       call print_type(ts_sp_uc)
    end if

#ifdef TRANSIESTA_DEBUG
    if(IONode)write(*,*)'Created TS-Global HS (100)'
    call sp_to_file(100+Node,ts_sp_uc)
#endif

    ! In order to ensure that the electrodes are in the
    ! tri-diagonal sparsity pattern, we can easily create
    ! the full sparsity pattern with the electrodes included
    ! and then recreate the tri-diagonal sparsity pattern
    ! This is probably the crudest way of doing it.
#ifdef MPI
    call newDistribution(nrows_g(ts_sp_uc),MPI_Comm_Self,dit, &
         name='TranSIESTA UC distribution')
#else    
    call newDistribution(nrows_g(ts_sp_uc),-1,dit, &
         name='TranSIESTA UC distribution')
#endif

    if ( bool ) then

       tsup_sp_uc = ts_sp_uc
       if ( IONode ) then
          write(*,'(/,a)') 'Created the TranSIESTA H,S sparsity pattern.'
          write(*,'(a)') 'TranSIESTA global update sparsity pattern same as H,S'
          call print_type(ts_sp_uc)
       end if

    else

       ! Create the update region (a direct subset of ts_sp_uc)
       call ts_Sparsity_Update(dit,ts_sp_uc, N_Elec, Elecs, &
            tsup_sp_uc)
       
       if ( IONode ) then
          write(*,'(/,a)') 'Created the TranSIESTA global update sparsity pattern:'
          call print_type(tsup_sp_uc)
       end if

    end if

    ! Read in the monitor lists...
    ! initialize the monitor list
!    if ( N_mon == 0 ) then
!       nullify(monitor_list)
!       call read_monitor('TS.DM.Monitor', &
!            dit, tsup_sp_uc, N_mon, monitor_list)
!       if ( N_mon > 0 .and. IONode ) then
!          if ( .not. IsVolt ) then
!             call re_alloc(iu_MON,1,1,1,N_mon)
!          else
!             call re_alloc(iu_MON,1,4,1,N_mon)
!          end if
!          do i = 1 , N_mon
!             ! open the files
!             if ( .not. IsVolt ) then
!                call io_assign(iu_MON(1,i))
!                open(iu_MON(1,i),file=fname_monitor( &
!                     monitor_list(i,1),monitor_list(i,2), &
!                     basename=trim(slabel)//'.TSMON'), &
!                     form='formatted',status='unknown')
!             else
!                call io_assign(iu_MON(1,i))
!                open(iu_MON(1,i),file=fname_monitor( &
!                     monitor_list(i,1),monitor_list(i,2), &
!                     basename=trim(slabel)//'.TSMONL'), &
!                     form='formatted',status='unknown')
!                call io_assign(iu_MON(2,i))
!                open(iu_MON(2,i),file=fname_monitor( &
!                     monitor_list(i,1),monitor_list(i,2), &
!                     basename=trim(slabel)//'.TSMONR'), &
!                     form='formatted',status='unknown')
!                call io_assign(iu_MON(3,i))
!                open(iu_MON(3,i),file=fname_monitor( &
!                     monitor_list(i,1),monitor_list(i,2), &
!                     basename=trim(slabel)//'.TSMONLN'), &
!                     form='formatted',status='unknown')
!                call io_assign(iu_MON(4,i))
!                open(iu_MON(4,i),file=fname_monitor( &
!                     monitor_list(i,1),monitor_list(i,2), &
!                     basename=trim(slabel)//'.TSMONRN'), &
!                     form='formatted',status='unknown')
!             end if
!          end do
!       end if
!    end if

#ifdef TRANSIESTA_DEBUG
    if(IONode)write(*,*)'Created TS-Global update (200)'
    call print_type(tsup_sp_uc)
    call sp_to_file(200+Node,tsup_sp_uc)
#endif

    ! Clean-up
    call delete(tmp_sp)

  end subroutine ts_sparse_init

! Returns the global sparsity pattern for the transiesta region
! Note that this will automatically detect whether there are
! cross connections from left-right, AND remove any
! z-connections (less orbitals to move about, and they should
! not exist).
  subroutine ts_Sparsity_Global(block_dist,s_sp, &
       N_Elec, Elecs, ts_sp)

    use parallel, only : IONode
#ifdef MPI
    use mpi_siesta
#endif
    use class_OrbitalDistribution
    use create_Sparsity_SC
    use m_ts_electype
    use m_ts_method
    use m_sparsity_handling

#ifdef TRANSIESTA_DEBUG
    use m_ts_debug
    use parallel, only: Node
#endif

! **********************
! * INPUT variables    *
! **********************
    ! The SIESTA distribution of the sparsity pattern
    type(OrbitalDistribution), intent(in) :: block_dist
    ! Sparsity patterns of SIESTA (local)
    type(Sparsity), intent(inout) :: s_sp
    ! All the electrodes
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    ! the returned update region.    
    type(Sparsity), intent(inout) :: ts_sp

! **********************
! * LOCAL variables    *
! **********************
    ! We need temporary sparsity patterns which will be deleted
    ! We need this to generate the Transiesta sparsity
    type(Sparsity) :: sp_global, sp_uc

    ! to globalize from the local sparsity pattern (SIESTA)
    ! and afterwards used as the looping mechanisms
    ! to create the mask for creating the UC transiesta pattern
    integer, pointer :: l_ncol(:) => null()
    integer, pointer :: l_ptr(:) => null()
    integer, pointer :: l_col(:) => null()

    ! Also used in non-MPI (to reduce dublicate code)
    integer :: no_l, no_u, uc_n_nzs, n_nzsg

    ! search logical to determine the update region...
    logical, allocatable :: l_HS(:)

    ! Loop-counters
    integer :: io, jo, ind
    integer :: iot, jot
    logical :: UseBulk

    ! Initialize
    call delete(ts_sp)

    call attach(s_sp,nrows=no_l,nrows_g=no_u)

    ! Create the (local) SIESTA-UC sparsity...
#ifdef MPI
    call crtSparsity_SC(s_sp,sp_global, UC=.TRUE.)
    uc_n_nzs = nnzs(sp_global)

    ! point to the local (SIESTA-UC) sparsity pattern arrays
    call Sp_to_Spglobal(block_dist,sp_global,sp_uc)

    ! Delete the local UC sparsity pattern
    call delete(sp_global)

#else
    call crtSparsity_SC(s_sp,sp_uc    , UC=.TRUE.)
    uc_n_nzs = nnzs(sp_uc)
#endif

#ifdef TRANSIESTA_DEBUG
    if(IONode)write(*,*)'Created UC SIESTA sparsity (400)'
    call print_type(sp_uc)
    call sp_to_file(400+Node, sp_uc)
#endif

    ! Write that we have created it
    if ( IONode ) call print_type(sp_uc)

    ! Now we have the globalized SIESTA Unit-cell pattern, make it 
    ! only to the Transiesta sparsity pattern

    ! Immediately point the global arrays to their respective parts
    call attach(sp_uc,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nnzs=n_nzsg)

    ! allocate space for the MASK to create the TranSIESTA GLOBAL region
    allocate(l_HS(n_nzsg))
    call memory('A','L',n_nzsg,'transiesta')

    ! Initialize
    l_HS(:) = .false.

    ! We do not need to check the buffer regions...
    ! We know they will do NOTHING! :)
    do io = 1 , no_u

       iot = orb_type(io)
       if ( iot == TYP_BUFFER ) cycle

       ! The index in the pointer array is retrieved
       do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
          
          ! The unit-cell column index (remember we are looping a UC SP)
          jo = l_col(ind)

          ! If we are in the buffer region, cycle (l_HS(ind) =.false. already)
          ! note, that we have already *checked* ic
          jot = orb_type(jo)
          if ( jot == TYP_BUFFER ) cycle

          if ( iot > 0 ) then
             ! In order to allow to have the update sparsity pattern
             ! as a subset, we require that the Hamiltonian also
             ! has the electrode interconnects
             if ( Elecs(iot)%DM_update == 2 ) then
                UseBulk = .false.
             else
                UseBulk = Elecs(iot)%Bulk
             end if
          else if ( jot > 0 ) then
             if ( Elecs(jot)%DM_update == 2 ) then
                UseBulk = .false.
             else
                UseBulk = Elecs(jot)%Bulk
             end if
          else
             ! we are definitely not in an electrode
             UseBulk = .true.
          end if
             
          if ( UseBulk ) then
             ! here we create the update-density matrix on these criterias:
             !  1) no electrode-electrode connections
             !  2) no only-electrode connections
             !  3) add cross-terms

             l_HS(ind) = any((/iot,jot/)==TYP_DEVICE)

          else
             ! If not usebulk we update everything that is not buffer
             ! We also do not add things that are from one electrode
             ! to another
             if ( all(0 < (/iot,jot/)) ) then
                l_HS(ind) = iot == jot
             else
                l_HS(ind) = .true.
             end if
          end if

       end do

    end do

    ! We now have a MASK of the actual needed TranSIESTA sparsity pattern
    ! We create the TranSIESTA sparsity pattern
    call crtSparsity_SC(sp_uc,ts_sp,MASK=l_HS)

    ! clean-up
    deallocate(l_HS)
    call memory('D','L',n_nzsg,'transiesta')
    
    ! Furthermore, we dont need the SIESTA UC sparsity global...
    call delete(sp_uc)

  end subroutine ts_Sparsity_Global


! Returns the global sparsity pattern for the transiesta region
! Note that this will automatically detect whether there are
! cross connections from left-right, AND remove any
! z-connections (less orbitals to move about, and they should
! not exist).
  subroutine ts_Sparsity_Update(dit,s_sp, N_Elec, Elecs, &
       tsup_sp)

    use parallel, only : IONode, Node
    use geom_helper, only : UCORB
    use create_Sparsity_SC
    use class_OrbitalDistribution
    use m_ts_method
    use m_ts_electype
! **********************
! * INPUT variables    *
! **********************
    type(OrbitalDistribution), intent(inout) :: dit
    type(Sparsity), intent(inout) :: s_sp
    ! the electrodes
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    type(Sparsity), intent(inout) :: tsup_sp

! **********************
! * LOCAL variables    *
! **********************
    ! to globalize from the local sparsity pattern (SIESTA)
    ! and afterwards used as the looping mechanisms
    ! to create the mask for creating the UC transiesta pattern
    integer, pointer :: l_ncol(:) => null()
    integer, pointer :: l_ptr(:) => null()
    integer, pointer :: l_col(:) => null()

    ! Also used in non-MPI (to reduce dublicate code)
    integer :: no_l, no_u, n_nzs

    ! search logical to determine the update region...
    logical, allocatable :: lup_DM(:)
    logical :: direct_LR
    logical :: UseBulk, DM_CrossTerms

    ! Loop-counters
    integer :: lio, io, ind
    integer :: ict, jct
    
    ! Logical for determining the region
    logical :: i_in_C, j_in_C

    ! Initialize the tsup
    call delete(tsup_sp)

    call attach(s_sp,nrows=no_l,nrows_g=no_u,nnzs=n_nzs)

    ! allocate space for the MASK to create the TranSIESTA UPDATE region
    allocate(lup_DM(n_nzs))
    call memory('A','L',n_nzs,'transiesta')

    ! Initialize
    lup_DM(:) = .false.
    direct_LR = .false.

    call attach(s_sp,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)

    ! We do not need to check the buffer regions...
    ! We know they will do NOTHING! :)
    do lio = 1 , no_l

       ! Shift out of the buffer region
       io = index_local_to_global(dit,lio,Node)

       ! If we are in the buffer region, cycle (lup_DM(ind) =.false. already)
       ict = orb_type(io)
       if ( ict == TYP_BUFFER ) cycle

       ! Loop the index of the pointer array
       do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

          ! If we are in the buffer region, cycle (lup_DM(ind) =.false. already)
          ! note, that we have already *checked* io
          jct = orb_type(l_col(ind))
          if ( jct == TYP_BUFFER ) cycle

          if ( ict > 0 ) then
             if ( Elecs(ict)%DM_update == 2 ) then ! update everything
                UseBulk = .false.
             else
                UseBulk = Elecs(ict)%Bulk
                DM_CrossTerms = Elecs(ict)%DM_update == 1
             end if
          else if ( jct > 0 ) then
             if ( Elecs(jct)%DM_update == 2 ) then ! update everything
                UseBulk = .false.
             else
                UseBulk = Elecs(jct)%Bulk
                DM_CrossTerms = Elecs(jct)%DM_update == 1
             end if
          else
             ! we are definitely not in an electrode
             ! just set it to be updated
             UseBulk = .false.
             DM_CrossTerms = .true.
          end if

          ! We check whether it is electrode-connections. 
          ! If, so, they are not used in transiesta:
          if      ( ict > 0 .and. jct > 0 ) then
             ! Remove connections between electrodes
             ! but maintain same electrode updates if not usebulk
             if ( ict == jct .and. .not. UseBulk ) then
                lup_DM(ind) = .true.
             else
                lup_DM(ind) = .false.
             end if

             ! This means that we have an INNER-cell connection
             !if ( jc == l_col(ind) ) then
                ! The first super-cell is the home-unit-cell
                ! hence it means a direct interaction in the unit-cell
                !print *,orb_type((/ic,jc/))
                !direct_LR = .true.
             !end if
          else if ( UseBulk ) then
             ! If UseBulk, we can safely update cross-terms
             ! between the central and the electrodes (the self-energies are the
             ! same)
             ! This will also remove any electrode terms
                
             i_in_C = ict == TYP_DEVICE
             j_in_C = jct == TYP_DEVICE
             
             if ( DM_CrossTerms ) then
                ! the user has requested to also update cross-terms
                lup_DM(ind) = i_in_C .or.  j_in_C
             else
                ! the user has requested to ONLY update the central region
                lup_DM(ind) = i_in_C .and. j_in_C
             end if

          else
             ! If not UseBulk the full Hamiltonian and the 
             ! self-energy terms are used.
             ! Hence, everything in the left-central-right
             ! is updated (note that we are sure to be
             ! in the range [noBufL+1;no_u-noBufR]
             ! Otherwise we needed to test that here.
             lup_DM(ind) = .true.
             ! It makes no sense to ususe UpdateDMCR in case of not using UseBulk
             ! as we update something that is not used...
             
          end if

       end do

    end do

    ! Tell the user about inner-cell connections
    if ( IONode .and. direct_LR ) then
       write(*,*) 'WARNING: Cross connections across the &
            &junction leads to tunneling.'
       write(*,*) 'WARNING: Transiesta will disregard these &
            &cross-terms'
       write(*,*) 'WARNING: Consider increasing the electrodes &
            &in the transport-direction.'
       write(0,*) 'WARNING: Cross connections across the &
            &junction leads to tunneling.'
       write(0,*) 'WARNING: Transiesta will disregard these &
            &cross-terms'
       write(0,*) 'WARNING: Consider increasing the electrodes &
            &in the transport-direction.'
    end if

    ! We now have a MASK of the actual needed TranSIESTA sparsity pattern
    ! We create the TranSIESTA sparsity pattern
    call crtSparsity_SC(s_sp,tsup_sp,MASK=lup_DM)

    call memory('D','L',n_nzs,'transiesta')
    deallocate(lup_DM)

  end subroutine ts_Sparsity_Update

  subroutine ts_Sparsity_Subset_pointer(dit,sp,sub_sp,ipnt)
    use class_OrbitalDistribution

! **********************
! * INPUT variables    *
! **********************
    ! The distribution pattern for everything
    type(OrbitalDistribution), intent(inout) :: dit
    ! The sparsity pattern we wish to point to
    type(Sparsity), intent(inout) :: sp
    ! The sparsity pattern we wish to point from
    type(Sparsity), intent(inout) :: sub_sp
    ! The pointer index
    type(iSpData1D), intent(inout) :: ipnt

! **********************
! * LOCAL variables    *
! **********************
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer, pointer :: sub_ncol(:), sub_ptr(:), sub_col(:)
    integer, pointer :: pnt(:)      => null()

    integer :: no_l, io, j, sub_ind, ind

    call attach(sp,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=no_l)
    call attach(sub_sp,n_col=sub_ncol,list_ptr=sub_ptr,list_col=sub_col, &
         nrows=io)
    if ( io /= no_l ) call die('Could not do index matching due to &
         &inconsistent sparsity patterns')

    ! Clear and create
    call delete(ipnt)
    call newiSpData1D(sub_sp,dit,ipnt,name='TS pointer')
    ! Point to the array
    pnt => val(ipnt)
    pnt(:) = 0

    ! Loop in the subset sparsity pattern
    do io = 1 , no_l

       ! Loop number of entries in the row...
       do j = 1 , sub_ncol(io)

          ! The index in the pointer array is retrieved
          sub_ind = sub_ptr(io) + j

          ! Loop in the super-set sparsity pattern
          super_idx: do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)

             ! If we have the same column index it must be
             ! the same entry they represent
             if ( sub_col(sub_ind) == l_col(ind) ) then
                pnt(sub_ind) = ind
                exit super_idx
             end if

          end do super_idx
       end do
    end do

    if ( any(pnt(:) == 0) ) then
       call die('An index could not be located in the super-set &
            &sparsity pattern. Are you surely having the correct &
            &sparsity?')
    end if

    if ( any(pnt(:) > nnzs(sp)) ) then
       call die('An index could not be located in the super-set &
            &sparsity pattern. Are you surely having the correct &
            &sparsity?')
    end if

  end subroutine ts_Sparsity_Subset_pointer

#ifdef NOT_USED
  ! Resets all degrees of freedom in the
  ! matrix for elements belonging to the central
  ! region AND the cross-terms with the electrodes
  subroutine ts_Reset_D_C(D2)

    use parallel, only : Node
    use class_OrbitalDistribution
    use class_dSpData2D
    use m_ts_method
! **********************
! * INPUT variables    *
! **********************
    type(dSpData2D), intent(inout) :: D2

! **********************
! * LOCAL variables    *
! **********************
    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: sp
    ! to globalize from the local sparsity pattern (SIESTA)
    ! and afterwards used as the looping mechanisms
    ! to create the mask for creating the UC transiesta pattern
    integer, pointer :: l_ncol(:) => null()
    integer, pointer :: l_ptr(:) => null()
    integer, pointer :: l_col(:) => null()
    real(dp), pointer :: a2(:,:)

    ! Also used in non-MPI (to reduce dublicate code)
    integer :: no_l, no_u

    ! Loop-counters
    integer :: io, ic, ind
    integer :: ict, jct

    dit => dist(D2)
    sp  => spar(D2)
    a2  => val(D2)

    call attach(sp,nrows=no_l,nrows_g=no_u)
    call attach(sp,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)

    do io = 1 , no_l

       ! Shift out of the buffer region
       ic = index_local_to_global(dit,io,Node)

       ! If we are in the buffer region, cycle (lup_DM(ind) =.false. already)
       ict = orb_type(ic)
       if ( ict == TYP_BUFFER ) cycle

       ! Loop the index of the pointer array
       do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)

          ! If we are in the buffer region, cycle (lup_DM(ind) =.false. already)
          ! note, that we have already *checked* ic
          jct = orb_type(l_col(ind))
          if ( jct == TYP_BUFFER ) cycle

          ! Do not reset diagonal and cross terms
          ! between electrodes
          if ( ict > 0 .and. jct > 0 ) cycle
          ! Do not reset in device 
          if ( ict == TYP_DEVICE .and. jct == TYP_DEVICE ) cycle

          ! Reset !
          a2(ind,:) = 0._dp

       end do

    end do

  end subroutine ts_Reset_D_C
#endif

  subroutine ts_Optimize(sp,N_Elec, Elecs,na_u,lasto,xa, algorithm)
    ! We return the pivoting indexes for the atoms
    use class_Sparsity
    use alloc, only : re_alloc, de_alloc
    use parallel, only : IONode
    use m_ts_tdir, only : ts_tdir
    use m_ts_method
    use m_ts_electype

    use m_bandwidth
    ! This sparsity pattern must be a full one
    type(Sparsity), intent(in out) :: sp
    ! Electrodes
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    ! number of atoms in the cell, also the last orbitals of each
    ! atom.
    integer, intent(in) :: na_u, lasto(0:na_u)
    ! The atomic coordinates
    ! This means that we can sort against the z-level by assigning 
    ! priority of atoms.
    real(dp), intent(in) :: xa(3,na_u)
    integer, intent(in) :: algorithm
    ! The pivoting array
    integer, pointer :: R(:)
    ! The priority array
    integer, pointer :: a_priority(:)

    integer, pointer :: a_mm(:,:)
    integer, pointer :: l_col(:), ncol(:), tmp(:)

    integer :: ptr
    integer :: na_b, nrg
    integer :: io, ia, iab, iab2, iac, co1, co2
    integer :: in_elec1, in_elec2

    ! Retrieve the lists
    call attach(sp,list_col=l_col,n_col=ncol,nrows_g=nrg)

    nullify(a_mm,R,tmp,a_priority)
#ifdef TS_PIVOT_ELEC
    na_b = na_u - na_Buf - sum(TotUsedAtoms(Elecs(:))) + N_Elec
#else
    na_b = na_u - na_Buf - sum(TotUsedAtoms(Elecs(:)))
#endif
    call re_alloc(a_mm,1,na_b,1,na_b)
    call re_alloc(R,1,na_b)
    call re_alloc(tmp,1,maxval(ncol))
    call re_alloc(a_priority,1,na_b)

    ! First we create the priority array
    ! We assign an integer that has an arbitrary range 
    ! dependent on the size of the system
    ! A high-z gets a high priority
    iab = 0
    in_elec1 = TYP_DEVICE
    do ia = 1 , na_u
       if ( a_isBuffer(ia) ) cycle ! the buffers don't matter
       ! check whether we have a new electrode
       call step_ia(iab,ia,in_elec1)
       if ( iab > na_b ) exit
       if ( 1 <= ts_tdir .and. ts_tdir <= 3 ) then
          ! TODO correct ts_tdir for rotated unit-cells!
          a_priority(iab) = nint(xa(ts_tdir,ia) * 10000._dp)
       else
          a_priority(iab) = 0
       end if
    end do

    ! Initialize the connections
    ! Notice that we build the matrix from backwards
    ! as we need to do the reversed Cuthill-Mckee algorithm
    a_mm(:,:) = 0
    a_mm(na_b,na_b) = 1
    iab = 0
    in_elec1 = TYP_DEVICE
    do ia = 1 , na_u - 1
       if ( a_isBuffer(ia) ) cycle
       ! this will only step if the electrode is changed to another
       ! electrode...
       ! Hence, one element for each electrode
       call step_ia(iab,ia,in_elec1)
       if ( iab > na_b ) exit

       ! Of course the atom connects to itself
       a_mm(iab,iab) = 1

       ! The current atoms orbitals
       co1 = lasto(ia-1) + 1
       co2 = lasto(ia) - co1

       iab2 = iab
       in_elec2 = in_elec1
       do iac = ia + 1 , na_u
          if ( a_isBuffer(iac) ) cycle
          call step_ia(iab2,iac,in_elec2)
          if ( iab2 > na_b ) exit

          ! search for connections
          connection: do io = lasto(iac-1) + 1 , lasto(iac)
             if ( ncol(io) == 0 ) cycle
             ! Retrieve the pointer
             ptr = list_ptr(sp,io)
             tmp(1:ncol(io)) = l_col(ptr+1:ptr+ncol(io)) - co1
             where ( tmp < 0 ) tmp = co2 + 1
             if ( any(tmp(1:ncol(io)) <= co2) ) then
                ! Create the connection to the other atoms
                ! here we do the reverse notation
                a_mm(iab,iab2) = 1
                a_mm(iab2,iab) = 1
                exit connection
             end if
          end do connection

       end do

    end do

    ! deallocate unused array
    call de_alloc(tmp)

    ! Calculate the pivoting...
    call BandWidth_pivoting(algorithm, na_b, a_mm, R, &
         priority = a_priority )

    ! We do not need the priority any-more
    call de_alloc(a_priority)
         
    if ( IONode ) then

       ! Create the pivoting array
       call re_alloc(tmp,1,na_b)
#ifdef TS_PIVOT_ELEC
       tmp(:)   = huge(1)
#else
       tmp(:)   = 0
#endif
       iab      = 0
       in_elec1 = TYP_DEVICE
       do ia = 1 , na_u
          ! notice that it is reversed
          if ( a_isBuffer(ia) ) cycle ! the buffers are not pivoted
          ! check whether we have a new electrode
          call step_ia(iab,ia,in_elec1)
          if ( iab > na_b ) exit
          ! This asserts that we point to the
          ! first index of each electrode
#ifdef TS_PIVOT_ELEC
          tmp(iab) = min(tmp(iab),ia) 
#else
          tmp(iab) = max(tmp(iab),ia) 
#endif
       end do

#ifdef TS_PIVOT_ELEC
       ! We need to correct for the pivoting of the electrodes
       ! Thus we need to check for entire electrode block
       ! pivoting, remember that the entire electrode is defined
       ! as 1 entity. This might not make sense, but is necessary
       ! for the transiesta algorithm
       iab = 0
       in_elec1 = TYP_DEVICE
       do ia = 1 , na_u
          ! notice that it is reversed
          if ( a_isBuffer(ia) ) cycle ! the buffers are never pivoted
          ! check whether we have a new electrode
          call step_ia(iab,ia,in_elec1)
          if ( iab > na_b ) exit
          ! if we are not in an electrode, simply skip...
          if ( in_elec1 /= TYP_DEVICE ) cycle
          
          ! we are not pivoting the electrode! YUPPI
          if ( R(iab) == iab ) cycle

          ! we will pivot an electrode.. :(
       end do
#endif       

       write(*,'(a)') 'transiesta: Tri-diagonal blocks can be &
            &optimized by the following pivoting'

       write(*,'(t5,a4,tr2,a2,tr2,a4)') 'From','->','To'

       iab = 0
       in_elec1 = TYP_DEVICE
       do ia = 1 , na_u

          ! we do not consider buffer atoms
          if ( a_isBuffer(ia) ) cycle

          ! step the R(...) counter according to the step
          call step_ia(iab,ia,in_elec1)
          if ( iab > na_b ) exit
#ifndef TS_PIVOT_ELEC
          ! Skip any electrodes
          if ( in_elec1 /= TYP_DEVICE ) cycle
#endif

          ! in case no pivoting is required
          if ( iab == R(iab) ) cycle
          
          if ( in_elec1 /= TYP_DEVICE ) then
             write(*,'(t5,i4,tr2,a2,tr2,i4,a)') &
                  tmp(R(iab)) , '->', tmp(iab) ,' *'
          else
             write(*,'(t5,i4,tr2,a2,tr2,i4)') &
                  tmp(R(iab)) , '->', tmp(iab)
          end if

       end do
       
       call de_alloc(tmp)

    end if
           
    call de_alloc(a_mm)
    call de_alloc(R)

  contains

    pure function offset(N_Elec,Elecs,ia)
      integer, intent(in) :: N_Elec
      type(Elec), intent(in) :: Elecs(N_Elec)
      integer, intent(in) :: ia
      integer :: offset
      offset = sum(TotUsedAtoms(Elecs(:)), MASK=Elecs(:)%idx_a < ia )
    end function offset

#ifdef MUMPS
    ! Add a routine for analysing the structure
    ! This will let users easily get access to the 
    ! different orderings in MUMPS and check them
    ! against each other...

#endif

    subroutine step_ia(iab,ia,in_elec)
      integer, intent(inout) :: iab
      integer, intent(in) :: ia
      integer, intent(inout) :: in_elec
      integer :: at
      at = atom_type(ia)
#ifdef TS_PIVOT_ELEC
      if ( at > TYP_DEVICE ) then
         if ( at == in_elec ) then
            ! do nothing
         else
            in_elec = at
            iab = iab + 1
         end if
      else
         in_elec = TYP_DEVICE
         iab = iab + 1
      end if
#else
      if ( at > TYP_DEVICE ) then
         ! we are in an electrode
         if ( in_elec == TYP_DEVICE ) then
            ! only once, increment it...
            iab = iab + 1
         end if
      else if ( in_elec > TYP_DEVICE ) then
         ! do nothing, just ensure no incrementing
      else
         iab = iab + 1
      end if
      in_elec = at
#endif
    end subroutine step_ia

  end subroutine ts_Optimize

#ifdef MPI_NEW
  subroutine ts_init_distribution(dit,sparse_pattern)

    use parallel, only : Node, Nodes
    use intrinsic_missing, only : UNIQC
    use geom_helper, only : UCORB

    use mpi_siesta, only : MPI_Integer, MPI_Sum
    use mpi_siesta, only : MPI_AllReduce, MPI_Comm_World

    use class_Sparsity
    use class_OrbitalDistribution

    type(OrbitalDistribution), intent(in) :: dit
    type(Sparsity), intent(inout) :: sparse_pattern

    ! pointers to the sparsity pattern
    integer, pointer :: l_ncol(:) => null()
    integer, pointer :: l_ptr(:) => null()
    integer, pointer :: l_col(:) => null()
    integer, allocatable :: tmp(:)

    integer :: io, jo, ind, no_l, no_u, N

    integer :: MPIerror

    ! Allocate all the nodes distribution sizes
    if ( allocated(d_dist) ) deallocate(d_dist)
    allocate(d_dist(0:Nodes-1,4))
    d_dist = 0

    call attach(sparse_pattern,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=no_l, nrows_g=no_u)
    ! Assign starting / ending
    d_dist(Node,1) = index_local_to_global(dit,1,Node)
    d_dist(Node,2) = index_local_to_global(dit,no_l,Node)

    ! allocate space to capture content
    N = maxval(l_ncol)
    allocate(tmp(N))
    ! find the local update size
    N = 0
    do io = 1 , no_l

       ind = l_ptr(io)
       jo  = l_ncol(io)
       tmp(1:jo) = UCORB(l_col(ind+1:ind+jo),no_u)
       N = N + UNIQC(tmp(1:jo))

    end do
    d_dist(Node,3) = N
    deallocate(tmp)

    allocate(tmp(3*Nodes))

    call MPI_AllReduce(d_dist(1,1),tmp,3*Nodes, &
         MPI_Integer, MPI_Sum, MPI_Comm_World, MPIerror)

    ! insert the values again
    ind = 0
    do jo = 1 , 3
       do io = 0 , Nodes - 1
          ind = ind + 1
          d_dist(io,jo) = tmp(ind)
       end do
    end do

    deallocate(tmp)
    
    ! create the stride in data
    do io = 1 , Nodes - 1
       d_dist(io,4) = d_dist(io-1,4) + d_dist(io-1,3)
    end do

  end subroutine ts_init_distribution

#endif
     
end module m_ts_sparse
