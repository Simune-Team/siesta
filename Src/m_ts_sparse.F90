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

! TODO.
!   I have a few ideas which could be brought into existance on top
!   of this version (quite easily, the reason for not doing it
!   now is that the code will change so much that no-one will be familiar
!   with the code). I have hence decided not to introduce the "full"
!   memory reduction scheme for this version.
!
!  1. Many electrode dependent arrays could be formed in a sparse format 
!  2. The codes for collecting charges in regions have not been
!     transferred. This should be do-able (and quite easily)
!     However, it would be MUCH better if this information
!     is only shown per-request.
!  3. Consider splitting Transiesta up in Gamma and non-Gamma
!     codes. This will probably dublicate some code, however,
!     I think it will help the unexperienced in delving into the 
!     wonderful art of Transiestaing... :)
!  
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
       Gamma,block_dist,sparse_pattern,na_u,lasto)

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
#endif
!    use m_monitor

! **********************
! * INPUT variables    *
! **********************
    character(len=*), intent(in) :: slabel
    logical, intent(in) :: IsVolt ! bias calculation
    integer, intent(in) :: N_Elec
    type(Elec), intent(inout) :: Elecs(N_Elec)
    ! A Gamma-calculation?
    logical, intent(in)  :: Gamma
    ! The distribution for the sparsity-pattern
    type(OrbitalDistribution), intent(inout) :: block_dist
    ! SIESTA local sparse pattern (not changed)
    type(Sparsity), intent(inout) :: sparse_pattern
    ! Number of atoms in the unit-cell
    integer, intent(in) :: na_u
    ! Last orbital of the equivalent unit-cell atom
    integer, intent(in) :: lasto(0:na_u)

! **********************
! * LOCAL variables    *
! **********************
    type(OrbitalDistribution) :: dit
    ! Temporary arrays for knowing the electrode size
    logical :: bool
    integer :: no_u_TS

    ! Number of orbitals in TranSIESTA
    no_u_TS = nrows_g(sparse_pattern) - no_Buf

    ! Do a crude check of the sizes
    if ( no_u_TS <= sum(TotUsedOrbs(Elecs)) ) then
       call die("The contact region size is &
            &smaller than the electrode size. Please correct.")
    end if

    if ( IsVolt ) then 

       ! Create the update region (a direct subset of the local sparsity pattern)
       ! Hence it is still a local sparsity pattern.
       call ts_Sparsity_Update(block_dist,sparse_pattern, Elecs, &
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
    call ts_Sparsity_Global(block_dist,sparse_pattern, &
         nnzs(sparse_pattern), Elecs, &
         ts_sp_uc)

    ! The update sparsity pattern can be simplied to the H,S sparsity
    ! pattern, if all electrodes have certain options to be the same.
    bool = all(Elecs(:)%Bulk) .and. all(Elecs(:)%DM_CrossTerms)
    bool = bool .or. all( .not. Elecs(:)%Bulk )

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
          write(*,'(/,a)') 'Created the TranSIESTA H,S sparsity pattern:'
          write(*,'(a)') 'TranSIESTA global update sparsity pattern same as H,S'
          call print_type(ts_sp_uc)
       end if
    else
       ! Create the update region (a direct subset of ts_sp_uc)
       call ts_Sparsity_Update(dit,ts_sp_uc, Elecs, &
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

  end subroutine ts_sparse_init


! Returns the global sparsity pattern for the transiesta region
! Note that this will automatically detect whether there are
! cross connections from left-right, AND remove any
! z-connections (less orbitals to move about, and they should
! not exist).
  subroutine ts_Sparsity_Global(block_dist,s_sp,n_nzs, &
       Elecs, &
       ts_sp)

    use parallel, only : IONode
#ifdef MPI
    use m_glob_sparse
    use mpi_siesta
#endif
    use class_OrbitalDistribution
    use create_Sparsity_SC
    use m_ts_electype
    use m_ts_method
#ifdef TRANSIESTA_DEBUG
    use m_ts_debug
#endif

! **********************
! * INPUT variables    *
! **********************
    ! The SIESTA distribution of the sparsity pattern
    type(OrbitalDistribution), intent(inout) :: block_dist
    ! Sparsity patterns of SIESTA (local)
    type(Sparsity), intent(inout) :: s_sp
    ! The number of non-zeroes in the sparsity pattern (local)
    integer, intent(in) :: n_nzs
    ! All the electrodes
    type(Elec), intent(in) :: Elecs(:)
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
#ifdef MPI
    ! to create the global sparsity pattern (TranSIESTA)
    integer, allocatable :: l_ncolg(:)
    integer, allocatable :: l_ptrg(:)
    integer, allocatable :: l_colg(:)
#endif

    ! Also used in non-MPI (to reduce dublicate code)
    integer :: no_local, no_u, uc_n_nzs, n_nzsg

    ! search logical to determine the update region...
    logical, allocatable :: l_HS(:)

    ! Loop-counters
    integer :: j ,io, jo, ind
    integer :: iot, jot
    logical :: UseBulk

    ! Initialize
    call delete(ts_sp)

    call attach(s_sp,nrows=no_local,nrows_g=no_u)

    ! Create the (local) SIESTA-UC sparsity...
#ifdef MPI
    call crtSparsity_SC(s_sp,sp_global, UC=.TRUE.)
    uc_n_nzs = nnzs(sp_global)
#else
    call crtSparsity_SC(s_sp,sp_uc    , UC=.TRUE.)
    uc_n_nzs = nnzs(sp_uc)
#endif

#ifdef MPI
    ! point to the local (SIESTA-UC) sparsity pattern arrays
    call attach(sp_global,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)
    call glob_sparse_numh(no_local,no_u,l_ncol,l_ncolg)
    call glob_sparse_listhptr(no_u,l_ncolg,l_ptrg)
    call glob_sparse_listh(no_local,no_u, uc_n_nzs, &
         l_ncol , l_ptr , l_col , &
         l_ncolg, l_ptrg, n_nzsg, l_colg)

    ! Delete the local UC sparsity pattern
    call delete(sp_global)
    
    ! Create the globalized UC sparsity pattern
    call newSparsity(sp_uc,no_u, no_u, &
         n_nzsg, l_ncolg, l_ptrg, l_colg, &
         name='SIESTA UC sparsity')
    
    ! Deallocate the arrays which we do not need
    deallocate(l_ncolg,l_ptrg,l_colg)
    call memory('D','I',no_u*2+n_nzsg,'globArrays')

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

       ! Loop number of entries in the row...
       do j = 1 , l_ncol(io)

          ! The index in the pointer array is retrieved
          ind = l_ptr(io) + j
          
          ! The unit-cell column index (remember we are looping a UC SP)
          jo = l_col(ind)

          ! If we are in the buffer region, cycle (l_HS(ind) =.false. already)
          ! note, that we have already *checked* ic
          jot = orb_type(jo)
          if ( jot == TYP_BUFFER ) cycle

          if ( iot > 0 ) then
             UseBulk = Elecs(iot)%Bulk
          else if ( jot > 0 ) then
             UseBulk = Elecs(jot)%Bulk
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
  subroutine ts_Sparsity_Update(dit,s_sp, Elecs, &
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
    type(Elec), intent(in) :: Elecs(:)
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
    integer :: no_local, no_u, n_nzs

    ! search logical to determine the update region...
    logical, allocatable :: lup_DM(:)
    logical :: direct_LR
    logical :: UseBulk, DM_CrossTerms

    ! Loop-counters
    integer :: io, ic, ind
    integer :: ict, jct
    
    ! Logical for determining the region
    logical :: i_in_C, j_in_C

    ! Initialize the tsup
    call delete(tsup_sp)

    call attach(s_sp,nrows=no_local,nrows_g=no_u,nnzs=n_nzs)

    ! allocate space for the MASK to create the TranSIESTA UPDATE region
    allocate(lup_DM(n_nzs))
    call memory('A','L',n_nzs,'transiesta')

    ! Initialize
    lup_DM(:) = .false.
    direct_LR = .false.

    call attach(s_sp,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)

    ! We do not need to check the buffer regions...
    ! We know they will do NOTHING! :)
    do io = 1 , no_local

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

          if ( ict > 0 ) then
             UseBulk = Elecs(ict)%Bulk
             DM_CrossTerms = Elecs(ict)%DM_CrossTerms
          else if ( jct > 0 ) then
             UseBulk = Elecs(jct)%Bulk
             DM_CrossTerms = Elecs(jct)%DM_CrossTerms
          else
             ! we are definitely not in an electrode
             ! just set it to be updated
             UseBulk = .false.
             DM_CrossTerms = .true.
          end if

          ! We check whether it is electrode-connections. 
          ! If, so, they are not used in transiesta:
          if      ( all((/ict,jct/) > 0) ) then
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
       write(*,*) 'WARNING: Consider increasing the junction &
            &in the z-direction.'
       write(0,*) 'WARNING: Cross connections across the &
            &junction leads to tunneling.'
       write(0,*) 'WARNING: Transiesta will disregard these &
            &cross-terms'
       write(0,*) 'WARNING: Consider increasing the junction &
            &in the z-direction.'
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
    integer, pointer :: l_ncol(:)   => null()
    integer, pointer :: l_ptr(:)    => null()
    integer, pointer :: l_col(:)    => null()
    integer, pointer :: sub_ncol(:) => null()
    integer, pointer :: sub_ptr(:)  => null()
    integer, pointer :: sub_col(:)  => null()
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
    integer :: in_elec1, in_elec2, N_buff
    integer :: afrom, ato

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
      offset = sum(TotUsedAtoms(Elecs(:)), MASK=Elecs(:)%idx_na < ia )
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

    integer :: io, jo, ind, no_local, no_u, N

    integer :: MPIerror

    ! Allocate all the nodes distribution sizes
    if ( allocated(d_dist) ) deallocate(d_dist)
    allocate(d_dist(0:Nodes-1,4))
    d_dist = 0

    call attach(sparse_pattern,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=no_local, nrows_g=no_u)
    ! Assign starting / ending
    d_dist(Node,1) = index_local_to_global(dit,1,Node)
    d_dist(Node,2) = index_local_to_global(dit,no_local,Node)

    ! allocate space to capture content
    N = maxval(l_ncol)
    allocate(tmp(N))
    ! find the local update size
    N = 0
    do io = 1 , no_local

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
