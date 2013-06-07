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
  type(Sparsity), save :: ts_sp_uc ! TS-GLOBAL only UC

  ! We will save the "update region sparsity"
  ! Note that this is NOT the same as the sparsity pattern
  ! provided by SIESTA!
  ! After using this sparsity pattern, we do not need the listud[g] arrays
  ! as this sparsity pattern is a reflection of the mask used by listud[g]
  ! Lastly the usage of a MASKED sparsity pattern, will reduce the clutter
  ! between MPI and non-MPI codes as it will be the same array in both 
  ! circumstances.
  type(Sparsity), save :: tsup_sp_uc ! TS-update-GLOBAL only UC

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

contains

! This routine setups what-ever is needed to do the
! memory reduced TranSIESTA code.
! This means collecting information about which region needs
! update, etc.
  subroutine ts_sparse_init(Gamma,block_dist,sparse_pattern,na_u,lasto)

    use class_OrbitalDistribution

#ifdef MPI
    use mpi_siesta, only : MPI_Comm_Self
#endif 
    use parallel, only: IONode, Node
    ! Used when creating the update region...
    use m_ts_options, only : UseBulk, UpdateDMCR

    use m_ts_electype
    use m_ts_options, only : IsVolt
    use m_ts_options, only : ElLeft, ElRight
    use m_ts_options, only : na_BufL, no_BufL
    use m_ts_options, only : na_BufR, no_BufR

! **********************
! * INPUT variables    *
! **********************
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
    integer :: no_L, no_R
    integer :: no_u_LCR
    ! Calculate the number of used atoms/orbitals in left/right
    no_L = TotUsedOrbs(ElLeft)
    no_R = TotUsedOrbs(ElRight)

    ! Number of orbitals in TranSIESTA
    no_u_LCR = nrows_g(sparse_pattern) - no_BufL - no_BufR

    ! Do a crude check of the sizes
    if ( no_u_LCR <= no_L + no_R ) then
       call die("The contact region size is &
            &smaller than the electrode size. Please correct.")
    end if

    if ( IsVolt ) then    
       ! Create the update region (a direct subset of the local sparsity pattern)
       ! Hence it is still a local sparsity pattern.
       call ts_Sparsity_Update(block_dist,sparse_pattern,UseBulk,UpdateDMCR, &
            no_BufL, no_BufR, no_L, no_R, &
            no_u_LCR, ltsup_sp_sc)
       
       if ( IONode ) then
          write(*,'(/,a)') 'Created the TranSIESTA local update sparsity pattern:'
          call print_type(ltsup_sp_sc)
       end if

#ifdef TRANSIESTA_DEBUG
       call sp_to_file(3000+Node,ltsup_sp_sc)
#endif

       ! Create the pointer from the local transiesta update sparsity 
       ! to the local siesta sparsity
       call ts_Sparsity_Subset_pointer(block_dist,sparse_pattern,ltsup_sp_sc, &
            ltsup_sc_pnt)
    end if

    ! Create the global transiesta H(k), S(k) sparsity pattern
    call ts_Sparsity_Global(block_dist,sparse_pattern, &
         nnzs(sparse_pattern), &
         UseBulk, UpdateDMCR, &
         no_BufL,no_BufR,no_L,no_R, &
         no_u_LCR, &
         ts_sp_uc)

    if ( IONode ) then
       write(*,'(/,a)') 'Created the TranSIESTA H,S sparsity pattern:'
       call print_type(ts_sp_uc)
    end if

#ifdef TRANSIESTA_DEBUG
    call sp_to_file(1000+Node,ts_sp_uc)
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

    ! Create the update region (a direct subset of ts_sp_uc)
    call ts_Sparsity_Update(dit,ts_sp_uc,UseBulk,UpdateDMCR, &
         no_BufL, no_BufR, no_L, no_R, &
         no_u_LCR, tsup_sp_uc)

    if ( IONode ) then
       write(*,'(/,a)') 'Created the TranSIESTA global update sparsity pattern:'
       call print_type(tsup_sp_uc)
    end if


#ifdef TRANSIESTA_DEBUG
    call print_type(tsup_sp_uc)
    call sp_to_file(2000+Node,tsup_sp_uc)
#endif


#ifdef TRANSIESTA_DEBUG
    call die('')
#endif
        
  end subroutine ts_sparse_init


! Returns the global sparsity pattern for the transiesta region
! Note that this will automatically detect whether there are
! cross connections from left-right, AND remove any
! z-connections (less orbitals to move about, and they should
! not exist).
  subroutine ts_Sparsity_Global(block_dist,s_sp,n_nzs, &
       UseBulk, UpdateDMCR, &
       no_BufL,no_BufR,no_L,no_R, &
       no_u_LCR, &
       ts_sp)

    use parallel, only : IONode, Node
#ifdef MPI
    use m_glob_sparse
    use mpi_siesta
#endif
    use class_OrbitalDistribution
    use create_Sparsity_SC

! **********************
! * INPUT variables    *
! **********************
    ! If we have a Gamma calculation
    logical, intent(in) :: UseBulk, UpdateDMCR
    ! Number of orbitals in each segment
    integer, intent(in) :: no_BufL, no_BufR, no_L, no_R
    ! Number of rows in the transiesta SP
    integer, intent(in) :: no_u_LCR
    ! The SIESTA distribution of the sparsity pattern
    type(OrbitalDistribution), intent(inout) :: block_dist
    ! Sparsity patterns of SIESTA and the returned update region.    
    ! SIESTA sparsity pattern is the local one...
    type(Sparsity), intent(inout) :: s_sp, ts_sp
    ! The number of non-zeroes in the sparsity pattern (local)
    integer, intent(in) :: n_nzs

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
    logical, allocatable :: lup_DM(:)
    logical :: direct_LR

    ! Loop-counters
    integer :: i,j ,io,jo, ic,jc, ind
    integer :: lio, ljc, lind, no_u_LC

    ! Logical for determining the region
    logical :: i_in_C, j_in_C

    ! Initialize
    call delete(ts_sp)

    call retrieve(s_sp,nrows=no_local,nrows_g=no_u)

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
    call retrieve(sp_global,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)
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
    write(*,*)'Created UC SIESTA sparsity'
    call print_type(sp_uc)
    call sp_to_file(400+Node, sp_uc)
#endif

    ! Write that we have created it
    if ( IONode ) call print_type(sp_uc)

    ! Now we have the globalized SIESTA Unit-cell pattern, make it 
    ! only to the Transiesta sparsity pattern

    ! Immediately point the global arrays to their respective parts
    call retrieve(sp_uc,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nnzs=n_nzsg)

    ! allocate space for the MASK to create the TranSIESTA GLOBAL region
    allocate(lup_DM(n_nzsg))
    call memory('A','L',n_nzsg,'transiesta')

    ! Initialize
    lup_DM(:) = .false.
    direct_LR = .false.

    ! Retrieve the ending point of LC region
    no_u_LC = no_u_LCR - no_R

    ! We do not need to check the buffer regions...
    ! We know they will do NOTHING! :)
    do io = no_BufL + 1 , no_u - no_BufR
       ! Shift out of the buffer region
       ic = io - no_BufL
       ! Loop number of entries in the row...
       do j = 1 , l_ncol(io)

          ! The index in the pointer array is retrieved
          ind = l_ptr(io) + j
          
          ! The unit-cell column index (remember we are looping a UC SP)
          jc = l_col(ind) - no_BufL

          ! If we are in the buffer region, cycle (lup_DM(ind) =.false. already)
          ! note, that we have already *checked* ic
          if ( jc < 1 ) cycle
          if ( no_u_LCR < jc ) cycle

          ! We check whether it is electrode-connections. 
          ! If, so, they are not used in transiesta:
          if      ( ic <= no_L .and. no_u_LC < jc ) then
             lup_DM(ind) = .false.

             ! This means that we have an INNER-cell connection
             ! TODO, do a check on the local nodes SP for this
             !direct_LR = direct_LR .or. ( jc == jo - no_BufL )
          else if ( jc <= no_L .and. no_u_LC < ic ) then
             lup_DM(ind) = .false.

             ! This means that we have an INNER-cell connection
             ! TODO, do a check on the local nodes SP for this
             !direct_LR = direct_LR .or. ( jc == jo - no_BufL )
          else if ( UseBulk .and. ic <= no_L .and. jc <= no_L ) then
             lup_DM(ind) = .false.

             ! in case of UseBulk we remove the left electrode from
             ! the Hamiltonian
          else if ( UseBulk .and. no_u_LC < ic .and. no_u_LC < jc ) then
             lup_DM(ind) = .false.

             ! in case of UseBulk we remove the right electrode from
             ! the Hamiltonian
          else
             lup_DM(ind) = .true.
          end if

       end do

    end do

    ! We now have a MASK of the actual needed TranSIESTA sparsity pattern
    ! We create the TranSIESTA sparsity pattern
    call crtSparsity_SC(sp_uc,ts_sp,MASK=lup_DM)

    ! clean-up
    deallocate(lup_DM)
    call memory('D','L',n_nzsg,'transiesta')
    
    ! Furthermore, we dont need the SIESTA UC sparsity global...
    ! TODO check that this object is in fact deleted...
    call delete(sp_uc)

  end subroutine ts_Sparsity_Global

#ifdef TRANSIESTA_DEBUG
  
  subroutine sp_to_file(u,sp)
    use geom_helper, only : UCORB
#ifdef MPI
    use mpi_siesta, only : MPI_Comm_World
#endif
    integer, intent(in) :: u
    type(Sparsity), intent(inout) :: sp
    integer :: io,jo,j,ind
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    call retrieve(sp,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)

    write(u,'(i5)') nrows(sp)

    do io = 1 , nrows(sp)
       if ( l_ncol(io) == 0 ) cycle
       do j = 1 , l_ncol(io)
          ind = l_ptr(io) + j
          jo = UCORB(l_col(ind),nrows(sp))
          write(u,'(3(tr1,i5))') io,jo,1
       end do
    end do
    if ( ind /= nnzs(sp) ) then
       call die('Have not looped through all things')
    end if
#ifdef MPI
    call MPI_Barrier(MPI_Comm_World,io)
#endif

  end subroutine sp_to_file

#endif


! Returns the global sparsity pattern for the transiesta region
! Note that this will automatically detect whether there are
! cross connections from left-right, AND remove any
! z-connections (less orbitals to move about, and they should
! not exist).
  subroutine ts_Sparsity_Update(dit,s_sp, &
       UseBulk, UpdateDMCR, &
       no_BufL,no_BufR,no_L,no_R, &
       no_u_LCR, &
       tsup_sp)

    use parallel, only : IONode, Node
    use geom_helper, only : UCORB
    use create_Sparsity_SC
    use class_OrbitalDistribution

! **********************
! * INPUT variables    *
! **********************
    type(OrbitalDistribution), intent(inout) :: dit
    type(Sparsity), intent(inout) :: s_sp
    ! Options for creating the update
    logical, intent(in) :: UseBulk, UpdateDMCR
    ! Number of orbitals in each segment
    integer, intent(in) :: no_BufL, no_BufR, no_L, no_R
    ! Number of rows in the transiesta SP
    integer, intent(in) :: no_u_LCR
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

    ! Loop-counters
    integer :: j ,io,jo, ic,jc, ind
    integer :: lio, ljc, lind, no_u_LC

    ! Logical for determining the region
    logical :: i_in_C, j_in_C

    ! Initialize the tsup
    call delete(tsup_sp)

    call retrieve(s_sp,nrows=no_local,nrows_g=no_u,nnzs=n_nzs)

    ! allocate space for the MASK to create the TranSIESTA UPDATE region
    allocate(lup_DM(n_nzs))
    call memory('A','L',n_nzs,'transiesta')

    ! Initialize
    lup_DM(:) = .false.
    direct_LR = .false.

    ! Retrieve the ending orbital of the LC region
    no_u_LC = no_u_LCR - no_R

    call retrieve(s_sp,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)

    ! We do not need to check the buffer regions...
    ! We know they will do NOTHING! :)
    do io = 1 , no_local
       ! Shift out of the buffer region
       ic = index_local_to_global(dit,io,Node) - no_BufL

       ! If we are in the buffer region, cycle (lup_DM(ind) =.false. already)
       if ( ic < 1 ) cycle
       if ( no_u_LCR < ic ) cycle

       ! Loop number of entries in the row...
       do j = 1 , l_ncol(io)

          ! The index in the pointer array is retrieved
          ind = l_ptr(io) + j
          
          ! The unit-cell column index, without buffer
          jc = UCORB(l_col(ind),no_u) - no_BufL

          ! If we are in the buffer region, cycle (lup_DM(ind) =.false. already)
          ! note, that we have already *checked* ic
          if ( jc < 1 ) cycle
          if ( no_u_LCR < jc ) cycle

          ! We check whether it is electrode-connections. 
          ! If, so, they are not used in transiesta:
          if      ( ic <= no_L .and. no_u_LC < jc ) then
             lup_DM(ind) = .false.

             ! This means that we have an INNER-cell connection
             ! TODO, do a check on the local nodes SP for this
             !direct_LR = direct_LR .or. ( jc == jo - no_BufL )
          else if ( jc <= no_L .and. no_u_LC < ic ) then
             lup_DM(ind) = .false.

             ! This means that we have an INNER-cell connection
             ! TODO, do a check on the local nodes SP for this
             !direct_LR = direct_LR .or. ( jc == jo - no_BufL )
          else if ( UseBulk ) then
             ! If UseBulk, we can safely update cross-terms
             ! between the central and the electrodes (the self-energies are the
             ! same)
             ! However, cross-term updates are controlled by the user 
             ! via the option UpdateDMCR
             ! This will also remove any electrode terms

             i_in_C = ( no_L < ic .and. ic <= no_u_LC )
             j_in_C = ( no_L < jc .and. jc <= no_u_LC )

             if ( UpdateDMCR ) then
                ! the user has requested to ONLY update the central region
                lup_DM(ind) = i_in_C .and. j_in_C
             else
                ! the user has requested to also update cross-terms
                lup_DM(ind) = i_in_C .or.  j_in_C
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

    integer :: no_l, io, j, jo, sub_ind, ind

    call retrieve(sp,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=no_l)
    call retrieve(sub_sp,n_col=sub_ncol,list_ptr=sub_ptr,list_col=sub_col, &
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
    
end module m_ts_sparse
