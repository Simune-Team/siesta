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
!  3. Adding a charge calculation when TranSIESTA ends might be
!     the best, and most effective way out of this (while still
!     retaining some information for the user).
!     Mainly because it is then performed in parallel.
!  4. Consider splitting Transiesta up in Gamma and non-Gamma
!     calculation. This will probably dublicate some code, however,
!     I think it will help the unexperienced in delving into the 
!     wonderful art of Transiestaing... :)
!  
module m_ts_mem_sparsity

  use class_Sparsity

  implicit none

  integer, parameter :: dp = selected_real_kind(10,100)

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

  ! Logical variable that describes the solution method on
  ! LEFT-RIGHT-EQUILIBRIUM contour points.
  ! This is realized by the fact that for:
  !    UseBulk .and. UpdateDMCR
  ! the needed part of the GF is only the C...C regions:
  !
  !  -------------------------------------
  !  | L...L | L...C   0     0   |   0   |
  !  | C...L | C...C C...C C...C |   0   |
  !  |   0   | C...C C...C C...C |   0   |
  !  |   0   | C...C C...C C...C | C...R |
  !  |   0   |   0     0   R...C | R...R |
  !  -------------------------------------
  ! 
  ! This means we can solve the following instead:
  ! G_F^{-1} G_F I_P = I \times I_P,
  ! where I_P:
  !  ---------------------
  !  |   0     0     0   |
  !  |   1     0     0   |
  !  |   0     1     0   |
  !  |   0     0     1   |
  !  |   0     0     0   |
  !  ---------------------
  ! Note, that this can ONLY be used in EQUI contour points.
  ! In principle we can obtain the EXACT size of the problem
  ! For very large electrodes. This could come in handy.
  logical, save :: GF_INV_EQUI_PART = .false.

contains

! This routine setups what-ever is needed to do the
! memory reduced TranSIESTA code.
! This means collecting information about which region needs
! update, etc.
  subroutine ts_mem_init(Gamma,sparse_pattern,na_u,lasto)

    ! Used when creating the update region...
    use m_ts_options, only : UseBulk, UpdateDMCR

    ! left stuff
    use m_ts_options, only : na_BufL => NBufAtL
    !use m_ts_options, only : na_L_HS => NUsedAtomsL
    use m_ts_options, only : no_L_HS => NUsedOrbsL
    use m_ts_options, only : NRepA1L, NRepA2L
    ! right stuff
    use m_ts_options, only : na_BufR => NBufAtR
    !use m_ts_options, only : na_R_HS => NUsedAtomsR
    use m_ts_options, only : no_R_HS => NUsedOrbsR
    use m_ts_options, only : NRepA1R, NRepA2R

! **********************
! * INPUT variables    *
! **********************
    ! A Gamma-calculation?
    logical, intent(in)  :: Gamma
    ! SIESTA local sparse pattern (not changed)
    type(Sparsity), intent(inout) :: sparse_pattern
    ! Number of atoms in the unit-cell
    integer, intent(in) :: na_u
    ! Last orbital of the equivalent unit-cell atom
    integer, intent(in) :: lasto(0:na_u)

! **********************
! * LOCAL variables    *
! **********************
    ! Temporary arrays for knowing the electrode size
    integer :: no_L, no_R
    integer :: no_BufL, no_BufR
    integer :: no_u_LCR

    ! Calculate the number of used atoms/orbitals in left/right
    !na_L = na_L_HS * NRepA1L * NRepA2L
    no_L = no_L_HS * NRepA1L * NRepA2L
    !na_R = na_R_HS * NRepA1R * NRepA2R
    no_R = no_R_HS * NRepA1R * NRepA2R

    ! Calculate the number of orbitals not used (i.e. those 
    ! in the buffer regions)
    ! Left has the first atoms
    no_BufL = lasto(na_BufL)
    ! Right has the last atoms
    no_BufR = lasto(na_u) - lasto(na_u - na_BufR)

    ! Number of orbitals in TranSIESTA
    no_u_LCR = nrows_g(sparse_pattern) - no_BufL - no_BufR

    ! Do a crude check of the sizes
    if ( no_u_LCR <= no_L + no_R ) then
       call die("The contact region size is &
            &smaller than the electrode size. Please correct.")
    end if
    
    ! Setup the correct handling of EQUILIBRIUM solution method:
    ! See above the global variable for its use.
    GF_INV_EQUI_PART = UseBulk .and. UpdateDMCR

    ! We make sure that the arrays are deleted, before entrance
    call delete(ts_sp_uc)
    call delete(tsup_sp_uc)

    ! Create the sparsity patterns we need...
    call ts_Sparsity_Global(Gamma,sparse_pattern, nnzs(sparse_pattern), &
         UseBulk, UpdateDMCR, &
         no_BufL,no_BufR,no_L,no_R, &
         no_u_LCR, &
         ts_sp_uc,tsup_sp_uc)
        
  end subroutine ts_mem_init


! Returns the global sparsity pattern for the transiesta region
! Note that this will automatically detect whether there are
! cross connections from left-right, AND remove any
! z-connections (less orbitals to move about, and they should
! not exist).
  subroutine ts_Sparsity_Global(Gamma,s_sp,maxnh, &
       UseBulk, UpdateDMCR, &
       no_BufL,no_BufR,no_L,no_R, &
       no_u_LCR, &
       ts_sp,tsup_sp)

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
    logical, intent(in) :: Gamma, UseBulk, UpdateDMCR
    ! Number of orbitals in each segment
    integer, intent(in) :: no_BufL, no_BufR, no_L, no_R
    ! Number of rows in the transiesta SP
    integer, intent(in) :: no_u_LCR
    ! Sparsity patterns of SIESTA and the returned update region.    
    ! SIESTA sparsity pattern is the local one...
    type(Sparsity), intent(inout) :: s_sp, ts_sp, tsup_sp
    ! The number of non-zeroes in the sparsity pattern (local)
    integer, intent(in) :: maxnh

! **********************
! * LOCAL variables    *
! **********************
    ! Helpers to generate the sparsity patterns
    type(OrbitalDistribution) :: dit
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
    integer :: no_local, no_u, maxnhg

    ! search logical to determine the update region...
    logical, allocatable :: lup_DM(:)
    logical :: direct_LR

    ! Loop-counters
    integer :: i,j ,io,jo, ic,jc, ind
    integer :: lio, ljc, lind

    ! Logical for determining the region
    logical :: i_in_C, j_in_C

    no_local = nrows  (s_sp)
    no_u     = nrows_g(s_sp)

    ! Create the Fake distribution
    ! The Block-size is the number of orbitals, i.e. all on the first processor
#ifdef MPI
    call newDistribution(no_u,MPI_COMM_WORLD,dit,name='SIESTA global distribution')
#else
    call newDistribution(no_u,-1,dit,name='SIESTA global distribution')
#endif

#ifdef TRANSIESTA_DEBUG
    write(*,*)'Starting new Distribution'
#endif

    ! point to the local sparsity pattern arrays
    l_ncol   => n_col   (s_sp)
    l_ptr    => list_ptr(s_sp)
    l_col    => list_col(s_sp)

#ifdef MPI
    call glob_sparse_numh(no_local,no_u,l_ncol,l_ncolg)
    call glob_sparse_listhptr(no_u,l_ncolg,l_ptrg)
    call glob_sparse_listh(no_local,no_u, maxnh, &
         l_ncol , l_ptr , l_col , &
         l_ncolg, l_ptrg, maxnhg, l_colg)

    call newSparsity(sp_global,no_u, no_u, &
         maxnhg, l_ncolg, l_ptrg, l_colg, &
         name='SIESTA-full sparsity')

    ! Deallocate the arrays which we do not need
    deallocate(l_ncolg,l_ptrg,l_colg)
    call memory('D','I',no_u*2+maxnhg,'globArrays')

#ifdef TRANSIESTA_DEBUG
    write(*,*)'Created FULL SIESTA sparsity'
    call print_type(sp_global)

    call sp_to_file(400+Node, sp_global)
#endif

    ! Create the SIESTA-UC sparsity...
    call crtSparsity_SC(sp_global,sp_uc, &
         UC=.TRUE.)

#ifdef TRANSIESTA_DEBUG
    write(*,*)'Created UC SIESTA sparsity'
    call print_type(sp_uc)

    call sp_to_file(500+Node, sp_uc)
#endif

    ! Delete the full sparsity pattern
    ! TODO check that it IS deleted
    call delete(sp_global)

    ! Write that we have created it
    if ( IONode ) call print_type(sp_uc)

#else 
    call crtSparsity_SC(s_sp,sp_uc, UC=.TRUE.)
#endif

    ! Immediately point the global arrays to their respective parts
    l_ncol => n_col   (sp_uc)
    l_ptr  => list_ptr(sp_uc)
    l_col  => list_col(sp_uc)

    ! allocate space for the MASK to create the TranSIESTA GLOBAL region
    maxnhg = nnzs(sp_uc)
    allocate(lup_DM(maxnhg))
    call memory('A','L',maxnhg,'transiesta')

    ! Initialize
    lup_DM(:) = .false.
    direct_LR = .false.

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
          if      ( ic <= no_L .and. no_u_LCR - no_R < jc ) then
             lup_DM(ind) = .false.

             ! This means that we have an INNER-cell connection
             ! TODO, do a check on the local nodes SP for this
             !direct_LR = direct_LR .or. ( jc == jo - no_BufL )
          else if ( jc <= no_L .and. no_u_LCR - no_R < ic ) then
             lup_DM(ind) = .false.

             ! This means that we have an INNER-cell connection
             ! TODO, do a check on the local nodes SP for this
             !direct_LR = direct_LR .or. ( jc == jo - no_BufL )
          else
             lup_DM(ind) = .true.
          end if

       end do

    end do

#ifdef TRANSIESTA_DEBUG
    write(*,*)'Created the update region'
    write(*,'(10000(tr1,l1))') lup_DM
#endif

    ! TODO, the above dissemation of the sparsity
    ! can be performed in parallel, however, it requires
    ! manual updates of the lup_DM array afterwards.
    ! This is easier, and will only be performed once!
    
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
    call crtSparsity_SC(sp_uc,ts_sp,MASK=lup_DM)

    if ( IONode ) then
       write(*,'(/,a)') 'Created the TranSIESTA sparsity pattern:'
       call print_type(ts_sp)
    end if

    ! clean-up
    deallocate(lup_DM)
    call memory('D','L',maxnhg,'transiesta')
    
    ! Furthermore, we dont need the SIESTA UC sparsity global...
    ! TODO check that this object is in fact deleted...
    call delete(sp_uc)

    ! obtain the TranSIESTA GLOBAL sparsity
    l_ncol => n_col   (ts_sp)
    l_ptr  => list_ptr(ts_sp)
    l_col  => list_col(ts_sp)

    ! Reallocate the MASK to figure out the update sparsity pattern
    maxnhg = nnzs(ts_sp)
    allocate(lup_DM(maxnhg))
    call memory('A','L',maxnhg,'transiesta')

    ! Now we can create the TranSIESTA update region
    ! As the Transiesta sparsity pattern is contained in the scheme
    ! of buffer-atoms, we still need to skip the buffer atoms.
    ! keep in mind that Transiesta sparsity pattern is a subset of SIESTA
    ! pattern, WITHOUT reindexing the rows and columns.
    ! Also note that this sparsity pattern is a direct subset of
    ! the transiesta sparsity pattern. 
    ! Also, it will ALWAYS be smaller as the intra-electrode cross-terms
    ! will never be updated.
    do io = no_BufL + 1 , no_u - no_BufR
       ! Shift out of the buffer region
       ic = io - no_BufL
       ! Loop number of entries in the row...
       do j = 1 , l_ncol(io)

          ! The index in the pointer array is retrieved
          ind = l_ptr(io) + j

          ! The unit-cell equivalent index (without buffer-region)
          jc = l_col(ind) - no_BufL

          if ( jc < 1 ) call die('Transiesta sparsity pattern is &
               &faulty. Please check with the developers.')
          if ( no_u_LCR < jc ) call die('Transiesta sparsity pattern is &
               &faulty. Please check with the developers.')
          
          if ( UseBulk ) then
             ! If UseBulk, we can safely update cross-terms
             ! between the central and the electrodes (the self-energies are the
             ! same)
             ! However, cross-term updates are controlled by the user 
             ! via the option UpdateDMCR

             i_in_C = ( no_L < ic .and. ic <= no_u_LCR - no_R )
             j_in_C = ( no_L < jc .and. jc <= no_u_LCR - no_R )

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
             ! TODO why can't we use UpdateDMCR in case of not using UseBulk?

          end if

          ! if it is a cross-junction term (or z-connection term)
          if ( ic <= no_L .and. no_u_LCR - no_R < jc ) then
             call die('Transiesta sparsity pattern is &
                  &faulty. Please check with the developers.')
          elseif ( jc <= no_L .and. no_u_LCR - no_R < ic ) then
             call die('Transiesta sparsity pattern is &
                  &faulty. Please check with the developers.')
          end if

       end do

    end do

    ! We need to create a unit-cell only sparsity pattern
    ! for the transiesta update region. 
    ! This may seem ludicris to have 2 sparsity patterns...
    ! However, this should save more space than they consume... 
    call crtSparsity_SC(ts_sp,tsup_sp,MASK=lup_DM)
    
    if ( IONode ) then
       write(*,'(a)') 'Created the TranSIESTA update sparsity pattern:'
       call print_type(tsup_sp)
       write(*,*) '' ! newline...
    end if

    ! Clean up
    call memory('D','L',maxnhg,'transiesta')
    deallocate(lup_DM)

#ifdef TRANSIESTA_DEBUG
    call sp_to_file(1000+Node,ts_sp)
    call sp_to_file(2000+Node,tsup_sp)

    call die('')
  contains
    
    subroutine sp_to_file(u,sp)
      use geom_helper, only : UCORB
      integer, intent(in) :: u
      type(Sparsity), intent(inout) :: sp
      integer :: io,jo,j,ind
      integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
      l_ncol => n_col   (sp)
      l_ptr  => list_ptr(sp)
      l_col  => list_col(sp)
      
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
  end subroutine ts_Sparsity_Global


  ! A subroutine for printing out the charge distribution in the cell
  ! it will currently only handle the full charge distribution, and
  ! not per k-point.
  subroutine ts_print_charges(dit, sparse_pattern, &
       na_u, lasto, &
       nspin, n_nzs, DM, S)
    ! left stuff
    use m_ts_options, only : na_BufL => NBufAtL
    use m_ts_options, only : no_L_HS => NUsedOrbsL
    use m_ts_options, only : NRepA1L, NRepA2L
    ! right stuff
    use m_ts_options, only : na_BufR => NBufAtR
    use m_ts_options, only : no_R_HS => NUsedOrbsR
    use m_ts_options, only : NRepA1R, NRepA2R
    use m_ts_mem_scat, only : get_scat_region
    use parallel, only : IONode, Node
#ifdef MPI
    use mpi_siesta
#endif
    use class_OrbitalDistribution
    use class_Sparsity
    use geom_helper, only : UCORB

! **********************
! * INPUT variables    *
! **********************
    type(OrbitalDistribution), intent(inout) :: dit
    ! SIESTA local sparse pattern (not changed)
    type(Sparsity), intent(inout) :: sparse_pattern
    ! Number of atoms in the unit-cell
    integer, intent(in) :: na_u
    ! Last orbital of the equivalent unit-cell atom
    integer, intent(in) :: lasto(0:na_u)
    ! Number of non-zero elements
    integer, intent(in) :: nspin, n_nzs
    ! The density matrix and overlap
    real(dp), intent(in) :: DM(n_nzs,nspin), S(n_nzs)

! **********************
! * LOCAL variables    *
! **********************
    real(dp) :: Q(0:9,nspin,2)
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer :: no_lo, no_u, lio, io, ind, jo, ispin, r
    integer :: no_u_TS
    integer :: no_BufL, no_BufR
    integer :: no_L, no_R
#ifdef MPI
    integer :: MPIerror
#endif

    ! Calculate the buffer region and electrode orbitals
    no_L = no_L_HS * NRepA1L * NRepA2L
    no_R = no_R_HS * NRepA1R * NRepA2R

    ! Calculate the number of orbitals not used (i.e. those 
    ! in the buffer regions)
    ! Left has the first atoms
    no_BufL = lasto(na_BufL)
    ! Right has the last atoms
    no_BufR = lasto(na_u) - lasto(na_u - na_BufR)


    no_lo = nrows  (sparse_pattern)
    no_u  = nrows_g(sparse_pattern)
    no_u_TS = no_u - no_BufL - no_BufR

    ! The sparse lists.
    l_ncol => n_col   (sparse_pattern)
    l_ptr  => list_ptr(sparse_pattern)
    l_col  => list_col(sparse_pattern)

    ! Initialize charges
    Q = 0._dp

    do ispin = 1 , nspin
       do lio = 1 , no_lo

          ! obtain the global index of the orbital.
          io = index_local_to_global(dit,lio,Node) - no_BufL

          ! Loop number of entries in the row... (index frame)
          do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
             
             ! as the local sparsity pattern is a super-cell pattern,
             ! we need to check the unit-cell orbital
             ! The unit-cell column index
             jo = UCORB(l_col(ind),no_u) - no_BufL
             
             r = get_scat_region(io,no_L,jo,no_R,no_u_TS)
             
             Q(r,ispin,2) = Q(r,ispin,2) + &
                  DM(ind,ispin) * S(ind)
          end do
       end do
    end do

#ifdef MPI
    call MPI_Reduce(Q(0,1,2),Q(0,1,1),10*nspin,MPI_Double_Precision,MPI_SUM, &
         0,MPI_Comm_World,MPIerror)
#endif

    if ( IONode ) then 
       write(*,'(/,a)') 'transiesta: Charge distribution:'
       if ( nspin > 1 ) then
          write(*,'(a,2(f12.5,tr1))') &
               'Total charge                  [Q]    :', &
               sum(Q(:,1,1)),sum(Q(:,2,1))
          if ( no_BufL > 0 ) write(*,'(a,2(f12.5,tr1),/,a,2(f12.5,tr1))') &
               'Left buffer                   [LB]   :',Q(1,1,1), Q(1,2,1), &
               'Left buffer/left electrode    [LB-L] :',Q(2,1,1), Q(2,2,1)
          write(*,'(a,2(f12.5,tr1),4(/,a,2(f12.5,tr1)))') &
               'Left electrode                [L]    :',Q(3,1,1), Q(3,2,1), &
               'Left electrode/device         [L-C]  :',Q(4,1,1), Q(4,2,1), &
               'Device                        [C]    :',Q(5,1,1), Q(5,2,1), &
               'Device/right electrode        [C-R]  :',Q(6,1,1), Q(6,2,1), &
               'Right electrode               [R]    :',Q(7,1,1), Q(7,2,1)
          if ( no_BufR > 0 ) write(*,'(a,2(f12.5,tr1),/,a,2(f12.5,tr1))') &
               'Right electrode/right buffer  [R-RB] :',Q(8,1,1), Q(8,2,1), &
               'Right buffer                  [RB]   :',Q(9,1,1), Q(9,2,1)
          write(*,'(a,2(f12.5,tr1),/)') &
               'Other                         [O]    :',Q(0,1,1), Q(0,2,1)

       else
          write(*,'(a,f12.5)') &
               'Total charge                  [Q]    :',sum(Q(:,1,1))
          if ( no_BufL > 0 ) write(*,'(a,f12.5,/,a,f12.5)') &
               'Left buffer                   [LB]   :',Q(1,1,1), &
               'Left buffer/left electrode    [LB-L] :',Q(2,1,1)
          write(*,'(a,f12.5,4(/,a,f12.5))') &
               'Left electrode                [L]    :',Q(3,1,1), &
               'Left electrode/device         [L-C]  :',Q(4,1,1), &
               'Device                        [C]    :',Q(5,1,1), &
               'Device/right electrode        [C-R]  :',Q(6,1,1), &
               'Right electrode               [R]    :',Q(7,1,1)
          if ( no_BufR > 0 ) write(*,'(a,f12.5,/,a,f12.5)') &
               'Right electrode/right buffer  [R-RB] :',Q(8,1,1), &
               'Right buffer                  [RB]   :',Q(9,1,1)
          write(*,'(a,f12.5,/)') &
               'Other                         [O]    :',Q(0,1,1)
       end if
    end if
    
  end subroutine ts_print_charges
  
end module m_ts_mem_sparsity
