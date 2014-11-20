! This module contains
! the different regions used in tbtrans

! Coded by Nick Papior Andersen
! 2014

module m_tbt_regions

  use m_region
  use class_OrbitalDistribution
  use class_Sparsity

  implicit none

  save

  ! The full tbtrans sparsity region in the unit-cell equivalent
  ! region, this is used to accomodate the Hamiltonian and overlap
  ! matrices. In principle this could be made obsolete.
  type(Sparsity) :: sp_uc ! TBT-GLOBAL (UC)

  ! the different regions that becomes the electrodes
  type(tRegion), allocatable, target :: r_aEl_alone(:), r_oEl_alone(:)

  ! the different regions that connects to the equivalent
  ! electrodes.
  ! I.e. it is the regions that goes from the electrode
  ! and down to the central region (without any central
  ! region overlap)
  type(tRegion), allocatable, target :: r_aEl(:)   , r_oEl(:)
  type(tRegion), allocatable, target :: r_aElpD(:) , r_oElpD(:)

  ! the device region (the calculated GF region)
  type(tRegion) :: r_aDev, r_oDev

  ! The buffer region, just for completeness
  type(tRegion) :: r_aBuf, r_oBuf

  ! k -> k' regions
  type(tRegion), allocatable :: r_ak(:)

contains

  subroutine tbt_init_regions(N_Elec, Elecs, cell, na_u, lasto, dit, sp, &
       nsc, isc_off)

    use fdf
    use fdf_extra
    use parallel, only : Node, Nodes
#ifdef MPI
    use mpi_siesta, only : MPI_Comm_World
#endif

    use intrinsic_missing, only : SPC_PROJ, VNORM, VEC_PROJ
    use m_ts_electype
    use m_ts_method, only : ts_init_regions
    use m_ts_method, only : atom_type, TYP_DEVICE, TYP_BUFFER

    use m_ts_sparse, only : ts_Sparsity_Global

#ifdef TRANSIESTA_DEBUG
    use m_ts_debug
#endif

    use m_sparsity_handling
    
    ! Number of electrodes
    integer, intent(in) :: N_Elec
    ! electrodes
    type(Elec), intent(inout) :: Elecs(N_Elec)
    ! The device region unit-cell
    real(dp), intent(in) :: cell(3,3)
    ! Last orbital of each atom
    integer, intent(in) :: na_u, lasto(0:na_u)
    ! The distribution for the sparsity pattern
    type(OrbitalDistribution), intent(in) :: dit
    ! The sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! The supercell information
    integer, intent(in) :: nsc, isc_off(3,nsc)

    integer :: iEl, na
    integer :: list_a(na_u)

    ! A temporary sparsity pattern
    type(Sparsity) :: sp_tmp

    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline => null()
    character(len=50) :: g
    integer :: i, ia, ia1, ia2, no_u
    type(tRegion) :: r_tmp, r_tmp2, r_tmp3, r_Dev
    real(dp) :: p(3), contrib

    no_u = lasto(na_u)
    
    ! Instantiate regions
    if ( fdf_defined('TBT.Atoms.Buffer') ) then
       call ts_init_regions('TBT',N_Elec,Elecs,na_u,lasto)
    else
       call ts_init_regions('TS',N_Elec,Elecs,na_u,lasto)
    end if

    ! populate the buffer region 
    na = 0
    do ia = 1 , na_u
       if ( atom_type(ia) == TYP_BUFFER ) then
          na = na + 1
          list_a(na) = ia
       end if
    end do
    call region_list(r_aBuf,na,list_a, name = '[A]-buffer')
    call region_Atom2Orb(r_aBuf,na_u,lasto,r_oBuf)
    r_oBuf%name = '[O]-buffer'

    ! Create the sparsity pattern and remove the buffer atoms...
    if ( r_oBuf%n > 0 ) then
       sp_tmp = sp
       call Sp_remove_region(dit,sp_tmp,r_oBuf,sp)
       call delete(sp_tmp)
    end if

#ifdef TRANSIESTA_DEBUG
    open(file='NO_BUF_SP',unit=1400,form='formatted')
    call sp_to_file(1400,sp)
    close(1400)
#endif

    ! Create all the "alone" electrode regions
    allocate(r_aEl_alone(N_Elec),r_oEl_alone(N_Elec))
    do iEl = 1 , N_Elec

       ! Create electrode region
       ia1 = Elecs(iEl)%idx_a
       ia2 = ia1 - 1 + TotUsedAtoms(Elecs(iEl))
       call region_range(r_aEl_alone(iEl), ia1, ia2)
       r_aEl_alone(iEl)%name = '[A]-'//trim(Elecs(iEl)%name)

       ia1 = Elecs(iEl)%idx_o
       ia2 = ia1 - 1 + TotUsedOrbs(Elecs(iEl))
       call region_range(r_oEl_alone(iEl),ia1,ia2)
       r_oEl_alone(iEl)%name = '[O]-'//trim(Elecs(iEl)%name)
       
       ! Check that we have a legal region
       if ( region_overlaps(r_aEl_alone(iEl),r_aBuf) ) then
          write(*,*)'Overlapping electrode: '//trim(Elecs(iEl)%name)
          call die('Buffer region overlaps with an electrode &
               &please correct your input!')
       end if

    end do

    ! Read in device region via the new block
    if ( fdf_block('TBT.Atoms.Device',bfdf) ) then

       ! read by line and set them to be buffer atoms
       do while ( fdf_bline(bfdf,pline) ) 
          ! empty line
          if ( fdf_bnnames(pline) == 0 ) cycle
       
          g = fdf_bnames(pline,1)
          if ( leqi(g,'atom') ) then
             ! We can read in a range
             call fdf_brange(pline,r_tmp,1,na_u)
             call region_copy(r_aDev,r_tmp2)
             call region_union(r_tmp2,r_tmp,r_aDev)
             
          end if
          
       end do
       call region_delete(r_tmp,r_tmp2)

    else
       
       ! populate the device region with all but the 
       ! electrodes and buffer atoms
       na = 0
       do ia = 1 , na_u
          if ( atom_type(ia) == TYP_DEVICE ) then
             na = na + 1
             list_a(na) = ia
          end if
       end do

       call region_list(r_aDev,na,list_a)

    end if

    if ( r_aDev%n == 0 ) then
       call die('Zero atoms are in the device region...???')
    end if
    ! list_a(1:na) now contains the device region

    ! Create device region
    call region_Atom2Orb(r_aDev,na_u,lasto,r_oDev)

    ! In case the user wants "a correct DOS"
    ! in this region, we extend it
    if ( fdf_get('TBT.Atoms.Device.Connect',.false.) ) then

       ! TBTrans will truncate connections at electrode interfaces.
       call region_connect(r_oDev, dit, sp, r_tmp)
       call region_union(r_oDev,r_tmp,r_tmp2)
       call region_copy(r_tmp2,r_oDev)
       call region_delete(r_tmp,r_tmp2)
       call region_Orb2Atom(r_oDev,na_u,lasto,r_aDev)
       ! Ensure that we do not have any atoms from the electrodes
       ! This will make it behave like the "old" tbtrans
       ! in that we do not correctly capture the DOS as S_C-El /= 0 
       ! We print to the user when this occurs...
       na = r_aDev%n
       ! We remove all "electrode" implicit regions
       do iEl = 1 , N_Elec
          call region_complement(r_aDev,r_aEl_alone(iEl),r_tmp)
          call region_copy(r_tmp,r_aDev)
       end do
       call region_delete(r_tmp)

       if ( na /= r_aDev%n .and. Node == 0 ) then
          write(*,'(a)')'tbtrans: Device regions connects directly with electrodes'
          write(*,'(a)')'tbtrans: If the overlap is large this might produce spurious effects in DOS calculations'
       end if

       ! In its current state we force the entire atoms
       ! to be in the orbital connection scheme (even though
       ! some orbitals might not connect...)
       call region_Atom2Orb(r_aDev,na_u,lasto,r_oDev)

    end if

#ifdef TRANSIESTA_DEBUG
    open(file='FULL_SP',unit=1400,form='formatted')
    call sp_to_file(1400,sp)
    close(1400)
#endif

    ! Allocate the different regions
    allocate(r_aEl(N_Elec),r_oEl(N_Elec))
    allocate(r_aElpD(N_Elec),r_oElpD(N_Elec))

    do iEl = 1 , N_Elec

       ! Remove the connections that cross the boundary
       ! starting from this electrode
       call region_connect(r_oEl_alone(iEl), dit, sp, r_tmp)
       call region_union(r_oEl_alone(iEl), r_tmp, r_tmp2)

       ! Calculate the transport direction in the device cell.
       p = SPC_PROJ(cell,Elecs(iEl)%ucell(:,Elecs(iEl)%t_dir))

       ! Loop over cell vectors
       do i = 1 , 3 

          ! project the unit-cell vector onto each cell component
          contrib = VNORM(VEC_PROJ(cell(:,i),p))

          ! If the contribution in this cell direction is too
          ! small we consider it not to be important.
          ! TODO this might in certain skewed examples be a bad choice.
          if ( contrib < 1.e-5_dp ) cycle

          ! Remove connections from this electrode across the boundary...
          call Sp_remove_crossterms(dit,sp,nsc,isc_off, &
               i, &
               sp, r = r_tmp2)

       end do

       sp%data%name = 'TBT-sparsity'

       ! Check that the device region does not overlap
       if ( region_overlaps(r_aEl_alone(iEl),r_aDev) ) then
          write(*,*)'Overlapping electrode: '//trim(Elecs(iEl)%name)
          call die('Device region overlaps with an electrode &
               &please correct your input!')
       end if

    end do

    do iEl = 1 , N_Elec - 1

       ! in order to get the correct connections
       ! i.e. without connections across the device region
       ! we need to remove the connections to the other 
       ! electrodes
       ! Step 1. build a unified region of all the following
       !         electrodes!
       call region_delete(r_tmp2)
       do i = iEl + 1 , N_Elec
          call region_copy(r_tmp2,r_tmp)
          call region_union(r_oEl_alone(i),r_tmp,r_tmp2)
       end do
       call region_delete(r_tmp)

       ! First we update the sparsity pattern to remove any connections
       ! between the electrode and the other ones
       ! It will NOT remove connections between the central region and
       ! the other electrodes!
       sp_tmp = sp
       call Sp_remove_region2region(dit,sp_tmp,r_oEl_alone(iEl),r_tmp2,sp)
       call delete(sp_tmp)

    end do

#ifdef TRANSIESTA_DEBUG
    open(file='NO_ELECTRODE_CONNECTIONS_SP',unit=1400,form='formatted')
    call sp_to_file(1400,sp)
    close(1400)
#endif

    ! Create the global transiesta H(k), S(k) sparsity pattern
    call ts_Sparsity_Global(dit,sp, N_Elec, Elecs, sp_uc)

    do iEl = 1 , N_Elec

       if ( mod(iEl-1,Nodes) /= Node ) cycle

       ! Sort the regions according to the connections to
       ! the device...
       call region_copy(r_oEl_alone(iEl), r_oEl(iEl))

       ! this is a sort of the electrode
       ! TODO, try without sorting the electrode...
       call region_sort(r_oEl(iEl), dit, sp, r_oEl_alone(iEl), R_SORT_MAX_FRONT )
       r_tmp%n = 1 ! simple force of do-loop

       ! Create the region that connects out to the device
       do while ( r_tmp%n /= 0 )

          ! Create the region that connects to the last part of the 
          ! added orbitals
          call region_connect(r_oEl(iEl), dit, sp, r_tmp, except = r_oDev)

          ! r_tmp contains the connecting region (except the device region)

          ! Append the newly found region that is connecting out to the
          ! full region
          call region_union(r_oEl(iEl), r_tmp, r_tmp2)
          call region_copy(r_tmp2,r_oEl(iEl))

          ! we sort the newly attached region
          call region_sort(r_oEl(iEl), dit, sp, r_tmp, R_SORT_MAX_BACK )
       end do

       call region_delete(r_tmp,r_tmp2)

       ! This aligns the atoms in the same way the orbitals 
       ! introduce the atoms.
       call region_Orb2Atom(r_oEl(iEl), na_u, lasto , r_aEl(iEl))

       ! Create the region that connects the electrode-followed
       ! region to the central region
       call region_connect(r_oEl(iEl), dit, sp, Elecs(iEl)%o_inD)
       ! Append the newly found region that is connecting out to the
       ! full region
       call region_union(r_oEl(iEl), Elecs(iEl)%o_inD, r_oElpD(iEl))

       ! We now know how many orbitals that we are down-folding the 
       ! electrode self-energy to
       if ( Elecs(iEl)%o_inD%n == 0 ) then
          call die('The electrode down-folding region is 0 in the device &
               &region. Please expand your device region.')
       end if

       ! Create the atom equivalent regions
       call region_Orb2Atom(r_oElpD(iEl) , na_u, lasto, r_aElpD(iEl) )

       ! Clean up the elements that we do not need
       call region_delete(r_tmp2)

    end do

    do iEl = 1 , N_Elec

#ifdef MPI
       ! Bcast the region
       call region_MPI_Bcast(r_oEl(iEl),mod(iEl-1,Nodes),MPI_Comm_World)
       call region_MPI_Bcast(r_aEl(iEl),mod(iEl-1,Nodes),MPI_Comm_World)
       call region_MPI_Bcast(Elecs(iEl)%o_inD,mod(iEl-1,Nodes),MPI_Comm_World)
       call region_MPI_Bcast(r_oElpD(iEl),mod(iEl-1,Nodes),MPI_Comm_World)
       call region_MPI_Bcast(r_aElpD(iEl),mod(iEl-1,Nodes),MPI_Comm_World)
#endif

       if ( iEl > 1 ) then
       ! Check that the region does not overlap with any previous 
       ! electrode region...
       do i = 1 , iEl - 1
          if ( region_overlaps(r_aEl(iEl),r_aEl(i)) ) then
             call die('Electrode regions connect across the device region, &
                  &please increase your device region!')
          end if
       end do
       end if

       ! Set the names
       r_oEl(iEl)%name    = '[O]-'//trim(Elecs(iEl)%name)//' folding region'
       r_aEl(iEl)%name    = '[A]-'//trim(Elecs(iEl)%name)//' folding region'
       r_oElpD(iEl)%name  = '[O]-'//trim(Elecs(iEl)%name)//' folding El + D'
       r_aElpD(iEl)%name  = '[A]-'//trim(Elecs(iEl)%name)//' folding El + D'
       
       ! Check that the electrode down-folded self-energy
       ! is fully contained
       ia1 = minval(region_pivot(r_oDev,Elecs(iEl)%o_inD%r))
       if ( ia1 <= 0 ) then
          call die('A downfolded region is not existing. Programming error')
       end if

    end do

    ! We need the complement of the device region
    call region_range(r_tmp,1,no_u)
    call region_complement(r_tmp,r_oDev,r_tmp3)
    
    ! We sort the device region based on the
    ! first orbital that the electrode connects the most too.
    ! This seems like a good choice as we know 
    ! it will be on the boundary between one electrode
    ! and the device region

    ! Re-create the device region from a sorting algorithm.
    ia1 = r_oDev%n
    call region_copy(r_oDev, r_Dev)
    call region_list(r_oDev,1,Elecs(1)%o_inD%r(1:1))

    do while ( r_tmp%n /= 0 )

       ! Create the region that connects to the last part of the 
       ! added orbitals
       call region_connect(r_oDev, dit, sp, r_tmp, except = r_tmp3)

       if ( r_tmp%n == 0 .and. ia1 /= r_oDev%n ) then
          ! In case the connecting region is empty,
          ! say for capacitors we need to force the next region.
          ! In this case, we take some "random" orbital.
          do i = 1 , r_Dev%n
             if ( in_region(r_oDev,r_Dev%r(i)) ) cycle
             call region_range(r_tmp,r_Dev%r(i),r_Dev%r(i))
             exit
          end do
       end if
       ! r_tmp contains the connecting region (except all dwn-folding regions)

       ! Append the newly found region that is connecting out to the
       ! full region
       call region_union(r_oDev, r_tmp, r_tmp2)
       call region_copy(r_tmp2,r_oDev)

       ! we sort the newly attached region
       call region_sort(r_oDev, dit, sp, r_tmp, R_SORT_MAX_BACK )
       
    end do
    r_oDev%name = '[O]-device'
    call region_delete(r_Dev)


    ! Check that we have correctly re-captured the device region.
    if ( region_overlaps(r_oDev,r_tmp3) ) then
       print *,'Overlapping regions...'
       call die('Error in programming')
    end if
    if ( r_oDev%n + r_tmp3%n /= no_u .or. &
         r_oDev%n /= ia1 ) then
       print *,r_oDev%n,r_tmp3%n,no_u,ia1
       call die('Error in number of orbitals...')
    end if

    call region_Orb2Atom(r_oDev,na_u,lasto,r_aDev)
    r_aDev%name = '[A]-device'

    ! The down-folded region can "at-will" be sorted
    ! in the same manner it is seen in the device region.
    ! We enforce this as it increases the chances of consecutive 
    ! memory layout.
    do iEl = 1 , N_Elec
       
       call region_copy(Elecs(iEl)%o_inD,r_tmp2)

       ! Loop on the device region and copy
       ! region, in order
       ia1 = 0
       do i = 1 , r_oDev%n
          if ( .not. in_region(r_tmp2,r_oDev%r(i)) ) cycle
          
          ia1 = ia1 + 1
          Elecs(iEl)%o_inD%r(ia1) = r_oDev%r(i)
          
       end do

       ! Copy this information to the ElpD
       i = r_oElpD(iEl)%n
       if ( ia1 /= Elecs(iEl)%o_inD%n ) &
            call die('Error programming, ia1')
       r_oElpD(iEl)%r(i-ia1+1:i) = Elecs(iEl)%o_inD%r(1:ia1)

    end do

    ! Clean-up
    call region_delete(r_tmp,r_tmp2,r_tmp3)

    ! Do a final check that all regions are correctly setup
    ! We know that the sum of each segment has to be the 
    ! total number of orbitals in the region.
    i = r_oDev%n
    do iEl = 1 , N_Elec
       i = i + r_oEl(iEl)%n
    end do
    i = i + r_oBuf%n
    if ( i /= no_u ) then
       call die('Something went wrong when asserting the &
            &total number of orbitals. Have you requested &
            &something not applicable?')
    end if
    
  end subroutine tbt_init_regions

  subroutine tbt_init_kregions(dit,sp, na_u, lasto, nsc, isc_off)
    
    ! the orbital distribution for the sparsity pattern
    type(OrbitalDistribution), intent(in) :: dit
    ! sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! Atomic orbital configuration
    integer, intent(in) :: na_u, lasto(0:na_u)
    ! the supercell information
    integer, intent(in) :: nsc(3), isc_off(3,product(nsc))

    ! ** local variables
    integer :: no_l, no_u
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)

    call attach(sp,nrows=no_l,nrows_g=no_u, &
         n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)
    if ( no_u /= no_l ) call die('Error in fold creation')

    ! here we need to find the regions
    ! that are split by a "Gamma"-region...
    ! possibly this requirement could be
    ! relaxed by imposing a specific k-point for certain
    ! regions by the user.

    print '(a)',"TBtrans, currently k->k' does not work..."

    call die('')

  end subroutine tbt_init_kregions

  subroutine tbt_print_regions(N_Elec, Elecs)

    use parallel, only : Node
    use m_ts_electype
    
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    integer :: iEl

    if ( Node /= 0 ) return

    write(*,*) ! new-line

    ! Print out the buffer regions
    if ( r_aBuf%n > 0 ) then
       call region_print(r_aBuf, seq_max = 10 )
       !call region_print(r_oBuf, seq_max = 8 )
    end if

    ! Print out the device region
    write(*,'(a,i0)')'tbtrans: # of device region orbitals: ',r_oDev%n
    call region_print(r_aDev, seq_max = 10 )
    !call region_print(r_oDev, seq_max = 8 )

    ! Print out all the electrodes + their projection region
    do iEl = 1 , N_Elec
       write(*,*) ! new-line
       write(*,'(3a,i0)')'tbtrans: # of ',trim(Elecs(iEl)%name), &
            ' scattering orbitals: ',Elecs(iEl)%o_inD%n
       write(*,'(3a,i0)')'tbtrans: # of ',trim(Elecs(iEl)%name), &
            ' down-folding orbitals: ',r_oElpD(iEl)%n
       call region_print(r_aEl  (iEl) , seq_max = 10 )
       !call region_print(r_oEl  (iEl) , seq_max = 8 )
       call region_print(r_aElpD(iEl) , seq_max = 10 )
       !call region_print(r_oElpD(iEl) , seq_max = 8 )
       !call region_print(r_aElinD(iEl), seq_max = 10 )
       !call region_print(r_oElinD(iEl), seq_max = 8 )
    end do

  end subroutine tbt_print_regions

end module m_tbt_regions
