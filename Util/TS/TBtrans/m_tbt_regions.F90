! This module contains
! the different regions used in tbtrans

! Coded by Nick Papior Andersen
! 2014

module m_tbt_regions

  use m_region
  use class_OrbitalDistribution
  use class_Sparsity
  ! To re-use as much from transiesta
  use m_ts_method, only : r_aBuf, r_oBuf
  use m_ts_method, only : r_aDev => r_aC, r_oDev => r_oC

  implicit none

  private
  save

  ! The full tbtrans sparsity region in the unit-cell equivalent
  ! region, this is used to accomodate the Hamiltonian and overlap
  ! matrices. In principle this could be made obsolete.
  type(Sparsity), public :: sp_uc ! TBT-GLOBAL (UC)

  ! The UC sparsity pattern in the device region
  type(Sparsity), public :: sp_dev

  ! the different regions that becomes the electrodes
  type(tRgn), allocatable, target, public :: r_aEl_alone(:), r_oEl_alone(:)

  ! the different regions that connects to the equivalent
  ! electrodes.
  ! I.e. it is the regions that goes from the electrode
  ! and down to the central region (without any central
  ! region overlap)
  type(tRgn), allocatable, target, public :: r_aEl(:)   , r_oEl(:)
  type(tRgn), allocatable, target, public :: r_aElpD(:) , r_oElpD(:)

  ! the device region (the calculated GF region)
  public :: r_aDev, r_oDev

  ! The buffer region, just for completeness
  public :: r_aBuf, r_oBuf

  public :: tbt_init_regions, tbt_read_regions
  public :: tbt_region_options
  public :: tbt_print_regions

contains

  subroutine tbt_init_regions(N_Elec, Elecs, cell, na_u, lasto, dit, sp, &
       nsc, isc_off)

    use fdf
    use fdf_extra
    use parallel, only : Node, Nodes
#ifdef MPI
    use mpi_siesta, only : MPI_Comm_World, MPI_Barrier
#endif

    use m_pivot

    use m_io_s, only : file_exist
    use geom_helper, only : iaorb
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
    character(len=50) :: g, csort
    integer :: i, ia, ia1, ia2, no_u, iu, n
    type(tRgn) :: r_tmp, r_tmp2, r_tmp3, r_Els, r_Dev, priority
    logical :: sort_orb
    real(dp) :: p(3), contrib

    no_u = lasto(na_u)
    
    ! Instantiate regions
    if ( fdf_defined('TBT.Atoms.Buffer') ) then
       call ts_init_regions('TBT',N_Elec,Elecs,na_u,lasto)
    else
       call ts_init_regions('TS',N_Elec,Elecs,na_u,lasto)
    end if

    ! Re=name the read-in information
    r_aBuf%name = '[A]-buffer'
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
       call rgn_range(r_aEl_alone(iEl), ia1, ia2)
       r_aEl_alone(iEl)%name = '[A]-'//trim(Elecs(iEl)%name)

       ia1 = Elecs(iEl)%idx_o
       ia2 = ia1 - 1 + TotUsedOrbs(Elecs(iEl))
       call rgn_range(r_oEl_alone(iEl),ia1,ia2)
       r_oEl_alone(iEl)%name = '[O]-'//trim(Elecs(iEl)%name)
       
       ! Check that we have a legal region
       if ( rgn_overlaps(r_aEl_alone(iEl),r_aBuf) ) then
          write(*,*)'Overlapping electrode: '//trim(Elecs(iEl)%name)
          call die('Buffer region overlaps with an electrode &
               &please correct your input!')
       end if

    end do

    ! Delete to be ready to populate the device
    call rgn_delete(r_aDev)

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
             if ( r_tmp%n == 0 ) &
                  call die('Could not read in any atoms &
                  &in line of TBT.Atoms.Device')
             call rgn_union(r_aDev,r_tmp,r_aDev)
             
          end if
          
       end do
       call rgn_delete(r_tmp)

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

       call rgn_list(r_aDev,na,list_a)

    end if

    if ( r_aDev%n == 0 ) then
       call die('Zero atoms are in the device region...???')
    end if

    ! Create device region
    call rgn_Atom2Orb(r_aDev,na_u,lasto,r_oDev)

    if ( Node == 0 ) then
       write(*,'(/,a)')'tbtrans: Analyzing electrode sparsity pattern to create optimal tri-diagonal blocks...'
    end if

    ! In case the user wants "a correct DOS"
    ! in this region, we extend it
    if ( fdf_get('TBT.Atoms.Device.Connect',.false.) ) then

       ! TBTrans will truncate connections at electrode interfaces.
       call rgn_sp_connect(r_oDev, dit, sp, r_tmp)
       if ( r_tmp%n == 0 ) &
            call die('No orbitals connect to the specified device &
            &region. This is not allowed.')
       call rgn_append(r_oDev,r_tmp,r_oDev)
       call rgn_delete(r_tmp)
       call rgn_Orb2Atom(r_oDev,na_u,lasto,r_aDev)
       ! Ensure that we do not have any atoms from the electrodes
       ! This will make it behave like the "old" tbtrans
       ! in that we do not correctly capture the DOS as S_C-El /= 0 
       ! We print to the user when this occurs...
       na = r_aDev%n
       ! We remove all "electrode" implicit regions
       do iEl = 1 , N_Elec
          call rgn_complement(r_aEl_alone(iEl),r_aDev,r_aDev)
       end do

       if ( na /= r_aDev%n .and. Node == 0 ) then
          write(*,'(a)')'tbtrans: Device regions connects directly with electrodes'
          write(*,'(a)')'tbtrans: If the overlap is large this might produce spurious effects in DOS calculations'
       end if

       ! In its current state we force the entire atoms
       ! to be in the orbital connection scheme (even though
       ! some orbitals might not connect...)
       call rgn_Atom2Orb(r_aDev,na_u,lasto,r_oDev)

    end if

    ! Makes searching a little faster
    call rgn_sort(r_aDev)
    call rgn_sort(r_oDev)

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
       call rgn_sp_connect(r_oEl_alone(iEl), dit, sp, r_tmp)
       call rgn_union(r_oEl_alone(iEl), r_tmp, r_tmp2)

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
       if ( rgn_overlaps(r_aEl_alone(iEl),r_aDev) ) then
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
       call rgn_delete(r_tmp)
       do i = iEl + 1 , N_Elec
          call rgn_append(r_tmp,r_oEl_alone(i),r_tmp)
       end do

       ! First we update the sparsity pattern to remove any connections
       ! between the electrode and the other ones
       ! It will NOT remove connections between the central region and
       ! the other electrodes!
       sp_tmp = sp
       call Sp_remove_region2region(dit,sp_tmp,r_oEl_alone(iEl),r_tmp,sp)
       call delete(sp_tmp)

    end do

#ifdef TRANSIESTA_DEBUG
    open(file='NO_ELECTRODE_CONNECTIONS_SP',unit=1400,form='formatted')
    call sp_to_file(1400,sp)
    close(1400)
#endif

    ! Create the global transiesta H(k), S(k) sparsity pattern
    call ts_Sparsity_Global(dit,sp, N_Elec, Elecs, sp_uc)

    ! At this point we have reduced the sparsity
    ! pattern to it's correct size.
    ! Try and read in the regions
    call tbt_read_regions(N_Elec,Elecs,na_u,lasto,no_u)
    if ( r_oEl(1)%n > 0 ) then
       ! We say that if the electrode orbitals are allocated
       ! we have read in the file-information
       call rgn_delete(r_tmp,r_tmp2,r_tmp3)
       return
    end if

    ! Create the electrode down-folding regions.
    ! Note that sorting according to the more advanced methods
    ! is not directly applicable as the methods involve non-stringent
    ! ending elements.

    do iEl = 1 , N_Elec

       if ( mod(iEl-1,Nodes) /= Node ) cycle

       ! Sort the regions according to the connections to
       ! the device...
       call rgn_copy(r_oEl_alone(iEl), r_oEl(iEl))

       ! this is a sort of the electrode
       ! TODO, try without sorting the electrode...
       call rgn_sp_sort(r_oEl(iEl), dit, sp, r_oEl_alone(iEl), R_SORT_MAX_FRONT )

       ! Create the region that connects out to the device
       do

          ! Create the region that connects to the last part of the 
          ! added orbitals
          call rgn_sp_connect(r_oEl(iEl), dit, sp, r_tmp, except = r_oDev)

          ! If no additional orbitals are found, exit
          if ( r_tmp%n == 0 ) exit

          ! r_tmp contains the connecting region (except the device region)

          ! Append the newly found region that is connecting out to the
          ! full region
          call rgn_append(r_oEl(iEl), r_tmp, r_oEl(iEl))

          ! we sort the newly attached region
          call rgn_sp_sort(r_oEl(iEl), dit, sp, r_tmp, R_SORT_MAX_BACK )

       end do

       call rgn_delete(r_tmp)

       ! This aligns the atoms in the same way the orbitals 
       ! introduce the atoms.
       call rgn_Orb2Atom(r_oEl(iEl), na_u, lasto , r_aEl(iEl))

       ! Create the region that connects the electrode-followed
       ! region to the central region
       call rgn_sp_connect(r_oEl(iEl), dit, sp, Elecs(iEl)%o_inD)
       ! Append the newly found region that is connecting out to the full region
       call rgn_append(r_oEl(iEl), Elecs(iEl)%o_inD, r_oElpD(iEl))

       ! We now know how many orbitals that we are down-folding the 
       ! electrode self-energy to
       if ( Elecs(iEl)%o_inD%n == 0 ) then
          call die('The electrode down-folding region is 0 in the device &
               &region. Please expand your device region.')
       end if

       ! Create the atom equivalent regions
       call rgn_Orb2Atom(r_oElpD(iEl) , na_u, lasto, r_aElpD(iEl) )

    end do

    do iEl = 1 , N_Elec

#ifdef MPI
       ! Bcast the region
       call rgn_MPI_Bcast(r_oEl(iEl),mod(iEl-1,Nodes),MPI_Comm_World)
       call rgn_MPI_Bcast(r_aEl(iEl),mod(iEl-1,Nodes),MPI_Comm_World)
       call rgn_MPI_Bcast(Elecs(iEl)%o_inD,mod(iEl-1,Nodes),MPI_Comm_World)
       call rgn_MPI_Bcast(r_oElpD(iEl),mod(iEl-1,Nodes),MPI_Comm_World)
       call rgn_MPI_Bcast(r_aElpD(iEl),mod(iEl-1,Nodes),MPI_Comm_World)
#endif

       if ( iEl > 1 ) then
       ! Check that the region does not overlap with any previous 
       ! electrode region...
       do i = 1 , iEl - 1
          if ( rgn_overlaps(r_aEl(iEl),r_aEl(i)) ) then
             call die('Electrode regions connect across the device region, &
                  &please increase your device region!')
          end if
       end do
       end if

       ! Set the names
       r_oEl(iEl)%name   = '[O]-'//trim(Elecs(iEl)%name)//' folding region'
       r_aEl(iEl)%name   = '[A]-'//trim(Elecs(iEl)%name)//' folding region'
       r_oElpD(iEl)%name = '[O]-'//trim(Elecs(iEl)%name)//' folding El + D'
       r_aElpD(iEl)%name = '[A]-'//trim(Elecs(iEl)%name)//' folding El + D'
       
       ! Check that the electrode down-folded self-energy
       ! is fully contained
       ia1 = minval(rgn_pivot(r_oDev,Elecs(iEl)%o_inD%r))
       if ( ia1 <= 0 ) then
          call die('A downfolded region is not existing. Programming error')
       end if

    end do

    if ( Node == 0 ) then
       write(*,'(a)')'tbtrans: Analyzing device sparsity pattern to &
            &create optimal tri-diagonal blocks...'
    end if
    
    ! We sort the device region based on the
    ! first orbital that the electrode connects the most too.
    ! This seems like a good choice as we know 
    ! it will be on the boundary between one electrode
    ! and the device region

    ! Re-create the device region from a sorting algorithm.
    ! However, in cases where the user is sure that
    ! the device region is correctly sorted 
    ! then we need not re-sort it.
    ! This seems like the best choice when looking
    ! at TB models where number of connections is the same

    ! Prepare total region with all electrodes
    call rgn_copy(Elecs(1)%o_inD,r_tmp)
    do i = 2 , N_Elec
       ! We have to use union as they may overlap
       call rgn_union(r_Els,Elecs(i)%o_inD,r_Els)
    end do

    csort = 'atom+'//trim(Elecs(1)%name)
    csort = fdf_get('TS.BTD.Pivot',csort)
    csort = fdf_get('TBT.BTD.Pivot',csort)
    csort = fdf_get('TBT.BTD.Pivot.Device',csort)
    g = csort

    ! we default to do atomic sorting, faster and very consistent
    sort_orb = index(g,'orb') > 0

    if ( sort_orb ) then

       ! the specification is orb + something
       i = index(g,'orb')
       g(i:i+3) = ' ' ! also remove +, we ADJUSTL below

       call rgn_copy(r_oDev,r_tmp)

       sp_tmp = sp_uc

    else

       ! in case the user supplies 'atom+something'
       i = index(g,'atom')
       if ( i > 0 ) then
          g(i:i+4) = ' ' ! also remove +
       else
          csort = 'atom+'//trim(csort)
       end if

       call rgn_copy(r_Els,r_tmp)
       call rgn_Orb2Atom(r_tmp,na_u,lasto,r_Els)
       call rgn_copy(r_aDev,r_tmp)

       call SpOrb_to_SpAtom(dit,sp_uc,na_u,lasto,sp_tmp)

    end if

    ! Print out what we found
    if ( Node == 0 ) then
       write(*,'(a)')'tbtrans: BTD pivoting scheme in device: '//trim(csort)
    end if

    ! The device region
    call rgn_copy(r_tmp,r_Dev)
    ! Sort the electrode region to make it faster
    call rgn_sort(r_Els)

    ! Now r_tmp contains the copy of the "sub" section
    ! of the wanted region, i.e. the device region
    n = nrows_g(sp_tmp)

    ! We need the complement of the device region
    call rgn_range(r_tmp2,1,n)
    call rgn_complement(r_Dev,r_tmp2,r_tmp3) ! r_3 == complement of a/oDev

    ! Create priority list
    ! Currently the priority is hard to create using 
    ! a down-folded region.
    ! We cannot assert the semi-infinite sorting, or at least much less
    call rgn_init(priority,n)
    call crt_El_priority(N_Elec,Elecs,priority,na_u,lasto, is_orb = sort_orb )

    ! Left adjust the string
    g = ADJUSTL(g)

    if ( leqi(g,'none') ) then

       ! do nothing, it should not be sorted.
       
    else if ( leqi(g,'CM') ) then

       call sp_pvt(n,sp_tmp,r_tmp2, PVT_CUTHILL_MCKEE, r_Dev, start = r_Els)

    else if ( leqi(g,'CM+priority') ) then

       call sp_pvt(n,sp_tmp,r_tmp2, PVT_CUTHILL_MCKEE, r_Dev, start = r_Els, &
            priority = priority%r )

    else if ( leqi(g,'rev-CM') ) then

       call sp_pvt(n,sp_tmp,r_tmp2, PVT_REV_CUTHILL_MCKEE, r_Dev, start = r_Els)

    else if ( leqi(g,'rev-CM+priority') ) then

       call sp_pvt(n,sp_tmp,r_tmp2, PVT_REV_CUTHILL_MCKEE, r_Dev, start = r_Els , &
            priority = priority%r )

    else if ( leqi(g,'GPS') ) then

       call sp_pvt(n,sp_tmp,r_tmp2, PVT_GPS, r_Dev)

    else if ( leqi(g,'GPS+priority') ) then

       call sp_pvt(n,sp_tmp,r_tmp2, PVT_GPS, r_Dev, priority = priority%r )

    else if ( leqi(g,'rev-GPS') ) then

       call sp_pvt(n,sp_tmp,r_tmp2, PVT_REV_GPS, r_Dev)
       
    else if ( leqi(g,'rev-GPS+priority') ) then
       
       call sp_pvt(n,sp_tmp,r_tmp2, PVT_REV_GPS, r_Dev, priority = priority%r )

    else if ( leqi(g,'PCG') ) then

       call sp_pvt(n,sp_tmp,r_tmp2, PVT_PCG, r_Dev)

    else if ( leqi(g,'PCG+priority') ) then

       call sp_pvt(n,sp_tmp,r_tmp2, PVT_PCG, r_Dev, priority = priority%r )

    else if ( leqi(g,'rev-PCG') ) then

       call sp_pvt(n,sp_tmp,r_tmp2, PVT_REV_PCG, r_Dev)
       
    else if ( leqi(g,'rev-PCG+priority') ) then
       
       call sp_pvt(n,sp_tmp,r_tmp2, PVT_REV_PCG, r_Dev, priority = priority%r )
       
    else if ( leqi(g,'GGPS') ) then

       call sp_pvt(n,sp_tmp,r_tmp2, PVT_GGPS, r_Dev)

    else if ( leqi(g,'GGPS+priority') ) then

       call sp_pvt(n,sp_tmp,r_tmp2, PVT_GGPS, r_Dev, priority = priority%r )

    else if ( leqi(g,'rev-GGPS') ) then

       call sp_pvt(n,sp_tmp,r_tmp2, PVT_REV_GGPS, r_Dev)

    else if ( leqi(g,'rev-GGPS+priority') ) then
       
       call sp_pvt(n,sp_tmp,r_tmp2, PVT_REV_GGPS, r_Dev, priority = priority%r )
       
    else
       

       ! We select all the atoms from the specified electrode
       ! that connect into the device
       ! This forces the initial search-pattern
       ! to use the "electrode-surface" as an initial SP guess
       ! Note that the connecting electrode orbitals are already
       ! sorted wrt. back-connectivity

       ! Figure out which electrode has been given
       iEl = 0
       do i = 1 , N_Elec
          if ( leqi(Elecs(i)%name,g) ) iEl = i
       end do
       if ( iEl == 0 ) then
          print *,csort
          call die('Could find the electrode in &
               &TBT.BTD.Pivot.Device in the list of electrodes, &
               &please correct sorting method.')
       end if

       if ( sort_orb ) then
          call rgn_copy(Elecs(iEl)%o_inD,r_tmp2)
       else
          call rgn_Orb2Atom(Elecs(iEl)%o_inD,na_u,lasto,r_tmp2)
       end if
       if ( r_tmp2%n == 0 ) then
          call die('Could not determine sorting electrode. &
               &Please provide a valid electrode name or none.')
       end if
       
       do

          ! Create the region that connects to the last part of the 
          ! added orbitals
          call rgn_sp_connect(r_tmp2, dit, sp_tmp, r_tmp, except = r_tmp3)
          
          if ( r_tmp%n == 0 .and. r_Dev%n /= r_tmp2%n ) then

             ! In case the connecting region is empty,
             ! say for capacitors we need to force the next region.
             ! In this case, we take some "random" orbital.
             do i = 1 , r_Dev%n
                if ( in_rgn(r_tmp2,r_Dev%r(i)) ) cycle
                call rgn_range(r_tmp,r_Dev%r(i),r_Dev%r(i))
                exit
             end do

             ! get the atomic index
             if ( sort_orb ) then
                i = iaorb(r_tmp%r(1),lasto)
             else
                i = r_tmp%r(1)
             end if
                
             if ( Node == 0 ) then
                write(*,'(/,a,i0,a)')'WARNING TBT: &
                     &Random atom ',i,' has been added &
                     &due to non-completeness of the connectivity &
                     &graph.'
                write(*,'(a)')'WARNING TBT: Expect sub-optimal &
                     &BTD format.'
             end if

             if ( sort_orb ) then
                call rgn_range(r_tmp,lasto(i-1)+1,lasto(i))
             end if

             ! Assert that none of the orbitals exist in the region
             if ( rgn_overlaps(r_tmp2,r_tmp) ) then
                call die('This is extremely difficult. &
                     &Please do not sort the BTD format according to &
                     &an electrode as it cannot figure out what to do.')
             end if
             
             call rgn_append(r_tmp2, r_tmp, r_tmp2)
             cycle
             
          end if
          ! r_tmp contains the connecting region (except all dwn-folding regions)

          ! If no additional orbitals are found, exit
          if ( r_tmp%n == 0 ) exit

          ! Append the newly found region that is connecting out to the
          ! full region
          call rgn_append(r_tmp2, r_tmp, r_tmp2)
          
          ! we sort the newly attached region
          call rgn_sp_sort(r_tmp2, dit, sp_tmp, r_tmp, R_SORT_MAX_BACK )
          
       end do

    end if

    ! Check that there is no overlap with the other regions
    if ( rgn_overlaps(r_tmp3,r_tmp2) ) then
       call rgn_print(r_Dev)
       call rgn_print(r_tmp3)
       call rgn_print(r_tmp2)
       print *,'Overlapping device, down-folding and/or buffer region...'
       call die('Error in programming')
    end if

    ! Check that we have correctly re-captured the device region.
    if ( .not. sort_orb ) then
       if ( r_tmp2%n + r_tmp3%n /= na_u .or. &
            r_tmp2%n /= r_Dev%n ) then
          print *,r_tmp2%n,r_tmp3%n,na_u,r_Dev%n
          call die('Error in number of atoms...')
       end if
    end if

    ! Create complement of device region
    call rgn_range(r_tmp,1,no_u)
    call rgn_complement(r_oDev,r_tmp,r_tmp3) ! r_3 == complement of oDev

    ! Copy over the sorted region
    if ( sort_orb ) then
       call rgn_copy(r_tmp2,r_oDev)
    else
       call rgn_Atom2Orb(r_tmp2,na_u,lasto,r_oDev)
    end if
    r_oDev%name = '[O]-device'

    ! Check that r_oDev contains all but the folded
    ! region and the buffer
    if ( r_oDev%n + r_tmp3%n /= no_u ) then
       print *,r_oDev%n,r_tmp3%n,no_u
       call die('Error in number of orbitals...')
    end if

    call rgn_Orb2Atom(r_oDev,na_u,lasto,r_aDev)
    r_aDev%name = '[A]-device'

    ! The down-folded region can "at-will" be sorted
    ! in the same manner it is seen in the device region.
    ! We enforce this as it increases the chances of consecutive 
    ! memory layout.
    do iEl = 1 , N_Elec
       
       call rgn_copy(Elecs(iEl)%o_inD,r_tmp2)
       call rgn_sort(r_tmp2)

       ! Loop on the device region and copy
       ! region, in order
       ia1 = 0
       do i = 1 , r_oDev%n
          if ( .not. in_rgn(r_tmp2,r_oDev%r(i)) ) cycle
          
          ia1 = ia1 + 1
          Elecs(iEl)%o_inD%r(ia1) = r_oDev%r(i)
          
       end do

       ! create the pivoting table
       call rgn_copy(Elecs(iEl)%o_inD,Elecs(iEl)%inDpvt)
       Elecs(iEl)%inDpvt%r(:) = rgn_pivot(r_oDev,Elecs(iEl)%o_inD%r(:))

       ! Copy this information to the ElpD
       i = r_oElpD(iEl)%n
       if ( ia1 /= Elecs(iEl)%o_inD%n ) &
            call die('Error programming, ia1')
       r_oElpD(iEl)%r(i-ia1+1:i) = Elecs(iEl)%o_inD%r(1:ia1)

    end do

    ! Do a final check that all regions are correctly setup
    ! We know that the sum of each segment has to be the 
    ! total number of orbitals in the region.
    i = r_oDev%n
    do iEl = 1 , N_Elec
       i = i + r_oEl(iEl)%n
    end do
    i = i + r_oBuf%n
    if ( i /= no_u ) then
       if ( Node == 0 ) then
          write(*,'(a,i0)')'Buffer orbitals: ',r_oBuf%n
          write(*,'(a,i0)')'Device orbitals: ',r_oDev%n
          do iEl = 1 , N_Elec
             write(*,'(a,i0)')trim(Elecs(iEl)%name)//' orbitals: ',r_oEl(iEl)%n
          end do
          if ( i > no_u ) then
             ! find the overlapping orbitals
             call rgn_union(r_oBuf,r_oDev,r_tmp)
             do iEl = 1 , N_Elec
                call rgn_union(r_tmp,r_oEl(iEl),r_tmp)
             end do
             r_tmp%name = 'Double counted orbitals'
             call rgn_print(r_tmp)
          else
             ! find the missing orbitals
             call rgn_range(r_tmp,1,no_u)
             if ( r_oBuf%n > 0 ) then
                call rgn_complement(r_oBuf,r_tmp,r_tmp)
             end if
             call rgn_complement(r_oDev,r_tmp,r_tmp)
             do iEl = 1 , N_Elec
                call rgn_complement(r_oEl(iEl),r_tmp,r_tmp)
             end do
             r_tmp%name = 'Missing orbitals'
             call rgn_print(r_tmp)
          end if
          write(*,'(a,2(tr1,i0))')'Total number of orbitals vs. counted:',no_u,i

          write(*,'(/,a)')'Missing/Excess orbitals can happen if your device &
               &region is ill-formatted.'
          write(*,'(a)')'Suppose you create a device region which disconnects &
               &certain non-device region orbitals from the electrode regions.'
          write(*,'(a)')'Then this will occur, please ensure that you have &
               &defined your device region such that the above does not occur.'
       end if
#ifdef MPI
       call MPI_Barrier(MPI_Comm_World,i)
#endif
       call die('Something went wrong when asserting the &
            &total number of orbitals. Have you requested &
            &something not applicable?')
    end if

    ! Clean-up
    call rgn_delete(r_tmp,r_tmp2,r_tmp3,r_Els,priority)
    call delete(sp_tmp)

    if ( Node == 0 ) then
       write(*,'(a)')'tbtrans: Done analyzing sparsity pattern...'

       ! Write the information 
       g = fdf_get('TBT.Region.File','DOESNOTEXIST')
       if ( g /= 'DOESNOTEXIST' ) then

       ! Create the file
       call io_assign(iu)
       open(file=trim(g),unit=iu,form='formatted')

       ! write information about the system size
       write(iu,'(4(i0,tr1))') na_u, no_u, N_Elec

       ! if buffer atoms exists, we write them
       write(iu,'(2(i0,tr1))') r_aBuf%n,r_oBuf%n
       if ( r_aBuf%n > 0 ) then
          ! write alone electrode, atoms and orbitals
          write(iu,'(12(i0,tr1))') (r_aBuf%r(i),i=1,r_aBuf%n)
          write(iu,'(12(i0,tr1))') (r_oBuf%r(i),i=1,r_oBuf%n)
       end if
             
       ! Write out each electrode individually
       do iEl = 1 , N_Elec

          ! Write name
          write(iu,'(a)') trim(Elecs(iEl)%name)
          ! write number of atoms and orbitals
          i = product(Elecs(iEl)%rep)
          write(iu,'(2(i0,tr1))') Elecs(iEl)%na_used*i,Elecs(iEl)%no_used*i
          ! write alone electrode, atoms and orbitals
          write(iu,'(12(i0,tr1))') (r_aEl_alone(iEl)%r(i),i=1,r_aEl_alone(iEl)%n)
          write(iu,'(12(i0,tr1))') (r_oEl_alone(iEl)%r(i),i=1,r_oEl_alone(iEl)%n)
          write(iu,'(2(i0,tr1))') r_oEl(iEl)%n
          write(iu,'(12(i0,tr1))') (r_oEl(iEl)%r(i),i=1,r_oEl(iEl)%n)
          write(iu,'(2(i0,tr1))') Elecs(iEl)%o_inD%n
          write(iu,'(12(i0,tr1))') (Elecs(iEl)%o_inD%r(i),i=1,Elecs(iEl)%o_inD%n)
       end do
             
       write(iu,'(2(i0,tr1))') r_aDev%n,r_oDev%n
       write(iu,'(12(i0,tr1))') (r_aDev%r(i),i=1,r_aDev%n)
       write(iu,'(12(i0,tr1))') (r_oDev%r(i),i=1,r_oDev%n)
       
       call io_close(iu)

       end if

    end if

  end subroutine tbt_init_regions

  subroutine tbt_read_regions(N_Elec, Elecs, na_u, lasto, no_u)
    use parallel, only : Node
    use fdf
#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD
#endif

    use m_ts_electype
    use m_io_s, only : file_exist
    use m_region

    integer, intent(in) :: N_Elec
    type(Elec), intent(inout) :: Elecs(N_Elec)
    integer, intent(in) :: na_u, lasto(0:na_u), no_u

    ! Local variables
    character(len=250) :: fname, g
    integer :: iu, fa_u, fo_u, f_Elec
    integer :: iEl, i, itmp(2)
    type(tRgn) :: r

    ! We try and see if a user file exists and then create all
    ! regions.
    fname = fdf_get('TBT.Region.File','DOESNOTEXIST')
    if ( .not. file_exist(fname, Bcast = .true.) ) return

    if ( Node == 0 ) then

    write(*,'(a)')'tbtrans: Re-reading previously created regions...'

    ! Read in the quantities
    call io_assign(iu)
    open(file=trim(fname),unit=iu,form='formatted')

    ! write information about the system size
    read(iu,*) fa_u, fo_u, f_Elec
    if ( fa_u /= na_u ) &
         call mdie('Error in number of atoms in full structure.')
    if ( fo_u /= no_u ) &
         call mdie('Error in number of orbitals in full structure.')
    if ( f_Elec /= N_Elec ) &
         call mdie('Error in number of electrodes.')

    read(iu,*) itmp
    if ( itmp(1) /= r_aBuf%n ) &
         call mdie('Error in number of buffer atoms.')
    if ( itmp(2) /= r_oBuf%n ) &
         call mdie('Error in number of buffer orbitals.')
    ! if buffer atoms exists, we write them
    if ( itmp(1) > 0 ) then
       call rgn_range(r,1,itmp(1))
       ! write alone electrode, atoms and orbitals
       read(iu,*) (r%r(i),i=1,r%n)
       if ( any(r%r /= r_aBuf%r) ) &
            call mdie('Buffer atoms not equivalent.')
       call rgn_range(r,1,itmp(2))
       read(iu,*) (r%r(i),i=1,r%n)
       if ( any(r%r /= r_oBuf%r) ) &
            call mdie('Buffer orbitals not equivalent.')
    end if

    ! Read each electrode individually
    do iEl = 1 , N_Elec

       ! read name
       read(iu,*) g
       if ( trim(g) /= trim(Elecs(iEl)%name) ) &
            call mdie('Electrode name not the same.', &
            'Read: '//trim(g)//' vs. fdf: '//trim(Elecs(iEl)%name))
       
       ! read number of atoms and orbitals
       read(iu,*) itmp
       i = product(Elecs(iEl)%rep)
       if ( itmp(1) /= Elecs(iEl)%na_used * i ) &
            call mdie('Electrode used atoms are not equivalent.', &
            'Electrode: '//trim(g))
       if ( itmp(2) /= Elecs(iEl)%no_used * i ) &
            call mdie('Electrode used orbitals are not equivalent.', &
            'Electrode: '//trim(g))

       ! read alone electrode, atoms and orbitals
       call rgn_range(r,1,r_aEl_alone(iEl)%n)
       read(iu,*) (r%r(i),i=1,r%n)
       if ( any(r%r /= r_aEl_alone(iEl)%r) ) &
            call mdie('Electrode atom positions are not equivalent.', &
            'Electrode: '//trim(g))
       call rgn_range(r,1,r_oEl_alone(iEl)%n)
       read(iu,*) (r%r(i),i=1,r%n)
       if ( any(r%r /= r_oEl_alone(iEl)%r) ) &
            call mdie('Electrode orbital positions are not equivalent.', &
            'Electrode: '//trim(g))

       read(iu,*) itmp(1)
       call rgn_range(r_oEl(iEl),1,itmp(1))
       read(iu,*) (r_oEl(iEl)%r(i),i=1,r_oEl(iEl)%n)

       read(iu,*) itmp(1)
       call rgn_range(Elecs(iEl)%o_inD,1,itmp(1))
       read(iu,*) (Elecs(iEl)%o_inD%r(i),i=1,itmp(1))

    end do

    read(iu,*) itmp
    ! Check that the device regions fit
    if ( itmp(1) /= r_aDev%n ) &
         call mdie('Device region atoms are not the same.', &
         'Please delete file: '//trim(fname)//' or revert the device region block')
    read(iu,*) (r_aDev%r(i),i=1,r_aDev%n)
    call rgn_range(r_oDev,1,itmp(2))
    read(iu,*) (r_oDev%r(i),i=1,r_oDev%n)

    call io_close(iu)

    call rgn_delete(r)

    end if

#ifdef MPI
    do iEl = 1 , N_Elec
       call rgn_MPI_Bcast(r_oEl(iEl),0,MPI_Comm_World)
       call rgn_MPI_Bcast(Elecs(iEl)%o_inD,0,MPI_Comm_World)
    end do
    call rgn_MPI_Bcast(r_aDev,0,MPI_Comm_World)
    call rgn_MPI_Bcast(r_oDev,0,MPI_Comm_World)
#endif
    r_aDev%name = '[A]-device'
    r_oDev%name = '[O]-device'

    do iEl = 1 , N_Elec

       call rgn_Orb2Atom(r_oEl(iEl), na_u, lasto , r_aEl(iEl))
       r_oEl(iEl)%name    = '[O]-'//trim(Elecs(iEl)%name)//' folding region'
       r_aEl(iEl)%name    = '[A]-'//trim(Elecs(iEl)%name)//' folding region'
       call rgn_append(r_oEl(iEl), Elecs(iEl)%o_inD, r_oElpD(iEl))
       r_oElpD(iEl)%name  = '[O]-'//trim(Elecs(iEl)%name)//' folding El + D'
       call rgn_Orb2Atom(r_oElpD(iEl) , na_u, lasto, r_aElpD(iEl) )
       r_aElpD(iEl)%name  = '[A]-'//trim(Elecs(iEl)%name)//' folding El + D'

       call rgn_copy(Elecs(iEl)%o_inD,Elecs(iEl)%inDpvt)
       Elecs(iEl)%inDpvt%r(:) = rgn_pivot(r_oDev,Elecs(iEl)%o_inD%r(:))

    end do

  contains

    subroutine mdie(msg,amsg)
      character(len=*), intent(in) :: msg
      character(len=*), intent(in), optional :: amsg

      write(*,'(a)') 'TBtrans will die due to erroneous reading of &
           &the '//trim(fname)//' file.'
      write(*,'(a)') 'This will happen if you change the device region &
           &or the electrode positions.'
      write(*,'(a)') 'The easiest thing is to delete the file: '//trim(fname)
      if ( present(amsg) ) then
         write(*,'(a)') amsg
      end if
      call die(msg)
    end subroutine mdie

  end subroutine tbt_read_regions

  subroutine tbt_region_options( save_DATA )
    use dictionary
#ifdef MPI
    use mpi_siesta, only : MPI_Comm_Self
#endif
    use m_sparsity_handling, only : Sp_retain_rgn

    type(dict), intent(in) :: save_DATA
    type(OrbitalDistribution) :: fdit

    integer :: no_u

    ! Make sure to initialize the device region
    ! sparsity pattern
    call delete(sp_dev)
#ifdef NCDF_4
    if ( 'orb-current' .in. save_DATA ) then

       call attach(sp_uc,nrows_g=no_u)
#ifdef MPI
       call newDistribution(no_u,MPI_Comm_Self,fdit,name='TBT-fake dist')
#else
       call newDistribution(no_u,-1           ,fdit,name='TBT-fake dist')
#endif

       call Sp_retain_rgn(fdit,sp_uc,r_oDev,sp_dev)
       call delete(fdit)

    end if
#endif

  end subroutine tbt_region_options

  subroutine tbt_print_regions(N_Elec, Elecs)

    use parallel, only : Node
    use m_verbosity, only : verbosity
    use m_ts_electype
    
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    integer :: i, verb
    type(tRgn) :: r

    if ( Node /= 0 ) return

    if ( verbosity < 3 ) return

    ! Print out the buffer regions
    if ( r_aBuf%n > 0 ) then
       call rgn_print(r_aBuf, seq_max = 12 )
       if ( verbosity > 7 ) &
            call rgn_print(r_oBuf, seq_max = 10 )
    end if

    ! Print out the device region
    write(*,'(a,i0)')'tbtrans: # of device region orbitals: ',r_oDev%n
    if ( verbosity > 4 ) &
         call rgn_print(r_aDev, seq_max = 12 )
    if ( verbosity > 7 ) &
         call rgn_print(r_oDev, seq_max = 10 )

    ! Print out all the electrodes + their projection region
    do i = 1 , N_Elec
       write(*,*) ! new-line
       write(*,'(3a,i0)')'tbtrans: # of ',trim(Elecs(i)%name), &
            ' scattering orbitals: ',Elecs(i)%o_inD%n
       write(*,'(3a,i0)')'tbtrans: # of ',trim(Elecs(i)%name), &
            ' down-folding orbitals: ',r_oElpD(i)%n
       if ( verbosity > 4 ) &
            call rgn_print(r_aEl  (i) , seq_max = 12 )
       if ( verbosity > 7 ) &
            call rgn_print(r_oEl  (i) , seq_max = 10 )
       if ( verbosity > 3 ) then
          call rgn_intersection(r_aElpD(i),r_aDev,r)
          r%name = '[A]-'//trim(Elecs(i)%name)//' folding in D'
          call rgn_print(r, seq_max = 12 )
       end if
       if ( verbosity > 7 ) then
          call rgn_intersection(r_oElpD(i),r_oDev,r)
          r%name = '[O]-'//trim(Elecs(i)%name)//' folding in D'
          call rgn_print(r, seq_max = 10 )
       end if
    end do

    ! Clean-up
    call rgn_delete(r)

  end subroutine tbt_print_regions

  subroutine crt_El_priority(N_Elec,Elecs,pr,na_u,lasto,is_orb)
    use m_ts_electype
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    type(tRgn), intent(inout) :: pr ! Needs to be pre-allocated
    integer, intent(in) :: na_u, lasto(0:na_u)
    logical, intent(in) :: is_orb
    integer :: i, iEl, j
    type(tRgn) :: rgn
    
    ! Initialize
    pr%r = 0
    do iEl = 1 , N_Elec
       if ( is_orb ) then
          call rgn_copy(Elecs(iEl)%o_inD,rgn)
       else
          call rgn_Orb2Atom(Elecs(iEl)%o_inD,na_u,lasto,rgn)
       end if
       
       i = rgn%n
       do j = 1 , rgn%n
          pr%r(rgn%r(j)) = i
          i = i - 1
       end do
    end do

    call rgn_delete(rgn)

  end subroutine crt_El_priority

end module m_tbt_regions
