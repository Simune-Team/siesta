! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

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

  public :: tbt_init_regions
  public :: tbt_region_options
  public :: tbt_print_regions

contains

  subroutine tbt_init_regions(N_Elec, Elecs, cell, na_u, xa, lasto, &
       dit, sp, &
       nsc, isc_off)

    use fdf
    use fdf_extra
    use parallel, only : Node, Nodes
#ifdef MPI
    use mpi_siesta, only : MPI_Comm_World, MPI_Barrier
#endif

    use m_pivot
#ifdef GRAPHVIZ
    use m_pivot_methods, only : sp2graphviz
#endif

    use geom_helper, only : iaorb
    use intrinsic_missing, only : SPC_PROJ, VNORM, VEC_PROJ
    use m_ts_electype
    use m_ts_method, only : ts_init_regions
    use m_ts_method, only : atom_type, TYP_DEVICE, TYP_BUFFER

    use m_ts_sparse, only : ts_Sparsity_Global
    use m_ts_pivot, only : ts_pivot

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
    ! The atomic coordinates
    real(dp), intent(in) :: xa(3,na_u)
    ! The distribution for the sparsity pattern
    type(OrbitalDistribution), intent(in) :: dit
    ! The sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! The supercell information
    integer, intent(in) :: nsc, isc_off(3,nsc)

    integer :: iEl, na

    ! A temporary sparsity pattern
    type(Sparsity) :: sp_tmp

    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline => null()
    character(len=50) :: g, csort
    integer :: i, ia, ia1, ia2, no_u
    type(tRgn) :: r_tmp, r_tmp2, r_tmp3, r_Els, priority
    real(dp) :: tmp

    no_u = lasto(na_u)
    
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
       call rgn_range(r_oEl_alone(iEl), ia1, ia2)
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
          if ( leqi(g,'atom') .or. leqi(g,'position') ) then
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
       call rgn_range(r_aDev,1,na_u)

    end if


    ! We immediately remove the buffer and electrode atoms
    if ( r_aBuf%n > 0 ) then
       call rgn_complement(r_aBuf,r_aDev,r_aDev)
    end if
    do iEl = 1 , N_Elec
       call rgn_complement(r_aEl_alone(iEl),r_aDev,r_aDev)
    end do


    if ( r_aDev%n == 0 ) then
       call die('Zero atoms are in the device region...???')
    end if

    ! Create device region
    call rgn_Atom2Orb(r_aDev,na_u,lasto,r_oDev)

    if ( Node == 0 ) then
       write(*,'(/,a)')'tbtrans: Analyzing electrode sparsity &
            &pattern to create optimal tri-diagonal blocks...'
    end if

    ! In case the user wants "a correct DOS"
    ! in this region, we extend it
    if ( fdf_get('TBT.Atoms.Device.Connect',.false.) ) then

       ! TBTrans will truncate connections at electrode interfaces.
       call rgn_sp_connect(r_oDev, dit, sp, r_tmp)
       if ( r_tmp%n == 0 ) &
            call die('No orbitals connect to the specified device &
            &region. This is not allowed.')

       ! Convert connecting region to atoms
       call rgn_Orb2Atom(r_tmp,na_u,lasto,r_oDev)
       call rgn_delete(r_tmp)
       ! Remove buffer atoms (in case the electrode is too small)
       if ( r_aBuf%n > 0 ) then
          call rgn_complement(r_aBuf,r_oDev,r_oDev)
       end if
       ! Remove all electrodes from the region
       do iEl = 1 , N_Elec
          call rgn_complement(r_aEl_alone(iEl),r_oDev,r_oDev)
       end do

       ! Append to device region
       if ( r_oDev%n > 0 ) then
          call rgn_append(r_aDev,r_oDev,r_aDev)
       else if ( Node == 0 ) then
          ! It only connects to electrodes
          write(*,'(a)')'tbtrans: Device regions &
               &connects directly with electrodes'
          write(*,'(a)')'tbtrans: If the overlap is large this might &
               &produce spurious effects in DOS calculations'
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
       i = Elecs(iEl)%pvt(Elecs(iEl)%t_dir)
       
       ! Remove connections from this electrode across the boundary...
       call Sp_remove_crossterms(dit,sp,nsc,isc_off, &
            i, &
            sp, r = r_tmp2)

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

    ! Create the electrode down-folding regions.
    ! Note that sorting according to the more advanced methods
    ! is not directly applicable as the methods involve non-stringent
    ! ending elements.

    do iEl = 1 , N_Elec

       if ( mod(iEl-1,Nodes) /= Node ) cycle

       ! Create pivoting region (except device)
       call rgn_range(r_tmp,1,no_u)
       if ( r_oBuf%n > 0 ) then
          call rgn_complement(r_oBuf,r_tmp,r_tmp)
       end if
       call rgn_complement(r_oDev,r_tmp,r_oEl(iEl))

       ! Ensure the initial electrode orbitals are allocated
       call rgn_copy(r_oEl_alone(iEl),Elecs(iEl)%o_inD)

       ! Sort according to the connectivity of the electrode
       ! This will also reduce the pivoting table (r_oEl) to
       ! _only_ have orbitals from the electrode up to the
       ! device region.
       ! First get the actual sub region that connects
       ! from the electrode to the device
       g = 'atom+'//trim(Elecs(iEl)%name)
       call ts_pivot( dit, sp, &
            1, Elecs(iEl:iEl), &
            cell, na_u, xa, lasto, &
            r_oEl(iEl), g, extend = .false.)
       csort = fdf_get('TBT.BTD.Pivot.Elecs',trim(g))
       csort = fdf_get('TBT.BTD.Pivot.Elec.'//&
            trim(Elecs(iEl)%name),trim(csort))
       if ( .not. leqi(g,csort) ) then
          call ts_pivot( dit, sp, &
               1, Elecs(iEl:iEl), &
               cell, na_u, xa, lasto, &
               r_oEl(iEl), csort)
          ! if the majority of the electrode pivoting
          ! elements are in the upper half, then reverse
          tmp = sum(rgn_pivot(r_oEl(iEl),r_oEl_alone(iEl)%r))
          ! average position in pivoting array
          tmp = tmp / r_oEl_alone(iEl)%n
          if ( tmp > 0.5_dp * r_oEl(iEl)%n ) then
             ! the actual electrode position
             ! is primarily placed closer to the device,
             ! hence we reverse the list to get a
             ! better sorting
             call rgn_reverse(r_oEl(iEl))
          end if
          
       end if
       
       ! This aligns the atoms in the same way the orbitals 
       ! introduce the atoms.
       call rgn_Orb2Atom(r_oEl(iEl), na_u, lasto , r_aEl(iEl))

       ! Create the region that connects the electrode-followed
       ! region to the central region
       call rgn_sp_connect(r_oEl(iEl), dit, sp, Elecs(iEl)%o_inD)
       ! To streamline the memory layout in the tri-diagonal
       ! blocks we might as well sort this
       call rgn_sort(Elecs(iEl)%o_inD)
       ! Append the found region that is connecting out to the device region
       call rgn_append(r_oEl(iEl), Elecs(iEl)%o_inD, r_oElpD(iEl))

       ! We now know how many orbitals that we are down-folding the 
       ! electrode self-energy to
       if ( Elecs(iEl)%o_inD%n == 0 ) then
          ! Print-out for debugging purposes
          call rgn_print(r_oDev)
          call rgn_print(r_oEl(iEl))
          call die('The electrode down-folding region is 0 in the device &
               &region. Please expand your device region.')
       end if

       ! Create the atom equivalent regions
       call rgn_Orb2Atom(r_oElpD(iEl) , na_u, lasto, r_aElpD(iEl) )

    end do

    ! Possibly deleting it
    call rgn_delete(r_tmp2)
    do iEl = 1 , N_Elec

#ifdef MPI
       i = mod(iEl-1,Nodes)
       ! Bcast the region
       call rgn_MPI_Bcast(r_oEl(iEl),i)
       call rgn_MPI_Bcast(r_aEl(iEl),i)
       call rgn_MPI_Bcast(Elecs(iEl)%o_inD,i)
       call rgn_MPI_Bcast(r_oElpD(iEl),i)
       call rgn_MPI_Bcast(r_aElpD(iEl),i)
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

       ! Prepare the inDpvt array (will be filled in tri_init
       ! as we sort each block individually)
       call rgn_copy(Elecs(iEl)%o_inD,Elecs(iEl)%inDpvt)
       
       ! Check that the electrode down-folded self-energy
       ! is fully contained
       ia1 = minval(rgn_pivot(r_oDev,Elecs(iEl)%o_inD%r))
       if ( ia1 <= 0 ) then
          call die('A downfolded region is not existing. Programming error')
       end if

       ! Collect all electrode down-fold regions into one
       call rgn_union(r_tmp2,r_oEl(iEl),r_tmp2)

#ifdef GRAPHVIZ
       ! If the user requests GRAPHVIZ output
       if ( fdf_get('TBT.BTD.Pivot.Graphviz',.false.) .and. Node == 0 ) then
          csort = trim(Elecs(iEl)%name) // '.pvt'
          call sp2graphviz(csort,sp,pvt=r_oEl(iEl))
       end if
#endif

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

    ! Prepare the check-regions

    csort = 'atom+'//trim(Elecs(1)%name)
    csort = fdf_get('TS.BTD.Pivot',trim(csort))
    csort = fdf_get('TBT.BTD.Pivot',trim(csort))
    csort = fdf_get('TBT.BTD.Pivot.Device',trim(csort))
    call ts_pivot( dit, sp, &
         N_Elec, Elecs, &
         cell, na_u, xa, lasto, &
         r_oDev, csort)

    ! Print out what we found
    if ( Node == 0 ) then
       write(*,'(a)')'tbtrans: BTD pivoting scheme in device: '//trim(csort)
    end if

    ! Check that there is no overlap with the other regions
    if ( rgn_overlaps(r_tmp2, r_oDev) ) then
       call rgn_print(r_oDev)
       call rgn_print(r_tmp2)
       print *,'Overlapping device, down-folding region(s)...'
       call die('tbt_regions: Error in programming, electrode down')
    end if

    if ( rgn_overlaps(r_oBuf, r_oDev) ) then
       call rgn_print(r_oDev)
       call rgn_print(r_oBuf)
       print *,'Overlapping device and buffer region...'
       call die('tbt_regions: Error in programming, buffer')
    end if

    ! Ensure that the number of device orbitals + electrode downfolding
    ! + buffer orbitals equal the full system
    if ( r_tmp2%n + r_oDev%n + r_oBuf%n /= no_u ) then
       r_tmp2%name = 'Electrodes'
       r_oDev%name = 'Device'
       call rgn_print(r_oBuf)
       call rgn_print(r_oDev)
       call rgn_print(r_tmp2)
       call die('tbt_regions: Error in programming, total')
    end if

    call rgn_Orb2Atom(r_oDev,na_u,lasto,r_aDev)
    r_oDev%name = '[O]-device'
    r_aDev%name = '[A]-device'
#ifdef GRAPHVIZ
    ! If the user requests GRAPHVIZ output
    if ( fdf_get('TBT.BTD.Pivot.Graphviz',.false.) .and. Node == 0 ) then
       csort = 'device.pvt'
       call sp2graphviz(csort,sp,pvt=r_oDev)
    end if
#endif

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

    end if

  contains

    function sort_contain(str,name) result(contain)
      use m_char, only : lcase
      character(len=*), intent(in) :: str, name
      logical :: contain

      character(len=len(str)) :: lstr
      character(len=len_trim(name)) :: lname

      integer :: i

      contain = .false.

      lstr = lcase(str)
      lname = lcase(trim(name))

      ! check whether it is in this stuff
      i = index(lstr,lname)
      if ( i > 1 ) then
         contain = scan(lstr(i-1:i-1),'+ ') == 1
         i = i + len_trim(lname)
         if ( i <= len(str) ) then
            contain = contain .and. scan(lstr(i:i),'+ ') == 1
         end if
      else if ( i == 1 ) then
         i = i + len_trim(lname)
         if ( i <= len(str) ) then
            contain = scan(lstr(i:i),'+ ') == 1
         else
            ! it was found and the string is too short
            contain = .true.
         end if
      end if
      
    end function sort_contain

  end subroutine tbt_init_regions


  subroutine tbt_region_options( save_DATA )
    use dictionary
#ifdef MPI
    use mpi_siesta, only : MPI_Comm_Self
#endif
    use m_sparsity_handling, only : Sp_retain_rgn

    type(dict), intent(in) :: save_DATA
#ifdef NCDF_4
    type(OrbitalDistribution) :: fdit

    integer :: no_u
#endif

    ! Make sure to initialize the device region
    ! sparsity pattern
    call delete(sp_dev)
#ifdef NCDF_4
    if ( ('orb-current' .in. save_DATA) .or. &
         ('proj-orb-current' .in. save_DATA) ) then

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
    integer :: i
    type(tRgn) :: r

    if ( Node /= 0 ) return

    if ( verbosity < 3 ) return

    ! Print out the buffer regions
    if ( r_aBuf%n > 0 ) then
       call local_print(r_aBuf,.false.)
       call local_print(r_oBuf,.true.)
    end if

    ! Print out the device region
    write(*,'(a,i0)')'tbtrans: # of device region orbitals: ',r_oDev%n

    call local_print(r_aDev,.false.)
    call local_print(r_oDev,.true.)

    ! Print out all the electrodes + their projection region
    do i = 1 , N_Elec
       write(*,*) ! new-line
       write(*,'(3a,i0)')'tbtrans: # of ',trim(Elecs(i)%name), &
            ' scattering orbitals: ',Elecs(i)%o_inD%n
       write(*,'(3a,i0)')'tbtrans: # of ',trim(Elecs(i)%name), &
            ' down-folding orbitals: ',r_oElpD(i)%n
       call local_print(r_aEl(i),.false.)
       call local_print(r_oEl(i),.true.)
       if ( verbosity > 3 ) then
          call rgn_intersection(r_aElpD(i),r_aDev,r)
          r%name = '[A]-'//trim(Elecs(i)%name)//' folding in D'
          call local_print(r, .false. )
       end if
       if ( verbosity > 7 ) then
          call rgn_intersection(r_oElpD(i),r_oDev,r)
          r%name = '[O]-'//trim(Elecs(i)%name)//' folding in D'
          call local_print(r,.true.)
       end if
    end do

    ! Clean-up
    call rgn_delete(r)

  contains
    
    subroutine local_print(r_in,is_orb)
      type(tRgn), intent(in) :: r_in
      logical, intent(in) :: is_orb
      type(tRgn) :: r
      integer :: seq, mid, high

      if ( is_orb ) then
         seq = 10
         mid = 7
         high = 9
      else
         seq = 12
         mid = 4
         high = 6
      end if

      if ( verbosity > high ) then
         ! print-unsorted
         call rgn_print(r_in, seq_max = seq )
      else if ( verbosity > mid ) then
         call rgn_copy(r_in,r)
         call rgn_sort(r)
         call rgn_print(r, seq_max = seq )
         call rgn_delete(r)
      end if

    end subroutine local_print

  end subroutine tbt_print_regions

end module m_tbt_regions
