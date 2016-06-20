!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2014, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

! This particular solution method relies on solving the GF
! with the tri-diagonalization routine.
! This will leverage memory usage and also the execution time.

module m_ts_tri_init

  use precision, only : dp, i8b
  use m_region

  use m_ts_electype

  use m_ts_pivot

  implicit none

  ! arrays for containing the tri-diagonal matrix part sizes
  type(tRgn) :: c_Tri

  public :: ts_tri_init
  public :: ts_tri_analyze
  public :: c_Tri

  private
  
contains

  subroutine ts_tri_init( dit, sparse_pattern , N_Elec, Elecs, &
       IsVolt, ucell, na_u, xa, lasto , nsc, isc_off, method )

    use alloc, only : re_alloc, de_alloc
    use parallel, only : IONode
    use fdf, only : fdf_get, leqi
    use fdf_extra
#ifdef MPI
    use mpi_siesta, only : MPI_Comm_Self
#endif

    use class_OrbitalDistribution
    use class_Sparsity
    use create_Sparsity_Union
    use create_Sparsity_SC
    use m_sparsity_handling

    use m_pivot

    use m_ts_electype
    use m_ts_sparse, only : ts_sp_calculation

    use m_ts_method

    use m_ts_tri_common
    use m_ts_rgn2trimat

#ifdef TRANSIESTA_DEBUG
    use m_ts_debug
#endif

    type(OrbitalDistribution), intent(inout) :: dit
    type(Sparsity), intent(inout) :: sparse_pattern ! the local sparse pattern
    integer, intent(in) :: N_Elec
    type(Elec), intent(inout) :: Elecs(N_Elec)
    logical, intent(in) :: IsVolt
    real(dp), intent(in) :: ucell(3,3)
    integer, intent(in) :: na_u, lasto(0:na_u)
    real(dp), intent(in) :: xa(3,na_u)
    integer, intent(in) :: nsc(3), isc_off(3,product(nsc))

    ! The method used to partition the BTD format
    integer, intent(in) :: method

    type(OrbitalDistribution) :: fdit
    type(Sparsity) :: tmpSp1, tmpSp2

    integer :: idx, no
    integer :: i, io, els, iEl, no_u_TS
    integer :: padding, worksize
    character(len=NAME_LEN) :: csort
    
    ! Regions used for sorting the device region
    type(tRgn) :: r_tmp
    
    no_u_TS = nrows_g(sparse_pattern) - no_Buf

    ! In order to ensure that the electrodes are in the
    ! tri-diagonal sparsity pattern, we can easily create
    ! the full sparsity pattern with the electrodes included
    ! and then recreate the tri-diagonal sparsity pattern
    ! This is probably the crudest way of doing it.
#ifdef MPI
    call newDistribution(nrows_g(sparse_pattern),MPI_Comm_Self,fdit, &
         name='TranSIESTA UC distribution')
#else    
    call newDistribution(nrows_g(sparse_pattern),-1,fdit, &
         name='TranSIESTA UC distribution')
#endif

    ! We need to regenerate the sparsity pattern without the 
    ! electrode connections etc. this is because ts_sp_uc
    ! has different size dependent on the electrodes bulk settings.
    call ts_Sp_calculation(dit,sparse_pattern,N_Elec,Elecs, &
         ucell, nsc, isc_off, tmpSp2)
    
    call crtSparsity_SC(tmpSp2, tmpSp1, UC = .true. )

    ! point to the local (SIESTA-UC) sparsity pattern arrays
    call Sp_to_Spglobal(dit,tmpSp1,tmpSp2)

    ! This works as creating a new sparsity deletes the previous
    ! and as it is referenced several times it will not be actually
    ! deleted...
    do iEl = 1 , N_Elec

       idx = Elecs(iEl)%idx_o
       no  = TotUsedOrbs(Elecs(iEl))

       ! we first create the super-set sparsity
       tmpSp1 = tmpSp2
       call crtSparsity_Union(fdit,tmpSp1, &
            idx,idx,no,no, tmpSp2)

       ! Create the o_inD
       call rgn_range(Elecs(iEl)%o_inD,idx,idx+no-1)

    end do
    call delete(tmpSp1) ! clean up

#ifdef TRANSIESTA_DEBUG
    if(IONode)write(*,*)'Created TS-tri + elecs (4000)'
    call sp_to_file(4000,tmpSp2)
#endif

    ! Get sorting method, we default to sort
    ! the BTD matrix according to the connection
    ! scheme of the first electrode.
    csort = fdf_get('TS.BTD.Pivot','atom+'//trim(Elecs(1)%name))
    call ts_pivot( fdit, tmpSp2, &
         N_Elec, Elecs, &
         ucell, na_u, xa, lasto, &
         r_pvt, csort)

    ! In transiesta we can immediately calculate the
    ! tri-diagonal matrix.
    ! Then we can re-arrange the pivot indices in each block
    ! to be as consecutive in each electrode as possible.
    ! This will greatly speed up complex Gf.G.Gf matrix
    ! products
    call rgn_delete(c_Tri)

    ! Get the current sorting method
    if ( IONode ) &
         write(*,'(/,2a)') 'transiesta: Determining an optimal &
         &tri-matrix using: ',trim(csort)

    ! Create a new tri-diagonal matrix, do it in parallel
    call ts_rgn2TriMat(N_Elec, Elecs, IsVolt, &
         fdit, tmpSp2, r_pvt, c_Tri%n, c_Tri%r, &
         method, 0, par = .true. )
    call delete(tmpSp2) ! clean-up
    call delete(fdit)

    ! Sort the pivoting table for the electrodes
    ! such that we reduce the Gf.Gamma.Gf
    ! However, this also makes it easier to
    ! insert the self-energy as they become consecutive
    ! in index
    call ts_pivot_tri_sort_El(r_pvt,N_Elec,Elecs,c_Tri)

    if ( c_Tri%n < 2 ) then
       call die('Erroneous transiesta BTD format. &
            &Check with the developers')
    end if

    ! Recalculate number of orbitals in TS
    no_u_TS = nrows_g(sparse_pattern) - no_Buf

    if ( r_pvt%n /= no_u_TS ) then
       call die('Error in size estimation, the sparse pattern &
            &removal is erroneous')
    end if
       
    ! Now r_pvt contains the sorted device region according to
    ! electrode (1). (if asked for)
    ! From this we can generate a "better" tri-diagonal matrix than
    ! from the non-sorted one (provided that the user have provided
    ! the atoms in sub-optimal order).

    do i = 1 , N_Elec

       idx = Elecs(i)%idx_o
       no = TotUsedOrbs(Elecs(i))

       ! Create the pivoting table for the electrodes
       call rgn_init(r_tmp,no)
       do io = 1 , no 
          ! equals the io'th orbital index in the 
          !           TS region.  io == ts_i
          r_tmp%r(io) = rgn_pivot(r_pvt,idx-1+io)
       end do
       
       ! Sort it to be able to gather the indices
       ! in the correct order
       call rgn_sort(r_tmp)
       ! pivot the o_inD back
       call rgn_init(Elecs(i)%o_inD,no)
       call rgn_init(Elecs(i)%inDpvt,no)
       do io = 1 , no 
          Elecs(i)%o_inD%r(io)  = r_pvt%r(r_tmp%r(io))

          ! Create the pivot table for the electrode scattering matrix
          !  (1) = an electrode orbital which is seen first in the device
          ! Note that this array's meaning is different in TBtrans
          Elecs(i)%inDpvt%r(io) = Elecs(i)%o_inD%r(io) - idx + 1

       end do

    end do
    call rgn_delete(r_tmp)


    ! Calculate size of the tri-diagonal matrix
    els = nnzs_tri(c_Tri%n,c_Tri%r)
    
    if ( IONode ) then
       write(*,'(a)') 'transiesta: Established a near-optimal partition &
            &for the tri-diagonal matrix.'

       call rgn_print(c_Tri, name = 'BTD partitions' , &
            seq_max = 10 , collapse = .false. )

       write(*,'(a,i0,a,i0)') 'transiesta: Matrix elements in tri / full: ', &
            els,' / ',no_u_TS**2
    end if

    if ( .not. IsVolt ) return

    if ( ts_A_method == TS_BTD_A_COLUMN ) then
       ! Get the padding for the array to hold the entire column
       call GFGGF_needed_worksize(c_Tri%n, c_Tri%r, &
            N_Elec, Elecs, padding, worksize)
       if ( IONode ) then
          write(*,'(a,i0)') 'transiesta: Padding + work elements: ', &
               padding + worksize
       end if
    end if
    do iEl = 1 , N_Elec
       i = consecutive_Elec_orb(Elecs(iEl),r_pvt)
       if ( IONode ) then
          write(*,'(3a,i0,a)') 'transiesta: Electrode ',trim(Elecs(iEl)%name), &
               ' splitting Gamma ',i-1,' times.'
       end if
    end do

  end subroutine ts_tri_init

  subroutine ts_tri_analyze( dit, sparse_pattern , N_Elec, Elecs, &
       ucell, na_u, lasto , nsc, isc_off , method )
    
    use parallel, only : IONode
    use fdf, only : fdf_get, leqi
#ifdef MPI
    use mpi_siesta, only : MPI_Comm_Self
#endif

    use class_OrbitalDistribution
    use class_Sparsity
    use create_Sparsity_Union
    use create_Sparsity_SC
    use m_sparsity_handling

    use m_pivot
#ifdef GRAPHVIZ
    use m_pivot_methods, only : sp2graphviz
#endif

    use m_ts_electype
    use m_ts_sparse, only : ts_sp_calculation

    use m_ts_method

    use m_ts_tri_common
    use m_ts_rgn2trimat

#ifdef TRANSIESTA_DEBUG
    use m_ts_debug
#endif

    type(OrbitalDistribution), intent(inout) :: dit
    type(Sparsity), intent(inout) :: sparse_pattern ! the local sparse pattern
    integer, intent(in) :: N_Elec
    type(Elec), intent(inout) :: Elecs(N_Elec)
    real(dp), intent(in) :: ucell(3,3)
    integer, intent(in) :: na_u, lasto(0:na_u)
    integer, intent(in) :: nsc(3), isc_off(3,product(nsc))

    ! The method used to partition the BTD format
    integer, intent(in) :: method

    type(OrbitalDistribution) :: fdit
    type(Sparsity) :: tmpSp1, tmpSp2

    integer :: no_u_TS, i, iEl, no

    integer :: n, n_nzs
    integer, pointer, contiguous :: ncol(:), l_ptr(:), l_col(:)

    character(len=40), parameter :: fmt = '(/,''TS.BTD.Pivot '',a,''+'',a)'
    character(len=50) :: fmethod
    character(len=4) :: corb
    
    ! Regions used for sorting the device region
    type(tRgn) :: r_tmp, r_Els, r_El, full, priority
    integer :: orb_atom
    
    no_u_TS = nrows_g(sparse_pattern) - no_Buf

    ! In order to ensure that the electrodes are in the
    ! tri-diagonal sparsity pattern, we can easily create
    ! the full sparsity pattern with the electrodes included
    ! and then recreate the tri-diagonal sparsity pattern
    ! This is probably the crudest way of doing it.
#ifdef MPI
    call newDistribution(nrows_g(sparse_pattern),MPI_Comm_Self,fdit, &
         name='TranSIESTA UC distribution')
#else    
    call newDistribution(nrows_g(sparse_pattern),-1,fdit, &
         name='TranSIESTA UC distribution')
#endif

    ! We need to regenerate the sparsity pattern without the 
    ! electrode connections etc. this is because ts_sp_uc
    ! has different size dependent on the electrodes bulk settings.
    call ts_Sp_calculation(dit,sparse_pattern,N_Elec,Elecs, &
         ucell, nsc, isc_off, tmpSp2)
    
    call crtSparsity_SC(tmpSp2,tmpSp1, UC = .TRUE. )

    ! point to the local (SIESTA-UC) sparsity pattern arrays
    call Sp_to_Spglobal(dit,tmpSp1,tmpSp2)

    call rgn_delete(r_Els)
    do iEl = 1 , N_Elec

       i  = Elecs(iEl)%idx_o
       no = TotUsedOrbs(Elecs(iEl))

       ! we first create the super-set sparsity
       tmpSp1 = tmpSp2
       call crtSparsity_Union(fdit,tmpSp1, &
            i,i,no,no, tmpSp2)

       ! Create the o_inD
       call rgn_range(Elecs(iEl)%o_inD,i,i+no-1)
       call rgn_append(r_Els,Elecs(iEl)%o_inD,r_Els)

    end do

    ! Write out all pivoting etc. analysis steps
    if ( IONode ) write(*,'(/,a)') 'transiesta: BTD analysis'

    ! Make a copy
    tmpSp1 = tmpSp2

    ! Attach the sparsity pattern of the orbitals
    ! later (tmpSp2) may be atom
    call attach(tmpSp1, n_col = ncol, list_ptr = l_ptr, &
         list_col = l_col , nrows_g = no , nnzs = n_nzs )

    fmethod = 'orb+none'
    if ( IONode ) write(*,fmt) 'orb','none'
    call tri(r_pvt)

    orb_atom_switch: do orb_atom = 1 , 2
    ! The user can skip the orbital analysis if it takes too long
    if ( orb_atom == 1 ) then
       ! We default to only looking at the atomic sparsity
       ! pattern. This is *much* faster and does provide
       ! a very near optimal sparse pattern. 
       ! The user can select to do both.
       if ( leqi(fdf_get('TS.BTD.Analyze','atom'),'atom') ) cycle
       corb = 'orb'

       call rgn_copy(r_pvt,full)

    else
       corb = 'atom'

       ! Convert the sparsity pattern to the atom
       call SpOrb_to_SpAtom(fdit,tmpSp1,na_u,lasto,tmpSp2)
       ! *** the distribution will always
       !     be bigger than for the atoms, hence we need
       !     not re-construct it ***

       ! Reduce the searching place of atoms
       call rgn_copy(r_aC,full)

       ! Prepare total region with all electrodes
       call rgn_orb2atom(r_Els,na_u,lasto,r_tmp)
       call rgn_copy(r_tmp,r_Els)

    end if

    ! Sort the device region
    call rgn_sort(full)

    n = nrows_g(tmpSp2)

    ! Create priority list
    call rgn_init(priority,no)
    call crt_El_priority(N_Elec,Elecs,priority, &
         na_u,lasto,is_orb = orb_atom == 1 )

#ifdef GRAPHVIZ
    ! Attach sparsity pattern to current designation ('atom' == atom sp)
    call attach(tmpSp2, n_col = ncol, list_ptr = l_ptr, &
         list_col = l_col , nnzs = n_nzs )
    call rgn_init(r_El,n)
    r_El%r(:) = 0
    do i = 1 , n
       if ( orb_atom == 1 ) then
          r_El%r(i) = orb_type(i)
       else
          r_El%r(i) = atom_type(i)
       end if
    end do
    if ( IONode ) &
         call sp2graphviz('GRAPHVIZ_'//trim(corb)//'.gv', &
         n,n_nzs,ncol,l_ptr,l_col, types = r_El%r )
#endif
  
    ! Attach the sparsity pattern of the orbitals
    call attach(tmpSp1, n_col = ncol, list_ptr = l_ptr, &
         list_col = l_col , nrows_g = no , nnzs = n_nzs )


    ! *** Start analyzing sparsity pattern


    do iEl = 1 , N_Elec
       fmethod = trim(corb)//'+'//trim(Elecs(iEl)%name)
       if ( IONode ) write(*,fmt) trim(corb),trim(Elecs(iEl)%name)
       if ( orb_atom == 1 ) then
          call rgn_copy(Elecs(iEl)%o_inD,r_El)
       else
          call rgn_copy(Elecs(iEl)%o_inD,r_tmp)
          call rgn_orb2atom(r_tmp,na_u,lasto,r_El)
       end if
       call rgn_copy(r_El,r_tmp)
       ! we sort the newly attached region
       call rgn_sp_sort(r_El, fdit, tmpSp2, r_tmp, R_SORT_MAX_FRONT )
       call rgn_copy(r_tmp,r_El)
       call rgn_delete(r_tmp)
       do 
          call rgn_sp_connect(r_El, fdit, tmpSp2, r_tmp)
          if ( r_tmp%n == 0 .and. r_El%n /= no_u_TS ) then
             do i = 1 , nrows_g(tmpSp2)
                if ( in_rgn(r_El,i) ) cycle
                if ( orb_atom == 1 ) then
                   if ( orb_type(i) /= TYP_BUFFER ) then
                      call rgn_range(r_tmp,i,i)
                      exit
                   end if
                else
                   if ( atom_type(i) /= TYP_BUFFER ) then
                      call rgn_range(r_tmp,i,i)
                      exit
                   end if
                end if
             end do
             if ( r_tmp%n /= 0 ) then
                do i = 1 , r_tmp%n
                   if ( in_rgn(r_El,r_tmp%r(i)) ) then
                      call die('This is extremely difficult. &
                           &Please do not sort the BTD format as &
                           it cannot figure out what to do.')
                   end if
                end do
                call rgn_append(r_El, r_tmp, r_El)
                cycle
             end if
          end if
          if ( r_tmp%n == 0 ) exit
          call rgn_append(r_El, r_tmp, r_El)
          call rgn_sp_sort(r_El, fdit, tmpSp2, r_tmp, R_SORT_MAX_BACK )
       end do
       if ( orb_atom == 1 ) then
          call tri(r_El)
       else
          call rgn_atom2orb(r_El,na_u,lasto,r_tmp)
          call tri(r_tmp)
       end if
    end do

    fmethod = trim(corb)//'+CM'
    if ( IONode ) write(*,fmt) trim(corb),'CM'
    call sp_pvt(n,tmpSp2,r_tmp, PVT_CUTHILL_MCKEE, sub = full, start = r_Els)
    if ( orb_atom == 1 ) then
       call tri(r_tmp)
    else
       call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
       call tri(r_El)
    end if

    fmethod = trim(corb)//'+rev-CM'
    if ( IONode ) write(*,fmt) trim(corb),'rev-CM'
    call rgn_reverse(r_tmp)
    if ( orb_atom == 1 ) then
       call tri(r_tmp)
    else
       call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
       call tri(r_El)
    end if

    fmethod = trim(corb)//'+CM+priority'
    if ( IONode ) write(*,fmt) trim(corb),'CM+priority'
    call sp_pvt(n,tmpSp2,r_tmp, PVT_CUTHILL_MCKEE, sub = full, start = r_Els, &
         priority = priority%r)
    if ( orb_atom == 1 ) then
       call tri(r_tmp)
    else
       call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
       call tri(r_El)
    end if

    fmethod = trim(corb)//'+rev-CM+priority'
    if ( IONode ) write(*,fmt) trim(corb),'rev-CM+priority'
    call rgn_reverse(r_tmp)
    if ( orb_atom == 1 ) then
       call tri(r_tmp)
    else
       call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
       call tri(r_El)
    end if
    
    fmethod = trim(corb)//'+GPS'
    if ( IONode ) write(*,fmt) trim(corb),'GPS'
    call sp_pvt(n,tmpSp2,r_tmp, PVT_GPS, sub = full)
    if ( orb_atom == 1 ) then
       call tri(r_tmp)
    else
       call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
       call tri(r_El)
    end if

    fmethod = trim(corb)//'+rev-GPS'
    if ( IONode ) write(*,fmt) trim(corb),'rev-GPS'
    call rgn_reverse(r_tmp)
    if ( orb_atom == 1 ) then
       call tri(r_tmp)
    else
       call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
       call tri(r_El)
    end if

    fmethod = trim(corb)//'+GPS+priority'
    if ( IONode ) write(*,fmt) trim(corb),'GPS+priority'
    call sp_pvt(n,tmpSp2,r_tmp, PVT_GPS, sub = full, priority = priority%r)
    if ( orb_atom == 1 ) then
       call tri(r_tmp)
    else
       call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
       call tri(r_El)
    end if

    fmethod = trim(corb)//'+rev-GPS+priority'
    if ( IONode ) write(*,fmt) trim(corb),'rev-GPS+priority'
    call rgn_reverse(r_tmp)
    if ( orb_atom == 1 ) then
       call tri(r_tmp)
    else
       call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
       call tri(r_El)
    end if

    fmethod = trim(corb)//'+PCG'
    if ( IONode ) write(*,fmt) trim(corb),'PCG'
    call sp_pvt(n,tmpSp2,r_tmp, PVT_PCG, sub = full)
    if ( orb_atom == 1 ) then
       call tri(r_tmp)
    else
       call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
       call tri(r_El)
    end if

    fmethod = trim(corb)//'+rev-PCG'
    if ( IONode ) write(*,fmt) trim(corb),'rev-PCG'
    call rgn_reverse(r_tmp)
    if ( orb_atom == 1 ) then
       call tri(r_tmp)
    else
       call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
       call tri(r_El)
    end if

    fmethod = trim(corb)//'+PCG+priority'
    if ( IONode ) write(*,fmt) trim(corb),'PCG+priority'
    call sp_pvt(n,tmpSp2,r_tmp, PVT_PCG, sub = full, priority = priority%r)
    if ( orb_atom == 1 ) then
       call tri(r_tmp)
    else
       call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
       call tri(r_El)
    end if

    fmethod = trim(corb)//'+rev-PCG+priority'
    if ( IONode ) write(*,fmt) trim(corb),'rev-PCG+priority'
    call rgn_reverse(r_tmp)
    if ( orb_atom == 1 ) then
       call tri(r_tmp)
    else
       call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
       call tri(r_El)
    end if

    fmethod = trim(corb)//'+GGPS'
    if ( IONode ) write(*,fmt) trim(corb),'GGPS'
    call sp_pvt(n,tmpSp2,r_tmp, PVT_GGPS, sub = full)
    if ( orb_atom == 1 ) then
       call tri(r_tmp)
    else
       call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
       call tri(r_El)
    end if

    fmethod = trim(corb)//'+rev-GGPS'
    if ( IONode ) write(*,fmt) trim(corb),'rev-GGPS'
    call rgn_reverse(r_tmp)
    if ( orb_atom == 1 ) then
       call tri(r_tmp)
    else
       call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
       call tri(r_El)
    end if

    fmethod = trim(corb)//'+GGPS+priority'
    if ( IONode ) write(*,fmt) trim(corb),'GGPS+priority'
    call sp_pvt(n,tmpSp2,r_tmp, PVT_GGPS, sub = full, priority = priority%r)
    if ( orb_atom == 1 ) then
       call tri(r_tmp)
    else
       call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
       call tri(r_El)
    end if

    fmethod = trim(corb)//'+rev-GGPS+priority'
    if ( IONode ) write(*,fmt) trim(corb),'rev-GGPS+priority'
    call rgn_reverse(r_tmp)
    if ( orb_atom == 1 ) then
       call tri(r_tmp)
    else
       call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
       call tri(r_El)
    end if

#ifdef SIESTA__METIS
    fmethod = trim(corb)//'+metis'
    if ( IONode ) write(*,fmt) trim(corb),'metis'
    call sp_pvt(n,tmpSp2,r_tmp, PVT_METIS, sub = full)
    if ( orb_atom == 1 ) then
       call tri(r_tmp)
    else
       call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
       call tri(r_El)
    end if

    fmethod = trim(corb)//'+metis+priority'
    if ( IONode ) write(*,fmt) trim(corb),'metis+priority'
    call sp_pvt(n,tmpSp2,r_tmp, PVT_METIS, sub = full, &
         priority = priority%r)
    if ( orb_atom == 1 ) then
       call tri(r_tmp)
    else
       call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
       call tri(r_El)
    end if
#endif

    end do orb_atom_switch

    call rgn_delete(r_tmp,r_Els,r_El,full,priority)

    call delete(tmpSp1) ! clean up
    call delete(tmpSp2)
    call delete(fdit)

    if ( IONode ) write(*,*) ! new-line
    
  contains

    ! Print out all relevant information for this
    ! pivoting scheme
    subroutine tri(r_pvt)
      use m_pivot_methods, only : bandwidth, profile
      use fdf, only : fdf_overwrite
      type(tRgn), intent(inout) :: r_pvt

      integer :: bw, els, pad, work, i
      integer(i8b) :: prof
      type(tRgn) :: ctri
      character(len=132) :: fname

      bw   = bandwidth(no,n_nzs,ncol,l_ptr,l_col,r_pvt)
      prof = profile(no,n_nzs,ncol,l_ptr,l_col,r_pvt)
      if ( IONode ) then
         write(*,'(tr3,a,t23,i10,/,tr3,a,t13,i20)') &
              'Bandwidth: ',bw,'Profile: ',prof
      end if

      ! Only if it is defined
      fname = fdf_get('TS.BTD.Output',' ')
      if ( len_trim(fname) > 0 ) then
         fname = 'TS.BTD.Output '//trim(fmethod)
         call fdf_overwrite(fname)
      end if

      ! Create a new tri-diagonal matrix, do it in parallel
      call ts_rgn2TriMat(N_Elec, Elecs, .true., &
           fdit, tmpSp1, r_pvt, ctri%n, ctri%r, &
           method, 0, par = .true. )

      ! Sort the pivoting table for the electrodes
      ! such that we reduce the Gf.Gamma.Gf
      ! However, this also makes it easier to
      ! insert the self-energy as they become consecutive
      ! in index, all-in-all, win-win!
      call ts_pivot_tri_sort_El(r_pvt,N_Elec,Elecs,ctri)
      
      ! Calculate size of the tri-diagonal matrix
      els = nnzs_tri(ctri%n,ctri%r)
    
      if ( IONode ) then
         call rgn_print(ctri, name = 'BTD partitions' , &
              seq_max = 10 , indent = 3 , collapse = .false. )
         
         pad = no_u_TS ** 2
         write(*,'(tr3,a,i0,'' / '',f9.5)') &
              'Matrix elements in tri / % of full: ', &
              els , real(els,dp)/real(pad,dp) * 100._dp
      end if

      if ( ts_A_method == TS_BTD_A_COLUMN ) then
         ! Get the padding for the array to hold the entire column
         call GFGGF_needed_worksize(ctri%n, ctri%r, &
              N_Elec, Elecs, pad, work)
      else
         pad = 0
         work = 0
      end if
      prof = pad + work + els * 2 ! total mem
      if ( IONode ) then
         write(*,'(tr3,a,t30,i20)') 'BTD x 2 + padding + work: ', prof
         write(*,'(tr3,a,t39,f8.2,a)') 'Rough estimation of MEMORY: ', &
              real(prof,dp) * 16._dp / 1024._dp ** 2,' MB'
      end if
      do i = 1 , N_Elec
         work = consecutive_Elec_orb(Elecs(i),r_pvt)
         if ( IONode ) then
            write(*,'(tr3,2a,t39,i0)') trim(Elecs(i)%name), &
                 ' splitting Gamma:',work-1
         end if
      end do

      call rgn_delete(ctri)
      
    end subroutine tri

  end subroutine ts_tri_analyze

end module m_ts_tri_init