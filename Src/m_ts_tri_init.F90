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
! Nick Papior Andersen, 2014, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

! This particular solution method relies on solving the GF
! with the tri-diagonalization routine.
! This will leverage memory usage and also the execution time.

module m_ts_tri_init

  use precision, only : dp
  use m_region

  implicit none

  ! arrays for containing the tri-diagonal matrix part sizes
  type(tRgn) :: c_Tri

  public :: ts_tri_init
  public :: c_Tri

  private
  
contains

  subroutine ts_tri_init( dit, sparse_pattern , N_Elec, Elecs, &
       IsVolt, ucell, na_u, lasto , nsc, isc_off, method )

    use alloc, only : re_alloc, de_alloc
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
    integer, intent(in) :: nsc(3), isc_off(3,product(nsc))

    ! The method used to partition the BTD format
    integer, intent(in) :: method

    type(OrbitalDistribution) :: fdit
    type(Sparsity) :: tmpSp1, tmpSp2

    integer :: idx, no
    integer :: i, io, els, iEl, no_u_TS
    integer :: padding, worksize
    character(len=NAME_LEN) :: ctmp
    
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
    
    call crtSparsity_SC(tmpSp2,tmpSp1, UC = .TRUE. )

    ! point to the local (SIESTA-UC) sparsity pattern arrays
    call Sp_to_Spglobal(dit,tmpSp1,tmpSp2)
    
    call delete(tmpSp1)

    ! Get sorting method, we default to sort
    ! the BTD matrix according to the connection
    ! scheme of the first electrode.
    ctmp = fdf_get('TS.BTD.Sort',trim(Elecs(1)%name))
    if ( leqi(ctmp,'none') ) then

       ! do nothing, the pivoting array
       ! is already arranged correctly

    else ! the user *must* have supplied an electrode

       ! Figure out which electrode has been given
       iEl = 0
       do i = 1 , N_Elec
          if ( Elecs(i)%name == ctmp ) iEl = i
       end do
       if ( iEl == 0 ) call die('Could find the electrode in &
            &TS.BTD.Sort in the list of electrodes, &
            &please correct sorting method.')

       ! First we sort our pivot-region to minimize the elements
       ! in the BTD format
       i = Elecs(iEl)%idx_o
       call rgn_range(r_pvt,i,i-1+TotUsedOrbs(Elecs(iEl)))
       call rgn_copy(r_pvt,r_tmp)
       ! we sort the newly attached region
       call rgn_sp_sort(r_pvt, fdit, tmpSp2, r_tmp, R_SORT_MAX_FRONT )
       call rgn_copy(r_tmp,r_pvt)
       call rgn_delete(r_tmp)
       do 
          
          ! Create attached region starting from electrode 1
          call rgn_sp_connect(r_pvt, fdit, tmpSp2, r_tmp)

          if ( r_tmp%n == 0 .and. r_pvt%n /= no_u_TS ) then

             ! we need to ensure that all (even non-connected)
             ! orbitals are taken into account
             do i = 1 , nrows_g(tmpSp2)
                ! if it already exists, skip it
                if ( in_rgn(r_pvt,i) ) cycle
                if ( orb_type(i) /= TYP_BUFFER ) then
                   call rgn_range(r_tmp,i,i)
                   exit
                end if
             end do

             if ( r_tmp%n /= 0 ) then

                ! get the atomic index
                i = iaorb(r_tmp%r(1),lasto)
                
                if ( IONode ) then
                   write(*,'(/,a,i0,a)')'WARNING TS: &
                        &Random atom ',i,' has been added &
                        &due to non-completeness of the connectivity &
                        &graph.'
                   write(*,'(a)')'WARNING TS: Expect sub-optimal &
                        &BTD format.'
                end if
                call rgn_range(r_tmp,lasto(i-1)+1,lasto(i))
                ! Assert that none of the orbitals exist in the, region
                do i = 1 , r_tmp%n
                   if ( in_rgn(r_pvt,r_tmp%r(i)) ) then
                      call die('This is extremely difficult. &
                           &Please do not sort the BTD format as &
                           it cannot figure out what to do.')
                   end if
                end do

                call rgn_append(r_pvt, r_tmp, r_pvt)
                cycle

             end if
             
          end if
          
          ! If no additional orbitals are found, exit
          if ( r_tmp%n == 0 ) exit

          ! Append the newly found region that is connecting out to the
          ! full region
          call rgn_append(r_pvt, r_tmp, r_pvt)
          
          ! we sort the newly attached region
          call rgn_sp_sort(r_pvt, fdit, tmpSp2, r_tmp, R_SORT_MAX_BACK )
       
       end do
       
       call rgn_delete(r_tmp)

    end if
    if ( r_pvt%n /= nrows_g(sparse_pattern) - no_Buf ) then
       call die('Error in size estimation, the sparse pattern &
            &removal is erroneous')
    end if
       
    ! Now r_pvt contains the sorted device region according to
    ! electrode (1). (if asked for)
    ! From this we can generate a "better" tri-diagonal matrix than
    ! from the non-sorted one (provided that the user have provided
    ! the atoms in sub-optimal order).

    
    ! This works as creating a new sparsity deletes the previous
    ! and as it is referenced several times it will not be actually
    ! deleted...
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
       do io = 1 , no 
          Elecs(i)%o_inD%r(io) = r_pvt%r(r_tmp%r(io))
       end do

       ! Create the pivot table for the electrode scattering
       ! matrix
       call rgn_copy(Elecs(i)%o_inD,Elecs(i)%inDpvt)
       Elecs(i)%inDpvt%r = Elecs(i)%inDpvt%r - idx + 1

       ! we first create the super-set sparsity
       tmpSp1 = tmpSp2
       call crtSparsity_Union(fdit,tmpSp1, &
            idx,idx,no,no, tmpSp2)
    end do
    call rgn_delete(r_tmp)
    call delete(tmpSp1)

    ! Initialize the tri-diagonal container
    call rgn_delete(c_Tri)

#ifdef TRANSIESTA_DEBUG
    if(IONode)write(*,*)'Created TS-tri + elecs (4000)'
    call sp_to_file(4000,tmpSp2)
#endif

    if ( IONode ) &
         write(*,'(/,a)') 'transiesta: Determining an optimal tri-matrix...'

    ! Create a new tri-diagonal matrix, do it in parallel
    call ts_rgn2TriMat(N_Elec, Elecs, IsVolt, &
         fdit, tmpSp2, r_pvt, c_Tri%n, c_Tri%r, &
         method, 0, par = .true. )
    call delete(tmpSp2) ! clean up
    call delete(fdit)

    c_Tri%name = ' '

    if ( c_Tri%n < 2 ) then
       call die('Erroneous transiesta BTD format. &
            &Check with the developers')
    end if

    ! Calculate size of the tri-diagonal matrix
    els = nnzs_tri(c_Tri%n,c_Tri%r)
    
    if ( IONode ) then
       write(*,'(a)') 'transiesta: Established a near-optimal partition &
            &for the tri-diagonal matrix.'

       call rgn_print(c_Tri, name = 'BTD partitions' , &
            seq_max = 10 )

       write(*,'(a,i0,a,i0)') 'transiesta: Matrix elements in tri / full: ', &
            els,' / ',nrows_g(sparse_pattern)**2
    end if

    if ( .not. IsVolt ) return

    ! Get the padding for the array to hold the entire column
    call GFGGF_needed_worksize(c_Tri%n, c_Tri%r, &
         N_Elec, Elecs, padding, worksize)
    if ( IONode ) then
       write(*,'(a,i0)') 'transiesta: Padding + work elements: ', &
            padding + worksize
    end if

    ! Check that we can contain the full column
    if ( maxval(TotUsedOrbs(Elecs)) * no_u_TS > els ) then
       call die('Transiesta tri-diagonal partitioning is too good to &
            &perform GFGGF calculation. Cannot sustain the full column. &
            &Use the full/MUMPS TS.SolutionMethod instead')
    end if
    
  end subroutine ts_tri_init

end module m_ts_tri_init
