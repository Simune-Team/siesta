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

  use precision, only : dp, i8b
  use m_region

  implicit none

  ! arrays for containing the tri-diagonal matrix part sizes
  type(tRgn) :: c_Tri

  public :: ts_tri_init
  public :: ts_tri_analyze
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
    integer, intent(in) :: nsc(3), isc_off(3,product(nsc))

    ! The method used to partition the BTD format
    integer, intent(in) :: method

    type(OrbitalDistribution) :: fdit
    type(Sparsity) :: tmpSp1, tmpSp2

    integer :: idx, no
    integer :: i, io, els, iEl, no_u_TS
    integer :: padding, worksize
    character(len=NAME_LEN) :: ctmp
    logical :: sort_orb
    
    ! Regions used for sorting the device region
    type(tRgn) :: r_tmp, r_Els
    
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
    end do
    
    tmpSp1 = tmpSp2

    ! Get sorting method, we default to sort
    ! the BTD matrix according to the connection
    ! scheme of the first electrode.
    ctmp = fdf_get('TS.BTD.Sort',trim(Elecs(1)%name))
    ! we default to do atomic sorting, faster and very consistent
    sort_orb = (index(ctmp,'orb') > 0)
    if ( sort_orb ) then
       ! the specification is orb something
       i = index(ctmp,'orb')
       ctmp(i:i+2) = ' '
       ! Prepare total region with all electrodes
       call rgn_range(r_Els,Elecs(1)%idx_o, &
            Elecs(1)%idx_o-1+TotUsedOrbs(Elecs(1)))
       do i = 2 , N_Elec
          call rgn_range(r_tmp,Elecs(i)%idx_o, &
               Elecs(i)%idx_o-1+TotUsedOrbs(Elecs(i)))
          call rgn_append(r_Els,r_tmp,r_Els)
       end do
       call rgn_copy(r_pvt,r_tmp)

       no_u_TS = nrows_g(tmpSp2) - no_Buf

    else
       ! in case the user supplies 'atom+something'
       i = index(ctmp,'atom')
       if ( i > 0 ) ctmp(i:i+3) = ' '
       ! We are doing atomic comparison
       call rgn_range(r_Els,Elecs(1)%idx_a, &
            Elecs(1)%idx_a-1+TotUsedAtoms(Elecs(1)))
       do i = 2 , N_Elec
          call rgn_range(r_tmp,Elecs(i)%idx_a, &
               Elecs(i)%idx_a-1+TotUsedAtoms(Elecs(i)))
          call rgn_append(r_Els,r_tmp,r_Els)
       end do
       call rgn_copy(r_aC,r_tmp)
       call SpOrb_to_SpAtom(fdit,tmpSp1,na_u,lasto,tmpSp2)
       ! *** the distribution will always
       !     be bigger than for the atoms, hence we need
       !     not re-construct it
       ! ***

       no_u_TS = nrows_g(tmpSp2) - na_Buf

    end if
    ! Remove all '+' in the string
    i = index(ctmp,'+')
    if ( i > 0 ) ctmp(i:i) = ' '
    ! Left adjust the string
    ctmp = ADJUSTL(ctmp)
    if ( leqi(ctmp,'none') ) then

       ! do nothing, the pivoting array
       ! is already arranged correctly

    else if ( leqi(ctmp,'CM') ) then

       call sp_pvt(tmpSp2,r_pvt, PVT_CUTHILL_MCKEE, r_tmp, start = r_Els)

    else if ( leqi(ctmp,'rev-CM') ) then

       call sp_pvt(tmpSp2,r_pvt, PVT_REV_CUTHILL_MCKEE, r_tmp, start = r_Els)

    else if ( leqi(ctmp,'GPS') ) then

       call sp_pvt(tmpSp2,r_pvt, PVT_GPS, r_tmp)

    else if ( leqi(ctmp,'rev-GPS') ) then

       call sp_pvt(tmpSp2,r_pvt, PVT_REV_GPS, r_tmp)

    else if ( leqi(ctmp,'GGPS') ) then

       call sp_pvt(tmpSp2,r_pvt, PVT_GGPS, r_tmp)

    else if ( leqi(ctmp,'rev-GGPS') ) then

       call sp_pvt(tmpSp2,r_pvt, PVT_REV_GGPS, r_tmp)

    else ! the user *must* have supplied an electrode

       ! Figure out which electrode has been given
       iEl = 0
       do i = 1 , N_Elec
          if ( Elecs(i)%name == ctmp ) iEl = i
       end do
       if ( iEl == 0 ) then
          print *,ctmp
          call die('Could find the electrode in &
               &TS.BTD.Sort in the list of electrodes, &
               &please correct sorting method.')
       end if

       ! First we sort our pivot-region to minimize the elements
       ! in the BTD format
       if ( sort_orb ) then
          i = Elecs(iEl)%idx_o
          call rgn_range(r_pvt,i,i-1+TotUsedOrbs(Elecs(iEl)))
       else
          i = Elecs(iEl)%idx_a
          call rgn_range(r_pvt,i,i-1+TotUsedAtoms(Elecs(iEl)))
       end if
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
                if ( sort_orb ) then
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

                if ( sort_orb ) then
                   ! get the atomic index
                   i = iaorb(r_tmp%r(1),lasto)
                   call rgn_range(r_tmp,lasto(i-1)+1,lasto(i))
                else
                   i = r_tmp%r(1)
                end if
                if ( IONode ) then
                   write(*,'(/,a,i0,a)')'WARNING TS: &
                        &Random atom ',i,' has been added &
                        &due to non-completeness of the connectivity &
                        &graph.'
                   write(*,'(a)')'WARNING TS: Expect sub-optimal &
                        &BTD format.'
                end if

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

    end if
    if ( .not. sort_orb ) then
       call rgn_atom2orb(r_pvt,na_u,lasto,r_tmp)
       call rgn_copy(r_tmp,r_pvt)

       tmpSp2 = tmpSp1

    end if
    call rgn_delete(r_tmp,r_Els)

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

    end do
    call rgn_delete(r_tmp)
    
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
            els,' / ',no_u_TS**2
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

  subroutine ts_tri_analyze( dit, sparse_pattern , N_Elec, Elecs, &
       ucell, na_u, lasto , nsc, isc_off , method )
    
    use parallel, only : IONode
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
    real(dp), intent(in) :: ucell(3,3)
    integer, intent(in) :: na_u, lasto(0:na_u)
    integer, intent(in) :: nsc(3), isc_off(3,product(nsc))

    ! The method used to partition the BTD format
    integer, intent(in) :: method

    type(OrbitalDistribution) :: fdit
    type(Sparsity) :: tmpSp1, tmpSp2

    integer :: no_u_TS, i, iEl, no

    integer :: n, n_nzs
    integer, pointer :: ncol(:), l_ptr(:), l_col(:)

    character(len=40), parameter :: fmt = '(/,''TS.BTD.Sort '',a,''+'',a)'
    character(len=4) :: corb
    
    ! Regions used for sorting the device region
    type(tRgn) :: r_tmp, r_Els, r_El, full
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

    do iEl = 1 , N_Elec

       i  = Elecs(iEl)%idx_o
       no = TotUsedOrbs(Elecs(iEl))

       ! we first create the super-set sparsity
       tmpSp1 = tmpSp2
       call crtSparsity_Union(fdit,tmpSp1, &
            i,i,no,no, tmpSp2)
    end do

    ! Write out all pivoting etc. analysis steps
    if ( IONode ) write(*,'(/,a)') 'transiesta: BTD analysis'

    ! Make a copy
    tmpSp1 = tmpSp2

    ! Attach the sparsity pattern
    call attach(tmpSp1, n_col = ncol, list_ptr = l_ptr, &
         list_col = l_col , nrows_g = n , nnzs = n_nzs )

    if ( IONode ) write(*,fmt) 'orb','none'
    call tri(r_pvt)

    orb_atom_switch: do orb_atom = 1 , 2
    if ( orb_atom == 1 ) then
       corb = 'orb'

       call rgn_copy(r_pvt,full)

       ! Prepare total region with all electrodes
       call rgn_range(r_Els,Elecs(1)%idx_o, &
            Elecs(1)%idx_o-1+TotUsedOrbs(Elecs(1)))
       do i = 2 , N_Elec
          call rgn_range(r_tmp,Elecs(i)%idx_o, &
               Elecs(i)%idx_o-1+TotUsedOrbs(Elecs(i)))
          call rgn_append(r_Els,r_tmp,r_Els)
       end do

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
       call rgn_range(r_Els,Elecs(1)%idx_a, &
            Elecs(1)%idx_a-1+TotUsedAtoms(Elecs(1)))
       do i = 2 , N_Elec
          call rgn_range(r_tmp,Elecs(i)%idx_a, &
               Elecs(i)%idx_a-1+TotUsedAtoms(Elecs(i)))
          call rgn_append(r_Els,r_tmp,r_Els)
       end do

    end if

    do iEl = 1 , N_Elec
       if ( IONode ) write(*,fmt) trim(corb),trim(Elecs(iEl)%name)
       if ( orb_atom == 1 ) then
          i = Elecs(iEl)%idx_o
          call rgn_range(r_El,i,i-1+TotUsedOrbs(Elecs(iEl)))
       else
          i = Elecs(iEl)%idx_a
          call rgn_range(r_El,i,i-1+TotUsedAtoms(Elecs(iEl)))
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

    if ( IONode ) write(*,fmt) trim(corb),'CM'

    call sp_pvt(tmpSp2,r_tmp, PVT_CUTHILL_MCKEE, sub = full, start = r_Els)
    if ( orb_atom == 1 ) then
       call tri(r_tmp)
    else
       call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
       call tri(r_El)
    end if

    if ( IONode ) write(*,fmt) trim(corb),'rev-CM'
    call sp_pvt(tmpSp2,r_tmp, PVT_REV_CUTHILL_MCKEE, sub = full, start = r_Els)
    if ( orb_atom == 1 ) then
       call tri(r_tmp)
    else
       call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
       call tri(r_El)
    end if
    
    if ( IONode ) write(*,fmt) trim(corb),'GPS'
    call sp_pvt(tmpSp2,r_tmp, PVT_GPS, sub = full)
    if ( orb_atom == 1 ) then
       call tri(r_tmp)
    else
       call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
       call tri(r_El)
    end if

    if ( IONode ) write(*,fmt) trim(corb),'rev-GPS'
    call sp_pvt(tmpSp2,r_tmp, PVT_REV_GPS, sub = full)
    if ( orb_atom == 1 ) then
       call tri(r_tmp)
    else
       call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
       call tri(r_El)
    end if

    if ( IONode ) write(*,fmt) trim(corb),'GGPS'
    call sp_pvt(tmpSp2,r_tmp, PVT_GGPS, sub = full)
    if ( orb_atom == 1 ) then
       call tri(r_tmp)
    else
       call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
       call tri(r_El)
    end if

    if ( IONode ) write(*,fmt) trim(corb),'rev-GGPS'
    call sp_pvt(tmpSp2,r_tmp, PVT_REV_GGPS, sub = full)
    if ( orb_atom == 1 ) then
       call tri(r_tmp)
    else
       call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
       call tri(r_El)
    end if

    end do orb_atom_switch

    call rgn_delete(r_tmp,r_Els,r_El)

    call delete(tmpSp1) ! clean up
    call delete(tmpSp2)
    call delete(fdit)

    if ( IONode ) write(*,*) ! new-line
    
  contains

    ! Print out all relevant information for this
    ! pivoting scheme
    subroutine tri(r_pvt)
      use m_pivot_methods, only : bandwidth, profile
      type(tRgn), intent(in) :: r_pvt

      integer :: bw, els, pad, work
      integer(i8b) :: prof
      type(tRgn) :: ctri


      bw   = bandwidth(n,n_nzs,ncol,l_ptr,l_col,r_pvt)
      prof = profile(n,n_nzs,ncol,l_ptr,l_col,r_pvt)
      if ( IONode ) then
         write(*,'(tr3,a,t23,i10,/,tr3,a,t13,i20)') &
              'Bandwidth: ',bw,'Profile: ',prof
      end if

      ! Create a new tri-diagonal matrix, do it in parallel
      call ts_rgn2TriMat(N_Elec, Elecs, .true., &
           fdit, tmpSp1, r_pvt, ctri%n, ctri%r, &
           method, 0, par = .true. )
      
      ! Calculate size of the tri-diagonal matrix
      els = nnzs_tri(ctri%n,ctri%r)
    
      if ( IONode ) then
         call rgn_print(ctri, name = 'BTD partitions' , &
              seq_max = 10 , indent = 3 )
         
         pad = no_u_TS ** 2
         write(*,'(tr3,a,2(i0,'' / ''),f8.6)') &
              'Matrix elements in tri / full / fraction: ', &
              els,pad,real(els,dp)/real(pad,dp)
      end if
      
      ! Get the padding for the array to hold the entire column
      call GFGGF_needed_worksize(ctri%n, ctri%r, &
           N_Elec, Elecs, pad, work)
      prof = pad + work + els * 2
      if ( IONode ) then
         write(*,'(tr3,a,t40,i10)') 'Bias padding + work elements: ', &
              pad + work
         write(*,'(tr3,a,t30,i20)') 'BTD x 2 + padding + work: ', prof
      end if

      call rgn_delete(ctri)
      
    end subroutine tri
      
  end subroutine ts_tri_analyze

end module m_ts_tri_init
