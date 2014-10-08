! Module to handle sparsity formats
! We supply routines for nullifying specific atoms
! in a sparsity format. 
! We also provide routines for removal of those entries in the equivalent
! data arrays.

! Fully created by Nick Papior Andersen, 2014

module m_sparsity_handling

  use precision, only : dp
  use class_OrbitalDistribution
  use class_Sparsity
  use class_dData1D
  use class_dSpData1D
  use class_dData2D
  use class_dSpData2D
  use geom_helper, only : iaorb, ucorb
  use m_region

  implicit none

contains

  ! Creates a sparsity pattern in the atomic subspace.
  ! I.e. we first reduce the sparsity pattern to the Gamma
  ! sparsity pattern and create the sparsity pattern for the 
  ! atoms (it will most likely be relatively dense, but...)
  subroutine SpOrb_to_SpAtom(dit,in,na_u,lasto,out)
    ! the sparsity pattern may be distributed
    type(OrbitalDistribution), intent(in) :: dit
    ! sparsity pattern to be reduced to a sparsity pattern for
    ! atomic connections
    type(Sparsity), intent(inout) :: in
    ! number of orbitals per atom
    integer, intent(in) :: na_u, lasto(0:na_u)
    ! The out-put sparsity is _never_ distributed
    ! downfolding to atoms is a huge decrease, and 
    ! it makes no sense to distribue it...
    type(Sparsity), intent(inout) :: out

    ! The new sparsity pattern
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer, allocatable :: num(:), listptr(:), list(:)
    integer :: ia, ja, io, lio, ind, indx, idx
    integer :: no_l, no_u, nnzs
    integer :: a_c(na_u), ic

    call attach(in,nrows=no_l,nrows_g=no_u, &
         n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)
    if ( no_u /= no_l ) call die('Error in conversion SpOrb2SpAt')

    allocate(num(na_u))
    a_c(:) = 0
    do ia = 1 , na_u
       ! we figure out the number of atom
       ! connections for each atom

       num(ia) = 0

       ! Initialize count
       ic = 0
       do io = lasto(ia-1) + 1 , lasto(ia)

          ! Check the local orbital
          lio = index_global_to_local(dit,io)
          if ( lio <= 0 ) cycle

          do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
          
             ! Get connecting atom
             ja = iaorb(l_col(ind),lasto)
             if ( ic == 0 ) then
                ic = ic + 1
                a_c(1) = ja
             else if ( .not. any(ja == a_c(1:ic)) ) then
                ic = ic + 1
                a_c(ic) = ja
             end if
             
          end do

       end do

       num(ia) = ic
       
    end do

    ! Create listptr
    allocate(listptr(na_u))
    listptr(1) = 0
    do ia = 2 , na_u
       listptr(ia) = listptr(ia-1) + num(ia)
    end do

    nnzs = listptr(na_u) + num(na_u)

    ! Create actual sparsity pattern 
    allocate(list(nnzs))
    indx = 1
    do ia = 1 , na_u

       ! Initialize count
       ic = 0
       do io = lasto(ia-1) + 1 , lasto(ia)

          ! Check the local orbital
          lio = index_global_to_local(dit,io)
          if ( lio <= 0 ) cycle

          do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

             ! Get connecting atom
             ja = iaorb(l_col(ind),lasto)
             if ( ic == 0 ) then
                ic = ic + 1
                a_c(1) = ja
             else if ( .not. any(ja == a_c(1:ic)) ) then
                ic = ic + 1
                a_c(ic) = ja
             end if
             
          end do
       
       end do

       num(idx:idx+ic-1) = a_c(1:ic)
       indx = indx + ic
       
    end do
    if ( nnzs /= indx ) then
       call die('Could not ensure sparsity pattern')
    end if

    ! Create new sparsity pattern and copy over
    call newSparsity(out,na_u,na_u,nnzs,num,listptr,list, &
         name='Atomic ('//name(in)//')', &
         ncols=na_u,ncols_g=na_u)

    ! Clean up
    deallocate(num,listptr,list)

  end subroutine SpOrb_to_SpAtom

  subroutine Sp_remove_crossterms(dit,in,nsc,isc_off,dir,out,r)
    ! The distribution this sparsity pattern lives in
    type(OrbitalDistribution), intent(in) :: dit
    ! sparsity pattern to be reduced
    type(Sparsity), intent(inout) :: in
    ! The supercell connections (in integers)
    integer, intent(in) :: nsc, isc_off(3,0:nsc-1)
    ! The direction in the unit-cell we wish to remove connections to
    ! All OUT OF unit-cell connections will be removed!
    integer, intent(in) :: dir
    ! The output sparsity pattern
    type(Sparsity), intent(inout) :: out
    ! The region that we restrict our limitation to
    type(tRegion), intent(in), optional :: r
    
    ! ** local variables
    ! tssp variables
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer :: no_l, no_u, nnzs
    integer, allocatable :: num(:), listptr(:), list(:)
    integer :: lio, io, ind, indx, is

    type(tRegion) :: rin

    ! set the connections to zero for the supplied atoms
    ! We do NOT reduce the sparsity row-size, ONLY the number of
    ! entries.

    call attach(in,nrows=no_l,nrows_g=no_u, &
         n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)

    if ( present(r) ) then
       call region_copy(r,rin)
    else
       call region_range(rin,1,no_u)
    end if

    allocate(num(no_l))
    do lio = 1 , no_l

       ! Initialize sparsity to 0 entries
       num(lio) = 0

       if ( l_ncol(lio) == 0 ) cycle

       io = index_local_to_global(dit,lio)
       if ( .not. in_region(rin,io) ) then
          ! we are not asked to remove these cross-terms
          num(lio) = l_ncol(lio)
          cycle
       end if

       do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

          if ( in_region(rin,ucorb(l_col(ind),no_u)) ) then

             is = (l_col(ind)-1)/no_u
             if ( isc_off(dir,is) /= 0 ) cycle

          end if
          
          ! The orbital exists on the atom
          num(lio) = num(lio) + 1
          
       end do
       
    end do

    ! Create listptr
    allocate(listptr(no_l))
    listptr(1) = 0
    do lio = 2 , no_l
       listptr(lio) = listptr(lio-1) + num(lio-1)
    end do

    nnzs = listptr(no_l) + num(no_l)

    ! Create actual sparsity pattern 
    allocate(list(nnzs))
    indx = 0
    do lio = 1 , no_l

       ! This will be zero if l_ncol(lio) is... :)
       if ( num(lio) == 0 ) cycle

       io = index_local_to_global(dit,lio)
       if ( .not. in_region(rin,io) ) then
          ! we are not asked to remove these cross-terms
          ind = l_ptr(lio)
          list(indx+1:indx+l_ncol(lio)) = l_col(ind+1:ind+l_ncol(lio))
          indx = indx + l_ncol(lio)
          cycle
       end if

       do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

          if ( in_region(rin,UCORB(l_col(ind),no_u)) ) then
             
             is = (l_col(ind)-1)/no_u
             if ( isc_off(dir,is) /= 0 ) cycle

          end if

          indx = indx + 1

          list(indx) = l_col(ind)
          
       end do
       
    end do
    if ( nnzs /= indx ) then
       call die('Could not ensure sparsity pattern')
    end if

    ! Create new sparsity pattern and copy over
    call newSparsity(out,no_l,no_u,nnzs,num,listptr,list, &
         name='Truncated '//trim(name(in)), &
         ncols=ncols(in),ncols_g=ncols_g(in))

    call region_delete(rin)
    
    ! Clean up
    deallocate(num,listptr,list)

  end subroutine Sp_remove_crossterms
  
  subroutine Sp_remove_region(dit,in,rr,out)
    ! The distribution this sparsity pattern lives in
    type(OrbitalDistribution), intent(in) :: dit
    ! sparsity pattern to be reduced
    type(Sparsity), intent(inout) :: in
    ! The region we wish to remove
    type(tRegion), intent(in) :: rr
    ! The region of orbitals that will be removed
    type(Sparsity), intent(inout) :: out
    
    ! ** local variables
    ! tssp variables
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer :: no_l, no_u, nnzs
    integer, allocatable :: num(:), listptr(:), list(:)
    integer :: lio, io, ind, jo, indx

    ! set the connections to zero for the supplied atoms
    ! We do NOT reduce the sparsity row-size, ONLY the number of
    ! entries.

    call attach(in,nrows=no_l,nrows_g=no_u, &
         n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)

    allocate(num(no_l))
    do lio = 1 , no_l

       ! Initialize sparsity to 0 entries
       num(lio) = 0

       if ( l_ncol(lio) == 0 ) cycle

       io = index_local_to_global(dit,lio)

       if ( in_region(rr,io) ) cycle

       do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
          
          jo = ucorb(l_col(ind),no_u)
          if ( in_region(rr,jo) ) cycle
       
          ! The orbital exists on the atom
          num(lio) = num(lio) + 1
          
       end do
       
    end do

    ! Create listptr
    allocate(listptr(no_l))
    listptr(1) = 0
    do lio = 2 , no_l
       listptr(lio) = listptr(lio-1) + num(lio-1)
    end do

    nnzs = listptr(no_l) + num(no_l)

    ! Create actual sparsity pattern 
    allocate(list(nnzs))
    indx = 0
    do lio = 1 , no_l

       io = index_local_to_global(dit,lio)
       if ( in_region(rr,io) ) cycle

       do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

          jo = ucorb(l_col(ind),no_u)
          if ( in_region(rr,jo) ) cycle

          indx = indx + 1

          list(indx) = l_col(ind)
          
       end do
       
    end do
    if ( nnzs /= indx ) then
       call die('Could not ensure sparsity pattern')
    end if

    ! Create new sparsity pattern and copy over
    call newSparsity(out,no_l,no_u,nnzs,num,listptr,list, &
         name='Truncated '//trim(name(in)), &
         ncols=ncols(in),ncols_g=ncols_g(in))

    ! Clean up
    deallocate(num,listptr,list)
    
  end subroutine Sp_remove_region

  subroutine Sp_retain_region(dit,in,rr,out)
    ! The distribution this sparsity pattern lives in
    type(OrbitalDistribution), intent(in) :: dit
    ! sparsity pattern to be reduced
    type(Sparsity), intent(inout) :: in
    ! The region to retain
    type(tRegion), intent(in) :: rr
    ! The region of orbitals that will be removed
    type(Sparsity), intent(inout) :: out
    
    type(tRegion) :: full, rem_r
    integer :: no_u

    call attach(in,nrows_g=no_u)

    ! instead of creating two different
    ! algorithms we create the complement region
    ! and then remove that!
    ! HOW genius! :)

    call region_range(full,1,no_u)
    
    ! create complement of rem_r
    call region_complement(full,rr,rem_r)
    ! clean-up
    call region_delete(full)

    call Sp_remove_region(dit,in,rem_r,out)

    call region_delete(rem_r)
    
  end subroutine Sp_retain_region

  subroutine Sp_remove_region2region(dit,in,r1,r2,out)
    ! The distribution this sparsity pattern lives in
    type(OrbitalDistribution), intent(in) :: dit
    ! sparsity pattern to be reduced
    type(Sparsity), intent(inout) :: in
    ! The cross-terms "from" and "to", or "to" and "from"
    type(tRegion), intent(in) :: r1, r2
    ! The region of orbitals that will be removed
    type(Sparsity), intent(inout) :: out
    
    ! ** local variables
    ! tssp variables
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer :: no_l, no_u, nnzs
    integer, allocatable :: num(:), listptr(:), list(:)
    integer :: lio, io, ind, jo, indx, ridx

    ! set the connections to zero for the supplied atoms
    ! We do NOT reduce the sparsity row-size, ONLY the number of
    ! entries.

    call attach(in,nrows=no_l,nrows_g=no_u, &
         n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)

    allocate(num(no_l))
    do lio = 1 , no_l

       ! Initialize sparsity to 0 entries
       num(lio) = 0

       if ( l_ncol(lio) == 0 ) cycle

       io = index_local_to_global(dit,lio)

       if ( in_region(r1,io) ) then
          ridx = 1
       else if ( in_region(r2,io) ) then
          ridx = 2
       else
          num(lio) = l_ncol(lio)
          cycle
       end if

       do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
          
          if ( ridx == 1 ) then
             jo = ucorb(l_col(ind),no_u)
             ! the i'th orbital is in region 1
             ! now if jo is in region 2 we have a match
             ! and remove that orbital connection
             if ( in_region(r2,jo) ) cycle
          else if ( ridx == 2 ) then
             jo = ucorb(l_col(ind),no_u)
             if ( in_region(r1,jo) ) cycle
          end if
       
          ! The orbital exists on the atom
          num(lio) = num(lio) + 1
          
       end do
       
    end do

    ! Create listptr
    allocate(listptr(no_l))
    listptr(1) = 0
    do lio = 2 , no_l
       listptr(lio) = listptr(lio-1) + num(lio-1)
    end do

    nnzs = listptr(no_l) + num(no_l)

    ! Create actual sparsity pattern 
    allocate(list(nnzs))
    indx = 0
    do lio = 1 , no_l

       io = index_local_to_global(dit,lio)
       if ( in_region(r1,io) ) then
          ridx = 1
       else if ( in_region(r2,io) ) then
          ridx = 2
       else
          ridx = 0
       end if

       do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

          if ( ridx == 1 ) then
             jo = ucorb(l_col(ind),no_u)
             ! the i'th orbital is in region 1
             ! now if jo is in region 2 we have a match
             ! and remove that orbital connection
             if ( in_region(r2,jo) ) cycle
          else if ( ridx == 2 ) then
             jo = ucorb(l_col(ind),no_u)
             if ( in_region(r1,jo) ) cycle
          end if

          indx = indx + 1
          list(indx) = l_col(ind)
          
       end do
       
    end do
    if ( nnzs /= indx ) then
       call die('Could not ensure sparsity pattern')
    end if

    ! Create new sparsity pattern and copy over
    call newSparsity(out,no_l,no_u,nnzs,num,listptr,list, &
         name='Truncated '//trim(name(in)), &
         ncols=ncols(in),ncols_g=ncols_g(in))

    ! Clean up
    deallocate(num,listptr,list)
    
  end subroutine Sp_remove_region2region

  subroutine dSpData1D_to_Sp(A,out,B)

    ! the sparse data that needs to be transferred
    ! to the new sparsity format (out)
    type(dSpData1D), intent(inout) :: A,B
    ! The sparsity it needs to be transferred to
    type(Sparsity), intent(inout) :: out

    type(dData1D) :: dat
    real(dp), pointer :: inA(:), outA(:)

    ! get lists of sparsity patterns
    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: sp
    integer, pointer :: l_ncol(:) , l_ptr(:) , l_col(:)
    integer, pointer :: n_ncol(:), n_ptr(:), n_col(:)
    integer :: lio, io, nr, lnr, lind, rind, nind

    ! get value and distribution
    inA => val(A)
    dit => dist(A)
    sp  => spar(A)

    ! Allocate the new data array
    io = nnzs(out)
    call newdData1D(dat,io,trim(name(A))//' reduced')
    outA => val(dat)

    call attach(sp ,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=lnr,nrows_g=nr)
    call attach(out,n_col=n_ncol,list_ptr=n_ptr,list_col=n_col, &
         nrows=lio,nrows_g=io)
    if ( lnr /= lio .or. nr /= io ) then
       call die('Error in reduced sparsity pattern')
    end if
    
    nind = 0
    do lio = 1 , lnr
       
       if ( l_ncol(lio) == 0 ) cycle
       if ( n_ncol(lio) == 0 ) cycle
       
       do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
          
          find: do rind = n_ptr(lio) + 1 , n_ptr(lio) + n_ncol(lio)
             if ( l_col(lind) == n_col(rind) ) then
                outA(rind) = inA(lind)
                nind = nind + 1
                exit find
             end if
          end do find

       end do

       ! Check that copying goes as planned
       if ( nind /= n_ptr(lio) + n_ncol(lio) ) then
          call die('Error in sparsity copying')
       end if
       
    end do

    ! copy over data
    call newdSpData1D(out,dat,dit,B,trim(name(A))//' reduced')

    ! Delete the data (we have copied it...)
    call delete(dat)
    
  end subroutine dSpData1D_to_Sp

  subroutine dSpData2D_to_Sp(A,out,B)

    ! the sparse data that needs to be transferred
    ! to the new sparsity format (out)
    type(dSpData2D), intent(inout) :: A, B
    ! The sparsity it needs to be transferred to
    type(Sparsity), intent(inout) :: out

    type(dData2D) :: dat
    real(dp), pointer :: inA(:,:), outA(:,:)

    ! get lists of sparsity patterns
    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: sp
    integer, pointer :: l_ncol(:) , l_ptr(:) , l_col(:)
    integer, pointer :: n_ncol(:), n_ptr(:), n_col(:)
    integer :: lio, io, nr, lnr, lind, rind, nind
    integer :: dim_sp
    integer :: d

    ! get value and distribution
    inA => val(A)
    dit => dist(A)
    sp  => spar(A)

    ! Create the data array for the reduced size
    dim_sp = spar_dim(A)
    io = nnzs(out)
    if ( dim_sp == 1 ) then
       d = size(inA,dim=2)
       call newdData2D(dat,io,d,trim(name(A))//' reduced')
    else
       d = size(inA,dim=1)
       call newdData2D(dat,d,io,trim(name(A))//' reduced')
    end if
    outA => val(dat)

    call attach(sp ,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=lnr,nrows_g=nr)
    call attach(out,n_col=n_ncol,list_ptr=n_ptr,list_col=n_col, &
         nrows=lio,nrows_g=io)
    if ( lnr /= lio .or. nr /= io ) then
       call die('Error in reduced sparsity pattern')
    end if
    
    nind = 0
    do lio = 1 , lnr
       
       if ( l_ncol(lio) == 0 ) cycle
       if ( n_ncol(lio) == 0 ) cycle
       
       select case ( d )
       case ( 1 )
          do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
             
             find1: do rind = n_ptr(lio) + 1 , n_ptr(lio) + n_ncol(lio)
                if ( l_col(lind) == n_col(rind) ) then
                   outA(rind,:) = inA(lind,:)
                   nind = nind + 1
                   exit find1
                end if
             end do find1
          end do
       case ( 2 )
          do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
             
             find2: do rind = n_ptr(lio) + 1 , n_ptr(lio) + n_ncol(lio)
                if ( l_col(lind) == n_col(rind) ) then
                   outA(:,rind) = inA(:,lind)
                   nind = nind + 1
                   exit find2
                end if
             end do find2
          end do
       end select

       ! Check that copying goes as planned
       if ( nind /= n_ptr(lio) + n_ncol(lio) ) then
          call die('Error in sparsity copying')
       end if
       
    end do

    ! copy over data
    call newdSpData2D(out,dat,dit,B,trim(name(A))//' reduced')

    ! Delete the data (we have copied it...)
    call delete(dat)
    
  end subroutine dSpData2D_to_Sp

  subroutine dSpData2D_interp(N,A_2Ds,x,x0)
    use m_interpolate
    integer, intent(in) :: N
    type(dSpData2D), intent(inout) :: A_2Ds(N)
    ! The real values that we wish to interpolate from
    real(dp), intent(in) :: x(N), x0 ! and to

    ! Container to do array assignments
    type :: A2D
       real(dp), pointer :: v(:,:)
    end type A2D
    type(A2D) :: array(N)

    integer :: io, i, j, n_nzs, sp_dim, dim2
    real(dp) :: y(N), y0


    ! prep the array assignments...
    do i = 1 , N
       array(i)%v => val(A_2Ds(i))
    end do

    ! Get the dimensionality and the sparsity dimension
    sp_dim = spar_dim(A_2Ds(1))
    if ( sp_dim == 1 ) then
       dim2 = size(array(1)%v(1,:))
    else
       call die('How on earth should we interpolate &
            dependent data...')
    end if
    n_nzs = nnzs(A_2Ds(1))

    ! figure out the sparsity dimension
    do j = 1 , dim2
       do io = 1 , n_nzs
          
          ! Copy over 
          do i = 1 , N
             y(i) = array(i)%v(io,j)
          end do

          call interp_spline(N,x,y,x0,y0)
             
          ! Copy over...
          array(1)%v(io,j) = y0
          
       end do
    end do
    
  end subroutine dSpData2D_interp

end module m_sparsity_handling
