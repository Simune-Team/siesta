
! This module provides handles for holding regions of orbitals/atoms.
! It is a convenient in-place type to easily handle
! different kinds of regions of any kind.
! It easily handles connections with regions to
! generate limited sets of regions (or expand them).

! Fully created by Nick Papior Andersen, 2014

module m_region
  
  use intrinsic_missing, only : uniqc
  use geom_helper, only : ucorb, iaorb
  use class_OrbitalDistribution
  use class_Sparsity

  implicit none

  private

  ! max length of region
  integer, parameter, public :: R_NAME_LEN = 50

  type :: tRegion
     ! The name of the region
     character(len=R_NAME_LEN) :: name = ' '
     ! The quantities in the region
     integer :: n = 0
     integer, pointer :: r(:) => null()
  end type tRegion

  public :: tRegion
  public :: region_delete
  public :: region_connect
  public :: region_intersection
  public :: region_union
  public :: region_complement
  public :: region_range
  public :: region_list
  public :: region_correct_atom
  public :: region_atom2orb
  public :: region_orb2atom
  public :: region_overlaps
  public :: region_sort
  public :: region_insert
  public :: region_print
  public :: region_copy
  public :: in_region, region_pivot
#ifdef MPI
  public :: region_MPI_union
  public :: region_MPI_Bcast
#endif

  ! Sorting methods for the regions
  integer, parameter, public :: R_SORT_MAX_FRONT = 1
  integer, parameter, public :: R_SORT_MAX_BACK = 2

  interface region_remove
     module procedure region_remove_region
     module procedure region_remove_list
  end interface region_remove

  interface region_insert
     module procedure region_insert_region
     module procedure region_insert_list
  end interface region_insert

contains

  ! Get the index of the element that equal 'i'
  ! We refer to this as the pivoting element.
  ! In general this is the same as doing:
  !   A(r%r == i)
  elemental function region_pivot(r,i) result(idx)
    type(tRegion), intent(in) :: r
    integer, intent(in) :: i
    integer :: idx
    do idx = 1 , r%n
       if ( r%r(idx) == i ) return
    end do
    idx = 0
  end function region_pivot

  ! Check whether an element exists in this region.
  ! Performs an easy check instead of doing it manually.
  elemental function in_region(r,i) result(in)
    type(tRegion), intent(in) :: r
    integer, intent(in) :: i
    logical :: in
    in = any(i == r%r)
  end function in_region

  ! Copies one region to another one
  ! This emplies that the element 'to' is deleted before
  ! copy is performed.
  subroutine region_copy(from,to)
    type(tRegion), intent(in) :: from
    type(tRegion), intent(inout) :: to
    call region_list(to,from%n,from%r,name=from%name)
  end subroutine region_copy

  ! Fully deletes a region (irrespective of it's current state)
  subroutine region_delete(r)
    type(tRegion), intent(inout) :: r
    r%name = ' '
    r%n = 0
    if ( associated(r%r) ) then
       deallocate(r%r)
    end if
    nullify(r%r)
  end subroutine region_delete

  ! Allows removing certain elements from a region.
  subroutine region_remove_list(r,n,list,rout)
    type(tRegion), intent(in) :: r
    type(tRegion), intent(inout) :: rout
    integer, intent(in) :: n, list(n)

    integer :: i, ni
    integer, allocatable :: tmp_r(:)

    if ( n == 0 ) then
       call region_copy(r,rout)
       return
    end if

    call region_delete(rout)
    
    allocate(tmp_r(r%n))
    ni = 0
    do i = 1 , r%n
       if ( all(r%r(i)/=list(:)) ) then
          ni = ni + 1
          tmp_r(ni) = r%r(i)
       end if
    end do

    call region_list(rout,ni,tmp_r(1:ni))

    deallocate(tmp_r)

  end subroutine region_remove_list
  subroutine region_remove_region(r,rr,rout)
    type(tRegion), intent(in) :: r
    type(tRegion), intent(in) :: rr
    type(tRegion), intent(inout) :: rout
    
    call region_remove_list(r,rr%n,rr%r,rout)

  end subroutine region_remove_region


  ! Generates a new region which connects to 'r'
  ! We have several options to control the output region
  !   1. Only do a "first touch" thereby only gathering the 
  !      elements that directly connects with 'r'
  !   2. 'follow' will iteratively follow the sparsity pattern
  !      and create a region which encompass everything
  !   3. The user can request that a certain region be not crossed
  !      I.e. if 'except' is provided any connections to that region
  !      will be disregarded.
  !      This is usefull if you want to create a limited sparsity pattern
  !      which connects from 1 orbital but ends before some orbital <idx>
  !   4. The user can request a region of orbitals which are in charge of 
  !      the connections, one could easily imagine 2 orbitals where
  !      one is connecting out, the other does not, in that case would
  !      'connect_from' only contain one orbital.
  ! NOTE: It DOES work in parallel
  subroutine region_connect(r,dit,sp,cr,except, connect_from, follow)

    ! the region we wish to find the connections to
    type(tRegion), intent(in) :: r
    ! The sparsity patterns distribution
    type(OrbitalDistribution), intent(in) :: dit
    ! The sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! the connecting region
    type(tRegion), intent(inout) :: cr
    ! a region which must not be crossed (exception region)
    type(tRegion), intent(in), optional :: except
    ! a region which contains the orbitals that are 
    ! responsible for the connection
    type(tRegion), intent(inout), optional :: connect_from
    ! Tell whether we should follow the connection
    ! or stop at the first iteration
    logical, intent(in), optional :: follow

    ! ** local variables
#ifdef MPI
    type(tRegion) :: tmp
#endif
    integer :: i, j, io, jo, ind, no_l, no_u, it, rt
    integer, allocatable :: ct(:), rr(:)
    integer, pointer :: err(:)
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)

    call attach(sp,nrows=no_l,nrows_g=no_u, &
         n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)

    ! ensure nullification of cr
    call region_delete(cr)

    allocate(rr(r%n))
    allocate(ct(no_u-r%n))

    ! Default the err region to be the r-region
    err => r%r
    if ( present(except) ) err => except%r

    rt = 0
    it = 0
    do i = 1 , r%n

       ! Orbital that should be folded from
       io = index_global_to_local(dit,r%r(i))
       if ( io <= 0 ) cycle ! the orbital does not exist on this node

       ! Count number of available orbitals to fold to
       j = 0
       do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)

          ! UC-orb
          jo = ucorb(l_col(ind),no_u)

          ! Ensure that it is not a folding to the same region
          if ( any(jo == r%r) ) cycle
          if ( any(jo == err) ) cycle ! in case er has been provided

          if ( it == 0 ) then
             it = it + 1
             ct(it) = jo
          else if ( .not. any(jo == ct(1:it)) ) then
             it = it + 1
             ct(it) = jo
          end if

          ! The connect_from region
          ! will only make sense if we do not follow the orbitals
          ! (else we still get where the connection starts)!
          if ( rt == 0 ) then
             rt = rt + 1
             rr(rt) = io
          else if ( .not. any(io == rr(1:rt)) ) then
             rt = rt + 1
             rr(rt) = io
          end if

       end do
    end do

#ifdef MPI
    ! We have found the first group of regions
    ! Now we need to align the regions for all the processors
    ! in this orbital distribution (before we move on 
    ! and check for "follow"

    ! First we do the ct region
    call region_list(tmp,it,ct)
    call region_MPI_union(dit,tmp)
    it = tmp%n
    if ( it > 0 ) then
       ct(1:it) = tmp%r
    end if
    call region_list(tmp,rt,rr)
    call region_MPI_union(dit,tmp)
    rt = tmp%n
    if ( rt > 0 ) then
       rr(1:rt) = tmp%r
    end if
    call region_delete(tmp) ! clean-up

#endif

    ! If we are supposed to follow
    ! then loop and extend ct
    if ( present(follow) ) then
    if ( follow ) then

       i = 1
       do while ( i <= it )

          ! Orbital that should be folded from
          io = ct(i)

          ! Count number of available orbitals to fold to
          j = 0
          do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)

             ! UC-orb
             jo = ucorb(l_col(ind),no_u)

             ! Ensure that it is not a folding to the same region
             if ( any(jo == r%r) ) cycle
             if ( any(jo == err) ) cycle ! in case er has been provided

             if ( it == 0 ) then
                it = it + 1
                ct(it) = jo
             else if ( .not. any(jo == ct(1:it)) ) then
                it = it + 1
                ct(it) = jo
             end if
             
          end do
          
          ! step the orbital that we now follow
          i = i + 1

       end do

    end if
    end if

    ! We now have a list of orbitals that needs to be folded to
    ! Copy the list over
    call region_list(cr,it,ct(1:it))

    deallocate(ct)

#ifdef MPI
    call region_MPI_union(dit,cr)
#endif

    if ( present(connect_from) ) then

       call region_list(connect_from,rt,rr(1:rt))

#ifdef MPI
       call region_MPI_union(dit,connect_from)
#endif

    end if

    deallocate(rr)

  end subroutine region_connect


  ! This routine WORKS IN PARALLEL
  ! Will sort a region of orbitals based on the sparsity
  ! pattern. In effect this can be used as a pivoting method for
  ! obtaing the smallest bandwidth of a matrix.
  ! It currently implements a method which
  !  -> takes a region 'r' and sort the elements 'sr' in that region
  !     All orbitals in sr MUST exist in 'r', it checks and dies if this
  !     is not fulfilled.
  !     We then allow for two different sorting methods:
  !  1. Sort according to the maximum connections in front of the region 'sr'
  !     I.e. we find number of connection orbitals for all in 'sr' which does
  !     does not connect into region 'r'.
  !     We then place the most connecting orbital in 'sr' last in region 'r'
  !  2. Sort according to the maximum connections back to the region 'r'
  !     I.e. we find the number of connection orbitals for all in 'sr' which 
  !     connects into region 'r'.
  !     We then place 
  subroutine region_sort(r,dit,sp, sr, method)

#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD, MPI_Integer
    use mpi_siesta, only : MPI_MAX, MPI_AllReduce
#endif

    ! the region we wish to find the connections to
    type(tRegion), intent(inout) :: r
    ! The sparsity patterns distribution
    type(OrbitalDistribution), intent(in) :: dit
    ! The sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! the sorting region (i.e. the orbitals that are allowed
    ! to be pivoted)
    type(tRegion), intent(inout) :: sr
    ! The method used for sorting
    integer, intent(in) :: method

    ! ** local variables
    type(tRegion) :: r_tmp
    integer :: no_u, i, io, ci, ind, jo, last_n, idx_max
    integer, allocatable :: n_c(:), cur_con(:)
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
#ifdef MPI
    integer :: comm
    integer :: MPIerror
#endif

    if ( r%n == 0 ) return

#ifdef MPI
    comm = dist_comm(dit)
#endif

    ! Attach to the sparsity pattern...
    call attach(sp,nrows_g=no_u, &
         n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)

    ! First we ensure that the lists are the same
    call region_intersection(r,sr,r_tmp)
    if ( r_tmp%n /= sr%n ) then
       call die('The regions are not well-defined')
    end if
    call region_delete(r_tmp)

    ! Prepare the collection arrays
    allocate(n_c(sr%n))
    n_c = 0 ! initialize to ensure MPI reduction
    allocate(cur_con(no_u))

    select case ( method ) 

    case ( R_SORT_MAX_FRONT ) 
       ! we sort and move the orbitals
       ! down in the rank according to their maximum
       ! range.
       ! A higher number of connections outside of the region
       ! will move it further down.

       ! Step 1. find the number of connections for
       ! each orbital in the sorting region
       
       do i = 1 , sr%n

          ! Orbital that should be folded from
          io = index_global_to_local(dit,sr%r(i))
          if ( io <= 0 ) cycle ! the orbital does not exist on this node
          
          ci = 0
          do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
             
             ! UC-orb
             jo = ucorb(l_col(ind),no_u)
             
             ! Ensure that it is not folding to the same region
             if ( in_region(r,jo) ) cycle

             if ( ci == 0 ) then
                ci = ci + 1
                cur_con(ci) = jo
             else if ( .not. any(jo == cur_con(1:ci)) ) then
                ci = ci + 1
                cur_con(ci) = jo
             end if
          end do

          ! As the region is created "from a frontal search"
          ! we will reverse the result (maxloc returns the first max)
          n_c(sr%n-i+1) = ci

       end do

#ifdef MPI
       if ( dist_nodes(dit) > 1 ) then
          call MPI_AllReduce(n_c,cur_con,sr%n,MPI_Integer, &
               MPI_MAX, comm, MPIerror)
          n_c(1:sr%n) = cur_con(1:sr%n)
       end if
#endif

    case ( R_SORT_MAX_BACK ) 

       ! we sort and move the orbitals
       ! back in the rank according to their maximum
       ! range backwards.
       ! A higher number of connections inside of the region
       ! will move it further back.

       ! Step 1. find the number of connections for
       ! each orbital in the sorting region
       
       do i = 1 , sr%n

          ! Orbital that should be folded from
          io = index_global_to_local(dit,sr%r(i))
          if ( io <= 0 ) cycle ! the orbital does not exist on this node

          if ( l_ncol(io) == 0 ) cycle

          ci = 0
          do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
             
             ! UC-orb
             jo = ucorb(l_col(ind),no_u)
             
             ! Ensure that it is folding to the same region
             if ( .not. in_region(r,jo) ) cycle

             if ( ci == 0 ) then
                ci = ci + 1
                cur_con(ci) = jo
             else if ( .not. any(jo == cur_con(1:ci)) ) then
                ci = ci + 1
                cur_con(ci) = jo
             end if
          end do

          ! As the region is created "from a frontal search"
          ! we will reverse the result (maxloc returns the first max)
          n_c(sr%n-i+1) = ci

       end do

#ifdef MPI
       if ( dist_nodes(dit) > 1 ) then
          call MPI_AllReduce(n_c,cur_con,sr%n,MPI_Integer, &
               MPI_MAX, comm, MPIerror)
          n_c(1:sr%n) = cur_con(1:sr%n)
       end if
#endif

       ! We now reverse the values as we then control that the one
       ! with the least connections back, will be placed last
       n_c(:) = maxval(n_c,dim=1) - n_c(:)
                 
    case default

       call die('not implemented')

    end select

    select case ( method ) 

    case ( R_SORT_MAX_FRONT, R_SORT_MAX_BACK )

       ! We now have how many orbitals they all connect to
       ! Sort them accordingly!
       last_n = r%n ! where to position the next orbital

       do jo = 1 , sr%n

          idx_max = maxloc(n_c,dim=1)
          
          ! The index is reversed, reverse back and retrieve
          ! the equivalent orbital
          io = sr%r(sr%n - idx_max + 1)
          do i = 1 , last_n
             if ( r%r(i) == io ) then
                ! swap the positions
                r%r(i) = r%r(last_n)
                r%r(last_n) = io
                last_n = last_n - 1
                exit ! we are not going to find it later...
             end if
          end do

          n_c(idx_max) = -1 ! reset connections, we already
                            ! processed the orbital

       end do

       ! We have now sorted the entries according to the maximum
       ! connections across the existing region.

    end select

    deallocate(n_c,cur_con)

  end subroutine region_sort

  ! Creates a region in terms of the INTERSECTION of two regions 
  subroutine region_intersection(r1,r2,ir)
    ! the regions we wish to operate on
    type(tRegion), intent(in) :: r1, r2
    ! the intersection region
    type(tRegion), intent(inout) :: ir

    ! ** local variables
    integer :: i, it
    integer, allocatable :: ct(:)

    if ( r1%n == 0 .or. r2%n == 0 ) then
       call region_delete(ir)
       return
    end if

    it = 0
    allocate(ct(min(r1%n,r2%n)))

    do i = 1 , r1%n
       
       if ( .not. any(r1%r(i) == r2%r) ) cycle
       it = it + 1
       ct(it) = r1%r(i)

    end do

    ! We now have a list of orbitals that needs to be folded to
    ! Copy the list over
    call region_list(ir,it,ct(1:it))

    deallocate(ct)
    
  end subroutine region_intersection

  ! Creates a region in terms of the UNION of two regions 
  subroutine region_union(r1,r2,ur)
    ! the regions we wish to operate on
    type(tRegion), intent(in) :: r1, r2
    ! the unified region
    type(tRegion), intent(inout) :: ur

    ! ** local variables
    integer :: i, it
    integer, allocatable :: ct(:)

    if ( r1%n == 0 ) then
       call region_copy(r2,ur)
       return
    else if ( r2%n == 0 ) then
       call region_copy(r1,ur)
       return
    end if

    allocate(ct(r1%n+r2%n))

    ! Copy over r1
    ct(1:r1%n) = r1%r(:)

    it = r1%n
    do i = 1 , r2%n
       
       if ( any(r2%r(i) == ct(1:it)) ) cycle

       it = it + 1
       ct(it) = r2%r(i)

    end do

    ! We now have a list of orbitals that needs to be folded to
    ! Copy the list over
    call region_list(ur,it,ct(1:it))

    deallocate(ct)

  end subroutine region_union

  subroutine region_insert_list(r,n,list,rout,idx)
    type(tRegion), intent(in) :: r
    integer, intent(in) :: n, list(n)
    type(tRegion), intent(inout) :: rout
    ! The place of insertion
    integer, intent(in) :: idx
    type(tRegion) :: rtmp
    call region_list(rtmp,n,list)
    call region_insert_region(r,rtmp,rout,idx)
    call region_delete(rtmp)
  end subroutine region_insert_list

  ! Inserts a region in another region.
  subroutine region_insert_region(r,rin,rout,idx)
    ! The region we wish to insert in
    type(tRegion), intent(in) :: r
    type(tRegion), intent(in) :: rin
    type(tRegion), intent(inout) :: rout
    ! The place of insertion
    !  [idx] at index in r%r ( after idx )
    ! hence, idx == 0 will be in the beginning.
    integer, intent(in) :: idx

    ! Local index
    integer :: lidx, i,j 
    integer, allocatable :: itmp(:)

    if ( idx > r%n ) then
       call die('Error, index does not exist')
    else if ( idx == r%n ) then
       call region_union(r,rin,rout)
       return
    else if ( idx <= 0 ) then
       call region_union(rin,r,rout)
       return
    end if

    ! We are placing it somewhere in the middle

    ! Insert list...
    ! This is not soo easy.
    ! First we need to find the first value before
    ! the list (in case it already exists)
    lidx = idx
    do while ( in_region(rin,r%r(lidx)) )
       lidx = lidx - 1
       if ( lidx == 0 ) then
          call region_union(rin,r,rout)
          return
       end if
    end do

    ! Grab value before removing potential duplicates
    i = r%r(lidx)
    call region_remove(r,rin,rout)
    ! Back-retrieve the new position
    i = region_pivot(rout,i)
    if ( i == 0 ) call die('Error in algorithm')
    allocate(itmp(rout%n+rin%n))
    itmp(1:i) = rout%r(1:i)
    itmp(i+1:i+rin%n) = rin%r(:)
    j = i + rin%n + 1
    itmp(j:) = rout%r(i+1:)
    call region_list(rout,size(itmp),itmp)
    deallocate(itmp)

  end subroutine region_insert_region

  ! Creates a region in terms of the COMPLEMENT of two regions 
  subroutine region_complement(r1,r2,cr)
    ! the regions we wish to operate on
    type(tRegion), intent(in) :: r1, r2
    ! the complementary region
    type(tRegion), intent(inout) :: cr

    ! ** local variables
    integer :: i, it
    integer, allocatable :: ct(:)

    allocate(ct(r1%n))

    it = 0
    do i = 1 , r1%n
       
       if ( any(r1%r(i) == r2%r) ) cycle

       it = it + 1
       ct(it) = r1%r(i)

    end do

    ! Copy the list over
    call region_list(cr,it,ct(1:it))

    deallocate(ct)

  end subroutine region_complement

  ! Easy setup of a consecutive orbital range
  subroutine region_range(r,o1,o2)
    ! The region containing the range [o1;o2]
    type(tRegion), intent(inout) :: r
    ! The limits on the range
    integer, intent(in) :: o1, o2
    
    integer :: io
    call region_delete(r)

    r%n = abs(o2 - o1) + 1
    allocate(r%r(r%n))
    if ( o1 <= o2 ) then
       do io = o1 , o2
          r%r(io-o1+1) = io
       end do
    else
       do io = o2 , o1
          r%r(io-o2+1) = io
       end do
    end if

  end subroutine region_range

  subroutine region_list(r,n,list,name)
    ! Region to put list in
    type(tRegion), intent(inout) :: r
    ! list to copy over
    integer, intent(in) :: n, list(n)
    character(len=*), intent(in), optional :: name

    call region_delete(r)
    r%n = n
    if ( n > 0 ) then
       allocate(r%r(n))
       r%r(1:n) = list(1:n)
    end if
    if ( present(name) ) r%name = name
    
  end subroutine region_list

  subroutine region_correct_atom(r,na_u,lasto)
    ! Region which we want to extend with the atoms
    ! that it already overlaps
    type(tRegion), intent(inout) :: r
    ! The last orbitals of each atom
    integer, intent(in) :: na_u, lasto(0:na_u)

    ! ** local variables
    integer :: a(na_u), i, j, na, no
    character(len=R_NAME_LEN) :: tmp

    na = 0
    a = 0
    do i = 1 , r%n
       
       ! Get the atom index of the orbital
       j = iaorb(r%r(i),lasto)

       if ( na == 0 ) then
          na = na + 1
          a(na) = j
       else if ( .not. any( j == a(1:na) ) ) then
          na = na + 1
          a(na) = j
       end if

    end do

    ! Count orbital size
    no = 0
    do i = 1 , na
       no = no + lasto(a(i)) - lasto(a(i)-1)
    end do

    ! to be sure we just depopulate the region
    ! and populate it with the correct orbitals
    tmp = r%name
    call region_delete(r)
    r%name = tmp
    r%n = no
    nullify(r%r)
    if ( r%n > 0 ) then
       allocate(r%r(no))
       no = 0
       do i = 1 , na
          do j = lasto(a(i)-1) + 1 , lasto(a(i))
             no = no + 1
             r%r(no) = j
          end do
       end do
    end if

  end subroutine region_correct_atom

  subroutine region_Atom2Orb(ar,na_u,lasto,or)
    ! Region which we want to extend with the atoms
    ! that it already overlaps
    type(tRegion), intent(in) :: ar
    ! The last orbitals of each atom
    integer, intent(in) :: na_u, lasto(0:na_u)
    type(tRegion), intent(inout) :: or

    ! ** local variables
    integer :: i, cr, no
    character(len=R_NAME_LEN) :: tmp

    ! First calculate size of region in orbital
    ! space
    no = 0
    do i = 1 , ar%n
       no = no + lasto(ar%r(i))-lasto(ar%r(i)-1)
    end do

    ! to be sure we just depopulate the region
    ! and populate it with the correct orbitals
    tmp = or%name
    call region_delete(or)
    or%name = tmp
    or%n = no
    nullify(or%r)
    if ( no > 0 ) then
       allocate(or%r(no))
       no = 0
       do i = 1 , ar%n
          do cr = lasto(ar%r(i)-1) + 1 , lasto(ar%r(i))
             no = no + 1
             or%r(no) = cr
          end do
       end do
    end if
    
  end subroutine region_Atom2Orb

  subroutine region_Orb2Atom(or,na_u,lasto,ar)
    ! The orbital region we wish to convert to an atomic
    ! region
    type(tRegion), intent(in) :: or
    ! The last orbitals of each atom
    integer, intent(in) :: na_u, lasto(0:na_u)
    type(tRegion), intent(inout) :: ar

    ! ** local variables
    integer :: io, ia, na, a
    integer, allocatable :: ca(:)
    character(len=R_NAME_LEN) :: tmp

    ! The maximum number of atoms in the orbital region
    na = or%n / minval(lasto(1:na_u)-lasto(0:na_u-1)) + 1
    allocate(ca(na))
    ia = 1
    ca(1) = iaorb(or%r(1),lasto)
    do io = 2 , or%n
       a = iaorb(or%r(io),lasto)
       if ( any(a == ca(1:ia)) ) cycle
       ia = ia + 1
       if ( ia > na ) call die('Error in program')
       ca(ia) = a
    end do
    na = ia

    ! to be sure we just depopulate the region
    ! and populate it with the correct atoms
    tmp = ar%name
    call region_list(ar,na,ca(1:na), name = tmp )

  end subroutine region_Orb2Atom

  ! Creates the folded to orbitals by considering
  function region_overlaps(r1,r2) result(overlap)
    ! the regions we wish to find the union of
    type(tRegion), intent(in) :: r1, r2
    logical :: overlap

    ! ** local variables
    integer :: i

    overlap = .false.

    if ( r1%n == 0 ) return
    if ( r2%n == 0 ) return

    do i = 1 , r1%n
       
       if ( any(r1%r(i) == r2%r(:)) ) then
          overlap = .true. 
          return
       end if

    end do
    overlap = .false.

  end function region_overlaps

  subroutine region_print(r, seq_max)
    type(tRegion), intent(in) :: r
    integer, intent(in), optional :: seq_max
 
    integer :: i, lseq_max

    lseq_max = 7
    if ( present(seq_max) ) lseq_max = seq_max

    write(*,'(a,i0,2a)') 'Region (',r%n,'): ',trim(r%name)
    write(*,'(tr2,a)',advance='no') '['
    do i = 1 , r%n - 1
       write(*,'(tr1,i0,a)',advance='no') r%r(i),', '
       if ( mod(i,lseq_max) == 0 ) then
          write(*,'(/,tr3)',advance='no')
       end if
    end do
    if ( r%n > 0 ) then
       write(*,'(tr1,i0,a)') r%r(r%n),' ]'
    else
       write(*,'(tr1,a)') ' ]'
    end if

  end subroutine region_print
    

#ifdef MPI
  subroutine region_MPI_union(dit,r)
    use mpi_siesta, only : MPI_AllReduce, MPI_Integer, MPI_Sum
    use mpi_siesta, only : MPI_Bcast
    use mpi_siesta, only : MPI_Recv, MPI_Send, MPI_STATUS_SIZE
    use mpi_siesta, only : MPI_Get_Count

    type(OrbitalDistribution), intent(in) :: dit
    type(tRegion), intent(inout) :: r
    
    ! Our temporary region
    integer :: nt, ct, iN, it
    integer, allocatable :: rd(:)
    character(len=R_NAME_LEN) :: tmp
    integer :: comm

    integer :: MPIerror, MPIstatus(MPI_STATUS_SIZE)

    ! Get the communicator assigned to this
    ! distribution
    comm = dist_comm(dit)
    if ( dist_nodes(dit) == 1 ) return
    call MPI_AllReduce(r%n,nt,1,MPI_Integer, MPI_Sum, &
         comm, MPIerror)

    ! Allocate space for the information
    allocate(rd(nt+1))

    ! We let the first node in the distribution collect
    ! the data, then we b-cast it...
    if ( dist_node(dit) == 0 ) then
       if ( r%n > 0 ) then
          rd(1:r%n) = r%r(:)
       end if
       ct = r%n + 1
       do iN = 1 , dist_nodes(dit) - 1
          call MPI_Recv(rd(ct),nt-ct+1,MPI_Integer, &
               iN,0, comm, MPIstatus, MPIerror)
          call MPI_Get_Count(MPIstatus, MPI_Integer, it, MPIerror)
          ct = ct + it
       end do
       ! Count the actual number of unique entries
       nt = uniqc(rd(1:nt-1))
       ! Sort them...
       ct = 1
       it = 1
       do while ( ct < nt )
          it = it + 1
          if ( all(rd(it) /= rd(1:ct)) ) then
             ct = ct + 1
             rd(ct) = rd(it)
          end if
       end do
    else
       if ( r%n == 0 ) then
          call MPI_Send(rd(1),0,MPI_Integer, &
               0, 0, comm, MPIerror)
       else
          call MPI_Send(r%r(1),r%n,MPI_Integer, &
               0, 0, comm, MPIerror)
       end if
    end if
    
    ! Reduce the size of the array to the 
    ! actual size...
    ! The problem with distributions is their
    ! ability to discontinue the natural order...
    ! We should probably use a single optimization of the 
    ! sparsity pattern to reduce the bandwidth.

    call MPI_Bcast(nt,1,MPI_Integer, 0, comm, MPIerror)
    call MPI_Bcast(rd(1),nt,MPI_Integer, 0, comm, MPIerror)

    ! Deallocate and copy the new array
    tmp = r%name
    call region_list(r,nt,rd(1:nt),name=tmp)
    
    ! Clean-up
    deallocate(rd)

  end subroutine region_MPI_union

  subroutine region_MPI_Bcast(r,Bnode,Comm)
    use mpi_siesta, only : MPI_Bcast, MPI_Integer
    type(tRegion), intent(inout) :: r
    integer, intent(in) :: Bnode, Comm
    
    integer :: Node
    integer :: MPIerror

    call MPI_Comm_Rank(Comm,Node,MPIerror)

    call MPI_Bcast(r%n,1,MPI_Integer, Bnode, comm, MPIerror)
    if ( Node /= Bnode ) then
       allocate(r%r(r%n))
    end if
    call MPI_Bcast(r%r(1),r%n,MPI_Integer, Bnode, comm, MPIerror)

  end subroutine region_MPI_Bcast

#endif
  
end module m_region
