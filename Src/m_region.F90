
! This module provides handles for holding regions of orbitals/atoms.
! It is a convenient in-place type to easily handle
! different kinds of regions of any kind.
! It easily handles connections with regions to
! generate limited sets of regions (or expand them).

! Fully created by Nick Papior Andersen, 2014

module m_region
  
  use intrinsic_missing, only : UNIQ, UNIQC, SFIND
  use intrinsic_missing, only : SORT_QUICK
  use geom_helper, only : UCORB, IAORB

  use class_OrbitalDistribution
  use class_Sparsity

  implicit none

  private

  ! max length of region
  integer, parameter, public :: R_NAME_LEN = 50

  type :: tRgn
     ! The name of the region
     character(len=R_NAME_LEN) :: name = ' '
     ! The quantities in the region
     integer :: n = 0
     integer, pointer :: r(:) => null()
     logical :: sorted = .false.
  end type tRgn

  ! A linked list of regions
  type :: tRgnLL
     type(tRgn) :: rgn
     type(tRgnLL), pointer :: next => null()
  end type tRgnLL

  public :: tRgn, tRgnLL
  public :: rgn_init, rgn_init_pvt
  public :: rgn_assoc
  public :: rgn_delete, rgn_nullify, rgnll_delete
  public :: rgn_intersection
  public :: rgn_union, rgn_append
  public :: rgn_union_complement, rgn_complement
  public :: rgn_range
  public :: rgn_list
  public :: rgn_overlaps
  public :: rgn_sort, rgn_uniq
  public :: rgn_insert, rgn_remove
  public :: rgn_print
  public :: rgn_copy
  public :: rgn_reverse
  public :: rgn_wrap
  public :: in_rgn, rgn_pivot
  interface rgn_push
     module procedure rgn_push_val
     module procedure rgn_push_rgn
  end interface rgn_push
  public :: rgn_push, rgn_pop
#ifdef MPI
  public :: rgn_MPI_union
  public :: rgn_MPI_Bcast
#endif

  ! Regions which has to do with sparsity patterns
  public :: rgn_sp_connect
  interface rgn_sp_sort
     module procedure rgn_sp_sort_type
     module procedure rgn_sp_sort_explicit
  end interface rgn_sp_sort
  public :: rgn_sp_sort

  ! Regions which has to do with atoms/orbitals
  public :: rgn_correct_atom
  public :: rgn_atom2orb
  public :: rgn_orb2atom

  interface sum
     module procedure rgn_sum
  end interface sum
  public :: rgn_sum, sum

  ! Sorting methods for the regions
  ! We sort by the number of connections forward
  ! the maximum number of connections to elements
  ! not already in the region will be put lastly
  integer, parameter, public :: R_SORT_MAX_FRONT = 1
  ! We sort by the number of connections backward back into
  ! the region behind the current region.
  ! The one with the most connections back will be put first
  integer, parameter, public :: R_SORT_MAX_BACK = 2
  ! The one that reaches the farthest back into the region
  integer, parameter, public :: R_SORT_LONGEST_BACK = 3

  interface rgn_remove
     module procedure rgn_remove_rgn
     module procedure rgn_remove_list
  end interface rgn_remove

  interface rgn_insert
     module procedure rgn_insert_rgn
     module procedure rgn_insert_list
  end interface rgn_insert

contains

  recursive subroutine rgnll_delete(rll)
    type(tRgnLL), intent(inout) :: rll

    ! clean the region
    call rgn_delete(rll%rgn)

    if ( associated(rll%next) ) then
       call rgnll_delete(rll%next)
       deallocate(rll%next)
       nullify(rll%next)
    end if
    
  end subroutine rgnll_delete


  ! Easy initialization a region with optional name
  subroutine rgn_init(r,n,name,val)
    ! The region containing of size n
    type(tRgn), intent(inout) :: r
    ! size of the range
    integer, intent(in) :: n
    ! The optional name of the range
    character(len=*), intent(in), optional :: name
    ! default value (none will be set if not present)
    integer, intent(in), optional :: val
    
    call rgn_delete(r)

    ! Quick return for empty region
    if ( n == 0 ) return

    ! Create with no default values
    r%n = n
    allocate(r%r(n))
    call memory('A','I',n,'rgn-list')
    if ( present(val) ) r%r = val

    if ( present(name) ) then
       r%name = name
    end if
    
  end subroutine rgn_init

  ! From return in pvt the pivoting
  ! table of r2 in r1.
  ! Hence:
  !   pvt%r(1) == rgn_pivot(r1,r2%r(1))
  subroutine rgn_init_pvt(r1,r2,pvt)
    type(tRgn), intent(in) :: r1, r2
    type(tRgn), intent(inout) :: pvt

    integer :: i

    ! Initialize to zero if the elements
    ! are not aligned
    call rgn_init(pvt,r2%n,val = 0)
    pvt%name = trim(r2%name)//' -> '//trim(r1%name)
    
    do i = 1 , r2%n
       pvt%r(i) = rgn_pivot(r1,r2%r(i))
    end do

  end subroutine rgn_init_pvt
    

  ! Get the index of the element that equal 'i'
  ! We refer to this as the pivoting element.
  ! In general this is the same as doing:
  !   A(r%r == i)
  elemental function rgn_pivot(r,i) result(idx)
    type(tRgn), intent(in) :: r
    integer, intent(in) :: i
    integer :: idx
    if ( r%sorted ) then
       idx = 0
       if ( r%n == 0 ) return
       idx = SFIND(r%r(1:r%n),i)
    else
       do idx = 1 , r%n
          if ( r%r(idx) == i ) return
       end do
       idx = 0
    end if
  end function rgn_pivot

  ! Check whether an element exists in this region.
  ! Performs an easy check instead of doing it manually.
  elemental function in_rgn(r,i) result(in)
    type(tRgn), intent(in) :: r
    integer, intent(in) :: i
    logical :: in
    if ( r%n == 0 ) then
       in = .false.
       return
    end if
    if ( r%sorted ) then
       in = SFIND(r%r(1:r%n),i) > 0
    else
       in = any(i == r%r(1:r%n))
    end if
  end function in_rgn

  ! Copies one region to another one
  ! This emplies that the element 'to' is deleted before
  ! copy is performed.
  subroutine rgn_copy(from,to)
    type(tRgn), intent(in) :: from
    type(tRgn), intent(inout) :: to
    if ( from%n == 0 ) then
       call rgn_delete(to)
       to%name = from%name
    else if ( .not. associated(to%r,from%r) ) then
       call rgn_list(to,from%n,from%r,name=from%name)
       to%sorted = from%sorted
!    else
       ! The to and from have the same memory address
       ! Hence we can skip copying anything
       !print*,'skipping copy, same rgn'
    end if
  end subroutine rgn_copy

  ! Fully deletes a region (irrespective of it's current state)
  recursive subroutine rgn_delete(r,r1,r2,r3,r4,r5)
    type(tRgn), intent(inout) :: r
    type(tRgn), intent(inout), optional :: r1,r2,r3,r4,r5
    r%name = ' '
    r%n = 0
    if ( associated(r%r) ) then
       call memory('D','I',size(r%r),'rgn-list')
       deallocate(r%r)
    end if
    nullify(r%r)
    r%sorted = .false.
    if ( present(r1) ) call rgn_delete(r1,r2,r3,r4,r5)
  end subroutine rgn_delete

  ! Fully nullifies a region (irrespective of it's current state)
  recursive subroutine rgn_nullify(r,r1,r2,r3,r4,r5)
    type(tRgn), intent(inout) :: r
    type(tRgn), intent(inout), optional :: r1,r2,r3,r4,r5
    nullify(r%r) ! this pro-hibits the deletion of the array
    call rgn_delete(r)
    if ( present(r1) ) call rgn_nullify(r1,r2,r3,r4,r5)
  end subroutine rgn_nullify


  ! Assigns the left region with the right region
  ! hence, deleting one will create a memory leak if
  ! not handled carefully
  subroutine rgn_assoc(lhs,rhs, dealloc)
    type(tRgn), intent(inout) :: lhs, rhs
    ! whether a pre-deallocation of the lhs should occur
    logical, intent(in), optional :: dealloc
    if ( present(dealloc) ) then
       if ( dealloc ) call rgn_delete(lhs)
    end if
    ! default to nullify
    call rgn_nullify(lhs)
    lhs%name = rhs%name
    lhs%n = rhs%n
    lhs%r => rhs%r
    lhs%sorted = rhs%sorted

  end subroutine rgn_assoc

  ! Allows removing certain elements from a region.
  subroutine rgn_remove_list(r,n,list,rout)
    type(tRgn), intent(in) :: r
    type(tRgn), intent(inout) :: rout
    integer, intent(in) :: n, list(n)

    integer :: i, ni
    integer, allocatable :: tmp_r(:)

    if ( n == 0 ) then
       call rgn_copy(r,rout)
       return
    end if

    allocate(tmp_r(r%n))
    ni = 0
    do i = 1 , r%n
       if ( all(r%r(i)/=list(:)) ) then
          ni = ni + 1
          tmp_r(ni) = r%r(i)
       end if
    end do

    call rgn_list(rout,ni,tmp_r)

    deallocate(tmp_r)

  end subroutine rgn_remove_list
  subroutine rgn_remove_rgn(r,rr,rout)
    type(tRgn), intent(in) :: r
    type(tRgn), intent(in) :: rr
    type(tRgn), intent(inout) :: rout
    
    call rgn_remove_list(r,rr%n,rr%r,rout)

  end subroutine rgn_remove_rgn


  ! Returns a list of elements which all
  ! are between [1;r%n]
  ! Hence they will be wrapped
  ! An optional size can be specified
  ! if the array is shorter
  subroutine rgn_wrap(r,n)
    type(tRgn), intent(inout) :: r
    integer, intent(in), optional :: n
    integer :: i, ln

    ln = r%n
    if ( ln <= 0 ) return
    if ( present(n) ) ln = n

    i = 1
    do while ( i <= r%n )
       if ( r%r(i) < 0 ) then
          r%r(i) = ln + r%r(i) + 1
       else if ( r%r(i) == 0 ) then
          r%r(i) = ln
       else if ( ln < r%r(i) ) then
          r%r(i) = r%r(i) - ln
       else
          i = i + 1
       end if
    end do

  end subroutine rgn_wrap


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
  subroutine rgn_sp_connect(r,dit,sp,cr,except, connect_from, follow)

    ! the region we wish to find the connections to
    type(tRgn), intent(in) :: r
    ! The sparsity patterns distribution
    type(OrbitalDistribution), intent(in) :: dit
    ! The sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! the connecting region
    type(tRgn), intent(inout) :: cr
    ! a region which must not be crossed (exception region)
    type(tRgn), intent(in), optional :: except
    ! a region which contains the orbitals that are 
    ! responsible for the connection
    type(tRgn), intent(inout), optional :: connect_from
    ! Tell whether we should follow the connection
    ! or stop at the first iteration
    logical, intent(in), optional :: follow

    ! ** local variables
    type(tRgn) :: tmp, tmp2
    integer :: i, j, io, jo, ind, no_l, no_u, it, rt
    integer, allocatable :: ct(:), rr(:)
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)

    call attach(sp,nrows=no_l,nrows_g=no_u, &
         n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)

    ! ensure nullification of cr
    call rgn_delete(cr)

    allocate(rr(r%n))
    allocate(ct(no_u-r%n))

    call rgn_copy(r, tmp)
    ! In case r%n is extremely big we can with benefit &
    ! search the region in a sorted array.
    ! We sort it once, and search that array instead
    ! This requires a slightly increased memory, but 
    ! drastically improves performance.
    call rgn_sort(tmp)
    
    if ( present(except) ) then
       ! Exclude it's own region, (no connection back)
       if ( except%n > 0 ) then
          call rgn_copy(tmp,tmp2)
          call rgn_insert(tmp2,except,tmp,0)
          call rgn_delete(tmp2)
          call rgn_sort(tmp)
       end if
    end if

    rt = 0
    ! A duplicate of ct for faster checks (sorted search)
    call rgn_init(tmp2, no_u-r%n)
    tmp2%n = 0
    
    do i = 1 , r%n

       ! Orbital that should be folded from
       io = index_global_to_local(dit, r%r(i))
       if ( io <= 0 ) cycle ! the orbital does not exist on this node

       ! Count number of orbitals that are connected
       j = tmp2%n
       do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)

          ! UC-orb
          jo = ucorb(l_col(ind), no_u)

          ! Ensure that it is not a folding to the same region
          ! tmp is sorted, hence it uses SFIND
          if ( in_rgn(tmp, jo) ) cycle

          if ( .not. in_rgn(tmp2, jo) ) then
             tmp2%n = tmp2%n + 1
             ct(tmp2%n) = jo
             ! Retain two lists (in case there
             ! are super-cell connections we must do
             ! the sorting here, besides)
             tmp2%r(tmp2%n) = jo
             call rgn_sort(tmp2)
          end if

       end do
       
       if ( j /= tmp2%n ) then
          ! The connect_from region
          ! will only make sense if we do not follow the orbitals
          ! (else we still get where the connection starts)!
          ! This will never double check as
          ! we already have each orbital only once per core...
          rt = rt + 1
          rr(rt) = io
       end if

    end do

    ! Store the current size of ct
    it = tmp2%n

    ! If we are supposed to follow
    ! then loop and extend ct
    if ( present(follow) ) then
    if ( follow ) then

#ifdef MPI
       ! We have found the first group of regions
       ! Now we need to align the regions for all the processors
       ! in this orbital distribution (before we move on 
       ! and check for "follow"

       if ( dist_nodes(dit) > 1 ) then
          ! b-cast the connected to region
          call rgn_list(tmp2, it, ct)
          call rgn_MPI_union(dit, tmp2)
          it = tmp2%n
          if ( it > 0 ) then
             do i = 1 , it
                ct(i) = tmp2%r(i)
             end do
          end if
          
          ! A duplicate of ct for faster checks (sorted search)
          call rgn_init(tmp2, no_u-r%n)
          do i = 1 , it
             tmp2%r(i) = ct(i)
          end do
          tmp2%n = it

          ! Ensure it is sorted.
          call rgn_sort(tmp2)
       end if
#else
       ! Do nothing, the tmp2 array is already the correct size
       ! and sorted.
#endif

       i = 1
       do while ( i <= it )

          ! Orbital that should be folded from
          io = ct(i)

          do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)

             ! UC-orb
             jo = ucorb(l_col(ind), no_u)

             ! Ensure that it is not a folding to the same region
             if ( in_rgn(tmp, jo) ) cycle

             if ( .not. in_rgn(tmp2, jo) ) then
                it = it + 1
                ct(it) = jo
                ! Retain two lists (in case there
                ! are super-cell connections we must do
                ! the sorting here, besides)
                tmp2%n = it
                tmp2%r(it) = jo
                call rgn_sort(tmp2)
             end if
             
          end do
          
          ! step the orbital that we now follow
          i = i + 1

       end do
       
    end if
    end if

    ! Clean-up
    call rgn_delete(tmp2)

    ! We now have a list of orbitals that needs to be folded to
    ! Copy the list over
    call rgn_list(cr, it, ct)

    deallocate(ct)

#ifdef MPI
    call rgn_MPI_union(dit, cr)
#endif

    if ( present(connect_from) ) then

       call rgn_list(connect_from, rt, rr)

#ifdef MPI
       call rgn_MPI_union(dit, connect_from)
#endif

    end if

    call rgn_delete(tmp)
    deallocate(rr)

  end subroutine rgn_sp_connect


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
  subroutine rgn_sp_sort_type(r,dit,sp, sr, method)

#ifdef MPI
    use mpi_siesta, only : MPI_Integer
    use mpi_siesta, only : MPI_MAX, MPI_MIN, MPI_AllReduce
#endif

    ! the region we wish to find the connections to
    type(tRgn), intent(inout) :: r
    ! The sparsity patterns distribution
    type(OrbitalDistribution), intent(in) :: dit
    ! The sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! the sorting region (i.e. the orbitals that are allowed to be pivoted)
    type(tRgn), intent(in) :: sr
    ! The method used for sorting
    integer, intent(in) :: method

    ! ** local variables
    integer :: n, nzs
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)

    if ( r%n == 0 ) return

    ! Attach to the sparsity pattern...
    call attach(sp,nrows_g=n, nnzs=nzs, &
         n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)

    call rgn_sp_sort_explicit(r,n,nzs,l_ncol,l_ptr,l_col, &
         sr,method,dit=dit)

  end subroutine rgn_sp_sort_type

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
  subroutine rgn_sp_sort_explicit(r,n,nnzs,n_col,l_ptr,l_col, sr, method,dit)

#ifdef MPI
    use mpi_siesta, only : MPI_Integer
    use mpi_siesta, only : MPI_MAX, MPI_MIN, MPI_AllReduce
#endif

    ! the region we wish to find the connections to
    type(tRgn), intent(inout) :: r
    ! The sparsity pattern
    integer, intent(in) :: n, nnzs, n_col(n), l_ptr(n), l_col(nnzs)
    ! the sorting region (i.e. the orbitals that are allowed to be pivoted)
    type(tRgn), intent(in) :: sr
    ! The method used for sorting
    integer, intent(in) :: method
    ! The distribution
    type(OrbitalDistribution), intent(in), optional :: dit

    ! ** local variables
    type(tRgn) :: r_tmp, cc
    integer :: i, io, ind, jo, ci, last_n, idx_max
    integer, allocatable :: n_c(:)
#ifdef MPI
    integer :: comm
    integer :: MPIerror
#endif

    if ( r%n == 0 ) return

#ifdef MPI
    if ( present(dit) ) comm = dist_comm(dit)
#endif

    ! First we ensure that the lists are the same
    call rgn_intersection(r,sr,r_tmp)
    if ( r_tmp%n /= sr%n ) then
       call die('The regions are not well-defined')
    end if
    call rgn_delete(r_tmp)

    ! Prepare the collection arrays
    allocate(n_c(sr%n))
    n_c = 0 ! initialize to ensure MPI reduction
    call rgn_init(cc, n)
    cc%n = 0

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
          if ( present(dit) ) then
             io = index_global_to_local(dit,sr%r(i))
             if ( io <= 0 ) cycle ! the orbital does not exist on this node
          else
             io = sr%r(i)
          end if
          
          if ( n_col(io) == 0 ) cycle

          cc%n = 0
          do ind = l_ptr(io) + 1 , l_ptr(io) + n_col(io)
             
             ! UC-orb
             jo = ucorb(l_col(ind),n)
             
             ! Ensure that it is not folding to the same region
             if ( in_rgn(r,jo) ) cycle

             if ( .not. in_rgn(cc, jo) ) then
                cc%n = cc%n + 1
                cc%r(cc%n) = jo
                call rgn_sort(cc)
             end if

          end do

          ! As the region is created "from a frontal search"
          ! we will reverse the result (maxloc returns the first max)
          n_c(sr%n-i+1) = cc%n

       end do

#ifdef MPI
       if ( present(dit) ) then
          if ( dist_nodes(dit) > 1 ) then
             call MPI_AllReduce(n_c,cc%r,sr%n,MPI_Integer, &
                  MPI_MAX, comm, MPIerror)
             do i = 1 , sr%n
                n_c(i) = cc%r(i)
             end do
          end if
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
          if ( present(dit) ) then
             io = index_global_to_local(dit,sr%r(i))
             if ( io <= 0 ) cycle ! the orbital does not exist on this node
          else
             io = sr%r(i)
          end if

          if ( n_col(io) == 0 ) cycle

          cc%n = 0
          do ind = l_ptr(io) + 1 , l_ptr(io) + n_col(io)
             
             ! UC-orb
             jo = ucorb(l_col(ind),n)
             
             ! Ensure that it is folding to the same region
             if ( .not. in_rgn(r,jo) ) cycle

             if ( .not. in_rgn(cc, jo) ) then
                cc%n = cc%n + 1
                cc%r(cc%n) = jo
                call rgn_sort(cc)
             end if

          end do

          ! As the region is created "from a frontal search"
          ! we will reverse the result (maxloc returns the first max)
          n_c(sr%n-i+1) = cc%n

       end do

#ifdef MPI
       if ( present(dit) ) then
          if ( dist_nodes(dit) > 1 ) then
             call MPI_AllReduce(n_c,cc%r,sr%n,MPI_Integer, &
                  MPI_MAX, comm, MPIerror)
             do i = 1 , sr%n
                n_c(i) = cc%r(i)
             end do
          end if
       end if
#endif

       ! We now reverse the values as we then control that the one
       ! with the least connections back, will be placed last
       jo = maxval(n_c, dim=1)
       do i = 1 , sr%n
          n_c(i) = jo - n_c(i)
       end do

    case ( R_SORT_LONGEST_BACK ) 

       ! we sort and move the orbitals
       ! back in the rank according to their longest
       ! range backwards.
       ! A longer connection back inside of the region
       ! will move it further back.

       ! Step 1. find the number of connections for
       ! each orbital in the sorting region
       
       do i = 1 , sr%n

          ! Orbital that should be folded from
          if ( present(dit) ) then
             io = index_global_to_local(dit,sr%r(i))
             if ( io <= 0 ) cycle ! the orbital does not exist on this node
          else
             io = sr%r(i)
          end if

          if ( n_col(io) == 0 ) cycle

          ci = r%n
          do ind = l_ptr(io) + 1 , l_ptr(io) + n_col(io)
             
             ! UC-orb
             jo = ucorb(l_col(ind),n)
             
             ! Ensure that it is folding to the same region
             ci = min(rgn_pivot(r,jo),ci)

          end do

          ! As the region is created "from a frontal search"
          ! we will reverse the result (maxloc returns the first max)
          n_c(i) = ci

       end do

#ifdef MPI
       if ( present(dit) ) then
          if ( dist_nodes(dit) > 1 ) then
             call MPI_AllReduce(n_c,cc%r,sr%n,MPI_Integer, &
                  MPI_MIN, comm, MPIerror)
             do i = 1 , sr%n
                n_c(i) = cc%r(i)
             end do
          end if
       end if
#endif

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
          do i = last_n , 1 , -1
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

    case ( R_SORT_LONGEST_BACK )

       last_n = r%n - sr%n + 1

       do jo = 1 , sr%n

          ! The farthest back
          idx_max = minloc(n_c,dim=1)
          
          io = sr%r(idx_max)
          do i = last_n , r%n
             if ( r%r(i) == io ) then
                ! swap the positions
                r%r(i) = r%r(last_n)
                r%r(last_n) = io
                last_n = last_n + 1
                exit ! we are not going to find it later...
             end if
          end do

          n_c(idx_max) = r%n+1 ! reset connections, we already
                               ! processed the orbital

       end do

       ! We have now sorted the entries according to the maximum
       ! connections across the existing region.

    end select

    call rgn_delete(cc)
    deallocate(n_c)

    if ( r%sorted ) then
       do i = 2 , r%n
          if ( r%r(i-1) > r%r(i) ) then
             r%sorted = .false.
             exit
          end if
       end do
    end if

  end subroutine rgn_sp_sort_explicit


  ! Creates a region in terms of the INTERSECTION of two regions 
  subroutine rgn_intersection(r1,r2,ir)
    ! the regions we wish to operate on
    type(tRgn), intent(in) :: r1, r2
    ! the intersection region
    type(tRgn), intent(inout) :: ir

    ! ** local variables
    type(tRgn) :: sr
    integer :: i, it
    integer, allocatable :: ct(:)

    if ( r1%n == 0 .or. r2%n == 0 ) then
       call rgn_delete(ir)
       return
    end if

    ! Having a sorted array will greatly increase performance
    ! at minimal memory cost
    if ( r2%sorted ) then
       sr%sorted = .true.
       sr%n = r2%n
       sr%r => r2%r
    else
       call rgn_copy(r2,sr)
       call rgn_sort(sr)
    end if

    it = 0
    allocate(ct(min(r1%n,r2%n)))

    do i = 1 , r1%n
       
       if ( in_rgn(sr,r1%r(i)) ) then
          it = it + 1
          ct(it) = r1%r(i)
       end if

    end do

    if ( .not. r2%sorted ) call rgn_delete(sr)

    ! We now have a list of orbitals that needs to be folded to
    ! Copy the list over
    call rgn_list(ir,it,ct)

    deallocate(ct)
    
  end subroutine rgn_intersection

  ! Creates a region in terms of the UNION of two regions 
  subroutine rgn_union(r1,r2,ur)
    ! the regions we wish to operate on
    type(tRgn), intent(in) :: r1, r2
    ! the unified region
    type(tRgn), intent(inout) :: ur

    ! ** local variables
    integer :: i, it
    type(tRgn) :: sr
    integer, allocatable :: ct(:)

    if ( r1%n == 0 ) then
       call rgn_copy(r2,ur)
       return
    else if ( r2%n == 0 ) then
       call rgn_copy(r1,ur)
       return
    end if

    allocate(ct(r1%n+r2%n))

    ! Having a sorted array will greatly increase performance
    ! at minimal memory cost
    if ( r1%sorted ) then
       sr%sorted = .true.
       sr%n = r1%n
       sr%r => r1%r
    else
       call rgn_copy(r1,sr)
       call rgn_sort(sr)
    end if

    ! Copy over r1
    do i = 1 , r1%n
       ct(i) = r1%r(i)
    end do

    it = r1%n
    do i = 1 , r2%n
       
       if ( in_rgn(sr,r2%r(i)) ) cycle
       ! the above does not remove dublicates in r2 :(
       !if ( any(r2%r(i) == ct(1:it)) ) cycle

       it = it + 1
       ct(it) = r2%r(i)

    end do

    if ( .not. r1%sorted ) call rgn_delete(sr)

    ! We now have a list of orbitals that needs to be folded to
    ! Copy the list over
    call rgn_list(ur,it,ct)

    deallocate(ct)

  end subroutine rgn_union

  ! Creates a region by appending them, no checks of dublicates
  subroutine rgn_append(r1,r2,r)
    ! the regions we wish to operate on
    type(tRgn), intent(in) :: r1, r2
    ! the total region
    type(tRgn), intent(inout) :: r

    ! ** local variables
    logical :: sorted
    integer :: i
    integer, allocatable :: ct(:)

    if ( r1%n == 0 ) then
       call rgn_copy(r2,r)
       return
    else if ( r2%n == 0 ) then
       call rgn_copy(r1,r)
       return
    end if

    ! To retain the sorted attribute of the new region
    sorted = r1%sorted .and. r2%sorted .and. r1%r(r1%n) <= r2%r(1)
    
    allocate(ct(r1%n+r2%n))

    ! Copy over r1
    do i = 1 , r1%n
       ct(i) = r1%r(i)
    end do
    ! Copy over r2
    do i = 1 , r2%n
       ct(r1%n+i) = r2%r(i)
    end do

    ! We now have a list of orbitals that needs to be folded to
    ! Copy the list over
    call rgn_list(r,size(ct),ct)

    r%sorted = sorted

    deallocate(ct)

  end subroutine rgn_append

  subroutine rgn_insert_list(r,n,list,rout,idx)
    type(tRgn), intent(in) :: r
    integer, intent(in) :: n, list(n)
    type(tRgn), intent(inout) :: rout
    ! The place of insertion
    integer, intent(in) :: idx
    type(tRgn) :: rtmp
    call rgn_list(rtmp,n,list)
    call rgn_insert_rgn(r,rtmp,rout,idx)
    call rgn_delete(rtmp)
  end subroutine rgn_insert_list

  ! Inserts a region in another region.
  subroutine rgn_insert_rgn(r,rin,rout,idx)
    ! The region we wish to insert in
    type(tRgn), intent(in) :: r
    type(tRgn), intent(in) :: rin
    type(tRgn), intent(inout) :: rout
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
       call rgn_union(r,rin,rout)
       return
    else if ( idx <= 0 ) then
       call rgn_union(rin,r,rout)
       return
    end if

    ! We are placing it somewhere in the middle

    ! Insert list...
    ! This is not soo easy.
    ! First we need to find the first value before
    ! the list (in case it already exists)
    lidx = idx
    do while ( in_rgn(rin,r%r(lidx)) )
       lidx = lidx - 1
       if ( lidx == 0 ) then
          call rgn_union(rin,r,rout)
          return
       end if
    end do

    ! Grab value before removing potential duplicates
    i = r%r(lidx)
    call rgn_remove(r,rin,rout)
    ! Back-retrieve the new position
    i = rgn_pivot(rout,i)
    if ( i == 0 ) call die('Error in algorithm')
    allocate(itmp(rout%n+rin%n))
    itmp(1:i) = rout%r(1:i)
    itmp(i+1:i+rin%n) = rin%r(1:rin%n)
    j = i + rin%n + 1
    itmp(j:) = rout%r(i+1:rout%n)
    call rgn_list(rout,size(itmp),itmp)
    deallocate(itmp)

  end subroutine rgn_insert_rgn

  ! Creates a region in terms of the COMPLEMENT of two regions 
  ! Hence it returns elements in r2 which are not in r1
  subroutine rgn_complement(r1,r2,cr)
    ! the regions we wish to operate on
    type(tRgn), intent(in) :: r1, r2
    ! the complementary region
    type(tRgn), intent(inout) :: cr

    ! ** local variables
    type(tRgn) :: sr
    integer :: i, it
    integer, allocatable :: ct(:)

    if ( r1%n == 0 ) then
       call rgn_copy(r2,cr)
       return
    end if
    if ( r2%n == 0 ) then
       call rgn_delete(cr)
       return
    end if

    allocate(ct(r2%n))

    ! Having a sorted array will greatly increase performance
    ! at minimal memory cost
    if ( r1%sorted ) then
       sr%sorted = .true.
       sr%n = r1%n
       sr%r => r1%r
    else
       call rgn_copy(r1,sr)
       call rgn_sort(sr)
    end if

    it = 0
    do i = 1 , r2%n
       
       if ( .not. in_rgn(sr,r2%r(i)) ) then
          it = it + 1
          ct(it) = r2%r(i)
       end if

    end do

    if ( .not. r1%sorted ) call rgn_delete(sr)

    ! Copy the list over
    call rgn_list(cr,it,ct)

    deallocate(ct)

  end subroutine rgn_complement

  ! Creates a region in terms of the COMPLEMENT of two regions.
  ! Hence it returns elements in r2 not in r1, AND the elements
  ! in r1 not in r2
  subroutine rgn_union_complement(r1,r2,cr)
    ! the regions we wish to operate on
    type(tRgn), intent(in) :: r1, r2
    ! the complementary region
    type(tRgn), intent(inout) :: cr

    ! ** local variables
    type(tRgn) :: sr
    integer :: i, it
    integer, allocatable :: ct(:)

    if ( r1%n == 0 ) then
       call rgn_copy(r2,cr)
       return
    end if
    if ( r2%n == 0 ) then
       call rgn_copy(r1,cr)
       return
    end if

    allocate(ct(r1%n+r2%n))

    if ( r2%sorted ) then
       sr%sorted = .true.
       sr%n = r2%n
       sr%r => r2%r
    else
       call rgn_copy(r2,sr)
       call rgn_sort(sr)
    end if

    it = 0
    do i = 1 , r1%n
       if ( in_rgn(sr,r1%r(i)) ) cycle

       it = it + 1
       ct(it) = r1%r(i)

    end do
    if ( r2%sorted ) then
       call rgn_nullify(sr)
    else
       call rgn_delete(sr)
    end if
    if ( r1%sorted ) then
       sr%sorted = .true.
       sr%n = r1%n
       sr%r => r1%r
    else
       call rgn_copy(r1,sr)
       call rgn_sort(sr)
    end if
    do i = 1 , r2%n
       if ( in_rgn(sr,r2%r(i)) ) cycle

       it = it + 1
       ct(it) = r2%r(i)

    end do

    if ( .not. r1%sorted ) call rgn_delete(sr)

    ! Copy the list over
    call rgn_list(cr,it,ct)

    deallocate(ct)

  end subroutine rgn_union_complement

  ! Easy setup of a consecutive orbital range
  subroutine rgn_range(r,o1,o2)
    ! The region containing the range [o1;o2]
    type(tRgn), intent(inout) :: r
    ! The limits on the range
    integer, intent(in) :: o1, o2
    integer :: io

    io = abs(o2 - o1) + 1
    call rgn_init(r,io)

    if ( o1 <= o2 ) then
       do io = o1 , o2
          r%r(io-o1+1) = io
       end do
       ! It is sorted
       r%sorted = .true.
    else
       do io = o1 , o2 , -1
          r%r(o1-io+1) = io
       end do
    end if

  end subroutine rgn_range

  subroutine rgn_list(r,n,list,name)
    ! Region to put list in
    type(tRgn), intent(inout) :: r
    ! list to copy over
    integer, intent(in) :: n, list(n)
    character(len=*), intent(in), optional :: name
    integer :: i

    call rgn_delete(r)
    r%n = n
    if ( n > 0 ) then
       allocate(r%r(n))
       call memory('A','I',n,'rgn-list')
       do i = 1 , n
          r%r(i) = list(i)
       end do
    end if
    if ( present(name) ) r%name = name
    
  end subroutine rgn_list

  subroutine rgn_sort(r)
    type(tRgn), intent(inout) :: r
    if ( r%n > 0 ) then
       call sort_quick(r%n, r%r)
       r%sorted = .true.
    end if
  end subroutine rgn_sort

  subroutine rgn_uniq(r)
    type(tRgn), intent(inout) :: r
    character(len=R_NAME_LEN) :: name
    integer, allocatable :: ct(:)
    integer :: n

    if ( r%n > 0 ) then
       name = r%name
       n = uniqc(r%r(1:r%n))
       allocate(ct(n))
       ct(1:n) = uniq(r%r(1:r%n))
       call rgn_list(r,n,ct,name)
       deallocate(ct)
       call rgn_sort(r)
    end if

  end subroutine rgn_uniq

  subroutine rgn_reverse(r)
    type(tRgn), intent(inout) :: r
    integer :: i, tmp
    ! Reverse the list
    do i = 1 , r%n / 2
       tmp = r%r(i)
       r%r(i) = r%r(r%n+1-i)
       r%r(r%n+1-i) = tmp
    end do
  end subroutine rgn_reverse

  subroutine rgn_correct_atom(r,na_u,lasto)
    ! Region which we want to extend with the atoms
    ! that it already overlaps
    type(tRgn), intent(inout) :: r
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
    call rgn_delete(r)
    r%name = tmp
    r%n = no
    nullify(r%r)
    if ( r%n > 0 ) then
       allocate(r%r(no))
       call memory('A','I',no,'rgn-list')
       no = 0
       do i = 1 , na
          do j = lasto(a(i)-1) + 1 , lasto(a(i))
             no = no + 1
             r%r(no) = j
          end do
       end do
    end if

  end subroutine rgn_correct_atom

  subroutine rgn_Atom2Orb(ar,na_u,lasto,or)
    ! Region which we want to extend with the atoms
    ! that it already overlaps
    type(tRgn), intent(in) :: ar
    ! The last orbitals of each atom
    integer, intent(in) :: na_u, lasto(0:na_u)
    type(tRgn), intent(inout) :: or

    ! ** local variables
    integer :: i, cr, no
    character(len=R_NAME_LEN) :: tmp

    if ( ar%n == 0 ) then
       call rgn_delete(or)
       return
    end if

    ! First calculate size of region in orbital
    ! space
    no = 0
    do i = 1 , ar%n
       no = no + lasto(ar%r(i))-lasto(ar%r(i)-1)
    end do

    ! to be sure we just depopulate the region
    ! and populate it with the correct orbitals
    tmp = or%name
    call rgn_delete(or)
    or%name = tmp
    or%n = no
    nullify(or%r)
    if ( no > 0 ) then
       allocate(or%r(no))
       call memory('A','I',no,'rgn-list')
       no = 0
       do i = 1 , ar%n
          do cr = lasto(ar%r(i)-1) + 1 , lasto(ar%r(i))
             no = no + 1
             or%r(no) = cr
          end do
       end do
    end if

  end subroutine rgn_Atom2Orb

  subroutine rgn_Orb2Atom(or,na_u,lasto,ar)
    ! The orbital region we wish to convert to an atomic
    ! region
    type(tRgn), intent(in) :: or
    ! The last orbitals of each atom
    integer, intent(in) :: na_u, lasto(0:na_u)
    type(tRgn), intent(inout) :: ar

    ! ** local variables
    integer :: io, ia, na, a
    integer :: ca(na_u)
    character(len=R_NAME_LEN) :: tmp

    if ( or%n == 0 ) then
       call rgn_delete(ar)
       return
    end if

    ! The maximum number of atoms in the orbital region
    ia = 1
    ca(1) = iaorb(or%r(1),lasto)
    do io = 2 , or%n
       a = iaorb(or%r(io),lasto)
       if ( any(a == ca(1:ia)) ) cycle
       ia = ia + 1
       if ( ia > na_u ) call die('Error in program')
       ca(ia) = a
    end do
    na = ia

    ! to be sure we just depopulate the region
    ! and populate it with the correct atoms
    tmp = ar%name
    call rgn_list(ar,na,ca, name = tmp )

  end subroutine rgn_Orb2Atom

  ! Whether any element exists in the other one.
  function rgn_overlaps(r1,r2) result(overlap)
    ! the regions we wish to find the union of
    type(tRgn), intent(in) :: r1, r2
    logical :: overlap

    ! ** local variables
    integer :: i

    overlap = .false.

    if ( r1%n == 0 ) return
    if ( r2%n == 0 ) return

    ! We loop over the long one
    if ( r1%n < r2%n ) then
       
       do i = 1 , r2%n
          if ( in_rgn(r1,r2%r(i)) ) then
             overlap = .true. 
             return
          end if
       end do

    else

       do i = 1 , r1%n
          if ( in_rgn(r2,r1%r(i)) ) then
             overlap = .true. 
             return
          end if
       end do

    end if

  end function rgn_overlaps

  subroutine rgn_print(r, name, seq_max, indent, collapse)
    type(tRgn), intent(in) :: r
    character(len=*), intent(in), optional :: name
    integer, intent(in), optional :: seq_max, indent
    logical, intent(in), optional :: collapse
 
    character(len=5) :: fmt
    integer :: i, j, k, c, lseq_max
    logical :: lcollapse

    lseq_max = 7
    if ( present(seq_max) ) lseq_max = seq_max
    lcollapse = .true.
    if ( present(collapse) ) lcollapse = collapse

    fmt = ' '
    if ( present(indent) ) write(fmt,'(a,i0,a)')'tr',indent,','
    
    if ( present(name) ) then
       write(*,'('//fmt//'a,i0,2a)') trim(name)//' (',r%n,'): ',trim(r%name)
    else
       write(*,'('//fmt//'a,i0,2a)') 'Region (',r%n,'): ',trim(r%name)
    end if
    write(*,'('//fmt//'tr2,a)',advance='no') '['
    if ( .not. lcollapse ) then
       do i = 1 , r%n - 1
          write(*,'(tr1,i0,a)',advance='no') r%r(i),','
          if ( mod(i,lseq_max) == 0 ) then
             write(*,'(/,'//fmt//'tr3)',advance='no')
          end if
       end do
    else

       j = 1
       c = 0
       do while ( j <= r%n ) 

          ! Count the consecutive numbers
          k = 0
          i = j
          do while ( r%r(i) - r%r(j) == k )
             k = k + 1
             i = i + 1
             ! Ensure we do not by pass the regions
             if ( i > r%n ) exit
          end do

          if ( k > 1 ) then
             ! This is a consecutive block
             if ( lseq_max - c < 2 ) then
                write(*,'(/,'//fmt//'tr3)',advance='no') !
                c = 0
             end if
             ! Write out the current line
             write(*,'(tr1,i0,a,i0)',advance='no') &
                  r%r(j),' -- ',r%r(j+k-1)
             i = c + 3
             if ( i >= lseq_max ) then
                c = 0
             else
                c = i
             end if
          else
             ! This is not a consecutive block
             write(*,'(tr1,i0)',advance='no') r%r(j)
             c = c + 1
             if ( mod(c,lseq_max) == 0 ) c = 0
          end if

          ! Write ending comma
          if ( j + k <= r%n ) write(*,'(a)',advance='no') ','

          if ( c == 0 .and. j + k < r%n ) then
             ! Should we start a new line
             write(*,'(/,'//fmt//'tr3)',advance='no')
          end if

          ! Update the currently reached element
          j = j + k

       end do
    end if
    if ( .not. lcollapse ) then
       if ( r%n > 0 ) then
          write(*,'(tr1,i0,a)') r%r(r%n),' ]'
       else
          write(*,'(tr1,a)') ' ]'
       end if
    else
       write(*,'(a)') ' ]'
    end if

  end subroutine rgn_print

  function rgn_push_val(r,val) result(good)
    type(tRgn), intent(inout) :: r
    integer, intent(in) :: val
    logical :: good

    good = ( size(r%r) > r%n ) 
    if ( .not. good ) return

    if ( r%sorted ) then
       r%sorted = r%r(r%n) <= val
    end if

    r%n = r%n + 1
    r%r(r%n) = val

  end function rgn_push_val

  function rgn_push_rgn(r,push) result(good)
    type(tRgn), intent(inout) :: r
    type(tRgn), intent(in) :: push
    integer :: i
    logical :: good

    good = ( size(r%r) >= r%n + push%n ) 
    if ( .not. good ) return
    
    if ( r%sorted .and. push%sorted ) then
       r%sorted = r%r(r%n) <= push%r(1)
    end if
    
    do i = 1 , push%n
       r%r(r%n+i) = push%r(i)
    end do
    r%n = r%n + push%n

  end function rgn_push_rgn

  pure function rgn_sum(r) result(s)
    type(tRgn), intent(in) :: r
    integer :: s

    if ( r%n == 0 ) then
       s = 0
    else
       s = sum(r%r)
    end if
    
  end function rgn_sum


  ! Popping of an index of a region
  function rgn_pop(r,idx,val) result(out)
    type(tRgn), intent(inout) :: r
    integer, intent(in), optional :: val, idx
    integer :: out, i, j

    if ( r%n == 0 ) then
       out = 0
       return
    end if

    if ( present(idx) ) then
       i = idx
    else if ( present(val) ) then
       i = rgn_pivot(r,val)
    else
       i = 1 
    end if

    ! get the value
    out = r%r(i)
    if ( r%n > 1 .and. i < r%n ) then
       ! remove it
       do j = i , r%n - 1
          r%r(j) = r%r(j+1)
       end do
    end if
    r%n = r%n - 1
    
  end function rgn_pop

#ifdef MPI
  subroutine rgn_MPI_union(dit,r)
    use mpi_siesta, only : MPI_AllReduce, MPI_Integer, MPI_Sum
    use mpi_siesta, only : MPI_Bcast
    use mpi_siesta, only : MPI_Recv, MPI_Send, MPI_STATUS_SIZE
    use mpi_siesta, only : MPI_Get_Count

    type(OrbitalDistribution), intent(in) :: dit
    type(tRgn), intent(inout) :: r
    
    ! Our temporary region
    integer :: i, nt, ct, iN, it
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
          do i = 1 , r%n
             rd(i) = r%r(i)
          end do
       end if
       ct = r%n + 1
       do iN = 1 , dist_nodes(dit) - 1
          call MPI_Recv(rd(ct),nt-ct+1,MPI_Integer, &
               iN, 0, comm, MPIstatus, MPIerror)
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
    call rgn_list(r,nt,rd, name=tmp)
    
    ! Clean-up
    deallocate(rd)

  end subroutine rgn_MPI_union

  subroutine rgn_MPI_Bcast(r,Bnode,Comm)
    use mpi_siesta, only : MPI_Bcast, MPI_Comm_World
    use mpi_siesta, only : MPI_Logical, MPI_Integer

    type(tRgn), intent(inout) :: r
    integer, intent(in) :: Bnode
    integer, intent(in), optional :: Comm

    integer :: Node, n, lComm
    integer :: MPIerror

    lComm = MPI_COMM_WORLD
    if ( present(Comm) ) lComm = Comm

    ! Get current node to 
    call MPI_Comm_Rank(lComm,Node,MPIerror)

    n = r%n
    call MPI_Bcast(n,1,MPI_Integer, Bnode, lComm, MPIerror)
    if ( Node /= Bnode ) then
       ! ensures that it is having the correct size
       call rgn_init(r,n)
    end if

    ! If empty-simply return without doing any B-cast
    if ( n == 0 ) return
    
    ! B-cast array
    call MPI_Bcast(r%r(1),r%n, MPI_Integer, Bnode, lComm, MPIerror)
    ! B-cast sorted token
    call MPI_Bcast(r%sorted,1, MPI_Logical, Bnode, lComm, MPIerror)

  end subroutine rgn_MPI_Bcast

#endif
  
end module m_region
