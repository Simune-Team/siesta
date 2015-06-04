!!@LICENSE
!
! ==================================================================
! Allocation, reallocation, and deallocation utility routines
! Written by J.M.Soler. May 2000.
! Re-organized by A. Garcia, June 2015.
! ==================================================================
! SUBROUTINE alloc_default( old, new, restore, &
!                           copy, shrink, imin, routine )
!   Sets defaults for allocation
! INPUT (optional):
!   type(allocDefaults) restore : default settings to be restored
!   logical             copy    : Copy old array to new array?
!   logical             shrink  : Reduce array size?
!   integer             imin    : First index (typically 1 in Fortan,
!                                              0 in C)
!   character(len=*)    routine : Name of calling routine
! OUTPUT (optional):
!   type(allocDefaults) old     : default settings before the call
!   type(allocDefaults) new     : default settings after the call
! BEHAVIOR:
!   All these defaults can be superseeded by optional arguments in
!     each call to re_alloc.
!   Initial default values: copy    = .true.
!                           shrink  = .true.
!                           imin    = 1
!                           routine = 'unknown'
!   If restore is present together with any of copy, shrink, imin, or
!   routine, these are applied AFTER resetting the restore defaults.
! USAGE:
!   In order to restore the allocation defaults possibly set by the
! calling routine, the suggested construction is:
!   use alloc_module
!   type(allocDefaults) oldDefaults
!   call alloc_default( old=oldDefaults, routine=..., &
!                       copy=..., shrink=... )
!   call re_alloc(...)
!   call alloc_default( restore=oldDefaults )
! Notice that, if the restore call is skipped, the new defaults will
! stay in effect until a new call to alloc_dafault is made.
! ==================================================================
! SUBROUTINE re_alloc( array, [i1min,] i1max,
!                      [[i2min,] i2max, [[i3min,] i3max]],
!                      name, routine, copy, shrink )
! INPUT:
!   integer       :: i1min        : Lower bound of first dimension
!                                   If not present, it is fixed by 
!                                   the last call to alloc_default.
!                                   If present and the rank is 2(3),
!                                   then i2min(&i3min) must also be
!                                   present
!   integer       :: i1max        : Upper bound of first dimension
!   integer       :: i2min,i2max  : Bounds of second dimension, if
!                                   applicable
!   integer       :: i3min,i3max  : Bounds of third dimension, if appl.
!
! INPUT (optional):
!   character*(*) :: name         : Actual array name or a label for it
!   character*(*) :: routine      : Name of the calling routine
!                                     or routine section
!   logical       :: copy         : Save (copy) contents of old array 
!                                   to new array?
!   logical       :: shrink       : Reallocate if the new array bounds 
!                                   are contained within the old ones? 
!                                   If not present, copy and/or shrink
!                                      are fixed by the last call to
!                                      alloc_default. 
! INPUT/OUTPUT:
!   TYPE, pointer :: array : Array to be allocated or reallocated.
!                            Implemented types and ranks are:
!                              integer,    rank 1, 2, 3
!                              integer*8,  rank 1
!                              real*4,     rank 1, 2, 3, 4
!                              real*8,     rank 1, 2, 3, 4
!                              complex*16, rank 1, 2
!                              logical,    rank 1, 2, 3
!                              character(len=*), rank 1
! BEHAVIOR:
!   Pointers MUST NOT enter in an undefined state. Before using them
! for the first time, they must be nullified explicitly. Alternatively,
! in f95, they can be initialized as null() upon declaration.
!   If argument array is not associated on input, it is just allocated.
!   If array is associated and has the same bounds (or smaller bonds
! and shrink is false) nothing is done. Thus, it is perfectly safe and
! efficient to call re_alloc repeatedly without deallocating the array.
! However, subroutine dealloc is provided to eliminate large arrays
! when they are not needed.
!   In order to save (copy) the contents of the old array, the new array
! needs to be allocated before deallocating the old one. Thus, if the
! contents are not needed, or if reducing memory is a must, calling
! re_alloc with copy=.false. makes it to deallocate before allocating.
!   The elements that are not copied (because copy=.false. or because
! they are outside the bounds of the input array) return with value
! zero (integer and real), .false. (logical), or blank (character).
!   If imin>imax for any dimension, the array pointer returns
! associated to a zero-size array.
!   Besides allocating or reallocating the array, re_alloc keeps a count
! of memory usage (and prints a report in a file previously fixed by a
! call to alloc_report). Thus, it is not recommended to call re_alloc
! to reallocate an array not allocated by it the first time.
!   If name is not present, the memory count associated to the
! allocation/deallocation is still made, but the allocation report
! will account the array under the generic name 'unknown'.
!   The routine argument is NOT used to classify the counting of
! memory usage. The classification uses only the name argument.
! This is because a pointer may be allocated in one routine and
! deallocated in a different one (i.e. when it is used to return an
! array whose size is not known in advance). However, the routine
! argument is reported when alloc_report=4, and it is also used to
! report in which routine the memory peack occurs. If you want the
! routine name to be used for classification, you should include it
! as part of the name argument, like in name='matInvert '//'aux'.
! ==================================================================---
! SUBROUTINE de_alloc( array, name, routine )
! INPUT (optional):
!   character*(*) :: name    : Actual array name or a label for it
!   character*(*) :: routine : Name of the calling routine
!                                or routine section
! INPUT/OUTPUT:
!   TYPE, pointer :: array : Array be deallocated (same types and
!                            kinds as in re_alloc).
! BEHAVIOR:
!   Besides deallocating the array, re_alloc decreases the count of
! memory usage previously counted by re_alloc. Thus, dealloc should 
! not be called to deallocate an array not allocated by re_alloc.
! Equally, arrays allocated or reallocated by re_alloc should be 
! deallocated by dealloc.
! ==================================================================---
MODULE alloc

  use memory_log, only: alloc_count

  implicit none

PUBLIC ::             &
  alloc_default,      &! Sets allocation defaults
  re_alloc,           &! Allocation/reallocation
  de_alloc,           &! Deallocation
  allocDefaults        ! Derived type to hold allocation defaults

PRIVATE      ! Nothing is declared public beyond this point

integer, parameter :: sp = selected_real_kind(5,10)
integer, parameter :: dp = selected_real_kind(10,100)

! Interfaces to external routines that must be provided
! by the calling program
!
interface
   ! Message and integer code
   ! If 'code' is 0, this is the last call in a series
   ! (see below for usage)
   subroutine alloc_error_report(str,code)
     character(len=*), intent(in) :: str
     integer, intent(in)          :: code
   end subroutine alloc_error_report
end interface

  interface de_alloc
    module procedure &
      dealloc_i1, dealloc_i2, dealloc_i3,             &
      dealloc_E1,                                     &
      dealloc_r1, dealloc_r2, dealloc_r3, dealloc_r4, &
      dealloc_d1, dealloc_d2, dealloc_d3, dealloc_d4, &
      dealloc_z1, dealloc_z2,                         &
      dealloc_l1, dealloc_l2, dealloc_l3,             &
      dealloc_s1
  end interface

  interface re_alloc
    module procedure &
      realloc_i1, realloc_i2, realloc_i3,             &
      realloc_E1,                                     &
      realloc_r1, realloc_r2, realloc_r3, realloc_r4, &
      realloc_d1, realloc_d2, realloc_d3, realloc_d4, &
      realloc_z1, realloc_z2,                         &
      realloc_l1, realloc_l2, realloc_l3,             &
      realloc_s1
!    module procedure & ! AG: Dangerous!!!
!      realloc_i1s, realloc_i2s, realloc_i3s,              &
!      realloc_r1s, realloc_r2s, realloc_r3s, realloc_r4s, &
!      realloc_d1s, realloc_d2s, realloc_d3s, realloc_d4s, &
!      realloc_l1s, realloc_l2s, realloc_l3s
  end interface

  ! Initial default values
  character(len=*), parameter :: &
    DEFAULT_NAME = 'unknown_name'         ! Array name default
  character(len=*), parameter :: &
    DEFAULT_ROUTINE = 'unknown_routine'   ! Routine name default
  integer, save ::               &
    REPORT_LEVEL = 0,            &! Level (detail) of allocation report
    REPORT_UNIT  = 0              ! Output file unit for report
  character(len=50), save ::     &
    REPORT_FILE = 'alloc_report'  ! Output file name for report
  real(dp), save ::              &
    REPORT_THRESHOLD = 0          ! Memory threshold (in bytes) to print
                                  ! the memory use of any given array 

  ! Derived type to hold allocation default options
  type allocDefaults
    private
    logical ::          copy = .true.              ! Copy default
    logical ::          shrink = .true.            ! Shrink default
    integer ::          imin = 1                   ! Imin default
    character(len=32):: routine = DEFAULT_ROUTINE  ! Routine name default
  end type allocDefaults

  ! Object to hold present allocation default options
  type(allocDefaults), save :: DEFAULT

  ! Internal auxiliary type for a binary tree
  type TREE
    character(len=80)   :: name  ! Name of an allocated array
    real(DP)            :: mem   ! Present memory use of the array
    real(DP)            :: max   ! Maximum memory use of the array
    real(DP)            :: peak  ! Memory use of the array during
                                 !   peak of total memory
    type(TREE), pointer :: left  ! Pointer to data of allocated arrays 
                                 !   preceeding in alphabetical order
    type(TREE), pointer :: right ! Pointer to data of allocated arrays 
                                 !   trailing in alphabetical order
  end type TREE

  ! Global variables used to store allocation data
  real(DP),   parameter     :: MBYTE = 1.e6_dp
  type(TREE), pointer, save :: REPORT_TREE
  real(DP),            save :: TOT_MEM  = 0._dp
  real(DP),            save :: PEAK_MEM = 0._dp
  character(len=80),   save :: PEAK_ARRAY = ' '
  character(len=32),   save :: PEAK_ROUTINE = ' '
  integer,             save :: MAX_LEN  = 0
  
  ! Other common variables
  integer :: IERR
  logical :: ASSOCIATED_ARRAY, NEEDS_ALLOC, NEEDS_COPY, NEEDS_DEALLOC

CONTAINS

! ==================================================================

SUBROUTINE alloc_default( old, new, restore,          &
                          routine, copy, shrink, imin )
implicit none
type(allocDefaults), optional, intent(out) :: old, new
type(allocDefaults), optional, intent(in)  :: restore
character(len=*),    optional, intent(in)  :: routine
logical,             optional, intent(in)  :: copy, shrink
integer,             optional, intent(in)  :: imin

if (present(old))     old = DEFAULT
if (present(restore)) DEFAULT = restore
if (present(copy))    DEFAULT%copy   = copy
if (present(shrink))  DEFAULT%shrink = shrink
if (present(imin))    DEFAULT%imin   = imin
if (present(routine)) DEFAULT%routine = routine
if (present(new))     new = DEFAULT

END SUBROUTINE alloc_default

! ==================================================================
! Integer array reallocs
! ==================================================================

SUBROUTINE realloc_i1( array, i1min, i1max, &
                       name, routine, copy, shrink )
! Arguments
implicit none
integer, dimension(:),      pointer    :: array
integer,                    intent(in) :: i1min
integer,                    intent(in) :: i1max
character(len=*), optional, intent(in) :: name
character(len=*), optional, intent(in) :: routine
logical,          optional, intent(in) :: copy
logical,          optional, intent(in) :: shrink

! Internal variables and arrays
character, parameter           :: type='I'
integer, parameter             :: rank=1
integer, dimension(:), pointer :: old_array
integer, dimension(2,rank)     :: b, c, new_bounds, old_bounds

! Get old array bounds
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array          ! Keep pointer to old array
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if

! Copy new requested array bounds
new_bounds(1,:) = (/ i1min /)
new_bounds(2,:) = (/ i1max /)

! Find if it is a new allocation or a true reallocation,
! and if the contents need to be copied (saved)
! Argument b returns common bounds
! Options routine also reads common variable ASSOCIATED_ARRAY,
! and it sets NEEDS_ALLOC, NEEDS_DEALLOC, and NEEDS_COPY
call options( b, c, old_bounds, new_bounds, copy, shrink )
! Deallocate old space
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array), type, name, routine )
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if

! Allocate new space
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1)), stat=IERR )
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = 0
end if

! Copy contents and deallocate old space
if (NEEDS_COPY) then
  array(c(1,1):c(2,1)) = old_array(c(1,1):c(2,1))
  call alloc_count( -size(old_array), type, name, routine )
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_i1

! ==================================================================
SUBROUTINE realloc_i2( array, i1min,i1max, i2min,i2max,       &
                       name, routine, copy, shrink )
implicit none
character, parameter                   :: type='I'
integer, parameter                     :: rank=2
integer, dimension(:,:),    pointer    :: array, old_array
integer,                    intent(in) :: i1min, i1max, i2min, i2max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min, i2min /)
new_bounds(2,:) = (/ i1max, i2max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1),b(1,2):b(2,2)), stat=IERR )
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = 0
end if
if (NEEDS_COPY) then
      array(c(1,1):c(2,1),c(1,2):c(2,2)) =  &
  old_array(c(1,1):c(2,1),c(1,2):c(2,2))
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_i2
! ==================================================================

SUBROUTINE realloc_i3( array, i1min,i1max, i2min,i2max, i3min,i3max, &
                       name, routine, copy, shrink )
implicit none
character, parameter                   :: type='I'
integer, parameter                     :: rank=3
integer, dimension(:,:,:),  pointer    :: array, old_array
integer,                    intent(in) :: i1min,i1max, i2min,i2max, &
                                          i3min,i3max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array 
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min, i2min, i3min /)
new_bounds(2,:) = (/ i1max, i2max, i3max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1),b(1,2):b(2,2),b(1,3):b(2,3)),stat=IERR)
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = 0
end if
if (NEEDS_COPY) then
      array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3)) =  &
  old_array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3))
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_i3
! ==================================================================
SUBROUTINE realloc_E1( array, i1min, i1max, &
                       name, routine, copy, shrink )
! Arguments
implicit none
integer*8, dimension(:),    pointer    :: array
integer,                    intent(in) :: i1min
integer,                    intent(in) :: i1max
character(len=*), optional, intent(in) :: name
character(len=*), optional, intent(in) :: routine
logical,          optional, intent(in) :: copy
logical,          optional, intent(in) :: shrink

! Internal variables and arrays
character, parameter                   :: type='I'
integer, parameter                     :: rank=1
integer*8, dimension(:), pointer       :: old_array
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds

! Get old array bounds
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array          ! Keep pointer to old array
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if

! Copy new requested array bounds
new_bounds(1,:) = (/ i1min /)
new_bounds(2,:) = (/ i1max /)

! Find if it is a new allocation or a true reallocation,
! and if the contents need to be copied (saved)
! Argument b returns common bounds
! Options routine also reads common variable ASSOCIATED_ARRAY,
! and it sets NEEDS_ALLOC, NEEDS_DEALLOC, and NEEDS_COPY
call options( b, c, old_bounds, new_bounds, copy, shrink )

! Deallocate old space
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if

! Allocate new space
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1)), stat=IERR )
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = 0
end if

! Copy contents and deallocate old space
if (NEEDS_COPY) then
  array(c(1,1):c(2,1)) = old_array(c(1,1):c(2,1))
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if

END SUBROUTINE realloc_E1

! ==================================================================
! Single precision real array reallocs
! ==================================================================
SUBROUTINE realloc_r1( array, i1min, i1max,        &
                       name, routine, copy, shrink )
implicit none
character, parameter                   :: type='R'
integer, parameter                     :: rank=1
real(SP), dimension(:),     pointer    :: array, old_array
integer,                    intent(in) :: i1min, i1max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array 
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min /)
new_bounds(2,:) = (/ i1max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1)), stat=IERR )
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = 0._sp
end if
if (NEEDS_COPY) then
  array(c(1,1):c(2,1)) = old_array(c(1,1):c(2,1))
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_r1
! ==================================================================
SUBROUTINE realloc_r2( array, i1min,i1max, i2min,i2max, &
                       name, routine, copy, shrink )
implicit none
character, parameter             :: type='R'
integer, parameter               :: rank=2
real(SP), dimension(:,:),   pointer    :: array, old_array
integer,                    intent(in) :: i1min, i1max, i2min, i2max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array 
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min, i2min /)
new_bounds(2,:) = (/ i1max, i2max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1),b(1,2):b(2,2)), stat=IERR )
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = 0._sp
end if
if (NEEDS_COPY) then
      array(c(1,1):c(2,1),c(1,2):c(2,2)) =  &
  old_array(c(1,1):c(2,1),c(1,2):c(2,2))
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_r2
! ==================================================================
SUBROUTINE realloc_r3( array, i1min,i1max, i2min,i2max, i3min,i3max, &
                       name, routine, copy, shrink )
implicit none
character, parameter                   :: type='R'
integer, parameter                     :: rank=3
real(SP), dimension(:,:,:), pointer    :: array, old_array
integer,                    intent(in) :: i1min,i1max, i2min,i2max, &
                                          i3min,i3max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array          ! Keep pointer to old array
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min, i2min, i3min /)
new_bounds(2,:) = (/ i1max, i2max, i3max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1),b(1,2):b(2,2),b(1,3):b(2,3)),stat=IERR)
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = 0._sp
end if
if (NEEDS_COPY) then
      array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3)) =  &
  old_array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3))
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_r3
! ==================================================================
SUBROUTINE realloc_r4( array, i1min,i1max, i2min,i2max, &
                              i3min,i3max, i4min,i4max, &
                       name, routine, copy, shrink )
implicit none
character, parameter                   :: type='R'
integer, parameter                     :: rank=4
real(SP), dimension(:,:,:,:), pointer  :: array, old_array
integer,                    intent(in) :: i1min,i1max, i2min,i2max, &
                                          i3min,i3max, i4min,i4max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array 
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min, i2min, i3min, i4min /)
new_bounds(2,:) = (/ i1max, i2max, i3max, i4max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1),b(1,2):b(2,2), &
                 b(1,3):b(2,3),b(1,4):b(2,4)),stat=IERR)
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = 0._sp
end if
if (NEEDS_COPY) then
      array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3),c(1,4):c(2,4))= &
  old_array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3),c(1,4):c(2,4))
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_r4
! ==================================================================
! Double precision real array reallocs
! ==================================================================
SUBROUTINE realloc_d1( array, i1min, i1max,        &
                       name, routine, copy, shrink )
implicit none
character, parameter                   :: type='D'
integer, parameter                     :: rank=1
real(DP), dimension(:),     pointer    :: array, old_array
integer,                    intent(in) :: i1min, i1max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min /)
new_bounds(2,:) = (/ i1max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1)), stat=IERR )
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = 0._dp
end if
if (NEEDS_COPY) then
  array(c(1,1):c(2,1)) = old_array(c(1,1):c(2,1))
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_d1
! ==================================================================
SUBROUTINE realloc_d2( array, i1min,i1max, i2min,i2max, &
                       name, routine, copy, shrink )
implicit none
character, parameter                   :: type='D'
integer, parameter                     :: rank=2
real(DP), dimension(:,:),   pointer    :: array, old_array
integer,                    intent(in) :: i1min, i1max, i2min, i2max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
integer                                :: i1, i2
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array 
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min, i2min /)
new_bounds(2,:) = (/ i1max, i2max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1),b(1,2):b(2,2)), stat=IERR )
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = 0._dp
end if
if (NEEDS_COPY) then
!      array(c(1,1):c(2,1),c(1,2):c(2,2)) =  &
!  old_array(c(1,1):c(2,1),c(1,2):c(2,2))
  do i2 = c(1,2),c(2,2)
  do i1 = c(1,1),c(2,1)
    array(i1,i2) = old_array(i1,i2)
  end do
  end do
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_d2
! ==================================================================
SUBROUTINE realloc_d3( array, i1min,i1max, i2min,i2max, i3min,i3max, &
                       name, routine, copy, shrink )
implicit none
character, parameter                   :: type='D'
integer, parameter                     :: rank=3
real(DP), dimension(:,:,:), pointer    :: array, old_array
integer,                    intent(in) :: i1min,i1max, i2min,i2max, &
                                          i3min,i3max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
integer                                :: i1, i2, i3
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array 
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min, i2min, i3min /)
new_bounds(2,:) = (/ i1max, i2max, i3max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1),b(1,2):b(2,2),b(1,3):b(2,3)),stat=IERR)
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = 0._dp
end if
if (NEEDS_COPY) then
!      array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3)) =  &
!  old_array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3))
  do i3 = c(1,3),c(2,3)
  do i2 = c(1,2),c(2,2)
  do i1 = c(1,1),c(2,1)
    array(i1,i2,i3) = old_array(i1,i2,i3)
  end do
  end do
  end do
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_d3
! ==================================================================
SUBROUTINE realloc_d4( array, i1min,i1max, i2min,i2max, &
                              i3min,i3max, i4min,i4max, &
                       name, routine, copy, shrink )
implicit none
character, parameter                   :: type='D'
integer, parameter                     :: rank=4
real(DP), dimension(:,:,:,:), pointer  :: array, old_array
integer,                    intent(in) :: i1min,i1max, i2min,i2max, &
                                          i3min,i3max, i4min,i4max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array 
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min, i2min, i3min, i4min /)
new_bounds(2,:) = (/ i1max, i2max, i3max, i4max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1),b(1,2):b(2,2), &
                 b(1,3):b(2,3),b(1,4):b(2,4)),stat=IERR)
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = 0._dp
end if
if (NEEDS_COPY) then
      array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3),c(1,4):c(2,4))= &
  old_array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3),c(1,4):c(2,4))
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_d4
! ==================================================================
! Double precision complex array reallocs
! ==================================================================
SUBROUTINE realloc_z1( array, i1min, i1max,        &
                       name, routine, copy, shrink )
implicit none
character, parameter                   :: type='D'
integer, parameter                     :: rank=1
complex(DP), dimension(:),  pointer    :: array, old_array
integer,                    intent(in) :: i1min, i1max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min /)
new_bounds(2,:) = (/ i1max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -2*size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1)), stat=IERR )
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( 2*size(array), type, name, routine )
  array = 0._dp
end if
if (NEEDS_COPY) then
  array(c(1,1):c(2,1)) = old_array(c(1,1):c(2,1))
  call alloc_count( -2*size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_z1
! ==================================================================
SUBROUTINE realloc_z2( array, i1min,i1max, i2min,i2max, &
                       name, routine, copy, shrink )
implicit none
character, parameter                   :: type='D'
integer, parameter                     :: rank=2
complex(DP), dimension(:,:),  pointer  :: array, old_array
integer,                    intent(in) :: i1min, i1max, i2min, i2max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
integer                                :: i1, i2
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array 
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min, i2min /)
new_bounds(2,:) = (/ i1max, i2max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -2*size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1),b(1,2):b(2,2)), stat=IERR )
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( 2*size(array), type, name, routine )
  array = 0._dp
end if
if (NEEDS_COPY) then
!      array(c(1,1):c(2,1),c(1,2):c(2,2)) =  &
!  old_array(c(1,1):c(2,1),c(1,2):c(2,2))
  do i2 = c(1,2),c(2,2)
  do i1 = c(1,1),c(2,1)
    array(i1,i2) = old_array(i1,i2)
  end do
  end do
  call alloc_count( -2*size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_z2
! ==================================================================
! Logical array reallocs
! ==================================================================
SUBROUTINE realloc_l1( array, i1min,i1max,  &
                       name, routine, copy, shrink )
implicit none
character, parameter                   :: type='L'
integer, parameter                     :: rank=1
logical, dimension(:),      pointer    :: array, old_array
integer,                    intent(in) :: i1min,i1max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array 
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min /)
new_bounds(2,:) = (/ i1max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1)), stat=IERR )
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = .false.
end if
if (NEEDS_COPY) then
  array(c(1,1):c(2,1)) = old_array(c(1,1):c(2,1))
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_l1
! ==================================================================
SUBROUTINE realloc_l2( array, i1min,i1max, i2min,i2max, &
                       name, routine, copy, shrink )
implicit none
character, parameter                   :: type='L'
integer, parameter                     :: rank=2
logical, dimension(:,:),    pointer    :: array, old_array
integer,                    intent(in) :: i1min,i1max, i2min,i2max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array 
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min, i2min /)
new_bounds(2,:) = (/ i1max, i2max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1),b(1,2):b(2,2)), stat=IERR )
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = .false.
end if
if (NEEDS_COPY) then
      array(c(1,1):c(2,1),c(1,2):c(2,2)) = &
  old_array(c(1,1):c(2,1),c(1,2):c(2,2))
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_l2
! ==================================================================
SUBROUTINE realloc_l3( array, i1min,i1max, i2min,i2max, i3min,i3max, &
                       name, routine, copy, shrink )
implicit none
character, parameter                   :: type='L'
integer, parameter                     :: rank=3
logical, dimension(:,:,:),  pointer    :: array, old_array
integer,                    intent(in) :: i1min,i1max, i2min,i2max, &
                                          i3min,i3max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array 
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min, i2min, i3min /)
new_bounds(2,:) = (/ i1max, i2max, i3max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1),b(1,2):b(2,2),b(1,3):b(2,3)),stat=IERR)
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = .false.
end if
if (NEEDS_COPY) then
      array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3)) = &
  old_array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3))
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_l3
! ==================================================================
! Realloc routines with assumed lower bound = 1
!AG: Extremely dangerous -- do not use.
! ==================================================================
SUBROUTINE realloc_i1s( array, i1max, &
                        name, routine, copy, shrink )
! Arguments
implicit none
integer, dimension(:),      pointer    :: array
integer,                    intent(in) :: i1max
character(len=*), optional, intent(in) :: name
character(len=*), optional, intent(in) :: routine
logical,          optional, intent(in) :: copy
logical,          optional, intent(in) :: shrink

call realloc_i1( array, DEFAULT%imin, i1max, &
                 name, routine, copy, shrink )

END SUBROUTINE realloc_i1s
! ==================================================================
SUBROUTINE realloc_i2s( array, i1max, i2max,  &
                        name, routine, copy, shrink )
implicit none
integer, dimension(:,:),    pointer    :: array
integer,                    intent(in) :: i1max, i2max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
call realloc_i2( array, DEFAULT%imin, i1max, DEFAULT%imin, i2max, &
                 name, routine, copy, shrink )
END SUBROUTINE realloc_i2s
! ==================================================================
SUBROUTINE realloc_i3s( array, i1max, i2max, i3max,  &
                        name, routine, copy, shrink )
implicit none
integer, dimension(:,:,:),  pointer    :: array
integer,                    intent(in) :: i1max, i2max, i3max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
call realloc_i3( array, DEFAULT%imin, i1max, DEFAULT%imin, i2max, &
                 DEFAULT%imin, i3max,                             &
                 name, routine, copy, shrink )
END SUBROUTINE realloc_i3s
! ==================================================================
SUBROUTINE realloc_r1s( array, i1max, &
                        name, routine, copy, shrink )
implicit none
real(SP), dimension(:),     pointer    :: array
integer,                    intent(in) :: i1max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
call realloc_r1( array, DEFAULT%imin, i1max, &
                 name, routine, copy, shrink )
END SUBROUTINE realloc_r1s
! ==================================================================
SUBROUTINE realloc_r2s( array, i1max, i2max, &
                        name, routine, copy, shrink )
implicit none
real(SP), dimension(:,:),   pointer    :: array
integer,                    intent(in) :: i1max, i2max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
call realloc_r2( array, DEFAULT%imin, i1max, DEFAULT%imin, i2max, &
                 name, routine, copy, shrink )
END SUBROUTINE realloc_r2s
! ==================================================================
SUBROUTINE realloc_r3s( array, i1max, i2max, i3max, &
                        name, routine, copy, shrink )
implicit none
real(SP), dimension(:,:,:), pointer    :: array
integer,                    intent(in) :: i1max, i2max, i3max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
call realloc_r3( array, DEFAULT%imin, i1max, DEFAULT%imin, i2max, &
                 DEFAULT%imin, i3max,                             &
                 name, routine, copy, shrink )
END SUBROUTINE realloc_r3s
! ==================================================================
SUBROUTINE realloc_r4s( array, i1max, i2max, i3max, i4max, &
                        name, routine, copy, shrink )
implicit none
real(SP), dimension(:,:,:,:), pointer  :: array
integer,                    intent(in) :: i1max, i2max, i3max, i4max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
call realloc_r4( array, DEFAULT%imin, i1max, DEFAULT%imin, i2max, &
                        DEFAULT%imin, i3max, DEFAULT%imin, i4max, &
                 name, routine, copy, shrink )
END SUBROUTINE realloc_r4s
! ==================================================================
SUBROUTINE realloc_d1s( array, i1max, &
                        name, routine, copy, shrink )
implicit none
real(DP), dimension(:),     pointer    :: array
integer,                    intent(in) :: i1max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
call realloc_d1( array, DEFAULT%imin, i1max, &
                 name, routine, copy, shrink )
END SUBROUTINE realloc_d1s
! ==================================================================
SUBROUTINE realloc_d2s( array, i1max, i2max, &
                        name, routine, copy, shrink )
implicit none
real(DP), dimension(:,:),   pointer    :: array
integer,                    intent(in) :: i1max, i2max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
call realloc_d2( array, DEFAULT%imin, i1max, DEFAULT%imin, i2max, &
                 name, routine, copy, shrink )
END SUBROUTINE realloc_d2s
! ==================================================================
SUBROUTINE realloc_d3s( array, i1max, i2max, i3max, &
                        name, routine, copy, shrink )
implicit none
real(DP), dimension(:,:,:), pointer    :: array
integer,                    intent(in) :: i1max, i2max, i3max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
call realloc_d3( array, DEFAULT%imin, i1max, DEFAULT%imin, i2max, &
                 DEFAULT%imin, i3max,                             &
                 name, routine, copy, shrink )
END SUBROUTINE realloc_d3s
! ==================================================================
SUBROUTINE realloc_d4s( array, i1max, i2max, i3max, i4max, &
                        name, routine, copy, shrink )
implicit none
real(DP), dimension(:,:,:,:), pointer  :: array
integer,                    intent(in) :: i1max, i2max, i3max, i4max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
call realloc_d4( array, DEFAULT%imin, i1max, DEFAULT%imin, i2max, &
                        DEFAULT%imin, i3max, DEFAULT%imin, i4max, &
                 name, routine, copy, shrink )
END SUBROUTINE realloc_d4s
! ==================================================================
SUBROUTINE realloc_l1s( array, i1max, &
                        name, routine, copy, shrink )
implicit none
logical, dimension(:),      pointer    :: array
integer,                    intent(in) :: i1max
character(len=*), optional, intent(in) :: name
character(len=*), optional, intent(in) :: routine
logical,          optional, intent(in) :: copy
logical,          optional, intent(in) :: shrink
call realloc_l1( array, DEFAULT%imin, i1max, &
                 name, routine, copy, shrink )
END SUBROUTINE realloc_l1s
! ==================================================================
SUBROUTINE realloc_l2s( array, i1max, i2max, &
                        name, routine, copy, shrink )
implicit none
logical, dimension(:,:),    pointer    :: array
integer,                    intent(in) :: i1max, i2max
character(len=*), optional, intent(in) :: name
character(len=*), optional, intent(in) :: routine
logical,          optional, intent(in) :: copy
logical,          optional, intent(in) :: shrink
call realloc_l2( array, DEFAULT%imin, i1max, DEFAULT%imin, i2max, &
                 name, routine, copy, shrink )
END SUBROUTINE realloc_l2s
! ==================================================================
SUBROUTINE realloc_l3s( array, i1max, i2max, i3max, &
                        name, routine, copy, shrink )
implicit none
logical, dimension(:,:,:),  pointer    :: array
integer,                    intent(in) :: i1max, i2max, i3max
character(len=*), optional, intent(in) :: name
character(len=*), optional, intent(in) :: routine
logical,          optional, intent(in) :: copy
logical,          optional, intent(in) :: shrink
call realloc_l3( array, DEFAULT%imin, i1max, DEFAULT%imin, i2max, &
                 DEFAULT%imin, i3max, name, routine, copy, shrink )
END SUBROUTINE realloc_l3s
! ==================================================================
! Character vector realloc
! ==================================================================
SUBROUTINE realloc_s1( array, i1min, i1max, &
                       name, routine, copy, shrink )
! Arguments
implicit none
character(len=*), dimension(:),      pointer    :: array
integer,                    intent(in) :: i1min
integer,                    intent(in) :: i1max
character(len=*), optional, intent(in) :: name
character(len=*), optional, intent(in) :: routine
logical,          optional, intent(in) :: copy
logical,          optional, intent(in) :: shrink

! Internal variables and arrays
character, parameter           :: type='S'
integer, parameter             :: rank=1
character(len=len(array)), dimension(:), pointer :: old_array
integer, dimension(2,rank)     :: b, c, new_bounds, old_bounds

! Get old array bounds
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array          ! Keep pointer to old array
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if

! Copy new requested array bounds
new_bounds(1,:) = (/ i1min /)
new_bounds(2,:) = (/ i1max /)

! Find if it is a new allocation or a true reallocation,
! and if the contents need to be copied (saved)
! Argument b returns common bounds
! Options routine also reads common variable ASSOCIATED_ARRAY,
! and it sets NEEDS_ALLOC, NEEDS_DEALLOC, and NEEDS_COPY
call options( b, c, old_bounds, new_bounds, copy, shrink )

! Deallocate old space
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array)*len(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if

! Allocate new space
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1)), stat=IERR )
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array)*len(array), type, name, routine )
  array = ''
end if

! Copy contents and deallocate old space
if (NEEDS_COPY) then
  array(c(1,1):c(2,1)) = old_array(c(1,1):c(2,1))
  call alloc_count( -size(old_array)*len(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if

END SUBROUTINE realloc_s1
! ==================================================================
! Dealloc routines
! ==================================================================
SUBROUTINE dealloc_i1( array, name, routine )

! Arguments
implicit none
integer, dimension(:),      pointer    :: array
character(len=*), optional, intent(in) :: name
character(len=*), optional, intent(in) :: routine

if (associated(array)) then
  call alloc_count( -size(array), 'I', name, routine )
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if

END SUBROUTINE dealloc_i1

! ==================================================================
SUBROUTINE dealloc_i2( array, name, routine )
implicit none
integer, dimension(:,:),    pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'I', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )

end if
END SUBROUTINE dealloc_i2
! ==================================================================
SUBROUTINE dealloc_i3( array, name, routine )
implicit none
integer, dimension(:,:,:),  pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'I', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )

end if
END SUBROUTINE dealloc_i3
! ==================================================================
SUBROUTINE dealloc_E1( array, name, routine )

! Arguments
implicit none
integer*8, dimension(:),      pointer    :: array
character(len=*), optional, intent(in) :: name
character(len=*), optional, intent(in) :: routine

if (associated(array)) then
  call alloc_count( -size(array), 'I', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )

end if

END SUBROUTINE dealloc_E1

! ==================================================================
SUBROUTINE dealloc_r1( array, name, routine )
implicit none
real(SP), dimension(:),     pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'R', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_r1
! ==================================================================
SUBROUTINE dealloc_r2( array, name, routine )
implicit none
real(SP), dimension(:,:),   pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'R', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_r2
! ==================================================================
SUBROUTINE dealloc_r3( array, name, routine )
implicit none
real(SP), dimension(:,:,:), pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'R', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_r3
! ==================================================================
SUBROUTINE dealloc_r4( array, name, routine )
implicit none
real(SP), dimension(:,:,:,:), pointer  :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'R', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_r4
! ==================================================================
SUBROUTINE dealloc_d1( array, name, routine )
implicit none
real(DP), dimension(:),     pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'D', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_d1
! ==================================================================
SUBROUTINE dealloc_d2( array, name, routine )
implicit none
real(DP), dimension(:,:),   pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'D', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_d2
! ==================================================================
SUBROUTINE dealloc_d3( array, name, routine )
implicit none
real(DP), dimension(:,:,:), pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'D', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_d3
! ==================================================================
SUBROUTINE dealloc_d4( array, name, routine )
implicit none
real(DP), dimension(:,:,:,:), pointer  :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'D', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_d4
! ==================================================================
! COMPLEX versions
!
SUBROUTINE dealloc_z1( array, name, routine )
implicit none
complex(DP), dimension(:),   pointer   :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -2*size(array), 'D', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_z1
! ==================================================================
SUBROUTINE dealloc_z2( array, name, routine )
implicit none
complex(DP), dimension(:,:),  pointer  :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -2*size(array), 'D', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_z2
! ==================================================================
SUBROUTINE dealloc_l1( array, name, routine )
implicit none
logical, dimension(:),      pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'L', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_l1
! ==================================================================
SUBROUTINE dealloc_l2( array, name, routine )
implicit none
logical, dimension(:,:),    pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'L', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_l2
! ==================================================================
SUBROUTINE dealloc_l3( array, name, routine )
implicit none
logical, dimension(:,:,:),  pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'L', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_l3
! ==================================================================
SUBROUTINE dealloc_s1( array, name, routine )
implicit none
character(len=*), dimension(:), pointer :: array
character(len=*), optional, intent(in)  :: name, routine
if (associated(array)) then
  call alloc_count( -size(array)*len(array), 'S', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_s1

! ==================================================================
! Internal subroutines
! ==================================================================

SUBROUTINE options( final_bounds, common_bounds, &
                    old_bounds, new_bounds, copy, shrink )
! Arguments
integer, dimension(:,:), intent(out) :: final_bounds
integer, dimension(:,:), intent(out) :: common_bounds
integer, dimension(:,:), intent(in)  :: old_bounds
integer, dimension(:,:), intent(in)  :: new_bounds
logical,       optional, intent(in)  :: copy
logical,       optional, intent(in)  :: shrink

! Internal variables and arrays
logical want_shrink


!! AG*****
!  It might be worthwhile to check whether the user
!  atttemps to use bounds which do not make sense,
!  such as zero, or with upper<lower...
!!***

! Find if it is a new allocation or a true reallocation,
! and if the contents need to be copied (saved)
if (ASSOCIATED_ARRAY) then

  ! Check if array bounds have changed
  if ( all(new_bounds==old_bounds) ) then
    ! Old and new arrays are equal. Nothing needs to be done
    NEEDS_ALLOC   = .false. 
    NEEDS_DEALLOC = .false.
    NEEDS_COPY    = .false.
  else 

    ! Want to shrink?
    if (present(shrink)) then
      want_shrink = shrink
    else
      want_shrink = DEFAULT%shrink
    end if

    if (.not. want_shrink  &
        .and. all(new_bounds(1,:)>=old_bounds(1,:)) &
        .and. all(new_bounds(2,:)<=old_bounds(2,:)) ) then
      ! Old array is already fine. Nothing needs to be done
      NEEDS_ALLOC   = .false. 
      NEEDS_DEALLOC = .false.
      NEEDS_COPY    = .false.
    else
      ! Old array needs to be substituted by a new array
      NEEDS_ALLOC   = .true.
      NEEDS_DEALLOC = .true.
      if (present(copy)) then
        NEEDS_COPY = copy
      else
        NEEDS_COPY = DEFAULT%copy
      end if

      ! Ensure that bounds shrink only if desired
      if (want_shrink) then
        final_bounds(1,:) = new_bounds(1,:)
        final_bounds(2,:) = new_bounds(2,:)
      else
        final_bounds(1,:) = min( old_bounds(1,:), new_bounds(1,:) )
        final_bounds(2,:) = max( old_bounds(2,:), new_bounds(2,:) )
      end if

      ! Find common section of old and new arrays
      common_bounds(1,:) = max( old_bounds(1,:), final_bounds(1,:) )
      common_bounds(2,:) = min( old_bounds(2,:), final_bounds(2,:) )
    end if

  end if

else
  ! Old array does not exist. Allocate new one
  NEEDS_ALLOC   = .true. 
  NEEDS_DEALLOC = .false.
  NEEDS_COPY    = .false.
  final_bounds(1,:) = new_bounds(1,:)
  final_bounds(2,:) = new_bounds(2,:)
end if

END SUBROUTINE options

! ==================================================================

SUBROUTINE alloc_err( ierr, name, routine, bounds )
implicit none

integer,                    intent(in) :: ierr
character(len=*), optional, intent(in) :: name
character(len=*), optional, intent(in) :: routine
integer, dimension(:,:), optional, intent(in) :: bounds

integer i
character(len=128) :: msg

if (ierr/=0) then
   write(msg,*) 'alloc_err: allocate status error', ierr
   call alloc_error_report(trim(msg),1)
  if (present(name).and.present(routine)) then
     write(msg,*) 'alloc_err: array ', name, &
               ' requested by ', routine
     call alloc_error_report(trim(msg),2)
  elseif (present(name)) then
     write(msg,*) 'alloc_err: array ', name, &
               ' requested by unknown'
     call alloc_error_report(trim(msg),3)
  elseif (present(routine)) then
     write(msg,*) 'alloc_err: array unknown', &
               ' requested by ', routine
     call alloc_error_report(trim(msg),4)
  endif
  if (present(bounds)) then
     write(msg,'(a,i3,2i10)') ('alloc_err: dim, lbound, ubound:',  &
          i,bounds(1,i),bounds(2,i),                         &
          i=1,size(bounds,dim=2))            
     call alloc_error_report(trim(msg),5)
  endif
end if
call alloc_error_report("alloc_err: end of error report",0)

END SUBROUTINE alloc_err

! ==================================================================

END MODULE alloc
