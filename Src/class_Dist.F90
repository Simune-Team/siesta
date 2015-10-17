!

  module class_Dist
#ifdef MPI
    use mpi
#endif

    use class_BlockCyclicDist
    use class_PEXSIDist

  implicit none

  character(len=*), parameter :: mod_name=__FILE__

#ifndef MPI
  ! To allow compilation in serial mode
  integer, parameter, private :: MPI_COMM_NULL = -huge(1)
  integer, parameter, private :: MPI_UNDEFINED = -huge(1)
#endif

!---------------------------------------------
  type, private :: distKind
     private
     integer :: i
  end type distKind

  type(distKind), parameter, public          :: TYPE_NULL = distKind(-1)
  type(distKind), parameter, public          :: TYPE_BLOCK_CYCLIC = distKind(1)
  type(distKind), parameter, public          :: TYPE_PEXSI = distKind(2)

  public operator(==)
  interface operator(==)
     module procedure isequal_
  end interface
!---------------------------------------------

  type Dist_
     integer              :: refCount = 0
     character(len=36)    :: id
     !------------------------
     character(len=256)   :: name = "null Dist"
     !------------------------
     type(distKind)        :: dist_type
     type(BlockCyclicDist) :: bdist
     type(PEXSIDist)       :: pdist
     !------------------------------------------
  end type Dist_

  type Dist
     type(Dist_), pointer :: data => null()
  end type Dist
  
  type(Dist_), pointer :: obj => null()

  public :: newDistribution, distType
  public :: ref_comm, ranks_in_ref_comm
  public :: num_local_elements, node_handling_element
  public :: index_local_to_global, index_global_to_local

  interface newDistribution
     module procedure newDistribution_
  end interface
  interface distType
     module procedure dist_type
  end interface

  interface ref_comm
     module procedure ref_comm_
  end interface ref_comm

  interface ranks_in_ref_comm
     module procedure ranks_in_ref_comm_
  end interface

  interface num_local_elements
     module procedure num_local_elements_
  end interface
  interface index_local_to_global
     module procedure index_local_to_global_
  end interface
  interface index_global_to_local
     module procedure index_global_to_local_
  end interface
  interface node_handling_element
     module procedure node_handling_element_
  end interface

!====================================    
#define TYPE_NAME Dist
#include "basic_type.inc"
!====================================    

     subroutine delete_Data(spdata)
      type(Dist_) :: spdata
        if (spdata%dist_type == TYPE_BLOCK_CYCLIC) then
           call delete(spdata%bdist)
        else if (spdata%dist_type == TYPE_PEXSI) then
           call delete(spdata%pdist)
        else if (spdata%dist_type == TYPE_NULL) then
           ! do nothing
        else
           ! do nothing
        end if

     end subroutine delete_Data

  function dist_type(this) result(dtype)
    type(Dist)  :: this
    type(distKind) :: dtype
  
    dtype = obj%dist_type

  end function dist_type

  function isequal_ (a,b) result(iseq)
    type(distKind), intent(in) :: a, b
    logical                 :: iseq
    
    iseq = (a%i == b%i)
  end function isequal_
!
  subroutine newDistribution_(this,Ref_Comm,Ranks_in_Ref_Comm, &
       dist_type,Blocksize,name)
     !........................................
     ! Constructor
     !........................................
     type (Dist), intent(inout) :: this
     integer, intent(in)                       :: ref_comm
     integer, intent(in)                       :: ranks_in_ref_comm(:)
     type(distKind), intent(in)                :: dist_type
     integer, intent(in)                       :: Blocksize
     character(len=*), intent(in), optional    :: name

     integer :: error

     call init(this)

     if (dist_type == TYPE_BLOCK_CYCLIC) then
        call newDistribution(obj%bdist,ref_comm,ranks_in_ref_comm,&
             BlockSize,name)
     else if (dist_type == TYPE_PEXSI) then
        call newDistribution(obj%pdist,ref_comm,ranks_in_ref_comm,&
             BlockSize,name)
     else if (dist_type == TYPE_NULL) then
        ! do nothing
     else
        ! do nothing
     end if

     obj%dist_type = dist_type

     if (present(name)) then
        obj%name = trim(name)
     else
        obj%name = "(Abstract Distribution from BS and Group)"
     endif
     call tag_new_object(this)

   end subroutine newDistribution_

!-----------------------------------------------------------
   function ranks_in_ref_comm_(this) result(ranks)
     type(Dist), intent(in)  :: this
     integer, allocatable               :: ranks(:)

     if (dist_type(this) == TYPE_BLOCK_CYCLIC) then
        ranks =  ranks_in_ref_comm(obj%bdist)
     else if (dist_type(this) == TYPE_PEXSI) then
        ranks =  ranks_in_ref_comm(obj%pdist)
     else if (dist_type(this) == TYPE_NULL) then
        allocate(ranks(0))
     else
        allocate(ranks(0))
     end if

   end function ranks_in_ref_comm_
!-----------------------------------------------------------
   function ref_comm_(this) result(comm)
     type(Dist), intent(in)  :: this
     integer                 :: comm

     if (dist_type(this) == TYPE_BLOCK_CYCLIC) then
        comm =  ref_comm(obj%bdist)
     else if (dist_type(this) == TYPE_PEXSI) then
        comm =  ref_comm(obj%pdist)
     else if (dist_type(this) == TYPE_NULL) then
        comm = MPI_COMM_NULL
     else
        comm = MPI_COMM_NULL
     end if

   end function ref_comm_
!-----------------------------------------------------------
   function num_local_elements_(this,nels,Node) result(nl)
     type(Dist), intent(in)  :: this
     integer, intent(in)                    :: nels
     integer, intent(in)                    :: Node
     integer                                :: nl

     if (dist_type(this) == TYPE_BLOCK_CYCLIC) then
        nl =  num_local_elements(obj%bdist,nels,Node)
     else if (dist_type(this) == TYPE_PEXSI) then
        nl =  num_local_elements(obj%pdist,nels,Node)
     else if (dist_type(this) == TYPE_NULL) then
        nl = 0
     else
        nl = 0
     end if

   end function num_local_elements_

   function index_local_to_global_(this,il,Node) result(ig)
     type(Dist), intent(in)  :: this
     integer, intent(in)                    :: il
     integer, intent(in)                    :: Node
     integer                                :: ig

     if (dist_type(this) == TYPE_BLOCK_CYCLIC) then
        ig = index_local_to_global(obj%bdist,il,Node)
     else if (dist_type(this) == TYPE_PEXSI) then
        ig = index_local_to_global(obj%pdist,il,Node)
     else if (dist_type(this) == TYPE_NULL) then
        ig = 0
     else
        ig = 0
     end if

   end function index_local_to_global_

   function index_global_to_local_(this,ig,Node) result(il)
     type(Dist), intent(in)  :: this
     integer, intent(in)                    :: ig
     integer, intent(in), optional          :: Node
     integer                                :: il

     if (dist_type(this) == TYPE_BLOCK_CYCLIC) then
        il = index_global_to_local(obj%bdist,ig,Node)
     else if (dist_type(this) == TYPE_PEXSI) then
        il = index_global_to_local(obj%pdist,ig,Node)
     else if (dist_type(this) == TYPE_NULL) then
        il = 0
     else
        il = 0
     end if
 
   end function index_global_to_local_

   function node_handling_element_(this,ig) result(proc)
     type(Dist), intent(in)  :: this
     integer, intent(in)                    :: ig
     integer                                :: proc

     if (dist_type(this) == TYPE_BLOCK_CYCLIC) then
        proc = node_handling_element(obj%bdist,ig)
     else if (dist_type(this) == TYPE_PEXSI) then
        proc = node_handling_element(obj%pdist,ig)
     else if (dist_type(this) == TYPE_NULL) then
        proc = MPI_UNDEFINED
     else
        proc = MPI_UNDEFINED
     end if

   end function node_handling_element_

end module class_Dist
