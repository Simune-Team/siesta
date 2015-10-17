!
  module class_PEXSIDist
#ifdef MPI
    use mpi
#else    
    integer, parameter, private :: MPI_GROUP_NULL = -huge(1)
    integer, parameter, private :: MPI_COMM_NULL = -huge(1)
    integer, parameter, private :: MPI_UNDEFINED = -huge(1)
#endif

    implicit none

  character(len=*), parameter :: mod_name=__FILE__

  type PEXSIDist_
     integer              :: refCount = 0
     character(len=36)    :: id
     !------------------------
     character(len=256)   :: name = "null PEXSIDist"
     !------------------------
     integer  :: ref_comm  = MPI_COMM_NULL
     integer  :: group     = MPI_GROUP_NULL      ! MPI group
     integer  :: node      = MPI_UNDEFINED       ! MPI rank in group  (my_proc)
     integer, allocatable :: ranks_in_ref_comm(:)
     integer  :: nodes = 0       ! MPI size of group  (nprocs)
     integer  :: node_io = -1    ! Node capable of IO
     !------------------------
     integer  :: blocksize = 0   
     !                   
     integer  :: isrcproc = 0    ! Processor holding the first element (unused)
     !------------------------------------------
     !
  end type PEXSIDist_

  type PEXSIDist
     type(PEXSIDist_), pointer :: data => null()
  end type PEXSIDist

  type(PEXSIDist_), pointer :: obj => null()

  public :: newDistribution
  public :: ref_comm, ranks_in_ref_comm
  public :: num_local_elements, node_handling_element
  public :: index_local_to_global, index_global_to_local

  interface newDistribution
     module procedure newPEXSIDistribution
  end interface

  interface ref_comm
     module procedure ref_comm_
  end interface ref_comm

  interface ranks_in_ref_comm
     module procedure ranks_in_ref_comm_
  end interface ranks_in_ref_comm

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
#define TYPE_NAME PEXSIDist
#include "basic_type.inc"
!====================================    

     subroutine delete_Data(spdata)
      type(PEXSIDist_) :: spdata
      ! do nothing
     end subroutine delete_Data

     subroutine newPEXSIDistribution(this, &
          ref_comm,ranks_in_ref_comm,BlockSize,name)
     !........................................
     ! Constructor
     !........................................
     type (PEXSIDist), intent(inout) :: this
     integer, intent(in)                       :: ref_comm
     integer, intent(in)                       :: ranks_in_ref_comm(:)
     integer, intent(in)                       :: Blocksize
     character(len=*), intent(in), optional :: name

     integer :: error, gsize, ref_group

     call init(this)
     obj => this%data
     
     obj%blocksize = Blocksize
     obj%ref_comm = ref_comm
     
     ! Do we need to allocate, or a simple assignment would suffice?
     gsize = size(ranks_in_ref_comm)
     if (allocated(obj%ranks_in_ref_comm)) deallocate(obj%ranks_in_ref_comm)
     allocate(obj%ranks_in_ref_comm(gsize))
     obj%ranks_in_ref_comm(:)  = ranks_in_ref_comm(:)

     ! Define a group internally only for the purposes of
     ! getting the ranks consistently

#ifdef MPI
     call MPI_Comm_Group(ref_comm,ref_group,error)
     call MPI_Group_Incl(ref_group,gsize,ranks_in_ref_comm,obj%group,error)
     call MPI_Group_Rank( obj%Group, obj%node, error )
     call MPI_Group_Size( obj%Group, obj%nodes, error )
#else
     obj%node = 0
     obj%nodes = 1
#endif
     obj%node_io = 0

     if (present(name)) then
        obj%name = trim(name)
     else
        obj%name = "(PEXSI Dist from BlockSize and Group)"
     endif
     call tag_new_object(this)

   end subroutine newPEXSIDistribution

!-----------------------------------------------------------
   function ref_comm_(this) result(comm)
     type(PEXSIDist), intent(in)  :: this
     integer                            :: comm

     obj => this%data
     comm = obj%ref_comm
   end function ref_comm_

!-----------------------------------------------------------
   function ranks_in_ref_comm_(this) result(ranks)
     type(PEXSIDist), intent(in)        :: this
     integer, allocatable               :: ranks(:)

     obj => this%data
     allocate(ranks(size(obj%ranks_in_ref_comm)))
     ranks(:) = obj%ranks_in_ref_comm
   end function ranks_in_ref_comm_
!-----------------------------------------------------------
   function num_local_elements_(this,nels,Node) result(nl)
     type(PEXSIDist), intent(in)  :: this
     integer, intent(in)                    :: nels
     integer, intent(in)                    :: Node
     integer                                :: nl

     integer :: remainder

     remainder = nels - obj%blocksize * obj%nodes
     if (Node == (obj%Nodes - 1)) then
        nl = obj%blocksize + remainder
     else if (Node >= obj%Nodes) then
        nl = 0
     else 
        nl = obj%blocksize
     endif

   end function num_local_elements_

   function index_local_to_global_(this,il,Node) result(ig)
     type(PEXSIDist), intent(in)  :: this
     integer, intent(in)                    :: il
     integer, intent(in)                    :: Node
     integer                                :: ig

     if (Node >= obj%Nodes) then
        ig = 0
     else
        ig = obj%blocksize*Node + il
     endif

   end function index_local_to_global_

   function index_global_to_local_(this,ig,Node) result(il)
     type(PEXSIDist), intent(in)  :: this
     integer, intent(in)                    :: ig
     integer, intent(in), optional          :: Node
     integer                                :: il

     integer :: owner, myrank, error

     if (present(Node)) then
        if (Node >= obj%nodes) then
           il = 0
        else
           il = ig - obj%blocksize*Node 
        endif
     else
        ! Assume that we only want a non-zero value if the orb
        ! is local to this node
        owner = node_handling_element_(this,ig)
        if (owner == obj%node) then
           il = ig - obj%blocksize * obj%node
        else
           il = 0
        endif
     endif
 
   end function index_global_to_local_

   function node_handling_element_(this,ig) result(proc)
     type(PEXSIDist), intent(in)  :: this
     integer, intent(in)                    :: ig
     integer                                :: proc

     ! Assume bs=2, norbs=13, nodes=5
     ! Then, the distribution is 2, 2, 2, 2, 5 for procs 0,1,2,3,4
     ! For ig=13, proc=6 according to the naive calculation (first line)
     ! We have to correct this to assign this orb to proc 4.
     ! Same for ig=10: proc=5
     ! In this case the load balancing is quite bad.

     proc = (ig-1) / obj%blocksize
     if (proc > obj%Nodes-1) proc = obj%Nodes - 1

   end function node_handling_element_


end module class_PEXSIDist
