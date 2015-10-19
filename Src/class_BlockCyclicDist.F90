!
  module class_BlockCyclicDist
#ifdef MPI
    use mpi
#else    
    integer, parameter, private :: MPI_GROUP_NULL = -huge(1)
    integer, parameter, private :: MPI_COMM_NULL = -huge(1)
    integer, parameter, private :: MPI_UNDEFINED = -huge(1)
#endif
  implicit none

  character(len=*), parameter :: mod_name=__FILE__

  type BlockCyclicDist_
     integer              :: refCount = 0
     character(len=36)    :: id
     !------------------------
     character(len=256)   :: name = "null BlockCyclicDist"
     !------------------------
     integer  :: ref_comm = MPI_COMM_NULL
     integer  :: group     = MPI_GROUP_NULL      ! MPI group
     integer  :: node      = MPI_UNDEFINED       ! MPI rank in group  (my_proc)
     integer, allocatable :: ranks_in_ref_comm(:)
     integer  :: nodes = 0       ! MPI size of group  (nprocs)
     integer  :: node_io = -1    ! Node capable of IO
     !------------------------
     integer  :: blocksize = 0   ! 
     !                   
     integer  :: isrcproc = 0    ! Processor holding the first element (unused)
     !------------------------------------------
  end type BlockCyclicDist_

  type BlockCyclicDist
     type(BlockCyclicDist_), pointer :: data => null()
  end type BlockCyclicDist

  type(BlockCyclicDist_), pointer :: obj => null()

  public :: newDistribution
!  public :: group
  public :: ref_comm, ranks_in_ref_comm
  public :: num_local_elements, node_handling_element
  public :: index_local_to_global, index_global_to_local

  interface newDistribution
     module procedure newBlockCyclicDistribution
  end interface newDistribution

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
#define TYPE_NAME BlockCyclicDist
#include "basic_type.inc"
!====================================    

     subroutine delete_Data(spdata)
      type(BlockCyclicDist_) :: spdata
      !
      integer :: ierr
      if (spdata%group /= MPI_GROUP_NULL) then
         call MPI_Group_free(spdata%group,ierr)
      endif
      if (allocated(spdata%ranks_in_ref_comm)) then
         deallocate(spdata%ranks_in_ref_comm)
      endif
     end subroutine delete_Data

     subroutine newBlockCyclicDistribution(this, &
          ref_comm,ranks_in_ref_comm,BlockSize,name)
     !........................................
     ! Constructor
     !........................................
     type (BlockCyclicDist), intent(inout)     :: this
     integer, intent(in)                       :: ref_comm
     integer, intent(in)                       :: ranks_in_ref_comm(:)
     integer, intent(in)                       :: Blocksize
     character(len=*), intent(in), optional :: name

     integer :: error, gsize, ref_group

     call init(this)
     obj => this%data

     obj%ref_comm = ref_comm
     obj%blocksize = Blocksize
     
     ! Do we need to allocate, or a simple assignment would suffice?
     gsize = size(ranks_in_ref_comm)
     if (allocated(obj%ranks_in_ref_comm)) deallocate(obj%ranks_in_ref_comm)
     allocate(obj%ranks_in_ref_comm(gsize))
     obj%ranks_in_ref_comm(:)  = ranks_in_ref_comm(:)

     ! Define a group internally only for the purposes of
     ! getting the ranks consistently
     
#ifdef MPI
     call MPI_Comm_Group(ref_comm,ref_group,error)
     call MPI_Group_Incl(ref_group,gsize,obj%ranks_in_ref_comm,obj%group,error)
     call MPI_Group_Rank( obj%Group, obj%node, error )
     call MPI_Group_Size( obj%Group, obj%nodes, error )
     call MPI_Group_free(ref_group, error)
#else
     obj%node = 0
     obj%nodes = 1
#endif
!     print *, "bc: node, ranks in ref: ", obj%node, obj%ranks_in_ref_comm(:)
     obj%node_io = 0

     if (present(name)) then
        obj%name = trim(name)
     else
        obj%name = "(Distribution from BlockSize and Group)"
     endif
     call tag_new_object(this)

   end subroutine newBlockCyclicDistribution

!-----------------------------------------------------------
   function ref_comm_(this) result(comm)
     type(BlockCyclicDist), intent(in)  :: this
     integer                            :: comm

     obj => this%data
     comm = obj%ref_comm
   end function ref_comm_
!-----------------------------------------------------------
   function ranks_in_ref_comm_(this) result(ranks)
     type(BlockCyclicDist), intent(in)  :: this
     integer, allocatable               :: ranks(:)

     obj => this%data
     allocate(ranks(size(obj%ranks_in_ref_comm)))
     ranks(:) = obj%ranks_in_ref_comm
   end function ranks_in_ref_comm_
!-----------------------------------------------------------
   function num_local_elements_(this,nels,Node) result(nl)
     type(BlockCyclicDist), intent(in)  :: this
     integer, intent(in)                    :: nels
     integer, intent(in)                    :: Node
     integer                                :: nl

     integer, external :: numroc
     
     obj => this%data
     if (Node >= obj%Nodes) then
        nl = 0
     else
        nl = numroc(nels,obj%blocksize,Node,  &
                 obj%isrcproc,obj%nodes)
     endif

   end function num_local_elements_

   function index_local_to_global_(this,il,Node) result(ig)
     type(BlockCyclicDist), intent(in)  :: this
     integer, intent(in)                    :: il
     integer, intent(in)                    :: Node
     integer                                :: ig

     integer, external :: indxl2g
     
     obj => this%data
     if (Node >= obj%Nodes) then
        ig = 0
     else
        ig = indxl2g(il,obj%blocksize,Node, &
                  obj%isrcproc,obj%nodes)
     endif

   end function index_local_to_global_

   function index_global_to_local_(this,ig,Node) result(il)
     type(BlockCyclicDist), intent(in)  :: this
     integer, intent(in)                    :: ig
     integer, intent(in), optional          :: Node
     integer                                :: il

     integer :: owner, myrank, error
     integer, external :: indxg2l

     obj => this%data
     if (present(Node)) then
        if (Node >= obj%nodes) then
           il = 0
        else
           il = indxg2l(ig,obj%blocksize,Node,  &
                     obj%isrcproc,obj%nodes)
        endif
     else
        ! Assume that we only want a non-zero value if the orb
        ! is local to this node
        owner = node_handling_element_(this,ig)
        if (owner == obj%node) then
           il = indxg2l(ig,obj%blocksize,myrank,  &
                        obj%isrcproc,obj%nodes)
        else
           il = 0
        endif
     endif
 
   end function index_global_to_local_

   function node_handling_element_(this,ig) result(proc)
     type(BlockCyclicDist), intent(in)  :: this
     integer, intent(in)                    :: ig
     integer                                :: proc

     integer :: dummy_Node = 0
     integer, external :: indxg2p

     obj => this%data
     proc = indxg2p(ig,obj%blocksize,dummy_Node,  &
                    obj%isrcproc,obj%nodes)

   end function node_handling_element_

end module class_BlockCyclicDist
