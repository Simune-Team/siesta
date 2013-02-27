!
  module class_PEXSIDist
  
  implicit none

  character(len=*), parameter :: mod_name=__FILE__

  type PEXSIDist_
     integer              :: refCount = 0
     character(len=36)    :: id
     !------------------------
     character(len=256)   :: name = "null PEXSIDist"
     !------------------------
     integer  :: comm = -1       ! MPI communicator
     integer  :: node = -1       ! MPI rank in comm  (my_proc)
     integer  :: nodes = 0       ! MPI size of comm  (nprocs)
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

  public :: newDistribution
  public :: num_local_elements, node_handling_element
  public :: index_local_to_global, index_global_to_local

  interface newDistribution
     module procedure newPEXSIDistribution
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
#define TYPE_NAME PEXSIDist
#include "basic_type.inc"
!====================================    

     subroutine delete_Data(spdata)
      type(PEXSIDist_) :: spdata
      ! do nothing
     end subroutine delete_Data

  subroutine newPEXSIDistribution(Blocksize,Comm,this,name)
     !........................................
     ! Constructor
     !........................................
     type (PEXSIDist), intent(inout) :: this
     integer, intent(in)                       :: Blocksize
     integer, intent(in)                       :: Comm
     character(len=*), intent(in), optional :: name

     integer :: error

     call init(this)

     this%data%blocksize = Blocksize
     this%data%comm      = Comm

#ifdef MPI
     call MPI_Comm_Rank( Comm, this%data%node, error )
     call MPI_Comm_Size( Comm, this%data%nodes, error )
#else
     this%data%node = 0
     this%data%nodes = 1
#endif
     this%data%node_io = 0

     if (present(name)) then
        this%data%name = trim(name)
     else
        this%data%name = "(PEXSI Dist from BlockSize and Comm)"
     endif
     call tag_new_object(this)

   end subroutine newPEXSIDistribution

!-----------------------------------------------------------
   function num_local_elements_(this,nels,Node) result(nl)
     type(PEXSIDist), intent(in)  :: this
     integer, intent(in)                    :: nels
     integer, intent(in)                    :: Node
     integer                                :: nl

     integer :: remainder

     remainder = nels - this%data%blocksize * this%data%nodes
     if (Node == (this%data%Nodes - 1)) then
        nl = this%data%blocksize + remainder
     else
        nl = this%data%blocksize
     endif

   end function num_local_elements_

   function index_local_to_global_(this,il,Node) result(ig)
     type(PEXSIDist), intent(in)  :: this
     integer, intent(in)                    :: il
     integer, intent(in)                    :: Node
     integer                                :: ig

     ig = this%data%blocksize*Node + il

   end function index_local_to_global_

   function index_global_to_local_(this,ig,Node) result(il)
     type(PEXSIDist), intent(in)  :: this
     integer, intent(in)                    :: ig
     integer, intent(in), optional          :: Node
     integer                                :: il

     integer :: owner, myrank, error

     if (present(Node)) then
        il = ig - this%data%blocksize*Node 
     else
        ! Assume that we only want a non-zero value if the orb
        ! is local to this node
        owner = node_handling_element_(this,ig)
        call MPI_Comm_Rank( this%data%Comm, myrank, error )
        if (owner == myrank) then
           il = ig - this%data%blocksize*myrank
        else
           il = 0
        endif
     endif
 
   end function index_global_to_local_

   function node_handling_element_(this,ig) result(proc)
     type(PEXSIDist), intent(in)  :: this
     integer, intent(in)                    :: ig
     integer                                :: proc

     proc = (ig-1) / this%data%blocksize
     if (proc == this%data%Nodes) proc = proc - 1

   end function node_handling_element_


end module class_PEXSIDist
