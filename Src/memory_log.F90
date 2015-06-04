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
! ==================================================================

! Written by J.M.Soler. May 2000.
! Re-organized by A. Garcia, June 2015
!
! ==================================================================
! SUBROUTINE memory_report( level, unit, file, printNow, threshold )
!   Sets the output file for the allocation report
! INPUT (optional):
!   integer      :: level     : Level (detail) of report
!   integer      :: unit      : Output file unit
!   character*(*):: file      : Output file name
!   logical      :: printNow  : If present & true => print report now
!   real(dp)     :: threshold : Memory threshold (in bytes) to print
!                               the memory use of any given array 
! BEHAVIOR:
!   The detail/extent of the report increses with the value of level:
! level=0 : no report at all (the default)
! level=1 : only total memory peak and where it occurred
! level=2 : detailed report created but printed only upon request
! level=3 : detailed report printed at every new memory peak
! level=4 : print every individual reallocation or deallocation
!   If unit is present, memory_report merely takes note of it for
! future use, assuming that it has been already open outside.
! In this case, file is not used.
!   If unit is absent, and file is present, a file with that
! name is open for future use.
!   If both arguments are absent, a file named 'memory_report'
! is open for future use.
!   If memory_report is called with printNow=.true. several times in
! a program, with the same unit or file argument, the subsequent 
! reports are written consecutively in the same file, each with a 
! time stamp header.
!   If threshold is not present, threshold=0 is assumed.
!   In parallel execution, the report sections that involve every
! reallocation (levels 1, 3, and 4) are written only by node 0.
! The section that is written upon request (level 2) is written
! only by the node with the highest peak of memory up to that time,
! but it contains a summary of the memory used by all other nodes.
!   In parallel execution, the nodes that share the same file
! system (e.g. different chip cores or NFS-connected nodes) write
! on the same file. Otherwise they write on files with the same name 
! in their local disks.
! ==================================================================---

MODULE memory_log

  use precision, only: dp        ! Double precision real type
  use parallel,  only: Node      ! My processor node index
  use parallel,  only: Nodes     ! Number of parallel processors
  use parallel,  only: ionode    ! Am I the I/O processor?
  use parallel,  only: parallel_init  ! Initialize parallel variables
  use m_io,      only: io_assign ! Get and reserve an available IO unit
#ifdef MPI
!  use mpi_siesta
  use mpi_siesta, only: MPI_AllGather
  use mpi_siesta, only: MPI_Barrier
  use mpi_siesta, only: MPI_Bcast
  use mpi_siesta, only: MPI_Comm_World
  use mpi_siesta, only: MPI_double_precision
  use mpi_siesta, only: MPI_integer
  use mpi_siesta, only: MPI_character
#endif

  implicit none

PUBLIC ::             &
  memory_report,       &! Sets log report defaults
  memory_event,        &! Memory counting for allocs
  type_mem              ! Converter

public :: memory   ! The old (re-furbished) routine

integer, public :: mem_stat    ! (legacy) For use in calls to allocate 
                               ! and deallocate

PRIVATE      ! Nothing is declared public beyond this point

  ! Initial default values
  character(len=*), parameter :: &
    DEFAULT_NAME = 'unknown_name'         ! Array name default
  character(len=*), parameter :: &
    DEFAULT_ROUTINE = 'unknown_routine'   ! Routine name default

  integer, save ::               &
    REPORT_LEVEL = 0,            &! Level (detail) of allocation report
    REPORT_UNIT  = 0              ! Output file unit for report

  character(len=50), save ::     &
    REPORT_FILE = 'memory_report'  ! Output file name for report
  real(dp), save ::              &
    REPORT_THRESHOLD = 0          ! Memory threshold (in bytes) to print
                                  ! the memory use of any given array 

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
  
CONTAINS

! ==================================================================

SUBROUTINE memory_report( level, unit, file, printNow, threshold )

implicit none

integer,          optional, intent(in) :: level, unit
character(len=*), optional, intent(in) :: file
logical,          optional, intent(in) :: printNow
real(dp),         optional, intent(in) :: threshold

logical open

#ifdef MPI
integer MPIerror
#endif

if (present(level)) then
  REPORT_LEVEL = level
end if

if (node == 0) then
  if (present(unit)) then  ! Assume that unit has been open outside
    if (unit > 0) then
      REPORT_UNIT = unit
      if (present(file)) then
        REPORT_FILE = file
      else
        REPORT_FILE = 'unknown'
      end if
    end if
  else if (present(file)) then    ! If file is the same, do nothing
    if (file /= REPORT_FILE) then ! Check if file was open outside
      REPORT_FILE = file
      inquire( file=REPORT_FILE, opened=open, number=REPORT_UNIT )
      if (.not.open) then         ! Open new file
        call io_assign(REPORT_UNIT)
        open( REPORT_UNIT, file=REPORT_FILE, status='unknown')
        write(REPORT_UNIT,*) ' '  ! Overwrite previous reports
      end if
    end if
  else if (REPORT_UNIT==0) then   ! No unit has been open yet
    REPORT_FILE = 'memory_report'
    call io_assign(REPORT_UNIT)
    open( REPORT_UNIT, file=REPORT_FILE, status='unknown')
    write(REPORT_UNIT,*) ' '      ! Overwrite previous reports
  end if
end if

#ifdef MPI
! Distribute information to other nodes and open REPORT_UNIT
call MPI_Bcast(REPORT_LEVEL,1,MPI_integer,0,MPI_Comm_World,MPIerror)
call MPI_Bcast(REPORT_UNIT,1,MPI_integer,0,MPI_Comm_World,MPIerror)
call MPI_Bcast(REPORT_FILE,50,MPI_character,0,MPI_Comm_World,MPIerror)
! JMS: open file only in node 0
!if (node > 0) then
!  open( REPORT_UNIT, file=REPORT_FILE )
!end if
#endif

if (present(threshold)) REPORT_THRESHOLD = threshold

if (present(printNow)) then
  if (printNow) call print_report( )
end if

END SUBROUTINE memory_report

! ==================================================================
! Internal subroutines
! ==================================================================

SUBROUTINE memory_event( bytes, aname )

implicit none

integer, intent(in)          :: bytes
character(len=*), intent(in) :: aname

character(len=1)    :: memType, task
real(DP)            :: delta_mem
logical             :: newPeak
logical,  save      :: header_written = .false.
logical,  save      :: tree_nullified = .false.
integer             :: memSize

if (REPORT_LEVEL <= 0) return

MAX_LEN = max( MAX_LEN, len(trim(aname)) )

! Find memory increment and total allocated memory
delta_mem = real(bytes,kind=dp)
TOT_MEM = TOT_MEM + delta_mem
if (TOT_MEM > PEAK_MEM+0.5_dp) then
  newPeak = .true.
  PEAK_MEM = TOT_MEM
  PEAK_ARRAY = aname
  PEAK_ROUTINE = '-'
!  print'(/,a,f18.6),a,/)',
!    'memory: Memory peak =', PEAK_MEM/MBYTE, ' Mbytes'
else
  newPeak = .false.
end if

! Add/subtract/classify array memory
if (REPORT_LEVEL > 1) then
  if (.not.tree_nullified) then
    nullify(report_tree)
    tree_nullified = .true.
  end if
  call tree_add( report_tree, aname, delta_mem )
  if (newPeak) call tree_peak( report_tree )
end if

! Print report, but only in node 0, as not all 
! processors may follow the same route here
!   The detail/extent of the report increses with the value of level:
! level=0 : no report at all (the default)
! level=1 : only total memory peak and where it occurred
! level=2 : detailed report created but printed only upon request
! level=3 : detailed report printed at every new memory peak
! level=4 : print every individual reallocation or deallocation


if (newPeak .and. (REPORT_LEVEL==1 .or. REPORT_LEVEL==3) .and. &
    node == 0) then
  call print_report
end if

if (REPORT_LEVEL == 4 .and. node == 0) then
  if (.not.header_written) then
    write(REPORT_UNIT,'(/,a7,9x,1x,a4,28x,1x,2a15)') &
     'Routine', 'Name', 'Incr. (MB)', 'Total (MB)'
    header_written = .true.
  end if
  write(REPORT_UNIT,'(a32,1x,2f15.6)') &
     aname, delta_mem/MBYTE, TOT_MEM/MBYTE
end if
END SUBROUTINE memory_event

! ==================================================================

INTEGER FUNCTION type_mem( var_type )
!
! It is not clear that the sizes assumed are universal for
! non-Cray machines...
!
implicit none
character, intent(in) :: var_type
character(len=40)     :: message

select case( var_type )
#ifdef OLD_CRAY
  case('I')
    type_mem = 8
  case('R')
    type_mem = 8
  case('L')
    type_mem = 8
#else
  case('I')
    type_mem = 4
  case('R')
    type_mem = 4
  case('L')
    type_mem = 4
#endif
case('E')
  type_mem = 8
case('D')
  type_mem = 8
case('S')
  type_mem = 1
case default
  write(message,"(2a)") &
    'memory_log: ERROR: unknown type = ', var_type
  call die(trim(message))
end select

END FUNCTION type_mem

! ==================================================================

RECURSIVE SUBROUTINE tree_add( t, name, delta_mem )

implicit none
type(TREE),       pointer    :: t
character(len=*), intent(in) :: name
real(DP),         intent(in) :: delta_mem

logical, save :: warn_negative = .true.

if (.not.associated(t)) then
  allocate( t )
  t%name = name
  t%mem  = delta_mem
  t%max  = delta_mem
  t%peak = 0._dp
  nullify( t%left, t%right )
else if (name == t%name) then
  t%mem = t%mem + delta_mem
  ! The abs is to handle the case of apparent de_allocs without re_allocs,
  ! caused by routine/name argument mismatches
  if (abs(t%mem) > abs(t%max)) t%max = t%mem
else if ( llt(name,t%name) ) then
  call tree_add( t%left, name, delta_mem )
else
  call tree_add( t%right, name, delta_mem )
end if

if (warn_negative .and. t%mem<0._dp) then
  call parallel_init()   ! Make sure that node and Nodes are initialized
  if (Node==0) then
   write(6,'(/,a,/,2a,/,a,f18.0,a)')  &
      'WARNING: alloc-realloc-dealloc name mismatch',  &
      '         Name: ', trim(name),                   &
      '         Size: ', t%mem, ' Bytes'
    if (Nodes>1) write(6,'(9x,a,i6)') 'Node:', Node
    write(6,'(9x,a)') 'Subsequent mismatches will not be reported'
    warn_negative = .false.  ! Print this warning only once
  end if
end if

END SUBROUTINE tree_add

! ==================================================================

RECURSIVE SUBROUTINE tree_peak( t )

implicit none
type(TREE), pointer :: t

if (.not.associated(t)) return

t%peak = t%mem
call tree_peak( t%left )
call tree_peak( t%right )

END SUBROUTINE tree_peak

! ==================================================================

RECURSIVE SUBROUTINE tree_print( t )

implicit none
type(TREE), pointer :: t

if (.not.associated(t)) return

call tree_print( t%left )

if (abs(t%max) >= REPORT_THRESHOLD) then
  write(REPORT_UNIT,'(a,1x,3f15.6,f9.2)') &
    t%name(1:MAX_LEN), t%mem/MBYTE, t%max/MBYTE, t%peak/MBYTE, &
    100._dp * t%peak / (PEAK_MEM + tiny(PEAK_MEM) )
end if

call tree_print( t%right )

END SUBROUTINE tree_print

! ==================================================================

SUBROUTINE print_report

implicit none

character(len=80)   :: string = 'Name'
character           :: date*8, time*10, zone*5
integer             :: iNode, peakNode
real(dp)            :: maxPeak
real(dp),allocatable:: nodeMem(:), nodePeak(:)

#ifdef MPI
integer           :: MPIerror
#endif

! Make sure that variables node and Nodes are initialized
call parallel_init()

! Allocate and initialize two small arrays
allocate( nodeMem(0:Nodes-1), nodePeak(0:Nodes-1) )

! Initializations for Nodes=1 (serial case)
nodeMem(node) = TOT_MEM
nodePeak(node) = PEAK_MEM
peakNode = node

! In parallel, find the memory values of all nodes
#ifdef MPI
if (Nodes > 1) then
  ! Gather the present and peak memories of all nodes
  call MPI_AllGather( TOT_MEM, 1, MPI_double_precision, &
                      nodeMem, 1, MPI_double_precision, &
                      MPI_COMM_WORLD, MPIerror )
  call MPI_AllGather( PEAK_MEM, 1, MPI_double_precision, &
                      nodePeak, 1, MPI_double_precision, &
                      MPI_COMM_WORLD, MPIerror )
  ! Find the node with the highest peak of memory
  maxPeak = 0
  do iNode = 0,Nodes-1
    if (nodePeak(iNode) > maxPeak) then
      peakNode = iNode
      maxPeak = nodePeak(iNode)
    end if
  end do ! iNode
  ! Change the writing node for the peak-node information
  if (node==0 .and. peakNode/=0) close( unit=REPORT_UNIT )
  call MPI_Barrier( MPI_COMM_WORLD, MPIerror )
  if (node==peakNode .and. peakNode/=0) &
    open( unit=REPORT_UNIT, file=REPORT_FILE, &
          status='unknown', position='append' )
end if ! (Nodes>1)
#endif

! The report is printed by the highest-peak node
if (node == peakNode) then

  ! AG: Commented out to allow multiple batches of information
  ! if (REPORT_LEVEL < 4) rewind(REPORT_UNIT)

  call date_and_time( date, time, zone )

  write(REPORT_UNIT,'(/,a,16a)')                &
    'Allocation summary at ',                   &
    date(1:4),'/',date(5:6),'/',date(7:8),' ',  &
    time(1:2),':',time(3:4),':',time(5:10),' ', &
    zone(1:3),':',zone(4:5)

  if (Nodes > 1) then
    write(REPORT_UNIT,'(/,(a,f18.6,a))')            &
      'Present memory all nodes : ', sum(nodeMem)/MBYTE,  ' MB', &
      'Added peak mem all nodes : ', sum(nodePeak)/MBYTE, ' MB', &
      'Min peak memory in a node: ', minval(nodePeak)/MBYTE, ' MB', &
      'Max peak memory in a node: ', maxval(nodePeak)/MBYTE, ' MB'
!   Impractical for many nodes:
!    write(REPORT_UNIT,'(/,a,/,(i6,f12.6))') &
!      'Memory peaks of individual nodes (Mb):', &
!      (iNode,nodePeak(iNode)/MBYTE,iNode=0,Nodes-1)
    write(REPORT_UNIT,'(/,a,i6)') &
      'Maximum peak of memory occurred in node:', peakNode
  end if

  write(REPORT_UNIT,'(2(/,a,f18.6,a),/,2a,/,2a)')            &
    'Present memory allocation: ', TOT_MEM/MBYTE,  ' MB', &
    'Maximum memory allocation: ', PEAK_MEM/MBYTE, ' MB', &
    'Occurred after allocating: ', trim(PEAK_ARRAY),         &
    'In routine:                ', trim(PEAK_ROUTINE)

  if (REPORT_LEVEL > 1) then
    if (REPORT_THRESHOLD > 0._dp) then
      write(REPORT_UNIT,'(/,a,f12.6,a,/,a,1x,3a15,a9)') &
        'Allocated sizes (in MByte) of arrays larger than ', &
        REPORT_THRESHOLD/MBYTE, ' MB:', &
        string(1:MAX_LEN), 'Present', 'Maximum', 'At peak', '%'
    else
      write(REPORT_UNIT,'(/,a,/,a,1x,3a15,a9)') &
        'Allocated array sizes (in MByte):', &
        string(1:MAX_LEN), 'Present', 'Maximum', 'At peak', '%'
    end if
    call tree_print( report_tree )
  end if

end if ! (node == peakNode)

! Change again the writing node for the rest of the report
#ifdef MPI
if (node==peakNode .and. peakNode/=0) close( unit=REPORT_UNIT )
call MPI_Barrier( MPI_COMM_WORLD, MPIerror )
if (node==0 .and. peakNode/=0) &
  open( unit=REPORT_UNIT, file=REPORT_FILE, &
        status='unknown', position='append' )
#endif

deallocate( nodeMem, nodePeak )

END SUBROUTINE print_report
!
! This routine is called "by hand" when allocations are not
! handled by the 'alloc' mechanism, which is only suitable
! for pointers.
!
subroutine memory( Task, Type, NElements, CallingRoutine, &
                         stat,id)
! 
! This subroutine keeps track of information relating to the use 
! of dynamic memory

! Input :
! character*1 Task  : job type = 'A' -> allocate
!                   :            'D' -> deallocate
! character*1 Type  : type of variable = 'I' = integer
!                   :                    'S' = single precision real
!                   :                    'D' = double precision real
!                   :                    'X' = grid precision real
!                   :                    'L' = logical
!                   :                    'C' = single precision complex
!                   :                    'Z' = double precision complex
!                   :                    'S' = character data (we assume takes one word)
!                   :                    'E' = double precision integer
! integer NElements : number of array elements being 
!                   : allocated/deallocated
! character         :
!   CallingRoutine  : string containing the name of the calling routine

! Created by J.D. Gale, October 1999

! Stat and ID keyword arguments added by Alberto Garcia, 2005

      implicit none

      integer, intent(in)                 :: NElements
      character(len=1), intent(in)        :: Task, Type
      character(len=*), intent(in)        :: CallingRoutine
      integer, intent(in), optional       :: stat
      character(len=*), intent(in), optional    :: id

! Local variables
      integer         :: allocSize, bytes
      character(len=1):: allocType

      if (present(stat)) then
         if (stat .ne. 0) then
            if (present(id)) then
               call die(Task // "-llocation failed in " // &
                        CallingRoutine // id)
            else
               call die(Task // "-llocation failed in " // &
                        CallingRoutine)
            endif
         endif
      endif

      select case(Type)
      case('S')
        allocType = 'R'
        allocSize = NElements
      case('C')
        allocType = 'R'
        allocSize = NElements*2
      case('Z')
        allocType = 'D'
        allocSize = NElements*2
      case('X')
#ifdef GRID_DP
        allocType = 'D'
#else
        allocType = 'R'
#endif
        allocSize = NElements
      case default
        allocType = Type
        allocSize = NElements
      end select
      if (Task=='D') allocSize = -allocSize
      bytes = allocSize*type_mem(allocType)
      call memory_event(bytes,   &
           aname=trim(CallingRoutine)//'@'//'unknown' )

   end subroutine memory

END MODULE memory_log
