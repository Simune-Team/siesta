module m_timer_tree

  !
  ! Timing framework that knows about the tree structure among
  ! timed sections: an enclosing section has sub-sections as
  ! 'children', and this relationship is preserved in the report.
  ! Future extensions include the generation of a (possibly sorted)
  ! flat report, and the accumulation of timings for specific classes
  ! of sections.
  !
  ! This module keeps a "global_section" which is the parent of all
  ! first-level sections. Users should *not* handle this section.
  !
  ! Exported routines, all accepting a string (section name) as argument:
  !
  !  timer_on    :   Flags the section as 'active' and (resumes) timing it.
  !  timer_off   :   Removes the 'active' flag and computes the elapsed time.
  !  timer_report:   Produces a report (as of now it is a full report)
  !
  !  The timing is based by default on calls to wall_time()
  !  This behavior can be changed by setting the variable 'use_walltime' to
  !  .false.
  ! 
  !  Implementation notes: 
  !    - The fixed-size array for the children is not robust enough, and
  !      should be replaced by a more flexible structure, or at least by
  !      a re-sizable pointer or allocatable array.
  !      
  !
  ! Alberto Garcia, January 2013, re-using some pieces of code by Jose Soler
  ! 
  implicit none

  integer, parameter  :: dp = selected_real_kind(10,100)
  integer, parameter  :: NMAX_CHILDREN = 10
  integer, parameter  :: maxLength = 40

  ! Derived type to hold time data
  type times_t
     private
     character(len=maxLength):: name=' '  ! Name of program or code section
     integer :: nCalls=0           ! Number of calls made to the program
     real(dp):: totTime=0          ! Total time for this program
     real(dp):: lastTime=0         ! Total time in last call
  end type times_t

  type section_t
     type(times_t) :: data
     type(section_t), pointer               :: caller => null()
     logical       :: active=.false.     ! Is program time being counted?
     integer       :: nchildren = 0
     type(section_t), dimension(:), pointer :: child
  end type section_t

  type(section_t), pointer :: global_section => null()
  type(section_t), pointer :: last_active => null()

  real(dp)       :: globaltime
  logical,public :: use_walltime = .true.

  type(section_t), pointer :: p
  type(times_t), pointer   :: pd

  real(dp)                 :: treal
  real(dp)                 :: timeNow, deltaTime

  public :: timer_on, timer_off, timer_all_off, timer_report
  private

  interface
     subroutine die(str)
       character(len=*), intent(in), optional :: str
     end subroutine die
  end interface

CONTAINS

  !------------------------------------------------
  subroutine timer_on(secname)
    character(len=*)    :: secname

    integer :: loc

    ! Use an extra enclosing level for everything,
    ! so that multiple user "trees" can be supported

    if (.not. associated(global_section)) then
       allocate(global_section)
       p => global_section
       p%active = .true.
       pd => p%data
       pd%name = "global_section"
       pd%nCalls = pd%nCalls + 1
       last_active => p
       call current_time( treal )
       pd%lastTime = treal
    endif

    ! Find proper place
    p => last_active
    if (p%nchildren==0) then
       allocate(p%child(NMAX_CHILDREN))
    endif
    call child_index(secname,p%child,loc)
    if (loc == 0) then
       ! New child
       if (p%nchildren == NMAX_CHILDREN) then
          ! re_allocate in future
          call die("too many children")
       endif
       p%nchildren = p%nchildren + 1
       !print *, "New child: " // trim(secname) // " of " // trim(p%data%name)
       p=>p%child(p%nchildren)
       p%caller => last_active
    else
       ! Existing child
       p => p%child(loc)
       p%caller => last_active
    endif

    p%active = .true.
    pd => p%data
    pd%name = secname
    pd%nCalls = pd%nCalls + 1
    last_active => p
    ! Find present time
    call current_time( treal )
    pd%lastTime = treal

  end subroutine timer_on

  !------------------------------------------------
  subroutine timer_off(secname)
    character(len=*)    :: secname

    character(len=256) :: msg

    p => last_active
    if (trim(p%data%name) /= trim(secname)) then
       write(msg,"(a,a,a,a)") "Wrong sequence in 'timer_off'. Last: ", &
                               trim(p%data%name), ". This: ", trim(secname)
       call die(msg)
    endif

    if (.not. p%active) then
       write(msg,"(a,a)") trim(p%data%name), " not active!"
       call die(msg)
    endif
    p%active = .false.
    pd => p%data

    call current_time( treal )  
    deltaTime = treal - pd%lastTime
    pd%lastTime = deltaTime
    pd%totTime  = pd%totTime + deltaTime

    if (associated(p%caller)) then
       last_active => p%caller
    else
       call die('timer_tree: attempt to close global_section?')
    endif
  end subroutine timer_off

  !------------------------------------------------
  subroutine timer_all_off()
    !
    !  Closes all outstanding active sections
    !

    if (.not. associated(last_active)) RETURN

    p => last_active
    do
       if (associated(p,global_section)) exit
       call timer_off(p%data%name)
       p => p%caller
    enddo
  end subroutine timer_all_off

  !------------------------------------------------
  subroutine timer_report(secname)
    character(len=*), optional    :: secname

    integer :: i
    type(times_t), pointer :: qd

    p => global_section
    ! Assign to the global section the sum of the times
    ! of its children
    globaltime = 0
    do i = 1, p%nchildren
       qd => p%child(i)%data
       globaltime = globaltime + qd%totTime
    enddo
    p%data%totTime = globaltime + 1.0e-6_dp

    write(*,"(/,a20,T30,a6,a12,a8)") "Section","Calls","Walltime","%"
    call walk_tree(p,0)

  end subroutine timer_report

  !------------------------------------------------
  recursive subroutine walk_tree(p,level)
    type(section_t), intent(in),target  :: p
    integer, intent(in)          :: level

    integer :: i
    character(len=40) fmtstr

    pd => p%data
    write(fmtstr,"(a,i0,a1,a)") "(", level+1, "x", ",a20,T30,i6,f12.3,f8.2)"
    write(*,fmtstr) pd%name, pd%nCalls, pd%totTime, 100*pd%totTime/globaltime
    if (p%nchildren /= 0) then
       do i=1,p%nchildren
          call walk_tree(p%child(i),level+1)
       enddo
    endif

  end subroutine walk_tree

  !------------------------------------------------
  subroutine child_index( child, childData, ichild ) 

    implicit none
    character(len=*),intent(in) :: child  ! Child name
    type(section_t), dimension(:) ,intent(in) :: childData
    integer,intent(out)                       :: ichild     ! 0 if not found

    integer :: kChild, nsize

    ichild = 0
    nsize = size(childData)

    if (nsize>0) then
       search: do kChild = 1,nsize
          if (trim(childData(kChild)%data%name) == trim(child)) then
             iChild = kChild
             exit search
          end if
       end do search
    end if

  end subroutine child_index

  !------------------------------------------------
  subroutine current_time(t)
    !
    ! CPU or walltime, depending on the setting of 'use_walltime'
    !
    use m_walltime, only: wall_time
    real(dp), intent(out) :: t

    if (use_walltime) then
       call wall_time(t)
    else
       call cpu_time(treal)
       t = treal
    endif
  end subroutine current_time

end module m_timer_tree
