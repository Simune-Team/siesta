module m_timer_tree

implicit none

integer, parameter  :: dp = selected_real_kind(10,100)
integer, parameter  :: NMAX_CHILDREN = 20
integer, parameter  :: maxLength = 40

! Derived type to hold time data
  type times_t
    private
    character(len=maxLength):: name=' '  ! Name of program or code section
    integer :: nCalls=0           ! Number of calls made to the program
    real(dp):: totTime=0          ! Total time for this program
    real(dp):: comTime=0          ! Communications time
    real(dp):: lastTime=0         ! Total time in last call
    real(dp):: lastComTime=0      ! Communications time in last call
  end type times_t

  type section_t
     type(times_t) :: data
     type(section_t), pointer               :: caller => null()
     logical       :: active=.false.     ! Is program time being counted?
     integer       :: nchildren = 0
     type(section_t), dimension(:), pointer :: child
  end type section_t
     
  type(section_t), pointer :: root_section => null()
  type(section_t), pointer :: last_active => null()

  type(section_t), pointer :: p
  type(times_t), pointer   :: pd

  real                     :: treal
  real(dp)                 :: timeNow, deltaTime

  public :: timer_on, timer_off, timer_report
  private

CONTAINS

  subroutine timer_on(secname)
    character(len=*)    :: secname

    integer :: loc

    ! Use an extra enclosing level for everything,
    ! so that multiple user "trees" can be supported

      if (.not. associated(root_section)) then
         allocate(root_section)
         p => root_section
         p%active = .true.
         pd => p%data
         pd%name = "root_section"
         pd%nCalls = pd%nCalls + 1
         last_active => p
         call cpu_time( treal )         ! Standard Fortran95
         pd%lastTime = treal
      endif

      ! Find proper place
         p => last_active
         if (p%nchildren==0) then
            allocate(p%child(NMAX_CHILDREN))
         endif
         call child_index(secname,p%child,loc)
         if (loc == 0) then
            if (p%nchildren == NMAX_CHILDREN) STOP "too many children"
            p%nchildren = p%nchildren + 1
!            print *, "New child: " // trim(secname) // " of " // trim(p%data%name)
            p=>p%child(p%nchildren)
            p%caller => last_active
         else
            p => p%child(loc)
            p%caller => last_active
         endif

               
      p%active = .true.
      pd => p%data
      pd%name = secname
      pd%nCalls = pd%nCalls + 1
      last_active => p
      ! Find present CPU time and convert it to double precision
      call cpu_time( treal )         ! Standard Fortran95
      pd%lastTime = treal

    end subroutine timer_on

  subroutine timer_off(secname)
    character(len=*)    :: secname

    p => last_active
    if (trim(p%data%name) /= trim(secname)) then
       print *, "on/off mismatch: ", trim(p%data%name), " ", trim(secname)
       stop
    endif

    if (.not. p%active) stop "not active!"
    p%active = .false.
    pd => p%data

    call cpu_time( treal )  
    deltaTime = treal - pd%lastTime
    pd%lastTime = deltaTime
    pd%totTime  = pd%totTime + deltaTime

!    print *, "stopping: " // trim(secname)

    if (associated(p%caller)) then
       last_active => p%caller
    else
!       print *, "reached top in off"
    endif
  end subroutine timer_off

  subroutine timer_report(secname)
    character(len=*), optional    :: secname

    p => root_section
    call walk_tree(p,0)
  end subroutine timer_report

  recursive subroutine walk_tree(p,level)
  type(section_t), intent(in),target  :: p
  integer, intent(in)          :: level

  integer :: i
  character(len=40) fmtstr

  pd => p%data
  write(fmtstr,"(a,i0,a1,a)") "(", level+1, "x", ",a20,T30,i6,f12.2)"
  write(*,fmtstr) pd%name, pd%nCalls, pd%totTime
  if (p%nchildren /= 0) then
     do i=1,p%nchildren
        call walk_tree(p%child(i),level+1)
     enddo
  endif

end subroutine walk_tree


subroutine child_index( child, childData, ichild )   ! Get index of child in data

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

end module m_timer_tree
