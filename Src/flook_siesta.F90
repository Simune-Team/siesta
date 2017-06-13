module flook_siesta

#ifdef SIESTA__FLOOK
  use flook

  use siesta_options
  use siesta_geom
  use parallel, only : IONode, Node, Nodes

#endif

  implicit none

  private

  ! Signals to LUA
  ! Right after reading initial options 
  integer, parameter, public :: LUA_INITIALIZE = 1
  ! Right before SCF step starts, but at each MD step
  integer, parameter, public :: LUA_INIT_MD = 2
  ! at the start of each SCF step
  integer, parameter, public :: LUA_SCF_LOOP = 3
  ! after each SCF has finished
  integer, parameter, public :: LUA_FORCES = 4
  ! when moving the atoms, right after the FORCES step
  integer, parameter, public :: LUA_MOVE = 5
  ! when SIESTA is complete, just before it exists
  integer, parameter, public :: LUA_ANALYSIS = 6

#ifdef SIESTA__FLOOK

  public :: slua_init, slua_call, slua_close

  ! Internal parameters
  logical, save :: slua_run = .false.
  character(len=512), save, public :: slua_file = ' '
  ! Debugging flag for both parallel and serial debugging
  logical, save, public :: slua_debug = .false.

contains

  subroutine slua_init(LUA)

    use fdf, only : fdf_get
    use m_os, only : file_exist
    
    type(luaState), intent(inout) :: LUA

    character(len=30) :: fortran_msg

    character(*), parameter :: fortran_static_lua = '&
siesta = { &
    Node = 1, &
    INITIALIZE = 1, &
    INIT_MD = 2, &
    SCF_LOOP = 3, &
    FORCES = 4, &
    MOVE = 5, &
    ANALYSIS = 6, &
    state = 0, &
    IOprint = function(self, ...) &
       if self.IONode then &
          print(...) &
       end &
    end, &
    print = function(self, ...) &
       print(...) &
    end, &
} &
IOprint = function(...) &
   siesta:IOprint(...) &
end &
siesta_comm = function(...) end'

    character(*), parameter :: unit_static_lua = '&
siesta.Units = { &
    Ang    = 1. / 0.529177, &
    eV     = 1. / 13.60580, &
    kBar   = 1. / 1.47108e5, &
    Debye  = 0.393430, &
    amu    = 2.133107, &
} &
siesta.Units.GPa = siesta.Units.kBar * 10 &
siesta.Units.Kelvin = siesta.Units.eV / 11604.45'

    ! For error-handling with lua
    integer :: err
    character(len=2048) :: err_msg
    
    ! First retrieve lua file
    slua_file = fdf_get('LUA.Script',' ')
    ! Immediately return if the file is not specified...
    if ( len_trim(slua_file) == 0 ) return

    ! Default debugging only on the io-node.
    slua_debug = fdf_get('LUA.Debug',.false.)
    slua_debug = slua_debug .and. IONode 
    if ( fdf_get('LUA.Debug.MPI',.false.) ) then
       ! Only if requesting parallel debug should all processors
       ! use the debugging.
       slua_debug = .true.
    end if

    ! Check that all sees the files
    ! Currently this is a limitation of this simple
    ! flook hook.
    ! This is because the calling of the fortran
    ! function from Lua is not following the
    ! same path. One could easily provide
    ! signals to the other routines by send-catch
    ! messages, but for now this is not the scope
    slua_run = file_exist(slua_file, all = .true. )

    ! If not all processors sees the file, 
    ! we return immediately...
    if ( .not. slua_run ) then
       ! Let the user know if their system
       ! does not allow running flook
       if ( IONode .and. file_exist(slua_file) ) then
          write(*,'(a)') 'siesta-lua: WARNING'
          write(*,'(a)') 'siesta-lua: The file-system does not &
               &support file-views on all running processors.'
          write(*,'(a)') 'siesta-lua: siesta-lua CANNOT be runned on &
               &this file-system. Try running with only one node.'
          write(*,'(a)') 'siesta-lua: WARNING'
       else if ( IONode ) then
          write(*,'(a)') 'siesta-lua: WARNING'
          write(*,'(3a)') 'siesta-lua: File ', trim(slua_file), &
               ' could not be found!'
          write(*,'(a)') 'siesta-lua: WARNING'
       end if
       return
    end if
    
    ! Initialize the Lua state
    call lua_init(LUA)

    ! Create LUA table for data container
    call lua_run(LUA, code = fortran_static_lua )
    ! Append the unit table for SIESTA unit conversion
    call lua_run(LUA, code = unit_static_lua )

    ! Register siesta calls to communicate to the lua layer
    ! Old names for backwards compatibility
    call lua_register(LUA,'siesta_get', slua_receive_siesta)
    call lua_register(LUA,'siesta_return', slua_send_siesta)

    call lua_register(LUA,'siesta_receive', slua_receive_siesta)
    call lua_register(LUA,'siesta_send', slua_send_siesta)
    ! Make local siesta.receive and siesta.send
    call lua_run(LUA, code = 'siesta.receive = siesta_receive' )
    call lua_run(LUA, code = 'siesta.send = siesta_send' )

    ! Only used for printing information about
    ! what can be retrieved
    ! This function will return different things
    ! dependen on where in the routine it is called
    call lua_register(LUA,'_internal_print_allowed', slua_siesta_print_objects)
    call lua_run(LUA, code = 'siesta.print_allowed = _internal_print_allowed' )

    write(fortran_msg,'(a,i0)') 'siesta.Node = ',Node + 1
    call lua_run(LUA, code = fortran_msg )
    write(fortran_msg,'(a,i0)') 'siesta.Nodes = ',Nodes
    call lua_run(LUA, code = fortran_msg )
    if ( IONode ) then
       call lua_run(LUA, code = 'siesta.IONode = true' )
    else
       call lua_run(LUA, code = 'siesta.IONode = false' )
    end if

    ! Run the requested lua-script
    err_msg = " "
    call lua_run(LUA, slua_file, error = err, message=err_msg)
    if ( err /= 0 ) then
       write(*,'(a)') trim(err_msg)
       call die('LUA initialization failed, please check your Lua script!!!')
    end if
    
  end subroutine slua_init

  subroutine slua_call(LUA, state)
    type(luaState), intent(inout) :: LUA
    integer, intent(in) :: state
    character(len=30) :: tmp

    ! For error-handling with lua
    integer :: err
    character(len=2048) :: err_msg

    ! Return immediately if we should not run
    if ( .not. slua_run ) return

    ! Transfer the state to the lua interpreter such
    ! that decisions can be made as to which steps
    ! to take
    write(tmp,'(a,i0)') 'siesta.state = ',state
    call lua_run(LUA, code = tmp )

    if ( slua_debug ) then
       write(*,'(a,i0)') 'siesta-lua: calling siesta_comm() @ ',state
    end if

    ! Call communicator
    call lua_run(LUA, code = 'siesta_comm()', error = err, message=err_msg )
    if ( err /= 0 ) then
       write(*,'(a)') trim(err_msg)
       call die('LUA could not run siesta_comm() without an error, please &
            &check your Lua script')
    end if

  end subroutine slua_call

  subroutine slua_close(LUA)
    type(luaState), intent(inout) :: LUA
    ! Return immediately if we should not run
    if ( .not. slua_run ) return
    call lua_close(LUA)
  end subroutine slua_close


  ! ! ! ! ! ! ! 
  ! The remaining functions/routines are private
  ! methods.
  ! ! ! ! ! ! !


  function slua_receive_siesta(state) result(nret) bind(c)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_int

    use dictionary
    use siesta_dicts

    ! Define the state
    type(c_ptr), value :: state
    ! Define the in/out
    integer(c_int) :: nret

    type(luaState) :: lua
    type(luaTbl) :: tbl

    type(dict) :: keys

    if ( slua_debug ) then
       write(*,'(a,i0)') '  lua: siesta_receive, Node = ',Node + 1
    end if

    call lua_init(LUA,state)

    ! Retrieve information
    call slua_get_tbl_to_dict(lua,keys)

    ! open global siesta table
    tbl = lua_table(LUA,'siesta')
    
    ! Expose the dictionary
    call slua_put_dict(tbl,options,keys)
    call slua_put_dict(tbl,variables,keys)

    call lua_close_tree(tbl)

    ! Clean-up
    call delete(keys)

    ! this function returns nothing
    nret = 0

  end function slua_receive_siesta

  function slua_send_siesta(state) result(nret) bind(c)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_int

    use dictionary
    use siesta_dicts

    ! Define the state
    type(c_ptr), value :: state
    ! Define the in/out
    integer(c_int) :: nret

    type(luaState) :: lua
    type(luaTbl) :: tbl
    type(dict) :: keys

    if ( slua_debug ) then
       write(*,'(a,i0)') '  lua: siesta_send, Node = ',Node + 1
    end if

    call lua_init(LUA,state)

    ! Retrieve information
    call slua_get_tbl_to_dict(lua,keys)

    ! open global siesta table
    tbl = lua_table(LUA,'siesta')
    
    ! Expose the dictionary
    call slua_get_dict(tbl,options,keys)
    call slua_get_dict(tbl,variables,keys)

    call lua_close_tree(tbl)

    ! Clean-up
    call delete(keys)

    ! this function returns nothing
    nret = 0

  end function slua_send_siesta

  subroutine slua_get_tbl_to_dict(lua,keys)

    use variable
    use dictionary

    type(luaState), intent(inout) :: lua
    type(dict), intent(inout) :: keys

    type(luaTbl) :: tbl
    character(len=255) :: name
    integer :: i, N

    ! Clean the dictionary
    call delete(keys)

    ! Retrieve the table @ the top
    tbl = lua_table(LUA)

    ! Traverse all elements in the table
    N = len(tbl)
    do i = 1 , N
       call lua_get(tbl,i,name)
       keys = keys // (trim(name).kv.1)
    end do
    ! Loop through all keys in it.
!    print *,'Number of elements passed: ',len(tbl),len(keys)
!    call print(keys)

    call lua_close(tbl)

  end subroutine slua_get_tbl_to_dict
  

  subroutine slua_put_dict(tbl,dic,keys)

    use variable
    use dictionary

    type(luaTbl), intent(inout) :: tbl
    type(dict), intent(inout) :: dic
    type(dict), intent(inout), optional :: keys

    ! Sadly we need a pointer for all variables that we might
    ! expect to handle.
    character(len=DICT_KEY_LENGTH) :: key
    type(dict) :: pd ! pointer to the dictionary
    type(var) :: v

    if ( present(keys) ) then

       ! Loop over all entries in the keys dictionary
       pd = .first. keys
       do while ( .not. (.empty. pd) )
          key = trim(.key. pd)
          if ( key .in. dic ) then
             call associate(v,dic,key)
             call slua_put_var(key)
             call nullify(v) ! do not delete, simply nullify
          end if
          pd = .next. pd
       end do
       
    else
       
       ! Loop over all entries
       pd = .first. dic
       do while ( .not. (.empty. pd) )
          key = .key. pd
          call associate(v,dic,trim(key))
          call slua_put_var(key)
          call nullify(v)
          pd = .next. pd
       end do

    end if

  contains

    subroutine slua_put_var(key)
      character(len=*), intent(in) :: key
      character(len=1), pointer :: a1(:)
      logical, pointer :: b0, b1(:), b2(:,:)
      integer, pointer :: i0, i1(:), i2(:,:)
      real(sp), pointer :: s0, s1(:), s2(:,:)
      real(dp), pointer :: d0, d1(:), d2(:,:)
      character(len=4) :: t
      character(len=255) :: lkey, rkey
      integer :: lvls
      lvls = 0
      t = which(v)
      if ( slua_debug ) then
         write(*,'(4a)') '    siesta2lua; dtype = ',t,', var = ',trim(key)
      end if
!      print *,'Attempt storing: ',trim(key), ' type= ',t
      select case ( t )
      case ( 'a1', 'b0' , 'i0', 's0', 'd0' )
         ! We need to handle the variables secondly
         call key_split(key,lkey,rkey)
         if ( len_trim(rkey) == 0 ) then
            lvls = 0
            rkey = lkey
         else
            call lua_open(tbl,lkey,lvls = lvls)
         end if
      end select
      select case ( t )
      case ( 'a1' )
         call associate(a1,v)
         lkey = cunpack(a1)
         call lua_set(tbl,rkey,lkey(1:len_trim(lkey)))
      case ( 'b0' )
         call associate(b0,v)
         call lua_set(tbl,rkey,b0)
      case ( 'b1' )
         call associate(b1,v)
         call lua_set(tbl,key,b1)
      case ( 'b2' ) 
         call associate(b2,v)
         call lua_set(tbl,key,b2)
      case ( 'i0' ) 
         call associate(i0,v)
         call lua_set(tbl,rkey,i0)
      case ( 'i1' ) 
         call associate(i1,v)
         call lua_set(tbl,key,i1)
      case ( 'i2' ) 
         call associate(i2,v)
         call lua_set(tbl,key,i2)
      case ( 's0' ) 
         call associate(s0,v)
         call lua_set(tbl,rkey,s0)
      case ( 's1' ) 
         call associate(s1,v)
         call lua_set(tbl,key,s1)
      case ( 's2' ) 
         call associate(s2,v)
         call lua_set(tbl,key,s2)
      case ( 'd0' ) 
         call associate(d0,v)
         call lua_set(tbl,rkey,d0)
!         print *,'setting: '//trim(key)//' to ',d0
      case ( 'd1' ) 
         call associate(d1,v)
         call lua_set(tbl,key,d1)
      case ( 'd2' ) 
         call associate(d2,v)
         call lua_set(tbl,key,d2)
!         print *,'setting: '//trim(key)//' to ',d2
      end select
!      print *,'Done storing: ',trim(key), ' type= ',t
      call lua_close(tbl,lvls = lvls)
    end subroutine slua_put_var

  end subroutine slua_put_dict

  subroutine slua_get_dict(tbl,dic,keys)

    use variable
    use dictionary

    type(luaTbl), intent(inout) :: tbl
    type(dict), intent(inout) :: dic
    type(dict), intent(in), optional :: keys

    ! Sadly we need a pointer for all variables that we might
    ! expect to handle.
    character(len=DICT_KEY_LENGTH) :: key
    type(dict) :: pd ! pointer to the dictionary
    type(var) :: v

    if ( present(keys) ) then

       ! Loop over all entries in the keys dictionary
       pd = .first. keys
       do while ( .not. (.empty. pd) )
          key = .key. pd
          if ( key .in. dic ) then
             call associate(v,dic,trim(key))
             call slua_get_var(key)
             call nullify(v)
          end if
          pd = .next. pd
       end do
       
    else
       
       ! Loop over all entries
       pd = .first. dic
       do while ( .not. (.empty. pd) )
          key = .key. pd
          call associate(v,dic,trim(key))
          call slua_get_var(key)
          call nullify(v)
          pd = .next. pd
       end do

    end if

  contains

    subroutine slua_get_var(key)
      character(len=*), intent(in) :: key
      character(len=256) :: V0
      character(len=1), pointer :: a1(:)
      logical, pointer :: b0, b1(:), b2(:,:)
      integer, pointer :: i0, i1(:), i2(:,:)
      real(sp), pointer :: s0, s1(:), s2(:,:)
      real(dp), pointer :: d0, d1(:), d2(:,:)
      character(len=4) :: t
      character(len=255) :: lkey, rkey
      integer :: na1
      integer :: lvls
      lvls = 0
      t = which(v)
      if ( slua_debug ) then
         write(*,'(4a)') '    lua2siesta; dtype = ',t,', var = ',trim(key)
      end if
      select case ( t )
      case ( 'a1', 'b0' , 'i0', 's0', 'd0' )
         ! We need to handle the variables secondly
         call key_split(key,lkey,rkey)
         if ( len_trim(rkey) == 0 ) then
            lvls = 0
            rkey = lkey
         else
            call lua_open(tbl,lkey,lvls = lvls)
         end if
      end select
!      print *,'Attempt retrieving: ',trim(key), ' type= ',t
      select case ( t )
      case ( 'a1' )
         call associate(a1,v)
         call lua_get(tbl,rkey,V0)
         na1 = len_trim(V0)
         a1 = ' '
         a1(1:na1) = cpack(V0(1:na1))
      case ( 'b0' ) 
         call associate(b0,v)
         call lua_get(tbl,rkey,b0)
      case ( 'b1' ) 
         call associate(b1,v)
         call lua_get(tbl,key,b1)
      case ( 'b2' ) 
         call associate(b2,v)
         call lua_get(tbl,key,b2)
      case ( 'i0' ) 
         call associate(i0,v)
         call lua_get(tbl,rkey,i0)
      case ( 'i1' ) 
         call associate(i1,v)
         call lua_get(tbl,key,i1)
      case ( 'i2' ) 
         call associate(i2,v)
         call lua_get(tbl,key,i2)
      case ( 's0' ) 
         call associate(s0,v)
         call lua_get(tbl,rkey,s0)
      case ( 's1' ) 
         call associate(s1,v)
         call lua_get(tbl,key,s1)
      case ( 's2' ) 
         call associate(s2,v)
         call lua_get(tbl,key,s2)
      case ( 'd0' ) 
         call associate(d0,v)
         call lua_get(tbl,rkey,d0)
!         print *,'getting: '//trim(key)//' as ',d0
      case ( 'd1' ) 
         call associate(d1,v)
         call lua_get(tbl,key,d1)
      case ( 'd2' ) 
         call associate(d2,v)
         call lua_get(tbl,key,d2)
!         print *,'getting: '//trim(key)//' as ',d2
      end select
!      print *,'Done retrieving: ',trim(key), ' type= ',t
      call lua_close(tbl,lvls = lvls)
    end subroutine slua_get_var

  end subroutine slua_get_dict
  
  subroutine key_split(key,lkey,rkey)
    character(len=*), intent(in) :: key
    character(len=*), intent(inout) :: lkey, rkey
    integer :: i
    i = index(key,'.',back=.true.)
    if ( i > 0 ) then
       lkey = trim(adjustl(key(1:i-1)))
       rkey = key(i+1:)
    else
       lkey = trim(adjustl(key))
       rkey = ' '
    end if
  end subroutine key_split


  function slua_siesta_print_objects(state) result(nret) bind(c)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_int

    use dictionary
    use siesta_dicts

    use parallel, only : Node

    ! Define the state
    type(c_ptr), value :: state
    ! Define the in/out
    integer(c_int) :: nret

    type(dict) :: et
    character(len=DICT_KEY_LENGTH) :: key
    character(len=12) :: fmt = '(tr2,a,'','')'

    ! Currently we only let the current io-node
    ! print out information.
    nret = 0
    if ( .not. slua_debug .and. .not. IONode ) return

    ! Print out information
    write(*,'(a)') '-- siesta table structure available in LUA'
    write(*,'(a)') 'siesta = {'
    write(*,'(tr2,a,i0,'','')') 'Node = ',Node + 1

    ! Loop across all keys in the dictionaries
    et = .first. options
    do while ( .not. (.empty. et) )
       key = .key. et
       write(*,fmt) trim(key)
       et = .next. et
    end do

    ! Loop across all keys in the dictionaries
    et = .first. variables
    do while ( .not. (.empty. et) )
       key = .key. et
       write(*,fmt) trim(key)
       et = .next. et
    end do
    
    write(*,'(a)') '}'

  end function slua_siesta_print_objects

#endif

end module flook_siesta
