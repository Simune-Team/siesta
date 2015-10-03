! Module for all mixing methods in a standard way

! This module implements mixing of the Pulay and Broyden
! type.

! The Pulay method is implemented in the fast calculation
! setup and in the stable method.
! The stable method is executed if the inversion fails.
!  - Stable: G.Kresse and J.Furthmuller, Comp. Mat. Sci. 6, 15, 1996
!  - gr (guarenteed-reduction) : http://arxiv.org/pdf/cond-mat/0005521.pdf

! The Broyden scheme is implemented

! All implemented methods employ a restart with variable
! history saving.

module m_mixing
  
  use precision, only: dp
  use parallel, only: IONode, Node
  use alloc, only: re_alloc, de_alloc

  ! MPI stuff
#ifdef MPI
  use mpi_siesta
#endif

  ! Intrinsic classes for retaining history
  use class_dData1D
  use class_Fstack_dData1D

  implicit none
  
  private

  save

  integer, parameter :: MIX_LINEAR = 1
  integer, parameter :: MIX_PULAY = 2
  integer, parameter :: MIX_BROYDEN = 3
  integer, parameter :: MIX_FIRE = 4

  type tMixer

     ! Name of mixer
     character(len=20) :: name

     ! The different saved variables per iteration
     ! and their respective stacks
     type(Fstack_dData1D), allocatable :: stack(:)
     
     ! The method of the mixer
     integer :: m = MIX_PULAY
     ! In case the mixing method has a variant
     ! this denote the variant
     ! This value is thus specific for each method
     integer :: v = 0

     ! Different mixers may have different histories
     integer :: n_hist = 2
     
     ! Number of iterations using this mixer
     ! There are a couple of signals here
     !   < 0 :
     !     switch mixer after convergence
     !  == 0 :
     !     only use this mixer until convergence
     !   > 0 :
     !     after having runned n_itt, goto
     !     index given by number
     integer :: n_itt = 0
     ! The currently reached iteration
     integer :: cur_itt = 0
     ! When mod(cur_itt,restart_itt) == 0 the history will
     ! be "killed"
     integer :: restart = 1000000000
     integer :: restart_save = 0

     ! The next mixing method following this method
     type(tMixer), pointer :: next => null()

     ! The mixing parameter used for this mixer
     real(dp) :: w = 0._dp

     ! linear array of real variables used specifically
     ! for this mixing type
     real(dp), pointer :: rv(:) => null()
     integer, pointer :: iv(:) => null()

  end type tMixer
  
  ! Debug mixing runs
  logical :: debug_mix = .false.
  ! In case of parallel mixing this also contains the node number
  character(len=20) :: debug_msg = 'mix:'
  
  public :: tMixer
  public :: mixing_init, mixing_history_clear
  public :: mixing, mixing_reset
  public :: mixing_print

  public :: MIX_LINEAR, MIX_FIRE, MIX_PULAY, MIX_BROYDEN

  interface mixing
     module procedure mixing_1d, mixing_2d
  end interface mixing

contains

  subroutine mixing_init( mixers , suffix, prefix )

    use fdf

    ! The array of mixers (has to be nullified upon entry)
    type(tMixer), allocatable, target :: mixers(:)
    character(len=*), intent(in), optional :: suffix, prefix

    ! Block constructs
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline

    ! number of history steps saved
    integer :: n_hist, n_kick, n_restart, n_save
    integer :: n_lin_before, n_lin_after
    real(dp) :: w, w_kick
    real(dp) :: w_lin_before, w_lin_after
    integer :: n_broy, n_broy_orig
    real(dp) :: w_broy, w_broy_p
    integer :: n_pulay, n_pulay_orig
    real(dp) :: w_pulay, w_pulay_damp

    logical :: bool
    type(tMixer), pointer :: m
    integer :: nm, im, im2
    character(len=10) :: lsuf, lpre
    character(len=70) :: method, variant, opt, opt2

    lpre = 'SCF.'
    if ( present(prefix) ) lpre = trim(prefix) // '.'
    lsuf = ''
    if ( present(suffix) ) lsuf = '.' // trim(suffix)

    ! ensure nullification
    if ( allocated(mixers) ) deallocate(mixers)
    
    ! General options
    if ( present(prefix) ) then
       w = 0.25_dp
    else
       w = fdf_get('DM.MixingWeight',0.25_dp)
    end if
    w = fdf_get('SCF.Mix.Weight',w)


    n_kick = fdf_get('DM.NumberKick',0)
    n_kick = fdf_get('SCF.Mix.Kick',n_kick)
    n_kick = fdf_get(trim(lpre)//'Mix.Kick',n_kick)

    ! Restart after this number of iterations
    n_restart = fdf_get('SCF.Mix.Restart',0)
    n_restart = fdf_get(trim(lpre)//'Mix.Restart',n_restart)
    n_save = fdf_get('SCF.Mix.Restart.Save',1)
    n_save = fdf_get(trim(lpre)//'Mix.Restart.Save',n_save)
    n_save = max(0,n_save)

    if ( present(prefix) ) then
       w_kick = w
    else
       w_kick = fdf_get('DM.KickMixingWeight',w)
    end if
    w_kick = fdf_get('SCF.Mix.Kick.Weight',w_kick)
    w_kick = fdf_get(trim(lpre)//'Mix.Kick.Weight',w_kick)

    n_lin_before = fdf_get('SCF.Mix.Linear.Before',0)
    n_lin_before = fdf_get(trim(lpre)//'Mix.Linear.Before',n_lin_before)
    w_lin_before = fdf_get('SCF.Mix.Linear.Before.Weight',w)
    w_lin_before = fdf_get(trim(lpre)//'Mix.Linear.Before.Weight',w)
    
    n_lin_after = fdf_get('SCF.Mix.Linear.After',0)
    n_lin_after = fdf_get(trim(lpre)//'Mix.Linear.After',n_lin_after)
    w_lin_after = fdf_get('SCF.Mix.Linear.After.Weight',w)
    w_lin_after = fdf_get(trim(lpre)//'Mix.Linear.After.Weight',w)

    ! Get original implementation quantities
    if ( present(prefix) ) then
       n_pulay_orig = n_hist
       n_broy_orig = n_hist
    else
       ! Let Mix.History decide if the others are not given
       n_hist = fdf_get('SCF.Mix.History',0)
       n_pulay_orig = fdf_get('DM.NumberPulay',max(n_hist,2))
       n_broy_orig = fdf_get('DM.NumberBroyden',max(n_hist,0))
    end if

    ! default to pulay history
    n_pulay = fdf_get('SCF.Mix.Pulay.History',n_pulay_orig)
    n_pulay = fdf_get(trim(lpre)//'Mix.Pulay.History',n_pulay)
    ! Linear mixing weight for pulay 
    w_pulay = fdf_get('SCF.Mix.Pulay.Weight',w)
    w_pulay = fdf_get(trim(lpre)//'Mix.Pulay.Weight',w_pulay)
    ! mixing weight of residual in pulay mixing
    w_pulay_damp = fdf_get('SCF.Mix.Pulay.Damping',w_pulay)
    w_pulay_damp = fdf_get(trim(lpre)//'Mix.Pulay.Damping',w_pulay_damp)

    ! default to pulay history
    n_hist = fdf_get('SCF.Mix.History',n_pulay_orig)
    n_hist = fdf_get(trim(lpre)//'Mix.History',n_hist)

    n_broy = fdf_get('DM.NumberBroyden',0)
    n_broy = fdf_get('SCF.Mix.Broyden.History',n_broy)
    n_broy = fdf_get(trim(lpre)//'Mix.Broyden.History',n_broy)
    w_broy = fdf_get('SCF.Mix.Broyden.Weight',w)
    ! Weight for prime
    w_broy_p = fdf_get('SCF.Mix.Broyden.WeightP',0.01_dp)
    w_broy_p = fdf_get(trim(lpre)//'Mix.Broyden.WeightP',w_broy_p)


    ! Retrieve default method for mixing
    if ( present(prefix) ) then
       method = fdf_get('SCF.Mix','pulay')
    else
       if ( n_broy /= 0 ) then
          method = fdf_get('SCF.Mix','broyden')
       else
          method = fdf_get('SCF.Mix','pulay')
       end if
    end if
    method = fdf_get(trim(lpre)//'Mix',trim(method))

    variant = fdf_get('SCF.Mix.Variant','original')
    variant = fdf_get(trim(lpre)//'Mix.Variant',trim(variant))

    ! Debug options
    bool = fdf_get('SCF.Mix.Debug',.false.)
    if ( fdf_get(trim(lpre)//'Mix.Debug',bool) ) then
       debug_mix = IONode
       debug_msg = 'mix:'
    end if
    bool = fdf_get('SCF.Mix.Debug.MPI',.false.)
    if ( fdf_get(trim(lpre)//'Mix.Debug.MPI',bool) ) then
       debug_mix = .true.
       write(debug_msg,'(a,i0,a)') 'mix (',Node,'):'
    end if

    ! Read in blocks of different mixers
    bool = fdf_block(trim(lpre)//'Mix',bfdf)
    if ( .not. bool ) then
       bool = fdf_block('SCF.Mix',bfdf)
    end if

    if ( bool ) then

       ! We have a block of mixers
       ! This list _only_ lists all mixing blocks that are to be defined

       nm = 0
       do while ( fdf_bline(bfdf,pline) ) 
          if ( fdf_bnnames(pline) == 0 ) cycle
          nm = nm + 1
       end do
       if ( nm == 0 ) then
          call die('mixing: No mixing schemes selected. &
               &Please at least add one mixer.')
       end if

       ! Allocate different mixers (we always add a linear
       ! mixer in case the user "forgot"
       call alloc_init(nm)

       ! Rewind to grab names.
       call fdf_brewind(bfdf)
       nm = 0
       do while ( fdf_bline(bfdf,pline) ) 

          if ( fdf_bnnames(pline) == 0 ) cycle

          nm = nm + 1
          mixers(nm)%name = fdf_bnames(pline,1)

       end do

       ! Now read all mixers for this segment and their options
       do im = 1 , nm

          m => mixers(im)

          call read_block(m, .true.)
          
       end do

    else

       nm = 1 ! the default mixer
       if ( n_lin_before > 0 ) nm = nm + 1
       if ( n_lin_after > 0 ) nm = nm + 1
       if ( n_kick > 0 ) nm = nm + 1
       
       call alloc_init(nm)

       ! Current processing index
       im = 0

       if ( n_lin_before > 0 ) then
          im = im + 1
          m => mixers(im)
          m%name = 'Linear-Before'
          m%m = MIX_LINEAR
          m%n_itt = n_lin_before
          m%w = w_lin_before
          m%next => mixers(im+1)
          call read_block(m, .false.)
       end if

       ! Setup the actual mixer
       im = im + 1
       im2 = im ! this index
       m => mixers(im)
       m%m = get_method(method)
       select case ( m%m )
       case ( MIX_LINEAR )
          m%name = 'Linear'
       case ( MIX_PULAY )
          m%name = 'Pulay'
       case ( MIX_BROYDEN )
          m%name = 'Broyden'
       case default
          call die('mix: Unknown mixing option, error.')
       end select
       m%v = get_variant(m%m,variant)
       call read_block(m, .false. )

       if ( n_lin_after > 0 ) then
          im = im + 1
          m => mixers(im)
          ! Signal to switch to this index after
          ! convergence
          mixers(im2)%n_itt = - im
          m%name = 'Linear-After'
          m%m = MIX_LINEAR
          m%w = w_lin_after
          m%n_itt = n_lin_after
          ! jump back to previous after having run a
          ! few iterations
          m%next => mixers(im2)
          call read_block(m, .false.)
       end if

       ! In case we have a kick, apply the kick here
       ! This overrides the "linear.after" option
       if ( n_kick > 0 ) then
          im = im + 1
          m => mixers(im)
          m%name = 'Linear-Kick'
          m%n_itt = 1
          m%m = MIX_LINEAR
          m%w = w_kick
          
          call read_block(m, .false.)
          m%next => mixers(im2)
          
          ! set the default mixer to kick
          mixers(im2)%n_itt = n_kick
          mixers(im2)%next => m
          
       end if

    end if

    ! Create history stack and associate correct 
    ! stack pointers
    call mixing_history_clear(mixers)

  contains
    
    function get_method(str) result(m)
      character(len=*), intent(in) :: str
      integer :: m
      
      if ( leqi(str,'linear') ) then
         m = MIX_LINEAR
      else if ( leqi(str,'pulay') ) then
         m = MIX_PULAY
      else if ( leqi(str,'broyden') ) then
         m = MIX_BROYDEN
      else if ( leqi(str,'fire') ) then
         m = MIX_FIRE
         call die('mixing: FIRE currently not supported.')
      else
         call die('mixing: Unknown mixing variant.')
      end if

    end function get_method

    function get_variant(m,str) result(v)
      integer, intent(in) :: m
      character(len=*), intent(in) :: str
      integer :: v

      v = 0
      select case ( m )
      case ( MIX_LINEAR )
         ! no variants
      case ( MIX_PULAY ) 
         v = 0
         ! We do not implement tho non-stable version
         ! There is no need to have an inferior Pulay mixer...
         if ( leqi(str,'original') .or. &
              leqi(str,'kresse') .or. leqi(str,'stable') ) then
            ! stable version, will nearly always succeed on inversion
            v = 0
         else if ( leqi(str,'gr') .or. &
              leqi(str,'guarenteed-reduction') .or. &
              leqi(str,'bowler-gillan') ) then
            ! Guarenteed reduction version
            v = 1
         end if
      case ( MIX_BROYDEN )
         v = 0
      case ( MIX_FIRE ) 
         ! no variants
      end select

    end function get_variant

    subroutine read_block(m, force)
      type(tMixer), intent(inout) :: m
      ! Force the block to exist
      logical, intent(in) :: force

      integer :: n
      
      ! create block string
      opt = trim(lpre)//'Mix'//trim(lsuf)//'.'//trim(m%name)

      ! Read the options for this mixer
      if ( fdf_block(opt,bfdf) ) then
         
         method = ' '
         variant = ' '
         
         ! First read read generic things
         
         ! read options
         do while ( fdf_bline(bfdf,pline) )
            if ( fdf_bnnames(pline) == 0 ) cycle
            
            opt = fdf_bnames(pline,1)
            
            if ( leqi(opt,'method') ) then
               
               ! setting the method
               method = fdf_bnames(pline,2)
               m%m = get_method(method)
               if ( len_trim(variant) > 0 ) then
                  m%v = get_variant(m%m,variant)
               end if

            else if ( leqi(opt,'next') ) then

               nullify(m%next)

               opt2 = fdf_bnames(pline,2)
               do im2 = 2 , nm
                  if ( im2 == im ) continue
                  if ( leqi(opt2,mixers(im2)%name) ) then
                     m%next => mixers(im2)
                     exit
                  end if
               end do

               if ( .not. associated(m%next) ) then
                  call die('mixing: Could not find next mixer. &
                       &Ensure all mixers exist and their names.')
               end if

            else if ( leqi(opt,'iterations') .or. &
                 leqi(opt,'nitt') ) then

               m%n_itt = fdf_bintegers(pline,1)

            else if ( leqi(opt,'history') ) then

               m%n_hist = fdf_bintegers(pline,1)

            else if ( leqi(opt,'restart') ) then

               m%restart = fdf_bintegers(pline,1)

            else if ( leqi(opt,'restart.save') ) then

               m%restart_save = fdf_bintegers(pline,1)
               m%restart_save = min(0,m%restart_save)

            else if ( leqi(opt,'variant') ) then

               variant = fdf_bnames(pline,2)
               m%v = get_variant(m%m,variant)

            end if

         end do
         
         call fdf_brewind(bfdf)

      else if ( force ) then

         write(*,*) 'Could not find block:'
         write(*,*) trim(opt)
         
         call die('mixing: Could not find block &
              &for mixing parameters')

      end if

      ! Initialize generic options for this
      ! mixing scheme
      select case ( m%m )

      case ( MIX_PULAY )

         m%n_hist = n_pulay
         m%w = w_pulay_damp

         allocate(m%rv(1))
         m%rv(1) = w_pulay

      case ( MIX_BROYDEN )

         ! step history as there is 1 extra 0
         ! index
         m%n_hist = m%n_hist + 1

         m%w = w_broy

         ! allocate temporary array
         n = 1 + m%n_hist * (m%n_hist + 1)
         allocate(m%rv(n))

         m%rv(1) = w_broy_p

      end select

      if ( 0 < m%restart .and. m%restart < m%n_hist ) then
         ! signal to never restart
         m%restart = 0
      else if ( m%m == MIX_PULAY .and. m%v == 1 ) then
         ! Ensure the restart is an even number
         m%restart = m%restart + mod(m%restart,2)
      end if

      m%restart_save = min(m%n_hist - 1,m%restart_save)

      ! Read the options for this mixer
      if ( fdf_block(opt,bfdf) ) then
         
         ! read options
         do while ( fdf_bline(bfdf,pline) )
            if ( fdf_bnnames(pline) == 0 ) cycle
            
            opt = fdf_bnames(pline,1)

            if ( leqi(opt,'mixing.weight') .or. &
                 leqi(opt,'alpha') .or. leqi(opt,'w') ) then
               
               m%w = fdf_breals(pline,1)
               
            else if ( leqi(opt,'weightP') .or. &
                 leqi(opt,'w.prime') .or. leqi(opt,'alpha.init') ) then

               m%rv(1) = fdf_breals(pline,1)

            else if ( leqi(opt,'damping') ) then

               m%w = fdf_breals(pline,1)
               
            end if

         end do
         
      end if

    end subroutine read_block

    subroutine alloc_init(n)
      integer :: n

      allocate(mixers(n))
      mixers(:)%w = w
      mixers(:)%n_hist = n_hist
      mixers(:)%restart = n_restart
      mixers(:)%restart_save = n_save
      
    end subroutine alloc_init

  end subroutine mixing_init


  subroutine mixing_1d( mix, iscf, n, oldF, newF )
    ! The current mixing method 
    type(tMixer), pointer :: mix
    integer, intent(in) :: iscf, n
    real(dp), intent(inout) :: oldF(n), newF(n)

    integer :: info

    ! Step iterator (so first mixing has cur_itt == 1)
    mix%cur_itt = mix%cur_itt + 1

    ! Select mixer
    select case ( mix%m )

    case ( MIX_LINEAR )

       if ( debug_mix ) &
            write(*,'(2a)') trim(debug_msg),' linear'

       call mixing_linear(mix, iscf, n, oldF, newF, info)

    case ( MIX_PULAY )
       
       if ( debug_mix ) then
          select case ( mix%v )
          case ( 0 )
             write(*,'(2a)') trim(debug_msg),' Pulay'
          case ( 1 )
             write(*,'(2a)') trim(debug_msg),' Pulay, GR'
          end select
       end if

       call mixing_pulay(mix, iscf, n, oldF, newF , info)

    case ( MIX_BROYDEN )

       if ( debug_mix ) &
            write(*,'(2a)') trim(debug_msg),' Broyden'

       call mixing_broyden(mix, iscf, n, oldF, newF, info)
       
    end select

    ! Copy over new function
!$OMP parallel workshare default(shared)
    oldF = newF
!$OMP end parallel workshare
    
    ! check whether we should change the mixer
    if ( mix%n_itt == 0 ) then
       ! do nothing, we continue, indefinetely
    else if ( mix%n_itt < 0 ) then
       ! this is checked outside
    else if ( mix%n_itt <= mix%cur_itt ) then

       call mixing_step( mix )

    end if

  end subroutine mixing_1d
  
  subroutine mixing_2d( mix, iscf, n1, n2, oldF, newF )
    type(tMixer), pointer :: mix
    integer, intent(in) :: iscf, n1, n2
    real(dp), intent(inout) :: oldF(n1,n2), newF(n1,n2)

    ! Simple wrapper for 1D
    call mixing_1d( mix, iscf, n1*n2 , oldF(1,1), newF(1,1) )
    
  end subroutine mixing_2d


  ! Returns the value array from the stack(:)
  ! Returns this array:
  !    mix%stack(sidx)(hidx) ! defaults to the last item
  function getstackval(mix,sidx,hidx) result(d1)
    type(tMixer), intent(in) :: mix
    integer, intent(in) :: sidx
    integer, intent(in), optional :: hidx

    real(dp), pointer :: d1(:)

    ! Local arrays
    type(dData1D), pointer :: dD1
    if ( present(hidx) ) then
       dD1 => get_pointer(mix%stack(sidx),hidx)
    else
       dD1 => get_pointer(mix%stack(sidx), &
            n_items(mix%stack(sidx)))
    end if
    
    d1 => val(dD1)
    
  end function getstackval


  ! Returns true if the following 
  ! "advanced" mixer is 'm'
  function is_next(mix,method,next) result(bool)
    type(tMixer), intent(in), target :: mix
    integer, intent(in) :: method
    type(tMixer), pointer, optional :: next

    logical :: bool

    type(tMixer), pointer :: m

    bool = .false.
    m => mix%next
    do while ( associated(m) )

       if ( m%m == MIX_LINEAR ) then
          m => m%next
       else if ( m%m == method ) then
          bool = .true.
          exit
       else
          ! Quit if it does not do anything
          exit
       end if

       ! this will prevent cyclic combinations
       if ( associated(m,mix) ) exit

    end do

    if ( present(next) ) then
       next => m
    end if

  end function is_next


  ! Step the mixing object and ensure that
  ! the old history is either copied over or freed
  subroutine mixing_step( mix )

    type(tMixer), pointer :: mix

    select case ( mix%m )
    case ( MIX_PULAY )

       call reset(mix%stack(1))
       call reset(mix%stack(2))
       call reset(mix%stack(3))
       
    case ( MIX_BROYDEN )
       
       ! delete all BROYDEN history
       call reset(mix%stack(1))
       call reset(mix%stack(2))
       
    end select

    if ( associated(mix%next) ) then
       mix => mix%next
       mix%cur_itt = 0
       if ( debug_mix ) write(*,'(2a)') &
            trim(debug_msg),' switching mixer...'
    end if

  end subroutine mixing_step
    

  ! Here the actual mixing methods will be employed
  subroutine mixing_linear( mix, iscf, n, oldF, newF, info)

    ! The current mixing method 
    type(tMixer), intent(inout) :: mix
    integer, intent(in) :: iscf

    integer, intent(in) :: n
    real(dp), intent(inout) :: oldF(n), newF(n)
    integer, intent(out) :: info

    type(tMixer), pointer :: next
    real(dp), pointer :: tmp1(:), tmp2(:)
    real(dp) :: alpha
    integer :: in

    info = 0

    if ( is_next(mix,MIX_PULAY,next=next) ) then

       ! If the following uses history, add that information
       ! to the history.
       call update_res(next%stack(1),n,oldF,newF)
       if ( next%v /= 1 ) then
          call update_F(next%stack(3),n,newF)
       end if
       
    else if ( is_next(mix,MIX_BROYDEN,next) ) then
       
       ! If the following uses history, add that information
       ! to the history.
       call update_res(next%stack(1),n,oldF,newF)

       tmp1 => getstackval(next,1)
       if ( n_items(next%stack(2)) > 0 ) then
          tmp2 => getstackval(next,2)
!$OMP parallel workshare default(shared)
          tmp2 = next%w * tmp1
!$OMP end parallel workshare
       else

          nullify(tmp2)
          allocate(tmp2(n))
!$OMP parallel workshare default(shared)
          tmp2 =  next%w * tmp1
!$OMP end parallel workshare
          call push_F(next%stack(2),n,tmp2)
          deallocate(tmp2)

       end if

    end if

    alpha = mix%w

    if ( debug_mix ) write(*,'(2a,e10.4)') &
         trim(debug_msg),' alpha = ',alpha

!$OMP parallel do default(shared), private(in)
    do in = 1 , n
       newF(in) = oldF(in) + alpha * (newF(in) - oldF(in))
    end do
!$OMP end parallel do

  end subroutine mixing_linear


  ! Pulay mixing
  ! This has a variant in case the first (fast)
  ! gets unstable in the inversion algorithm
  subroutine mixing_pulay( mix, iscf, n, oldF, newF, info)
    ! The current mixing method 
    type(tMixer), intent(inout) :: mix
    integer, intent(in) :: iscf

    ! Input/output arrays
    integer, intent(in) :: n
    real(dp), intent(inout) :: oldF(n), newF(n)

    ! run-time information
    integer, intent(out) :: info

    type(dData1D) :: dD1

    ! Temporary arrays for local data structures
    real(dp), dimension(:), pointer :: res, rres
    real(dp), dimension(:), pointer :: rres1, rres2

    integer :: nh, ns
    integer :: i, j
    ! Used arrays
    real(dp), allocatable :: alpha(:), b(:,:), bi(:,:)
    real(dp) :: G, ssum

#ifdef MPI
    integer :: MPIerror
#endif

    real(dp), external :: ddot

    ! The Pulay mixing variable is called G
    G = mix%w

    ! The Pulay stacks has this data layout:

    ! - stack(1)
    !   The residuals for each iteration (n_hist)
    ! - stack(2)
    !   The residuals of the residuals for each iteration (n_hist - 1)

    nh = max_size(mix%stack(1))
    ns = n_items(mix%stack(1))

    ! Add the residual to the stack
    call push_res(mix%stack(1), n, oldF, newF)

    select case ( mix%v ) 
    case ( 0 ) ! stable pulay mixing

       if ( debug_mix ) &
            write(*,'(a,2(a,i0))') trim(debug_msg), &
            ' n_hist = ',min(nh,ns+1), ' / ',nh

       ns = n_items(mix%stack(1))
       
       ! Add the residuals of the residuals if applicable
       if ( ns >= 2 ) then

          call push_rres(mix%stack(2),mix%stack(1))

       else 

          ! Store output F
          call push_F(mix%stack(3), n, newF)

          ! The first Pulay step will do linear mixing
          res => getstackval(mix,1)

!$OMP parallel workshare default(shared)
          newF = oldF + mix%rv(1) * res
!$OMP end parallel workshare

          return

       end if

       ! Update the residual to reflect the input residual
       rres1 => getstackval(mix,1,ns-1)
       rres2 => getstackval(mix,3)

!$OMP parallel workshare default(shared)
       rres1 = rres1 - rres2 + oldF
!$OMP end parallel workshare
       
       ! Store the output F (we have used it now)
       call push_F(mix%stack(3), n, newF)

    case ( 1 ) ! Guaranteed reduction Pulay

       ! The history in this scheme is a little obscure... :)
       if ( mod(mix%cur_itt,2) == 1 ) then

          if ( debug_mix ) &
               write(*,'(2a)') trim(debug_msg), &
               ' Direct mixing'
          
          ! this will happen on the:
          !   1, 3, 5, ... iterations

          if ( n_items(mix%stack(2)) > 0 ) then

             ! Update RRes[i-2] to its value that should be used
             ! subsequently.
             ! From the previous iteration this array is now
             ! -Res[i-2]
             res => getstackval(mix,1)
             rres => getstackval(mix,2)
!$OMP parallel workshare default(shared)
             rres = rres + res
!$OMP end parallel workshare

          end if

          ! we continue with the input newF, no things to do

          return

       else
          
          if ( debug_mix ) &
               write(*,'(a,2(a,i0))') trim(debug_msg), &
               ' n_hist = ',min(nh,ns+1), ' / ',nh

          ! now we can calculate the RRes[i]
          call push_rres(mix%stack(2),mix%stack(1))

       end if

    case default
       
       call die('mixing: Unknown Pulay variant.')

    end select

    ! Number of saved residuals, so far
    ns = n_items(mix%stack(1))

    ! Number of history steps for the double residual
    nh = n_items(mix%stack(2))

    ! Allocate arrays
    allocate(b(nh,nh),bi(nh,nh),alpha(nh))

    ! Solve the Pulay mixing problem
    call pulay_stable()

    ! Deallocate local arrays
    deallocate(alpha,b,bi)

    if ( info /= 0 ) then ! failed inversion, goto linear...

       if ( IONode ) &
            write(*,'(2a)') trim(debug_msg), &
            ' Pulay -- failed, > linear'

       res => getstackval(mix,1)

!$OMP parallel workshare default(shared)
       newF = oldF + mix%rv(1) * res
!$OMP end parallel workshare

    end if

    ! For the Guaranteed Reduction method we need to 
    ! update the history...
    select case ( mix%v )
    case ( 1 )

       ! In the subsequent iteration we
       ! do not need Res[i], as we have to update the Res with the
       ! minimum Res path.

       ! Get the current RRes[i-1] and Res[i]
       rres => getstackval(mix,2)
       res => getstackval(mix,1)

       ! Resubtract res to get -Res[i-1]
       ! the RRes[i-1] will be updated in the next loop
!$OMP parallel workshare default(shared)
       rres = rres - res
!$OMP end parallel workshare

       ! delete latest residual
       call pop(mix%stack(1),dD1)
       call delete(dD1)

       ! Update the current residual to reflect the
       ! used residual in the algorithm
       res => getstackval(mix,1)
!$OMP parallel workshare default(shared)
       res = res - oldF + newF
!$OMP end parallel workshare

    end select

    if ( mix%restart > 0 ) then
    if ( mod(mix%cur_itt,mix%restart) == 0 ) then

       if ( IONode ) then
          write(*,'(a)')'mix: Pulay -- resetting history'
       end if
       ! The user has requested to restart the
       ! mixing scheme now
       j = mix%restart_save
       if ( j == 0 ) then
          call reset(mix%stack(1))
          call reset(mix%stack(2))
          call reset(mix%stack(3))
       else
          call reset(mix%stack(1),-j)
          call reset(mix%stack(2),-j+1)
       end if
          
    end if
    end if

  contains

    subroutine pulay_stable()

      info = 0
      
      ! Create A_ij coefficients for inversion
      do i = 1 , nh
         
         ! Get RRes[i] array
         rres1 => getstackval(mix,2,i)

         do j = 1 , i

            ! Get RRes[j] array
            rres2 => getstackval(mix,2,j)

            ! B(i,j) = B(j,i) = dot_product(RRes[i],RRes[j])
            b(i,j) = ddot(n,rres1,1,rres2,1)
            b(j,i) = b(i,j)

         end do
         
      end do
      
#ifdef MPI
      ! Global operations, but only for the non-extended entries
      call MPI_AllReduce(b(1,1),bi(1,1),nh*nh, &
           MPI_double_precision, MPI_Sum, &
           MPI_Comm_World,MPIerror)
      ! copy over reduced arrays
      b = bi
#endif
      
      ! Get inverse of matrix
      call inverse(nh, b, bi, info)
      if ( info /= 0 ) then
         
         ! return immediately, the algorithm
         ! will try the stable one, then the linear one
         
         return
         
      else ! inversion succeeded
         
         ! Get current residual
         res => getstackval(mix,1)
         
         do i = 1 , nh
            
            ! Calculate the coefficients on all processors
            alpha(i) = 0._dp
            
            do j = 1 , nh
               
               ! Get j'th residual array
               rres => getstackval(mix,2,j)
               
               ssum = ddot(n,rres,1,res,1)
               
               alpha(i) = alpha(i) - bi(i,j) * ssum
               
            end do
            
         end do
         
#ifdef MPI
         ! Reduce the alpha
         call MPI_AllReduce(alpha(1),b(1,1),nh, &
              MPI_double_precision, MPI_Sum, &
              MPI_Comm_World,MPIerror)
         alpha(:) = b(:,1)
#endif
      end if
      
      ! if debugging print out the different variables
      if ( debug_mix ) then
         write(*,'(2a,f10.6,a,100(tr1,e10.4))') &
              trim(debug_msg),' G = ',G,', alpha = ',alpha
      end if
      
      ! Read former matrices for mixing .........
      res => getstackval(mix,1,ns)
      
      ! Copy over input dm, and add the linear mixing
!$OMP parallel workshare default(shared)
      newF = oldF + G * res
!$OMP end parallel workshare
      
      do i = 1 , nh
         
         ! Get Res[i] and RRes[i]
         res => getstackval(mix,1,i)
         rres => getstackval(mix,2,i)
         
!$OMP parallel workshare default(shared)
         newF = newF + alpha(i) * ( res + G * rres )
!$OMP end parallel workshare
         
      end do
      
    end subroutine pulay_stable

  end subroutine mixing_pulay


  ! Broyden mixing
  subroutine mixing_broyden( mix, iscf, n, oldF, newF, info)
    ! The current mixing method 
    type(tMixer), intent(inout) :: mix
    integer, intent(in) :: iscf

    ! Input/output arrays
    integer, intent(in) :: n
    real(dp), intent(inout) :: oldF(n), newF(n)
    ! information about algorithm
    integer, intent(out) :: info

    ! Temporary arrays for local data structures
    real(dp), dimension(:), pointer :: tmp1, tmp2, res

    integer :: is, ns, nm
    integer :: i, j, k, l, in
    ! Used arrays
    real(dp), allocatable, dimension(:,:) :: dFdF, a, ao, b, bb
    real(dp), allocatable, dimension(:) :: F, dFF
    real(dp), pointer :: l_w(:), l_dFdF(:), l_wp
    real(dp) :: jinv0, norm, Fnorm, dMax, rtmp
    type(dData1D) :: a1D, b1D

    real(dp), external :: ddot
#ifdef MPI
    integer :: MPIerror
#endif

    info = 0

    ! Get the initial Jacobian
    jinv0 = mix%w
    
    ! Arrays in the stack
    !  s(1) == Res[i]
    !  s(2) == J_i * Res[i]) / norm(RRes[i]) ...
    !  rv   == w' // dFdF // w //

    ! Number of saved histories, so far
    ns = n_items(mix%stack(1))
    nm = max_size(mix%stack(1))

    if ( debug_mix ) &
         write(*,'(a,2(a,i0))') trim(debug_msg), &
         ' n_hist = ',min(nm,ns+1), ' / ',nm

    if ( ns == 0 ) then

       ! Add to the history of the residual
       call push_res(mix%stack(1),n,oldF,newF)

       ! Add new array
       allocate(F(n))
       tmp1 => getstackval(mix,1,1)
!$OMP parallel workshare default(shared)
       F = tmp1 * jinv0
!$OMP end parallel workshare
       call push_F(mix%stack(2),n,F)
       deallocate(F)

       info = 0

       ! Do linear interpolation
!$OMP parallel workshare default(shared)
       newF = oldF + jinv0 * (newF - oldF)
!$OMP end parallel workshare

       return
       
    end if

    ! Broyden implementation

    ! Retrieve all data segments
    l_wp => mix%rv(1) ! w'
    i = 1
    l_dFdF => mix%rv(i+1:i+nm**2) ! dFdF matrix
    i = i + nm ** 2
    l_w => mix%rv(i+1:i+nm) ! previous weights

    ! Get the two temporary arrays of the latest iteration
    tmp1 => getstackval(mix,1,ns)
    tmp2 => getstackval(mix,2,ns)

    ! Allocate all work arrays
    allocate(F(n),dFdF(nm,nm),dFF(ns))
    allocate(a(ns,ns),ao(ns,ns),b(ns,ns),bb(ns,ns))
    
    ! Calculate Res[i+1] - Res[i]
    ! From this point on 'newF' is Res[i]
    dMax = 0._dp
    norm = 0._dp
    Fnorm = 0._dp
!$OMP parallel do default(shared), private(in), &
!$OMP&reduction(max:dMax), reduction(+:norm,Fnorm)
    do in = 1 , n
       F(in) = newF(in) - oldF(in)
       dMax = max(dMax,abs(F(in)))
       Fnorm = Fnorm + F(in) * F(in)
       ! Create the RRes[i]
       tmp1(in) = F(in) - tmp1(in)
       norm = norm + tmp1(in) * tmp1(in)
    end do
!$OMP end parallel do

#ifdef MPI
    call MPI_AllReduce(dMax,rtmp,1, &
         MPI_double_precision, MPI_Max, &
         MPI_Comm_World,MPIerror)
    dMax = rtmp
    call MPI_AllReduce(norm,rtmp,1, &
         MPI_double_precision, MPI_Sum, &
         MPI_Comm_World,MPIerror)
    norm = sqrt(rtmp)
    call MPI_AllReduce(Fnorm,rtmp,1, &
         MPI_double_precision, MPI_Sum, &
         MPI_Comm_World,MPIerror)
    Fnorm = sqrt(rtmp)
#else
    norm = sqrt(norm)
    Fnorm = sqrt(Fnorm)
#endif

    ! Store weight of this iteration
    ! This is variable weight
    l_w(ns) = exp( 1._dp / (dMax + 0.2_dp) )
    if ( debug_mix ) &
         write(*,'(2(a,e10.4))') &
         trim(debug_msg)//' weight = ',l_w(ns), &
         ' , norm = ',Fnorm

    ! normalize RRes[i] and the Jacobian
    rtmp = 1._dp / norm
!$OMP parallel workshare default(shared)
    tmp1 = tmp1 * rtmp
    tmp2 = jinv0 * tmp1 + tmp2 * rtmp
!$OMP end parallel workshare

    ! Copy over previous data
    dFdF = reshape(l_dFdF,(/nm,nm/))

    ! Calculate scalar product between <RRes[is]|RRes[ns]>
    ! and between <RRes[is]|Res[ns]>
    do is = 1 , ns - 1
       res => getstackval(mix,1,is)
       dFdF(is,ns) = ddot(n,res(1),1,tmp1(1),1)
       dFF(is) = ddot(n,res(1),1,F(1),1)
    end do
    ! NOTE that for is == ns we get 1., <RRes[ns]|RRes[ns]> / norm ** 2
    dFdF(ns,ns) = 1._dp
    ! create last scalar producet <RRes[i]|Res[ns]>
    dFF(ns) = ddot(n,tmp1(1),1,F(1),1)

#ifdef MPI
    if ( ns > 1 ) then
       call MPI_AllReduce(dFdF(1,ns),a(1,1),ns-1, &
            MPI_double_precision, MPI_Sum, &
            MPI_Comm_World,MPIerror)
       dFdF(1:ns-1,ns) = a(1:ns-1,1)
    end if
    call MPI_AllReduce(dFF(1),a(1,1),ns, &
         MPI_double_precision, MPI_Sum, &
         MPI_Comm_World,MPIerror)
    dFF(:) = a(:,1)
#endif

    ! Create symmetry scalar product
    do is = 1 , ns - 1
       dFdF(ns,is) = dFdF(is,ns)
    end do
!    if ( debug_mix ) &
!         write(*,'(2a,100(tr1,e10.4))') trim(debug_msg), &
!         ' dFdF(:,1) = ',dFdF(1:ns,1)

    ! Calculate the weights by first calculating
    ! the inverse matrix

    ! Prepare to invert matrix
    do j = 1 , ns
       do i = 1 , ns
          a(i,j) = l_w(i) * l_w(j) * dFdF(i,j)
       end do
       a(j,j) = a(j,j) + l_wp ** 2
    end do

    ! Invert the matrix
    call inverse(ns, a, b, info)
    if ( info /= 0 ) then
       
       deallocate(F,dFF,dFdF,a,ao,b,bb)
       
       return
       
    end if
    
    ! Calculate beta bar
    do j = 1 , ns

       do i = 1 , ns

          bb(i,j) = 0._dp
          
          do is = 1 , ns
             bb(i,j) = bb(i,j) - l_w(i) * l_w(is) * &
                  b(i,is) * dFdF(j,is)
          end do
          
       end do
       
       bb(j,j) = bb(j,j) + 1._dp
       
    end do

    ! Compute all coefficients of the Zk | ui
    ! This is done recursively
    ! NOTE: We are done using b, so reuse it
    do l = 1 , ns

       ! Re-calculate alpha matrix
       do k = 1 , l ! ell
          do i = 1 , l ! ell
             a(k,i) = b(k,i) * l_w(k) * l_w(i)
          end do

          if ( l > 1 ) then

          ! Add matrix
          do j = 1 , l - 1 ! ell - 1
             do i = 1 , l - 1
                a(k,j) = a(k,j) + bb(k,i) * ao(i,j)
             end do
          end do

          end if
       
       end do

       ! prep for next loop
       ao(1:l,1:l) = a(1:l,1:l)

    end do

!    if ( debug_mix ) then
!       do is = 1 , ns
!          write(*,'(2a,i0,a,100(tr1,e10.4))') &
!               trim(debug_msg),' alpha(:,',is,') = ',a(:,is)
!       end do
!    end if

    ! Update the output F
    do is = 1 , ns
       
       res => getstackval(mix,2,is)
       
       b(is,1) = sum( a(:,is) * dFF )

       if ( is == 1 ) then
!$OMP parallel workshare default(shared)
          newF = oldF + jinv0*F - b(is,1) * res
!$OMP end parallel workshare
       else
!$OMP parallel workshare default(shared)
          newF = newF - b(is,1) * res
!$OMP end parallel workshare
       end if
       
    end do
    if ( debug_mix ) then
       write(*,'(2a,100(tr1,e10.4))') &
            trim(debug_msg),' gamma = ', b(:,1)
    end if

    ! Step Broyden
    if ( ns == nm ) then
       ! ensure that we reuse old data
       dFdF(1:nm-1,1:nm-1) = dFdF(2:nm,2:nm)
       l_w(1:nm-1) = l_w(2:nm)
    end if
       
    ! Copy back result for following iteration
    l_dFdF = reshape(dFdF,(/nm*nm/))

    deallocate(dFdF,a,ao,b,bb,dFF)

    ! Now add history segment to Broyden
    if ( ns == nm ) then

       ! get the data
       call get(mix%stack(1),1,a1D)
       call get(mix%stack(2),1,b1D)

    else

       ! Create new history
       call newdData1D(a1D, n, '(res)')
       call newdData1D(b1D, n, '(Jacobian-res)')

       ns = ns + 1

    end if

    ! Get residual and Jacobian residual
    tmp1 => val(a1D)
    tmp2 => val(b1D)

!$OMP parallel workshare default(shared)
    tmp1 = F
    ! newF is the "updated" residual
    tmp2 = newF - oldF
!$OMP end parallel workshare

    ! push to stack
    call push(mix%stack(1), a1D)
    call delete(a1D)
    call push(mix%stack(2), b1D)
    call delete(b1D)

    deallocate(F)

    if ( mix%restart > 0 ) then
    if ( mod(mix%cur_itt,mix%restart) == 0 ) then

       if ( IONode ) then
          write(*,'(a)')'mix: Broyden -- resetting history'
       end if

       ! The user has requested to restart the
       ! mixing scheme now
       j = mix%restart_save
       call reset(mix%stack(1),-j)
       call reset(mix%stack(2),-j)
       
    end if
    end if

  end subroutine mixing_broyden

  ! Deletes all history, 
  ! Possibly reassign number of new histories saved
  subroutine mixing_history_clear( mixers )
    type(tMixer), intent(inout), target :: mixers(:)

    type(tMixer), pointer :: m
    integer :: im, is, ns

    ! Clean up all arrays and reference counted
    ! objects
    do im = 1 , size(mixers)

       m => mixers(im)

       ! reset history track
       m%cur_itt = 0

       ! do not try and de-allocate something not
       ! allocated
       if ( .not. allocated(m%stack) ) cycle

       ns = size(m%stack)
       do is = 1 , ns
          call delete(m%stack(is))
       end do

       ! clean-up
       deallocate(m%stack)

    end do

    ! First we allocate the "advanced" mixing methods
    do im = 1 , size(mixers)
       m => mixers(im)

       select case ( m%m ) 

       case ( MIX_PULAY )

          allocate(m%stack(3))

          ! allocate Res'[i], Pulay
          call new(m%stack(1), m%n_hist)
          ! allocate Res[i+1] - Res[i], Pulay
          call new(m%stack(2), m%n_hist-1)
          if ( m%v == 0 ) then
             ! The out of the latest iteration
             call new(m%stack(3), 1)
          end if

       case ( MIX_BROYDEN )

          allocate(m%stack(2))

          ! allocate Res[i], Broyden
          call new(m%stack(1), m%n_hist)
          ! allocate 'u', Broyden
          call new(m%stack(2), m%n_hist)

       end select

    end do

  end subroutine mixing_history_clear

  subroutine mixing_reset( mixs )
    type(tMixer), allocatable, target :: mixs(:)
    
    type(tMixer), pointer :: m

    integer :: im, is, ns

    if ( .not. allocated(mixs) ) return

    do im = 1 , size(mixs)
       m => mixs(im)
       if ( allocated(m%stack) ) then
          ns = size(m%stack)
          do is = 1 , ns
             call delete(m%stack(is))
          end do
          deallocate(m%stack)
       end if
       if ( associated(m%rv) ) then
          deallocate(m%rv)
          nullify(m%rv)
       end if
       if ( associated(m%iv) ) then
          deallocate(m%iv)
          nullify(m%iv)
       end if
    end do

    deallocate(mixs)

  end subroutine mixing_reset

  subroutine mixing_print( mixers , prefix )
    
    type(tMixer), intent(in), target :: mixers(:)
    character(len=*), intent(in), optional :: prefix

    type(tMixer), pointer :: m
    character(len=50) :: lpre, fmt

    integer :: i

    if ( .not. IONode ) return

    lpre = 'SCF'
    if ( present(prefix) ) lpre = trim(prefix)

    fmt = 'mix.'//trim(lpre)//':'

    ! Print out options for all mixers
    do i = 1 , size(mixers)

       m => mixers(i)

       select case ( m%m )
          
       case ( MIX_LINEAR )
          
          write(*,'(2a,t50,''= '',a)') trim(fmt), &
               ' Linear mixing',trim(m%name)
          write(*,'(2a,t50,''= '',f12.6)') trim(fmt), &
               '    Mixing weight',m%w

       case ( MIX_PULAY )
          
          write(*,'(2a,t50,''= '',a)') trim(fmt), &
               ' Pulay mixing',trim(m%name)

          if ( m%v == 0 ) then
             write(*,'(2a,t50,''= '',a)') trim(fmt), &
                  '    Variant','stable'
          else if ( m%v == 1 ) then
             write(*,'(2a,t50,''= '',a)') trim(fmt), &
                  '    Variant','GR'
          end if

          write(*,'(2a,t50,''= '',i0)') trim(fmt), &
               '    History steps',m%n_hist
          write(*,'(2a,t50,''= '',f12.6)') trim(fmt), &
               '    Initial linear mixing weight',m%rv(1)
          write(*,'(2a,t50,''= '',f12.6)') trim(fmt), &
               '    Damping',m%w
          if ( m%restart > 0 ) then
             write(*,'(2a,t50,''= '',i0)') trim(fmt), &
                  '    Restart steps',m%restart
             write(*,'(2a,t50,''= '',i0)') trim(fmt), &
                  '    Restart save steps',m%restart_save
          end if

       case ( MIX_BROYDEN )
          
          write(*,'(2a,t50,''= '',a)') trim(fmt), &
               ' Broyden mixing',trim(m%name)

          !write(*,'(2a,t50,''= '',a)') trim(fmt), &
          !     '    Variant','original'

          write(*,'(2a,t50,''= '',i0)') trim(fmt), &
               '    History steps',m%n_hist - 1
          write(*,'(2a,t50,''= '',f12.6)') trim(fmt), &
               '    Jacobian weight',m%w
          write(*,'(2a,t50,''= '',f12.6)') trim(fmt), &
               '    Weight prime',m%rv(1)
          if ( m%restart > 0 ) then
             write(*,'(2a,t50,''= '',i0)') trim(fmt), &
                  '    Restart steps',m%restart
             write(*,'(2a,t50,''= '',i0)') trim(fmt), &
                  '    Restart save steps',m%restart_save
          end if
          
       case ( MIX_FIRE )
          
          write(*,'(2a,t50,''= '',a)') trim(fmt), &
               ' Fire mixing',trim(m%name)

       end select

       if ( m%n_itt < 0 ) then
          write(*,'(2a,t50,''= '',a)') trim(fmt), &
               '    After convergence mixer',trim(mixers(-m%n_itt)%name)
       else if ( m%n_itt > 0 ) then
          write(*,'(2a,t50,''= '',i0)') trim(fmt), &
               '    Number of mixing iterations',m%n_itt
       end if

       if ( associated(m%next) .and. m%n_itt > 0 ) then
          write(*,'(2a,t50,''= '',a)') trim(fmt), &
               '    Following mixing method',trim(m%next%name)
       end if
          
    end do

  end subroutine mixing_print


  ! Calculate the inverse of a matrix
  subroutine inverse(n, A, B, info )

    integer, intent(in) :: n
    real(dp), intent(in)  :: A(n,n)
    real(dp), intent(out) :: B(n,n)
    integer, intent(out) :: info

    integer :: i, j
    ! Local arrays
    real(dp) :: pm(n,n), work(n*4), err
    ! Relative tolerance dependent on the magnitude
    real(dp), parameter :: etol = 1.e-6_dp
    integer :: ipiv(n)

    ! initialize info
    info = 0

    ! simple check and fast return
    if ( n == 1 ) then

       B(1,1) = 1._dp / A(1,1)

       return

    end if

    ! Copy over to that we can check the accuracy afterwards
    B = A

    call dgetrf(n,n,B,n,ipiv,info)
    if ( info /= 0 ) return

    call dgetri(n,B,n,ipiv,work,n*4,info)
    if ( info /= 0 ) return

    ! Check correcteness
    pm = matmul(A,B)

    do j = 1 , n 
       do i = 1 , n
          if ( i == j ) then
             err = pm(i,j) - 1._dp
          else
             err = pm(i,j)
          end if

          ! This is pretty strict tolerance!
          if ( abs(err) > etol ) then
             ! Signal failure in inversion
             info = - n - 1
             return
          end if

       end do
    end do

  end subroutine inverse


  ! Stack handling routines

  function stack_check(stack,n) result(check)
    type(Fstack_dData1D), intent(inout) :: stack
    integer, intent(in) :: n
    logical :: check

    ! Local arrays
    type(dData1D), pointer :: dD1

    if ( n_items(stack) == 0 ) then
       check = .true.
    else

       ! Check that the stack stored arrays are
       ! of same size...
       
       dD1 => get_pointer(stack,1)
       check = n == size(dD1)

    end if

  end function stack_check
    

  subroutine push_F(s_F,n,F)
    type(Fstack_dData1D), intent(inout) :: s_F
    integer, intent(in) :: n
    real(dp), intent(in) :: F(n)

    type(dData1D) :: dD1
    real(dp), pointer :: sF(:)
    integer :: in, ns

    if ( .not. stack_check(s_F,n) ) then
       call die('mixing: history has changed size...')
    end if

    in = n_items(s_F)
    ns = max_size(s_F)

    if ( in == ns ) then
       
       ! we have to cycle the storage
       call get(s_F,1,dD1)

    else

       call newdData1D(dD1, n, '(F)')

    end if

    sF => val(dD1)

    call dcopy(n,F,1,sF,1)

    ! Push the data to the stack
    call push(s_F,dD1)

    ! Delete double reference
    call delete(dD1)

  end subroutine push_F

  subroutine update_F(s_F,n,F)
    type(Fstack_dData1D), intent(inout) :: s_F
    integer, intent(in) :: n
    real(dp), intent(in) :: F(n)

    type(dData1D), pointer :: dD1
    real(dp), pointer :: FF(:)
    integer :: in

    if ( .not. stack_check(s_F,n) ) then
       call die('mixing: history has changed size...')
    end if

    in = n_items(s_F)

    if ( in == 0 ) then
       
       ! We need to add it as it does not exist
       call push_F(s_F,n,F)

    else

       ! we have an entry, update the latest
       dD1 => get_pointer(s_F,in)

       FF => val(dD1)

!$OMP parallel workshare default(shared)
       FF = F
!$OMP end parallel workshare

    end if

  end subroutine update_F

  subroutine update_res(s_res,n,oldF,newF)
    type(Fstack_dData1D), intent(inout) :: s_res
    integer, intent(in) :: n
    real(dp), intent(in) :: oldF(n), newF(n)

    type(dData1D), pointer :: dD1
    real(dp), pointer :: res(:)
    integer :: in

    if ( .not. stack_check(s_res,n) ) then
       call die('mixing: history has changed size...')
    end if

    in = n_items(s_res)

    if ( in == 0 ) then
       
       ! We need to add it as it does not exist
       call push_res(s_res,n,oldF,newF)

    else

       ! we have an entry, update the latest
       dD1 => get_pointer(s_res,in)

       res => val(dD1)

!$OMP parallel workshare default(shared)
       res = newF - oldF
!$OMP end parallel workshare

    end if

  end subroutine update_res

  subroutine push_res(s_res,n,oldF,newF)
    type(Fstack_dData1D), intent(inout) :: s_res
    integer, intent(in) :: n
    real(dp), intent(in) :: oldF(n), newF(n)

    type(dData1D) :: dD1
    real(dp), pointer :: res(:)
    integer :: in, ns

    if ( .not. stack_check(s_res,n) ) then
       call die('mixing: history has changed size...')
    end if

    in = n_items(s_res)
    ns = max_size(s_res)

    if ( in == ns ) then
       
       ! we have to cycle the storage
       call get(s_res,1,dD1)

    else

       call newdData1D(dD1, n, '(res)')

    end if

    res => val(dD1)

!$OMP parallel workshare default(shared)
    res = newF - oldF
!$OMP end parallel workshare

    ! Push the data to the stack
    call push(s_res,dD1)

    ! Delete double reference
    call delete(dD1)

  end subroutine push_res

  subroutine push_rres(s_rres,s_res)
    type(Fstack_dData1D), intent(inout) :: s_rres
    type(Fstack_dData1D), intent(in) :: s_res

    type(dData1D) :: dD1
    type(dData1D), pointer :: pD1
    real(dp), pointer :: res1(:), res2(:), rres(:)
    integer :: in, ns

    if ( n_items(s_res) < 2 ) then
       call die('mixing: Residual residuals cannot be calculated, &
            &inferior residual size.')
    end if

    in = n_items(s_res)

    ! First get the value of in
    pD1 => get_pointer(s_res,in-1)
    res1 => val(pD1)
    ! get the value of in
    pD1 => get_pointer(s_res,in)
    res2 => val(pD1)

    in = n_items(s_rres)
    ns = max_size(s_rres)

    if ( in == ns ) then
       
       ! we have to cycle the storage
       call get(s_rres,1,dD1)

    else

       call newdData1D(dD1, size(res1), '(res)')

    end if

    ! Get the residual of the residual
    rres => val(dD1)

!$OMP parallel workshare default(shared)
    rres = res2 - res1
!$OMP end parallel workshare

    ! Push the data to the stack
    call push(s_rres,dD1)

    ! Delete double reference
    call delete(dD1)

  end subroutine push_rres

end module m_mixing


! Also the mixing container
module m_mixing_scf

  use class_Fstack_dData1D
  use m_mixing, only: tMixer
  use m_mixing, only: mixing_reset, mixing_history_clear

  implicit none

  private
  save

  type(tMixer), allocatable, target :: scf_mixs(:)
  type(tMixer), pointer :: scf_mix => null()

  public :: scf_mixs, scf_mix

  public :: reset_mixing_scf
  public :: mixing_scf_converged
  public :: mixing_scf_history_clear

contains

  subroutine mixing_scf_converged( SCFconverged )
    use parallel, only: IONode
    logical, intent(inout) :: SCFconverged
    integer :: i

    ! Return if no convergence
    if ( .not. SCFconverged ) return

    if ( scf_mix%n_itt < 0 ) then
       ! this means that we skip to the 
       ! following algorithm
       scf_mix => scf_mixs(-scf_mix%n_itt)
       SCFconverged = .false.

       if ( allocated(scf_mix%stack) ) then

          do i = 1 , size(scf_mix%stack)
             ! delete all but one history
             ! This should be fine
             call reset(scf_mix%stack(i), -1)
          end do

       end if

       if ( IONode ) then
         write(*,'(a)') ':!: SCF cycle continued for different mixer'
       end if

    end if

  end subroutine mixing_scf_converged

  subroutine reset_mixing_scf()

    nullify(scf_mix)
    call mixing_reset( scf_mixs )

  end subroutine reset_mixing_scf

  subroutine mixing_scf_history_clear( )
    
    call mixing_history_clear( scf_mixs )
    scf_mix => scf_mixs(1)
    
  end subroutine mixing_scf_history_clear

end module m_mixing_scf
  
