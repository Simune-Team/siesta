! Module for all mixing methods in a standard way

! This module implements mixing of the Pulay and Broyden
! type.

! The Pulay method is implemented in the fast calculation
! setup and in the stable method.
! The stable method is executed if the inversion fails.

! The Broyden scheme is implemented

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
    integer :: n_hist, n_kick
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
       n_pulay_orig = fdf_get('DM.NumberPulay',2)
       n_broy_orig = fdf_get('DM.NumberBroyden',0)
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
       allocate(mixers(nm))
       mixers(:)%w = w
       mixers(:)%n_hist = n_hist

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
       
       allocate(mixers(nm))
       mixers(:)%w = w
       mixers(:)%n_hist = n_hist

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
         if ( leqi(str,'pulay') .or. leqi(str,'original') ) then
            ! "non" stable version, might 
            ! error on inversion
            v = 0
         else if ( leqi(str,'kresse') .or. leqi(str,'stable') ) then
            ! stable version, will nearly always succeed on inversion
            v = 1
         end if
      case ( MIX_BROYDEN ) 
         if ( leqi(str,'broyden') .or. leqi(str,'original') ) then
            ! cycle the history
            v = 0
         else if ( leqi(str,'restart') ) then
            v = 1
         end if
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

  end subroutine mixing_init


  subroutine mixing_1d( mix, iscf, n, oldF, newF )
    ! The current mixing method 
    type(tMixer), pointer :: mix
    integer, intent(in) :: iscf, n
    real(dp), intent(inout) :: oldF(n), newF(n)

    integer :: info, o_ns
    real(dp) :: rtmp

    if ( mix%m == MIX_BROYDEN ) then
       ! we need to check for first linear thing
       o_ns = n_items(mix%stack(1))
    end if

    ! Add history (if needed)
    call mix_add_history( mix, iscf, n, oldF, newF )

    ! Select mixer
    select case ( mix%m )

    case ( MIX_LINEAR )

       if ( debug_mix ) write(*,'(2a)') trim(debug_msg),' linear'

       call mixing_linear(mix, n, oldF, newF )

    case ( MIX_PULAY )
       
       if ( debug_mix ) then
          if ( mix%v == 0 ) then
             write(*,'(2a)') trim(debug_msg),' Pulay'
          else
             write(*,'(2a)') trim(debug_msg),' Pulay, stable'
          end if
       end if

       call mixing_pulay(mix, n, oldF, newF , info )

       if ( info /= 0 .and. mix%v == 0 ) then

          ! Switch to the stable Pulay mixer
          mix%v = 1

          if ( debug_mix ) &
               write(*,'(2a)') trim(debug_msg), &
               ' Pulay -- failed, trying stable'

          ! Do stable Pulay mixing (if this fails it will revert
          ! to linear mixing)
          call mixing_pulay(mix, n, oldF, newF , info )

          ! swap back to fast method
          mix%v = 0

       end if

       if ( info /= 0 ) then

          ! All Pulay mixers failed, run linear
          if ( IONode ) then
             write(*,'(a)') 'mix: Pulay stable -- failed, > linear'
          end if

          rtmp = mix%w
          mix%w = mix%rv(1)
          call mixing_linear(mix, n, oldF, newF )
          mix%w = rtmp

       end if

    case ( MIX_BROYDEN )

       if ( o_ns == 0 ) then

          if ( debug_mix ) write(*,'(2a)') &
               trim(debug_msg),' Broyden -- init'

          call mixing_linear(mix, n, oldF, newF)

          ! Signal no error
          info = 0

       else

          if ( debug_mix ) then
            if ( mix%v == 0 ) then
              write(*,'(2a)') trim(debug_msg),' Broyden'
            else
              write(*,'(2a)') trim(debug_msg),' Broyden, restart'
            end if
          end if

          call mixing_broyden(mix, n, oldF, newF, info)

       end if

       if ( info /= 0 ) then

          if ( IONode ) then
             write(*,'(a)')'mix: Broyden -- failed, > linear'
          end if
          
          call mixing_linear(mix, n, oldF, newF)
          
       end if
       
    end select

    ! Copy over new function
!$OMP parallel workshare default(shared)
    oldF = newF
!$OMP end parallel workshare

    ! Step iterator
    mix%cur_itt = mix%cur_itt + 1
    
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

  subroutine mix_add_history( mix, iscf, n, oldF, newF )
    ! The current mixing scheme
    type(tMixer), pointer :: mix
    integer, intent(in) :: iscf
    integer, intent(in) :: n
    real(dp), intent(inout) :: oldF(n), newF(n)

    select case ( mix%m ) 
    case ( MIX_LINEAR )
       if ( is_next(mix,MIX_PULAY) ) then
          call mix_add_history_pulay(mix, iscf, n, oldF, newF)
       else if ( is_next(mix,MIX_BROYDEN) ) then
          call mix_add_history_broyden(mix, iscf, n, oldF, newF, &
               .true.)
       end if
    case ( MIX_PULAY )
       call mix_add_history_pulay(mix, iscf, n, oldF, newF)
    case ( MIX_BROYDEN )
       call mix_add_history_broyden(mix, iscf, n, oldF, newF, .false.)
    end select

  end subroutine mix_add_history

  subroutine mix_add_history_pulay( mix, iscf, n, oldF, newF )
    ! The current mixing scheme
    type(tMixer), pointer :: mix
    integer, intent(in) :: iscf
    integer, intent(in) :: n
    real(dp), intent(inout) :: oldF(n), newF(n)

    ! Temporary fields
    type(dData1D) :: a1D, b1D
    real(dp), pointer :: told(:) , tres(:)
    
    integer :: in, ni, nm

    ni = n_items(mix%stack(1))
    nm = max_size(mix%stack(1))

    ! if the history is 0 we need to
    ! initialize the size of the array
    if ( ni > 0 ) then
       in = getnsize(mix)
       if ( in /= n ) then
          call die('mixing: Number of elements has changed &
               &before resetting mixing arrays.')
       end if
    end if

    ! If the history is full, simply request the first one
    if ( ni == nm ) then
       
       ! Retrieve first item (copy it so that when we
       ! push we do not actually delete it)
       
       call get(mix%stack(1),1,a1D)
       call get(mix%stack(2),1,b1D)

    else

       ! Create new history
       call newdData1D(a1D, n, '(in)')
       call newdData1D(b1D, n, '(res)')

       ni = ni + 1

    end if

    ! Get old and residue
    told => val(a1D)
    tres => val(b1D)

    ! push to stack
    call push(mix%stack(1), a1D)
    call delete(a1D)
    call push(mix%stack(2), b1D)
    call delete(b1D)

    if ( debug_mix ) &
         write(*,'(a,2(a,i0))') trim(debug_msg), &
         ' n_hist = ',ni, ' / ',nm

    ! Store input and residual
!$OMP parallel do default(shared), private(in)
    do in = 1 , n
       told(in) = oldF(in)
       tres(in) = newF(in) - oldF(in)
    end do
!$OMP end parallel do

  end subroutine mix_add_history_pulay

  subroutine mix_add_history_broyden( mix, iscf, n, oldF, newF , &
       overwrite )
    ! The current mixing scheme
    type(tMixer), pointer :: mix
    integer, intent(in) :: iscf
    integer, intent(in) :: n
    real(dp), intent(inout) :: oldF(n), newF(n)
    ! In case we do linear interpolation before Broyden
    ! we can still re-use the last iteration from the linear
    ! mixing. In this case we allow overwriting for n_items == 1 and overwrite
    logical, intent(in) :: overwrite

    ! Temporary fields
    type(dData1D) :: a1D, b1D
    real(dp), pointer :: tres(:) , tJac(:)
    
    integer :: ni, nm
    real(dp) :: jinv0

    jinv0 = mix%w
    ni = n_items(mix%stack(1))
    nm = max_size(mix%stack(1))

    if ( ni > 0 .and. .not. overwrite ) return

    if ( ni == 0 ) then
    
       ! Create new history
       call newdData1D(a1D, n, '(res)')
       call newdData1D(b1D, n, '(Jacobian-res)')

       ! push to stack
       call push(mix%stack(1), a1D)
       call push(mix%stack(2), b1D)
       
       ni = ni + 1

    else

       ! get pointers
       call get(mix%stack(1),1,a1D)
       call get(mix%stack(2),1,b1D)

    end if

    ! Get residual and Jacobian residual
    tres => val(a1D)
    tJac => val(b1D)

    ! push to stack
    call delete(a1D)
    call delete(b1D)

    if ( debug_mix ) &
         write(*,'(a,2(a,i0))') trim(debug_msg), &
         ' n_hist = ',ni, ' / ',nm
    
    ! Store input and residual
!$OMP parallel workshare default(shared)
    tres = newF - oldF
    tJac = jinv0 * tres
!$OMP end parallel workshare

    ! The Broyden scheme updates the history
    ! subsequently...

  end subroutine mix_add_history_broyden


  ! Returns the value array from the stack(:)
  ! Returns this array:
  !    mix%stack(sidx)(hidx)
  function getstackval(mix,sidx,hidx) result(d1)
    type(tMixer), intent(in) :: mix
    integer, intent(in) :: sidx, hidx
    real(dp), pointer :: d1(:)

    ! Local arrays
    type(dData1D), pointer :: dD1

    dD1 => get_pointer(mix%stack(sidx),hidx)
    d1 => val(dD1)
    
  end function getstackval


  function getnsize(mix) result(n)
    type(tMixer), intent(in) :: mix
    integer :: n
    
    ! Local arrays
    type(dData1D), pointer :: dD1

    dD1 => get_pointer(mix%stack(1),1)
    n = size(dD1)

  end function getnsize

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
       ! if the following one is broyden
       ! clean history.
       if ( is_next(mix,MIX_BROYDEN) ) then
          ! delete all pulay history
          call reset(mix%stack(1))
          call reset(mix%stack(2))
       end if
       
    case ( MIX_BROYDEN )
       
       ! if the following one is broyden
       ! clean history.
       if ( is_next(mix,MIX_LINEAR) ) then
          ! delete all but one BROYDEN history
          call reset(mix%stack(1), -1)
          call reset(mix%stack(2), -1)
       else if ( is_next(mix,MIX_PULAY) ) then
          ! delete all BROYDEN history
          call reset(mix%stack(1))
          call reset(mix%stack(2))
       end if
       
    end select

    if ( associated(mix%next) ) then
       mix => mix%next
       mix%cur_itt = 0
       if ( debug_mix ) write(*,'(2a)') &
            trim(debug_msg),' switching mixer...'
    end if

  end subroutine mixing_step
    

  ! Here the actual mixing methods will be employed
  subroutine mixing_linear( mix, n, oldF, newF)

    ! The current mixing method 
    type(tMixer), intent(inout) :: mix

    integer, intent(in) :: n
    real(dp), intent(inout) :: oldF(n), newF(n)

    real(dp) :: alpha
    integer :: in

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
  subroutine mixing_pulay( mix, n, oldF, newF, info)
    ! The current mixing method 
    type(tMixer), intent(inout) :: mix

    ! Input/output arrays
    integer, intent(in) :: n
    real(dp), intent(inout) :: oldF(n), newF(n)

    ! run-time information
    integer, intent(out) :: info

    ! Temporary arrays for local data structures
    real(dp), dimension(:), pointer :: F
    real(dp), dimension(:), pointer :: res1, res2
    real(dp), dimension(:), pointer :: rres1, rres2

    integer :: ns, nb
    integer :: i, j, in
    ! Used arrays
    real(dp), allocatable :: alpha(:), b(:,:), bi(:,:)
    real(dp) :: G, ssum

#ifdef MPI
    integer :: MPIerror
#endif

    real(dp), external :: ddot

    ! Denote no warning
    info = 0

    ! The Pulay mixing variable is called G
    G = mix%w

    ! Number of saved histories, so far
    ns = n_items(mix%stack(1))

    if ( ns == 1 ) then

       ssum = mix%w
       mix%w = mix%rv(1)
       call mixing_linear( mix, n, oldF, newF )
       mix%w = ssum

       return

    end if

    select case ( mix%v )
    case ( 0 ) ! Default (fast) pulay mixing
       nb = ns + 1
       i  = ns
    case ( 1 ) ! Stable (slow) pulay mixing
       nb = ns - 1
       i  = ns - 1
    case default
       call die('mixing: Pulay unknown variant')
    end select

    ! Allocate arrays
    allocate(b(nb,nb),bi(nb,nb),alpha(i))

    ! Create A_ij coefficients for inversion
    select case ( mix%v )
    case ( 0 ) ! Default (fast) pulay mixing

       do i = 1 , ns

          ! Get i'th residual array
          res1 => getstackval(mix,2,i)

          do j = 1 , i

             ! Get j'th residual array
             res2 => getstackval(mix,2,j)

             ! B(i,j) = B(j,i) = dot_product(Res(i)*Res(j))
             ssum = ddot(n,res1(1),1,res2(1),1)

             b(i,j) = ssum
             b(j,i) = ssum

          end do

          ! Now extend the matrix with ones in an extra colum
          ! and row ...
          b(i,nb) = 1.0_dp
          b(nb,i) = 1.0_dp

       end do

       ! ... except in the extra diagonal entry
       b(nb,nb) = 0.0_dp

    case ( 1 ) ! Stable (slow) pulay mixing

       do i = 1 , nb

          ! Get i'th residual array
          res1 => getstackval(mix,2,i)
          rres1 => getstackval(mix,2,i+1)

          do j = 1 , i

             ! Get j'th residual array
             res2 => getstackval(mix,2,j)
             rres2 => getstackval(mix,2,j+1)

             ! B(i,j) = B(j,i) = dot_product(Delta Res(i)*Delta Res(j))
             ssum = 0._dp
!$OMP parallel do default(shared), private(in), reduction(+:ssum)
             do in = 1 , n
                ssum = ssum + (rres1(in)-res1(in)) * &
                     (rres2(in)-res2(in))
             end do
!$OMP end parallel do

             b(i,j) = ssum
             b(j,i) = ssum
             
          end do
          
       end do
       
    end select

#ifdef MPI
    ! Global operations, but only for the non-extended entries
    call MPI_AllReduce(b(1,1),bi(1,1),nb*nb, &
         MPI_double_precision, MPI_Sum, &
         MPI_Comm_World,MPIerror)
    ! copy over reduced arrays
    b = bi
    if ( mix%v == 0 ) then
       b(:,nb) = 1._dp
       b(nb,:) = 1._dp
       b(nb,nb) = 0._dp
    end if
#endif

    ! Get inverse of matrix
    call inverse(nb, b, bi, info)
    if ( info /= 0 ) then
     
       ! Deallocate local arrays
       deallocate(alpha,b,bi)

       ! return immediately, the algorithm
       ! will try the stable one, then the linear one

       return
       
    else ! inversion succeeded


    select case ( mix%v )
       
    case ( 0 ) ! Default (fast) pulay mixing
       
       do i = 1 , ns
          alpha(i) = bi(i,nb)
       end do
       
    case ( 1 ) ! Stable (slow) pulay mixing

       ! Get current residual
       res2 => getstackval(mix,2,ns)

       do i = 1 , nb
          
          ! Calculate the coefficients on all processors
          alpha(i) = 0._dp
             
          do j = 1 , nb
             
             ! Get j'th residual array
             res1 => getstackval(mix,2,j)
             rres1 => getstackval(mix,2,j+1)
             
             ssum = 0._dp
!$OMP parallel do default(shared), &
!$OMP&private(in), reduction(+:ssum)
             do in = 1 , n
                ssum = ssum + (rres1(in)-res1(in)) * res2(in)
             end do
!$OMP end parallel do
             
             alpha(i) = alpha(i) - bi(i,j) * ssum
             
          end do
          
       end do
       
#ifdef MPI
       ! Reduce the alpha
       call MPI_AllReduce(alpha(1),b(1,1),nb, &
            MPI_double_precision, MPI_Sum, &
            MPI_Comm_World,MPIerror)
       alpha(:) = b(:,1)
#endif

    end select

    end if

    ! if debugging print out the different variables
    if ( debug_mix ) then
       write(*,'(2a,f10.6,a,100(tr1,e10.4))') &
            trim(debug_msg),' G = ',G,', alpha = ',alpha
    end if

    ! Calculate new function
    select case ( mix%v )

    case ( 0 ) ! Default (fast) pulay mixing

       F => getstackval(mix,1,1)
       res1 => getstackval(mix,2,1)

!$OMP parallel workshare default(shared)
       newF(:) = alpha(1) * ( F + G * res1 )
!$OMP end parallel workshare
       
       do i = 2 , ns
          
          F => getstackval(mix,1,i)
          res1 => getstackval(mix,2,i)
          
!$OMP parallel workshare default(shared)
          newF = newF + alpha(i) * ( F + G * res1 )
!$OMP end parallel workshare
       
       end do
       
    case ( 1 ) ! Stable (slow) pulay mixing

       ! Read former matrices for mixing .........
       res1 => getstackval(mix,2,ns)

       ! Copy over input dm, and add the linear mixing
!$OMP parallel workshare default(shared)
       newF = oldF + G * res1
!$OMP end parallel workshare
       
       do i = 1 , nb

          ! two consecutive functions and residuals
          res1 => getstackval(mix,1,i)
          rres1 => getstackval(mix,2,i)
          res2 => getstackval(mix,1,i+1)
          rres2 => getstackval(mix,2,i+1)

!$OMP parallel workshare default(shared)
          newF = newF + alpha(i) * &
               ( res2 - res1 + G * ( rres2 - rres1 ) )
!$OMP end parallel workshare
          
       end do

    end select
    
    ! Deallocate local arrays
    deallocate(alpha,b,bi)
    
  end subroutine mixing_pulay


  ! Broyden mixing
  subroutine mixing_broyden( mix, n, oldF, newF, info)
    ! The current mixing method 
    type(tMixer), intent(inout) :: mix

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

    if ( ns == 0 ) then

       ! Signal an error
       info = -1001

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

       if ( mix%v == 0 ) then ! Cycling history

          ! get the data
          call get(mix%stack(1),1,a1D)
          call get(mix%stack(2),1,b1D)

       else if ( mix%v == 1 ) then ! Restart history
          
          call reset(mix%stack(1))
          call reset(mix%stack(2))

          ! Create new history
          call newdData1D(a1D, n, '(res)')
          call newdData1D(b1D, n, '(Jacobian-res)')

          ns = 1

       end if

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

    if ( debug_mix ) &
         write(*,'(a,2(a,i0))') trim(debug_msg), &
         ' n_hist = ',ns, ' / ',nm

    deallocate(F)

  end subroutine mixing_broyden

  ! Deletes all history, 
  ! Possibly reassign number of new histories saved
  subroutine mixing_history_clear( mixers )
    type(tMixer), intent(inout), target :: mixers(:)

    type(tMixer), pointer :: m, next
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

       case ( MIX_PULAY, MIX_BROYDEN )

          allocate(m%stack(2))

          ! allocate Fin, Pulay
          ! allocate Res[i], Broyden
          call new(m%stack(1), m%n_hist)
          ! allocate Res[i], Pulay
          ! allocate 'u', Broyden
          call new(m%stack(2), m%n_hist)

          if ( m%n_itt < 0 ) then

             next => mixers(-m%n_itt)
             ! Also allocate the next one
             if ( .not. allocated(next%stack) ) then
                allocate(next%stack(2))
                
                next%stack(1) = m%stack(1)
                next%stack(2) = m%stack(2)

             end if
             
          end if

       end select

    end do

    ! Point the simpler methods to support the later
    ! versions
    do im = 1 , size(mixers)
       m => mixers(im)

       ! In certain cases this can already be
       ! allocated n_itt < 0
       if ( allocated(m%stack) ) cycle
       
       select case ( m%m ) 
          
       case ( MIX_LINEAR )
          
          if ( is_next(m,MIX_PULAY,next=next) ) then
             
             allocate(m%stack(2))

             m%stack(1) = next%stack(1)
             m%stack(2) = next%stack(2)

          else if ( is_next(m,MIX_BROYDEN,next=next) ) then
             
             allocate(m%stack(2))
             
             m%stack(1) = next%stack(1)
             m%stack(2) = next%stack(2)
             
          end if

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
                  '    Variant','original'
          else if ( m%v == 1 ) then
             write(*,'(2a,t50,''= '',a)') trim(fmt), &
                  '    Variant','stable'
          end if

          write(*,'(2a,t50,''= '',i0)') trim(fmt), &
               '    History steps',m%n_hist
          write(*,'(2a,t50,''= '',f12.6)') trim(fmt), &
               '    Initial linear mixing weight',m%rv(1)
          write(*,'(2a,t50,''= '',f12.6)') trim(fmt), &
               '    Damping',m%w

       case ( MIX_BROYDEN )
          
          write(*,'(2a,t50,''= '',a)') trim(fmt), &
               ' Broyden mixing',trim(m%name)

          if ( m%v == 0 ) then
             write(*,'(2a,t50,''= '',a)') trim(fmt), &
                  '    Variant','original'
          else if ( m%v == 1 ) then
             write(*,'(2a,t50,''= '',a)') trim(fmt), &
                  '    Variant','restart'
          end if

          write(*,'(2a,t50,''= '',i0)') trim(fmt), &
               '    History steps',m%n_hist - 1
          write(*,'(2a,t50,''= '',f12.6)') trim(fmt), &
               '    Jacobian weight',m%w
          write(*,'(2a,t50,''= '',f12.6)') trim(fmt), &
               '    Weight prime',m%rv(1)
          
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
          write(*,'(2a,t50,''= '',i0)') trim(fmt), &
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
  
