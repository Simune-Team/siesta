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
! This code segment has been fully created by:
! Nick Papior Andersen, 2012, nickpapior@gmail.com
!

module m_ts_io_contour
!
! Routines that are used to read in and print out the contour for integration of the GFs
! 
  use precision, only : dp

  implicit none

  character(len=200), parameter :: OPT_N   = '(''ts_options: '',a)'
  character(len=200), parameter :: OPT_C   = '(''ts_options: '',a,t53,''=    '',a)'
  character(len=200), parameter :: OPT_F   = '(''ts_options: '',a,t53,''='',f10.4)'
  character(len=200), parameter :: OPT_INT = '(''ts_options: '',a,t53,''='',i5)'
  character(len=200), parameter :: OPT_F_U = '(''ts_options: '',a,t53,''='',f10.4,tr1,a)'
  character(len=200), parameter :: OPT_G_U = '(''ts_options: '',a,t53,''='',g11.4,tr1,a)'

contains

!  subroutine ts_read_contour_options(prefix, Elecs, N_poles, kT, IsVolt, V)
!
!    use parallel, only : IONode, Nodes, operator(.parcount.)
!    use fdf
!
!    character(len=*), intent(in) :: prefix
!    type(Elec), intent(in) :: Elecs(:)
!    integer, intent(in) :: N_poles
!    real(dp), intent(in) :: kT, V ! in Ry
!    logical :: IsVolt
!    character(len=C_N_NAME_LEN), allocatable :: tmp(:), nContours(:)
!
!    ! local variables
!    integer :: i,j,k,N, N_bias
!    integer :: different_poles
!
!    ! figure out how many different poles we need for the "fake" c_io types which contains the poles
!    different_poles = 1 ! this is the last electrodes mu
!    do i = 1 , size(Elecs) - 1
!       ! we only count the different mu-levels for the last electrode
!       if ( count(abs(Elecs(i+1:)%mu - Elecs(i)%mu) < 1.e-10_dp) == 0 ) then
!          different_poles = different_poles + 1
!       end if
!    end do
!
!    ! Count the number of equilibrium contour segments
!    N = 0
!    do i = 1 , size(Elecs)
!       N = N + Eq_segs(Elecs(i))
!    end do
!    if ( N == 0 ) then
!       call die('You cannot have assigned any contours for the electrodes.')
!    end if
!
!    ! collect all equilibrium names
!    allocate(tmp(N))
!    tmp(:) = ' '
!    N = 0
!    do i = 1 , size(Elecs)
!       do j = 1 , Eq_segs(Elecs(i))
!          N = N + 1
!          tmp(N) = Elecs(i)%Eq_seg(j)
!       end do
!    end do
!
!    ! find all unique equilibrium names
!    N = 0
!    uniq_names: do i = 2 , size(tmp)
!       j = 0
!       do 
!          j = j + 1
!          if ( i <= j ) exit
!          if ( tmp(i) == tmp(j) ) cycle uniq_names
!       end do
!       N = N + 1
!    end do uniq_names
!
!    ! Allocate space for the equilibrium contours
!    N_eq_io = N + different_poles
!    allocate(eq_io(N_eq_io))
!    allocate(nContours(N_eq_io))
!    
!    ! populate unique equilibrium names
!    N = 0
!    uniq_names_pop: do i = 1 , size(tmp)
!       j = 0
!       do 
!          j = j + 1
!          if ( j > N ) exit
!          if ( tmp(i) == nContours(j) ) cycle uniq_names_pop
!       end do
!       N = N + 1
!       if ( N > size(nContours) ) call die('Unique contour setup went wrong, contact devs')
!       nContours(N) = tmp(i)
!    end do uniq_names_pop
!    if ( N /= size(nContours) ) then
!       call die('ERROR: We have not populated enough contours')
!    end if
!    deallocate(tmp)
!
!    ! read in the equilibrium contours
!    do i = 1 , N
!       
!       ! read in the contour
!       call ts_read_contour_block(prefix,'',nContours(i),eq_io(i), kT,V)
!
!    end do
!    deallocate(nContours)
!
!    ! We here create the "fake" pole contours
!    j = 0
!    k = 0
!    do i = N + 1, N + different_poles
!       ! assign name to the eq_io
!       write(eq_io(i)%name,'(a,i0)') 'pole-',i-N
!       eq_io(i)%N = N_poles
!       eq_io(i)%type = 'eq'
!       eq_io(i)%part = 'pole'
!       c_io(i)%method = 'residual'
!       ! find the next mu level
!       do 
!          j = j + 1
!          if ( j > different_poles ) call die('Pole setup creation went wrong, contact devs')
!          ! we only count the different mu-levels for the last electrode
!          if ( j == different_poles ) then
!             k = k + 1
!             c_io(i)%a = Elecs(j)%mu
!             c_io(i)%b = Elecs(j)%mu
!             c_io(i)%d = Elecs(j)%mu
!             exit
!          else if ( count(abs(Elecs(j+1:)%mu - Elecs(j)%mu) < 1.e-10_dp) == 0 ) then
!             k = k + 1
!             c_io(i)%a = Elecs(j)%mu
!             c_io(i)%b = Elecs(j)%mu
!             c_io(i)%d = Elecs(j)%mu
!             exit
!          end if
!       end do
!       
!    end do
!
!    if ( k < different_poles ) then
!       call die('Pole setup went wrong, please contact the developers')
!    end if
!
!    ! we have now read in all information regarding the equilbrium contour
!
!    ! find number of integrals in the non-equilibrium contour
!    if ( IsVolt ) then
!
!       
!    end if
!
!    write(*,*) 'TODO suggest for the user the "best" offset in the fermi tails based &
!         &on the two closests fermi-levels'
!    
!    ! At this point we have read in all information regarding
!    ! the contours
!
!    do i = 1 , size(Elecs)
!
!       j = 1
!       cur = ts_c_io_get_index(eq_io,Elecs(i)%Eq_seg(j))
!       if ( Eq_segs(Elecs(i)) > 1 ) then
!          next = ts_c_io_get_index(eq_io,Elecs(i)%Eq_seg(j+1))
!          call ts_fix_contour(IsVolt, eq_io(cur),next=eq_io(next))
!       end if
!          
!       do j = 2 , Eq_segs(Elecs(i)) - 1
!          prev = ts_c_io_get_index(eq_io,Elecs(i)%Eq_seg(j-1))
!          cur = ts_c_io_get_index(eq_io,Elecs(i)%Eq_seg(j))
!          next = ts_c_io_get_index(eq_io,Elecs(i)%Eq_seg(j+1))
!          call ts_fix_contour(IsVolt, eq_io(cur), &
!               prev=eq_io(prev), next=eq_io(next))
!       end do
!
!       j = Eq_segs(Elecs(i))
!       cur = ts_c_io_get_index(eq_io,Elecs(i)%Eq_seg(j))
!       if ( Eq_segs(Elecs(i)) > 1 ) then
!          prev = ts_c_io_get_index(eq_io,Elecs(i)%Eq_seg(j-1))
!          call ts_fix_contour(IsVolt, eq_io(cur),prev=eq_io(prev))
!       end if
!
!       write(*,*)'TODO need to check the circle-line-tail sequence'
!    end do
!
!    if ( IsVolt ) then
!       if ( N_neq_io == 1 ) then
!          call ts_fix_contour(IsVolt, neq_io(1))
!       end if
!       do i = 1 , N_neq_io - 1
!          if ( 1 < i ) then
!             call ts_fix_contour(IsVolt, neq_io(i),prev=neq_io(i-1), next=neq_io(i+1))
!          else
!             call ts_fix_contour(IsVolt, neq_io(i),next=neq_io(i+1))
!          end if
!       end do
!       call ts_fix_contour(IsVolt, neq_io(N_neq_io),prev=neq_io(N_neq_io-1))
!
!       if ( N_neq_tail_io == 1 ) then
!          call ts_fix_contour(IsVolt, neq_tail_io(1))
!       end if
!       do i = 1 , N_neq_tail_io - 1
!          if ( 1 < i ) then
!             call ts_fix_contour(IsVolt, neq_tail_io(i), &
!                  prev=neq_tail_io(i-1), next=neq_tail_io(i+1))
!          else
!             call ts_fix_contour(IsVolt, neq_tail_io(i), &
!                  next=neq_tail_io(i+1))
!          end if
!       end do
!       call ts_fix_contour(IsVolt, neq_tail_io(N_neq_tail_io), &
!            prev=neq_tail_io(N_neq_tail_io-1))
!       
!    end if
!
!    ! at this point we have created all values in the contours
!    ! and made them available for creation
!
!  end subroutine ts_read_contour_options


  ! *****
  ! This routine fixes the inputs for the contours according to those given by 
  ! the input electrode
  ! 1.) It fixes the bounds next to each other if they have designated
  !     'next' or 'previous'
  ! 2.) Set the method to be equilibrium/non-equilibrium/transport
  ! 3.) Calculate number of points if dE specified
  ! 4.) Checks whether the contours are connected
  subroutine ts_fix_contour(cur,next,prev)
    use m_ts_io_ctype
    use fdf, only : leqi
    type(ts_c_io), intent(inout) :: cur
    type(ts_c_io), intent(inout), optional :: next, prev

    ! local variables
    real(dp) :: val

    ! we need this to "look ahead"
    if ( present(next) ) then
       if ( leqi(next%ca,'prev') .or. leqi(next%ca,'previous') ) then
          next%a = cur%b
          if ( leqi(cur%cb,'next') ) then
             call die('Connecting two contours by next and prev is invalid. &
                  &An explicit value is needed.')
          end if
       end if
    end if

    if ( leqi(cur%cb,'next') ) then
       if ( present(next) ) then
          cur%b = next%a
       else
          call die('The contour segment is not &
               &attached to a following segment.')
       end if
    end if

    if ( leqi(cur%ca,'prev') .or. leqi(cur%ca,'previous') ) then
       if ( present(prev) ) then
          cur%a = prev%b
       else
          call die('The contour is not fully connected')
       end if
    end if

    ! we can compare bounds
    if ( present(prev) ) then
       if (abs(cur%a - prev%b) > 1.e-8_dp ) then
          call die('Contour: '//trim(prev%name)//' and '//trim(cur%name)// &
               ' are not connected.')
       end if
    end if
    
    if ( present(next) ) then
       if (abs(cur%b - next%a) > 1.e-8_dp ) then
          call die('Contour: '//trim(cur%name)//' and '//trim(next%name)// &
               ' are not connected.')
       end if
    end if
    
    ! at this point both boundaries MUST exist
    if ( len_trim(cur%cd) > 0 .and. len_trim(cur%cN) == 0 ) then
       write(*,*)'TODO check limits for infinity'
       cur%N = nint(abs(cur%b - cur%a)/cur%d)
    end if
         
  end subroutine ts_fix_contour

!  ! *****
!  ! This routine only asserts that the electrodes contour parts are:
!  !  1.) Connected
!  !  2.) Equilibrium starts with circle contours
!  !  3.) The contours have the assigned part as said by the user
!  !  4.) The connected contours cannot both have "next" "prev" for a connecting part
!  !  5.) The tail parts are checked
!  subroutine ts_assert_contour_Elec(El,IsVolt)
!    use m_ts_electype
!    use fdf, only : leqi
!    type(Elec), intent(in) :: El
!    logical, intent(in) :: IsVolt
!
!    ! local variables
!    integer :: i, idx
!    logical :: is_circle
!    ! we need to check that the contours here match
!
!    if ( .not. contour_connected(El%Eq_seg) ) then
!       call die('Equilibrium contours in electrode: '//trim(Name(El))//' cannot be connected &
!            &correctly, please check for inconsistencies in the input.')
!    end if
!    do i = 1 , Eq_segs(El)
!
!       idx = contour_exists(El%Eq_seg(i),'Equilibrium contour could not be found: &
!            &'//trim(El%Eq_seg(i)))
!
!       ! check that it really is an equilibrium contour
!       if ( .not. leqi(c_io(idx)%type,'eq') ) then
!          call die('Contour '//trim(c_io(idx)%name)//' is not an equilibrium &
!               &contour as has been assigned in electrode: '//trim(name(El)))
!       end if
!       
!       ! the first contours must be circles
!       if ( i == 1 .and. .not. leqi(c_io(idx)%part,'circle') ) then
!          call die('First equilibrium contour is not a circle contour.')
!       else if ( c_io(idx)%part == 'circle' .and. .not. is_circle ) then
!          call die('You cannot move back to a circle contour.')
!       end if
!       is_circle = leqi(c_io(idx)%part,'circle')
!       if ( is_circle ) then
!          ! ensure that the circle ends before the fermi level
!          if ( .not. c_io(idx)%b < El%mu ) then
!             call die('The circle contour ends after the Fermi-level. &
!                  &This will constitute an error-prone numerical integration. Please correct')
!          end if
!       end if
!
!       if ( i /= size(El%Eq_seg) .and. leqi(c_io(idx)%part,'tail') ) then
!          call die('Only the last part of an equilibrium contour integration &
!               &must be assigned as a tail integral.')
!       else if ( i == size(El%Eq_seg) .and. .not. leqi(c_io(idx)%part,'tail') ) then
!          call die('Last part of the equilibrium contour integration must &
!               &be a tail part.')
!       end if
!       
!    end do
!
!
!    if ( IsVolt ) then
!       if ( .not. contour_connected(El%nEq_seg) ) then
!          call die('Non-equilibrium contours in electrode: '//trim(Name(El))//' cannot be connected &
!               &correctly, please check for inconsistencies in the input.')
!       end if
!       do i = 1 , nEq_segs(El)
!          
!          idx = contour_exists(El%nEq_seg(i),'Non-equilibrium contour could not be found: &
!               &'//trim(El%nEq_seg(i))
!
!          ! check that it really is a non-equilibrium contour
!          if ( .not. leqi(c_io(idx)%type,'neq') ) then
!             call die('Contour '//trim(c_io(idx)%name)//' is not a non-equilibrium &
!                  &contour as has been assigned in electrode: '//trim(name(El)))
!          end if
!          
!          if ( (i /= 1 .and. i /= size(El%nEq_seg)) .and. leqi(c_io(idx)%part,'tail') ) then
!             call die('Only the first/last part of a non-equilibrium contour integration &
!                  &must be assigned as a tail integral.')
!          else if ( i == 1 .and. .not. leqi(c_io(idx)%part,'tail') ) then
!             call die('First part of the non-equilibrium contour integration must &
!                  &be a tail part.')
!          else if ( i == size(El%nEq_seg) .and. .not. leqi(c_io(idx)%part,'tail') ) then
!             call die('Last part of the non-equilibrium contour integration must &
!                  &be a tail part.')
!          end if
!
!       end do
!    end if
!
!    if ( .not. contour_connected(El%t_seg) ) then
!       call die('Transport contours in electrode: '//trim(Name(El))//' cannot be connected &
!            &correctly, please check for inconsistencies in the input.')
!    end if
!    do i = 1 , t_segs(El)
!
!       idx = contour_exists(El%t_seg(i),'Transport contour could not be found: &
!            &'//trim(El%t_seg(i))
!
!       ! check that it really is a transport contour
!       if ( .not. leqi(c_io(idx)%type,'t') ) then
!          call die('Contour '//trim(c_io(idx)%name)//' is not a transport &
!               &contour as has been assigned in electrode: '//trim(name(El)))
!       end if
!       
!    end do
!
!  contains 
!
!    function contour_exists(name,err) result(idx)
!      character(len=*), intent(in) :: name, err
!      integer :: idx
!      idx = ts_c_io_get_index(c_io,name)
!      if ( idx < 1 ) then
!         call die(err)
!      end if
!    end function contour_exists
!
!    function contour_connected(seg) result(val)
!      character(len=C_N_NAME_LEN), allocatable, intent(in) :: seg(:)
!      logical :: val
!      integer :: i, idx1, idx2
!      val = .true.
!      if ( .not. allocated(seg) ) return
!
!      idx1 = ts_c_io_get_index(c_io,seg(1))
!      if ( leqi(c_io(idx1)%ca,'prev') .or. &
!           leqi(c_io(idx1)%ca,'previous') ) then
!         call die('First contour cannot start with previous: '//trim(seg(1)))
!      end if
!
!      do i = 2 , size(seg)
!         idx1 = ts_c_io_get_index(c_io,seg(i-1))
!         idx2 = ts_c_io_get_index(c_io,seg(i))
!         ! we can compare bounds
!         val = val .and. ( abs(c_io(idx1)%b - c_io(idx2)%a) > 1.e-8_dp ) 
!
!         ! we check "next" "prev"
!         val = val .and. .not. ( leqi(c_io(idx1)%cb,'next') .and. ( &
!              leqi(c_io(idx2)%ca,'prev') .or. leqi(c_io(idx2)%ca,'previous') ) )
!         
!      end do
!      
!    end function contour_connected
!
!  end subroutine ts_assert_contour_Elec

  subroutine write_e(str,val)

    use units, only : eV
    character(len=*), intent(in) :: str
    real(dp), intent(in) :: val
    ! This is our definition of infinity....
    if ( abs(val) > 10000._dp ) then
       if ( val < 0._dp ) then
          write(*,opt_c) trim(str),'-Infinity'
       else
          write(*,opt_c) trim(str),' Infinity'
       end if
    else
       write(*,opt_f_u) trim(str),val / eV,'eV'
    end if
  end subroutine write_e



!  subroutine ts_print_contour_warnings(cEq,cnEq,kT,Eq_Eta, nEq_Eta, N_poles, IsVolt)
!
!    use parallel, only : IONode, Nodes, operator(.parcount.)
!    use units, only : Pi, Kelvin
!
!    type(ts_io_c), intent(in) :: cEq(:), cnEq(:)
!    real(dp), intent(in) :: kT, Eq_Eta, nEq_Eta
!    integer, intent(in) :: N_poles
!    logical, intent(in) :: IsVolt
!    integer :: i
!    real(dp) :: diff
!
!    if ( .not. IONode ) return
!
!    call calc_Eta_Pole_Diff_Kelvin(Eq_Eta,diff)
!
!    ! say we warn if we are 3 kT away from the pole
!    if ( abs(diff) <= 10._dp ) then
!       write(*,'(a)') 'NOTICE : Equilibrium Eta value is &
!            &less than 10 Kelvin away from a pole.'
!       write(0,'(a)') 'NOTICE : Equilibrium Eta value is &
!            &less than 10 Kelvin away from a pole.'
!    end if
!    if ( abs(diff) <= 5._dp ) then
!       write(*,'(a)') 'WARNING: Equilibrium Eta value is &
!            &less than  5 Kelvin away from a pole.'
!       write(0,'(a)') 'WARNING: Equilibrium Eta value is &
!            &less than  5 Kelvin away from a pole.'
!    end if
!
!    if ( IsVolt ) then
!
!       call calc_Eta_Pole_Diff_Kelvin(nEq_Eta,diff)
!
!       ! say we warn if we are 3 kT away from the pole
!       if ( abs(diff) <= 10._dp ) then
!          write(*,'(a)') 'NOTICE : non-Equilibrium Eta value is &
!               &less than 10 Kelvin away from a pole.'
!          write(0,'(a)') 'NOTICE : non-Equilibrium Eta value is &
!               &less than 10 Kelvin away from a pole.'
!       end if
!       if ( abs(diff) <= 5._dp ) then
!          write(*,'(a)') 'WARNING: non-Equilibrium Eta value is &
!               &less than  5 Kelvin away from a pole.'
!          write(0,'(a)') 'WARNING: non-Equilibrium Eta value is &
!               &less than  5 Kelvin away from a pole.'
!       end if
!
!       i = 2 * ( sum(cEq(:)%N) + N_poles )
!       if ( mod(i,Nodes) /= 0 ) then
!          write(*,*) "NOTICE: Equilibrium energy contour points are not"
!          write(*,*) "        divisable by the number of nodes."
!          write(*,*) "        Better scalability is achived by changing:"
!          write(*,*) "          - TS.Contour.Eq.Circle.N"
!          write(*,*) "          - TS.Contour.Eq.Line.N"
!          write(*,*) "          - TS.Contour.Eq.Pole.N"
!          write(*,*) "          - %block TS.Contour.Eq"
!
!          ! Calculate optimal number of energy points
!          write(*,'(t10,a,i4)') "Used equilibrium # of energy points   : ",i
!          i = Nodes .PARCOUNT. i
!          write(*,'(t10,a,i4,tr1,a4,i3,/)') &
!               "Optimal equilibrium # of energy points: ",i, &
!               "+- i*",Nodes
!       end if
!       
!       i = sum(cnEq(:)%N)
!       if ( mod(i,Nodes) /= 0 ) then
!          write(*,*) "NOTICE: Non-equilibrium energy contour points are not"
!          write(*,*) "        divisable by the number of nodes."
!          write(*,*) "        Better scalability is achieved by changing:"
!          write(*,*) "          - TS.Contour.nEq.N"
!          write(*,*) "          - %block TS.Contour.nEq"
!          
!          ! Calculate optimal number of energy points
!          write(*,'(t10,a,i4)') "Used non-equilibrium # of energy points   : ",i
!          i = Nodes .PARCOUNT. i
!          write(*,'(t10,a,i4,tr1,a4,i3,/)') &
!               "Optimal non-equilibrium # of energy points: ",i, &
!               "+- i*",Nodes
!       end if
!       
!       i = 2 * ( sum(cEq(:)%N) + N_poles ) + sum(cnEq(:)%N)
!       if ( mod(i,Nodes) /= 0 ) then
!          write(*,*) "NOTICE: Total energy contour points are not"
!          write(*,*) "        divisable by the number of nodes."
!          
!          ! Calculate optimal number of energy points
!          write(*,'(t10,a,i4)') "Used # of energy points   : ",i
!          i = Nodes .PARCOUNT. i
!          write(*,'(t10,a,i4,tr1,a4,i3,/)') &
!               "Optimal # of energy points: ",i,"+- i*",Nodes
!       end if
!    else
!       i = sum(cEq(:)%N) + N_poles
!       
!       !   - The equilibrium parts are the same computational cost
!       !   * Solution make the equi contours divisible by Nodes
!       if ( mod(i,Nodes) /= 0 ) then
!          write(*,*) "NOTICE: Equilibrium energy contour points are not"
!          write(*,*) "        divisable by the number of nodes."
!          write(*,*) "        Better scalability is achived by changing:"
!          write(*,*) "          - TS.Contour.Eq.Circle.N"
!          write(*,*) "          - TS.Contour.Eq.Line.N"
!          write(*,*) "          - TS.Contour.Eq.Pole.N"
!          write(*,*) "          - %block TS.Contour.Eq"
!
!          ! Calculate optimal number of energy points
!          write(*,'(t10,a,i4)') "Used # of energy points   : ",i
!          i = Nodes .PARCOUNT. i
!          write(*,'(t10,a,i4)') "Optimal # of energy points: ",i
!       end if
!
!    end if
!
!  contains
!
!    subroutine calc_Eta_Pole_Diff_Kelvin(Eta,diff)
!      real(dp), intent(in) :: Eta
!      real(dp), intent(out) :: diff
!      integer :: iPole
!
!      ! We check whether the equilibrium shift into the imaginary plane
!      ! is too close to a Fermi-pole
!      ! We define this to be the case if 
!      iPole = nint(Eta / (Pi * kT) )
!      if ( iPole < 0 ) call die('Error in the eta value')
!      if ( mod(iPole,2) == 0 ) then
!         if ( iPole == 0 ) then
!            ! we have the closest pole at kT * Pi
!            diff = kT * Pi - Eta
!         else
!            ! we have the closest pole at either
!            ! iPole - 1 or iPole + 1
!            diff = Eta / ( kT * Pi ) 
!            
!            if ( abs(diff - real(iPole-1,dp)) < &
!                 abs(diff - real(iPole+1,dp)) ) then
!               diff = Eta - kT * Pi * real(iPole-1,dp)
!            else
!               diff = Eta - kT * Pi * real(iPole+1,dp)
!            end if
!         end if
!         diff = abs(diff) / Kelvin
!      else
!         ! The eta value lies close to a pole
!         ! check how close...
!         diff = abs( Eta - kT * Pi * real(iPole,dp) ) / Kelvin
!      end if
!
!    end subroutine calc_Eta_Pole_Diff_Kelvin
!
!  end subroutine ts_print_contour_warnings


!  ! Routine for "pretty" printing the contour points out in the out file
!  subroutine ts_print_contour(contour)
!    use parallel, only : IONode
!    use units,    only : eV
!
!    type(ts_ccontour), pointer :: contour(:)
! **********************
! * LOCAL variables    *
! **********************
!    type(ts_ccontour), pointer :: c
!    character(len=CC_TYPE_LEN) :: ctype
!    character(len=3) :: f_type
!    integer :: i, part, type
!
!    if ( .not. IONode ) return
!    
!    ! Initialize variables
!    part = -1
!    type = -1
!    nullify(c)
!
!    write(*,'(a)') "transiesta: contour integration path:"
!    write(f_type,'(a,i2)') 'a',CC_TYPE_LEN
!    write(*,'(1x,'//trim(f_type)//',''   '',2(tr1,a12),2(tr1,a14))') &
!         "Type  ","Re(c)[eV]","Im(c)[eV]","Re(weight)","Im(weight)"
!    do i = 1 , size(contour)
!       ! loop !
!       c => contour(i)
!       if ( part /= c%part ) then
!          write(*,'(1x,a)') part2str(c)
!          part = c%part
!       end if
!       if ( type /= c%type ) then
!          ctype = type2str(c)
!          type = c%type
!       end if
!       
!       ! Write out the contour information:
!       write(*,'(1x,'//trim(f_type)//','' : '',tr1,2(f12.5,tr1),2(f14.9,tr1))') &
!            ctype,c%c/eV,c%w
!    end do
!    write(*,*) ! New line
!
!  end subroutine ts_print_contour


! Write out the contour to a contour file
!  subroutine ts_io_contour(contour,slabel,suffix)
!    use parallel, only : IONode
!    use units, only : eV
!    type(ts_ccontour), intent(in) :: contour(:)
!    character(len=*), intent(in) :: slabel
!    character(len=*), intent(in), optional :: suffix
!
! *********************
! * LOCAL variables   *
! *********************
!    character(len=len_trim(slabel)+9) :: fname
!    integer :: i, unit, part
!
!    if ( .not. IONode ) return
!
!    if ( present(suffix) ) then
!       fname = trim(trim(slabel)//trim(suffix))
!    else
!       fname = trim(trim(slabel)//'.TSCC')
!    end if
!    call io_assign( unit )
!    open( unit, file=fname, status='unknown' )
!    write(unit,'(a)') '# Complex contour path'
!    write(unit,'(a,a12,3(tr1,a13))') '#','Re(c)[eV]','Im(c)[eV]','Re(w)','Im(w)'
!    part = contour(1)%part
!    do i = 1 , size(contour)
!       if ( part /= contour(i)%part ) then
!          write(unit,*)
!          part = contour(i)%part
!       end if
!       write(unit,'(4(e13.6,tr1))') contour(i)%c/eV,contour(i)%w
!    end do
!    call io_close( unit )
!
!  end subroutine ts_io_contour

end module m_ts_io_contour
