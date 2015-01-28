
module m_ts_global_vars
  
  implicit none
  
  save

  ! Whether transiesta is the solver
  logical :: TSmode = .false.

  ! Controls the change from diagon to transiesta solver
  logical :: TSinit = .false. , TSrun = .false.

  integer :: ts_istep ! FC step in phonon calculation

contains

  subroutine ts_method_init( start )

    use parallel, only : IONode
    
    logical, intent(in) :: start

    ! If we are not in transiesta mode, return
    if ( .not. TSmode ) return

    if ( start ) then

       ! We will immediately start Transiesta
       TSinit = .false.
       TSrun  = .true.
       
       if ( IONode ) then
          write(*,'(/a,/)')'transiesta: Starting immediately'
          write(*,'(a)')   '                     ************************'
          write(*,'(a)')   '                     *   TRANSIESTA BEGIN   *'
          write(*,'(a,/)') '                     ************************'
       end if

    else

       ! Tell transiesta to initialize the Hamiltonian with siesta

       TSinit = .true.
       TSrun  = .false.
       
       if ( IONode ) then
          write(*,'(/,a,/)') 'transiesta: Initialization run using siesta'
       end if

    end if
    
  end subroutine ts_method_init

end module m_ts_global_vars
