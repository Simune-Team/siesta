!==========================================================================*
!                                                                          *
!  TRANSIESTA MODULE m_ts_global_vars : Declaration of the TS variables    *
!  that are accessed in different parts of the code by using a:            *
!      use m_ts_global_vars                                                *
!  declaration, instead of passing as dummy arguments                      *  
!                                                                          *
!  Written by F.D.Novaes, Apr'10                                           *
!==========================================================================*

module m_ts_global_vars
  
  implicit none
  
  save

  ! Whether transiesta is the solver
  logical :: TSmode = .false.

  ! The current iteration in the SCF
  integer :: TSiscf = 1

  ! Controls the change from diagon to transiesta solver
  logical :: TSinit = .false. , TSrun = .false.

  integer :: ts_istep ! FC step in phonon calculation

contains

  subroutine ts_method_init( start )

    use parallel, only : Node
    
    logical, intent(in) :: start

    if ( start ) then

       ! We will immediately start Transiesta
       TSinit = .false.
       TSrun  = .true.
       
       if ( Node == 0 ) then
          write(*,'(a,/)') 'transiesta: Starting immediately'
          write(*,'(a)') '                     ************************'
          write(*,'(a)') '                     *   TRANSIESTA BEGIN   *'
          write(*,'(a)') '                     ************************'
       end if

    else

       ! Tell transiesta to initialize the Hamiltonian with siesta

       TSinit = .true.
       TSrun  = .false.
       
       if ( Node == 0 ) then
          write(*,'(a)') 'transiesta: Initialization run using siesta'
       end if

    end if
    
  end subroutine ts_method_init

end module m_ts_global_vars
