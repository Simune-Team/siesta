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
  
  use precision, only : dp
  
  save

  ! Whether transiesta is the solver
  logical :: TSmode = .false.

  ! The current iteration in the SCF
  integer :: TSiscf = 1

  ! Controls the change from diagon to transiesta solver
  logical :: TSinit = .false. , TSrun = .false.

  integer :: ts_istep ! FC step in phonon calculation

end module m_ts_global_vars
