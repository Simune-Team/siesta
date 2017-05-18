! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_diagon_opt
  
  implicit none
  
  public

  save

  logical, private :: called = .false.
  
  logical :: AllInOne = .false.
  logical :: DivideConquer = .true.
  logical :: DoRead = .true.
  logical :: FirstCall = .true.
  logical :: NoExpert = .false.
  logical :: Serial = .true.
  logical :: PreRotate = .false.
  logical :: SaveEigenvectors = .false.
  logical :: Use2D = .true.

  integer :: ictxt, i2d_ctxt
  integer :: maxclustersize = 12
  
#ifdef MPI
  logical :: ParallelOverK = .false.
#endif

contains

  subroutine init_diagon_opt(Gamma, nspin)
    use parallel, only: Nodes
    use fdf, only: fdf_get
    logical, intent(in) :: Gamma
    integer, intent(in) :: nspin

    ! if already read, return immediately
    if ( called ) return

#ifdef MPI
    if ( .not. Gamma ) then
       ParallelOverK = fdf_get( 'Diag.ParallelOverK', .false. )
       if ( nspin > 2 ) ParallelOverK = .false.
    end if
#endif

    called = .true.
    
  end subroutine init_diagon_opt
  
end module m_diagon_opt
