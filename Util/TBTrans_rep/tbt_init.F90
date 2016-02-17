! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
! Handling all the initializations
! This opens MPI channels, ensures the options are read in etc.
!
subroutine tbt_init()
  use fdf
  use sys, only : die
  use precision, only : dp
  use parallel, only : parallel_init, Nodes, IONode
  use parallel, only : ResetFirstCall, ParallelOverK
  use m_tbt_kpoints, only : setup_tbt_kscell
  use m_tbt_options
  use m_timer, only : timer_report
  use m_ts_contour, only : setup_contour, print_contour
  use alloc, only   : alloc_report
  use files, only   : slabel
  use m_timestamp, only : timestamp
  use m_wallclock, only : wallclock
  use m_spin, only : init_spin

  implicit none

  integer :: level
  real(dp) :: threshold
  real(dp) :: CCEmin
#ifdef MPI
  integer :: MPIerror
#endif

  ! Initialise MPI and set processor number
#ifdef MPI
  call MPI_Init(MPIerror)
#endif
  call parallel_init()
  ResetFirstCall =.true.
  ParallelOverK = .true.
#ifdef MPI
  if (.not. fdf_parallel()) then
     call die('tbt_init: ERROR: FDF module doesn''t have parallel support')
  endif
#endif

! Print version information ...........................................
  if (IOnode) then
#ifdef MPI
     write(6,'(a)') 'PARALLEL version'
#else
     write(6,'(a)') 'SERIAL version'
#endif
#ifdef CDF
     write(6,'(a)') 'NetCDF support'
#endif

#ifdef MPI
     if (Nodes.gt.1) then
        write(6,'(/,a,i4,a)') '* Running on ', Nodes, ' nodes in parallel'
     else
        write(6,'(/,a,i4,a)') '* Running in serial mode with MPI'
     endif
#else
     write(6,'(/,a,i4,a)') '* Running in serial mode'
#endif
     call timestamp('Start of run')
     call wallclock('Start of run')
  endif

! Start timer .....................................................
  call timer('tbtrans', 0)
  call timer('tbtrans', 1)
  
! Initialise read .................................................
  call tbt_reinit( sname )

! Set timer report file and threshold .............................
  threshold = fdf_get('timer_report_threshold', 0._dp)
  call timer_report( file=trim(slabel)//'_tbt.times', &
       threshold=threshold )

! Set allocation report level .........................................
! variables level and threshold imported from module siesta_options
  level = fdf_get('alloc_report_level', 0)
  threshold = fdf_get('alloc_report_threshold', 0._dp)
  call alloc_report( level=level, file=trim(slabel)//'_tbt.alloc', &
       threshold=threshold, printNow=.false. )

! initialize spin in the system
  call init_spin()

! Read tbtrans options
  call read_tbt_options()

! Read in k-point cell
  call setup_tbt_kscell()
  
  CCEmin = 0.0_dp ! The contour "lowest" energy does not make sense in a 
! transport calculation... It is not used.
! Create the contour lines
  call setup_contour(IsVolt,0,VoltL,0.0d0,VoltR, &
       0,0,0,0, &
       Emin, Emax, NPoints, &
       CCEmin, GFEta, kt)

! Print out the contour path
  call print_contour()

! Initialization now complete. Flush stdout.
  if (IOnode) call pxfflush( 6 )

end subroutine tbt_init
