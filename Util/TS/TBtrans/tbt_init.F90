! Handling all the initializations
! This opens MPI channels, ensures the options are read in etc.
!
subroutine tbt_init()

  use fdf
  use sys, only : die
  use precision, only : dp
  use parallel, only : parallel_init, Node, Nodes, IONode
  use parallel, only : ResetFirstCall, ParallelOverK
  use m_timer, only : timer_report
  use alloc, only   : alloc_report
  use files, only   : slabel
  use m_timestamp, only : timestamp
  use m_wallclock, only : wallclock
  use m_spin
#ifdef NCDF_4
  use nf_ncdf, only : ncdf_IOnode
#endif
#ifdef MPI
    use mpi_siesta, only : MPI_Barrier, MPI_Comm_World
#endif

  use class_Sparsity
  use class_dSpData1D
  use class_dSpData2D

  use m_ts_electrode, only : init_Electrode_HS
  use m_ts_electype

  use m_tbt_kpoint
  use m_tbt_regions
  use m_tbt_kregions, only : tbt_init_kregions, tbt_print_kRegions
  use m_tbt_hs
  use m_tbt_options
  use m_tbt_contour
  use m_tbt_gf
  use m_tbt_save
  use m_tbt_proj
#ifdef NCDF_4
  use m_tbt_dH, only : init_dH_options
#endif

  use m_sparsity_handling

#ifdef _OPENMP
  use omp_lib, only : omp_get_num_threads, omp_get_schedule
  use omp_lib, only : OMP_SCHED_STATIC, OMP_SCHED_DYNAMIC
  use omp_lib, only : OMP_SCHED_GUIDED, OMP_SCHED_AUTO
#else
!$ use omp_lib, only : omp_get_num_threads, omp_get_schedule
!$ use omp_lib, only : OMP_SCHED_STATIC, OMP_SCHED_DYNAMIC
!$ use omp_lib, only : OMP_SCHED_GUIDED, OMP_SCHED_AUTO
#endif

  implicit none

  integer :: level
  real(dp) :: threshold
#ifdef MPI
  integer :: MPIerror
#endif

!$ integer :: i, is
  integer :: iEl
  type(Sparsity) :: tmp_sp
  type(dSpData1D) :: tmp_1D
  type(dSpData2D) :: tmp_2D
  character(len=300) :: sname

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
     
     call prversion()

#ifdef MPI
     if (Nodes > 1) then
        write(*,'(/,a,i0,tr1,a)') '* Running on ', Nodes, &
             'nodes in parallel'
     else
        write(*,'(/,a)') '* Running in serial mode with MPI'
     endif
#else
     write(*,'(/,a)') '* Running in serial mode'
#endif
!$OMP parallel
!$OMP master
!$    i = omp_get_num_threads()
!$    write(*,'(a,i0,a)') '* Running ',i,' OpenMP threads.'
!$    write(*,'(a,i0,a)') '* Running ',Nodes*i,' processes.'
#ifdef _OPENMP
!$    write(*,'(a,i0,a)') '* OpenMP version ', _OPENMP
#endif
!$    call omp_get_schedule(i,is)
!$    select case ( i )
!$    case ( OMP_SCHED_STATIC ) 
!$    write(*,'(a,i0)') '* OpenMP runtime schedule STATIC, chunks ',is
!$    case ( OMP_SCHED_DYNAMIC ) 
!$    write(*,'(a,i0)') '* OpenMP runtime schedule DYNAMIC, chunks ',is
!$    case ( OMP_SCHED_GUIDED ) 
!$    write(*,'(a,i0)') '* OpenMP runtime schedule GUIDED, chunks ',is
!$    case ( OMP_SCHED_AUTO ) 
!$    write(*,'(a,i0)') '* OpenMP runtime schedule AUTO, chunks ',is
!$    case default
!$    write(*,'(a,i0)') '* OpenMP runtime schedule UNKNOWN, chunks ',is
!$    end select
!$OMP end master
!$OMP end parallel
     call timestamp('Start of run')
     call wallclock('Start of run')
  endif

  ! Start timer .....................................................
  call timer('tbtrans', 0)
  call timer('tbtrans', 1)
  
  ! Initialise read .................................................
  call tbt_reinit( sname , slabel )

  ! Set timer report file and threshold .............................
  threshold = fdf_get('timer_report_threshold', 0._dp)
  call timer_report( file=trim(slabel)//'.TBT.times', &
       threshold=threshold )

  ! Set allocation report level .........................................
  ! variables level and threshold imported from module siesta_options
  level = fdf_get('alloc_report_level', 0)
  threshold = fdf_get('alloc_report_threshold', 0._dp)
  call alloc_report( level=level, file=trim(slabel)//'.TBT.alloc', &
       threshold=threshold, printNow=.false. )

#ifdef NCDF_4
  ! In case the user wants to utilize the ncdf library
  call ncdf_IOnode(IONode)
#endif

  ! initialize spin in the system
  ! We read in information regarding spin
  call init_spin()

  ! Initialization now complete. Flush stdout.
  if (IOnode) call pxfflush( 6 )

  ! Initialize the HSfiles
  ! This will read in the HSfile and determine whether we should
  ! do interpolation due to bias not matching any TSHS files
  ! passed to the program.
  ! This will also read in the required information about the system
  call tbt_init_HSfile( nspin )

  ! Read in the options
  ! All generic options regarding the electrodes, etc. are read in
  ! here.
  call tbt_options( spin_idx, TSHS%na_u, TSHS%xa, TSHS%lasto )

  ! Setup the k-points
  call setup_kpoint_grid( TSHS%cell, N_Elec, Elecs )

  ! Check the electrode coordinates etc.
  do iEl = 1 , N_Elec
     call check_Elec(Elecs(iEl), TSHS%nspin, TSHS%cell, &
          TSHS%na_u, TSHS%xa, TSHS%lasto, &
          Elecs_xa_EPS)
  end do

  ! We have the contour now, so we can create the GF files
  do iEl = 1 , N_Elec

     if ( IONode ) write(*,*) ! newline

     ! initialize the electrode for Green's function calculation
     call init_Electrode_HS(Elecs(iEl))

     call do_Green(Elecs(iEl), &
          TSHS%cell,nkpnt,kpoint,kweight, &
          Elecs_xa_Eps, .false. )
     
     ! clean-up
     call delete(Elecs(iEl))
     
  end do

  if ( IONode ) write(*,*) ! newline

  if ( stop_after_GS ) then
     if ( IONode ) then
        write(*,'(a)')'tbtrans: Stopping program per user request.'
        write(*,'(a)')'tbtrans: Done creating all GF files.'
     end if
#ifdef MPI
     call MPI_Barrier(MPI_Comm_World,iEl)
#endif

     call tbt_end()

  end if

  ! Initialize tbtrans regions
  ! We pass a copy of the sparsity pattern as the sparsity pattern
  ! returned is the sparsity pattern minus the buffer atoms!
  ! Hence, in order to change the sparsity patterns of the data
  ! we need to retain both!
  tmp_sp = TSHS%sp
  call tbt_init_regions(N_Elec,Elecs,TSHS%cell,TSHS%na_u,TSHS%lasto, &
       TSHS%dit,tmp_sp, &
       product(TSHS%nsc),TSHS%isc_off)
  call tbt_init_kregions(r_aBuf,N_Elec,Elecs,TSHS%cell, &
       TSHS%dit,tmp_sp,TSHS%na_u,TSHS%xa,TSHS%lasto, &
       TSHS%nsc,TSHS%isc_off)

  call tbt_region_options( save_DATA )

  if ( Node == 0 ) then
     write(*,'(/,a,/)')'tbtrans: Reducing sparsity pattern...'
  end if

  ! Change the data 
  call dSpData1D_to_Sp(TSHS%S_1D,tmp_sp,tmp_1D)
  TSHS%S_1D = tmp_1D
  call delete(tmp_1D)
  call dSpData2D_to_Sp(TSHS%H_2D,tmp_sp,tmp_2D)
  TSHS%H_2D = tmp_2D
  call delete(tmp_2D)
  TSHS%sp = tmp_sp
  call delete(tmp_sp)

  call tbt_print_regions(N_Elec, Elecs)

  call tbt_print_kRegions( TSHS%cell )

#ifdef NCDF_4

  call init_dH_options( )

  ! Initialize the projections here.
  call init_proj( TSHS%na_u , TSHS%lasto , r_aDev , r_oDev, save_DATA )
  call init_proj_T( N_Elec, Elecs , save_DATA )

  call proj_print( N_Elec, Elecs )

#endif

  ! Now we have the sparsity patterns in the correct sparsity
  ! information and we have deleted all un-needed information.

end subroutine tbt_init
                       
