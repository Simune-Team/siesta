! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
#ifndef MPI
! Make sure ELPA is only used for MPI
# undef SIESTA__ELPA
#endif
module m_diag_option

  use precision, only: dp
  
  implicit none
  
  public

  save

  !> Whether diagonalization calls are made via LAPACK (true) or ScaLAPACK (false)
  logical :: Serial = .true.

  !> Whether ScaLAPACK uses a 2D distribution
  logical :: Use2D = .true.
  !> Number of processors on the columns for the 2D distribution
  integer :: ProcessorY = 1
  !> Whether we should use the upper or lower part of the Hermitian/symmetric matrices
  character :: UpperLower = 'L'

  ! Different choices of algorithms.

  ! Choose the algorithm
  ! Note that if any of these are true, the others should be false

  ! The 2stage solvers in LAPACK are (as of 3.7.1) not able to
  ! calculate the eigenvectors.
  ! Hence, they will not be documented.

  !> Use the divide-and-conquer algorithm
  integer, parameter :: DivideConquer = 1
  !> Use the 2-stage divide-and-conquer algorithm (only LAPACK)
  integer, parameter :: DivideConquer_2stage = 2
  !> Use the MRRR algorithm (RRR for LAPACK)
  integer, parameter :: MRRR = 3
  !> Use the 2-stage MRRR algorithm (only LAPACK)
  integer, parameter :: MRRR_2stage = 4
  !> Use the expert driver (if the above are all false)
  integer, parameter :: Expert = 5
  !> Use the 2-stage expert driver (only LAPACK ;\)
  integer, parameter :: Expert_2stage = 6
  !> Use the regular driver
  integer, parameter :: NoExpert = 7
  !> Use the 2-stage regular driver (only LAPACK)
  integer, parameter :: NoExpert_2stage = 8
  !> Use the ELPA driver
  integer, parameter :: ELPA_1stage = 9
  !> Use the 2-stage ELPA driver
  integer, parameter :: ELPA_2stage = 10

  integer :: algorithm = DivideConquer

  !> Tolerance for MRRR (LAPACK) and expert drivers
  real(dp) :: abstol = 1.e-8_dp
  !> Tolerance for expert ScaLAPACK driver
  real(dp) :: orfac = 1.e-3_dp

  !> Memory factor for the real work arrays
  real(dp) :: mem_factor = 1._dp

  logical :: ParallelOverK = .false.

contains

  subroutine read_diag(Gamma, nspin)

    use parallel, only: Nodes
    use fdf, only: fdf_get, leqi

    logical, intent(in) :: Gamma
    integer, intent(in) :: nspin

    logical :: dc_default

    character(len=32) :: algo

    ! No matter what we always re-read them.
    ! In case one will change the algorithm
    ! we will allow that
    ! (However, note that rdiag and cdiag haven't implemented this yet)

#ifdef MPI
    if ( Nodes > 1 .and. .not. Gamma ) then
       ParallelOverK = fdf_get( 'Diag.ParallelOverK', .false. )
       if ( nspin > 2 ) ParallelOverK = .false.
    end if

    if ( Nodes == 1 ) then
       Serial = .true.
       ParallelOverK = .false.
    else if ( ParallelOverK ) then
       Serial = .true.
    else
       Serial = .false.
    end if
#endif

    
    Use2D = fdf_get('Diag.Use2D',.true.)
    ProcessorY = fdf_get('Diag.ProcessorY', max(1, Nodes / 4))
    ProcessorY = max(1, ProcessorY)

    algo = fdf_get('Diag.UpperLower', 'lower')
    if ( leqi(algo, 'lower') .or. leqi(algo, 'l') ) then
       UpperLower = 'L'
    else if ( leqi(algo, 'upper') .or. leqi(algo, 'u') ) then
       UpperLower = 'U'
    else
       call die('diag: Unknown argument to Diag.UpperLower')
    end if

#ifdef SIESTA__MRRR
    dc_default = .false.
#else
    dc_default = .true.
#endif
    
    ! Decide the algorithm
    if ( fdf_get('Diag.DivideAndConquer',dc_default) ) then
       algorithm = DivideConquer
    end if
#ifdef SIESTA__MRRR
    if ( fdf_get('Diag.MRRR',.not. dc_default) ) then
       algorithm = MRRR
    end if
#endif
#ifdef SIESTA__ELPA
    if ( fdf_get('Diag.ELPA',.false.) ) then
       algorithm = ELPA_1stage
    end if
#endif
    if ( fdf_get('Diag.NoExpert',.false.) ) then
       algorithm = NoExpert
    end if

    if ( algorithm == DivideConquer ) then
       algo = fdf_get('Diag.Algorithm', 'D&C')
#ifdef SIESTA__MRRR
    else if ( algorithm == MRRR ) then
       algo = fdf_get('Diag.Algorithm', 'MRRR')
#endif
#ifdef SIESTA__ELPA
    else if ( algorithm == ELPA_1stage ) then
       algo = fdf_get('Diag.Algorithm', 'ELPA')
#endif
    else if ( algorithm == Expert ) then
       algo = fdf_get('Diag.Algorithm', 'expert')
    else
       algo = fdf_get('Diag.Algorithm', 'noexpert')
    end if


    ! Determine the global method
    if ( leqi(algo, 'D&C') .or. leqi(algo, 'divide-and-conquer') .or. &
         leqi(algo, 'DandC') .or. leqi(algo, 'vd') ) then
       algorithm = DivideConquer

    else if ( leqi(algo, 'D&C-2') .or. leqi(algo, 'D&C-2stage') .or. &
         leqi(algo, 'divide-and-conquer-2stage') .or. leqi(algo, 'DandC-2stage') .or. &
         leqi(algo, 'DandC-2') .or. leqi(algo, 'vd_2stage') ) then
#ifdef SIESTA__DIAG_2STAGE
       if ( Serial ) then
          algorithm = DivideConquer_2stage
       else
          algorithm = DivideConquer
       end if
#else
       algorithm = DivideConquer
#endif

#ifdef SIESTA__ELPA
    else if ( leqi(algo, 'elpa') .or. leqi(algo, 'elpa-1stage') ) then
       algorithm = ELPA_1stage
       
       ! ELPA requires 2D and non-serial
       Serial = .false.
       ParallelOverK = .false.
       Use2D = .true.

    else if ( leqi(algo, 'elpa-2stage') .or. leqi(algo, 'elpa-2') ) then
       algorithm = ELPA_2stage

       ! ELPA requires 2D distribution and non-serial
       Serial = .false.
       ParallelOverK = .false.
       Use2D = .true.
#endif
       
#ifdef SIESTA__MRRR
    else if ( leqi(algo, 'MRRR') .or. leqi(algo, 'RRR') .or. &
         leqi(algo, 'vr') ) then
       algorithm = MRRR

    else if ( leqi(algo, 'MRRR-2stage') .or. leqi(algo, 'RRR-2stage') .or. &
         leqi(algo, 'MRRR-2') .or. leqi(algo, 'RRR-2') .or. &
         leqi(algo, 'vr_2stage') ) then
# ifdef SIESTA__DIAG_2STAGE
       if ( Serial ) then
          algorithm = MRRR_2stage
       else
          algorithm = MRRR
       end if
# else
       algorithm = MRRR
# endif
#endif

    else if ( leqi(algo, 'expert') .or. leqi(algo, 'vx') ) then
       algorithm = Expert
       
    else if ( leqi(algo, 'expert-2stage') .or. leqi(algo, 'expert-2') .or. &
         leqi(algo, 'vx_2stage') ) then
#ifdef SIESTA__DIAG_2STAGE
       if ( Serial ) then
          algorithm = Expert_2stage
       else
          algorithm = Expert
       end if
#else
       algorithm = Expert
#endif
       
    else if ( leqi(algo, 'noexpert') .or. leqi(algo, 'regular') .or. &
         leqi(algo, 'v') ) then
       algorithm = NoExpert

    else if ( leqi(algo, 'noexpert-2stage') .or. leqi(algo, 'noexpert-2') .or. &
         leqi(algo, 'regular-2stage') .or. leqi(algo, 'regular-2') .or. &
         leqi(algo, 'v_2stage') ) then
#ifdef SIESTA__DIAG_2STAGE
       if ( Serial ) then
          algorithm = NoExpert_2stage
       else
          algorithm = NoExpert
       end if
#else
       algorithm = NoExpert
#endif

    else

       call die('diag: Unknown routine requested for the diagonalization')

    end if


    ! Retrieve tolerance
    abstol = fdf_get('Diag.AbsTol', 1.e-16_dp)
    orfac = fdf_get('Diag.OrFac', 1.e-3_dp)

    ! Currently this is not used (it shouldn't be needed)
    mem_factor = fdf_get('Diag.Memory', 1.0_dp)
    mem_factor = max(mem_factor, 1.0_dp)

  end subroutine read_diag

  subroutine print_diag()
    use parallel, only: IONode, Nodes

    if ( .not. IONode ) return

    write(*,*) ! new-line
    
    select case ( algorithm )
    case ( DivideConquer ) 
       write(*,'(a,t53,''= '',a)') 'diag: Algorithm', 'D&C'
    case ( DivideConquer_2stage ) 
       write(*,'(a,t53,''= '',a)') 'diag: Algorithm', 'D&C-2stage'
    case ( MRRR )
       write(*,'(a,t53,''= '',a)') 'diag: Algorithm', 'MRRR'
    case ( MRRR_2stage )
       write(*,'(a,t53,''= '',a)') 'diag: Algorithm', 'MRRR-2stage'
    case ( ELPA_1stage )
       write(*,'(a,t53,''= '',a)') 'diag: Algorithm', 'ELPA'
    case ( ELPA_2stage )
       write(*,'(a,t53,''= '',a)') 'diag: Algorithm', 'ELPA-2stage'
    case ( Expert )
       write(*,'(a,t53,''= '',a)') 'diag: Algorithm', 'expert'
    case ( Expert_2stage )
       write(*,'(a,t53,''= '',a)') 'diag: Algorithm', 'expert-2stage'
    case ( NoExpert )
       write(*,'(a,t53,''= '',a)') 'diag: Algorithm', 'noexpert'
    case ( NoExpert_2stage )
       write(*,'(a,t53,''= '',a)') 'diag: Algorithm', 'noexpert-2stage'
    end select


#ifdef MPI
    write(*,'(a,t53,''= '',tr2,l1)') 'diag: Use parallel 2D distribution', Use2D
    write(*,'(a,t53,''= '',i4,'' x '',i4)') 'diag: Parallel 2D distribution', &
         max(1,Nodes / ProcessorY), ProcessorY
    write(*,'(a,t53,''= '',tr2,l1)') 'diag: Parallel over k', ParallelOverK
#endif

    if ( UpperLower == 'L' ) then
       write(*,'(a,t53,''= '',a)') 'diag: Used triangular part', 'Lower'
    else
       write(*,'(a,t53,''= '',a)') 'diag: Used triangular part', 'Upper'
    end if
    
    write(*,'(a,t53,''= '', e10.3)') 'diag: Absolute tolerance', abstol
    write(*,'(a,t53,''= '', e10.3)') 'diag: Orthogonalization factor', orfac

    write(*,'(a,t53,''= '',f7.4)') 'diag: Memory factor', mem_factor

  end subroutine print_diag
  
end module m_diag_option
