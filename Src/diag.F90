! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.

! Module for performing diagonalization of both real symmetric
! and complex hermitian matrices.
!
! Its current structure is heavily based on the original
!   rdiag and cdiag
! routines and has been optimized and shortened.
!
! The basis of this routine is to provide diagonalization
! and it differs from the original routines in these respects:
!  1. It does NOT implement the generalized algorithms.
!     Every diagonalization is performed in a 2-step manner,
!     a) convert to generalized form (*potrf, *gst),
!     b) solve standard problem
!     This leverages the requirement of the external library
!     to implement the generalized forms (which does
!     exactly the same).
!  2. All pre-rotate and save-eigenvector options has been removed.
!     There is no point in having them if they were not used.
!     Using them should require changes to this.
!  3. In 2D parallel the original implementations *always*
!     allocated the equivalent 2D distribution arrays (H2D, S2D, Z2D).
!     However, when the node's distributed 2D number of elements
!     is <= size(H) then we do not need to allocate room for these
!     nodes. I.e. we can save 3 * size(H2D) elements in memory!
!     This will typically be so, and it may even happen that only a few
!     of the nodes actually require to allocate additional memory.
!  4. Change the default symmetric part from 'Upper' to 'Lower'.
!     This is because the pzhengst (in ScaLAPACK) performs better
!     for 'Lower' (see pdyengst/pzhengst documentations).
!  5. Added more solvers
!     - MRRR
!     - 2stage solvers (only in LAPACK, and currently only for jobz='N')
!     - no-expert drivers (regular ones)

module m_diag

  use precision
  use parallel, only : Node, Nodes, IONode

  implicit none

  private

  save

  ! The BLACS context
  integer :: iCTXT = -1
  integer :: iCTXT2D = -1

#ifdef _DIAG_WORK
  logical :: diag_work_r(2) = .true.
  logical :: diag_work_c(2) = .true.
#endif

  ! Initialization and exit routines for the blacs communicators
  ! Note that the real and complex algorithms share the blacs communicators.
  public :: diag_init, diag_exit

  ! Diagonalization routines (real and complex)
  public :: diag_r, diag_c

contains

  subroutine diag_init()

#ifdef MPI
    use m_diag_option, only: Use2D
    use parallel, only: ProcessorY
    use mpi_siesta, only: MPI_Comm_World
#endif
    
    integer :: nr, nc

    call diag_exit()

#ifdef _DIAG_WORK
    diag_work_r = .true.
    diag_work_c = .true.
#endif

#ifdef MPI
    
    ! Create a new context
    iCTXT = MPI_Comm_World
    nr = 1
    nc = Nodes
    call blacs_gridinit( iCTXT, 'C', nr, nc )

    if ( Use2D ) then

       ! Setup secondary grid for the 2D distribution
       nr = ProcessorY
       nc = Nodes/ProcessorY
       call blacs_get(iCTXT, 10, iCTXT2D)
       call blacs_gridinit(iCTXT2D, 'R', nr, nc)

    end if
#endif

  end subroutine diag_init

  subroutine diag_exit()

#ifdef MPI
    if ( iCTXT >= 0 ) then
       call blacs_gridexit(iCTXT)
       iCTXT = -1
    end if

    if ( iCTXT2D >= 0 ) then
       call blacs_gridexit(iCTXT2D)
       iCTXT2D = -1
    end if
#endif

  end subroutine diag_exit

  subroutine diag_correct_input(algo, jobz, range, uplo, trans, n, neigvec)

    use m_diag_option
    
    integer, intent(inout) :: algo
    character, intent(inout) :: jobz, range, uplo, trans
    integer, intent(in) :: n, neigvec

    ! Use lower part
    uplo = UpperLower

    if ( uplo == 'U' ) then
       trans = 'N'
    else
       trans = 'C'
    end if

    ! Set general Lapack/Scalapack parameters
    if ( neigvec > 0 ) then
       jobz = 'V'
       if ( neigvec == n ) then
          range = 'A'
       else
          range = 'I'
       end if
    else
       ! Note that we will query whether the routines allow
       ! retrieving the eigenvalues only
       jobz = 'N'
       range = 'A'
    end if

    ! Correct for special routines
    if ( Serial ) then

       if ( neigvec > 0 ) then
          select case ( algo )
          case ( DivideConquer_2stage )
             algo = DivideConquer
          case ( MRRR_2stage )
             algo = MRRR
          case ( Expert_2stage )
             algo = Expert
          case ( NoExpert_2stage )
             algo = NoExpert
          end select
       end if

#ifdef MPI
    else

       ! This is only in the case where it starts serial
       ! but then goes to parallel
       select case ( algo )
       case ( DivideConquer_2stage )
          algo = DivideConquer
       case ( MRRR_2stage )
          algo = MRRR
       case ( Expert_2stage )
          algo = Expert
       case ( NoExpert_2stage )
          algo = NoExpert
       end select

       if ( algo == DivideConquer ) then
          ! regardless of neigvec
          jobz = 'V'
          range = 'A'
       end if
#endif
    end if

  end subroutine diag_correct_input

  subroutine diag_c( H, S, n, nm, nml, w, Z, neigvec, iscf, ierror, BlockSize)
! ***************************************************************************
! Subroutine to solve all eigenvalues and eigenvectors of the
! complex general eigenvalue problem  H z = w S z,  with H and S
! complex hermitian matrices.
! Written by G.Fabricius and J.Soler, March 1998
! Rewritten by Julian Gale, August 2004
! Rewritten by Nick R. Papior, July 2017
! ************************** INPUT ******************************************
! complex*16 H(nml,nm)             : Hermitian H matrix
! complex*16 S(nml,nm)             : Hermitian S matrix
! integer n                        : Order of the generalized  system
! integer nm                       : Right hand dimension of H and S matrices
! integer nml                      : Left hand dimension of H and S matrices
!                                    which is greater than or equal to nm
! integer neigvec                  : No. of eigenvectors to calculate
! integer iscf                     : SCF cycle
! integer BlockSize                : Effective parallel block size
! ************************** OUTPUT *****************************************
! real*8 w(nml)                    : Eigenvalues
! complex*16 Z(nml,nm)             : Eigenvectors
! integer ierror                   : Flag indicating success code for routine
!                                  :  0 = success
!                                  : -1 = repeat call as memory is increased
!                                  :  1 = fatal error
! ************************* PARALLEL ****************************************
! When running in parallel this routine now uses Scalapack to perform a
! parallel matrix diagonalisation. This requires Scalapack and Blacs to
! be installed first. Although globally a 1-D block cyclic data distribution
! is employed, locally 1 or 2-D distributions are allowed for.
! The blocksize is now explicitly passed to the routine (A. Garcia, 2016)      
! The routine allows access to all the phases of diagonalisation for fuller
! control, and allows for parallel divide and conquer with reduced memory.
! The presence of eigenvalue clusters is checked for and the memory adjusted
! accordingly to try to guarantee convergence.
! Note that the blacs grid is only initialised on the first call to avoid
! exceeding MPI limits for split/group creation.
!
! When running in parallel and using a 2D distribution it is crucial that the
! Z array is also allocated.
! When using a 2D distribution the algorithm tries to reuse memory without
! double allocation. This is allowed if the 2D distribution # of elements
! is less than or equal to the number of elements in the local arrays.
! (NOTE this is typically so).
! 
!
! It allows for the following routines:
!   LAPACK (Serial or ParallelOverK):
!    - zheevd, zheevd_2stage
!    - zheevr, zheevr_2stage
!    - zheevx, zheevx_2stage
!    - zheev, zheev_2stage
!   ScaLAPACK (Parallel):
!    - pzheevd
!    - pzheevr
!    - pzheevx
!    - pzheev
! To enable the *vr* routines this file needs to be compiled with:
!   -DSIESTA__MRRR
! To enable the *v*_2stage routines this file needs to be compiled with:
!   -DSIESTA__DIAG_2STAGE
! Note that the LAPACK implementation of the 2stage solvers does (in 3.7.1) not
! implement the jobz='V' stage.
! ***************************************************************************

    ! Modules
    use m_diag_option ! just all
    
    use alloc
    use sys, only : die

    ! Passed variables
    integer :: ierror
    integer :: iscf
    integer :: n
    integer :: neigvec
    integer :: nm
    integer :: nml
    integer :: BlockSize
    real(dp) :: w(nml)
    complex(dp), target :: H(nml,nm)
    complex(dp), target :: S(nml,nm)
    complex(dp), target :: Z(nml,nm)

    ! Local variables
    type(allocDefaults) :: oldDefaults
    integer :: algo
#ifdef MPI

    ! Scale when transforming to generalized eigenvalue problem
    real(dp) :: scale

    ! Expert and MRRR drivers
    integer :: nz
    integer,  pointer :: iclustr(:) => null()
    real(dp), pointer :: gap(:) => null()

    ! 1D proc grid
    integer, target :: desch(9)

    ! Additional variables for a 2D proc grid
    integer :: np2d(2), my2d(2), mat_2d(2)
    integer, target :: desc_h2d(9)

    ! Pointers of H, S and Z
    complex(dp), pointer :: Hp(:,:) => null()
    complex(dp), pointer :: Sp(:,:) => null()
    complex(dp), pointer :: Zp(:,:) => null()

#endif

    ! Arguments to routines
    character :: jobz, range, uplo, trans

    ! Lapack info
    integer :: info

    ! Ok-eigenvalues
    integer :: neigok

    integer, pointer :: ifail(:) => null()
    integer, pointer :: isuppz(:) => null()
    integer :: il, iu
    real(dp) :: vl, vu
    
    ! Work sizes
    integer :: liwork, lrwork, lwork
#ifdef MPI
    integer, save :: lrwork_add = 0
#endif

    complex(dp), pointer :: work(:) => null()
    real(dp), pointer :: rwork(:) => null()
    integer, pointer :: iwork(:) => null()

    integer, pointer :: desc(:)

    ! Local variables for loops etc.
    integer :: i, nprc, mclustr

    integer, external :: numroc

    ! Start time count
    call timer('cdiag',1)

    ! Only re-initialize if the routine hasn't been setup.
    if ( iCTXT < 0 ) call diag_init()

!*******************************************************************************
! Setup                                                                        *
!*******************************************************************************
      
    ! Initialise error flag
    ierror = 0

    ! Trap n=1 case, which is not handled correctly otherwise (JMS 2011/07/19)
    if ( n == 1 ) then

       w(:) = 0._dp
       w(1) = real(H(1,1), dp) / real(S(1,1), dp)
       Z(:,:) = 0._dp
       Z(1,1) = 1._dp / sqrt( real(S(1,1), dp) )

       call timer('cdiag', 2)

       return

    end if

    ! Get old allocation defaults and set new ones
    call alloc_default( old=oldDefaults, &
         copy=.false., shrink=.true., &
         imin=1, routine='cdiag' )

    ! vl/il and vu/iu are not currently used, but they must be initialized
    vl = -huge(0._dp)
    vu = huge(0._dp)
    il = 1
    iu = neigvec

#ifdef MPI
    if ( .not. Serial) then

       ! Set up blacs descriptors for parallel case
       call descinit( desch, n, n, BlockSize, BlockSize, 0, 0, &
            iCTXT, n, info )
       if ( info /= 0 ) then
          call die('cdiag: Blacs setup has failed!')
       end if

       if ( Use2D ) then

          ! retrieve secondary grid
          call blacs_gridinfo(iCTXT2D, &
               np2d(1), np2d(2), my2d(1), my2d(2))
          
          ! Enquire size of local part of 2D matices
          mat_2d(1) = numroc(n, BlockSize, my2d(1), 0, np2d(1))
          mat_2d(2) = numroc(n, BlockSize, my2d(2), 0, np2d(2))

          ! Set up blacs descriptors for 2D case
          call descinit(desc_h2d, n, n, Blocksize, BlockSize, 0, 0, &
               iCTXT2D, mat_2d(1), info)
          if ( info /= 0 ) then
             call die('cdiag: Blacs setup has failed!')
          end if

          desc => desc_h2d(:)

       else
          
          desc => desch(:)
          
       end if

    end if
#endif


    algo = algorithm
    call diag_correct_input(algo, jobz, range, uplo, trans, n, neigvec)
    

    ! Initialize the variables for the different routines
    if ( Serial ) then

#ifdef SIESTA__MRRR
       if ( algo == MRRR .or. &
            algo == MRRR_2stage ) then
          call re_alloc(isuppz, 1, 2*n, name='isuppz')
       end if
#endif

#ifdef MPI
    else

       if ( algo == Expert .or. &
            algo == Expert_2stage ) then

          ! We will add a max-cluster-size of 12 to the
          ! lrwork array
          lrwork_add = max(lrwork_add, (12-1) * n)

          nprc = max(product(np2d), Nodes)
          call re_alloc(gap, 1, nprc, name='gap')
          call re_alloc(iclustr, 1, 2*nprc, name='iclustr')

       end if

       if ( Use2D ) then

          if ( product(mat_2d) > nml*nm ) then

             call re_alloc(Hp, 1, mat_2d(1), 1, mat_2d(2), name='H2D')
             call re_alloc(Sp, 1, mat_2d(1), 1, mat_2d(2), name='S2D')
             call re_alloc(Zp, 1, mat_2d(1), 1, mat_2d(2), name='Z2D')

          else
             
             Hp => S
             Sp => Z
             Zp => H
             
          end if
          
       else
          
          Hp => H
          Sp => S
          Zp => Z

       end if

#endif
    end if

    if ( algo == Expert .or. &
         algo == Expert_2stage ) then
       call re_alloc(ifail, 1, n, name='ifail')
    end if

    ! Perform work-size query
    lwork = 1
    lrwork = 1
    liwork = 1
    call re_alloc(work, 1, lwork, name='work')
    call re_alloc(rwork, 1, lrwork, name='rwork')
    call re_alloc(iwork, 1, liwork, name='iwork')

    ! Get memory requirements
    call work_query()
    if ( info /= 0 ) then
       write(*, *) 'cdiag: work-query info ', info
       call die('cdiag: work-query error')
    end if

    ! Add lrwork_add
    lrwork = lrwork + lrwork_add
    
#ifdef _DIAG_WORK
    if ( jobz == 'N' ) then
       if ( diag_work_c(1) ) then
          write(*,'(3a,i2,a,i5,3(a,i12))') &
               'cdiag-debug: jobz=',jobz,', algo=', algo, ', Node=',Node, &
               ', work=', lwork, ', rwork=', lrwork, ', iwork=', liwork
          diag_work_c(1) = .false.
       end if
    else
       if ( diag_work_c(2) ) then
          write(*,'(3a,i2,a,i5,3(a,i12))') &
               'cdiag-debug: jobz=',jobz,', algo=', algo, ', Node=',Node, &
               ', work=', lwork, ', rwork=', lrwork, ', iwork=', liwork
          diag_work_c(2) = .false.
       end if
    end if
#endif

    call re_alloc(work, 1, lwork, name='work')
    call re_alloc(rwork, 1, lrwork, name='rwork')
    call re_alloc(iwork, 1, liwork, name='iwork')

    ! Begin calculation
    neigok = n ! default to all (will only be overwritten for expert drivers)


#ifdef MPI
    if ( Use2D .and. .not. Serial ) then
       ! (re)-Distribute to new 2D layout
       ! Note that it has to be in this order
       call pzgemr2d(n, n, S, 1, 1, desch, Sp, 1, 1, desc_h2d, iCTXT)
       call pzgemr2d(n, n, H, 1, 1, desch, Hp, 1, 1, desc_h2d, iCTXT)
    end if
#endif
    
!*******************************************************************************
! Factorise overlap matrix                                                     *
!*******************************************************************************
    call timer('cdiag1',1)
    if ( Serial ) then
       call zpotrf(uplo,n,S,n,info)
#ifdef MPI
    else
       call pzpotrf(uplo,n,Sp,1,1,desc,info)
#endif
    end if
    if ( info /= 0 ) then
       print *, info
       call die('cdiag: Error in Cholesky factorisation')
    end if
    call timer('cdiag1',2)

!*******************************************************************************
! Transform problem to standard eigenvalue problem                             *
!*******************************************************************************
    call timer('cdiag2',1)
    if ( Serial ) then
       call zhegst(1,uplo,n,H,n,S,n,info)
#ifdef MPI
    else
       call pzhengst(1,uplo,n,Hp,1,1,desc,Sp,1,1, &
            desc,scale,work,lwork,info)
#endif
    end if
    if ( info /= 0 ) then
       print *, info
       call die('cdiag: Error in forward transformation')
    end if
    call timer('cdiag2',2)

!*******************************************************************************
! Solve standard eigenvalue problem                                            *
!*******************************************************************************
    call timer('cdiag3',1)
    if ( Serial ) then

       select case ( algo )
       case ( DivideConquer ) 
          call zheevd(jobz,uplo,n,H,n,&
               w,work,lwork,rwork,lrwork,iwork,liwork, &
               info)
          if ( neigvec > 0 ) then
             call zcopy(n*neigvec,H,1,Z,1)
          end if
          
#ifdef SIESTA__MRRR
       case ( MRRR )
          call zheevr(jobz,range,uplo, &
               n,H,n,vl,vu,il,iu,abstol, &
               neigok,w,Z,n, &
               isuppz, &
               work,lwork,rwork,lrwork,iwork,liwork, &
               info)
#endif
          
       case ( Expert )
          call zheevx(jobz,range,uplo,n,H,n,vl,vu,il,iu,abstol, &
               neigok,w,Z,n, &
               work,lwork,rwork,iwork,ifail, &
               info)
          
       case ( NoExpert )
          call zheev(jobz,uplo,n,H,n,w, &
               work,lwork,rwork, &
               info)
          if ( neigvec > 0 ) then
             call zcopy(n*neigvec,H,1,Z,1)
          end if
          
#ifdef SIESTA__DIAG_2STAGE
       case ( DivideConquer_2stage ) 
          call zheevd_2stage(jobz,uplo,n,H,n,&
               w,work,lwork,rwork,lrwork,iwork,liwork, &
               info)
          if ( neigvec > 0 ) then
             call zcopy(n*neigvec,H,1,Z,1)
          end if
          
# ifdef SIESTA__MRRR
       case ( MRRR_2stage )
          call zheevr_2stage(jobz,range,uplo, &
               n,H,n,vl,vu,il,iu,abstol, &
               neigok,w,Z,n, &
               isuppz, &
               work,lwork,rwork,lrwork,iwork,liwork, &
               info)
# endif

       case ( Expert_2stage )
          call zheevx_2stage(jobz,range,uplo,n,H,n,vl,vu,il,iu,abstol, &
               neigok,w,Z,n, &
               work,lwork,rwork,iwork,ifail, &
               info)

       case ( NoExpert_2stage )
          call zheev_2stage(jobz,uplo,n,H,n,w, &
               work,lwork,rwork, &
               info)
          if ( neigvec > 0 ) then
             call zcopy(n*neigvec,H,1,Z,1)
          end if
#endif

       end select
       
#ifdef MPI
    else

       select case ( algo )
       case ( DivideConquer )
          call pzheevd(jobz,uplo,n,Hp,1,1,desc, &
               w,Zp,1,1,desc, &
               work,lwork,rwork,lrwork,iwork,liwork, &
               info)

#ifdef SIESTA__MRRR
       case ( MRRR )
          call pzheevr(jobz,range,uplo,n,Hp,1,1,desc, &
               vl,vu,il,iu,neigok,nz,w, &
               Zp,1,1,desc, &
               work,lwork,rwork,lrwork,iwork,liwork, &
               info )
#endif

       case ( Expert ) 
          call pzheevx(jobz,range,uplo,n,Hp,1,1,desc,vl,vu,il,iu, &
               abstol,neigok,nz,w,orfac,Zp,1,1,desc, &
               work,lwork,rwork,lrwork,iwork,liwork, &
               ifail,iclustr,gap, &
               info)

          mclustr = 0
          do i = 1, nprc
             mclustr = max(iclustr(2*i) - iclustr(2*i-1), mclustr)
          end do

          if ( info == -25 ) then
             ! LRWORK is too small to compute all the eigenvectors
             ! However, I do not know by how much... ???
             call die('cdiag: Requires bigger rwork')

          else if ( mod(info,2) /= 0 .or. mod(info/8,2) /= 0 ) then
             
             ! One or more than one eigenvector failed to
             ! converge, we should warn the user to decrease
             ! the tolerance.
             if ( IONode ) then
                write(*,*) "cdiag: Decrease the absolute tolerance "//&
                     "due to insufficient eigenvector convergence..."
             end if
             call die('cdiag: Decrease the absolute tolerance!')
             
          else if ( mod(info/2, 2) /= 0 ) then

             ! We need to signal an increase in workspace
             if ( IONode ) then
                write(*,*) "cdiag: Increasing memory and trying diagonalization again"
             end if
             ierror = -1

             i = (mclustr-1) * n
             if ( lrwork_add < i ) then
                
                ! Try to increase the work-size
                lrwork_add = i
                call clean_memory()
                
                call timer('cdiag3', 2)
                call timer('cdiag', 2)

                return
             end if
             
          end if

       case ( NoExpert )
          call pzheev(jobz,uplo,n,Hp,1,1,desc, &
               w,Zp,1,1,desc, &
               work,lwork,rwork,lrwork, &
               info)

       end select

#endif
    end if

    ! Check error flag
    if ( info /= 0 ) then
       ierror = 1
       if ( info < 0 ) then
          call die('cdiag: Illegal argument to standard eigensolver')
       else
          call die('cdiag: Failure to converge standard eigenproblem')
       end if
       if ( neigok < neigvec ) then
          call die('cdiag: Insufficient eigenvalues converged')
       end if
    end if
    ! Ensure that the eigenvalues that haven't been calculated
    ! are "extreme" and hence not applicable
    if ( neigok < n ) then
       do i = neigok + 1 , n
          w(i) = huge(1._dp)
       end do
    end if
    call timer('cdiag3',2)

    
!*******************************************************************************
! Back transformation                                                          *
!*******************************************************************************
    if ( neigvec > 0 ) then
       call timer('cdiag4',1)
       if ( Serial ) then
          call ztrsm('L',uplo,trans,'N', &
               n, neigvec, dcmplx(1._dp, 0._dp), S, n, Z, n)
#ifdef MPI
       else
          call pztrsm('L',uplo,trans,'N',n,neigvec, &
               dcmplx(1._dp, 0._dp),Sp,1,1,desc,Zp,1,1,desc)
          if ( Use2D ) then
             call pzgemr2d(n,n,Zp,1,1,desc,Z,1,1,desc,iCTXT)
          end if
#endif
       end if
       call timer('cdiag4',2)
    end if
#ifdef MPI
    ! Rescale the eigenvalues
    if ( scale /= 1.0_dp .and. .not. Serial ) then
       call dscal(neigok,scale,w,1)
    end if
#endif
    if ( info /= 0 ) then
       call die('cdiag: Error in back transformation')
    end if

    call clean_memory()

    ! Stop time count
    call timer('cdiag',2)
    
  contains

    subroutine clean_memory()
      
!*******************************************************************************
! Clean up                                                                     *
!*******************************************************************************
      
      ! Deallocate workspace arrays
      if ( Serial ) then
#ifdef SIESTA__MRRR
         call de_alloc(isuppz, name='isuppz')
#endif
#ifdef MPI
      else

         if ( algo == Expert .or. &
              algo == Expert_2stage ) then
            call de_alloc(gap, name='gap')
            call de_alloc(iclustr, name='iclustr')
         end if

         if ( Use2D .and. product(mat_2d) > nml*nm ) then
            call de_alloc(Hp, name='H2D')
            call de_alloc(Sp, name='S2D')
            call de_alloc(Zp, name='Z2D')
         end if
#endif
      end if

      if ( algo == Expert .or. &
           algo == Expert_2stage ) then
         call de_alloc(ifail, name='ifail')
      end if

      call de_alloc(work, name='work')
      call de_alloc(iwork, name='iwork')
      call de_alloc(rwork, name='rwork')


      !  Restore old allocation defaults
      call alloc_default( restore=oldDefaults )

    end subroutine clean_memory
    
    subroutine work_query()
      integer :: l_lwork

      ! Initialize
      l_lwork = 0
      
      work(1) = 1
      rwork(1) = 1
      iwork(1) = 1

      lwork = -1
      lrwork = -1
      liwork = -1

      ! Get memory requirements
      if ( Serial ) then
         
         select case ( algo )
         case ( DivideConquer ) 
            call zheevd(jobz,uplo,n,H,n,&
                 w,work,lwork,rwork,lrwork,iwork,liwork, &
                 info)

#ifdef SIESTA__MRRR
         case ( MRRR )
            call zheevr(jobz,range,uplo, &
                 n,H,n,vl,vu,il,iu,abstol, &
                 neigok,w,Z,n, &
                 isuppz, &
                 work,lwork,rwork,lrwork,iwork,liwork, &
                 info)
#endif

         case ( Expert )
            call zheevx(jobz,range,uplo,n,H,n,vl,vu,il,iu,abstol, &
                 neigok,w,Z,n, &
                 work,lwork,rwork,iwork,ifail, &
                 info)
            ! The API does not state that it writes the work-query in
            ! rwork or iwork
            rwork(1) = max(nint(rwork(1)), 7*n)
            iwork(1) = max(iwork(1), 5*n)

         case ( NoExpert )
            call zheev(jobz,uplo,n,H,n,w, &
                 work,lwork,rwork, &
                 info)
            ! The API does not state that it writes the work-query in
            ! rwork
            rwork(1) = max(nint(rwork(1)), 3*n)

#ifdef SIESTA__DIAG_2STAGE
         case ( DivideConquer_2stage ) 
            call zheevd_2stage(jobz,uplo,n,H,n,&
                 w,work,lwork,rwork,lrwork,iwork,liwork, &
                 info)

# ifdef SIESTA__MRRR
         case ( MRRR_2stage )
            call zheevr_2stage(jobz,range,uplo, &
                 n,H,n,vl,vu,il,iu,abstol, &
                 neigok,w,Z,n, &
                 isuppz, &
                 work,lwork,rwork,lrwork,iwork,liwork, &
                 info)
# endif

         case ( Expert_2stage )
            call zheevx_2stage(jobz,range,uplo,n,H,n,vl,vu,il,iu,abstol, &
                 neigok,w,Z,n, &
                 work,lwork,rwork,iwork,ifail, &
                 info)
            ! The API does not state that it writes the work-query in
            ! rwork or iwork
            rwork(1) = max(nint(rwork(1)), 7*n)
            iwork(1) = max(iwork(1), 5*n)

         case ( NoExpert_2stage )
            call zheev_2stage(jobz,uplo,n,H,n,w, &
                 work,lwork,rwork, &
                 info)
            ! The API does not state that it writes the work-query in
            ! rwork
            rwork(1) = max(nint(rwork(1)), 3*n)
#endif

         case default
                        
            call die('cdiag: error in work_query')

         end select

#ifdef MPI
      else

         ! They all require pzhengst
         ! Well pzheevx exists in the generalized form, but
         ! we currently do not use it
         call pzhengst(1, uplo, &
              n, Hp, 1, 1, desc, Sp, 1, 1, desc, &
              scale, &
              work, lwork, &
              info)
         l_lwork = work(1)
         lwork = -1

         select case ( algo )
         case ( DivideConquer )
            call pzheevd(jobz, uplo, n, Hp, 1, 1, desc, &
                 w, &
                 Zp, 1, 1, desc, &
                 work, lwork, rwork, lrwork, iwork, liwork, &
                 info )
            
# ifdef SIESTA__MRRR
         case ( MRRR )
            call pzheevr(jobz, range, uplo, n, Hp, 1, 1, desc, &
                 vl, vu, il, iu, neigok, nz, w, &
                 Zp, 1, 1, desc, &
                 work, lwork, rwork, lrwork, iwork, liwork, &
                 info )
# endif
            
         case ( Expert ) 
            call pzheevx(jobz, range, uplo, n, Hp, 1, 1, desc, &
                 vl, vu, il, iu, abstol, neigok, nz, w, &
                 orfac, &
                 Zp, 1, 1, desc, &
                 work, lwork, rwork, lrwork, iwork, liwork, &
                 ifail, iclustr, gap, &
                 info )

         case ( NoExpert )
            call pzheev(jobz, uplo, n, Hp, 1, 1, desc, &
                 w, Zp, 1, 1, desc, &
                 work, lwork, rwork, lrwork, &
                 info )
            ! Possible bug in scalapack
            ! At least this makes it work!
#ifdef _DIAG_WORK
            if ( jobz == 'N' ) then
               if ( diag_work_c(1) ) then
                  if ( nint(rwork(1)) < 2*n ) then
                     write(*,'(3a,i2,a,i5,a)') &
                          'cdiag-debug: jobz=',jobz,', algo=', algo, ', Node=',Node, &
                          ' BUG in ScaLAPACK work-query'
                  end if
               end if
            else
               if ( diag_work_c(2) ) then
                  if ( nint(rwork(1)) < 4*n-2 ) then
                     write(*,'(3a,i2,a,i5,a)') &
                          'cdiag-debug: jobz=',jobz,', algo=', algo, ', Node=',Node, &
                          ' BUG in ScaLAPACK work-query'
                  end if
               end if
            end if
#endif
            if ( jobz == 'N' ) then
               rwork(1) = max(nint(rwork(1)), 2*n)
            else
               rwork(1) = max(nint(rwork(1)), 4*n-2)
            end if

         case default

            call die('cdiag: error in work_query')

         end select
#endif         
      end if

      lwork = nint(max(nint(real(work(1), dp)), l_lwork) * mem_factor)
      lrwork = nint(rwork(1) * mem_factor)
      liwork = iwork(1)

    end subroutine work_query
    
  end subroutine diag_c

  subroutine diag_r( H, S, n, nm, nml, w, Z, neigvec, iscf, ierror, BlockSize)
! ***************************************************************************
! Subroutine to solve all eigenvalues and eigenvectors of the
! real general eigenvalue problem  H z = w S z,  with H and S
! real symmetry matrices.
! Written by G.Fabricius and J.Soler, March 1998
! Rewritten by Julian Gale, August 2004
! Rewritten by Nick R. Papior, July 2017
! ************************** INPUT ******************************************
! real*8 H(nml,nm)                 : Symmetric H matrix
! real*8 S(nml,nm)                 : Symmetric S matrix
! integer n                        : Order of the generalized  system
! integer nm                       : Right hand dimension of H and S matrices
! integer nml                      : Left hand dimension of H and S matrices
!                                    which is greater than or equal to nm
! integer neigvec                  : No. of eigenvectors to calculate
! integer iscf                     : SCF cycle
! integer BlockSize                : Effective parallel block size
! ************************** OUTPUT *****************************************
! real*8 w(nml)                    : Eigenvalues
! real*8 Z(nml,nm)                 : Eigenvectors
! integer ierror                   : Flag indicating success code for routine
!                                  :  0 = success
!                                  : -1 = repeat call as memory is increased
!                                  :  1 = fatal error
! ************************* PARALLEL ****************************************
! When running in parallel this routine now uses Scalapack to perform a
! parallel matrix diagonalisation. This requires Scalapack and Blacs to
! be installed first. Although globally a 1-D block cyclic data distribution
! is employed, locally 1 or 2-D distributions are allowed for.
! The blocksize is now explicitly passed to the routine (A. Garcia, 2016)      
! The routine allows access to all the phases of diagonalisation for fuller
! control, and allows for parallel divide and conquer with reduced memory.
! The presence of eigenvalue clusters is checked for and the memory adjusted
! accordingly to try to guarantee convergence.
! Note that the blacs grid is only initialised on the first call to avoid
! exceeding MPI limits for split/group creation.
!
! When running in parallel and using a 2D distribution it is crucial that the
! Z array is also allocated.
! When using a 2D distribution the algorithm tries to reuse memory without
! double allocation. This is allowed if the 2D distribution # of elements
! is less than or equal to the number of elements in the local arrays.
! (NOTE this is typically so).
!
! It allows for the following routines:
!   LAPACK (Serial or ParallelOverK):
!    - dsyevd, dsyevd_2stage
!    - dsyevr, dsyevr_2stage
!    - dsyevx, dsyevx_2stage
!    - dsyev, dsyev_2stage
!   ScaLAPACK (Parallel):
!    - pdsyevd
!    - pdsyevr
!    - pdsyevx
!    - pdsyev
! To enable the *vr* routines this file needs to be compiled with:
!   -DSIESTA__MRRR
! To enable the *v*_2stage routines this file needs to be compiled with:
!   -DSIESTA__DIAG_2STAGE
! Note that the LAPACK implementation of the 2stage solvers does (in 3.7.1) not
! implement the jobz='V' stage.
! ***************************************************************************

    ! Modules
    use m_diag_option ! just all
    
    use alloc
    use sys, only : die

    ! Passed variables
    integer :: ierror
    integer :: iscf
    integer :: n
    integer :: neigvec
    integer :: nm
    integer :: nml
    integer :: BlockSize
    real(dp) :: w(nml)
    real(dp), target :: H(nml,nm)
    real(dp), target :: S(nml,nm)
    real(dp), target :: Z(nml,nm)

    ! Local variables
    type(allocDefaults) :: oldDefaults
    integer :: algo
#ifdef MPI

    ! Scale when transforming to generalized eigenvalue problem
    real(dp) :: scale

    ! Expert and MRRR drivers
    integer :: nz
    integer,  pointer :: iclustr(:) => null()
    real(dp), pointer :: gap(:) => null()

    integer, save :: lwork_add = 0

    ! 1D proc grid
    integer, target :: desch(9)

    ! Additional variables for a 2D proc grid
    integer :: np2d(2), my2d(2), mat_2d(2)
    integer, target :: desc_h2d(9)

    ! Pointers of H, S and Z
    real(dp), pointer :: Hp(:,:) => null()
    real(dp), pointer :: Sp(:,:) => null()
    real(dp), pointer :: Zp(:,:) => null()

#endif
    ! Arguments to routines
    character :: jobz, range, uplo, trans

    ! Lapack info
    integer :: info

    ! Ok-eigenvalues
    integer :: neigok

    integer, pointer :: ifail(:) => null()
    integer, pointer :: isuppz(:) => null()
    integer :: il, iu
    real(dp) :: vl, vu
    
    ! Work sizes
    integer :: liwork, lwork

    real(dp), pointer :: work(:) => null()
    integer, pointer :: iwork(:) => null()

    integer, pointer :: desc(:)

    ! Local variables for loops etc.
    integer :: i, nprc, mclustr

    integer, external :: numroc

    ! Start time count
    call timer('rdiag',1)

    ! Only re-initialize if the routine hasn't been setup.
    if ( iCTXT < 0 ) call diag_init()

!*******************************************************************************
! Setup                                                                        *
!*******************************************************************************
      
    ! Initialise error flag
    ierror = 0

    ! Trap n=1 case, which is not handled correctly otherwise (JMS 2011/07/19)
    if ( n == 1 ) then

       w(:) = 0._dp
       w(1) = H(1,1) / S(1,1)
       Z(:,:) = 0._dp
       Z(1,1) = 1._dp / sqrt( S(1,1) )

       call timer('rdiag', 2)

       return

    end if

    ! Get old allocation defaults and set new ones
    call alloc_default( old=oldDefaults, &
         copy=.false., shrink=.true., &
         imin=1, routine='rdiag' )

    ! vl/il and vu/iu are not currently used, but they must be initialized
    vl = -huge(0._dp)
    vu = huge(0._dp)
    il = 1
    iu = neigvec

#ifdef MPI
    if ( .not. Serial) then

       ! Set up blacs descriptors for parallel case
       call descinit( desch, n, n, BlockSize, BlockSize, 0, 0, &
            iCTXT, n, info )
       if ( info /= 0 ) then
          call die('rdiag: Blacs setup has failed!')
       end if

       if ( Use2D ) then

          ! retrieve secondary grid
          call blacs_gridinfo(iCTXT2D, &
               np2d(1), np2d(2), my2d(1), my2d(2))
          
          ! Enquire size of local part of 2D matices
          mat_2d(1) = numroc(n, BlockSize, my2d(1), 0, np2d(1))
          mat_2d(2) = numroc(n, BlockSize, my2d(2), 0, np2d(2))

          ! Set up blacs descriptors for 2D case
          call descinit(desc_h2d, n, n, Blocksize, BlockSize, 0, 0, &
               iCTXT2D, mat_2d(1), info)
          if ( info /= 0 ) then
             call die('rdiag: Blacs setup has failed!')
          end if

          desc => desc_h2d(:)

       else
          
          desc => desch(:)
          
       end if

    end if
#endif

    algo = algorithm
    call diag_correct_input(algo, jobz, range, uplo, trans, n, neigvec)
    

    ! Initialize the variables for the different routines
    if ( Serial ) then

#ifdef SIESTA__MRRR
       if ( algo == MRRR .or. &
            algo == MRRR_2stage ) then
          call re_alloc(isuppz, 1, 2*n, name='isuppz')
       end if
#endif

#ifdef MPI
    else

       if ( algo == Expert .or. &
            algo == Expert_2stage ) then

          ! We will add a max-cluster-size of 12 to the
          ! lwork array
          lwork_add = max(lwork_add, (12-1) * n)

          nprc = max(product(np2d), Nodes)
          call re_alloc(gap, 1, nprc, name='gap')
          call re_alloc(iclustr, 1, 2*nprc, name='iclustr')

       end if

       if ( Use2D ) then

          if ( product(mat_2d) > nml*nm ) then

             call re_alloc(Hp, 1, mat_2d(1), 1, mat_2d(2), name='H2D')
             call re_alloc(Sp, 1, mat_2d(1), 1, mat_2d(2), name='S2D')
             call re_alloc(Zp, 1, mat_2d(1), 1, mat_2d(2), name='Z2D')

          else
             
             Hp => S
             Sp => Z
             Zp => H
             
          end if
          
       else
          
          Hp => H
          Sp => S
          Zp => Z

       end if

#endif
    end if

    if ( algo == Expert .or. &
         algo == Expert_2stage ) then
       call re_alloc(ifail, 1, n, name='ifail')
    end if

    ! Perform work-size query
    lwork = 1
    liwork = 1
    call re_alloc(work, 1, lwork, name='work')
    call re_alloc(iwork, 1, liwork, name='iwork')

    ! Get memory requirements
    call work_query()
    if ( info /= 0 ) then
       write(*, *) 'rdiag: work-query info ', info
       call die('rdiag: work-query error')
    end if

    ! Add lwork_add
    lwork = lwork + lwork_add

#ifdef _DIAG_WORK
    if ( jobz == 'N' ) then
       if ( diag_work_r(1) ) then
          write(*,'(3a,i2,a,i5,2(a,i12))') &
               'rdiag-debug: jobz=',jobz,', algo=', algo, ', Node=',Node, &
               ', work=', lwork, ', iwork=', liwork
          diag_work_r(1) = .false.
       end if
    else
       if ( diag_work_r(2) ) then
          write(*,'(3a,i2,a,i5,2(a,i12))') &
               'rdiag-debug: jobz=',jobz,', algo=', algo, ', Node=',Node, &
               ', work=', lwork, ', iwork=', liwork
          diag_work_r(2) = .false.
       end if
    end if
#endif

    call re_alloc(work, 1, lwork, name='work')
    call re_alloc(iwork, 1, liwork, name='iwork')

    ! Begin calculation
    neigok = n ! default to all (will only be overwritten for expert drivers)


#ifdef MPI
    if ( Use2D .and. .not. Serial ) then
       ! (re)-Distribute to new 2D layout
       ! Note that it has to be in this order
       call pdgemr2d(n, n, S, 1, 1, desch, Sp, 1, 1, desc_h2d, iCTXT)
       call pdgemr2d(n, n, H, 1, 1, desch, Hp, 1, 1, desc_h2d, iCTXT)
    end if
#endif
    
!*******************************************************************************
! Factorise overlap matrix                                                     *
!*******************************************************************************
    call timer('rdiag1',1)
    if ( Serial ) then
       call dpotrf(uplo,n,S,n,info)
#ifdef MPI
    else
       call pdpotrf(uplo,n,Sp,1,1,desc,info)
#endif
    end if
    if ( info /= 0 ) then
       print *, info
       call die('rdiag: Error in Cholesky factorisation')
    end if
    call timer('rdiag1',2)

!*******************************************************************************
! Transform problem to standard eigenvalue problem                             *
!*******************************************************************************
    call timer('rdiag2',1)
    if ( Serial ) then
       call dsygst(1,uplo,n,H,n,S,n,info)
#ifdef MPI
    else
       call pdsyngst(1,uplo,n,Hp,1,1,desc,Sp,1,1, &
            desc,scale,work,lwork,info)
#endif
    end if
    if ( info /= 0 ) then
       print *, info
       call die('rdiag: Error in forward transformation')
    end if
    call timer('rdiag2',2)

!*******************************************************************************
! Solve standard eigenvalue problem                                            *
!*******************************************************************************
    call timer('rdiag3',1)
    if ( Serial ) then

       select case ( algo )
       case ( DivideConquer ) 
          call dsyevd(jobz,uplo,n,H,n,&
               w,work,lwork,iwork,liwork, &
               info)
          if ( neigvec > 0 ) then
             call dcopy(n*neigvec,H,1,Z,1)
          end if
          
#ifdef SIESTA__MRRR
       case ( MRRR )
          call dsyevr(jobz,range,uplo, &
               n,H,n,vl,vu,il,iu,abstol, &
               neigok,w,Z,n, &
               isuppz, &
               work,lwork,iwork,liwork, &
               info)
#endif
          
       case ( Expert )
          call dsyevx(jobz,range,uplo,n,H,n,vl,vu,il,iu,abstol, &
               neigok,w,Z,n, &
               work,lwork,iwork,ifail, &
               info)
          
       case ( NoExpert )
          call dsyev(jobz,uplo,n,H,n,w, &
               work,lwork, &
               info)
          if ( neigvec > 0 ) then
             call dcopy(n*neigvec,H,1,Z,1)
          end if
          
#ifdef SIESTA__DIAG_2STAGE
       case ( DivideConquer_2stage ) 
          call dsyevd_2stage(jobz,uplo,n,H,n,&
               w,work,lwork,iwork,liwork, &
               info)
          if ( neigvec > 0 ) then
             call dcopy(n*neigvec,H,1,Z,1)
          end if
          
# ifdef SIESTA__MRRR
       case ( MRRR_2stage )
          call dsyevr_2stage(jobz,range,uplo, &
               n,H,n,vl,vu,il,iu,abstol, &
               neigok,w,Z,n, &
               isuppz, &
               work,lwork,iwork,liwork, &
               info)
# endif

       case ( Expert_2stage )
          call dsyevx_2stage(jobz,range,uplo,n,H,n,vl,vu,il,iu,abstol, &
               neigok,w,Z,n, &
               work,lwork,iwork,ifail, &
               info)

       case ( NoExpert_2stage )
          call dsyev_2stage(jobz,uplo,n,H,n,w, &
               work,lwork, &
               info)
          if ( neigvec > 0 ) then
             call dcopy(n*neigvec,H,1,Z,1)
          end if
#endif

       end select
       
#ifdef MPI
    else

       select case ( algo )
       case ( DivideConquer )
          call pdsyevd(jobz,uplo,n,Hp,1,1,desc, &
               w,Zp,1,1,desc, &
               work,lwork,iwork,liwork, &
               info)

#ifdef SIESTA__MRRR
       case ( MRRR )
          call pdsyevr(jobz,range,uplo,n,Hp,1,1,desc, &
               vl,vu,il,iu,neigok,nz,w, &
               Zp,1,1,desc, &
               work,lwork,iwork,liwork, &
               info )
#endif

       case ( Expert ) 
          call pdsyevx(jobz,range,uplo,n,Hp,1,1,desc,vl,vu,il,iu, &
               abstol,neigok,nz,w,orfac,Zp,1,1,desc, &
               work,lwork,iwork,liwork, &
               ifail,iclustr,gap, &
               info)

          mclustr = 0
          do i = 1, nprc
             mclustr = max(iclustr(2*i) - iclustr(2*i-1), mclustr)
          end do

          if ( info == -25 ) then
             ! However, I do not know by how much... ???
             call die('rdiag: Requires bigger work')

          else if ( mod(info,2) /= 0 .or. mod(info/8,2) /= 0 ) then
             
             ! One or more than one eigenvector failed to
             ! converge, we should warn the user to decrease
             ! the tolerance.
             if ( IONode ) then
                write(*,*) "rdiag: Decrease the absolute tolerance "//&
                     "due to insufficient eigenvector convergence..."
             end if
             call die('rdiag: Decrease the absolute tolerance!')
             
          else if ( mod(info/2, 2) /= 0 ) then

             ! We need to signal an increase in workspace
             if ( IONode ) then
                write(*,*) "rdiag: Increasing memory and trying diagonalization again"
             end if
             ierror = -1

             i = (mclustr-1) * n
             if ( lwork_add < i ) then
                
                ! Try to increase the work-size
                lwork_add = i
                call clean_memory()
                
                call timer('rdiag3', 2)
                call timer('rdiag', 2)
                
                return
             end if
             
          end if

       case ( NoExpert )
          call pdsyev(jobz,uplo,n,Hp,1,1,desc, &
               w,Zp,1,1,desc, &
               work,lwork, &
               info)

       end select

#endif
    end if

    ! Check error flag
    if ( info /= 0 ) then
       ierror = 1
       if ( info < 0 ) then
          call die('rdiag: Illegal argument to standard eigensolver')
       else
          call die('rdiag: Failure to converge standard eigenproblem')
       end if
       if ( neigok < neigvec ) then
          call die('rdiag: Insufficient eigenvalues converged')
       end if
    end if
    ! Ensure that the eigenvalues that haven't been calculated
    ! are "extreme" and hence not applicable
    if ( neigok < n ) then
       do i = neigok + 1 , n
          w(i) = huge(1._dp)
       end do
    end if
    call timer('rdiag3',2)

    
!*******************************************************************************
! Back transformation                                                          *
!*******************************************************************************
    if ( neigvec > 0 ) then
       call timer('rdiag4',1)
       if ( Serial ) then
          call dtrsm('L',uplo,trans,'N', &
               n, neigvec, 1._dp, S, n, Z, n)
#ifdef MPI
       else
          call pdtrsm('L',uplo,trans,'N',n,neigvec, &
               1._dp,Sp,1,1,desc,Zp,1,1,desc)
          if ( Use2D ) then
             call pdgemr2d(n,n,Zp,1,1,desc,Z,1,1,desc,iCTXT)
          end if
#endif
       end if
       call timer('rdiag4',2)
    end if
#ifdef MPI
    ! Rescale the eigenvalues
    if ( scale /= 1.0_dp .and. .not. Serial ) then
       call dscal(neigok,scale,w,1)
    end if
#endif
    if ( info /= 0 ) then
       call die('rdiag: Error in back transformation')
    end if

    call clean_memory()

    ! Stop time count
    call timer('rdiag',2)
    
  contains

    subroutine clean_memory()
      
!*******************************************************************************
! Clean up                                                                     *
!*******************************************************************************
      
      ! Deallocate workspace arrays
      if ( Serial ) then
#ifdef SIESTA__MRRR
         call de_alloc(isuppz, name='isuppz')
#endif
#ifdef MPI
      else

         if ( algo == Expert .or. &
              algo == Expert_2stage ) then
            call de_alloc(gap, name='gap')
            call de_alloc(iclustr, name='iclustr')
         end if

         if ( Use2D .and. product(mat_2d) > nml*nm ) then
            call de_alloc(Hp, name='H2D')
            call de_alloc(Sp, name='S2D')
            call de_alloc(Zp, name='Z2D')
         end if
#endif
      end if

      if ( algo == Expert .or. &
           algo == Expert_2stage ) then
         call de_alloc(ifail, name='ifail')
      end if

      call de_alloc(work, name='work')
      call de_alloc(iwork, name='iwork')


      !  Restore old allocation defaults
      call alloc_default( restore=oldDefaults )

    end subroutine clean_memory
    
    subroutine work_query()
      integer :: l_lwork

      ! Initialize
      l_lwork = 0
      
      work(1) = 1
      iwork(1) = 1

      lwork = -1
      liwork = -1

      ! Get memory requirements
      if ( Serial ) then
         
         select case ( algo )
         case ( DivideConquer ) 
            call dsyevd(jobz,uplo,n,H,n,&
                 w,work,lwork,iwork,liwork, &
                 info)

#ifdef SIESTA__MRRR
         case ( MRRR )
            call dsyevr(jobz,range,uplo, &
                 n,H,n,vl,vu,il,iu,abstol, &
                 neigok,w,Z,n, &
                 isuppz, &
                 work,lwork,iwork,liwork, &
                 info)
#endif

         case ( Expert )
            call dsyevx(jobz,range,uplo,n,H,n,vl,vu,il,iu,abstol, &
                 neigok,w,Z,n, &
                 work,lwork,iwork,ifail, &
                 info)
            ! The API does not state that it writes the work-query in
            ! iwork
            iwork(1) = max(iwork(1), 5*n)

         case ( NoExpert )
            call dsyev(jobz,uplo,n,H,n,w, &
                 work,lwork, &
                 info)

#ifdef SIESTA__DIAG_2STAGE
         case ( DivideConquer_2stage ) 
            call dsyevd_2stage(jobz,uplo,n,H,n,&
                 w,work,lwork,iwork,liwork, &
                 info)

# ifdef SIESTA__MRRR
         case ( MRRR_2stage )
            call dsyevr_2stage(jobz,range,uplo, &
                 n,H,n,vl,vu,il,iu,abstol, &
                 neigok,w,Z,n, &
                 isuppz, &
                 work,lwork,iwork,liwork, &
                 info)
# endif

         case ( Expert_2stage )
            call dsyevx_2stage(jobz,range,uplo,n,H,n,vl,vu,il,iu,abstol, &
                 neigok,w,Z,n, &
                 work,lwork,iwork,ifail, &
                 info)
            ! The API does not state that it writes the work-query in
            ! iwork
            iwork(1) = max(iwork(1), 5*n)

         case ( NoExpert_2stage )
            call dsyev_2stage(jobz,uplo,n,H,n,w, &
                 work,lwork, &
                 info)
#endif

         case default
                        
            call die('rdiag: error in work_query')

         end select

#ifdef MPI
      else

         ! They all require pzhengst
         ! Well pzheevx exists in the generalized form, but
         ! we currently do not use it
         call pdsyngst(1, uplo, &
              n, Hp, 1, 1, desc, Sp, 1, 1, desc, &
              scale, &
              work, lwork, &
              info)
         l_lwork = work(1)
         lwork = -1

         select case ( algo )
         case ( DivideConquer )
            call pdsyevd(jobz, uplo, n, Hp, 1, 1, desc, &
                 w, &
                 Zp, 1, 1, desc, &
                 work, lwork, iwork, liwork, &
                 info )
            
# ifdef SIESTA__MRRR
         case ( MRRR )
            call pdsyevr(jobz, range, uplo, n, Hp, 1, 1, desc, &
                 vl, vu, il, iu, neigok, nz, w, &
                 Zp, 1, 1, desc, &
                 work, lwork, iwork, liwork, &
                 info )
# endif
            
         case ( Expert ) 
            call pdsyevx(jobz, range, uplo, n, Hp, 1, 1, desc, &
                 vl, vu, il, iu, abstol, neigok, nz, w, &
                 orfac, &
                 Zp, 1, 1, desc, &
                 work, lwork, iwork, liwork, &
                 ifail, iclustr, gap, &
                 info )

         case ( NoExpert )
            call pdsyev(jobz, uplo, n, Hp, 1, 1, desc, &
                 w, Zp, 1, 1, desc, &
                 work, lwork, &
                 info )

         case default

            call die('rdiag: error in work_query')

         end select
#endif         
      end if

      lwork = nint(max(nint(real(work(1), dp)), l_lwork) * mem_factor)
      liwork = iwork(1)

    end subroutine work_query
    
  end subroutine diag_r
  
end module m_diag


subroutine cdiag( H, S, n, nm, nml, w, Z, neigvec, iscf, ierror, BlockSize)
  use precision, only: dp
  use m_diag, only: diag_c
  integer :: nml, nm
  complex(dp), target :: H(nml,nm)
  complex(dp), target :: S(nml,nm)
  real(dp) :: w(nml)
  complex(dp), target :: Z(nml,nm)
  integer :: n
  integer :: neigvec
  integer :: iscf
  integer :: ierror
  integer :: BlockSize
  
  call diag_c(H, S, N, nm, nml, w, Z, neigvec, iscf, ierror, BlockSize)
  
end subroutine cdiag

subroutine rdiag( H, S, n, nm, nml, w, Z, neigvec, iscf, ierror, BlockSize)
  use precision, only: dp
  use m_diag, only: diag_r
  integer :: nml, nm
  real(dp), target :: H(nml,nm)
  real(dp), target :: S(nml,nm)
  real(dp) :: w(nml)
  real(dp), target :: Z(nml,nm)
  integer :: n
  integer :: neigvec
  integer :: iscf
  integer :: ierror
  integer :: BlockSize
  
  call diag_r(H, S, N, nm, nml, w, Z, neigvec, iscf, ierror, BlockSize)
  
end subroutine rdiag
