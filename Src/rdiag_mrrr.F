! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      module m_rdiag_mrrr
      implicit none

#ifdef SIESTA__NO_MRRR
  ! Make sure if the user-specifies the NO_MRRR pre-processor
  ! we disable MRRR
# ifdef SIESTA__MRRR
#  undef SIESTA__MRRR
# endif
#endif

#ifdef SIESTA__MRRR
      public :: rdiag_mrrr

      CONTAINS
      subroutine rdiag_mrrr(H,S,n,nm,nml,w,Z,neigvec,ierror)
! ***************************************************************************
! Subroutine to find a subset of eigenvalues and eigenvectors of the
! real general eigenvalue problem  H z = w S z,  with H and S
! real symmetric matrices. It uses the MRRR algorithm.

! Written by A. Garcia (Jan 2014), based on rdiag.

! ************************** INPUT ******************************************
! real*8 H(nml,nm)                 : Symmetric H matrix
! real*8 S(nml,nm)                 : Symmetric S matrix
! integer n                        : Order of the generalized  system
! integer nm                       : Right hand dimension of H and S matrices
! integer nml                      : Left hand dimension of H and S matrices
!                                    which is greater than or equal to nm
! integer neigvec                  : No. of eigenvectors to calculate
! ************************** OUTPUT *****************************************
! real*8 w(nml)                    : Eigenvalues
! real*8 Z(nml,nm)                 : Eigenvectors
! integer ierror                   : Flag indicating success code for routine
!                                  :  0 = success
!                                  :  1 = fatal error

!  Modules
      use precision
      use parallel,   only : Node, Nodes, BlockSize
#ifdef MPI
      use mpi_siesta, only : mpi_comm_world
#endif
      use alloc
      use sys,        only : die

      implicit          none

! Passed variables
      integer, intent(out)                 :: ierror
      integer, intent(in)                  :: n
      integer, intent(in)                  :: neigvec
      integer, intent(in)                  :: nm
      integer, intent(in)                  :: nml
      real(dp), intent(in)                 :: H(nml,nm)
      real(dp), intent(in)                 :: S(nml,nm)

      ! It should be possible to dimension Z to use the
      ! actual number of eigenvectors needed...
      real(dp), intent(out)                :: w(nml)
      real(dp), intent(out)                :: Z(nml,nm)

#ifdef MPI
      integer                 :: desch(9)
      integer                 :: MPIerror, mpi_comm_rows, mpi_comm_cols
      integer                 :: ictxt
      integer                 :: np0
      integer                 :: np_cols_1d
      integer                 :: np_rows_1d
      real(dp)                :: dscale
      logical                 :: BlacsOK

! Additional variables for a 2D proc grid
      integer                 :: np_rows, np_cols
      integer                 :: my_prow, my_pcol
      integer                 :: i2d_ctxt
      integer, dimension(9)   :: desc_h2d

      integer :: neigok, nz, vl, vu
      character(len=256) :: msg
      
! Matrices for 2D
      real(dp), pointer, dimension(:,:) ::  h2d =>null()
      real(dp), pointer, dimension(:,:) ::  s2d =>null()
      real(dp), pointer, dimension(:,:) ::  z2d =>null()
      integer ::  nh_rows, nh_cols

      integer                 :: info, i, n_col, n_row
      integer                 :: lwork, liwork
      real(dp), pointer       :: work(:) => null()
      integer , pointer       :: iwork(:) => null()

      real(dp), parameter     :: zero = 0.0_dp
      real(dp), parameter     :: one = 1.0_dp

      integer, external       :: numroc, indxl2g

!***********************************************************************
! Setup                                                                *
!***********************************************************************
      
!     Initialise error flag
      ierror = 0
      
      if ( n == 1 ) then
         w(:) = zero
         w(1) = H(1,1) / S(1,1)
         Z(:,:) = zero
         Z(1,1) = one / sqrt(S(1,1))
         return
      end if
      
!     Start time count
      call timer('rdiag_mrrr',1)
      
      
! Get Blacs context and initialise Blacs grid for main grid
      
      ! First, declare explicitly a context for the original
      ! distribution, which is block-cyclic over the columns,
      ! since our matrices are dimensioned as H(n,no_l) in diagg,
      ! i.e., each processor contains all the rows but just no_l
      ! columns.
      np_rows_1d = 1
      np_cols_1d = Nodes
      call blacs_get( -1, 0, ictxt )
      ! The design of the BLACS is weird.
      ! The above call is equivalent to setting:
      ictxt = mpi_comm_world    ! or whatever comm we are using
      ! The following call will OVERWRITE ictxt with a new context
      ! "C" means "column-major" ordering, which is not surprising,
      ! since we have np_rows_1d = 1
      call blacs_gridinit( ictxt, 'C', np_rows_1d, np_cols_1d )
      
      
      ! Set up blacs descriptors for parallel case
      BlacsOK = .true.
      ! This is wrong in the BSs
      ! since the effective row blocksize is n. It should be
      !                     n, n, n        , BlockSize
      call descinit( desch, n, n, BlockSize, BlockSize, 0, 0, 
     .     ictxt, n, info )
      if (info.ne.0) BlacsOK = .false.
      if (.not.BlacsOK) then
         call die('ERROR : Blacs setup has failed in rdiag!')
      endif
      
      ! Setup secondary grid for 2D distribution
      
      ! Selection of number of processor rows/columns
      ! We try to set up the grid square-like, i.e. start the search for possible
      ! divisors of Nodes with a number next to the square root of Nodes
      ! and decrement it until a divisor is found.
      
      do np_rows = NINT(SQRT(REAL(Nodes))),2,-1
         if(mod(Nodes,np_rows) == 0 ) exit
      end do
      ! at the end of the above loop, Nodes is always divisible by np_rows
      
      np_cols = Nodes/np_rows
      
      !if (Node == 0) print *, "np_rows, col: ", np_rows, np_cols
      
      !call blacs_get(ictxt, 10, i2d_ctxt)
      i2d_ctxt = mpi_comm_world ! or whatever comm we are using
      ! after the above call, i2d_ctxt will contain the communicator
      ! handle we used to initialize ictxt
      ! Now we set up the process grid in row-major, or "natural" order.
      !    0   1   2   3
      !    4   5   6   7
      call blacs_gridinit(i2d_ctxt, 'R', np_rows, np_cols)
      
      ! We need this call to know who we are (my_prow, my_pcol)
      call blacs_gridinfo(i2d_ctxt, np_rows, np_cols,
     .     my_prow, my_pcol)
      !if (Node == 0) print *, "ictxt, i2_ctxt:", ictxt, i2d_ctxt
      
      ! Enquire size of local part of 2D matrices
      ! since we are going to change and in principle we do not know.
      
      ! We might want to change the value of BlockSize...
      
      nh_rows = numroc(n, BlockSize, my_prow, 0, np_rows)
      nh_cols = numroc(n, BlockSize, my_pcol, 0, np_cols)
      
      ! Set up blacs descriptors for 2D case
      ! It is only necessary to give the leading (row) dimension of the array.
      call descinit(desc_h2d, n, n, Blocksize, BlockSize, 0, 0,
     .     i2d_ctxt,  nh_rows, info)
      if (info.ne.0) BlacsOK = .false.
      if (.not.BlacsOK) then
         call die('ERROR : Blacs setup has failed in rdiag!')
      endif
      
      ! Set up workspace arrays for 2D versions of 1D arrays
      call re_alloc( h2d, 1, nh_rows, 1, nh_cols, name='h2d')
      call re_alloc( s2d, 1, nh_rows, 1, nh_cols, name='s2d')
      call re_alloc( z2d, 1, nh_rows, 1, nh_cols, name='z2d')
      
      ! Copy arrays to new 2D distribution
      call pdgemr2d(n, n, H, 1, 1, desch, h2d, 1, 1, desc_h2d, ictxt)
      call pdgemr2d(n, n, S, 1, 1, desch, s2d, 1, 1, desc_h2d, ictxt)
      
      np0 = numroc(n,BlockSize,0,0,np_rows)
      ! lwork is used also in the forward transform.
      ! It will be reset below
      lwork = max(BlockSize*(np0+1),3*BlockSize)
      liwork = 3                ! for now: see query below

!     Allocate workspace arrays
      call re_alloc( work,    1,lwork,   name='work'  )
      call re_alloc( iwork,   1,liwork,  name='iwork'  )
      
!*************************************************************************
! Transform problem to standard eigenvalue problem                             
! by first factorizing the overlap matrix by Cholesky
!*************************************************************************
      call timer('rdiag1',1)
      
      call pdpotrf('U',n,s2d,1,1,desc_h2d,info)
      if (info.ne.0) then
         call die('Error in Cholesky factorisation in rdiag')
      endif
      call timer('rdiag1',2)
      
      !
      ! Now transform
      !
      call timer('rdiag2',1)
      
      call pdsyngst(1,'U',n,h2d,1,1,desc_h2d,s2d,1,1,
     .     desc_h2d,dscale,work,lwork,info)
      if (info.ne.0) then
         call die('Error in forward transformation in rdiag')
      endif
      call timer('rdiag2',2)
      
!     Solve standard eigenvalue problem 
      call timer('rdiag3',1)

      ! Find the optimal work-space sizes
      lwork = -1
      liwork = -1
      call pdsyevr('V','I','U',n,h2d,1,1,desc_h2d,vl,vu,1,
     .     neigvec,neigok,nz,w,z2d,1,1,
     .     desc_h2d,work,lwork,iwork,liwork,info)
      if (info /= 0) then
         write(msg,"(a,i0)") "Error in pdsyevr work estim. Info: ",
     $        info
         call die(msg)
      endif
      lwork = nint(work(1))
      liwork = iwork(1)
      call re_alloc(work, 1, lwork)
      call re_alloc(iwork, 1, liwork)

      call pdsyevr('V','I','U',n,h2d,1,1,desc_h2d,vl,vu,1,
     .     neigvec,neigok,nz,w,z2d,1,1,
     .     desc_h2d,work,lwork,iwork,liwork,info)
      
      if (info /= 0) then
         write(msg,"(a,i0)") "Error in pdsyevr. Info: ", info
         call die(msg)
      endif
      if (node==0) then
         if ((neigok /= neigvec)) then
            write(6,"(a,2i7)") "*** Warning: eigenpair mismatch:"
            write(6,"(a,3i7)")
     $           "Eigenvalues requested and found (w,z): ",
     $           neigvec , neigok, nz
            write(0,"(a,2i7)") "*** Warning: eigenpair mismatch:"
            write(0,"(a,3i7)")
     $           "Eigenvalues requested and found (w,z): ",
     $           neigvec , neigok, nz
         endif
      endif

      call timer('rdiag3',2)
      
!********************************************************************
! Back transformation                                               *
!********************************************************************
      call timer('rdiag4',1)
      
      call pdtrsm('Left','U','N','Non-unit',n,neigvec,one,
     .     s2d,1,1,desc_h2d,z2d,1,1,desc_h2d)
      ! Return to 1-D block-cyclic distribution
      call pdgemr2d(n,n,z2d,1,1,desc_h2d,Z,1,1,desch,ictxt)
      if (dscale.ne.one) then
         call dscal(n,dscale,w,1)
      endif
      
      if (info.ne.0) then
         call die('Error in back transformation in rdiag')
      endif
      call timer('rdiag4',2)
      
      CALL BLACS_GRIDEXIT( ICTXT )
      CALL BLACS_GRIDEXIT( I2D_CTXT )

!     Deallocate workspace arrays
      call de_alloc( work,    name='work'   )
      call de_alloc( iwork,   name='iwork'  )
      
      
      call de_alloc( h2d,   name='h2d'  )
      call de_alloc( s2d,   name='s2d'  )
      call de_alloc( z2d,   name='z2d'  )
      
      call timer('rdiag_mrrr',2)
      
#endif MPI
      end subroutine rdiag_mrrr
#else
      ! Some compilers do not allow empty modules
      CONTAINS
      subroutine rdiag_mrrr_dummy()
      end subroutine
#endif      
      end module m_rdiag_mrrr
