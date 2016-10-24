! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

      module m_rdiag_elpa

      implicit none

      public :: rdiag_elpa

      CONTAINS
      subroutine rdiag_elpa(H,S,n,nm,nml,w,Z,neigvec,ierror)
! ***************************************************************************
! Subroutine  to solve all eigenvalues and eigenvectors of the
! real general eigenvalue problem  H z = w S z,  with H and S
! real symmetric matrices.

! This version uses the ELPA library, using matrices in 2D block-cyclic
! distribution, and performing the Cholesky and Back-transform steps.


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
      use elpa1, only: get_elpa_row_col_comms
      use elpa2, only: solve_evp_real_2stage
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
      
! Matrices for 2D
      real(dp), pointer, dimension(:,:) ::  h2d =>null()
      real(dp), pointer, dimension(:,:) ::  s2d =>null()
      real(dp), pointer, dimension(:,:) ::  z2d =>null()
      integer ::  nh_rows, nh_cols

      integer                 :: info, i, n_col, n_row
      integer                 :: lwork
      real(dp), pointer       :: work(:) => null()

      real(dp), parameter     :: zero = 0.0_dp
      real(dp), parameter     :: one = 1.0_dp

      integer, external       :: numroc, indxl2g

!*******************************************************************************
! Setup                                                                        *
!*******************************************************************************

! Initialise error flag
      ierror = 0


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
          ictxt = mpi_comm_world   ! or whatever comm we are using
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
     .                 ictxt, n, info )
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
        enddo
        ! at the end of the above loop, Nodes is always divisible by np_rows

          np_cols = Nodes/np_rows

          if (Node == 0) print *, "np_rows, col: ", np_rows, np_cols

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
     .       my_prow, my_pcol)
          if (Node == 0) print *, "ictxt, i2_ctxt:", ictxt, i2d_ctxt

          ! Enquire size of local part of 2D matrices
          ! since we are going to change and in principle we do not know.

          ! We might want to change the value of BlockSize...

          nh_rows = numroc(n, BlockSize, my_prow, 0, np_rows)
          nh_cols = numroc(n, BlockSize, my_pcol, 0, np_cols)

          ! Set up blacs descriptors for 2D case
          ! It is only necessary to give the leading (row) dimension of the array.
          call descinit(desc_h2d, n, n, Blocksize, BlockSize, 0, 0,
     .                  i2d_ctxt,  nh_rows, info)
          if (info.ne.0) BlacsOK = .false.
          if (.not.BlacsOK) then
            call die('ERROR : Blacs setup has failed in rdiag!')
          endif

! Set up workspace arrays for 2D versions of 1D arrays
        call re_alloc( h2d, 1, nh_rows, 1, nh_cols, name='h2d')
        call re_alloc( s2d, 1, nh_rows, 1, nh_cols, name='s2d')
        call re_alloc( z2d, 1, nh_rows, 1, nh_cols, name='z2d')

        call mpi_barrier(mpi_comm_world, mpierror) ! for correct timings only
        if (node==0) print *, "About to re-distribute arrays..."
        call mpi_barrier(mpi_comm_world, mpierror) ! for correct timings only
        
        ! Copy arrays to new 2D distribution
        call pdgemr2d(n, n, H, 1, 1, desch, h2d, 1, 1, desc_h2d, ictxt)
        call pdgemr2d(n, n, S, 1, 1, desch, s2d, 1, 1, desc_h2d, ictxt)

        np0 = numroc(n,BlockSize,0,0,np_rows)
        lwork = max(BlockSize*(np0+1),3*BlockSize)

! Allocate workspace arrays
      call re_alloc( work,    1,lwork,   name='work'  )

        call mpi_barrier(mpi_comm_world, mpierror) ! for correct timings only
        if (node==0) print *, "About to factorize..."
        call mpi_barrier(mpi_comm_world, mpierror) ! for correct timings only
!*******************************************************************************
! Transform problem to standard eigenvalue problem                             
! by first factorizing the overlap matrix by Cholesky
!*******************************************************************************
        call pdpotrf('U',n,s2d,1,1,desc_h2d,info)
        if (info.ne.0) then
          call die('Error in Cholesky factorisation in rdiag')
        endif
        !
        ! Now transform
        !
        call mpi_barrier(mpi_comm_world, mpierror) ! for correct timings only
        if (node==0) print *, "About to transform..."
        call mpi_barrier(mpi_comm_world, mpierror) ! for correct timings only

        call pdsyngst(1,'U',n,h2d,1,1,desc_h2d,s2d,1,1,
     .       desc_h2d,dscale,work,lwork,info)
        if (info.ne.0) then
          call die('Error in forward transformation in rdiag')
        endif
!
!  Now call ELPA solver
!
   ! All ELPA routines need MPI communicators for communicating within
   ! rows or columns of processes, these are set in get_elpa_row_col_comms.

        call mpi_barrier(mpi_comm_world, mpierror) ! for correct timings only
        if (node==0) print *, "About to get comms..."
        call mpi_barrier(mpi_comm_world, mpierror) ! for correct timings only

        call get_elpa_row_col_comms(mpi_comm_world, my_prow, my_pcol,
     $                                  mpi_comm_rows, mpi_comm_cols)

        call mpi_barrier(mpi_comm_world, mpierror) ! for correct timings only
        if (node==0) print *, "About to solve..."
        call mpi_barrier(mpi_comm_world, mpierror) ! for correct timings only

        ! NOTE: Need the full matrix, not just a triangle...
        ! First, use z2d as scratch to hold the transpose of h2d
        call pdtran(n,n,1.d0,h2d,1,1,desc_h2d,0.d0,z2d,1,1,desc_h2d)
        ! Now z2d has in its lower triangle the symmetric part of h2d
        ! So we fill the lower triangle of h2d with the elements of z2d

        do i=1,nh_cols
          ! Get global column corresponding to i and number of local rows up to            
          ! and including the diagonal, these are unchanged in h2d

           n_col = indxl2g(i,     BlockSize, my_pcol, 0, np_cols)
           n_row = numroc (n_col, BlockSize, my_prow, 0, np_rows)
           h2d(n_row+1:nh_rows,i) = z2d(n_row+1:nh_rows,i)
        enddo

        !
        call solve_evp_real_2stage(n, neigvec, h2d, nh_rows, w,
     $              z2d, nh_rows, BlockSize,
     $              mpi_comm_rows, mpi_comm_cols, mpi_comm_world)

! This was the equivalent scalapack call...
!              call pdsyevx(jobz,range,'U',n,h2d,1,1,desc_h2d,vl,vu,1,
!     .                     neigvec,abstol,neigok,nz,w,orfac,z2d,1,1,
!     .                     desc_h2d,work,lwork,iwork,liwork,ifail,
!     .                     iclustr,gap,info)

!*******************************************************************************
! Back transformation                                                          *
!*******************************************************************************

        call mpi_barrier(mpi_comm_world, mpierror) ! for correct timings only
        if (node==0) print *, "About to back-transform..."
        call mpi_barrier(mpi_comm_world, mpierror) ! for correct timings only

              call pdtrsm('Left','U','N','Non-unit',n,neigvec,one,
     .                     s2d,1,1,desc_h2d,z2d,1,1,desc_h2d)

        call mpi_barrier(mpi_comm_world, mpierror) ! for correct timings only
        if (node==0) print *, "About to convert z to 1d..."
        call mpi_barrier(mpi_comm_world, mpierror) ! for correct timings only
              call pdgemr2d(n,n,z2d,1,1,desc_h2d,Z,1,1,desch,ictxt)
              if (dscale.ne.one) then
                 call dscal(n,dscale,w,1)
              endif
              
              if (info.ne.0) then
                 call die('Error in back transformation in rdiag')
              endif

        call mpi_barrier(mpi_comm_world, mpierror) ! for correct timings only
        if (node==0) print *, "About to exit"
        call mpi_barrier(mpi_comm_world, mpierror) ! for correct timings only


        CALL BLACS_GRIDEXIT( ICTXT )
        CALL BLACS_GRIDEXIT( I2D_CTXT )
              ! Do not call blacs_exit, as it seems to be 
              ! shutting down mpi...
              !!!!! CALL BLACS_EXIT( 0 )

! Deallocate workspace arrays
      call de_alloc( work,    name='work'   )
#endif MPI
      end subroutine rdiag_elpa
      end module m_rdiag_elpa

