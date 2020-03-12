!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2013, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

! Module for converting a sparsity pattern to an
! "optimal" tri-diagonal matrix...

! Routine for converting a TranSiesta sparsity pattern to a tri-diagonal form
! It requires that the sparsity pattern is fully contained in the current
! processor.

! This module has been fully developed by Nick Papior Andersen, 2013
! nickpapior@gmail.com

! Please contact the author before utilization in other routines.

module m_ts_rgn2trimat

  ! Use regions...
  use precision, only : dp, i8b
  use m_region

  use m_ts_electype
  use m_ts_tri_common, only : GFGGF_needed_worksize
  use m_ts_tri_common, only : nnzs_tri_i8b

  ! method of BTD matrix
  use m_ts_method, only: TS_BTD_A_PROPAGATION, TS_BTD_A_COLUMN
  use m_ts_method, only: ts_A_method

  implicit none

  private

  public :: ts_rgn2trimat

  integer, parameter :: VALID = 0
  integer, parameter :: NONVALID_SIZE = 1
  integer, parameter :: NONVALID_ELEMENT_CONTAIN = 2
  integer, parameter :: NONVALID_TS_ELECTRODE = 3
  
contains

  ! IF nblocks == 0 will create new partition
  subroutine ts_rgn2TriMat(N_Elec, Elecs, IsVolt, &
       dit, sp, r, nblocks, blocks, method, last_block, par)

    use class_OrbitalDistribution
    use class_Sparsity
    use parallel, only : IONode, Node, Nodes
    use fdf, only : fdf_get
#ifdef MPI
    use mpi_siesta
#endif
    use alloc, only : re_alloc, de_alloc
    use geom_helper, only: ucorb

    ! electrodes
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    ! Whether an entire column should be calculated
    logical, intent(in) :: IsVolt
    ! the distribution
    type(OrbitalDistribution), intent(inout) :: dit
    ! The sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! The region that we will create a tri-diagonal matrix on.
    type(tRgn), intent(in) :: r
    ! The sizes of the blocks in the tri-diagonal matrix
    integer, intent(inout) :: nblocks
    integer, pointer :: blocks(:)
    ! Which kind of method should be used to create the tri-diagonal
    integer, intent(in) :: method
    ! Whether we should retain the last block to a fixed size.
    integer, intent(in) :: last_block
    ! Whether the search should be performed in parallel or not
    logical, intent(in), optional :: par

    ! Local variables
    integer, pointer :: blocks_guess(:) => null()
    integer, pointer :: mm_col(:,:) => null()
    integer, pointer :: iwork(:,:) => null()

    integer :: ind
    integer :: i, nblock_guess
    ! In case of parallel
    integer :: init_min_b, init_max_b
    integer :: guess_start, guess_step, guess_end
    logical :: lpar
#ifdef MPI
    integer :: MPIerror
#endif

    call timer('TS-rgn2tri',1)

    lpar = .true.
    if ( present(par) ) lpar = par

    if ( r%n == 1 ) then
       
       ! Simple case, return immediately with default block-size
       
       call re_alloc(blocks, 1, 1, 'blocks', 'tsR2TM')
       blocks(1) = 1
       nblocks = 1

       call timer('TS-rgn2tri', 2)
       return
       
    end if

    ! Establish a guess on the partition of the tri-diagonal 
    ! matrix...
    call re_alloc(blocks_guess, 1, r%n, 'blocks_guess', 'tsR2TM')
    call re_alloc(blocks, 1, r%n, 'blocks', 'tsR2TM')
    call re_alloc(iwork, 1, r%n, 1, 2, 'iwork', 'tsR2TM')
    blocks_guess(:) = 0

    ! create array containing max-min for each ts-orbital
    call re_alloc(mm_col, 1, 2, 1, r%n, 'mm_col', 'tsR2TM')

    ! Set the min/max column indices in the pivoted matrix
    call set_minmax_col(sp, r, mm_col)

    nblocks = 2
    blocks(1) = r%n / 2
    blocks(2) = r%n / 2 + mod(r%n,2)
    if ( last_block > 0 ) then
       nblocks = 2
       blocks(2) = last_block
       blocks(1) = r%n - last_block
       ! Initialize the guess for the 
       call guess_TriMat_last(r%n,mm_col,nblock_guess,blocks_guess,last_block)
       if ( valid_tri(r%n,mm_col,nblock_guess, blocks_guess,last_block) == VALID ) then
          nblocks = nblock_guess
          blocks(1:nblocks) = blocks_guess(1:nblocks)
       end if
    end if

    if ( lpar ) guess_step = Nodes

    ! If the blocks are known by the user to not exceed a certain
    ! size, then we can greatly reduce the guessing step
    ! for huge systems: mm_col(2, 1)

    ! Now figure out the min/max connections.
    ! We do this by cheking the first connection-block and the
    ! min/max ranges

    ! Here we find the orbital in the first range
    ! that connects to the fewest amount of parent orbitals.
    ! This means that the block *could* potentially be as small
    ! as the one found.
    init_min_b = r%n
    init_max_b = 0
    do i = 1, mm_col(2, 1)
      ! bandwidth at row i
      init_min_b = min(init_min_b, mm_col(2, i) - mm_col(1, i))
      ! End guess at end
      init_max_b = max(init_max_b, mm_col(2, i))
    end do
    ! Allow 5% of the minimum block size +/- to search
    ! This *is* too much but to be on the safe side...
    i = max(nint(init_min_b * 0.05), 2)
    ! We allow the splitting of blocks to:
    !   [1 , 4]
    ! but this does not allow any
    !   [1 , >4]
    ! splittings.
    guess_start = max(1, init_min_b / 4 - i)
    guess_start = fdf_get('TS.BTD.Guess1.Min',guess_start)
#ifdef TBTRANS
    guess_start = fdf_get('TBT.BTD.Guess1.Min',guess_start)
#endif
    ! Stepping default to be 1% of minimium bandwith
    guess_step = max(1, nint(init_min_b * 0.01))
    guess_step = fdf_get('TS.BTD.Guess1.Step',guess_step)
#ifdef TBTRANS
    guess_step = fdf_get('TBT.BTD.Guess1.Step',guess_step)
#endif
    guess_end = max(1, min( r%n / 4, init_max_b + i))
    guess_end = fdf_get('TS.BTD.Guess1.Max',guess_end)
#ifdef TBTRANS
    guess_end = fdf_get('TBT.BTD.Guess1.Max',guess_end)
#endif

    ! In case the orbitals of this region is much smaller than
    ! max-block, then use the half 'no'
    guess_end = max(guess_end, guess_start + guess_step)
    guess_end = min(guess_end, r%n / 2)
    guess_start = max(1, min(guess_start, guess_end))

    ! Correct starting guess for the node
    if ( lpar ) guess_start = min(guess_start + Node, guess_end)

    ! We loop over all possibilities from the first part having size
    ! 2 up to and including total number of orbitals in the 
    ! In cases of MPI we do it distributed (however, the collection routine
    ! below could be optimized)
    do i = guess_start , guess_end , guess_step

      ! Make new guess...
      call guess_btd_blocks(r%n,mm_col,i,nblock_guess,blocks_guess,last_block)

      ! Quick escape if the first part is zero (signal from guess_btd_blocks)
      if ( blocks_guess(1) == 0 ) cycle

      if ( valid_tri(r%n,mm_col,nblock_guess,blocks_guess,last_block) == VALID ) then

        ! Try and even out the different blocks
        call smoothen_blocks(method,r%n,mm_col,nblock_guess,blocks_guess,last_block,iwork)

        call select_better(method, nblocks, blocks, nblock_guess, blocks_guess)

      end if

    end do

#ifdef MPI
    if ( lpar ) then
       ! Select the most optimal partition scheme...
       ! Only check up-till the largest block that is actually searched
       do i = 0 , Nodes - 1
          if ( i == Node ) then
             call MPI_Bcast(nblocks, 1, MPI_Integer, i, &
                  MPI_Comm_World, MPIerror)
             call MPI_Bcast(blocks(1), nblocks, MPI_Integer, i, &
                  MPI_Comm_World, MPIerror)
          else
             call MPI_Bcast(nblock_guess, 1, MPI_Integer, i, &
                  MPI_Comm_World, MPIerror)
             call MPI_Bcast(blocks_guess(1), nblock_guess, MPI_Integer, i, &
                  MPI_Comm_World, MPIerror)
             ! Only all the other nodes are allowed to check...
             ! TODO this could be made to a communication tree to limit communication
             call select_better(method, &
                  nblocks,blocks, nblock_guess, blocks_guess)
          end if
       end do
    end if
#endif

    call de_alloc(blocks_guess, 'blocks_guess', 'tsR2TM')
    ! Shrink to the found blocks
    call re_alloc(blocks,1, nblocks, 'blocks', 'tsR2TM', copy=.true., shrink=.true.)

    if ( nblocks < 2 ) then
       
       if ( IONode ) then 
          write(*,'(a)') 'Could not determine an optimal tri-diagonalization &
               &partition'
          write(*,'(2a)') 'Running on region: ',trim(r%name)
          if ( nblocks > 0 ) then
             write(*,'(a,i0)') 'Found: ',nblocks
             write(*,'(1000000(tr1,i0))') blocks
          else
             write(*,'(a)') 'None found...'
          end if
       end if
       call re_alloc(blocks, 1, 3, 'blocks', 'tsR2TM')
       call die('Not yet implemented')

    end if

    ! The blocks now have almost the same size and we will check that it
    ! is a valid thing, if not, we will revert to the other method of
    ! creating the tri-diagonal sparsity pattern
    ! We do not expect this to fail. After all we check that we can
    ! even out the partitions.
    ! The most probable thing is that the electrodes are not
    ! contained in the first two blocks.
    i = valid_tri(r%n, mm_col, nblocks, blocks, last_block)
    if ( i /= VALID ) then
       write(*,'(2a)') 'Running on region: ',trim(r%name)
       write(*,'(a,i0)') 'TranSiesta system size: ',r%n
       write(*,'(a,i0)') 'Current number of blocks: ',nblocks
       write(*,'(10000000(tr1,i0))') blocks
       write(*,'(a,i0)') 'Current part size: ',sum(blocks(:))
       select case ( i )
       case ( NONVALID_SIZE )
          write(*,'(a)') 'The size is not valid.'
       case ( NONVALID_ELEMENT_CONTAIN ) 
          write(*,'(a)') 'Some elements are not contained.'
       case ( NONVALID_TS_ELECTRODE ) 
          write(*,'(a)') 'The electrode is not fully encompassed.'
       case default
          write(*,'(a,i0,a)') 'Row ',-i,' not encompassed in the tri-matrix'
       end select
       call die('Contact the developers. (missing implementation). &
            &You appear to have a special form of electrode.')
    end if

    call de_alloc(mm_col, 'mm_col', 'tsR2TM')
    call de_alloc(iwork, 'iwork', 'tsR2TM')

    call timer('TS-rgn2tri',2)

    ! Write data (if user requests it, otherwise return immediately!)
    call write_pvt_btd(sp, r, nblocks, blocks)

  contains 

    subroutine select_better(method, nblock_cur, blocks_cur, nblock_guess, blocks_guess)

      integer, intent(in) :: method
      integer, intent(inout) :: nblock_cur
      integer, intent(in) :: nblock_guess
      integer, intent(inout) :: blocks_cur(max(nblock_cur,nblock_guess))
      integer, intent(in) :: blocks_guess(nblock_guess)
      
      logical :: copy
      integer :: cur_work, guess_work
      integer :: cur_pad, guess_pad

      ! We check whether the number of elements is smaller
      ! or that the number of blocks is greater (however, this should
      ! in principle always go together)
      ! If the method of optimization is memory:
      if ( method == 0 ) then
        
        copy = faster_blocks(nblock_cur,blocks_cur,nblock_guess,blocks_guess)
        
      else if ( method == 1 ) then
        
        ! We optimize for memory, i.e. we check for number of elements
        ! in this regard we also check whether we should allocate
        ! a work-array in case of bias calculations.
        if ( IsVolt .and. ts_A_method == TS_BTD_A_COLUMN ) then
          call GFGGF_needed_worksize(nblock_guess, blocks_guess, &
              N_Elec, Elecs, guess_pad, guess_work)
          call GFGGF_needed_worksize(nblock_cur, blocks_cur, &
              N_Elec, Elecs, cur_pad, cur_work)

          ! Update values
          guess_work = guess_work + guess_pad
          cur_work = cur_work + cur_pad

          ! Calculate difference between old and new
          cur_pad = cur_work - guess_work
        else
          ! No default values
          cur_pad = 0
        end if

        ! Difference between the two BTD matrices
        guess_pad = nnzs_tri_i8b(nblock_cur, blocks_cur) - &
            nnzs_tri_i8b(nblock_guess, blocks_guess)

        ! total difference in number of elements
        ! If this is positive the guessed BTD matrix has fewer
        ! memory elements allocated than the currently selected one
        cur_pad = guess_pad + cur_pad

        copy = cur_pad > 0
        if ( .not. copy ) then
          ! in case the work-size is the same we fall-back to the fastests method
          if ( cur_pad == 0 ) then
            copy = faster_blocks(nblock_cur,blocks_cur,nblock_guess,blocks_guess)
          end if
        end if

      else
        call die('Unknown optimization scheme for the tri-mat')
      end if

      if ( copy ) then
        nblock_cur = nblock_guess
        blocks_cur(1:nblock_cur) = blocks_guess(1:nblock_cur)
      end if

    end subroutine select_better

    subroutine write_pvt_btd(sp, r, nblocks, blocks)
      type(Sparsity), intent(inout) :: sp
      type(tRgn), intent(in) :: r
      integer, intent(in) :: nblocks
      integer, intent(in) :: blocks(nblocks)
      character(len=64) :: fname

      type(tRgn) :: pvt
      integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
      integer :: i, io, j, jo, ind, no_u, iu

      if ( .not. IONode ) return

      fname = fdf_get('TS.BTD.Output',' ')
      if ( len_trim(fname) == 0 ) return

      call attach(sp,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
          nrows_g=no_u)

      call rgn_init(pvt, no_u, val=0)
      ! Create the back-pivoting array
      do i = 1 , r%n
        pvt%r(r%r(i)) = i
      end do

      ! Write out the BTD format in a file to easily be processed
      ! by python, this is the pivoted sparsity pattern
      call io_assign(iu)
      open(iu, file=trim(fname)//'.sp',action='write')
      write(iu,'(i0)') r%n
      do i = 1 , r%n
        io = r%r(i)
        if ( l_ncol(io) == 0 ) cycle
        do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
          jo = UCORB(l_col(ind),no_u)
          j = pvt%r(jo)
          if ( j > i ) cycle ! only print upper half (it is hermitian)
          write(iu,'(2(i0,tr1),i1)') i, j, 1
        end do
      end do
      call io_close(iu)

      call io_assign(iu)
      open(iu, file=trim(fname)//'.pvt',action='write')
      write(iu,'(i0)') r%n
      do i = 1 , r%n
        ! store the pivoting table such that sp[i,j] = orig[ pvt[i] , pvt[j] ]
        write(iu,'(i0,tr1,i0)') i, r%r(i)
      end do
      call io_close(iu)

      call io_assign(iu)
      open(iu, file=trim(fname)//'.btd',action='write')
      ! write the tri-mat blocks
      write(iu,'(i0)') nblocks
      do i = 1 , nblocks
        write(iu,'(i0)') blocks(i)
      end do
      call io_close(iu)

      call rgn_delete(pvt)

    end subroutine write_pvt_btd

  end subroutine ts_rgn2TriMat

  subroutine guess_btd_blocks(no,mm_col,first_block,nblocks,blocks,last_block)

    integer, intent(in) :: no, mm_col(2,no) ! number of orbitals, max,min
    integer, intent(in) :: first_block
    integer, intent(inout) :: nblocks
    integer, intent(inout) :: blocks(no)
    integer, intent(in) :: last_block

    ! Local variables
    integer :: N

    if ( first_block > no ) &
        call die('Not allowed to do 1 tri-diagonal part')

    nblocks = 1
    blocks(1) = first_block
    N = blocks(1)
    
    do while ( N < no )
      
      if ( nblocks > no ) then
        print *,'Error',nblocks,no
        call die('Size error when guessing the tri-mat size')
      end if

      nblocks = nblocks + 1
      call set_next_block_size(no, mm_col, N-blocks(nblocks-1)+1, N, blocks(nblocks))
      
      N = N + blocks(nblocks)
      ! if a last-part was "forced" we do this here...
      if ( N + last_block > no ) then
        ! We need to add the former part with the "too many"
        ! orbitals
        blocks(nblocks) = blocks(nblocks) + no - N
        N = no
      end if

    end do

    if ( last_block > 0 ) then
      
      ! Correct so that we actually do contain the last_block
      ! in the last one
      N = blocks(nblocks) - last_block
      blocks(nblocks) = last_block
      blocks(nblocks-1) = blocks(nblocks-1) + N
      
      ! Signal that this is a faulty TRIMAT
      if ( N < 0 ) blocks(1) = 0
      
    end if

  end subroutine guess_btd_blocks

  subroutine guess_TriMat_last(no,mm_col,nblocks,blocks,last_block)

    integer, intent(in) :: no, mm_col(2,no) ! number of orbitals, max,min
    integer, intent(inout) :: nblocks
    integer, intent(inout) :: blocks(:)
    integer, intent(in) :: last_block

    ! Local variables
    integer :: N

    nblocks = 1
    blocks(1) = last_block
    N = blocks(1)
    do while ( N < no )
       nblocks = nblocks + 1
       if ( nblocks > size(blocks) ) then
          print *,'Error',nblocks,size(blocks)
          call die('Size error when guessing the tri-mat size')
       end if
       call guess_prev_block_size(no, mm_col, nblocks, nblocks, blocks)
       N = N + blocks(nblocks)
    end do

    ! Reverse
    blocks(1:nblocks) = blocks(nblocks:1:-1)

  end subroutine guess_TriMat_last

  function faster_blocks(n1,blocks1,n2,blocks2) result(faster)
    integer, intent(in) :: n1, blocks1(n1)
    integer, intent(in) :: n2, blocks2(n2)
    logical :: faster

    integer :: i, n
    real(dp) :: p1, p2, diff

    ! We estimate the fastest algorithm
    ! by the number of operations the matrices make
    ! inv = O^3, mm = O^3

    p1 = (R(blocks1(1)) + R(blocks1(2))) * R(blocks1(1)) * R(blocks1(1))
    p2 = (R(blocks2(1)) + R(blocks2(2))) * R(blocks2(1)) * R(blocks2(1))

    diff = p1 - p2

    n = min(n1, n2)
    do i = 2, n - 1

      p1 = (R(blocks1(i)) + R(blocks1(i-1)) + R(blocks1(i+1))) * R(blocks1(i)) * R(blocks1(i))
      p2 = (R(blocks2(i)) + R(blocks2(i-1)) + R(blocks2(i+1))) * R(blocks2(i)) * R(blocks2(i))

      diff = diff + p1 - p2

    end do

    ! The last entry
    p1 = (R(blocks1(n1)) + R(blocks1(n1-1))) * R(blocks1(n1)) * R(blocks1(n1))
    p2 = (R(blocks2(n2)) + R(blocks2(n2-1))) * R(blocks2(n2)) * R(blocks2(n2))

    diff = diff + p1 - p2

    if ( n1 > n2 ) then
      faster = .true.

      do i = n2, n1 - 1

        p1 = (R(blocks1(i)) + R(blocks1(i-1)) + R(blocks1(i+1))) * R(blocks1(i)) * R(blocks1(i))

        diff = diff + p1

        if ( diff > 0._dp ) return

      end do

    else
      faster = .false.

      do i = n1, n2 - 1

        p2 = (R(blocks2(i)) + R(blocks2(i-1)) + R(blocks2(i+1))) * R(blocks2(i)) * R(blocks2(i))

        diff = diff - p2

        if ( diff < 0._dp ) return

      end do

    end if

    faster = diff > 0._dp

  contains

    elemental function R(i) result(o)
      integer, intent(in) :: i
      real(dp) :: o
      o = real(i, dp)
    end function R
      
  end function faster_blocks


  ! We will guess the size of this (part) tri-diagonal part by
  ! searching for the size of the matrix that matches that of the previous 
  ! part.
  ! We save it in blocks(part)
  subroutine set_next_block_size(no,mm_col,sRow,eRow,block_size)
    integer, intent(in) :: no, mm_col(2,no)
    ! the part we are going to create
    ! eRow is sum(blocks(:-1))
    integer, intent(in) :: sRow, eRow
    integer, intent(inout) :: block_size
    ! Local variables
    integer :: i, mcol
    
    ! We will check in between the above selected rows and find the 
    ! difference in size...
    mcol = 0
    do i = sRow, eRow
      ! this is the # of elements from the RHS of the 'part-1'
      ! part of the tridiagonal matrix and out to the last element of
      ! this row...
      mcol = max(mcol, mm_col(2, i))
    end do

    ! In case there is no connection, we should
    ! force the next-part to be 1!
    block_size = max(1, mcol - eRow)

  end subroutine set_next_block_size

  ! We will guess the size of this (part) tri-diagonal part by
  ! searching for the size of the matrix that matches that of the previous 
  ! part.
  ! We save it in nblocks(part)
  subroutine guess_prev_block_size(no,mm_col,block,nblocks,blocks)
    integer, intent(in) :: no, mm_col(2,no)
    ! the part we are going to create
    integer, intent(in) :: block, nblocks
    integer, intent(inout) :: blocks(nblocks)
    ! Local variables
    integer :: i, sRow, eRow, mcol
    
    ! We are now checking a future block
    ! Hence we must ensure that the size is what
    ! is up to the last blocks size, and thats it...
    eRow = no
    if ( block > 2 ) then
       do i = 1 , block - 2
          eRow = eRow - blocks(i)
       end do
    end if
    sRow = eRow - blocks(block-1) + 1

    ! We will check in between the above selected rows and find the 
    ! difference in size...
    blocks(block) = 0
    do i = sRow, eRow
       ! this is the # of elements from the RHS of the 'block-1'
       ! block of the tridiagonal matrix and out to the last element of
       ! this row...
       mcol = sRow - mm_col(1,i)
       if ( blocks(block) < mcol ) then
          blocks(block) = mcol
       end if
    end do

    ! In case there is actually no connection, we should
    ! force the next-block to be 1!
    if ( blocks(block) == 0 ) then
       blocks(block) = 1
    end if

  end subroutine guess_prev_block_size

  subroutine smoothen_blocks(method,no,mm_col,nblocks,blocks, last_block, iwork)

    integer, intent(in) :: method ! the method used for creating the blocks
    integer, intent(in) :: no, mm_col(2,no)
    ! the block we are going to create
    integer, intent(in) :: nblocks
    integer, intent(inout) :: blocks(nblocks)
    integer, intent(in) :: last_block
    integer, intent(inout) :: iwork(nblocks,2)
    ! Local variables
    integer :: n, d
    logical :: changed

    if ( nblocks == 1 ) return

    ! Setup the cumultative blocks
    iwork(1,1) = blocks(1)
    iwork(1,2) = blocks(1)
    do n = 2, nblocks
      iwork(n,1) = blocks(n)
      iwork(n,2) = blocks(n) + iwork(n-1,2)
    end do

    select case ( method )
    case ( 0 ) ! performance

      do
        changed = .false.
        do n = 2 , nblocks - 1
          call smoothen_block(no, mm_col, nblocks, blocks, iwork(1,2), n, last_block)
          call diff_perf(n, nblocks, blocks, iwork(1,1), d)
          call select(d)
        end do
        n = 1
        call smoothen_block(no, mm_col, nblocks, blocks, iwork(1,2), n, last_block)
        call diff_perf(n, nblocks, blocks, iwork(1,1), d)
        call select(d)
        n = nblocks
        call smoothen_block(no, mm_col, nblocks, blocks, iwork(1,2), n, last_block)
        call diff_perf(n, nblocks, blocks, iwork(1,1), d)
        call select(d)
        
        if ( .not. changed ) exit
      end do

    case ( 1 ) ! memory
      
      do
        changed = .false.
        do n = 2 , nblocks - 1
          call smoothen_block(no, mm_col, nblocks, blocks, iwork(1,2), n, last_block)
          call diff_mem(n, nblocks, blocks, iwork(1,1), d)
          call select(d)
        end do
        n = 1
        call smoothen_block(no, mm_col, nblocks, blocks, iwork(1,2), n, last_block)
        call diff_mem(n, nblocks, blocks, iwork(1,1), d)
        call select(d)
        n = nblocks
        call smoothen_block(no, mm_col, nblocks, blocks, iwork(1,2), n, last_block)
        call diff_mem(n, nblocks, blocks, iwork(1,1), d)
        call select(d)
        
        if ( .not. changed ) exit
      end do
      
    case default

      call die('not implemented')

    end select

  contains

    subroutine select(d)
      integer, intent(in) :: d
      if ( d == 0 ) then
        ! Simply store. It could be that we swapped a few things
        call store_block(nblocks, blocks, iwork(1,1), n)
        call store_cum_block(nblocks, blocks, iwork(1,2), n)
      else if ( d > 0 ) then
        ! We have that the just tested BTD matrix is not as performant
        ! as the current
        call store_block(nblocks, iwork(1,1), blocks, n)
        ! the cumultative blocks are not changed since we restore them
      else
        ! the penalty is negative which means we have found a new candidate
        call store_block(nblocks, blocks, iwork(1,1), n)
        call store_cum_block(nblocks, blocks, iwork(1,2), n)
        changed = .true.
      end if
    end subroutine select

    subroutine store_block(nblocks, blocks, store_blocks, n)
      integer, intent(in) :: n, nblocks
      integer, intent(in) :: blocks(nblocks)
      integer, intent(inout) :: store_blocks(nblocks)
      
      store_blocks(n) = blocks(n)
      if ( n == 1 ) then
        store_blocks(n+1) = blocks(n+1)
      else if ( n == nblocks ) then
        store_blocks(n-1) = blocks(n-1)
      else
        store_blocks(n-1) = blocks(n-1)
        store_blocks(n+1) = blocks(n+1)
      end if
      
    end subroutine store_block

    subroutine store_cum_block(nblocks, blocks, cum_blocks, n)
      integer, intent(in) :: n, nblocks
      integer, intent(in) :: blocks(nblocks)
      integer, intent(inout) :: cum_blocks(nblocks)

      if ( n == 1 ) then
        cum_blocks(n) = blocks(n)
        cum_blocks(n+1) = cum_blocks(n) + blocks(n+1)
      else if ( n == nblocks ) then
        if ( nblocks == 2 ) then
          cum_blocks(n-1) = blocks(n-1)
        else
          cum_blocks(n-1) = cum_blocks(n-2) + blocks(n-1)
        end if
        cum_blocks(n) = cum_blocks(n-1) + blocks(n)
      else if ( n == 2 ) then
        cum_blocks(n-1) = blocks(n-1)
        cum_blocks(n) = cum_blocks(n-1) + blocks(n)
        cum_blocks(n+1) = cum_blocks(n) + blocks(n+1)
      else
        cum_blocks(n-1) = cum_blocks(n-2) + blocks(n-1)
        cum_blocks(n) = cum_blocks(n-1) + blocks(n)
        cum_blocks(n+1) = cum_blocks(n) + blocks(n+1)
      end if

    end subroutine store_cum_block
    
    function changed_block(nblocks, blocks, store_blocks, n) result(changed)
      integer, intent(in) :: n, nblocks
      integer, intent(in) :: blocks(nblocks)
      integer, intent(inout) :: store_blocks(nblocks)
      logical :: changed
      
      changed = store_blocks(n) /= blocks(n)
      if ( n == 1 ) then
        changed = changed .or. store_blocks(n+1) /= blocks(n+1)
      else if ( n == nblocks ) then
        changed = changed .or. store_blocks(n-1) /= blocks(n-1)
      else
        changed = changed .or. store_blocks(n-1) /= blocks(n-1)
        changed = changed .or. store_blocks(n+1) /= blocks(n+1)
      end if
      
    end function changed_block

  end subroutine smoothen_blocks

  subroutine smoothen_block(no,mm_col,nblocks,blocks, cum_blocks, n, last_block)
    integer, intent(in) :: no, mm_col(2,no)
    ! the block we are going to create
    integer, intent(in) :: nblocks, cum_blocks(nblocks)
    integer, intent(inout) :: blocks(nblocks)
    integer, intent(in) :: n, last_block
    
    ! Local variables
    integer :: copy_block
    integer :: sRow, eRow

    if ( last_block > 0 ) then
      
      ! We need the last one to be of a certain size.
      if ( n >= nblocks - 1 ) return
      
    end if

    select case ( nblocks )
    case ( 1 )
      call die('You cannot use tri-diagonalization &
          &without having at least 2 blocks')
    case ( 2 )
      ! For only two blocks there is no gain/loss regardless of matrix
      ! So set them equal
      blocks(1) = no / 2
      blocks(2) = blocks(1) + mod(no, 2)
      return
    end select

    ! Initialize
    copy_block = 0

    ! We do not allow to diminish the edges
    ! This is because we need additional checks of their
    ! regions, hence, we let the central blocks partition out
    ! the regions
    if ( n == 1 ) then
      
      ! We can always align the two blocks (sign is pointless here)
      sRow = 0
      do while ( blocks(1) /= copy_block )
        copy_block = blocks(1)
        call even_if_larger(sRow,blocks(1),blocks(2),sign=-1)
      end do

      return

    else if ( n == nblocks ) then

      ! We can always align the two blocks (sign is pointless here)
      sRow = 0
      do while ( blocks(nblocks) /= copy_block )
        copy_block = blocks(nblocks)
        call even_if_larger(sRow,blocks(nblocks),blocks(nblocks-1),sign=1)
      end do

      return

    end if

    ! Find min/max columns checked
    sRow = cum_blocks(n-1) + 1
    eRow = cum_blocks(n)

    ! We will continue to shift columns around
    ! until we do not shift columns any more...
    copy_block = 0
    do while ( blocks(n) /= copy_block )
        
       ! Copy the current blockition so that we can check in the
       ! next iteration...
       copy_block = blocks(n)

       ! TODO - consider adding a flag for memory reduced utilization of TRI
       ! this will require the two electrodes blocks to be larger
       ! than the other blocks...

       ! 1. if you wish to shrink it left, then:
       !    the first row must not have any elements
       !    extending into the right block
       if ( mm_col(2,sRow) <= eRow ) then
          call even_if_larger(sRow,blocks(n),blocks(n-1),sign= 1)
       end if

       ! 2. if you wish to shrink it right, then:
       !    the last row must not have any elements
       !    extending into the left block
       if ( sRow <= mm_col(1,eRow) ) then
          call even_if_larger(eRow,blocks(n),blocks(n+1),sign=-1)
       end if

    end do
     
  contains

    pure subroutine even_if_larger(Row,p1,p2,sign)
      integer , intent(inout) :: Row, p1, p2
      integer, intent(in) :: sign
      ! If we don't have p2 + 1 we could end up
      ! in 12 -- 11 swaps
      if ( p1 > p2 + 1) then
         p1 = p1 - 1
         p2 = p2 + 1
         Row = Row + sign
      end if
    end subroutine even_if_larger

  end subroutine smoothen_block

! Min and max column requires that the sparsity pattern
! supplied has already stripped off the buffer orbitals.
! Otherwise this will fail
  subroutine set_minmax_col(sp, r, mm_col)
    use class_Sparsity
    use geom_helper, only : UCORB
    ! The sparsity pattern
    type(Sparsity), intent(inout) :: sp
    type(tRgn), intent(in) :: r
    integer, intent(inout) :: mm_col(2,r%n)
    ! The results
    type(tRgn) :: pvt
    integer :: ir, row, ptr, nr, j
    integer, pointer :: l_col(:), l_ptr(:), ncol(:)

    call attach(sp,n_col=ncol,list_ptr=l_ptr,list_col=l_col,nrows_g=nr)

    ! Using a pivoting table reduces overhead
    ! of performing rgn_pivot on a non-sorted
    ! region! SUBSTANTIALLY!
    call rgn_init(pvt, nr, val=0)

    ! Create the back-pivoting array
    do ir = 1 , r%n
      pvt%r(r%r(ir)) = ir
    end do
    
    do ir = 1 , r%n

       ! Get original sparse matrix row
       row = r%r(ir)
       
       ! initialize to region row
       mm_col(1,ir) = ir
       mm_col(2,ir) = ir

       ! Loop on the sparse entries
       do ptr = l_ptr(row) + 1 , l_ptr(row) + ncol(row)
          j = pvt%r( ucorb(l_col(ptr),nr) )
          if ( j > 0 ) then
             if ( j < mm_col(1,ir) ) mm_col(1,ir) = j
             if ( j > mm_col(2,ir) ) mm_col(2,ir) = j
          end if
       end do

    end do

    call rgn_delete(pvt)

  end subroutine set_minmax_col

  !< Calculate the difference in memory requirement for two BTD
  !!  mem(block1) - mem(block2) = d
  !! If d > 0 block1 uses more memory than block2
  subroutine diff_mem(p, n, block1, block2, mem)
    integer, intent(in) :: p
    integer, intent(in) :: n
    integer, intent(in) :: block1(n), block2(n)
    integer, intent(inout) :: mem

    if ( p == 1 ) then

      mem = memB(block1(1), block1(2)) - memB(block2(1), block2(2))
      if ( n > 2 ) then
        mem = mem +  memB(block1(2), block1(3)) - memB(block2(2), block2(3))
      end if

    else if ( p == n ) then

      mem = memB(block1(n), block1(n-1)) - memB(block2(n), block2(n-1))
      if ( n > 2 ) then
        mem = mem +  memB(block1(n-1), block1(n-2)) - memB(block2(n-1), block2(n-2))
      end if

    else

      mem = memC(block1(p-1), block1(p), block1(p+1)) - memC(block2(p-1), block2(p), block2(p+1))
      if ( p > 2 ) then
        mem = mem + memB(block1(p-1), block1(p-2)) - memB(block2(p-1), block2(p-2))
      end if
      if ( p + 1 < n ) then
        mem = mem + memB(block1(p+1), block1(p+2)) - memB(block2(p+1), block2(p+2))
      end if

    end if

  contains

    !< Memory for boundary blocks
    pure function memB(p1, p2) result(p)
      integer, intent(in) :: p1, p2
      integer :: p
      p = p1 * (p1 + p2 * 2)
    end function memB

    !< Memory for center blocks
    pure function memC(p0, p1, p2) result(p)
      integer, intent(in) :: p0, p1, p2
      integer :: p
      p = p1 * (p1 + p0 * 2 + p2 * 2)
    end function memC

  end subroutine diff_mem

  !< Calculate the difference in approximate flops for two BTD
  !!   perf(block1) - perf(block2) = d
  !! If d > 0 block1 is slower than block2
  subroutine diff_perf(p, n, block1, block2, perf)
    use precision, only: dp
    integer, intent(in) :: p
    integer, intent(in) :: n
    integer, intent(in) :: block1(n), block2(n)
    integer, intent(inout) :: perf
    real(dp) :: one_third = 0.33333333333333333333333_dp
    real(dp) :: pf

    if ( p == 1 ) then

      pf = p_inv(block1(1)) - p_inv(block2(1)) + &
          p_inv(block1(2)) - p_inv(block2(2)) + &
          p_mm(block1(1), block1(2)) - p_mm(block2(1), block2(2))

      if ( 2 < n ) then
        pf = pf + p_mm(block1(2), block1(3)) - p_mm(block2(2), block2(3))
      end if

    else if ( p == n ) then

      pf = p_inv(block1(n)) - p_inv(block2(n)) + &
          p_inv(block1(n-1)) - p_inv(block2(n-1)) + &
          p_mm(block1(n), block1(n-1)) - p_mm(block2(n), block2(n-1))

      if ( 2 < n ) then
        pf = pf + p_mm(block1(n-1), block1(n-2)) - p_mm(block2(n-1), block2(n-2))
      end if
      
    else

      pf = p_inv(block1(p-1)) - p_inv(block2(p-1)) + &
          p_mm(block1(p), block1(p-1)) - p_mm(block2(p), block2(p-1)) + &
          p_inv(block1(p)) - p_inv(block2(p)) + & ! diagonal
          p_inv(block1(p+1)) - p_inv(block2(p+1)) + &
          p_mm(block1(p), block1(p+1)) - p_mm(block2(p), block2(p+1))

      if ( 2 < p ) then
        pf = pf + p_mm(block1(p-1), block1(p-2)) - p_mm(block2(p-1), block2(p-2))
      end if

      if ( p + 1 < n ) then
        pf = pf + p_mm(block1(p+1), block1(p+2)) - p_mm(block2(p+1), block2(p+2))
      end if
      
    end if

    ! This should remove possible overflows
    ! We could essentially also just return +1/0/-1
    ! But perhaps we can use the value to something useful?
    perf = nint(sign(abs(pf) ** one_third, pf))
      
  contains

    !< Order of inversion algorithm
    pure function p_inv(n) result(p)
      use precision, only: dp
      integer, intent(in) :: n
      real(dp) :: p
      p = real(n, dp) * real(n, dp) * real(n, dp)
    end function p_inv

    !< Order of matrix multiplication, here the matrices are [m,m] * [m,n] * 2 (upper and lower)
    pure function p_mm(m, n) result(p)
      use precision, only: dp
      integer, intent(in) :: m, n
      real(dp) :: p
      p = real(m, dp) * real(m, dp) * real(n * 2, dp)
    end function p_mm

  end subroutine diff_perf
  
  function valid_tri(no,mm_col,nblocks,blocks,last_block) result(val)
    integer, intent(in) :: no, mm_col(2,no)
    integer, intent(in) :: nblocks, blocks(nblocks), last_block
    integer :: val
    ! Local variables
    integer :: i, N, Nm1, Np1, col_min, col_max
    logical :: first

    if ( last_block > 0 ) then
      
      ! We do not allow something to not end in the requested number of orbitals.
      if ( blocks(nblocks) /= last_block ) then
        val = NONVALID_SIZE
        return
      end if
      
      if ( blocks(1) == 0 ) then
        val = NONVALID_SIZE
        return
      end if
      
    end if

    ! Default to a valid BTD
    val = VALID

    ! Check that every element is contained in the 
    ! tri-diagonal matrix...
    first = .true.
    Nm1 = 1
    Np1 = blocks(1) + blocks(2)

    ! Find min/max for first block
    call find_min_max(no, mm_col, 1, blocks(1), col_min, col_max)
    if ( col_min < Nm1 ) then
      ! If this ever occur it suggests that the 
      ! sparsity pattern is not fully symmetric !
      i = 1
      call print_error(first)
      val = - 1
    else if ( Np1 < col_max ) then
      i = 1
      call print_error(first)
      val = - 1
    end if

    ! Initialize loop counters
    N = blocks(1) + 1
    do i = 2 , nblocks - 1

      ! Find min/max for this block
      call find_min_max(no, mm_col, N, Np1, col_min, col_max)

      ! Update the size of the block after this
      Np1 = Np1 + blocks(i+1)

      if ( col_min < Nm1 ) then
        ! If this ever occur it suggests that the 
        ! sparsity pattern is not fully symmetric !
        call print_error(first)
        val = - i
        
      else if ( Np1 < col_max ) then
        call print_error(first)
        val = - i
        
      end if

      ! Update loop
      N = N + blocks(i)
      
      ! Update the previous block
      Nm1 = Nm1 + blocks(i-1)

    end do

    ! Find min/max for the last block
    call find_min_max(no, mm_col, N, Np1, col_min, col_max)
    if ( col_min < Nm1 ) then
      ! If this ever occur it suggests that the 
      ! sparsity pattern is not fully symmetric !
      i = nblocks
      call print_error(first)
      val = - nblocks
      
    else if ( Np1 < col_max ) then
      i = nblocks
      call print_error(first)
      val = - nblocks
      
    end if
    if ( val /= valid ) return

    ! Update total number of orbitals
    N = N + blocks(nblocks) - 1

    ! Size of tri-matrix is already calculated
    ! Easy check, if the number of rows
    ! does not sum up to the total number of rows.
    ! Then it must be invalid...
    if ( N /= no ) then
      val = NONVALID_SIZE
    end if
    
  contains

    pure subroutine find_min_max(no, mm_col, N1, N2, col_min, col_max)
      integer, intent(in) :: no, mm_col(2,no), N1, N2
      integer, intent(inout) :: col_min, col_max
      integer :: ir
      col_min = mm_col(1,N1)
      col_max = mm_col(2,N1)
      do ir = N1 + 1, N2
        col_min = min(col_min, mm_col(1,ir))
        col_max = max(col_max, mm_col(2,ir))
      end do
    end subroutine find_min_max

    subroutine print_error(first)
      use parallel, only : IONode
      logical, intent(inout) :: first
      if ( IONode ) then
        if ( first ) then
          first = .false.
          write(*,'(a)') 'BTD: Found non-symmetric matrix!'
          write(*,'(a)') '     If you are not using delta methods you are probably doing something wrong!'
          write(*,'(a)') ' block: block row is located in'
          write(*,'(a)') ' N    : size of block'
          write(*,'(a)') ' min_C: minimum element row connects to'
          write(*,'(a)') ' min_B: minimum element in the block'
          write(*,'(a)') ' max_B: maximum element in the block'
          write(*,'(a)') ' max_C: maximum element row connects to'
          write(*,'(6a8)') 'block', 'N', 'min_C', 'min_B', 'max_B', 'max_C'
        end if
        write(*,'(6i8)') i, blocks(i), col_min, Nm1, Np1, col_max
     end if
    end subroutine print_error

  end function valid_tri

end module m_ts_rgn2trimat

