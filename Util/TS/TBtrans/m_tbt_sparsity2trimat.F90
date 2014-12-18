!
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
!
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2013, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

! Module for converting a sparsity pattern to an
! "optimal" tri-diagonal matrix...

! Routine for converting a TranSIESTA sparsity pattern to a tri-diagonal form
! It requires that the sparsity pattern is fully contained in the current
! processor.

! This module has been fully developed by Nick Papior Andersen, 2013
! nickpapior@gmail.com

! Please contact the author before utilization in other routines.

module m_tbt_region2TriMat

  ! Use regions...
  use precision, only : dp, i8b
  use m_region

  implicit none

  private

  public :: tbt_region2TriMat

  integer, parameter :: VALID = 0
  integer, parameter :: NONVALID_SIZE = 1
  integer, parameter :: NONVALID_ELEMENT_CONTAIN = 2
  integer, parameter :: NONVALID_TS_ELECTRODE = 3
  
contains

  subroutine GFGGF_needed_worksize(N_tri_part, tri_parts, &
       N_Elec, Elecs, padding, worksize)
    use m_ts_electype
    integer, intent(in) :: N_tri_part
    integer, intent(in) :: tri_parts(N_tri_part)
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    integer, intent(out) :: padding, worksize

    integer :: els, n, tn, io
    integer :: no_max, cur_n

    ! calculate the maximum electrode size
    no_max = maxval(TotUsedOrbs(Elecs))

    ! We just need to find the maximum overlap of
    ! two regions.
    ! This will give the "pushed" number of elements
    ! that is required to prevent overwriting two
    ! quantities.
    ! At minimum this will most likely be the size
    ! of the last two parts due to that being calculated
    ! last.

    els = tri_parts(N_tri_part)**2
    worksize = tri_parts(N_tri_part)
    do n = 1 , N_tri_part - 1
       els = els + tri_parts(n)*( tri_parts(n) + 2 * tri_parts(n+1) )
       worksize = max(worksize,tri_parts(n))
    end do
    worksize = worksize * no_max

    ! subtract total column size
    ! to get the first matrix element of the current processing
    ! block (with an index shift of 1, so actually previous element
    ! of what is needed)
    tn = els - sum(tri_parts(1:n) * no_max)

    cur_n = 0
    io = 1
    do n = 1 , N_tri_part

       if ( 1 < n ) &
            cur_n = cur_n + tri_parts(n-1) * tri_parts(n)
       cur_n = cur_n + tri_parts(n) ** 2
       if ( n < N_tri_part ) &
            cur_n = cur_n + tri_parts(n) * tri_parts(n+1)
       
       if ( cur_n > tn ) then
          ! we have an overlap, calculate overlap
          ! and correct tn
          ! With ">" we do not need to correct the tn initialization
          ! of element - 1 as noted above
          padding = cur_n - tn
          ! We correct the starting index of tn
          tn = tn + padding
       end if
       
       ! we need to retain the column block
       ! for the next block...
       ! in that way we can still multiply the previous
       ! block with the current block.
       if ( n > 1 ) then
          tn = tn + tri_parts(n-1) * no_max
       end if

    end do
    tn = tn + tri_parts(N_tri_part) * no_max

    ! the padding must be the excess size we have appended to the matrix
    padding = tn - els

  end subroutine GFGGF_needed_worksize

  subroutine needed_mem(N_tri_part, tri_parts, worksize)
    use m_ts_electype
    use m_tbt_options, only : N_Elec, Elecs
    integer, intent(in) :: N_tri_part
    integer, intent(in) :: tri_parts(N_tri_part)
    integer, intent(out) :: worksize

    integer :: pad, n

    ! find at which point they will cross...
    worksize = tri_parts(N_tri_part)**2
    do n = 1 , N_tri_part - 1
       worksize = worksize + &
            tri_parts(n)*( tri_parts(n) + 2 * tri_parts(n+1) )
    end do
    ! for the two arrays
    worksize = worksize * 2

    call GFGGF_needed_worksize(N_tri_part, tri_parts, &
         N_Elec, Elecs, pad, n)
    worksize = worksize + pad + n
    
  end subroutine needed_mem

  ! IF parts == 0 will create new partition
  subroutine tbt_region2TriMat(dit,sp,r,parts,n_part, last_eq, par)

    use class_OrbitalDistribution
    use class_Sparsity
    use create_Sparsity_Union
    use parallel, only : IONode, Node, Nodes
    use fdf, only : fdf_get
#ifdef MPI
    use mpi_siesta
#endif
    use alloc, only : re_alloc, de_alloc
    use m_ts_electype
    
    use m_tbt_options, only : opt_TriMat_method

    ! the distribution
    type(OrbitalDistribution), intent(inout) :: dit
    ! The sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! The region that we will create a tri-diagonal matrix on.
    type(tRegion), intent(in) :: r
    ! The sizes of the parts in the tri-diagonal matrix
    integer, intent(out) :: parts
    integer, pointer :: n_part(:)
    ! Whether we should retain the last partition to a fixed size.
    integer, intent(in) :: last_eq
    ! Whether the search should be performed in parallel or not
    logical, intent(in), optional :: par

    ! Local variables
    integer, pointer :: guess_part(:) => null()
    integer, pointer :: mm_col(:,:) => null()
    integer :: i, no, guess_parts, max_block
    ! In case of parallel
    integer :: guess_start, guess_step
    logical :: copy_first, lpar
#ifdef MPI
    integer :: MPIerror
#endif

    lpar = .true.
    if ( present(par) ) lpar = par

    ! This is the size of the regional 3-diagonal matrix
    no = r%n
    if ( no <= 3 ) then
       call die('Erroneous sparsity pattern, only 3 orbitals')
    end if

    ! Establish a guess on the partition of the tri-diagonal 
    ! matrix...
    call re_alloc(guess_part, 1, no, &
         routine='tbtR2TM', name='guess_part')
    call re_alloc(n_part    , 1, no, &
         routine='tbtR2TM', name='n_part')
    guess_part(:) = 0

    ! create array containing max-min for each ts-orbital
    call re_alloc(mm_col  , 1, 2, 1, no, &
         routine='tbtR2TM', name='mm_col')
    do i = 1 , no
       mm_col(1,i) = min_col(sp,r,r%r(i))
       mm_col(2,i) = max_col(sp,r,r%r(i))
       !print '(a,tr1,i5,tr3,2(tr1,i5),tr4,i0)','Orb: ',i,mm_col(:,i),no
    end do

    parts = 2
    n_part(1) = no / 2
    n_part(2) = no / 2 + mod(no,2)
    if ( last_eq > 0 ) then
       parts = 2
       n_part(2) = last_eq
       n_part(1) = no - last_eq
       ! Initialize the guess for the 
       call guess_TriMat_last(no,mm_col,guess_parts,guess_part,last_eq)
       if ( valid_tri(no,r,mm_col,guess_parts, guess_part,last_eq) == VALID ) then
          parts = guess_parts
          n_part(1:parts) = guess_part(1:parts)
       end if
    end if

    ! We default to parallel
    guess_start = 2 + Node
    guess_step = Nodes
    if ( .not. lpar ) then
       guess_start = 2
       guess_step = 1
    end if

    ! If the first one happens to be the best partition, 
    ! but non-valid, we need to make sure to overwrite it
    copy_first = .false. ! currently TODO THIS COULD BE A PROBLEM

    ! If the blocks are known by the user to not exceed a certain
    ! size, then we can greatly reduce the guessing step
    ! for huge systems
    max_block = no / 4
    max_block = fdf_get('TS.TriMat.Block.Max',max_block)
    max_block = fdf_get('TBT.TriMat.Block.Max',max_block)
    
    ! We loop over all possibilities from the first part having size
    ! 2 up to and including total number of orbitals in the 
    ! In cases of MPI we do it distributed (however, the collection routine
    ! below could be optimized)
    do i = guess_start , max_block , guess_step

       ! Make new guess...
       call guess_TriMat(no,mm_col,i,guess_parts,guess_part,last_eq)

       ! If not valid tri-pattern, simply jump...
       if ( valid_tri(no,r,mm_col,guess_parts, guess_part,last_eq) /= VALID ) then
          !print*,guess_part(1:guess_parts)
          cycle
       end if

       call full_even_out_parts(opt_TriMat_method, &
            no,mm_col,guess_parts,guess_part,last_eq)

       if ( last_eq > 0 .and. guess_part(guess_parts) /= last_eq ) then
          call die('Something went terribly wrong...')
       end if

       if ( copy_first ) then
          ! ensure to copy it over (the initial one was not valid)
          copy_first = .false.
          parts = guess_parts
          n_part(1:parts) = guess_part(1:parts)
          cycle
       end if
       
       call select_better(opt_TriMat_method, &
            parts,n_part, guess_parts, guess_part)
       
    end do

#ifdef MPI
    if ( lpar ) then
       ! Select the most optimal partition scheme...
       do i = 0 , Nodes - 1
          if ( i == Node ) then
             call MPI_Bcast(parts, 1, MPI_Integer, i, &
                  MPI_Comm_World, MPIerror)
             call MPI_Bcast(n_part(1), parts, MPI_Integer, i, &
                  MPI_Comm_World, MPIerror)
          else
             call MPI_Bcast(guess_parts, 1, MPI_Integer, i, &
                  MPI_Comm_World, MPIerror)
             call MPI_Bcast(guess_part(1), guess_parts, MPI_Integer, i, &
                  MPI_Comm_World, MPIerror)
             ! Only all the other nodes are allowed to check...
             ! TODO this could be made to a communication tree to limit communication
             call select_better(opt_TriMat_method, &
                  parts,n_part, guess_parts, guess_part)
          end if
       end do
    end if
#endif

    ! Shrink to the found parts
    call re_alloc(n_part,1, parts, copy=.true., shrink=.true., &
         routine='tsSp2TM', name='n_part')
    call de_alloc(guess_part,routine='tsSp2TM',name='guess_part')

    if ( parts < 2 ) then
       if ( IONode ) then 
          write(*,'(a)') 'Could not determine an optimal tri-diagonalization &
               &partition'
          write(*,'(2a)') 'Running on region: ',trim(r%name)
          if ( parts > 0 ) then
             write(*,'(a,i0)') 'Found: ',parts
             write(*,'(1000000(tr1,i0))') n_part
          else
             write(*,'(a)') 'None found...'
          end if
       end if
       call re_alloc(n_part, 1, 3, routine='tsSp2TM',name='n_part')
       call die('Not yet implemented')

    end if

    ! The parts now have almost the same size and we will check that it
    ! is a valid thing, if not, we will revert to the other method of
    ! creating the tri-diagonal sparsity pattern
    ! We do not expect this to fail. After all we check that we can
    ! even out the partitions.
    ! The most probable thing is that the electrodes are not
    ! contained in the first two parts.
    i = valid_tri(no,r,mm_col,parts, n_part,last_eq)
    if ( i /= VALID ) then
       write(*,'(2a)') 'Running on region: ',trim(r%name)
       write(*,'(a,i0)') 'TranSIESTA system size: ',no
       write(*,'(a,i0)') 'Current parts: ',parts
       write(*,'(10000000(tr1,i0))') n_part
       write(*,'(a,i0)') 'Current part size: ',sum(n_part(:))
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

    call de_alloc(mm_col,routine='tsSp2TM',name='mm_col')

  contains 

    recursive subroutine select_better(method, parts,n_part, &
         guess_parts, guess_part)

      integer, intent(in)    :: method
      integer, intent(inout) :: parts
      integer, intent(in)    :: guess_parts
      integer, intent(inout) :: n_part(max(parts,guess_parts))
      integer, intent(in)    :: guess_part(guess_parts)
      logical :: copy
      integer :: part_work, guess_work

      copy = .false.

      ! We check whether the number of elements is smaller
      ! or that the number of parts is greater (however, this should
      ! in principle always go together)
      ! If the method of optimization is memory:
      if ( method == 0 ) then

         copy = faster_parts(parts,n_part,guess_parts,guess_part)

      else if ( method == 1 ) then

         ! We optimize for memory, i.e. we check for number of elements
         ! in this regard we also check whether we should allocate
         ! a work-array in case of bias calculations.
         call needed_mem(guess_parts,guess_part, guess_work)
         call needed_mem(parts, n_part, part_work)
         
         copy = part_work > guess_work
         if ( .not. copy ) then
            ! in case the work-size is the same...
            if ( part_work == guess_work ) then
               call select_better(0, parts,n_part, guess_parts, guess_part)
            end if
         end if

      else
         call die('Unknown optimization scheme for the tri-mat')
      end if

      if ( copy ) then
         parts = guess_parts
         n_part(1:parts) = guess_part(1:parts)
      end if

    end subroutine select_better

  end subroutine tbt_region2TriMat

  subroutine guess_TriMat(no,mm_col,first_part,parts,n_part,last_eq)

    integer, intent(in) :: no, mm_col(2,no) ! number of orbitals, max,min
    integer, intent(in) :: first_part
    integer, intent(out) :: parts
    integer, intent(out) :: n_part(:)
    integer, intent(in) :: last_eq

    ! Local variables
    integer :: N

    if ( first_part > no ) &
         call die('Not allowed to do 1 tri-diagonal part')

    parts = 1
    n_part(1) = first_part
    N = n_part(1)
    do while ( N < no )
       parts = parts + 1
       if ( parts > size(n_part) ) then
          print *,'Error',parts,size(n_part)
          call die('Size error when guessing the tri-mat size')
       end if
       call guess_next_part_size(no, mm_col, parts, parts, n_part)
       N = N + n_part(parts)
       if ( last_eq > 0 ) then
          ! if a last-part was "forced" we do this here...
          if ( N + last_eq > no ) then
             ! We need to add the former part with the "too many"
             ! orbitals
             n_part(parts) = n_part(parts) + no - N
             N = no
          end if
       end if
             
    end do

    if ( last_eq > 0 ) then
       ! Correct so that we actually do contain the last_eq
       ! in the last one
       n_part(parts) = last_eq
       n_part(parts-1) = n_part(parts-1) + no - sum(n_part(1:parts))
    end if

  end subroutine guess_TriMat

  subroutine guess_TriMat_last(no,mm_col,parts,n_part,last_eq)

    integer, intent(in) :: no, mm_col(2,no) ! number of orbitals, max,min
    integer, intent(out) :: parts
    integer, intent(out) :: n_part(:)
    integer, intent(in) :: last_eq

    ! Local variables
    integer :: N

    parts = 1
    n_part(1) = last_eq
    N = n_part(1)
    do while ( N < no )
       parts = parts + 1
       if ( parts > size(n_part) ) then
          print *,'Error',parts,size(n_part)
          call die('Size error when guessing the tri-mat size')
       end if
       call guess_prev_part_size(no, mm_col, parts, parts, n_part)
       N = N + n_part(parts)
    end do

    ! Reverse
    n_part(1:parts) = n_part(parts:1:-1)

  end subroutine guess_TriMat_last

  function calc_nnzs(parts,n_part) result(nnzs)
    integer, intent(in) :: parts
    integer, intent(in) :: n_part(parts)
    integer(i8b) :: i8_part(parts)
    integer :: nnzs
    i8_part(:) = n_part(:)
    ! Calculate size of the tri-diagonal matrix
    nnzs = sum(i8_part(:)**2) + sum(i8_part(1:parts-1)*i8_part(2:parts))
  end function calc_nnzs

  function faster_parts(np,i4n_part,ng,i4g_part) result(faster)
    integer, intent(in) :: np, i4n_part(np)
    integer, intent(in) :: ng, i4g_part(ng)
    logical :: faster

    integer :: i
    integer(i8b) :: guess_N, part_N, diff
    integer(i8b) :: n_part(np), guess_part(ng)

    n_part(:) = i4n_part(:)
    guess_part(:) = i4g_part(:)

    ! We estimate the fastest algorithm
    ! by the number of operations the matrices make

    diff = 0
    do i = 1 , max(np,ng)
       part_N = 0
       if ( i < np ) then
          part_N = 2 * 5 * n_part(i) ** 2 / 3 + 4 * n_part(i+1)
          part_N = part_N * n_part(i)
       end if
       guess_N = 0
       if ( i < ng ) then
          guess_N = 2 * 5 * guess_part(i) ** 2 / 3 + 4 * guess_part(i+1) 
          guess_N = guess_N * guess_part(i)
       end if
       diff = diff + part_N - guess_N

       if ( i == 1 ) cycle

       part_N = 0
       if ( i <= np ) then
          part_N = 2 * 5 * n_part(i) ** 2 / 3 + 4 * n_part(i-1) 
          part_N = part_N * n_part(i)
       end if
       guess_N = 0
       if ( i <= ng ) then
          guess_N = 2 * 5 * guess_part(i) ** 2 / 3 + 4 * guess_part(i-1)
          guess_N = guess_N * guess_part(i)
       end if
       diff = diff + part_N - guess_N
    end do

    faster = (diff > 0)

  end function faster_parts


  ! We will guess the size of this (part) tri-diagonal part by
  ! searching for the size of the matrix that matches that of the previous 
  ! part.
  ! We save it in n_part(part)
  subroutine guess_next_part_size(no,mm_col,part,parts,n_part)
    integer, intent(in) :: no, mm_col(2,no)
    ! the part we are going to create
    integer, intent(in) :: part, parts
    integer, intent(inout) :: n_part(parts)
    ! Local variables
    integer :: i, sRow, eRow, mcol
    
    ! We are now checking a future part
    ! Hence we must ensure that the size is what
    ! is up to the last parts size, and thats it...
    eRow = 0
    do i = 1 , part - 1
       eRow = eRow + n_part(i)
    end do
    sRow = eRow - n_part(part-1) + 1

    ! We will check in between the above selected rows and find the 
    ! difference in size...
    n_part(part) = 0
    do i = sRow, eRow
       ! this is the # of elements from the RHS of the 'part-1'
       ! part of the tridiagonal matrix and out to the last element of
       ! this row...
       mcol = mm_col(2,i) - eRow
       if ( n_part(part) < mcol ) then
          n_part(part) = mcol
       end if
    end do

    ! In case there is actually no connection, we should
    ! force the next-part to be 1!
    if ( n_part(part) == 0 ) then
       n_part(part) = 1
    end if

  end subroutine guess_next_part_size

  ! We will guess the size of this (part) tri-diagonal part by
  ! searching for the size of the matrix that matches that of the previous 
  ! part.
  ! We save it in n_part(part)
  subroutine guess_prev_part_size(no,mm_col,part,parts,n_part)
    integer, intent(in) :: no, mm_col(2,no)
    ! the part we are going to create
    integer, intent(in) :: part, parts
    integer, intent(inout) :: n_part(parts)
    ! Local variables
    integer :: i, sRow, eRow, mcol
    
    ! We are now checking a future part
    ! Hence we must ensure that the size is what
    ! is up to the last parts size, and thats it...
    eRow = no
    if ( part > 2 ) then
       do i = 1 , part - 2
          eRow = eRow - n_part(i)
       end do
    end if
    sRow = eRow - n_part(part-1) + 1

    ! We will check in between the above selected rows and find the 
    ! difference in size...
    n_part(part) = 0
    do i = sRow, eRow
       ! this is the # of elements from the RHS of the 'part-1'
       ! part of the tridiagonal matrix and out to the last element of
       ! this row...
       mcol = sRow - mm_col(1,i)
       if ( n_part(part) < mcol ) then
          n_part(part) = mcol
       end if
    end do

    ! In case there is actually no connection, we should
    ! force the next-part to be 1!
    if ( n_part(part) == 0 ) then
       n_part(part) = 1
    end if

  end subroutine guess_prev_part_size

  subroutine full_even_out_parts(method,no,mm_col,parts,n_part, last_eq)
    integer, intent(in) :: method ! the method used for creating the parts
    integer, intent(in) :: no, mm_col(2,no)
    ! the part we are going to create
    integer, intent(in) :: parts
    integer, intent(inout) :: n_part(parts)
    integer, intent(in) :: last_eq
    ! Local variables
    integer :: o_part(parts), mem_part(parts) , i, o_mem, n_mem, idx

    if ( method == 1 ) then
       ! we have a memory determining thing
       
       do
          o_part(:) = n_part(:)
          call needed_mem(parts,n_part,o_mem)
          do i = 1 , parts
             mem_part(:) = n_part(:)
             call even_out_parts(no, mm_col, parts, n_part, i, last_eq)
             call needed_mem(parts,n_part,n_mem)
             if ( n_mem > o_mem ) then
                ! copy back
                n_part(:) = mem_part(:)
             end if
          end do
          if ( maxval(abs(o_part-n_part)) == 0 ) exit
       end do
       
    else

       do
          o_part(:)   = n_part(:)
          mem_part(:) = n_part(:)
          ! Even out from the largest one first
          do i = 1 , parts
             idx = maxloc(mem_part,dim=1)
             mem_part(idx) = 0
             call even_out_parts(no, mm_col, parts, n_part, idx, last_eq)
          end do
          if ( maxval(abs(o_part-n_part)) == 0 ) exit
       end do

    end if

  end subroutine full_even_out_parts

  subroutine even_out_parts(no,mm_col,parts,n_part, n, last_eq)
    integer, intent(in) :: no, mm_col(2,no)
    ! the part we are going to create
    integer, intent(in) :: parts
    integer, intent(in out) :: n_part(parts)
    integer, intent(in) :: n, last_eq
    ! Local variables
    integer :: copy_n_part
    integer :: sRow, eRow
    integer :: i

    if ( parts < 2 ) call die('You cannot use tri-diagonalization &
         &without having at least 2 parts')
    
    if ( last_eq > 0 .and. n >= parts - 1 ) then
       ! We need the last one to be of a certain
       ! size.
       return

    end if
    
    if ( parts == 2 ) then

       i = 0
       copy_n_part = 0
       do while ( n_part(n) - copy_n_part /= 0 )
          
          ! Copy the current partition so that we can check in the
          ! next iteration...
          copy_n_part = n_part(n)

          ! TODO - consider adding a flag for memory reduced utilization of TRI
          ! this will require the two electrodes parts to be larger
          ! than the other parts...

          if ( n == 1 ) then
             ! If we have the first part we can always shrink it
             call even_if_larger(i,n_part(n),n_part(n+1),sign= 1)
          else if ( n == parts ) then
             ! If we have the last part we can always shrink it
             call even_if_larger(i,n_part(n),n_part(n-1),sign=-1)
          end if
       end do

       return

    end if

    ! We do not allow to diminish the edges
    ! This is because we need additional checks of their
    ! regions, hence, we let the central parts partition out
    ! the regions
    if ( n == 1 .or. n == parts ) return

    sRow = 1
    do i = 1 , n - 1
       sRow = sRow + n_part(i)
    end do
    eRow = sRow + n_part(n) - 1

    ! We will continue to shift columns around
    ! until we do not shift columns any more...
    copy_n_part = 0
    do while ( n_part(n) - copy_n_part /= 0 )

       ! Copy the current partition so that we can check in the
       ! next iteration...
       copy_n_part = n_part(n)

       ! TODO - consider adding a flag for memory reduced utilization of TRI
       ! this will require the two electrodes parts to be larger
       ! than the other parts...

       ! 1. if you wish to shrink it left, then:
       !    the first row must not have any elements
       !    extending into the right part
       if ( mm_col(2,sRow) <= eRow ) then
          call even_if_larger(sRow,n_part(n),n_part(n-1),sign= 1)
       end if

       ! 2. if you wish to shrink it right, then:
       !    the last row must not have any elements
       !    extending into the left part
       if ( sRow <= mm_col(1,eRow) ) then
          call even_if_larger(eRow,n_part(n),n_part(n+1),sign=-1)
       end if

    end do

  contains

    subroutine even_if_larger(Row,p1,p2,sign)
      integer, intent(in) :: sign
      integer , intent(inout) :: Row, p1, p2
      if ( p1 > p2 ) then
         p1 = p1 - 1
         p2 = p2 + 1
         Row = Row + sign * 1
      end if
    end subroutine even_if_larger

  end subroutine even_out_parts

! Min and max column requires that the sparsity pattern
! supplied has already stripped off the buffer orbitals.
! Otherwise this will fail
  function max_col(sp,r,row)
    use class_Sparsity
    use geom_helper, only : UCORB
    ! The sparsity pattern
    type(Sparsity), intent(inout) :: sp
    type(tRegion), intent(in) :: r
    ! the row which we will check for (in TranSIESTA counting)
    integer, intent(in) :: row 
    ! The result
    integer :: max_col, ptr, nr, j
    integer, pointer :: l_col(:), l_ptr(:), ncol(:)
    call attach(sp,n_col=ncol,list_ptr=l_ptr,list_col=l_col,nrows_g=nr)
    max_col = region_pivot(r,row)
    do ptr = l_ptr(row) + 1 , l_ptr(row) + ncol(row)
       j = region_pivot(r,ucorb(l_col(ptr),nr))
       if ( j > 0 ) then
          max_col = max(max_col,j)
       end if
    end do
  end function max_col

  function min_col(sp,r,row)
    use class_Sparsity
    use geom_helper, only : UCORB
    ! The sparsity pattern
    type(Sparsity), intent(inout) :: sp
    type(tRegion), intent(in) :: r
    ! the row which we will check for (in TranSIESTA counting)
    integer, intent(in) :: row
    ! The result
    integer :: min_col, ptr, nr, j
    integer, pointer :: l_col(:), l_ptr(:), ncol(:)
    call attach(sp,n_col=ncol,list_ptr=l_ptr,list_col=l_col,nrows_g=nr)
    min_col = region_pivot(r,row)
    do ptr = l_ptr(row) + 1 , l_ptr(row) + ncol(row)
       j = region_pivot(r,ucorb(l_col(ptr),nr))
       if ( j > 0 ) then
          min_col = min(min_col,j)
       end if
    end do
  end function min_col

  function valid_tri(no,r,mm_col,parts,n_part,last_eq) result(val) 
    integer, intent(in) :: no, mm_col(2,no)
    type(tRegion), intent(in) :: r
    integer, intent(in) :: parts, n_part(parts), last_eq
    integer :: val
    ! Local variables
    integer :: i, ir, N, Nm1, Np1

    val = VALID
    ! Calculate the size of the tri-matrix
    N = sum(n_part)

    ! Easy check, if the number of rows
    ! does not sum up to the total number of rows.
    ! Then it must be invalid...
    if ( N /= no ) then
       val = NONVALID_SIZE
       return
    end if

    ! check that all parts are at least size 2
    if ( any(n_part < 2 ) ) then
       val = NONVALID_SIZE
       return
    end if

    ! Check that every element is contained in the 
    ! tri-diagonal matrix...
    N = 1
    Nm1 = 1
    Np1 = n_part(1)
    do i = 1 , parts
       
       if ( i < parts ) then
          ! Update the size of the part after this
          Np1 = Np1 + n_part(i+1)
       end if
       
       do ir = N , N + n_part(i) - 1
          if ( mm_col(1,ir) < Nm1 .or. &
               mm_col(2,ir) > Np1 ) then
             print *,mm_col(1,ir),Nm1,Np1,mm_col(2,ir)
             val = - ir 
             return
          end if
       end do
       
       ! Update loop
       N = N + n_part(i)
       
       if ( i > 1 ) then
          ! Update the previous part
          Nm1 = Nm1 + n_part(i-1)
       end if
    end do

    if ( last_eq > 0 ) then
       ! We do not allow something to not end in the requested number 
       ! of orbitals.
       if ( n_part(parts) /= last_eq ) then
          val = NONVALID_SIZE
          return
       end if
    end if


  end function valid_tri

end module m_tbt_region2TriMat

