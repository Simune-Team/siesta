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

module m_ts_Sparsity2TriMat

  ! we heavily utilise the routines in this module
  use m_ts_method
  
  implicit none

  private

  public :: ts_Sparsity2TriMat

  integer, parameter :: VALID = 0
  integer, parameter :: NONVALID_SIZE = 1
  integer, parameter :: NONVALID_ELEMENT_CONTAIN = 2
  integer, parameter :: NONVALID_TS_ELECTRODE = 3
  
contains

  ! IF parts == 0 will create new partition
  subroutine ts_Sparsity2TriMat(dit,sp,parts,n_part)

    use class_OrbitalDistribution
    use class_Sparsity
    use create_Sparsity_Union
    use parallel, only : IONode, Node, Nodes
#ifdef MPI
    use mpi_siesta
#endif
    use alloc, only : re_alloc, de_alloc
    use m_ts_electype
    use m_ts_options, only : N_Elec, Elecs
    use m_ts_options, only : opt_TriMat_method
    use m_ts_options, only : IsVolt

    ! the distribution
    type(OrbitalDistribution), intent(inout) :: dit
    ! The sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! The sizes of the parts in the tri-diagonal matrix
    integer, intent(out) :: parts
    integer, pointer :: n_part(:)
    ! Local variables
    integer, pointer :: guess_part(:) => null()
    integer :: i, guess_parts
    logical :: copy_first
#ifdef MPI
    integer :: MPIerror
#endif

    if ( nrows_g(sp) <= 3 ) then
       call die('Erroneous sparsity pattern, only 3 orbitals')
    end if

    ! Establish a guess on the partition of the tri-diagonal 
    ! matrix...
    call re_alloc(guess_part, 1, nrows_g(sp), &
         routine='tsSp2TM', name='guess_part')
    call re_alloc(n_part    , 1, nrows_g(sp), &
         routine='tsSp2TM', name='n_part')
    guess_part(:) = 0

    ! We initialize to the standard 3-tri-diagonal matrix
    call set_3TriMat(nrows_g(sp),parts,n_part)
    ! If the first one happens to be the best partition, 
    ! but non-valid, we need to make sure to overwrite it
    copy_first = ( ts_valid_tri(sp,parts, n_part) /= VALID ) 

    ! TODO
    ! create array containing max-min for each orbital
    ! this will speed up this routine greatly!!!!
    
    ! We loop over all possibilities from the first part having size
    ! 2 up to and including total number of orbitals in the 
    ! In cases of MPI we do it distributed (however, the collection routine
    ! below could be optimized)
    do i = 2 , nrows_g(sp) / 10 , Nodes

       ! Make new guess...
       call guess_TriMat(sp,i,guess_parts,guess_part)

       ! If not valid tri-pattern, simply jump...
       if ( ts_valid_tri(sp,guess_parts, guess_part) /= VALID ) cycle

       call full_even_out_parts(opt_TriMat_method, &
            sp,guess_parts,guess_part)

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
#endif

    ! Shrink to the found parts
    call re_alloc(n_part,1, parts, copy=.true., shrink=.true., &
         routine='tsSp2TM', name='n_part')
    call de_alloc(guess_part,routine='tsSp2TM',name='guess_part')

    if ( parts < 3 ) then
       if ( IONode ) then 
          write(*,'(a)') 'Could not determine an optimal tri-diagonalization &
               &partition'
          if ( parts > 0 ) then
             write(*,'(a,i0)') 'Found: ',parts
             write(*,'(1000000(tr1,i0))') n_part
          else
             write(*,'(a)') 'None found...'
          end if
       end if
       call re_alloc(n_part, 1, 3, routine='tsSp2TM',name='n_part')
       call die('Not yet implemented')
       call set_3TriMat(nrows_g(sp),parts,n_part)

    end if

    ! The parts now have almost the same size and we will check that it
    ! is a valid thing, if not, we will revert to the other method of
    ! creating the tri-diagonal sparsity pattern
    ! We do not expect this to fail. After all we check that we can
    ! even out the partitions.
    ! The most probable thing is that the electrodes are not
    ! contained in the first two parts.
    if ( ts_valid_tri(sp, parts, n_part) /= VALID ) then
       i = nrows_g(sp) - no_Buf
       write(*,'(a,i0)') 'TranSIESTA system size: ',i
       write(*,'(a,i0)') 'Current parts: ',parts
       write(*,'(10000000(tr1,i0))') n_part
       i = ts_valid_tri(sp, parts, n_part)
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

  contains 

    recursive subroutine select_better(method, parts,n_part, guess_parts, guess_part)
      use m_ts_tri_scat, only : ts_needed_mem

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

         if ( guess_parts == parts ) then
            copy = faster_parts(parts,n_part,guess_part)
         end if
         ! This will take the correct value of true for the above check
         ! and this
         copy = copy .or. guess_parts > parts

      else if ( method == 1 ) then

         ! We optimize for memory, i.e. we check for number of elements
         ! in this regard we also check whether we should allocate
         ! a work-array in case of bias calculations.
         call ts_needed_mem(guess_parts,guess_part, guess_work)
         call ts_needed_mem(parts, n_part, part_work)
         
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

  end subroutine ts_Sparsity2TriMat

  subroutine guess_TriMat(sp,first_part,parts,n_part)
    use class_Sparsity
    use alloc, only: re_alloc
    use m_ts_method, only : no_Buf

    type(Sparsity), intent(inout) :: sp
    integer, intent(in) :: first_part
    integer, intent(out) :: parts
    integer, pointer :: n_part(:)

    ! Local variables
    integer :: N

    if ( first_part >= nrows_g(sp) - no_Buf ) &
         call die('Not allowed to do 1 tri-diagonal part')

    parts = 1
    n_part(1) = first_part
    N = n_part(1)
    do while ( N < nrows_g(sp) - no_Buf )
       parts = parts + 1
       if ( parts > size(n_part) ) then
          call die('Size error when guessing the tri-mat size')
       end if
       call guess_next_part_size(sp, parts, parts, n_part)
       N = N + n_part(parts)
    end do

  end subroutine guess_TriMat

  function calc_nnzs(parts,n_part) result(nnzs)
    integer, intent(in) :: parts
    integer, intent(in) :: n_part(parts)
    integer :: nnzs, i
    ! Calculate size of the tri-diagonal matrix
    nnzs = n_part(parts)**2
    do i = 1 , parts - 1
       nnzs = nnzs + n_part(i)*( n_part(i) + 2 * n_part(i+1) )
    end do
  end function calc_nnzs

  function faster_parts(parts,n_part,guess_part) result(faster)
    use precision, only: dp, i8b
    integer, intent(in) :: parts
    integer, intent(in) :: n_part(parts)
    integer, intent(in) :: guess_part(parts)

    integer :: i, diff
    integer(i8b) :: guess_N, part_N
    logical :: faster

    ! the faster must be the one with lowest total elemental
    ! operations
    guess_N = guess_part(1) ** 3
    part_N  = n_part(1)     ** 3
    diff = part_N - guess_N
    do i = 2 , parts - 1

       ! we suspect that inverting the matrix scales
       ! as the minor matrix inversions/solves
       ! I.e. we do not count the matrix multiplications

       guess_N =           guess_part(i) * guess_part(i-1) ** 2
       guess_N = guess_N + guess_part(i) ** 3
       guess_N = guess_N + guess_part(i) * guess_part(i+1) ** 2

       part_N  =           n_part(i)     * n_part(i-1)     ** 2
       part_N  = part_N  + n_part(i)     ** 3
       part_N  = part_N  + n_part(i)     * n_part(i+1)     ** 2

       diff = diff + part_N - guess_N

    end do
    guess_N = guess_part(parts) ** 3
    part_N  = n_part(parts)     ** 3
    
    diff = diff + part_N - guess_N

    faster = diff > 0

  end function faster_parts


  ! We will guess the size of this (part) tri-diagonal part by
  ! searching for the size of the matrix that matches that of the previous 
  ! part.
  ! We save it in n_part(part)
  subroutine guess_next_part_size(sp,part,parts,n_part)
    use class_Sparsity
    ! The sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! the part we are going to create
    integer, intent(in) :: part, parts
    integer, intent(inout) :: n_part(parts)
    ! Local variables
    integer :: i, sRow, eRow, mcol
    
    ! We are now checking a future part
    ! Hence we must ensure that the size is what
    ! is up to the last parts size, and thats it...
    sRow = 1
    do i = 1 , part - 2
       sRow = sRow + n_part(i)
    end do
    eRow = sRow + n_part(part-1) - 1
    ! We will check in between the above selected rows and find the 
    ! difference in size...
    n_part(part) = 0
    do i = sRow, eRow
       ! this is the # of elements from the RHS of the 'part-1'
       ! part of the tridiagonal matrix and out to the last element of
       ! this row...
       mcol = max_col(sp,tsorb2sorb(i)) - eRow
       if ( n_part(part) < mcol ) then
          n_part(part) = mcol
       end if
    end do

  end subroutine guess_next_part_size

  subroutine full_even_out_parts(method, sp,parts,n_part)
    use class_Sparsity
    use m_ts_tri_scat, only : ts_needed_mem
    integer, intent(in) :: method ! the method used for creating the parts
    ! The sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! the part we are going to create
    integer, intent(in) :: parts
    integer, intent(in out) :: n_part(parts)
    ! Local variables
    integer :: o_part(parts), mem_part(parts) , i, o_mem, n_mem, idx

    if ( method == 1 ) then
       ! we have a memory determining thing
       
       do
          o_part(:) = n_part(:)
          call ts_needed_mem(parts,n_part,o_mem)
          do i = 1 , parts
             mem_part(:) = n_part(:)
             call even_out_parts(sp, parts, n_part, i)
             call ts_needed_mem(parts,n_part,n_mem)
             if ( n_mem > o_mem ) then
                ! copy back
                n_part(:) = mem_part(:)
             end if
          end do
          if ( maxval(abs(o_part-n_part)) == 0 ) exit
       end do
       
    else

       do
          o_part(:) = n_part(:)
          mem_part(:) = n_part(:)
          ! Even out from the largest one first
          do i = 1 , parts
             idx = maxloc(mem_part,dim=1)
             mem_part(idx) = 0
             call even_out_parts(sp, parts, n_part, idx)
          end do
          if ( maxval(abs(o_part-n_part)) == 0 ) exit
       end do

    end if

  end subroutine full_even_out_parts

  subroutine even_out_parts(sp,parts,n_part, n)
    use class_Sparsity
    ! The sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! the part we are going to create
    integer, intent(in) :: parts
    integer, intent(in out) :: n_part(parts)
    integer, intent(in) :: n
    ! Local variables
    integer :: copy_n_part
    integer :: sRow, eRow
    integer :: i

    if ( parts < 2 ) call die('You cannot use tri-diagonalization &
         &without having at least 2 parts')
    
    copy_n_part = n_part(n) + 1

    sRow = 1
    do i = 1 , n - 1
       sRow = sRow + n_part(i)
    end do
    eRow = sRow + n_part(n) - 1

    ! We will continue to shift columns around
    ! until we do not shift columns any more...
    do while ( n_part(n) - copy_n_part /= 0 )

       ! Copy the current partition so that we can check in the
       ! next iteration...
       copy_n_part = n_part(n)

       ! TODO - consider adding a flag for memory reduced utilization of TRI
       ! this will require the two electrodes parts to be larger
       ! than the other parts...

       if ( n == 1 ) then
          ! If we have the first part we can always shrink it
          call even_if_larger(sRow,n_part(n),n_part(n+1))
       else if ( n == parts ) then
          ! If we have the last part we can always shrink it
          call even_if_larger(eRow,n_part(n),n_part(n-1))
       else
          ! all middle parts have some requirements:
          
          ! 1. if you wish to shrink it left, then:
          !    the first row must not have any elements
          !    extending into the right part
          if ( max_col(sp,tsorb2sorb(sRow)) <= eRow ) then
             call even_if_larger(sRow,n_part(n),n_part(n-1))
          end if

          ! 2. if you wish to shrink it right, then:
          !    the last row must not have any elements
          !    extending into the left part
          if ( sRow <= min_col(sp,tsorb2sorb(eRow)) ) then
             call even_if_larger(eRow,n_part(n),n_part(n+1))
          end if

       end if

    end do

  contains

    subroutine even_if_larger(Row,p1,p2)
      integer , intent(inout) :: Row, p1, p2
      if ( p1 > p2 ) then
         p1 = p1 - 1
         p2 = p2 + 1
         Row = Row + 1
      end if
    end subroutine even_if_larger

  end subroutine even_out_parts

! Min and max column requires that the sparsity pattern
! supplied has already stripped off the buffer orbitals.
! Otherwise this will fail
  function max_col(sp,row)
    use class_Sparsity
    use geom_helper, only : UCORB
    ! The sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! the row which we will check for (in SIESTA counting)
    integer, intent(in) :: row 
    ! The result
    integer :: max_col, ptr, nr, i, j
    integer, pointer :: l_col(:)
    call attach(sp,list_col=l_col,nrows_g=nr)
    ptr     =  list_ptr(sp,row)
    max_col = -1
    do i = 1 , n_col(sp,row)
       j = ucorb(l_col(ptr+i),nr)
       if ( orb_type(j) == TYP_BUFFER ) cycle
       max_col = max(max_col,j - orb_offset(j))
    end do
    ! Check the ts-region
    if ( max_col < 1 .or. nr - no_Buf < max_col ) &
         call die('Error in TS-sparsity pattern')
  end function max_col

  function min_col(sp,row)
    use class_Sparsity
    use geom_helper, only : UCORB
    use m_ts_method, only : no_Buf
    ! The sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! the row which we will check for (in SIESTA counting)
    integer, intent(in) :: row
    ! The result
    integer :: min_col, ptr, nr, i, j
    integer, pointer :: l_col(:)
    call attach(sp,list_col=l_col,nrows_g=nr)
    ptr     =  list_ptr(sp,row)
    min_col = huge(1)
    do i = 1 , n_col(sp,row)
       j = ucorb(l_col(ptr+i),nr)
       if ( orb_type(j) == TYP_BUFFER ) cycle
       min_col = min(min_col,j - orb_offset(j))
    end do
    if ( min_col < 1 .or. nr - no_Buf < min_col ) &
         call die('Error in TS-sparsity pattern')
  end function min_col

  function valid_tri(sp,parts,n_part) result(val) 
    use class_Sparsity
    type(Sparsity), intent(inout) :: sp
    integer, intent(in) :: parts, n_part(parts)
    integer :: val
    ! Local variables
    integer :: i, ir, N, Nm1, Np1

    val = VALID
    ! Calculate the size of the tri-matrix
    N = sum(n_part)

    ! Easy check, if the number of rows
    ! does not sum up to the total number of rows.
    ! Then it must be invalid...
    if ( N /= nrows_g(sp) - no_Buf ) then
       val = NONVALID_SIZE
       return
    end if

    ! check that all parts are at least size 2
    if ( any(n_part < 2) ) then
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
          if ( Nm1 > min_col(sp,tsorb2sorb(ir)) .or. &
               max_col(sp,tsorb2sorb(ir)) > Np1 ) then
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

  end function valid_tri

  function tsorb2sorb(tsorb) result(sorb)
    integer, intent(in) :: tsorb
    integer :: sorb

    ! start with assuming that the orb is the same orb
    sorb = tsorb
    if ( orb_offset(sorb) == 0 ) return

    ! we know that we have an offset
    ! find it by easy increments
    sorb = sorb + orb_offset(sorb) ! this is just a quick-jump
    do 
       if ( sorb - orb_offset(sorb) > tsorb ) then
          sorb = sorb - 1
       else if ( sorb - orb_offset(sorb) < tsorb ) then
          sorb = sorb + 1
       else
          exit ! we have found the correct index
       end if
    end do

  end function tsorb2sorb

  ! Validation routine for the tri-diagonal splitting
  ! with a Transiesta tri-diagonal matrix.
  ! It will first check for the electrode size
  function ts_valid_tri(sp,parts,n_part) result(val)
    use class_Sparsity
    use m_ts_electype
    use m_ts_options, only: N_Elec, Elecs
    type(Sparsity), intent(inout) :: sp
    integer, intent(in) :: parts, n_part(parts)
    integer :: val
    integer :: i, idx1, idx2, j

    do i = 1 , N_Elec
       j = Elecs(i)%idx_no
       idx1 = j - orb_offset(j)

       do j = 1 , parts - 1
          idx2 = idx1 - sum(n_part(1:j))
          if ( 0 <= idx2 ) then
             if ( idx2 + TotUsedOrbs(Elecs(i)) < n_part(j+1) ) then
                val = NONVALID_TS_ELECTRODE
                return
             end if
             exit
          end if
       end do
    end do

    val = valid_tri(sp,parts,n_part)
    if ( val /= VALID ) return
    
  end function ts_valid_tri


  subroutine set_3TriMat(no_u,parts,n_part)
    use m_ts_electype
    use m_ts_options, only : N_Elec, Elecs
    integer, intent(in) :: no_u
    integer, intent(out) :: parts
    integer, intent(out) :: n_part(3)
    integer :: noTS

    noTS = no_u - no_Buf
    parts = 3
    ! we need special handling if we only have one electrode
    if ( N_Elec == 1 ) then
       n_part(1) = TotUsedOrbs(Elecs(1))
       n_part(2) = noTS - n_part(1)
       n_part(2:3) = n_part(2) / 2
       if ( sum(n_part(1:3)) /= noTS ) then
          n_part(2) = n_part(2) + noTS - sum(n_part(1:3))
       end if
    else
       n_part(1) = sum(TotUsedOrbs(Elecs(1:N_Elec-1)))
       n_part(3) = TotUsedOrbs(Elecs(N_Elec))
       n_part(2) = noTS - n_part(1) - n_part(3)
    end if

  end subroutine set_3TriMat

end module m_ts_Sparsity2TriMat

