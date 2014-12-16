module m_redist_spmatrix

! Declares an auxiliary structure and functions to facilitate MPI transfers
! of sparse matrices in SIESTA
!
  integer, parameter, private :: dp = selected_real_kind(10,100)

  type, public :: dp_pointer
     ! boxed array pointer type
     real(dp), pointer :: data(:) => null()
  end type dp_pointer

  type, public ::  aux_matrix

     integer :: norbs = -1
     integer :: no_l  = -1
     integer :: nnzl  = -1
     integer, pointer :: numcols(:) => null()
     integer, pointer :: cols(:)    => null()
     ! array of 1D pointers
     ! to handle multiple value fields (e.g. spin components, or S, H; DM, EDM...)
     type(dp_pointer), dimension(:), pointer :: vals(:) => null()

  end type aux_matrix

public :: redistribute_spmatrix

CONTAINS
  subroutine redistribute_spmatrix(norbs,m1,dist1,m2,dist2,mpi_comm)

      use mpi, only: mpi_group_rank, mpi_comm_rank
      use mpi, only: MPI_UNDEFINED, mpi_integer
      use class_Dist
      use m_comm,      only: comm_t
      use m_transfers, only: do_transfers
      use alloc,       only: re_alloc, de_alloc

      implicit none

    integer, intent(in)       :: norbs   ! Overall number of rows
    type(aux_matrix) :: m1               ! Source matrix
    type(aux_matrix) :: m2               ! Destination matrix -- it is allocated
    type(Dist) :: dist1, dist2           ! Distributions
    integer, intent(in)   :: mpi_comm    ! Umbrella Communicator

    type(comm_t), dimension(:), allocatable, target :: comms
    type(comm_t), dimension(:), allocatable, target :: commsnnz
    type(comm_t), pointer :: c, cnnz

    integer ::  myrank1, myrank2, myid
    logical ::  proc_in_set1, proc_in_set2
    integer ::  ierr

    integer ::  i, io, g1, g2, j, nvals
    integer, parameter :: dp = selected_real_kind(10,100)
    real(dp), dimension(:), pointer  :: data1 => null(), data2 => null()

      g1 = group(dist1)
      g2 = group(dist2)
      call mpi_group_rank(g1,myrank1,ierr)
      call mpi_group_rank(g2,myrank2,ierr)
      proc_in_set1 = (myrank1 /= MPI_UNDEFINED)
      proc_in_set2 = (myrank2 /= MPI_UNDEFINED)

      call mpi_comm_rank(mpi_comm,myid,ierr)
!      print *, "rank, ing1?, ing2?", myid, proc_in_set1, proc_in_set2
      
      ! Figure out the communication needs
      call analyze_comms()

      ! In preparation for the transfer, we allocate
      ! storage for the second group of processors
      ! Note that m2%numcols (and, in general, any of the 2nd set 
      ! of arrays), will not be allocated by those processors
      ! not in the second set.


      if (proc_in_set2) then
         m2%norbs = norbs
         m2%no_l = num_local_elements(dist2,norbs,myrank2)
         call re_alloc(m2%numcols,1,m2%no_l,"m2%numcols","redistribute_spmatrix")
      endif

!      print *, "About to transfer numcols..."
      call do_transfers(comms,m1%numcols,m2%numcols, &
                        g1,g2,mpi_comm)

      if (proc_in_set1) then
         if (associated(m1%vals)) then
            nvals = size(m1%vals)
         else
            nvals = 0
         endif
      endif
      call MPI_Bcast(nvals,1,MPI_Integer,0,mpi_comm,ierr)
!      print *, "rank, nvals: ", myid, nvals
      
      ! Now we can figure out how many non-zeros there are
      if (proc_in_set2) then
         m2%nnzl = sum(m2%numcols(1:m2%no_l))
         call re_alloc(m2%cols,1,m2%nnzl,"m2%cols","redistribute_spmatrix")

         if (nvals > 0) then
            allocate(m2%vals(nvals))
            do j=1,nvals
               call re_alloc(m2%vals(j)%data,1,m2%nnzl,"m2%vals(j)%data","redistribute_spmatrix")
            enddo
         endif

      endif

      ! Generate a new comms-structure with new start/count indexes

      allocate(commsnnz(size(comms)))
      do i = 1, size(comms)
         c => comms(i)
         cnnz => commsnnz(i)

         cnnz%src = c%src
         cnnz%dst = c%dst
         if (myrank1 == c%src) then
            ! Starting position at source: previous cols plus 1
            cnnz%i1 = sum(m1%numcols(1:(c%i1-1))) + 1
            ! Number of items transmitted: total number of cols
            cnnz%nitems = sum(m1%numcols(c%i1 : c%i1 + c%nitems -1))
         endif
         if (myrank2 == c%dst) then
            ! Starting position at destination: previous cols plus 1
            cnnz%i2 = sum(m2%numcols(1 : (c%i2-1))) + 1
            ! Number of items transmitted: total number of cols
            cnnz%nitems = sum(m2%numcols(c%i2 : c%i2 + c%nitems -1))
         endif
      end do

!!$         do i = 1, size(comms)
!!$            c => commsnnz(i)
!!$            if (myrank1 == c%src) then
!!$               print "(a,i5,a,2i5,2i7,i5)", &
!!$                 "commnnz(src): ", i, " src, dst, i1, (), n:", &
!!$                 c%src, c%dst, c%i1, -1, c%nitems
!!$            endif
!!$            if (myrank2 == c%dst) then
!!$               print "(a,i5,a,2i5,2i7,i5)", &
!!$                 "commnnz(dst): ", i, " src, dst, (), i2, n:", &
!!$                 c%src, c%dst, -1, c%i2, c%nitems
!!$            endif
!!$         enddo

!      print *, "About to transfer cols..."
      ! Transfer the cols arrays
      call do_transfers(commsnnz,m1%cols,m2%cols, &
                        g1, g2, mpi_comm)

!      print *, "About to transfer values..."
      ! Transfer the values arrays
      do j=1, nvals
	 if (proc_in_set1) data1 => m1%vals(j)%data
	 if (proc_in_set2) data2 => m2%vals(j)%data
         call do_transfers(commsnnz,data1,data2, &
              g1,g2,mpi_comm)
      enddo
      nullify(data1,data2)
!      print *, "Done transfers."

      deallocate(commsnnz)
      deallocate(comms)

      CONTAINS

!-----------------------------------------------------
   subroutine analyze_comms()

      integer, allocatable, dimension(:) :: p1, p2, isrc, idst
      integer :: ncomms

      ! To turn on debug printing, set this to .true.
      logical, save :: comms_not_printed = .false. 

      ! Find the communication needs for each orbital
      ! This information is replicated in every processor
      ! (Note that the indexing functions are able to find
      !  out the information for any processor. For the
      ! block-cyclic and "pexsi" distributions, this is quite
      ! easy. For others, the underlying indexing arrays might
      ! be large...)

      ! It might not be necessary to have this in memory. It 
      ! can be done on the fly
      allocate(p1(norbs),p2(norbs),isrc(norbs),idst(norbs))

!      if (myid == 0) then
!         write(6,"(5a10)") "Orb", "p1", "i1", "p2", "i2"
!      endif
      do io = 1, norbs
         p1(io) = node_handling_element(dist1,io)
         p2(io) = node_handling_element(dist2,io)
         isrc(io) = index_global_to_local(dist1,io,p1(io))
         idst(io) = index_global_to_local(dist2,io,p2(io))
!         if (myid == 0) then
!            if ((norbs < 1000) .or. (mod(io,12) == 0)) then
!               write(6,"(5i10)") io, p1(io), isrc(io), p2(io), idst(io)
!            endif
!        endif
      enddo

      ! Aggregate communications
      ! First pass: find out how many there are, on the basis
      ! of groups of orbitals that share the same source and
      ! destination. Due to the form of the distributions, the
      ! local indexes are also correlative in that case, so we
      ! only need to check for p1 and p2. (Check whether this
      ! applies to every possible distribution...)

      ncomms = 1
      do io = 2, norbs
         if ((p1(io) /= p1(io-1)) .or. (p2(io) /= p2(io-1))) then
            ncomms = ncomms + 1
         else
            !
         endif
      enddo

      allocate(comms(ncomms))

      ! Second pass: Fill in the data structures
      ncomms = 1
      c => comms(ncomms)
      io = 1
      c%src = p1(io)
      c%dst = p2(io)
      c%i1  = isrc(io)
      c%i2  = idst(io)
      c%nitems = 1
      do io = 2, norbs
         if ((p1(io) /= p1(io-1)) .or. (p2(io) /= p2(io-1))) then
            ! end of group -- new communication
            ncomms = ncomms + 1
            c => comms(ncomms)
            c%src = p1(io)
            c%dst = p2(io)
            c%i1  = isrc(io)
            c%i2  = idst(io)
            c%nitems = 1
         else
            ! we stay in the same communication
            c%nitems = c%nitems + 1
         endif
      enddo

      if (myid == 0 .and. comms_not_printed) then
         do i = 1, ncomms
            c => comms(i)
            write(6,"(a,i5,a,2i5,2i7,i5)"), &
                 "comm: ", i, " src, dst, i1, i2, n:", &
                 c%src, c%dst, c%i1, c%i2, c%nitems
         enddo
         comms_not_printed = .false.
      endif

      deallocate(p1,p2,isrc,idst)

    end subroutine analyze_comms

  end subroutine redistribute_spmatrix

end module m_redist_spmatrix
