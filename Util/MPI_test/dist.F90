module m_matrix
  type, public ::  matrix
     integer :: norbs
     integer :: no_l
     integer :: nnzl
     integer, pointer :: numcols(:) => null()
     integer, pointer :: cols(:)    => null()
     !         real, pointer    :: vals(:)    => null()
  end type matrix
end module m_matrix



module m_redist
public :: redistribute
CONTAINS
  subroutine redistribute(norbs,m1,bcdist,m2,pxdist,mpi_comm)

      use mpi
      use class_BlockCyclicDist
      use class_PEXSIDist
      use m_matrix, only: matrix
      
      implicit none

    integer, intent(in)       :: norbs
    type(matrix) :: m1
    type(matrix) :: m2
    type(BlockCyclicDist) :: bcdist
    type(PEXSIDist)       :: pxdist
    integer, intent(in)   :: mpi_comm

    type comm
       integer :: src, dst, i1, i2, nitems
    end type comm
    type(comm), dimension(:), allocatable, target :: comms
    type(comm), dimension(:), allocatable, target :: commsnnz
    type(comm), pointer :: c, cnnz

    integer ::  myrank1, myrank2, myid
    logical ::  proc_in_set1, proc_in_set2
    integer ::  ierr

    integer ::  i, io

      call mpi_group_rank(bcdist%data%group,myrank1,ierr)
      call mpi_group_rank(pxdist%data%group,myrank2,ierr)
      proc_in_set1 = (myrank1 /= MPI_UNDEFINED)
      proc_in_set2 = (myrank2 /= MPI_UNDEFINED)

      call mpi_comm_rank(mpi_comm,myid,ierr)

      ! Figure out the communication needs
      call analyze_comms()

      ! In preparation for the transfer, we allocate
      ! storage for the second group of processors
      ! Note that m2%numcols (and, in general, any of the 2nd set 
      ! of arrays), will not be allocated by those processors
      ! not in the second set.

      if (proc_in_set2) then
         m2%no_l = num_local_elements(pxdist,norbs,myrank2)
         allocate(m2%numcols(m2%no_l))
      endif

      call do_transfers(comms,m1%numcols,m2%numcols)

      ! Now we can figure out how many non-zeros there are
      if (proc_in_set2) then
         m2%nnzl = sum(m2%numcols(1:m2%no_l))
         allocate(m2%cols(m2%nnzl))
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

      ! Transfer the cols arrays
      call do_transfers(commsnnz,m1%cols,m2%cols)

      ! Transfer the values arrays
!!!!      call do_transfers_real(commsnnz,m1%vals,m2%vals)

      CONTAINS

!-----------------------------------------------------
   subroutine analyze_comms()

      integer, allocatable, dimension(:) :: p1, p2, isrc, idst
      integer :: ncomms

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

      if (myid == 0) then
         write(6,"(5a10)") "Orb", "p1", "i1", "p2", "i2"
      endif
      do io = 1, norbs
         p1(io) = node_handling_element(bcdist,io)
         p2(io) = node_handling_element(pxdist,io)
         isrc(io) = index_global_to_local(bcdist,io,p1(io))
         idst(io) = index_global_to_local(pxdist,io,p2(io))
         if (myid == 0) then
            if ((norbs < 1000) .or. (mod(io,12) == 0)) then
               write(6,"(5i10)") io, p1(io), isrc(io), p2(io), idst(io)
            endif
         endif
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

      if (myid == 0) then
         do i = 1, ncomms
            c => comms(i)
            print "(a,i5,a,2i5,2i7,i5)", &
                 "comm: ", i, " src, dst, i1, i2, n:", &
                 c%src, c%dst, c%i1, c%i2, c%nitems
         enddo
      endif
    end subroutine analyze_comms

!--------------------------------------------------
   subroutine do_transfers(comms,data1,data2)
     type(comm), intent(in), target :: comms(:)
     integer, dimension(:), intent(in)  :: data1
     integer, dimension(:), intent(out) :: data2


      ! Do the actual transfers. 
      ! This version with non-blocking communications

     integer :: ncomms
     integer :: i
     integer :: nrecvs_local, nsends_local
     integer, allocatable :: statuses(:,:), local_reqR(:), local_reqS(:)


     ncomms = size(comms)

      ! Some bookkeeping for the requests
      nrecvs_local = 0
      nsends_local = 0
      do i=1,ncomms
         c => comms(i)
         if (myrank2 == c%dst) then
            nrecvs_local = nrecvs_local + 1
         endif
         if (myrank1 == c%src) then
            nsends_local = nsends_local + 1
         endif
      enddo
      allocate(local_reqR(nrecvs_local))
      allocate(local_reqS(nsends_local))
      allocate(statuses(mpi_status_size,nrecvs_local))

      ! First, post the receives
      nrecvs_local = 0
      do i=1,ncomms
         c => comms(i)
         if (myrank2 == c%dst) then
            nrecvs_local = nrecvs_local + 1
            call MPI_irecv(data2(c%i2),c%nitems,MPI_integer,c%src, &
                           i,mpi_comm,local_reqR(nrecvs_local),ierr)
         endif
      enddo

      ! Post the sends
      nsends_local = 0
      do i=1,ncomms
         c => comms(i)
         if (myrank1 == c%src) then
            nsends_local = nsends_local + 1
            call MPI_isend(data1(c%i1),c%nitems,MPI_integer,c%dst, &
                        i,mpi_comm,local_reqS(nsends_local),ierr)
         endif
      enddo

      ! A former loop of waits can be substituted by a "waitall",
      ! with every processor keeping track of the actual number of 
      ! requests in which it is involved.

      ! Should we wait also on the sends?

      call MPI_waitall(nrecvs_local, local_reqR, statuses, ierr)


      ! This barrier is needed, I think
      call MPI_Barrier(mpi_comm,ierr)

      deallocate(local_reqR, local_reqS, statuses)

    end subroutine do_transfers

  end subroutine redistribute

end module m_redist

      program dist

        ! Redistribution of orbital data.
        ! Two different distributions: 
        !   bcdist: block-cyclic (as in Siesta)
        !   pxdist: one block per processor, with fat last block (as in PEXSI)

      use mpi
      use class_BlockCyclicDist
      use class_PEXSIDist
      use m_matrix, only: matrix
      use m_redist, only: redistribute

      implicit none

      integer nprocs, i, j, ierr

      integer :: group_world, group1, group2
      integer :: nprocs1, nprocs2, bs, pbs, norbs, io, ncomms
      integer :: ig, ibeg, iend

      real    :: x


      integer ::  myrank1, myrank2, myid
      logical :: proc_in_set2, proc_in_set1
      integer, allocatable :: ranks(:)

      type(BlockCyclicDist) :: bcdist
      type(PEXSIDist)       :: pxdist

      type(matrix) :: m1, m2
         
!--------------------------------
      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

      if (myid  == 0) then
         print *, "Using ", nprocs, " procs in World."
         write(*,fmt="(a)",advance="no") "Enter number of procs in first set: "
         read *, nprocs1
         write(*,fmt="(a)",advance="no") "Enter number of procs in 2nd set: "
         read *, nprocs2
         write(*,fmt="(a)",advance="no") "Enter number of orbitals: "
         read *, norbs
         write(*,fmt="(a)",advance="no") "Enter blocksize for default dist: "
         read *, bs
      endif
      call MPI_Bcast(nprocs1,1,mpi_integer,0,MPI_COMM_WORLD,ierr)     
      call MPI_Bcast(nprocs2,1,mpi_integer,0,MPI_COMM_WORLD,ierr)     
      call MPI_Bcast(norbs,1,mpi_integer,0,MPI_COMM_WORLD,ierr)     
      call MPI_Bcast(bs,1,mpi_integer,0,MPI_COMM_WORLD,ierr)     

      call MPI_COMM_GROUP(MPI_COMM_WORLD, group_world, ierr) 

      ! New group, just for cleanliness
      allocate(ranks(nprocs1))
      do i = 1, nprocs1
         ranks(i) = i-1
      end do
      call MPI_Group_incl(group_world, nprocs1, ranks, group1, ierr)
      call newDistribution(bs,group1,bcdist,"bc dist")
      deallocate(ranks)

      ! New group, just for cleanliness
      allocate(ranks(nprocs2))
      do i = 1, nprocs2
         ranks(i) = i-1
      end do
      call MPI_Group_incl(group_world, nprocs2, ranks, group2, ierr)
      pbs = norbs/nprocs2
      call newDistribution(pbs,group2,pxdist,"px dist")
      deallocate(ranks)

      call mpi_group_rank(bcdist%data%group,myrank1,ierr)
      call mpi_group_rank(pxdist%data%group,myrank2,ierr)
      proc_in_set1 = (myrank1 /= MPI_UNDEFINED)
      proc_in_set2 = (myrank2 /= MPI_UNDEFINED)

      ! Create source matrix
      if (proc_in_set1) then
         m1%norbs = norbs
         m1%no_l = num_local_elements(bcdist,norbs,myrank1)
         allocate(m1%numcols(m1%no_l))
         do io = 1, m1%no_l
            call random_number(x)
            m1%numcols(io) = x*(0.4*norbs) + 1
         enddo
         m1%nnzl = sum(m1%numcols(1:m1%no_l))
         allocate(m1%cols(m1%nnzl))
         do j = 1, m1%nnzl
            call random_number(x)
            m1%cols(j) = x*norbs + 1
         enddo
      endif

      call redistribute(norbs,m1,bcdist,m2,pxdist,mpi_comm_world)


      if (proc_in_set1) then
         ibeg = 1
         do io = 1, m1%no_l
            iend = ibeg + m1%numcols(io) - 1
            ig = index_local_to_global(bcdist,io,myrank1)
            print "(a,i4,a,2i5,2x,i5,2x,10i3)", "Src: ", myrank1, &
                 " il, ig, ncols, cols: ", io, ig, &
                 m1%numcols(io), m1%cols(ibeg:iend)
            ibeg = iend + 1
         enddo
      endif

      call MPI_Barrier(mpi_comm_world,ierr)

      if (proc_in_set2) then
         ibeg = 1
         do io = 1, m2%no_l
            iend = ibeg + m2%numcols(io) - 1
            ig = index_local_to_global(pxdist,io,myrank2)
            print "(a,i4,a,2i5,2x,i5,2x,10i3)", "Dst: ", myrank2, &
              " il, ig, ncols, cols: ", io, ig, &
              m2%numcols(io), m2%cols(ibeg:iend)
            ibeg = iend + 1
         enddo
      endif

      call MPI_FINALIZE(ierr)

  end program dist




