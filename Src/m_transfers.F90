module m_transfers
!
! Transfer integer and real arrays 
!

  use mpi
  use m_comm, only: comm_t

public :: do_transfers

interface do_transfers
   module procedure do_transfers_int
   module procedure do_transfers_dp
end interface


CONTAINS

!--------------------------------------------------
   subroutine do_transfers_int(comms,data1,data2,g1,g2,mpi_comm)

     type(comm_t), intent(in), target     :: comms(:)
     integer, dimension(:), intent(in)  :: data1
     integer, dimension(:), intent(out) :: data2
     integer, intent(in)                :: g1
     integer, intent(in)                :: g2
     integer, intent(in)                :: mpi_comm

     integer                 :: basegroup, nsize1, nsize2, ierr
     integer, allocatable    :: comm_rank1(:), comm_rank2(:)


     integer :: ncomms
     integer :: i
     integer :: nrecvs_local, nsends_local
     integer, allocatable :: statuses(:,:), local_reqR(:), local_reqS(:)
     integer :: src_in_comm, dst_in_comm
     integer :: myrank1, myrank2, myrank
     type(comm_t), pointer :: c


      ! Find the rank correspondences, in case
      ! there is implicit renumbering at the time of group creation

      call  MPI_Comm_group( mpi_comm, basegroup, ierr )
      call  MPI_Comm_Rank( mpi_comm, myrank, ierr )

      call  MPI_Group_Size( g1, nsize1, ierr )
      call  MPI_Group_Size( g2, nsize2, ierr )

      allocate(comm_rank1(0:nsize1-1))
      call MPI_Group_translate_ranks( g1, nsize1, (/ (i,i=0,nsize1-1) /), &
                                      basegroup, comm_rank1, ierr )
!      print "(i4,a,10i3)", myrank, ":Ranks of g1 in base group:", comm_rank1

      allocate(comm_rank2(0:nsize2-1))
      call MPI_Group_translate_ranks( g2, nsize2, (/ (i,i=0,nsize2-1) /), &
                                      basegroup, comm_rank2, ierr )
!      print "(i4,a,10i3)", myrank,":Ranks of g2 in base group:", comm_rank2

      call mpi_group_rank(g1,myrank1,ierr)
!      print "(i4,a,2i3)", myrank,": ierr in rank1: ", ierr
      call mpi_group_rank(g2,myrank2,ierr)
!      print "(i4,a,2i3)", myrank,": ierr in rank2: ", ierr
      
!      print "(i4,a,2i3)", myrank,": Ranks in g1 and g2: ", myrank1, myrank2
!      print "(i4,a,2i3)", myrank,": g1 and g2: ", g1, g2


      ! Do the actual transfers. 
      ! This version with non-blocking communications

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
            src_in_comm = comm_rank1(c%src)
            call MPI_irecv(data2(c%i2),c%nitems,MPI_integer,src_in_comm, &
                           i,mpi_comm,local_reqR(nrecvs_local),ierr)
         endif
      enddo

      ! Post the sends
      nsends_local = 0
      do i=1,ncomms
         c => comms(i)
         if (myrank1 == c%src) then
            nsends_local = nsends_local + 1
            dst_in_comm = comm_rank2(c%dst)
            call MPI_isend(data1(c%i1),c%nitems,MPI_integer,dst_in_comm, &
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

    end subroutine do_transfers_int

!--------------------------------------------------
   subroutine do_transfers_dp(comms,data1,data2,g1,g2,mpi_comm)

     integer, parameter :: dp = selected_real_kind(10,100)

     type(comm_t), intent(in), target     :: comms(:)
     real(dp), dimension(:), pointer :: data1
     real(dp), dimension(:), pointer :: data2
     integer, intent(in)                :: g1
     integer, intent(in)                :: g2
     integer, intent(in)                :: mpi_comm

     integer                 :: basegroup, nsize1, nsize2, ierr
     integer, allocatable    :: comm_rank1(:), comm_rank2(:)


     integer :: ncomms
     integer :: i
     integer :: nrecvs_local, nsends_local
     integer, allocatable :: statuses(:,:), local_reqR(:), local_reqS(:)
     integer :: src_in_comm, dst_in_comm
     integer :: myrank1, myrank2
     type(comm_t), pointer :: c


      ! Find the rank correspondences, in case
      ! there is implicit renumbering at the time of group creation

      call  MPI_Comm_group( mpi_comm, basegroup, ierr )
      call  MPI_Group_Size( g1, nsize1, ierr )
      call  MPI_Group_Size( g2, nsize2, ierr )
      allocate(comm_rank1(0:nsize1-1))
      call MPI_Group_translate_ranks( g1, nsize1, (/ (i,i=0,nsize1-1) /), &
                                      basegroup, comm_rank1, ierr )
!      print "(a,10i3)", "Ranks of g1 in base group:", comm_rank1
      allocate(comm_rank2(0:nsize2-1))
      call MPI_Group_translate_ranks( g2, nsize2, (/ (i,i=0,nsize2-1) /), &
                                      basegroup, comm_rank2, ierr )
!      print "(a,10i3)", "Ranks of g2 in base group:", comm_rank2

      call mpi_group_rank(g1,myrank1,ierr)
      call mpi_group_rank(g2,myrank2,ierr)

      ! Do the actual transfers. 
      ! This version with non-blocking communications

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
            src_in_comm = comm_rank1(c%src)
            call MPI_irecv(data2(c%i2),c%nitems,MPI_Double_Precision,src_in_comm, &
                           i,mpi_comm,local_reqR(nrecvs_local),ierr)
         endif
      enddo

      ! Post the sends
      nsends_local = 0
      do i=1,ncomms
         c => comms(i)
         if (myrank1 == c%src) then
            nsends_local = nsends_local + 1
            dst_in_comm = comm_rank2(c%dst)
            call MPI_isend(data1(c%i1),c%nitems,MPI_Double_Precision,dst_in_comm, &
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

    end subroutine do_transfers_dp

  end module m_transfers
