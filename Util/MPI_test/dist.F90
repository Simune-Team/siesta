
      program dist

      use mpi
      use class_BlockCyclicDist
      use class_PEXSIDist

      implicit none

      integer n, myid, numprocs, i, rc, ntotal, ierr

      integer :: comm_sub, group_world, group_sub, myid_sub
      integer :: nsub, bs, pbs, norbs, io, ncomms
      integer :: no_l, no_l2, ig

      integer :: nrecvs_local, nsends_local
      integer, allocatable :: statuses(:,:), local_reqR(:), local_reqS(:)

      logical :: worker
      integer, allocatable :: ranks(:)
      integer, allocatable :: data1(:)
      integer, allocatable :: data2(:)

      type(BlockCyclicDist) :: bcdist
      type(PEXSIDist)       :: pxdist

      integer, allocatable, dimension(:) :: p1, p2, isrc, idst
      integer, allocatable, dimension(:) :: src, dst, i1, i2, nitems

      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      if (myid  == 0) then
         print *, "Using ", numprocs, " procs in World."
         print *, "Enter number of subset procs:"
         read *, nsub
         print *, "Enter number of orbitals:"
         read *, norbs
         print *, "Enter blocksize for default dist:"
         read *, bs
      endif
      call MPI_Bcast(nsub,1,mpi_integer,0,MPI_COMM_WORLD,ierr)     
      call MPI_Bcast(norbs,1,mpi_integer,0,MPI_COMM_WORLD,ierr)     
      call MPI_Bcast(bs,1,mpi_integer,0,MPI_COMM_WORLD,ierr)     

      call newDistribution(bs,MPI_Comm_World,bcdist,"bc dist")

      ! New communicator, just for cleanliness
      allocate(ranks(nsub))
      do i = 1, nsub
         ranks(i) = i-1
      end do

      worker = (myid < nsub)
!
      call MPI_COMM_GROUP(MPI_COMM_WORLD, group_world, ierr) 
      call MPI_Group_incl(group_world, nsub, ranks, group_sub, ierr)
      call MPI_Comm_create(MPI_COMM_WORLD, group_sub, comm_sub, ierr)
!
      pbs = norbs/nsub
      call newDistribution(pbs,group_sub,pxdist,"px dist")

      ! As "data", we simply store the global index of the orbital
      ! in the array "data1" in the first group of processors
      no_l = num_local_elements(bcdist,norbs,myid)
      allocate(data1(no_l))
      do io = 1, no_l
         ig = index_local_to_global(bcdist,io,myid)
         data1(io) = ig
      enddo

      ! In preparation for the transfer, we allocate
      ! storage for the second group of processors
      if (worker) then
         no_l2 = num_local_elements(pxdist,norbs,myid)
         allocate(data2(no_l2))
      endif

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

      allocate(src(ncomms),dst(ncomms),i1(ncomms),i2(ncomms),nitems(ncomms))

      ! Second pass: Fill in the data structures
      ncomms = 1
      src(ncomms) = p1(1)
      dst(ncomms) = p2(1)
      i1(ncomms) = isrc(1)
      i2(ncomms) = idst(1)
      nitems(ncomms) = 1
      do io = 2, norbs
         if ((p1(io) /= p1(io-1)) .or. (p2(io) /= p2(io-1))) then
            ! end of group -- new communication
            ncomms = ncomms + 1
            src(ncomms) = p1(io)
            dst(ncomms) = p2(io)
            i1(ncomms) = isrc(io)
            i2(ncomms) = idst(io)
            nitems(ncomms) = 1
         else
            ! we stay in the same communication
            nitems(ncomms) = nitems(ncomms) + 1
         endif
      enddo

      if (myid == 0) then
         do i = 1, ncomms
            print "(a,i5,a,2i5,2i7,i5)", &
                 "comm: ", i, " src, dst, i1, i2, n:", &
                 src(i), dst(i), i1(i), i2(i), nitems(i)
         enddo
      endif

      ! Do the actual transfers, with non-blocking communications

      nrecvs_local = 0
      nsends_local = 0
      do i=1,ncomms
         if (myid == dst(i)) then
            nrecvs_local = nrecvs_local + 1
         endif
         if (myid == src(i)) then
            nsends_local = nsends_local + 1
         endif
      enddo
      allocate(local_reqR(nrecvs_local))
      allocate(local_reqS(nsends_local))
      allocate(statuses(mpi_status_size,nrecvs_local))

      ! First, post the receives
      nrecvs_local = 0
      do i=1,ncomms
         if (myid == dst(i)) then
            nrecvs_local = nrecvs_local + 1
            call MPI_irecv(data2(i2(i)),nitems(i),MPI_integer,src(i), &
                           i,mpi_comm_world,local_reqR(nrecvs_local),ierr)
         endif
      enddo
      ! Post the sends
      nsends_local = 0
      do i=1,ncomms
         if (myid == src(i)) then
            nsends_local = nsends_local + 1
            call MPI_isend(data1(i1(i)),nitems(i),MPI_integer,dst(i), &
                        i,mpi_comm_world,local_reqS(nsends_local),ierr)
         endif
      enddo

      ! A former loop of waits can be substituted by a "waitall",
      ! with every processor keeping track of the actual number of 
      ! requests in which it is involved.

      ! Should we wait also on the sends?

      call MPI_waitall(nrecvs_local, local_reqR, statuses, ierr)


      ! This barrier is needed, I think
      call MPI_Barrier(mpi_comm_world,ierr)

      if (myid == 0) then
         print *, "Checking... (no extra output if correct)..."
      endif

      if (worker) then
         do io = 1, no_l2
            call MPI_Group_RANK(group_sub , myid, ierr )
            !print *, "Node: ", myid, "il, ig: ", io, data2(io)
            ig = index_local_to_global(pxdist,io,myid)
            if (ig /= data2(io)) then
               print "(a,i4,a,2i7,a,i7)", "Node: ", myid, &
                        "il, ig: ", io, data2(io), &
                        " should be ", ig
            endif
         enddo
      endif

 30   call MPI_FINALIZE(rc)

    end program dist




