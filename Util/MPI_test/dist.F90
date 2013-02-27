
      program dist

      use mpi
      use class_BlockCyclicDist
      use class_PEXSIDist

      implicit none

      integer n, myid, numprocs, i, rc, ntotal, ierr

      integer :: comm_sub, group_world, group_sub, myid_sub
      integer :: nsub, bs, pbs, norbs, io
      integer :: no_l, no_l2, ig
      integer :: status(mpi_status_size)

      logical :: worker
      integer, allocatable :: ranks(:)
      integer, allocatable :: data1(:)
      integer, allocatable :: data2(:)
      type(BlockCyclicDist) :: bcdist
      type(PEXSIDist)       :: pxdist

      integer, allocatable, dimension(:) :: p1, p2, isrc, idst
      integer, allocatable, dimension(:) :: requestR, requestS

      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      if (myid  == 0) then
         print *, "Using ", numprocs, " procs in World."
         print *, "Enter number of subset procs:"
         read *, nsub
         ! Avoid complications for now
!         nsub = numprocs
         print *, "Enter number of orbitals:"
         read *, norbs
         print *, "Enter blocksize for default dist:"
         read *, bs
      endif
      call MPI_Bcast(nsub,1,mpi_integer,0,MPI_COMM_WORLD,ierr)     
      call MPI_Bcast(norbs,1,mpi_integer,0,MPI_COMM_WORLD,ierr)     
      call MPI_Bcast(bs,1,mpi_integer,0,MPI_COMM_WORLD,ierr)     

      allocate(p1(norbs),p2(norbs),isrc(norbs),idst(norbs))
      allocate(requestS(norbs),requestR(norbs))

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

      no_l = num_local_elements(bcdist,norbs,myid)
      allocate(data1(no_l))
      do io = 1, no_l
         ig = index_local_to_global(bcdist,io,myid)
         data1(io) = ig
      enddo

      if (worker) then
         no_l2 = num_local_elements(pxdist,norbs,myid)
         allocate(data2(no_l2))
      endif

      if (myid == 0) then
         write(6,"(5a10)") "Orb", "p1", "i1", "p2", "i2"
      endif
      do io = 1, norbs
         p1(io) = node_handling_element(bcdist,io)
         p2(io) = node_handling_element(pxdist,io)
         isrc(io) = index_global_to_local(bcdist,io,p1(io))
         idst(io) = index_global_to_local(pxdist,io,p2(io))
         if (myid == 0) then
            write(6,"(5i10)") io, p1(io), isrc(io), p2(io), idst(io)
         endif
      enddo

      do io=1,norbs
         if (myid == p2(io)) then
            call MPI_irecv(data2(idst(io)),1,MPI_integer,p1(io), &
                           io,mpi_comm_world,requestR(io),ierr)
         endif
      enddo
      do io=1,norbs
         if (myid == p1(io)) then
            call MPI_isend(data1(isrc(io)),1,MPI_integer,p2(io), &
                        io,mpi_comm_world,requestS(io),ierr)
         endif
      enddo
      do io=1,norbs
         if (myid == p2(io)) then
            call MPI_wait(requestR(io),status,ierr)
         endif
      enddo

      call MPI_Barrier(mpi_comm_world,ierr)

      if (worker) then
         do io = 1, no_l2
            print *, "Node: ", myid, "il, ig: ", io, data2(io)
         enddo
      endif

 30   call MPI_FINALIZE(rc)

    end program dist




