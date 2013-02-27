
      program dist

      use mpi_siesta
      use class_BlockCyclicDist
      use class_PEXSIDist

      implicit none

      integer n, myid, numprocs, i, rc, ntotal, ierr

      integer :: comm_sub, group_world, group_sub, myid_sub
      integer :: nsub, bs, pbs, norbs, io, p1, p2, isrc, idst

      logical :: worker
      integer, allocatable :: ranks(:)
      type(BlockCyclicDist) :: bcdist
      type(PEXSIDist)       :: pxdist


      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      if (myid  == 0) then
         print *, "Using ", numprocs, " procs in World."
!         print *, "Enter number of subset procs:"
!         read *, nsub
         ! Avoid complications for now
         nsub = numprocs
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
      call newDistribution(pbs,Comm_sub,pxdist,"px dist")
      
      if (myid == 0) then
         write(6,"(5a10)") "Orb", "p1", "i1", "p2", "i2"
      endif
      do io = 1, norbs
         p1 = node_handling_element(bcdist,io)
         p2 = node_handling_element(pxdist,io)
         isrc = index_global_to_local(bcdist,io,p1)
         idst = index_global_to_local(pxdist,io,p2)
         if (myid == 0) then
            write(6,"(5i10)") io, p1, isrc, p2, idst
         endif
      enddo

 30   call MPI_FINALIZE(rc)

    end program dist




