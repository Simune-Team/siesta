! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!


C ==================================================================
C Computes a new mesh distribution in order to balance the load
C among the processes. The load depends on numphi^2.
C ==================================================================
C SUBROUTINE recomputeLMesh( nm, nmpl, nml, numphi )
C
C INPUT:
C integer nm(3)         : Number of Mesh divisions of each cell vector
C integer nmpl          : Local number of mesh points
C integer nml(3)        : Local number of Mesh div. of each cell vector
C
C INOUT:
C integer numphi(nmpl)  : Number of orbitals that intersects every mesh
C                         point. Its value change inside but it's not
C                         relevant for the calling routine.
C
C OUTPUT:
C The output values are in the module MeshComm (New limits of the mesh
C and the information needed to transfer data between data
C distributions)
C
C BEHAVIOR:
C The amount of computation in every point of the mesh depends on 
C numphi^2. We want to distribute numphi among a 2D grid of processors,
C of dimension processorY * ProcessorZ.
C
C First, we reduce the 3D vector numphi to a 2D vector numphiYZ.
C
C We use recomputeLimZ to compute the limits in dimension Z of numphiYZ.
C We split numphiYZ in ProcessorZ partitions using Nested disection. Then
C we compute the communications needed to move the data to the new
C distribution. 
C
C For every partition in dimension Z, we should redistribute the data in
C dimension Y among ProcessorY partitions using the Nested disection
C algorithm. This gives the limits of the new data distribution. Finally
C we compute the communications needed to move data from the uniform data
C distribution to the "Nested disection" data distribution.
C
C ==================================================================
      subroutine recomputeLMesh( nm, nmpl, nml, numphi )
C     Modules
      use precision, only: i8b
      use MeshComm
      use parallel
      use debugMpi
      implicit none
C     Input variables
      integer, intent(in)     :: nm(3), nmpl, nml(3)
      integer, intent(inout)  :: numphi(nmpl)
C     Local variables
      integer(i8b)               :: tt
      integer                 :: ii, i1, i2, i3, ProcessorZ, PY, PZ
      integer(i8b),  allocatable :: numphiYZ(:,:)
      integer,    allocatable :: newLimz(:), oriLimz(:)
      integer,    allocatable :: newLimy(:), oriLimy(:)

!------------------------------------------------------------------------- BEGIN
      do ii= 1, nmpl
        numphi(ii) = numphi(ii)*numphi(ii)
C        if (numphi(ii).gt.1) then
C          numphi(ii) = 1
C        endif
      enddo

C     Number of processors in dimension Z
      ProcessorZ = Nodes/ProcessorY

C     Compute Process ID in dimension Y and Z
      Py    = (Node/ProcessorZ) + 1
      Pz    = Node - (Py - 1)*ProcessorZ + 1

C     Compute the limits in dimensions Z and Y
      allocate(oriLimz(ProcessorZ+1))
      call memory('A','I',ProcessorZ+1,'recomputeLMesh')
      allocate(newLimz(ProcessorZ+1))
      call memory('A','I',ProcessorZ+1,'recomputeLMesh')
      allocate(oriLimy(ProcessorY+1))
      call memory('A','I',ProcessorY+1,'recomputeLMesh')
      allocate(newLimy(ProcessorY+1))
      call memory('A','I',ProcessorY+1,'recomputeLMesh')

      call computeLimits( ProcessorZ, oriLimZ, nm(3) )
      call computeLimits( ProcessorY, oriLimY, nm(2) )

C     Compute the acumulation of numphi in dimension YZ
      allocate(numphiYZ(nml(2),nml(3)))
      call memory('A','E',nml(2)*nml(3),'recomputeLMesh')
      ii    = 1
      do i1=1, nml(3)
        do i2=1, nml(2)
          tt=0
          do i3=1, nml(1)
            tt = tt + numphi(ii)
            ii = ii + 1
          enddo
          numphiYZ(i2,i1) = tt
        enddo
      enddo

      call recomputeLimZ( nml(2), nml(3), numphiYZ, nm(3), PY, PZ,
     &                    ProcessorZ, oriLimZ, newLimZ )

      call computeCommZ( PZ, ProcessorZ, oriLimZ, newLimZ )

      call recomputeLimY( nml(2), nml(3), numphiYZ, nm(2), PY, PZ,
     &                    ProcessorZ, oriLimZ, oriLimY, newLimY )

      call computeCommY( PZ, ProcessorZ, PY, ProcessorY,
     &                   oriLimZ, newLimZ, oriLimY, newLimY )

      udMeshLim(1,1) = 1
      udMeshLim(2,1) = nml(1)
      udMeshLim(1,2) = oriLimY(PY)
      udMeshLim(2,2) = oriLimY(PY+1)-1
      udMeshLim(1,3) = oriLimZ(PZ)
      udMeshLim(2,3) = oriLimZ(PZ+1)-1

      nddMeshLim(1,1) = 1
      nddMeshLim(2,1) = nml(1)
      nddMeshLim(1,2) = newLimY(PY)
      nddMeshLim(2,2) = newLimY(PY+1)-1
      nddMeshLim(1,3) = newLimZ(PZ)
      nddMeshLim(2,3) = newLimZ(PZ+1)-1

      deallocate(numphiYZ)
      call memory('D','E',nml(2)*nml(3),'recomputeLMesh')
      deallocate(newLimY)
      call memory('D','I',ProcessorY+1,'recomputeLMesh')
      deallocate(oriLimY)
      call memory('D','I',ProcessorY+1,'recomputeLMesh')
      deallocate(newLimz)
      call memory('D','I',ProcessorZ+1,'recomputeLMesh')
      deallocate(oriLimz)
      call memory('D','I',ProcessorZ+1,'recomputeLMesh')
!--------------------------------------------------------------------------- END
      end subroutine recomputeLMesh

C ==================================================================
C Computes the limits of one of the dimensions of the mesh.
C The mesh is uniformly distributed among nproc processes.
C ==================================================================
C SUBROUTINE computeLimits( nproc, limits, length )
C
C INPUT:
C integer nproc           : Number of processes in the current dimension.
C integer length          : Lenght of the grid in the current dimension.
C
C OUTPUT:
C integer limits(nproc+1) : limits of the grid between the several processes.
C
C BEHAVIOR:
C Given one of the dimensions of the grid, to be splitted among nproc
C parts, we compute the limits of everyone (Uniform distribution).
C
C ==================================================================
      subroutine computeLimits( nproc, limits, length )
      implicit none
C     Input variables
      integer, intent(in)  :: nproc, length
C     Output variables
      integer, intent(out) :: limits(nproc+1)
C     Local variables
      integer              :: BSize, Nrem, pid, newlim
!------------------------------------------------------------------------- BEGIN
      BSize = length/nproc
      Nrem  = length - BSize*nproc

      limits(1) = 1
      do pid= 1, nproc
        newlim = limits(pid) + BSize
        if (pid.le.NRem) newlim = newlim + 1
        limits(pid+1) = newlim
      enddo
!--------------------------------------------------------------------------- END
      end subroutine computeLimits

C ==================================================================
C Distributes numphiYZ in ProcessorZ parts and computes the new
C limits in dimension Z.
C ==================================================================
C SUBROUTINE recomputeLimZ( nmlY, nmlZ, numphiYZ, nmZ, PY, PZ,
C     &                     ProcessorZ, oriLimZ, newLimZ )
C
C INPUT:
C integer nmlY                  : Local number of points in dimension Y.
C integer nmlY                  : Local number of points in dimension Z.
C integer numphiYZ(nmlY,nmlZ)   : Load of the mesh reduced to a 2D array.
C integer nmZ                   : Global number of points in dimension Z.
C integer PY                    : process ID in dimension Y.
C integer PZ                    : process ID in dimension Z.
C integer ProcessorZ            : Number of processors in dimension Z.
C integer oriLimZ(ProcessorZ+1) : Original limits of the grid in dimension Z.
C
C OUTPUT:
C integer newLimZ(ProcessorZ+1) : New limits of the grid in dimension Z.
C
C BEHAVIOR:
C Every process reduces the 2D array numphiYZ to its dimension Z, in the
C vector numphiZ. All processes with PY=1 acumulates its numphiZ, with the
C rest of processes that have the same PZ (same column of processors).
C
C Process 0 (PY=1,PZ=1), concatenates the several numphiZ that recieves from
C the processors of the first row (PY=1) to form GnumphiZ.
C
C This process usses nested disection to split GnumphiZ among ProcessorZ
C parts. Finally, it communicates this information to all the processes.
C
C ==================================================================
      subroutine recomputeLimZ( nmlY, nmlZ, numphiYZ, nmZ, PY, PZ,
     &                          ProcessorZ, oriLimZ, newLimZ )
#ifdef MPI
      use precision, only: i8b
      use parallel
      use mpi_siesta
      implicit none
C     Input variables
      integer                :: nmlY, nmlZ, nmZ, PY, PZ, ProcessorZ,
     &                          oriLimZ(ProcessorZ+1),
     &                          newLimZ(ProcessorZ+1)
      integer(i8b)              :: numphiYZ(nmlY,nmlZ)

C     Local variables
      integer(i8b)              :: tt, nlist
      integer                :: ii, jj, P1, P2, PP, blkZ
      integer(i8b), allocatable :: numphiZ(:), sereBUF(:), GnumphiZ(:)
      integer                :: MPIerror, Status(MPI_Status_Size)
!------------------------------------------------------------------------- BEGIN
C     Compute the acumulation of numphiYZ in dimension Z
      allocate(numphiZ(nmlZ))
      call memory('A','E',nmlZ,'recomputeLimZ')
      numphiZ(1:nmlZ) = 0
      do jj=1, nmlZ
        tt=0
        do ii=1, nmlY
          tt = tt + numphiYZ(ii,jj)
        enddo
        numphiZ(jj) = tt
      enddo

C     Gather numphiZ from all processors in two steps:
C       1) All Processes with PY=1 (first row of processors) recieve
C          from those with the same PZ as them (same column).
C       2) Process 0, recieve from all processes with PY=1
      if (PY.eq.1) then
        allocate(sereBUF(nmlZ))
        call memory('A','E',nmlZ,'recomputeLimZ')
        PP = PZ - 1
        do P1= 2, ProcessorY
          PP = PP + ProcessorZ
          call mpi_recv( sereBUF, nmlZ, MPI_INTEGER8, PP,
     &                   1, MPI_Comm_world, Status, MPIerror )
          numphiZ = numphiZ + sereBUF
        enddo
        deallocate(sereBUF)
        call memory('D','E',nmlZ,'recomputeLimZ')

        if (PZ.eq.1) then
          allocate(GnumphiZ(nmZ))
          call memory('A','E',nmZ,'recomputeLimZ')
C         Copy the part that is in the current Process (Node=0)
          GnumphiZ(1:nmlZ) = numphiZ(1:nmlZ)

C         Recieve from those processes that are in the same row (PY=1)
          do P1= 2, ProcessorZ
            PP    = P1 - 1
            ii    = oriLimZ(P1)
            blkZ  = oriLimZ(P1+1) - oriLimZ(P1)
            call mpi_recv( GnumphiZ(ii), blkZ, MPI_INTEGER8, PP,
     &                     1, MPI_Comm_world, Status, MPIerror )
          enddo

          nlist = 0
          do ii= 1, nmZ
            nlist = nlist + GnumphiZ(ii)
          enddo

C         Recompute the limits in dimension Z
          newLimz(1)            = 1
          newLimz(ProcessorZ+1) = nmZ + 1
          call vecBisection( nmZ, GnumphiZ, ProcessorZ,
     &                       newLimz, nlist )

          deallocate(GnumphiZ)
          call memory('D','E',nmZ,'recomputeLimZ')
        else
          PP = 0
          call MPI_Send( numphiZ, nmlZ, MPI_INTEGER8, PP, 1,
     &                   MPI_Comm_World, MPIerror )
        endif
      else
        PP = PZ - 1
        call MPI_Send( numphiZ, nmlZ, MPI_INTEGER8, PP, 1,
     &                 MPI_Comm_World, MPIerror )
      endif

      deallocate(numphiZ)
      call memory('D','E',nmlZ,'recomputeLimZ')

      call MPI_Bcast( newLimz, ProcessorZ+1, MPI_integer, 0,
     &                MPI_Comm_World, MPIerror )
#endif
!--------------------------------------------------------------------------- END
      end subroutine recomputeLimZ

C ==================================================================
C numphiYZ was distributed in ProcessorZ parts. Each part is now
C distributed in ProcessorY parts. We compute the new limits in
C dimension Y. There are ProcessorZ different newLimY vectors.
C ==================================================================
C SUBROUTINE recomputeLimY( nmlY, nmlZ, numphiYZ, nmY, PY, PZ,
C     &                     ProcessorZ, oriLimZ, oriLimY,
C     &                     newLimY )
C
C INPUT:
C integer nmlY                  : Local number of points in dimension Y.
C integer nmlY                  : Local number of points in dimension Z.
C integer numphiYZ(nmlY,nmlZ)   : Load of the mesh reduced to a 2D array.
C integer nmY                   : Global number of points in dimension Y.
C integer PY                    : process ID in dimension Y.
C integer PZ                    : process ID in dimension Z.
C integer ProcessorZ            : Number of processors in dimension Z.
C integer oriLimZ(ProcessorZ+1) : Original limits of the grid in dimension Z.
C integer oriLimY(ProcessorY+1) : Original limits of the grid in dimension Y.
C
C OUTPUT:
C integer newLimY(ProcessorY+1) : New limits of the grid in dimension Y.
C
C BEHAVIOR:
C We should reduce the 2D array numphiYZ to its dimension Y, according to
C the new limits on dimension Z. Every column of processors should create
C a different GnumphiY.
C
C First, we use the precomputed communication struct to send, recieve or
C keep numphiY. Then all the processes from the first row concatenates
C the several numphiY that recieves from the rest of its column to form
C GnumphiY.
C
C Finally, the processes of the first row calls vecBisection to split
C GnumphiY among ProcessorY parts (limits in dimension Y) and send this
C information to the processes of its column.
C
C ==================================================================
      subroutine recomputeLimY( nmlY, nmlZ, numphiYZ, nmY, PY, PZ,
     &                          ProcessorZ, oriLimZ, oriLimY,
     &                          newLimY )
#ifdef MPI
      use precision, only: i8b
      use parallel
      use MeshComm
      use mpi_siesta
      implicit none
C     Input variables
      integer                :: nmlY, nmlZ, nmY, PY, PZ, ProcessorZ,
     &                          oriLimZ(ProcessorZ+1),
     &                          oriLimY(ProcessorY+1),
     &                          newLimY(ProcessorY+1)
      integer(i8b)              :: numphiYZ(nmlY,nmlZ)
C     Local variables
      integer                :: ii, i1, i2, P1, P2, PP, blkY, com
      integer(i8b)              :: nlist
      integer(i8b), allocatable :: numphiY(:), sereBUF(:), GnumphiY(:)
      integer                :: MPIerror, Status(MPI_Status_Size)
!------------------------------------------------------------------------- BEGIN

C     Compute the acumulation of numphi in dimension Y,
C     with the new ordering in dimension Z.
      allocate(numphiY(nmlY))
      call memory('A','E',nmlY,'recomputeLimY')
      allocate(sereBUF(nmlY))
      call memory('A','E',nmlY,'recomputeLimY')
      numphiY = 0

C     Loop for all the communications that we need to do
C     in dimension Z
      do com= 1, MCsize
        if (MCsrc(com).eq.Pz) then
          if (MCdst(com).eq.Pz) then
C           If source and destiny are the current process then
C           just accumulate numphiYZ in numphiY
            ii = MCinZ(com) - oriLimZ(PZ) + 1
            do i1= MCinZ(com), MCfiZ(com)
              do i2= 1, nmlY
                numphiY(i2) = numphiY(i2) + numphiYZ(i2,ii)
              enddo
              ii = ii + 1
            enddo
          else
C           Fill a buffer send it to the destiny node
            PP      = MCdst(com)-1 + (PY - 1)*ProcessorZ
            sereBUF = 0
            ii      = MCinZ(com) - oriLimZ(PZ) + 1
            do i1= MCinZ(com), MCfiZ(com)
              do i2= 1, nmlY
                sereBUF(i2) = sereBUF(i2) + numphiYZ(i2,ii)
              enddo
              ii = ii + 1
            enddo
            call MPI_Send( sereBUF, nmlY, MPI_INTEGER8, PP, 1,
     &                     MPI_Comm_World, MPIerror )
          endif
        else
C         Recive a buffer and accumulate it to numphiY
          PP = MCsrc(com)-1 + (PY - 1)*ProcessorZ
          call mpi_recv( sereBUF, nmlY, MPI_INTEGER8, PP,
     &                   1, MPI_Comm_world, Status, MPIerror )
          do i2= 1, nmlY
            numphiY(i2) = numphiY(i2) + sereBUF(i2)
          enddo
        endif
      enddo

      deallocate(MCBuffer)
      call memory('D','I',nmlY,'recomputeLimY')

      deallocate(sereBUF)
      call memory('D','E',nmlY,'recomputeLimY')

C     Gather numphiY from all processes to the firsts nodes of dimension Y
C     and recompute the limits in the dimension Y.
      if (PY.eq.1) then
        allocate(GnumphiY(nmY))
        call memory('A','E',nmY,'recomputeLMesh')

        GnumphiY(1:nmlY) = numphiY(1:nmlY)

        PP = Node
        do P1= 2, ProcessorY
          PP   = PP + ProcessorZ
          blkY = oriLimY(P1+1) - oriLimY(P1)
          ii   = oriLimY(P1)
          call mpi_recv( GnumphiY(ii), blkY, MPI_INTEGER8, PP,
     &                   1, MPI_Comm_world, Status, MPIerror )
        enddo

C       Compute the acumulation of numphi in dimension Y
        nlist = 0
        do ii= 1, nmY
          nlist = nlist + GnumphiY(ii)
        enddo

C       Recompute the limits in dimension Y
        newLimy(1)            = 1
        newLimy(ProcessorY+1) = nmY + 1
        call vecBisection( nmY, GnumphiY, ProcessorY,
     &                     newLimY, nlist )
        deallocate(GnumphiY)
        call memory('D','E',nmY,'recomputeLMesh')

        PP = Node
        do P1= 2, ProcessorY
          PP = PP + ProcessorZ
          call mpi_Send( newLimY, ProcessorY+1, MPI_integer, PP,
     &                   1, MPI_Comm_world, MPIerror )
        enddo
      else
        PP = PZ - 1
        call MPI_Send( numphiY, nmlY, MPI_INTEGER8, PP, 1,
     &                 MPI_Comm_World, MPIerror )
        call mpi_recv( newLimY, ProcessorY+1, MPI_integer, PP,
     &                 1, MPI_Comm_world, Status, MPIerror )
      endif

      deallocate(numphiY)
      call memory('D','E',nmlY,'recomputeLimY')
#endif
!--------------------------------------------------------------------------- END
      end subroutine recomputeLimY

C ==================================================================
C Given the new limits in dimension Z, computes the needed
C communications to converts data from the original distribution
C to the new one. We use scheduleComm to optimize the commmunication
C time and avoid death locks. The information is stored in the module 
C MeshComm.
C ==================================================================
C SUBROUTINE computeCommZ( PZ, ProcessorZ, oriLimZ, newLimZ )
C
C INPUT:
C integer PZ                    : process ID in dimension Z.
C integer ProcessorZ            : Number of processors in dimension Z.
C integer oriLimZ(ProcessorZ+1) : Original limits of the grid in dimension Z.
C integer newLimZ(ProcessorZ+1) : New limits of the grid in dimension Z.
C
C OUTPUT:
C This subroutine modifies values of module MeshComm.
C integer MCsize        : Number of communications.
C integer MCsrc(MCsize) : Source process.
C integer MCdst(MCsize) : Destiny process.
C integer MCinZ(MCsize) : Initial limit in dimension Z.
C integer MCfiZ(MCsize) : Final limit in dimension Z.
C
C BEHAVIOR:
C Every intersection between the original data distribution an the new
C one represents a new communication. First we count the number of
C communications. Then we make an ordered list of communications and we
C call scheduleComm to optimize the commmunication time and avoid death
C locks. Finally, every process stracts its own communications and
C it saves it in the module MeshComm.
C
C ==================================================================
      subroutine computeCommZ( PZ, ProcessorZ, oriLimZ, newLimZ )
      use MeshComm
      use scheComm
      implicit none
C     Input variables
      integer              :: PZ, ProcessorZ,
     &                        oriLimZ(ProcessorZ+1),
     &                        newLimZ(ProcessorZ+1)
C     Local variables
      integer              :: P1, P2, Ncom, com, pos, col
      integer, allocatable :: src(:), dst(:), ini(:), fin(:)
      type(COMM_T)          :: commZ
!------------------------------------------------------------------------- BEGIN
C     Count the number of communications
      P1  = 1
      P2  = 1
      com = 0
      do while(P1.ne.ProcessorZ .or. P2.ne.ProcessorZ)
        com = com + 1
        if (oriLimZ(P1+1).lt.newLimZ(P2+1)) then
          P1 = P1 + 1
        else if (oriLimZ(P1+1).gt.newLimZ(P2+1)) then
          P2 = P2 + 1
        else
          P1 = P1 + 1
          P2 = P2 + 1
        endif
      enddo
      NCom = com + 1

      allocate(src(NCom))
      call memory('A','I',NCom,'computeCommZ')
      allocate(dst(NCom))
      call memory('A','I',NCom,'computeCommZ')
      allocate(ini(NCom))
      call memory('A','I',NCom,'computeCommZ')
      allocate(fin(NCom))
      call memory('A','I',NCom,'computeCommZ')

      P1  = 1
      P2  = 1
      com = 0
      pos = 1
C     Compute the communications
      do while(P1.ne.ProcessorZ .or. P2.ne.ProcessorZ)
        com = com + 1
        src(com) = P1
        dst(com) = P2
        ini(com) = pos
        fin(com) = min(oriLimZ(P1+1),newLimZ(P2+1))-1
        if (oriLimZ(P1+1).lt.newLimZ(P2+1)) then
          P1  = P1 + 1
          pos = oriLimZ(P1)
        else if (oriLimZ(P1+1).gt.newLimZ(P2+1)) then
          P2  = P2 + 1
          pos = newLimZ(P2)
        else
          P1  = P1 + 1
          P2  = P2 + 1
          pos = oriLimZ(P1)
        endif
      enddo
      src(Ncom) = P1
      dst(Ncom) = P2
      ini(Ncom) = pos
      fin(Ncom) = oriLimZ(P1+1)-1

      commZ%np = ProcessorZ
C     reschedule the communications in order to minimize the time
      call scheduleComm( Ncom, src, dst, commZ )

C     Extract only the communications of the current process.
C        MCsize = Number of comunications of the current process
      com = 0
      do col= 1, commZ%Ncol
        if (commZ%ind(col,PZ).ne.0) com = com + 1
      enddo
      MCsize = com

      allocate(MCBuffer(MCsize*4))
      call memory('A','I',MCsize*4,'recomputeLMesh')

      MCsrc => MCBuffer(         1:  MCsize)
      MCdst => MCBuffer(  MCsize+1:2*MCsize)
      MCinZ => MCBuffer(2*MCsize+1:3*MCsize)
      MCfiZ => MCBuffer(3*MCsize+1:4*MCsize)
      com = 0
      do col= 1, commZ%Ncol
        pos = commZ%ind(col,PZ)
        if (pos.ne.0) then
          com = com + 1
          MCsrc(com) = src(pos)
          MCdst(com) = dst(pos)
          MCinZ(com) = ini(pos)
          MCfiZ(com) = fin(pos)
        endif
      enddo

      deallocate(commZ%ind)
      call memory('D','I',commZ%Ncol*ProcessorZ,'computeCommZ')
      deallocate(fin)
      call memory('D','I',NCom,'computeCommZ')
      deallocate(ini)
      call memory('D','I',NCom,'computeCommZ')
      deallocate(dst)
      call memory('D','I',NCom,'computeCommZ')
      deallocate(src)
      call memory('D','I',NCom,'computeCommZ')
!--------------------------------------------------------------------------- END
      end subroutine computeCommZ

C ==================================================================
C Given the new limits in dimension Y and Z, computes the needed
C communications to converts data from the original distribution
C to the new one. We use scheduleComm to optimize the commmunication
C time and avoid death locks. The information is stored in the module 
C MeshComm.
C ==================================================================
C SUBROUTINE computeCommY( PZ, ProcessorZ, PY, ProcessorY,
C     &                    oriLimZ, newLimZ, oriLimY, newLimY )
C
C INPUT:
C integer PZ                    : process ID in dimension Z.
C integer ProcessorZ            : Number of processors in dimension Z.
C integer PY                    : process ID in dimension Y.
C integer ProcessorY            : Number of processors in dimension Y.
C integer oriLimZ(ProcessorZ+1) : Original limits of the grid in dimension Z.
C integer newLimZ(ProcessorZ+1) : New limits of the grid in dimension Z.
C integer oriLimY(ProcessorY+1) : Original limits of the grid in dimension Y.
C integer newLimY(ProcessorY+1) : New limits of the grid in dimension Y.
C
C OUTPUT:
C This subroutine modifies values of module MeshComm.
C integer MCsize        : Number of communications.
C integer MCsrc(MCsize) : Source process.
C integer MCdst(MCsize) : Destiny process.
C integer MCinZ(MCsize) : Initial limit in dimension Z.
C integer MCfiZ(MCsize) : Final limit in dimension Z.
C integer MCinY(MCsize) : Initial limit in dimension Y.
C integer MCfiY(MCsize) : Final limit in dimension Y.
C
C BEHAVIOR:
C Every intersection between the original data distribution an the new
C one represents a new communication. First, the process 0 gather
C all the different newLimY from the rest of the processes. It counts
C the number of communications, make a list and calls scheduleComm to
C optimize the commmunication time and avoid death locks.
C
C Finally, every process recieve the information computed for the
C process 0 and saves it in the module MeshComm.
C
C ==================================================================
      subroutine computeCommY( PZ, ProcessorZ, PY, ProcessorY,
     &                         oriLimZ, newLimZ, oriLimY, newLimY )
#ifdef MPI
      use parallel, only: Node, Nodes
      use MeshComm
      use scheComm
      use mpi_siesta
      implicit none
C     Input variables
      integer              :: PZ, ProcessorZ, PY, ProcessorY,
     &                        oriLimZ(ProcessorZ+1),
     &                        newLimZ(ProcessorZ+1),
     &                        oriLimY(ProcessorY+1),
     &                        newLimY(ProcessorY+1)
C     Local variables
      integer, allocatable :: GnewLimY(:,:)
      integer              :: PP, PZ1, PZ2, PY1, PY2, Ncom, com,
     &                        poZ, poY, pos, col
      integer, allocatable :: src(:), dst(:), inZ(:), fiZ(:),
     &                        inY(:), fiY(:), sereBUF(:)
      type(COMM_T)         :: comm
      integer              :: MPIerror, Status(MPI_Status_Size)
!------------------------------------------------------------------------- BEGIN
      if (Node.eq.0) then
        allocate(GnewLimY(ProcessorY+1,ProcessorZ))
        call memory('A','I',(ProcessorY+1)*ProcessorZ,'computeCommY')

        GnewLimY(1:ProcessorY+1,1) = newLimY(1:ProcessorY+1)
        do PZ1= 2, ProcessorZ
          call mpi_recv( GnewLimY(1,PZ1), ProcessorY+1, MPI_integer,
     &                   PZ1-1, 1, MPI_Comm_world, Status, MPIerror )
        enddo

C       Count the number of communications
        PZ1 = 1
        PZ2 = 1
        com = 0
        do while(PZ1.ne.ProcessorZ .or. PZ2.ne.ProcessorZ)
          PY1 = 1
          PY2 = 1
          do while(PY1.ne.ProcessorY .or. PY2.ne.ProcessorY)
            com = com + 1
            if (oriLimY(PY1+1).lt.GnewLimY(PY2+1,PZ2)) then
              PY1 = PY1 + 1
            else if (oriLimY(PY1+1).gt.GnewLimY(PY2+1,PZ2)) then
              PY2 = PY2 + 1
            else
              PY1 = PY1 + 1
              PY2 = PY2 + 1
            endif
          enddo
          com = com + 1

          if (oriLimZ(PZ1+1).lt.newLimZ(PZ2+1)) then
            PZ1 = PZ1 + 1
          else if (oriLimZ(PZ1+1).gt.newLimZ(PZ2+1)) then
            PZ2 = PZ2 + 1
          else
            PZ1 = PZ1 + 1
            PZ2 = PZ2 + 1
          endif
        enddo
C       Last column
        PY1 = 1
        PY2 = 1
        do while(PY1.ne.ProcessorY .or. PY2.ne.ProcessorY)
          com = com + 1
          if (oriLimY(PY1+1).lt.GnewLimY(PY2+1,PZ2)) then
            PY1 = PY1 + 1
          else if (oriLimY(PY1+1).gt.GnewLimY(PY2+1,PZ2)) then
            PY2 = PY2 + 1
          else
            PY1 = PY1 + 1
            PY2 = PY2 + 1
          endif
        enddo
        NCom = com + 1

        allocate(src(NCom))
        call memory('A','I',NCom,'computeCommZ')
        allocate(dst(NCom))
        call memory('A','I',NCom,'computeCommZ')
        allocate(inZ(NCom))
        call memory('A','I',NCom,'computeCommZ')
        allocate(fiZ(NCom))
        call memory('A','I',NCom,'computeCommZ')
        allocate(inY(NCom))
        call memory('A','I',NCom,'computeCommZ')
        allocate(fiY(NCom))
        call memory('A','I',NCom,'computeCommZ')

C       Compute the communications
        PZ1 = 1
        PZ2 = 1
        com = 0
        poZ = 1
        do while(PZ1.ne.ProcessorZ .or. PZ2.ne.ProcessorZ)
          PY1 = 1
          PY2 = 1
          poY = 1
          do while(PY1.ne.ProcessorY .or. PY2.ne.ProcessorY)
            com      = com + 1
            src(com) = PZ1 + (PY1-1)*ProcessorZ
            dst(com) = PZ2 + (PY2-1)*ProcessorZ
            inZ(com) = poZ
            inY(com) = poY
            fiZ(com) = min(oriLimZ(PZ1+1),newLimZ(PZ2+1))-1
            fiY(com) = min(oriLimY(PY1+1),GnewLimY(PY2+1,PZ2))-1
            if (oriLimY(PY1+1).lt.GnewLimY(PY2+1,PZ2)) then
              PY1 = PY1 + 1
              poY = oriLimY(PY1)
            else if (oriLimY(PY1+1).gt.GnewLimY(PY2+1,PZ2)) then
              PY2 = PY2 + 1
              poY = GnewLimY(PY2,PZ2)
            else
              PY1 = PY1 + 1
              PY2 = PY2 + 1
              poY = oriLimY(PY1)
            endif
          enddo
          com = com + 1
          src(com) = PZ1 + (PY1-1)*ProcessorZ
          dst(com) = PZ2 + (PY2-1)*ProcessorZ
          inZ(com) = poZ
          inY(com) = poY
          fiZ(com) = min(oriLimZ(PZ1+1),newLimZ(PZ2+1))-1
          fiY(com) = min(oriLimY(PY1+1),GnewLimY(PY2+1,PZ2))-1

          if (oriLimZ(PZ1+1).lt.newLimZ(PZ2+1)) then
            PZ1 = PZ1 + 1
            poZ = oriLimZ(PZ1)
          else if (oriLimZ(PZ1+1).gt.newLimZ(PZ2+1)) then
            PZ2 = PZ2 + 1
            poZ = newLimZ(PZ2)
          else
            PZ1 = PZ1 + 1
            PZ2 = PZ2 + 1
            poZ = oriLimZ(PZ1)
          endif
        enddo
C       Last column
        PY1 = 1
        PY2 = 1
        poY = 1
        do while(PY1.ne.ProcessorY .or. PY2.ne.ProcessorY)
          com = com + 1
          src(com) = PZ1 + (PY1-1)*ProcessorZ
          dst(com) = PZ2 + (PY2-1)*ProcessorZ
          inZ(com) = poZ
          inY(com) = poY
          fiZ(com) = min(oriLimZ(PZ1+1),newLimZ(PZ2+1))-1
          fiY(com) = min(oriLimY(PY1+1),GnewLimY(PY2+1,PZ2))-1
          if (oriLimY(PY1+1).lt.GnewLimY(PY2+1,PZ2)) then
            PY1 = PY1 + 1
            poY = oriLimY(PY1)
          else if (oriLimY(PY1+1).gt.GnewLimY(PY2+1,PZ2)) then
            PY2 = PY2 + 1
            poY = GnewLimY(PY2,PZ2)
          else
            PY1 = PY1 + 1
            PY2 = PY2 + 1
            poY = oriLimY(PY1)
          endif
        enddo
        com = com + 1
        src(com) = PZ1 + (PY1-1)*ProcessorZ
        dst(com) = PZ2 + (PY2-1)*ProcessorZ
        inZ(com) = poZ
        inY(com) = poY
        fiZ(com) = min(oriLimZ(PZ1+1),newLimZ(PZ2+1))-1
        fiY(com) = min(oriLimY(PY1+1),GnewLimY(PY2+1,PZ2))-1

        deallocate(GnewLimY)
        call memory('D','I',(ProcessorY+1)*ProcessorZ,'computeCommY')

        comm%np = Nodes
C       reschedule the communications in order to minimize the time
        call scheduleComm( Ncom, src, dst, comm )

        allocate(sereBUF(comm%ncol*6))
        call memory('A','I',comm%ncol*6,'computeCommY')

        do PP= 2, Nodes
          com = 0
          do col= 1, comm%Ncol
           if (comm%ind(col,PP).ne.0) com = com + 1
          enddo
          Ncom = com

          com = 0
          do col= 1, comm%Ncol
            pos = comm%ind(col,PP)
            if (pos.ne.0) then
              com = com + 1
              sereBUF(com       ) = src(pos)
              sereBUF(com+  Ncom) = dst(pos)
              sereBUF(com+2*Ncom) = inZ(pos)
              sereBUF(com+3*Ncom) = fiZ(pos)
              sereBUF(com+4*Ncom) = inY(pos)
              sereBUF(com+5*Ncom) = fiY(pos)
            endif
          enddo
          call MPI_Send( Ncom, 1, MPI_integer,
     &                   PP-1, 1, MPI_Comm_World, MPIerror )
          call MPI_Send( sereBUF, Ncom*6, MPI_integer,
     &                   PP-1, 1, MPI_Comm_World, MPIerror )

        enddo

        deallocate(sereBUF)
        call memory('D','I',comm%ncol*6,'computeCommY')

        PP  = 1
        com = 0
        do col= 1, comm%Ncol
         if (comm%ind(col,PP).ne.0) com = com + 1
        enddo
        MCsize = com

        allocate(MCBuffer(MCsize*6))
        call memory('A','I',MCsize*6,'computeCommY')

        MCsrc => MCBuffer(         1:  MCsize)
        MCdst => MCBuffer(  MCsize+1:2*MCsize)
        MCinZ => MCBuffer(2*MCsize+1:3*MCsize)
        MCfiZ => MCBuffer(3*MCsize+1:4*MCsize)
        MCinY => MCBuffer(4*MCsize+1:5*MCsize)
        MCfiY => MCBuffer(5*MCsize+1:6*MCsize)

        com = 0
        do col= 1, comm%Ncol
          pos = comm%ind(col,PP)
          if (pos.ne.0) then
            com = com + 1
            MCsrc(com) = src(pos)
            MCdst(com) = dst(pos)
            MCinZ(com) = inZ(pos)
            MCfiZ(com) = fiZ(pos)
            MCinY(com) = inY(pos)
            MCfiY(com) = fiY(pos)
          endif
        enddo

        deallocate(src)
        call memory('D','I',NCom,'computeCommZ')
        deallocate(dst)
        call memory('D','I',NCom,'computeCommZ')
        deallocate(inZ)
        call memory('D','I',NCom,'computeCommZ')
        deallocate(fiZ)
        call memory('D','I',NCom,'computeCommZ')
        deallocate(inY)
        call memory('D','I',NCom,'computeCommZ')
        deallocate(fiY)
        call memory('D','I',NCom,'computeCommZ')

      else if (PY.eq.1) then
        call MPI_Send( newLimY, ProcessorY+1, MPI_integer,
     &                 0, 1, MPI_Comm_World, MPIerror )
      endif

      if (Node.ne.0) then
        call mpi_recv( MCsize, 1, MPI_integer, 0,  1,
     &                 MPI_Comm_world, Status, MPIerror )

        allocate(MCBuffer(MCsize*6))
        call memory('A','I',MCsize*6,'computeCommY')
 
        call mpi_recv( MCBuffer, MCsize*6, MPI_integer, 0,  1,
     &                 MPI_Comm_world, Status, MPIerror )
        MCsrc => MCBuffer(         1:  MCsize)
        MCdst => MCBuffer(  MCsize+1:2*MCsize)
        MCinZ => MCBuffer(2*MCsize+1:3*MCsize)
        MCfiZ => MCBuffer(3*MCsize+1:4*MCsize)
        MCinY => MCBuffer(4*MCsize+1:5*MCsize)
        MCfiY => MCBuffer(5*MCsize+1:6*MCsize)
      endif
      MCsrc(1:MCsize) = MCsrc(1:MCsize)-1
      MCdst(1:MCsize) = MCdst(1:MCsize)-1
#endif
!--------------------------------------------------------------------------- END
      end subroutine computeCommY

C ==================================================================
C Try to split the vector 'values' in 'nparts', trying to balance
C the values of the vector among the several parts. We use the
C nested disection algorithm to split the vector.
C ==================================================================
C subroutine vecBisection( nval, values, nparts, parts, total )
C
C INPUT:
C integer nval            : vector size.
C integer values          : vector to be splitted
C integer nparts          : number of parts that we want to calculate
C integer total           : Total load of vector 'values'
C
C OUTPUT:
C integer parts(nparts+1) : Limits of the partitions.
C
C BEHAVIOR:
C This is a recursive subroutine. Every time that we call vecBisection
C we try to split the vector in two parts with the same load. If nparts is
C an odd number the loads should be ponderated. For every half, we call
C vecBisection if the remaining number of parts is greater than 1.
C
C ==================================================================
      recursive subroutine vecBisection( nval, values,
     &                                   nparts, parts, total )
      use precision, only: i8b
      implicit none
C     Input variables
      integer               :: nval,  nparts, parts(nparts+1)
      integer(i8b)             :: values(nval), total

C     Local variables
      integer               :: half1, half2, ii
      integer(i8b)             :: LoadL, LoadT
!------------------------------------------------------------------------- BEGIN
      if (nparts.gt.1) then
        half1 = nparts/2
        half2 = nparts - half1
        LoadL = (total*half1)/nparts
        LoadT = 0
        ii    = 0
        do while(LoadT.lt.LoadL)
          ii    = ii + 1
          if (ii.eq.nval) STOP 'ERROR in vecBisection'
          LoadT = LoadT + values(ii)
        enddo
        if ((LoadT-values(ii)/2).gt.LoadL) then
          LoadT = LoadT - values(ii)
          ii    = ii - 1
        endif
        parts(half1+1) = parts(1) + ii
        if (half1.gt.1) then
          call vecBisection( ii, values, half1, parts, LoadT )
        endif
        if (half2.gt.1) then
          call vecBisection( nval-ii, values(ii+1), half2,
     &                       parts(half1+1), total-LoadT )
        endif
      endif
!--------------------------------------------------------------------------- END
      end subroutine vecBisection
