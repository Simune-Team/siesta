      subroutine initparallel(no_u,na_u,lasto,xa,ucell,rmaxh,rcoor,
     .                        isolve)
C
C  Initialises parallel parameters
C
C  Julian Gale, NRI, Curtin Uni, March 2004
C
      use precision
      use fdf
      use parallel,     only : Node, Nodes, BlockSize, ProcessorY
      use parallelsubs, only : set_processorY, set_BlockSizeDefault
      use parallelsubs, only : set_processorYdefault
      use spatial
      use sys
      use alloc
#ifdef MPI
      use mpi_siesta, only: mpi_comm_world, mpi_integer
#ifdef NODAT
      use mpi_siesta, only: mpi_double_precision
#else
      use mpi_siesta, only: DAT_double
#endif
#endif
      implicit none
C
C  Passed variables
C
      integer       :: no_u
      integer       :: na_u
      integer       :: isolve
      integer       :: lasto(0:na_u)
      real(dp)      :: xa(3,na_u)
      real(dp)      :: ucell(3,3)
      real(dp)      :: rcoor
      real(dp)      :: rmaxh
C
C  Local variables
C
#ifdef MPI
      integer       :: blocksizedefault
      integer       :: procYdefault
      integer       :: procYval
      integer       :: MPIerror
#endif
      integer       :: ii
      integer       :: maxorb
      logical       :: lspatialok
      logical, save :: first = .true.
      real(dp)      :: acell
      real(dp)      :: bcell
      real(dp)      :: ccell
      real(dp)      :: cellmin
      real(dp)      :: alpha
      real(dp)      :: beta
      real(dp)      :: gamma
      real(dp)      :: degtorad
      real(dp)      :: rcmax
      real(dp)      :: rcmaxopt

C  Initialise on first call
      if (first) then
        nullify(natomsNode)
        nullify(natomsL2G)
        nullify(natomsG2L)
        nullify(ncellnodeptr)
        nullify(lbuffercell)
        nullify(nL2G)
        nullify(nG2L)
        nullify(nNode)
        nullify(nOrbPerNode)
        first = .false.
      endif

C Set spatial decomposition flag according to SCF method
      if (isolve.eq.1) then
        lspatial = .true.
      else
        lspatial = .false.
      endif

C Spatial decomposition flag
      if (Node.eq.0) then
        rspatial = fdf_physical('RcSpatial',0.0d0,'Bohr')
      endif
#ifdef MPI
#ifdef NODAT
      call MPI_Bcast(rspatial,1,MPI_double_precision,0,MPI_Comm_World,
     .  MPIerror)
#else
      call MPI_Bcast(rspatial,1,DAT_double,0,MPI_Comm_World,MPIerror)
#endif
#endif

C Set ProcessorY  
#ifdef MPI
      if (Node.eq.0) then
        call set_processorYdefault(Nodes,procYdefault)
        procYval = fdf_integer('processorY',procYdefault)
      endif
      call MPI_Bcast(procYval,1,MPI_integer,0,MPI_Comm_World,MPIerror)
      call set_processorY(procYval)
#else
      call set_processorY(1)
#endif

C Override spatial cut-off
      if (lspatial.and.rspatial.gt.0.0d0) then
        rcmax = rspatial
      else
        rcmax = max(rmaxh,rcoor)
        rcmaxopt = rcmax
      endif

C Find cell parameters
      degtorad = 4.0d0*atan(1.0d0)/180.0d0
      call uncell(ucell,acell,bcell,ccell,alpha,beta,gamma,degtorad)

C Check that cut-off is less than cell parameter
      cellmin = acell
      if (bcell.lt.cellmin) cellmin = bcell
      if (ccell.lt.cellmin) cellmin = ccell
      if (cellmin.lt.2.0d0*rcmax) rcmax = 0.5d0*cellmin/dble(Nodes)

      if (.not.lspatial) then
C Set BlockSize
#ifdef MPI
        if (Node.eq.0) then
          call set_blocksizedefault(Nodes,no_u,blocksizedefault)
          BlockSize = fdf_integer('blocksize',blocksizedefault)
        endif
        call MPI_Bcast(BlockSize,1,MPI_integer,0,MPI_Comm_World,
     .    MPIerror)
#else
        call set_processorY(1)
        BlockSize = 8
#endif
      endif

#ifdef MPI
C Output indication of parallel parameters
      if (Node.eq.0) then
        if (lspatial) then
          write(6,'(/,a,f8.4,/)')
     .        '* Spatial decomposition: Cutoff = ',rcmax
        else
          write(6,'(/,a,2i4,/)')
     .        '* ProcessorY, Blocksize: ', procYval, Blocksize
        endif
      endif
#endif

C If spatial perform allocation of atoms to processors
      if (lspatial) then
        call setspatial(na_u,xa,ucell,rcmax,lspatialok)
        if (.not.lspatialok) then
          call die('Spatial decomposition failed')
        endif
        call setatomnodes(na_u,lasto,Node,Nodes)
        if (Nodes.gt.1) then
          ProcessorY = npgrid(2)
        else
          ProcessorY = 1
        endif
      else
C
C  Build dummy lists linking orbitals with parallel structure
C
        maxorb = lasto(na_u)
        call re_alloc(nL2G,1,maxorb,1,Nodes,name='nL2G')
        call re_alloc(nG2L,1,maxorb,name='nG2L')
        call re_alloc(nNode,1,maxorb,name='nNode')
        call re_alloc(nOrbPerNode,1,Nodes,name='nOrbPerNode')
C
        nOrbPerNode(1) = maxorb
        do ii = 1,maxorb
          nL2G(ii,1) = ii
          nG2L(ii) = ii
          nNode(ii) = 0
        enddo
      endif

      return
      end
