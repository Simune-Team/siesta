      subroutine iorho( task, fname, cell, mesh, nsm, maxp, nspin, rho,
     .                  found )

C *********************************************************************
C Saves the electron density at the mesh points.
C Writen by J.Soler July 1997.
C Parallel modifications added, while maintaining independence
C of density matrix on disk from parallel distribution. Uses a
C block distribution for density matrix. It is important to
C remember that the density matrix is divided so that all
C sub-points reside on the same node. Modified by J.D.Gale March 1999.
C NOTE : in order to achieve a consistent format of the disk file
C each record in the unformatted file corresponds to one pair of
C Y and Z values. Hence there will be a total of mesh(2) x mesh(3)
C records.
C *************************** INPUT **********************************
C character*(*) task      : 'read'/'READ' or 'write'/'WRITE'
C character*(*) fname     : File name for input or output
C integer nsm             : Number of sub-mesh points per mesh point
C integer maxp            : First dimension of array rho
C integer nspin           : Number of spin polarizations (1 or 2)
C ******************** INPUT or OUTPUT (depending on task) ***********
C real*8  cell(3,3)       : Lattice vectors
C integer mesh(3)         : Number of mesh divisions of each
C                           lattice vector
C real    rho(maxp,nspin) : Electron density
C                           Notice single precision in this version
C ******************** OUTPUT *****************************************
C logical found           : Were data found? (only when task='read')
C *************************** UNITS ***********************************
C Units should be consistent between task='read' and 'write'
C *********************************************************************

C
C  Modules
C
      use precision
      use parallel
      use fdf
#ifdef MPI
      use mpi
#endif

      implicit          none

      character*(*)     fname, task
      integer           maxp, mesh(3), nspin, nsm
      real              rho(maxp,nspin)
      double precision  cell(3,3)
      external          io_assign, io_close, paste, memory

c Internal variables and arrays
      character  paste*33, ffor*9, fform*(*)
      integer    i, ip, iu, is, j, ju, np, ns, Node, 
     .           npmax, BlockSizeY, RemainderY, ProcessorZ,
     .           BlockSizeZ, npl,
     .           NZ, iy, iz, izm, Ind, ir

#ifdef MPI
      integer    MPIerror, meshnsm(3), Request, meshl(3),
     .           Status(MPI_Status_Size), Nodes, BNode, NBlock
      logical    ltmp
      real, dimension(:), allocatable, save :: bdens
#endif

      logical    baddim, found

c Fix whether formatted or unformatted files wil be used
      parameter ( fform = 'unformatted' )

C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
#else
      Node = 0
#endif

#ifdef MPI
C Work out density block dimensions
      if (mod(Nodes,ProcessorY).gt.0) then
        write(6,'(''ERROR: ProcessorY must be a factor of the'',
     .    '' number of processors!'')')
        stop
      endif
      ProcessorZ = Nodes/ProcessorY
      BlockSizeY = ((((mesh(2)/nsm)-1)/ProcessorY) + 1)*nsm
      BlockSizeZ = ((((mesh(3)/nsm)-1)/ProcessorZ) + 1)*nsm
      allocate(bdens(BlockSizeY*mesh(3)))
      call memory('A','S',BlockSizeY*mesh(3),'iorho')
#else
      ProcessorZ = 1
      BlockSizeY = mesh(2)
      BlockSizeZ = mesh(3)
#endif

c Choose between read or write
      if (task.eq.'read' .or. task.eq.'READ') then

c       Check if input file exists
        if (Node.eq.0) then
          inquire( file=fname, exist=found )
        endif
#ifdef MPI
        call MPI_Bcast(found,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
        if (found) then

c         Open file
          if (Node.eq.0) then
            call io_assign( iu )
            open( iu, file=fname, form=fform, status='old' )      

c         Read cell vectors and number of points
            ffor = 'formatted'
            if (fform .eq. ffor) then
              read(iu) cell
              read(iu) mesh, ns
            else
              read(iu,*) cell
              read(iu,*) mesh, ns
            endif
          endif
#ifdef MPI
          call MPI_Bcast(cell(1,1),9,DAT_double,0,
     .      MPI_Comm_World,MPIerror)
          call MPI_Bcast(mesh,3,MPI_integer,0,MPI_Comm_World,MPIerror)
          call MPI_Bcast(ns,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif
          np = mesh(1) * mesh(2) * mesh(3)

#ifdef MPI
C  Get local dimensions
          meshnsm(1) = mesh(1)/nsm
          meshnsm(2) = mesh(2)/nsm
          meshnsm(3) = mesh(3)/nsm
          call HowManyMeshPerNode(meshnsm,Node,Nodes,npl,meshl)
#else
          npl = np
#endif

c  Check dimensions
          baddim = .false.
          if (ns .ne. nspin) baddim = .true.
          if (npl .gt. maxp) baddim = .true.

#ifdef MPI
C  Globalise dimension check
          call MPI_AllReduce(baddim,ltmp,1,MPI_logical,MPI_Lor,
     .      MPI_Comm_World,MPIerror)
          baddim = ltmp
#endif

          if (baddim) then
#ifdef MPI
C  Find largest value of npl
            call MPI_AllReduce(npl,npmax,1,MPI_integer,MPI_Max,
     .        MPI_Comm_World,MPIerror)
#else
            npmax = np
#endif
            if (Node.eq.0) then
              call io_assign( ju )
              open( ju, file='iorho.h', status='unknown' )
              write(ju,'(a)') 'C Dimensions for input to iorho'
              write(ju,'(6x,a,i8,a)') 'parameter ( nspin =', ns, ' )'
              write(ju,'(6x,a,i8,a)') 'parameter ( maxp  =', npmax,')'
              call io_close( ju )
              write(6,'(a)') 'iorho: ERROR: BAD DIMENSIONS'
              stop           'iorho: ERROR: BAD DIMENSIONS'
            else
              stop
            endif
          endif

C  Outer loop over spins
          do is = 1,ns

          Ind = 0

C  Loop over Z mesh direction
          do iz = 1,ProcessorZ

            if (iz.eq.ProcessorZ) then
              NZ = mesh(3) - (ProcessorZ-1)*BlockSizeZ
            else
              NZ = BlockSizeZ
            endif

C  Loop over local Z mesh points
            do izm = 1,NZ

C  Loop over blocks in Y mesh direction
              do iy = 1,ProcessorY

C  Work out size of density sub-matrix to be transfered
                if (iy.eq.ProcessorY) then
                  RemainderY = mesh(2) - (ProcessorY-1)*BlockSizeY
                else
                  RemainderY = BlockSizeY
                endif

#ifdef MPI
                NBlock = RemainderY*mesh(3)
C  Work out which node block is stored on
                BNode = (iy-1)*ProcessorZ + iz - 1

                if (BNode.eq.0.and.Node.eq.BNode) then
#endif
C  If density sub-matrix is local Node 0 then just read it in
                  if (fform .eq. ffor) then
                    do ir = 1,RemainderY
                      read(iu,*) (rho(Ind+ip,is),ip=1,mesh(3))
                      Ind = Ind + mesh(3)
                    enddo
                  else
                    do ir = 1,RemainderY
                      read(iu) (rho(Ind+ip,is),ip=1,mesh(3))
                      Ind = Ind + mesh(3)
                    enddo
                  endif

#ifdef MPI
                elseif (Node.eq.0) then
C  If this is Node 0 then read and send density sub-matrix
                  Ind = 0
                  if (fform .eq. ffor) then
                    do ir = 1,RemainderY
                      read(iu,*) (bdens(Ind+ip),ip=1,mesh(3))
                      Ind = Ind + mesh(3)
                    enddo
                  else
                    do ir = 1,RemainderY
                      read(iu) (bdens(Ind+ip),ip=1,mesh(3))
                      Ind = Ind + mesh(3)
                    enddo
                  endif
                  call MPI_ISend(bdens,NBlock,DAT_single,BNode,1,
     .              MPI_Comm_World,Request,MPIerror)
                  call MPI_Wait(Request,Status,MPIerror)

                elseif (Node.eq.BNode) then
C  If this is the Node where the density sub-matrix is, then receive
                  call MPI_IRecv(rho(Ind+1,is),NBlock,DAT_single,
     .              0,1,MPI_Comm_World,Request,MPIerror)
                  call MPI_Wait(Request,Status,MPIerror)
                  Ind = Ind + NBlock

                endif

                if (BNode.ne.0) then
                  call MPI_Barrier(MPI_Comm_World,MPIerror)
                endif
#endif

              enddo

            enddo

          enddo

          enddo

        endif

c  Close file
        if (Node.eq.0) then
          call io_close( iu )
        endif

      elseif (task.eq.'write' .or. task.eq.'WRITE') then

c  Open file
        if (Node.eq.0) then
          call io_assign( iu )
          open( iu, file=fname, form=fform, status='unknown' )      
        endif

        np = mesh(1) * mesh(2) * mesh(3)

c       Write data
        ffor = 'formatted'
        if (Node.eq.0) then
          if (fform .eq. ffor) then
            do i = 1,3
              write(iu,*) (cell(j,i),j=1,3)
            enddo
            write(iu,*) mesh, nspin
          else
            write(iu) cell
            write(iu) mesh, nspin
          endif
        endif

C  Outer loop over spins
        do is = 1,nspin

          Ind = 0

C  Loop over Z dimension of processor grid
          do iz = 1,ProcessorZ

            if (iz.eq.ProcessorZ) then
              NZ = mesh(3) - (ProcessorZ-1)*BlockSizeZ
            else
              NZ = BlockSizeZ
            endif

C  Loop over local Z mesh points
            do izm = 1,NZ

C  Loop over blocks in Y mesh direction
              do iy = 1,ProcessorY

C  Work out size of density sub-matrix to be transfered
                if (iy.eq.ProcessorY) then
                  RemainderY = mesh(2) - (ProcessorY-1)*BlockSizeY
                else
                  RemainderY = BlockSizeY
                endif

#ifdef MPI
                NBlock = RemainderY*mesh(3)
C  Work out which node block is stored on
                BNode = (iy-1)*ProcessorZ + iz - 1

                if (BNode.eq.0.and.Node.eq.BNode) then
#endif
C  If density sub-matrix is local Node 0 then just write it out
                  if (fform .eq. ffor) then
                    do ir = 1,RemainderY
                      write(iu,'(e15.6)') (rho(Ind+ip,is),
     .                  ip=1,mesh(3))
                      Ind = Ind + mesh(3)
                    enddo
                  else
                    do ir = 1,RemainderY
                      write(iu) (rho(Ind+ip,is),ip=1,mesh(3))
                      Ind = Ind + mesh(3)
                    enddo
                  endif

#ifdef MPI
                elseif (Node.eq.0) then
C  If this is Node 0 then recv and write density sub-matrix
                  call MPI_IRecv(bdens,NBlock,DAT_single,BNode,1,
     .              MPI_Comm_World,Request,MPIerror)
                  call MPI_Wait(Request,Status,MPIerror)

                elseif (Node.eq.BNode) then
C  If this is the Node where the density sub-matrix is, then send
                  call MPI_ISend(rho(Ind+1,is),NBlock,DAT_single,0,1,
     .              MPI_Comm_World,Request,MPIerror)
                  call MPI_Wait(Request,Status,MPIerror)
                  Ind = Ind + NBlock

                endif

                if (BNode.ne.0) then
                  call MPI_Barrier(MPI_Comm_World,MPIerror)
                  if (Node.eq.0) then
                    Ind = 0
                    if (fform .eq. ffor) then
                      do ir = 1,RemainderY
                        write(iu,'(e15.6)') (bdens(Ind+ip),ip=1,
     .                    mesh(3))
                        Ind = Ind + mesh(3)
                      enddo
                    else
                      do ir = 1,RemainderY
                        write(iu) (bdens(Ind+ip),ip=1,mesh(3))
                        Ind = Ind + mesh(3)
                      enddo
                    endif
                  endif
                endif
#endif

              enddo

            enddo

          enddo

        enddo

        if (Node.eq.0) then
c       Close file
          call io_close( iu )
        endif

      endif

#ifdef MPI
C Deallocate density buffer memory
      call memory('D','S',size(bdens),'iorho')
      deallocate(bdens)
#endif
      end


