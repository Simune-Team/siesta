      module meshdscf
C
C Stores quantities that are connected with Dscf in mesh local
C form when data is distributed for parallel execution
C
      implicit none

C ----------------------------------------------------------------------
C Dscf related variables for parallel distributed form
C ----------------------------------------------------------------------
C integer listdl(sum(numdl))           : List of non-zero elements in
C                                      : a row of DscfL
C integer listdlptr(nrowsDscfL)        : Pointer to row in listdl
C integer NeedDscfL(nuotot)            : Pointer as to whether a row of
C                                      : Dscf is needed in DscfL
C integer nrowsDscfL                   : Number of rows of DscfL
C integer numdl(nrowsDscfL)            : Number of non-zero elements in
C                                      : a row of DscfL
C real*8  DscfL(maxndl,nrowsDscfL)     : Local copy of Dscf elements
C                                      : needed for the local mesh
C ----------------------------------------------------------------------

      integer, save :: nrowsDscfL

      integer, dimension(:), allocatable, save :: listdl
      integer, dimension(:), allocatable, save :: listdlptr
      integer, dimension(:), allocatable, save :: NeedDscfL
      integer, dimension(:), allocatable, save :: numdl

      real*8,  dimension(:,:), allocatable, save :: DscfL

      end module meshdscf

      subroutine CreateLocalDscfPointers( nmpl, nuotot, indxuo, numd,
     .                                    listdptr, listd, Node, Nodes)
C
C Calculates the values of the orbitals at the mesh points
C
C ----------------------------------------------------------------------
C Input :
C ----------------------------------------------------------------------
C integer nmpl          : Number of mesh points in unit cell locally
C integer nuotot        : Total number of basis orbitals in unit cell
C integer indxuo(nuotot): Pointer from orbital in supercell to unit cell
C integer numd(nuo)     : Number of nonzero density-matrix
C                       : elements for each matrix row
C integer listdptr(nuo) : Pointer to start of rows of density-matrix
C integer listd(maxnh)  : Nonzero-density-matrix-element column
C                       : indexes for each matrix row
C integer Node          : Local node number
C integer Nodes         : Total number of nodes
C ----------------------------------------------------------------------
C Output :
C ----------------------------------------------------------------------
C All output quantities are in the module meshdscf
C ----------------------------------------------------------------------

C
C Modules
C
      use meshdscf, only: listdl, listdlptr, NeedDscfL, nrowsDscfL, 
     .                    numdl
      use meshphi,  only: endpht, lstpht
      use parallel
      use precision
#ifdef MPI
      use mpi
#endif

      implicit none

C
C Passed arguments
C
      integer, intent(in) :: nmpl
      integer, intent(in) :: nuotot
      integer, intent(in) :: indxuo(nuotot)
      integer, intent(in) :: numd(*)
      integer, intent(in) :: listdptr(*)
      integer, intent(in) :: listd(*)
      integer, intent(in) :: Node
      integer, intent(in) :: Nodes

C
C Local variables
C
      integer
     .  BNode, i, ii, io, iio, ip, imp, iu, numdele, maxndmax,
     .  nsize
#ifdef MPI
      integer
     .  MPIerror
#endif

      integer, dimension(:), allocatable :: ibuffer

C Create pointer as to whether a given row of DscfL is needed in NeedDscfL
      if (.not.allocated(NeedDscfL)) then
        allocate(NeedDscfL(nuotot))
        call memory('A','I',nuotot,'CreateLocalDscfPointers')
      endif
      NeedDscfL(1:nuotot) = 0
      do ip = 1,nmpl
        do imp = 1+endpht(ip-1), endpht(ip)
          i = lstpht(imp)
          iu = indxuo(i)
          NeedDscfL(iu) = 1
        enddo
      enddo
      nrowsDscfL = 0
      do i = 1,nuotot
        if (NeedDscfL(i).eq.1) then
          nrowsDscfL = nrowsDscfL + 1
          NeedDscfL(i) = nrowsDscfL
        endif
      enddo

C Allocate/reallocate memory for numdl and listdlptr
      if (.not.allocated(numdl).or.size(numdl).ne.nrowsDscfL) then
        if (allocated(numdl)) then
          call memory('D','I',size(numdl),'CreateLocalDscfPointers')
          deallocate(numdl)
          call memory('D','I',size(listdlptr),'CreateLocalDscfPointers')
          deallocate(listdlptr)
        endif
        allocate(numdl(max(1,nrowsDscfL)))
        call memory('A','I',max(1,nrowsDscfL),'CreateLocalDscfPointers')
        allocate(listdlptr(max(1,nrowsDscfL)))
        call memory('A','I',max(1,nrowsDscfL),'CreateLocalDscfPointers')
      endif

C Distribute information about numd globally
      maxndmax = 0
      do io = 1,nuotot
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          numdele = numd(iio)
        endif
#ifdef MPI
        call MPI_Bcast(numdele,1,MPI_integer,BNode,
     .    MPI_Comm_World,MPIerror)
#endif
        if (NeedDscfL(io).gt.0) numdl(NeedDscfL(io)) = numdele
        if (numdele.gt.maxndmax) maxndmax = numdele
      enddo

C Create listdlptr using numdl
      listdlptr(1) = 0
      do io = 2,nrowsDscfL
        listdlptr(io) = listdlptr(io-1) + numdl(io-1)
      enddo

C Allocate/reallocate listdl
      if (allocated(listdl)) then
        call memory('D','I',size(listdl),'CreateLocalDscfPointers')
        deallocate(listdl)
      endif
      if (nrowsDscfL.gt.0) then
        nsize = listdlptr(nrowsDscfL)+numdl(nrowsDscfL)
      else
        nsize = 1
      endif
      allocate(listdl(nsize))
      call memory('A','I',nsize,'CreateLocalDscfPointers')

C Distribute information about listd globally
      allocate(ibuffer(maxndmax))
      call memory('A','I',maxndmax,'CreateLocalDscfPointers')
      ibuffer(1:maxndmax) = 0
      do io = 1,nuotot
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          do ii = 1,numd(iio)
            ibuffer(ii) = listd(listdptr(iio)+ii)
          enddo
        endif
#ifdef MPI
        call MPI_Bcast(ibuffer,maxndmax,MPI_integer,BNode,
     .    MPI_Comm_World,MPIerror)
#endif
        iio = NeedDscfL(io)
        if (iio.gt.0) then
          do ii = 1,numdl(iio)
            listdl(listdlptr(iio)+ii) = ibuffer(ii)
          enddo
        endif
      enddo
      call memory('D','I',size(ibuffer),'CreateLocalDscfPointers')
      deallocate(ibuffer)

      end subroutine CreateLocalDscfPointers

      subroutine matrixOtoM( maxnd, numd, listdptr, maxndl, nuo, 
     .                       nuotot, nspin, Dscf, DscfL )
C ********************************************************************
C Transforms a matrix which is distributed by block cyclic 
C distribution of orbitals to a matrix that contains all
C the orbital rows needed for a mesh point distribution 
C over the nodes.
C Created by J.D.Gale, February 2000
C *********************** INPUT **************************************
C integer maxnd         : First dimension of Dscf
C integer numd(nuo)     : Number of non-zero elements in row of Dscf
C integer listdptr(nuo) : Pointer to start of rows in Dscf
C integer maxndl        : First dimension of DscfL
C integer nuo           : Local no. of orbitals in unit cell
C integer nuotot        : Total no. of orbitals in unit cell
C integer nspin         : Number of spin components
C real*8  Dscf(maxnd,nspin) : Matrix in orbital distributed form
C *********************** OUTPUT *************************************
C real*8  DscfL(maxndl,nspin) : Matrix in mesh distributed form
C ********************************************************************

C  Modules
      use precision
      use meshdscf, only: nrowsDscfL, numdl, listdlptr, NeedDscfL
#ifdef MPI
      use mpi
      use parallel
#endif

      implicit none

C Argument types and dimensions
      integer
     .   maxnd, maxndl, nspin, nuo, nuotot, numd(nuo), 
     .   listdptr(nuo)

      double precision
     .   Dscf(maxnd,nspin), DscfL(maxndl,nspin)

      external
     .   memory

C Internal variables and arrays
      integer
     .  ii, io, iio, il, ispin

#ifdef MPI
      integer 
     .  BNode, maxno, maxnog, MPIerror, Node, Nodes
      double precision, dimension(:), allocatable, save :: 
     .  buffer
#endif

C***********************
C  Parallel execution  *
C***********************
#ifdef MPI
C Get Node number
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)

C Find local and global size of maxno
      maxno = 0
      do io = 1,nrowsDscfL
        maxno = max(maxno,numdl(io))
      enddo
      call MPI_AllReduce(maxno,maxnog,1,MPI_integer,MPI_max,
     .  MPI_Comm_World,MPIerror)

C Allocate local Dscf storage array
      allocate(buffer(maxnog*nspin))
      call memory('A','D',maxnog*nspin,'matrixOtoM')

C Zero buffer as we will be passing the complete array
C as the number of explicitly non-zero elements will
C not be known on all nodes
      buffer(1:maxnog*nspin) = 0.0d0

C Loop over rows of Dscf broadcasting to all other nodes
      do io = 1,nuotot
        call WhichNodeOrb(io,Nodes,BNode)

C If this row is local to this node then copy into buffer
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          do ispin = 1,nspin
            ii = (ispin-1)*maxnog
            buffer(ii+1:ii+numd(iio)) = 
     .        Dscf(listdptr(iio)+1:listdptr(iio)+numd(iio),ispin)
          enddo
        endif

C Broadcast buffer array
        call MPI_Bcast(buffer,maxnog*nspin,DAT_double,BNode,
     .                 MPI_Comm_World,MPIerror)

C Get pointer for this row of Dscf and see if it is needed for DscfL
        il = NeedDscfL(io)
        if (il.gt.0) then
          do ispin = 1,nspin
            ii = (ispin-1)*maxnog
            DscfL(listdlptr(il)+1:listdlptr(il)+numdl(il),ispin) = 
     .        buffer(ii+1:ii+numdl(il))
          enddo
        endif
      enddo

C Deallocate buffer memory
      call memory('D','D',size(buffer),'matrixOtoM')
      deallocate(buffer)
#else
C*********************
C  Serial execution  *
C*********************
C Loop over rows of Dscf checking to see if they are in DscfL
      do ispin = 1,nspin
        do io = 1,nuotot

C Get pointer for this row of Dscf and see if it is needed for DscfL
          il = NeedDscfL(io)
          if (il.gt.0) then
            DscfL(listdlptr(il)+1:listdlptr(il)+numdl(il),ispin) = 
     .        Dscf(listdptr(io)+1:listdptr(io)+numdl(il),ispin)
          endif

        enddo
      enddo
#endif

      end subroutine matrixOtoM

      subroutine matrixMtoO( maxnvl, maxnv, numVs, listVsptr, nuo, 
     .                       nuotot, nspin, VsL, Vs )

C ********************************************************************
C Transforms a matrix which is distributed by mesh points to a matrix
C that is distributed by a block cyclic distribution over the orbitals
C and adds it to an existing array of this form.
C Created by J.D.Gale, February 2000
C *********************** INPUT **************************************
C integer maxnvl          : First dimension of VsL and maximum number
C                           of nonzero elements in VsL
C integer maxnv           : First dimension of Vs and maximum number
C                           of nonzero elements in Vs
C integer numVs(nuo)      : Number of non-zero elements in row of Vs
C integer listVsptr(nuo)  : Pointer to start of rows in Vs
C integer nuo             : Local no. of orbitals in unit cell
C integer nuotot          : Total no. of orbitals in unit cell
C integer nspin           : Number of spin components
C real*8  VsL(maxnvl,nspin) : Mesh contribution to be added to Vs
C ******************** INPUT AND OUTPUT *******************************
C real*8  Vs(maxnv,nspin) : Value of nonzero elements in each row of Vs
C                           to which the potential matrix elements are
C                           summed up
C *********************************************************************

C  Modules
      use precision
      use meshdscf
#ifdef MPI
      use mpi
      use parallel
#endif

      implicit none

C Argument types and dimensions
      integer
     .   maxnv, maxnvl, nspin, nuo, nuotot, numVs(nuo), 
     .   listVsptr(nuo)
      double precision
     .   Vs(maxnv,nspin), VsL(maxnvl,nspin)

C Internal variables and arrays
      integer
     .  i, iu, ii, ispin

#ifdef MPI
      integer 
     .  in, iul, maxnvg, MPIerror, nVsL, nVsLmax, Node, Nodes
      double precision, dimension(:), allocatable, save :: 
     .  Vi
      integer, dimension(:), allocatable, save :: 
     .  nVsLPtr, listViptr
#endif
      
C***********************
C  Parallel execution  *
C***********************
#ifdef MPI
C Get Node number
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)

C Find the maximum number of rows on anyone Node
      call MPI_AllReduce(nrowsDscfL,nVsLmax,1,MPI_integer,
     .  MPI_max,MPI_Comm_World,MPIerror)

C Find the maximum value of maxno on any Node and globally
      call MPI_AllReduce(maxnvl,maxnvg,1,MPI_integer,
     .  MPI_max,MPI_Comm_World,MPIerror)

C Allocate buffer memory for transfer of matrix and pointer
      allocate(Vi(maxnvg*nspin))
      call memory('A','D',maxnvg*nspin,'matrixMtoO')
      allocate(nVsLPtr(nuotot))
      call memory('A','I',nuotot,'matrixMtoO')
      allocate(listViptr(nuotot))
      call memory('A','I',nuotot,'matrixMtoO')

C Initialise buffers to avoid number trapping errors during transfer
      Vi(1:maxnvg) = 0.0d0
      listViptr(1:nuotot) = 0

C Loop over Nodes for broadcasting of local data
      do in = 0,Nodes-1

C Broadcast the number of rows stored locally
        if (in.eq.Node) nVsL = nrowsDscfL
        call MPI_Bcast(nVsL,1,MPI_integer,in,MPI_Comm_World,MPIerror)

C If this is the broadcasting Node then copy the pointer to the buffer
        if (in.eq.Node) then
          nVsLPtr(1:nuotot) = NeedDscfL(1:nuotot)
        endif

C Broadcast the pointer information
        call MPI_Bcast(nVsLPtr,nuotot,MPI_integer,in,MPI_Comm_World,
     .                 MPIerror)

C Need to broadcast the row pointer info for Vi (listdlptr)
        if (in.eq.Node.and.nrowsDscfL.gt.0) then
          listViptr(1:nrowsDscfL) = listdlptr(1:nrowsDscfL)
        endif
        call MPI_Bcast(listViptr,nuotot,MPI_integer,in,MPI_Comm_World,
     .                 MPIerror)

C If this is the broadcasting Node then copy the data to the buffer
        if (in.eq.Node) then
          do ispin = 1,nspin
            ii = (ispin-1)*maxnvg
            Vi(ii+1:ii+maxnvl) = VsL(1:maxnvl,ispin)
          enddo
        endif

C Broadcast the buffer information
        call MPI_Bcast(Vi,maxnvg*nspin,DAT_double,in,
     .                 MPI_Comm_World,MPIerror)

C Add those elements that are needed locally to the values already
C stored in the orbital oriented array
        do ispin = 1,nspin
          ii = (ispin-1)*maxnvg
          do i = 1,nuo
            call LocalToGlobalOrb(i,Node,Nodes,iu)
            iul = nVsLPtr(iu)
            if (iul.gt.0) then
              Vs(listVsptr(i)+1:listVsptr(i)+numVs(i),ispin) = 
     .          Vs(listVsptr(i)+1:listVsptr(i)+numVs(i),ispin) + 
     .          Vi(ii+listViptr(iul)+1:ii+listViptr(iul)+numVs(i))
            endif
          enddo
        enddo

      enddo

C Deallocate buffer and pointer memory
      call memory('D','I',size(listViptr),'matrixMtoO')
      deallocate(listViptr)
      call memory('D','I',size(nVsLPtr),'matrixMtoO')
      deallocate(nVsLPtr)
      call memory('D','D',size(Vi),'matrixMtoO')
      deallocate(Vi)
#else
C*********************
C  Serial execution  *
C*********************
C Add those elements that are needed locally to the values already
C stored in the orbital oriented array
      do ispin = 1,nspin
        do i = 1,nuotot
          iu = NeedDscfL(i)
          if (iu.gt.0) then
            Vs(listVsptr(i)+1:listVsptr(i)+numVs(i),ispin) = 
     .        Vs(listVsptr(i)+1:listVsptr(i)+numVs(i),ispin) + 
     .        VsL(listdlptr(iu)+1:listdlptr(iu)+numVs(i),ispin)
          endif
        enddo
      enddo
#endif

      end
