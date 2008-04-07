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
!******************************************************************************
! MODULE m_fft3d
! 3-D fast complex Fourier transform
! Written by J.D.Gale (July 1999), and J.M.Soler (Feb.2008).
!******************************************************************************
!
!   PUBLIC procedures available from this module:
! fft   : 3-D complex FFT. Old version of J.D.Gale
! fft3d : 3-D complex FFT. New version of J.M.Soler
!
!   PUBLIC parameters, types, and variables available from this module:
! none
!
!******************************************************************************
!
!   USED module procedures:
! use sys,    only: die           ! Termination routine
! use alloc,  only: de_alloc      ! De-allocation routine
! use alloc,  only: re_alloc      ! Re-allocation routine
! use mesh3D, only: copyMeshData  ! Copies data in a box of mesh points
! use mesh3D, only: fftMeshDistr  ! Sets a mesh distr. for FFTs
! use mesh3D, only: myMeshBox     ! Returns the mesh box of my processor
! use mesh3D, only: redistributeMeshData ! Changes data distribution
!
!   USED module parameters:
! use precision,  only: dp              ! Real double precision type
! use precision,  only: gp=>grid_p      ! Real precision type of mesh arrays
!
!   EXTERNAL procedures used:
! gpfa    : 1D complex FFT
! setgpfa : Initializes gpfa
! timer   : CPU time counter
!
!******************************************************************************
!
! SUBROUTINE fft3d( dat, meshDistr, nMesh, r2k )
!
! 3D complex FFT
!------------------------------ INPUT -----------------------------------------
! integer meshDistr : Mesh distribution ID of dat array
! integer nMesh(3)  : Mesh divisions of each cell axis
! integer r2k       : Direction of Fourier transform:
!                       r2k=+1 => r->k; r2k=-1 => k->r
!----------------------- INPUT and OUTPUT -------------------------------------
! real(gp) dat(:,:,:,:) : Complex data to be Fourier transformed:
!                           dat(:,:,:,1) = Real part
!                           dat(:,:,:,2) = Imaginary part
!                         The dat array must be declared target.
!----------------------------- UNITS ------------------------------------------
! Units of dat are arbitrary
!----------------------------- USAGE ------------------------------------------
! use precision, only: gp=>grid_p
! use mesh3D,    only: fftMeshDistr, myMeshBox, setMeshDistr
! use m_fft3d,   only: fft3d
! integer:: myBox(2,3), myDistr, nMesh(3)
! real(gp),allocatable:: dat(:,:,:,:)
! ...Find nMesh
! call setMeshDistr( myDistr, nMesh )
!    Or, alternatively
! call fftMeshDistr( nMesh, myDistr )
! call myMeshBox( nMesh, myDistr, myBox )
! allocate( dat(myBox(1,1):myBox(2,1),   &
!               myBox(1,2):myBox(2,2),   &
!               myBox(1,3):myBox(2,3),2) )
! ...Find dat in myBox
! call fft3d( dat, myDistr, nMesh, +1 )
!---------------------------- BEHAVIOUR ---------------------------------------
! It stops with an error message in the following cases:
! - The dimension of dat array in indexes 1-3 is not enough to contain myBox 
!   in meshDistr
! - The dimension of dat in index 4 is different from 2
!--------------------------- ALGORITHMS ---------------------------------------
! 1) dat array is redistributed into kDistr, as returned by fftMeshDistr. 
! 2) dat array is redistributed again into axisDistr(1), as returned from 
!    fftMeshDistr. In this mesh distribution, each processor's box spans the
!    whole x axis (or, more properly, the whole first cell vector).
! 3) dat is Fourier transformed in the x axis
! 4) dat is redistributed again into kDistr
! 5) Steps 2-4 are repeated for the y and z axes.
! 6) dat is redistributed bak into the input meshDistr.
! In serial execution, none of these data redistributions really occur. 
!   Rather, pointers are used to point the 'redistributed' arrays to the
!   original ones.
!
!******************************************************************************

      MODULE m_fft3d

      ! Used module procedures
      use sys,    only: die          ! Terminates execution
      use alloc,  only: de_alloc     ! Deallocates arrays
      use alloc,  only: re_alloc     ! (Re)allocates arrays
      use mesh3D, only: associateMeshTask ! Associates commun. tasks to distr.
      use mesh3D, only: copyMeshData ! Copies a box of mesh data
      use mesh3D, only: fftMeshDistr ! Sets mesh distributions for FFTs
      use mesh3D, only: myMeshBox    ! Returns the mesh box of my node
      use mesh3D, only: redistributeMeshData ! Changes data distribution

      ! Used module parameters
      use precision,    only: dp           ! Real double precision type
      use precision,    only: gp=>grid_p   ! Real type of mesh array data
! DEBUG
      use m_debug,      only: udebug    ! File unit for debug output
! END DEBUG

      PUBLIC:: 
     .  fft,    ! 3D complex FFT. Old JDG version.
     .  fft3d   ! 3D complex FFT. New JMS version.

      PRIVATE ! Nothing is declared public beyond this point

      external:: 
     .  gpfa,    ! 1D complex FFT
     .  setgpfa, ! Initializes gpfa
     .  timer    ! CPU time counter

      CONTAINS

!******************************************************************************

      subroutine fft3d( dat, meshDistr, nMesh, r2k )

      implicit none

      ! Passed arguments
      real(gp),target,intent(inout):: dat(0:,0:,0:,:) ! Data to be Fourier 
                                     ! transformed (one complex function):
                                     !   dat(:,:,:,1) = Real part
                                     !   dat(:,:,:,2) = Imaginary part
      integer,intent(in):: meshDistr ! Mesh distribution ID of dat array
      integer,intent(in):: nMesh(3)  ! Mesh divisions of each cell axis
      integer,intent(in):: r2k       ! Direction of Fourier transform:
                                     ! r2k=+1 => r->k; r2k=-1 => k->r

      ! Internal parameters, variables, and arrays
      character(len=*),parameter:: myName = 'fft3d '
      character(len=*),parameter:: errHead = myName//'ERROR: '
      integer,save:: a2my(3)=-1, aDistr(3)=-1, io2my=-1, my2a(3)=-1, 
     .               myDistr=-1, my2io=-1, oldMesh(3)=0
      logical:: finished, newMesh
      integer:: aBox(2,3), aMesh(3), datShape(4), i1, i2, i3, ia, iter,  
     .          maxTrigs, meshBox(2,3), myBox(2,3), nTrigs(3)
      real(dp),pointer,save:: trigs(:,:)=>null()
      real(gp),pointer:: myDat(:,:,:,:)=>null(), aDat(:,:,:,:)=>null()

      ! Start time counter
      call timer( 'fft3d', 1 )

      ! Get my node's mesh box in input mesh distribution
      call myMeshBox( nMesh, meshDistr, meshBox )

      ! Check dat array shape
      datShape = shape(dat)
      if (any(datShape(1:3) < meshBox(2,:)-meshBox(1,:)+1))
     .  call die(errHead//'dat shape inconsistent with meshDistr')
      if (datShape(4)/=2) 
     .  call die(errHead//'incorrect shape of dat array')

      ! Check if this is a new mesh and store it for later calls
      newMesh = any(nMesh/=oldMesh)
      oldMesh = nMesh

      ! Set trigs array, used by 1D-FFT routine gpfa
      if (newMesh) then
        maxTrigs = 256
        do iter = 1,2
          call re_alloc( trigs, 1,maxTrigs, 1,3, myName//'trigs' )
          do ia=1,3  ! Loop on cell axes
            call setgpfa( trigs(:,ia), maxTrigs, nTrigs(ia), nMesh(ia) )
          end do ! ia
          if (all(nTrigs<=maxTrigs)) then
            finished = .true.
            exit ! iter loop
          else
            finished = .false.
            maxTrigs = maxval(nTrigs)
          end if
        end do ! iter
        if (.not.finished) 
     .    call die(errHead//'trigs iteration not converged')
      end if ! (newMesh)

      ! (Re)set FFT mesh distributions myDistr and aDistr.
      ! myDistr is a homogeneous 3D mesh distribution.
      ! aDistr(1:3) are three homogen. 2D distribs. In each one,
      ! all the node's mesh boxes span one of the whole cell axes
      if (newMesh) call fftMeshDistr( nMesh, myDistr, aDistr )

! BEGIN DEBUG
      ! Check that aDistr mesh boxes span whole cell axes
      if (newMesh) then
        do ia = 1,3
          call myMeshBox( nMesh, aDistr(ia), aBox )
          aMesh(:) = aBox(2,:) - aBox(1,:) + 1
          if (aMesh(ia)/=nMesh(ia)) call die(errHead//'incorrect aBox')
        end do ! ia
      end if ! (newMesh)
! END DEBUG

      ! Associate communication tasks to distributions
      call associateMeshTask( io2my, meshDistr, myDistr )
      call associateMeshTask( my2io, meshDistr, myDistr )
      do ia = 1,3
        call associateMeshTask( my2a(ia), myDistr, aDistr(ia) )
        call associateMeshTask( a2my(ia), myDistr, aDistr(ia) )
      end do ! ia

      ! Redistribute data over myDistr
      ! Notice: if (meshDistr==myDistr) redistributeMeshData may
      ! just direct pointer myDat=>dat, without reallocating it
      call redistributeMeshData( meshDistr, dat, myDistr, myDat, io2my)

      ! FFT in first axis. Arguments of gpfa are:
      ! gpfa( realData(*), imagData(*), trigs(*),  
      !       dataSpan, vectorSpan, vectorSize, numVectors, fftDir )
      call myMeshBox( nMesh, aDistr(1), aBox )
      aMesh(:) = aBox(2,:) - aBox(1,:) + 1
      call redistributeMeshData( myDistr, myDat, aDistr(1), aDat,
     .                           my2a(1) )
      call gpfa( aDat(:,:,:,1), aDat(:,:,:,2), trigs(:,1),
     .           1, aMesh(1), aMesh(1), aMesh(2)*aMesh(3), -r2k )
      call redistributeMeshData( aDistr(1), aDat, myDistr, myDat,
     .                           a2my(1) )
      
      ! FFT in second axis
      call myMeshBox( nMesh, aDistr(2), aBox )
      aMesh(:) = aBox(2,:) - aBox(1,:) + 1
      call redistributeMeshData( myDistr, myDat, aDistr(2), aDat,
     .                           my2a(2) )
      do i3 = 0,aMesh(3)-1
        call gpfa( aDat(:,:,i3,1), aDat(:,:,i3,2), trigs(:,2), 
     .             aMesh(1), 1, aMesh(2), aMesh(1), -r2k )
      end do ! i3
      call redistributeMeshData( aDistr(2), aDat, myDistr, myDat,
     .                           a2my(2) )
      
      ! FFT in third axis
      call myMeshBox( nMesh, aDistr(3), aBox )
      aMesh(:) = aBox(2,:) - aBox(1,:) + 1
      call redistributeMeshData( myDistr, myDat, aDistr(3), aDat,
     .                           my2a(3) )
      call gpfa( aDat(:,:,:,1), aDat(:,:,:,2), trigs(:,3),
     .           aMesh(1)*aMesh(2), 1, aMesh(3), aMesh(1)*aMesh(2), 
     .           -r2k )
      call redistributeMeshData( aDistr(3), aDat, myDistr, myDat,
     .                           a2my(3) )
      
      ! Copy data to output distr. (unless myDat and dat are the same array)
      if (.not.associated(myDat,dat))
     .  call copyMeshData( nMesh, myDistr, myDat, meshBox, dat, my2io )

      ! Divide by number of points
      if (r2k>0) dat = dat / product(nMesh)

      ! Deallocate internal arrays (possibly allocated within
      ! redistributeMeshData)
      if (associated(aDat,myDat)) then
        nullify(aDat)
      else
        call de_alloc( aDat,  routine='redistributeMeshData' )
      end if
      if (associated(myDat,dat)) then
        nullify(myDat)
      else
        call de_alloc( myDat, routine='redistributeMeshData' )
      end if
        
      ! Stop time counter
      call timer( 'fft3d', 2 )

      end subroutine fft3d

!******************************************************************************

      subroutine fft( f, nMesh, isn )
C
C  Parallel FFT based on FFT routine of Templeton
C
C  On input :
C
C  real*8 f()       : contains data to be Fourier transformed
C  integer nMesh(3) : contains global dimensions of grid
C  integer isn      : indicates the direction of the transform
C
C  On exit :
C
C  real*8 f()      : contains data Fourier transformed
C
C  Julian Gale, July 1999
C

C
C  Modules
C
      use precision,    only : dp, grid_p
      use parallel,     only : Node, Nodes, ProcessorY
      use parallelsubs, only : HowManyMeshPerNode
      use sys,          only : die
      use alloc,        only : re_alloc, de_alloc
      use mesh,         only : nsm
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

C
C  Passed arguments
C
      real(grid_p)
     .  f(*)
      integer
     .  nMesh(3), isn
C
C  Local variables
C 
      integer
     .  maxmaxtrigs, ntrigs,
     .  n1, n2, n3, n2l, n3l, CMeshG(3), CMeshL(3), i, n,
     .  ng, nmeshl, ProcessorZ, IOffSet, OldMesh(3)
      integer, save ::
     .  maxtrigs
      logical
     .  LocalPoints, lredimension
      real(dp)
     .  scale
      real(dp), dimension(:,:), pointer, save ::   trigs
#ifdef MPI
      integer
     .  n1lf, nrem, Py, Pz
      real(grid_p), dimension(:), pointer  ::   ft
#endif

      save OldMesh
      data OldMesh/0,0,0/

C
C  Start time counter
C
      call timer( 'fft', 1 )

C
C  Set mesh size variables
C
      n1 = nMesh(1)
      n2 = nMesh(2)
      n3 = nMesh(3)
      CMeshG(1) = n1/nsm
      CMeshG(2) = n2/nsm
      CMeshG(3) = n3/nsm
      call HowManyMeshPerNode(CMeshG, Node, Nodes, nmeshl, CMeshL)
      n2l = CMeshL(2)*nsm
      n3l = CMeshL(3)*nsm

      n = n1*n2l*n3l
      ng = n1*n2*n3

C
C  Set logical as to whether there are any locally stored points
C
      LocalPoints = (n.gt.0) 

C
C  Work out processor grid size
C
      ProcessorZ = Nodes/ProcessorY
      if (ProcessorY*ProcessorZ.ne.Nodes) 
     $     call die('ERROR: ProcessorY must be a factor of the' //
     $     ' number of processors!')

C
C  Allocate trigs array
C
      if (.not.associated(trigs)) then
        maxtrigs = 256
        nullify( trigs )
        call re_alloc( trigs, 1, maxtrigs, 1, 3, name='trigs',
     &                 routine='fft' )
      endif
C
C  Initialise the tables for the FFT if the mesh has changed
C
 10   lredimension = .false.
      maxmaxtrigs = maxtrigs
      do i=1,3
        if (OldMesh(i).ne.nMesh(i)) then
          call setgpfa(trigs(1,i), maxtrigs, ntrigs, nMesh(i))
          if (ntrigs.gt.maxmaxtrigs) then
            lredimension = .true.
            maxmaxtrigs = ntrigs
          else
            OldMesh(i) = nMesh(i)
          endif
        endif
      enddo
      if (lredimension) then
C
C  Resize FFT array for trig factors and set OldMesh to 0 to force recalculation
C
        maxtrigs = maxmaxtrigs
        call re_alloc( trigs, 1, maxtrigs, 1, 3, name='trigs',
     &                 routine='fft' )
        OldMesh(1:3) = 0
        goto 10
      endif

C
C  FFT in X direction
C
      if (LocalPoints) then
        call gpfa(f,f(2),trigs(1,1),2,2*n1,n1,n2l*n3l,-isn)
      endif
#ifdef MPI
C***********************
C  2-D Processor Grid  *
C***********************
      n1lf = n1/ProcessorY
      nrem = n1 - n1lf*ProcessorY
      Py = (Node/ProcessorZ) + 1
      Pz = Node - (Py - 1)*ProcessorZ + 1
      if (Py.le.nrem) n1lf = n1lf + 1
C
C  Allocate local memory
C
      nullify( ft )
      call re_alloc( ft, 1, 2*n1lf*n2*n3l, name='ft', routine='fft' )
C
C  Redistribute data to be distributed by X and Z
C
      call redistribXZ(f,n1,n2l,n3l,ft,n1lf,n2,1,nsm,Node,Nodes)
C
C  FFT in Y direction
C
      if (LocalPoints) then
        do i=0,n3l-1
          IOffSet=2*n1lf*n2*i
          call gpfa(ft(IOffSet+1),ft(IOffSet+2),trigs(1,2),
     .            2*n1lf,2,n2,n1lf,-isn)
        enddo
      endif
C
C  Redistribute data back to original form
C
      call redistribXZ(f,n1,n2l,n3l,ft,n1lf,n2,-1,nsm,Node,Nodes)
C
C  Free local memory ready for re-use
C
      call de_alloc( ft, name='ft' )
C
C  Find new distributed x dimension
C
      n1lf = n1/ProcessorZ
      nrem = n1 - n1lf*ProcessorZ
      if (Pz.le.nrem) n1lf = n1lf + 1
C
C  Allocate local memory
C
      nullify( ft )   ! AG
      call re_alloc( ft, 1, 2*n1lf*n2l*n3, name='ft', routine='fft' )
C
C  Redistribute data to be distributed by X and Y
C
      call redistribXY(f,n1,n2l,n3l,ft,n1lf,n3,1,nsm,Node,Nodes)
C
C  FFT in Z direction
C
      if (LocalPoints) then
        call gpfa(ft,ft(2),trigs(1,3),2*n1lf*n2l,2,n3,n1lf*n2l,-isn)
      endif
C
C  Redistribute data to be distributed by Z again
C
      call redistribXY(f,n1,n2l,n3l,ft,n1lf,n3,-1,nsm,Node,Nodes)
C
C  Free local memory
C
      call de_alloc( ft, name='ft' )
#else
C
C  FFT in Y direction
C
      do i=0,n3-1
         IOffSet=2*n1*n2*i
         call gpfa(f(IOffSet+1),f(IOffSet+2),trigs(1,2),
     .             2*n1,2,n2,n1,-isn)
      enddo
C
C  FFT in Z direction
C
      call gpfa(f,f(2),trigs(1,3),2*n1*n2,2,n3,n1*n2,-isn)
#endif

C
C  Scale values
C
      if (LocalPoints.and.isn.gt.0) then
        scale=1.0_dp/dble(ng)
        do i=1,2*n
          f(i)=f(i)*scale
        enddo
      endif

C
C  Stop time counter
C
      call timer( 'fft', 2 )

      return
      end subroutine fft


!******************************************************************************

#ifdef MPI
      subroutine redistribXZ(f,n1,n2l,n3l,ft,n1lf,n2,idir,nsm,
     .  Node,Nodes)
C
C  This routine redistributes the data over the Nodes as needed
C  for the FFT routines between the arrays f and ft
C
C  Array f is distributed in the Y/Z direction while ft is distributed
C  in the X/Z direction
C
C  idir = direction of redistribution : > 0 => f->ft
C                                       < 0 => ft->f
C
C  Use the processor grid to divide communicator according to Z
C
C  Julian Gale, July 1999
C
      use precision,   only : dp, grid_p
      use mpi_siesta
      use parallel,    only : ProcessorY
      use sys,         only : die
      use alloc,       only : re_alloc, de_alloc

      implicit none

C
C  Passed arguments
C
      integer
     .  n1, n2, n2l, n1lf, n3l, Node, Nodes, idir, nsm
      real(grid_p)
     .  f(2,n1,n2l,n3l), ft(2,n1lf,n2,n3l)
C
C  Local variables
C
      integer
     .  BlockSizeY, BlockSizeYMax, jmin, jmax, jloc, n1lmax, NRemY,
     .  i, j, k, jl, kl, BNode, INode, SNode, jminS, jmaxS

      integer
     .  MPIerror, RequestR, RequestS, Status(MPI_Status_Size)

      integer, save ::
     .  ProcessorZ, Py, Pz, ZCommunicator, ZNode, ZNodes

      logical, save :: firsttimeZ = .true.

      real(grid_p), dimension(:,:,:,:), pointer  :: ftmp,ftmp2

      if (firsttimeZ) then
C
C  Determine processor grid coordinates
C
        ProcessorZ = Nodes/ProcessorY
        if (ProcessorY*ProcessorZ.ne.Nodes) 
     $     call die('ERROR: ProcessorY must be a factor of the' //
     $     ' number of processors!')
        Py = (Node/ProcessorZ) + 1
        Pz = Node - (Py - 1)*ProcessorZ + 1
C
C  Group processors into subsets by Z
C
        call MPI_Comm_Split(MPI_Comm_World, Pz, Py, ZCommunicator,
     .    MPIerror)
        call MPI_Comm_Rank(ZCommunicator,ZNode,MPIerror)
        call MPI_Comm_Size(ZCommunicator,ZNodes,MPIerror)
        firsttimeZ = .false.
      endif

      BlockSizeY = ((n2/nsm)/ProcessorY)*nsm
      NRemY = (n2 - ProcessorY*BlockSizeY)/nsm
      if (NRemY.gt.0) then
        BlockSizeYMax = BlockSizeY + nsm
      else
        BlockSizeYMax = BlockSizeY
      endif
      n1lmax = ((n1-1)/ProcessorY) + 1
C
C  Work out local dimensions
C
      jmin = (Py-1)*BlockSizeY + nsm*min(Py-1,NRemY) + 1
      jmax = jmin + BlockSizeY - 1
      if (Py-1.lt.NRemY) jmax = jmax + nsm
      jmax = min(jmax,n2)
      jloc = jmax - jmin + 1
C
C  Allocate local memory and initialise
C
      nullify( ftmp )
      call re_alloc( ftmp, 1, 2, 1, n1lmax, 1, BlockSizeYMax, 1, n3l,
     &               name='ftmp', routine='redistribXZ' )
      nullify( ftmp2 )
      call re_alloc( ftmp2, 1, 2, 1, n1lmax, 1, BlockSizeYMax, 1, n3l,
     &               name='ftmp2', routine='redistribXZ' )

      do i = 1,n3l
        do j = 1,BlockSizeYMax
          do k = 1,n1lmax
            ftmp(1,k,j,i) = 0.0_grid_p
            ftmp(2,k,j,i) = 0.0_grid_p
            ftmp2(1,k,j,i) = 0.0_grid_p
            ftmp2(2,k,j,i) = 0.0_grid_p
          enddo
        enddo
      enddo

      if (idir.gt.0) then
C***********
C  F -> FT *
C***********
C
C  Handle transfer of terms which are purely local
C
        do i = 1,n3l
          do j = jmin,jmax
            jl = j - jmin + 1
            kl = 0
            do k = 1+ZNode, n1, ZNodes
              kl = kl + 1
              ft(1,kl,j,i) = f(1,k,jl,i)
              ft(2,kl,j,i) = f(2,k,jl,i)
            enddo
          enddo
        enddo
C
C  Loop over all Node-Node vectors exchanging local data
C
        do INode = 1,ProcessorY-1
          BNode = (ZNode+INode)
          BNode = mod(BNode,ProcessorY)
          SNode = (ZNode-INode)
          SNode = mod(SNode+ProcessorY,ProcessorY)
C
C  Collect data to send
C
          do i = 1,n3l
            do jl = 1,jloc
              kl = 0
              do k = 1+BNode, n1, ZNodes
                kl = kl + 1
                ftmp(1,kl,jl,i) = f(1,k,jl,i)
                ftmp(2,kl,jl,i) = f(2,k,jl,i)
              enddo
            enddo
          enddo
C
C  Exchange data - send to right and receive from left
C
          call MPI_IRecv(ftmp2(1,1,1,1),2*n1lmax*BlockSizeYMax*n3l,
     .      MPI_grid_real,SNode,1,ZCommunicator,RequestR,MPIerror)
          call MPI_ISend(ftmp(1,1,1,1),2*n1lmax*BlockSizeYMax*n3l,
     .      MPI_grid_real,BNode,1,ZCommunicator,RequestS,MPIerror)
C
C  Wait for receive to complete
C
          call MPI_Wait(RequestR,Status,MPIerror)
C
C  Place received data into correct array
C
          jminS = SNode*BlockSizeY + nsm*min(SNode,NRemY) + 1
          jmaxS = jminS + BlockSizeY - 1
          if (SNode.lt.NRemY) jmaxS = jmaxS + nsm
          jmaxS = min(jmaxS,n2)
          do i = 1,n3l
            do j = jminS,jmaxS
              jl = j - jminS + 1
              kl = 0
              do k = 1+ZNode, n1, ZNodes
                kl = kl + 1
                ft(1,kl,j,i) = ftmp2(1,kl,jl,i)
                ft(2,kl,j,i) = ftmp2(2,kl,jl,i)
              enddo
            enddo
          enddo
C
C  Wait for send to complete
C
          call MPI_Wait(RequestS,Status,MPIerror)
        enddo
      elseif (idir.lt.0) then
C***********
C  FT -> F *
C***********
C
C  Handle transfer of terms which are purely local
C
        do i = 1,n3l
          do j = jmin,jmax
            jl = j - jmin + 1
            kl = 0
            do k = 1+ZNode, n1, ZNodes
              kl = kl + 1
              f(1,k,jl,i) = ft(1,kl,j,i) 
              f(2,k,jl,i) = ft(2,kl,j,i) 
            enddo
          enddo
        enddo
C
C  Loop over all Node-Node vectors exchanging local data
C
        do INode = 1,ProcessorY-1
          BNode = (ZNode+INode)
          BNode = mod(BNode,ProcessorY)
          SNode = (ZNode-INode)
          SNode = mod(SNode+ProcessorY,ProcessorY)
C
C  Collect data to send
C
          jminS = SNode*BlockSizeY + nsm*min(SNode,NRemY) + 1
          jmaxS = jminS + BlockSizeY - 1
          if (SNode.lt.NRemY) jmaxS = jmaxS + nsm
          jmaxS = min(jmaxS,n2)
          do i = 1,n3l
            do j = jminS,jmaxS
              jl = j - jminS + 1
              kl = 0
              do k = 1+ZNode, n1, ZNodes
                kl = kl + 1
                ftmp(1,kl,jl,i) = ft(1,kl,j,i) 
                ftmp(2,kl,jl,i) = ft(2,kl,j,i) 
              enddo
            enddo
          enddo
C
C  Exchange data - send to right and receive from left
C
          call MPI_IRecv(ftmp2(1,1,1,1),2*n1lmax*BlockSizeYMax*n3l,
     .      MPI_grid_real,BNode,1,ZCommunicator,RequestR,MPIerror)
          call MPI_ISend(ftmp(1,1,1,1),2*n1lmax*BlockSizeYMax*n3l,
     .      MPI_grid_real,SNode,1,ZCommunicator,RequestS,MPIerror)
C
C  Wait for receive to complete
C
          call MPI_Wait(RequestR,Status,MPIerror)
C
C  Place received data into correct array
C
          do i = 1,n3l
            do jl = 1,jloc
              kl = 0
              do k = 1+BNode, n1, ZNodes
                kl = kl + 1
                f(1,k,jl,i) = ftmp2(1,kl,jl,i) 
                f(2,k,jl,i) = ftmp2(2,kl,jl,i) 
              enddo
            enddo
          enddo
C
C  Wait for send to complete
C
          call MPI_Wait(RequestS,Status,MPIerror)
        enddo
      endif
C
C  Free local memory
C
      call de_alloc( ftmp2, name='ftmp2' )
      call de_alloc( ftmp, name='ftmp' )

      end subroutine redistribXZ

!******************************************************************************

      subroutine redistribXY(f,n1,n2l,n3l,ft,n1lf,n3,idir,nsm,
     .  Node,Nodes)
C
C  This routine redistributes the data over the Nodes as needed
C  for the FFT routines between the arrays f and ft
C
C  Array f is distributed in the Y/Z direction while ft is distributed
C  in the X/Y direction
C
C  idir = direction of redistribution : > 0 => f->ft
C                                       < 0 => ft->f
C
C  Use the processor grid to divide communicator according to Y
C
C  Currently written in not the most efficient but simple way! 
C  Need to improve communication later.
C
C  Julian Gale, July 1999
C
      use precision,   only: dp, grid_p
      use mpi_siesta
      use parallel,    only : ProcessorY
      use sys,         only : die
      use alloc,       only : re_alloc, de_alloc

      implicit none

C
C  Passed arguments
C
      integer
     .  n1, n3, n2l, n1lf, n3l, nsm, Node, Nodes, idir
      real(grid_p)
     .  f(2,n1,n2l,n3l), ft(2,n1lf,n2l,n3)
C
C  Local variables
C
      integer
     .  i, j, k, il, kl, BNode, BlockSizeZ, imin, imax, iloc,
     .  INode, n1lmax, SNode, iminS, imaxS, NRemZ, BlockSizeZMax

      integer
     .  MPIerror, RequestR, RequestS, Status(MPI_Status_Size)

      integer, save ::
     .  ProcessorZ, Py, Pz, YCommunicator, YNode, YNodes

      logical, save :: firsttimeY = .true.

      real(grid_p), dimension(:,:,:,:), pointer :: ftmp,ftmp2

      if (firsttimeY) then
C
C  Determine processor grid coordinates
C
        ProcessorZ = Nodes/ProcessorY
        if (ProcessorY*ProcessorZ.ne.Nodes) 
     $     call die('ERROR: ProcessorY must be a factor of the' //
     $     ' number of processors!')
        Py = (Node/ProcessorZ) + 1
        Pz = Node - (Py - 1)*ProcessorZ + 1
C
C  Group processors into subsets by Y
C
        call MPI_Comm_Split(MPI_Comm_World, Py, Pz, YCommunicator,
     .    MPIerror)
        call MPI_Comm_Rank(YCommunicator,YNode,MPIerror)
        call MPI_Comm_Size(YCommunicator,YNodes,MPIerror)
        firsttimeY = .false.
      endif

      BlockSizeZ = ((n3/nsm)/ProcessorZ)*nsm
      NRemZ = (n3 - ProcessorZ*BlockSizeZ)/nsm
      if (NRemZ.gt.0) then
        BlockSizeZMax = BlockSizeZ + nsm
      else
        BlockSizeZMax = BlockSizeZ
      endif
      n1lmax = ((n1-1)/ProcessorZ) + 1
C
C  Allocate local memory and initialise
C
      nullify( ftmp )
      call re_alloc( ftmp, 1, 2, 1, n1lmax, 1, n2l, 1, BlockSizeZMax,
     &               name='ftmp', routine='redistribXY' )
      nullify( ftmp2 )
      call re_alloc( ftmp2, 1, 2, 1, n1lmax, 1, n2l, 1, BlockSizeZMax,
     &               name='ftmp2', routine='redistribXY' )

      do i = 1,BlockSizeZMax
        do j = 1,n2l
          do k = 1,n1lmax
            ftmp(1,k,j,i) = 0.0_grid_p
            ftmp(2,k,j,i) = 0.0_grid_p
            ftmp2(1,k,j,i) = 0.0_grid_p
            ftmp2(2,k,j,i) = 0.0_grid_p
          enddo
        enddo
      enddo
C
C  Work out local dimensions
C
      imin = (Pz-1)*BlockSizeZ + nsm*min(Pz-1,NRemZ) + 1
      imax = imin + BlockSizeZ - 1
      if (Pz-1.lt.NRemZ) imax = imax + nsm
      imax = min(imax,n3)
      iloc = imax - imin + 1

      if (idir.gt.0) then
C***********
C  F -> FT *
C***********
C
C  Handle transfer of terms which are purely local
C
        do i = imin,imax
          il = i - imin + 1
          do j = 1,n2l
            kl = 0
            do k = 1+YNode, n1, YNodes
              kl = kl + 1
              ft(1,kl,j,i) = f(1,k,j,il)
              ft(2,kl,j,i) = f(2,k,j,il)
            enddo
          enddo
        enddo
C
C  Loop over all Node-Node vectors exchanging local data
C
        do INode = 1,ProcessorZ-1
          BNode = (YNode+INode)
          BNode = mod(BNode,ProcessorZ)
          SNode = (YNode-INode)
          SNode = mod(SNode+ProcessorZ,ProcessorZ)
C
C  Collect data to send
C
          do il = 1,iloc
            do j = 1,n2l
              kl = 0
              do k = 1+BNode, n1, YNodes
                kl = kl + 1
                ftmp(1,kl,j,il) = f(1,k,j,il)
                ftmp(2,kl,j,il) = f(2,k,j,il)
              enddo
            enddo
          enddo
C
C  Exchange data - send to right and receive from left
C
          call MPI_IRecv(ftmp2(1,1,1,1),2*n1lmax*n2l*BlockSizeZMax,
     .      MPI_grid_real,SNode,1,YCommunicator,RequestR,MPIerror)
          call MPI_ISend(ftmp(1,1,1,1),2*n1lmax*n2l*BlockSizeZMax,
     .      MPI_grid_real,BNode,1,YCommunicator,RequestS,MPIerror)
C
C  Wait for receive to complete
C
          call MPI_Wait(RequestR,Status,MPIerror)
C
C  Place received data into correct array
C
          iminS = SNode*BlockSizeZ + nsm*min(SNode,NRemZ) + 1
          imaxS = iminS + BlockSizeZ - 1
          if (SNode.lt.NRemZ) imaxS = imaxS + nsm
          imaxS = min(imaxS,n3)
          do i = iminS,imaxS
            il = i - iminS + 1
            do j = 1,n2l
              kl = 0
              do k = 1+YNode, n1, YNodes
                kl = kl + 1
                ft(1,kl,j,i) = ftmp2(1,kl,j,il)
                ft(2,kl,j,i) = ftmp2(2,kl,j,il)
              enddo
            enddo
          enddo
C
C  Wait for send to complete
C
          call MPI_Wait(RequestS,Status,MPIerror)
        enddo
      elseif (idir.lt.0) then
C***********
C  FT -> F *
C***********
C
C  Handle transfer of terms which are purely local
C
        do i = imin,imax
          il = i - imin + 1
          do j = 1,n2l
            kl = 0
            do k = 1+YNode, n1, YNodes
              kl = kl + 1
              f(1,k,j,il) = ft(1,kl,j,i)
              f(2,k,j,il) = ft(2,kl,j,i)
            enddo
          enddo
        enddo
C
C  Loop over all Node-Node vectors exchanging local data
C
        do INode = 1,ProcessorZ-1
          BNode = (YNode+INode)
          BNode = mod(BNode,ProcessorZ)
          SNode = (YNode-INode)
          SNode = mod(SNode+ProcessorZ,ProcessorZ)
C
C  Collect data to send
C
          iminS = SNode*BlockSizeZ + nsm*min(SNode,NRemZ) + 1
          imaxS = iminS + BlockSizeZ - 1
          if (SNode.lt.NRemZ) imaxS = imaxS + nsm
          imaxS = min(imaxS,n3)
          do i = iminS,imaxS
            il = i - iminS + 1
            do j = 1,n2l
              kl = 0
              do k = 1+YNode, n1, YNodes
                kl = kl + 1
                ftmp(1,kl,j,il) = ft(1,kl,j,i)
                ftmp(2,kl,j,il) = ft(2,kl,j,i)
              enddo
            enddo
          enddo
C
C  Exchange data - send to right and receive from left
C
          call MPI_IRecv(ftmp2(1,1,1,1),2*n1lmax*n2l*BlockSizeZMax,
     .      MPI_grid_real,BNode,1,YCommunicator,RequestR,MPIerror)
          call MPI_ISend(ftmp(1,1,1,1),2*n1lmax*n2l*BlockSizeZMax,
     .      MPI_grid_real,SNode,1,YCommunicator,RequestS,MPIerror)
C
C  Wait for receive to complete
C
          call MPI_Wait(RequestR,Status,MPIerror)
C
C  Place received data into correct array
C
          do il = 1,iloc
            do j = 1,n2l
              kl = 0
              do k = 1+BNode, n1, YNodes
                kl = kl + 1
                f(1,k,j,il) = ftmp2(1,kl,j,il)
                f(2,k,j,il) = ftmp2(2,kl,j,il)
              enddo
            enddo
          enddo
C
C  Wait for send to complete
C
          call MPI_Wait(RequestS,Status,MPIerror)
        enddo
      endif
C
C  Free local memory
C
      call de_alloc( ftmp2, name='ftmp2' )
      call de_alloc( ftmp, name='ftmp' )

      return
      end subroutine redistribXY
#endif

!******************************************************************************

      END MODULE m_fft3d
