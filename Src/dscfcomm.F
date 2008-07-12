 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      module m_dscfComm
C ----------------------------------------------------------------------
C Stores quantities used to transfer data between a matrix that is
C distributed by orbitals and one distributed by mesh
C ----------------------------------------------------------------------
C integer DCsize        : Number of communications
C integer DCself        : Number of orbitals that the current process should keep
C integer DCtotal       : Total Number of orbitals of the local matrix
C integer DCmax         : Maximum number of orbitals that we should send/receive
C                         from/to another process
C integer DCmaxnd       : Maximum size of matrix that we should send/receive
C                         from/to another process
C integer DCBuffer(:)   : Buffer that contains all the data
C integer DCsrc(DCsize) : Process ID of the source.
C integer DCdst(DCsize) : Process ID of the destiny.
C integer DCsic(DCsize) : Number of orbitals to send/receive
C integer DCinvp(DCtotal) : Permutation vector used to transform the received
C                           data to the local ordering
C ----------------------------------------------------------------------
C
C     Modules
      use precision
      implicit none
      integer,          save :: DCsize, DCself, DCtotal, DCmax, DCmaxnd
      integer, pointer, save :: DCBuffer(:)
      integer, pointer, save :: DCsrc(:), DCdst(:), DCsic(:), DCinvp(:)
      end module m_dscfComm


C ==================================================================
C Precomputes the needed communications to transform a matrix
C that is distributed by orbitals to one distributed by mesh
C ==================================================================
C subroutine dscfComm( nuotot, nrowsDscfL, NeedDscfL )
C
C INPUT:
C integer nuotot            : Total number of basis orbitals in unit cell
C integer nrowsDscfL        : Local number of rows of matrix DSCF
C integer NeedDscfL(nuotot) : Permutation vector of orbitals from global to local.
C
C OUTPUT:
C The output values are in the module m_dscfComm
C
C BEHAVIOR:
C Count the number of orbitals that we should receive from every process.
C Create a list of communications that every process should do. The process
C 0 receive all the communications and calls scheduleComm to planify
C the communication schedule. The rest of processes receive the schedule
C from process 0.
C
C ==================================================================
      subroutine dscfComm( nuotot, nrowsDscfL, NeedDscfL )
      use precision
      use parallel,     only : NODE, NODES
      use parallelsubs, only : GlobalToLocalOrb, WhichNodeOrb
      use alloc
      use scheComm
      use m_dscfComm
#ifdef MPI
      use mpi_siesta
#endif
      implicit none
C     Input variables
      integer               :: nuotot, nrowsDscfL
      integer               :: NeedDscfL(nuotot)

C     Local variables
      integer               :: ii, io, ANode, BNode, PP, ncom, com, pos,
     &                         tmp, nn
      integer, pointer      :: neights(:), Gneights(:,:), src(:),
     &                         dst(:), sereBUF(:), procSrc(:),
     &                         invpTmp(:), invp(:), xneights(:)
      type(COMM_T)          :: comm
#ifdef MPI
      integer               :: MPIerror, Status(MPI_Status_Size)
#endif
!------------------------------------------------------------------------- BEGIN
      call resetdscfComm( )

      nullify( neights, xneights, invpTmp, invp, procSrc )
      call re_alloc( xneights, 0, NODES, name='xneights',
     &               routine='dscfComm' )
      call re_alloc( invp, 1, nrowsDscfL, name='invp',
     &               routine='dscfComm' )
      call re_alloc( neights, 0, NODES-1, name='neights',
     &                 routine='dscfComm' )
      call re_alloc( invpTmp, 1, nrowsDscfL, name='invpTmp',
     &               routine='dscfComm' )
      call re_alloc( procSrc, 1, nrowsDscfL, name='procSrc',
     &               routine='dscfComm' )
      neights = 0
      ii      = 0
C     Count the number of rows of DSCFL that belongs to each domain
C     Create a permutation vector in order to recover the original
C     ID of every orbital. Associated to every orbital there is a PID.
      do io= 1, nuotot
        if (NeedDscfL(io).ne.0) then
          call WhichNodeOrb(io,Nodes,BNode)
          neights(BNode) = neights(BNode) + 1
          ii             = ii + 1
          invpTmp(ii)    = io
          procSrc(ii)    = BNode
        endif
      enddo

C     Reorder the permutation vector giving its PID
      xneights(0) = 1
      do BNode=1, NODES
        xneights(BNode) = xneights(BNode-1) + neights(BNode-1)
      enddo
      do io= 1, nrowsDscfL
        BNode           = procSrc(io)
        ii              = xneights(BNode)
        invp(ii)        = invpTmp(io)
        xneights(BNode) = xneights(BNode) + 1
      enddo
      xneights(0) = 1
      do BNode=1, NODES
        xneights(BNode) = xneights(BNode-1) + neights(BNode-1)
      enddo

      call de_alloc( invpTmp, name='invpTmp', routine='dscfComm' )
      call de_alloc( procSrc, name='procSrc', routine='dscfComm' )

C     Compute the communications needed. The master process will gather
C     the vector "neights" from all the other processes
      nullify( Gneights )
      if (Node.eq.0) then
        call re_alloc( Gneights, 0, NODES-1, 0, NODES-1,
     &                 name='Gneights', routine='dscfComm' )
      else
        call re_alloc( Gneights, 0, 0, 0, 0,
     &                 name='Gneights', routine='dscfComm' )
      endif

#ifdef MPI
      call MPI_Gather( neights, Nodes, MPI_INTEGER,
     &                 Gneights(0,0), Nodes, MPI_INTEGER,
     &                 0, MPI_COMM_WORLD, MPIerror )
#else
      Gneights(0,0) = neights(0)
#endif
      call de_alloc( neights, name='neights', routine='dscfComm' )

      if (Node.eq.0) then
C       Count the number of communications
        ncom = 0
        do ANode=0, Nodes-1
          do BNode=0, Nodes-1
            if (Anode.ne.BNode .and. Gneights(BNode,ANode).ne.0 ) then
              ncom = ncom + 1
            endif
          enddo
        enddo

        nullify( src, dst )
        call re_alloc( src, 1, ncom, name='src', routine='dscfComm' )
        call re_alloc( dst, 1, ncom, name='dst', routine='dscfComm' )

C       Compute the communications
        ncom = 0
        do ANode=0, Nodes-1
          do BNode=0, Nodes-1
            if (Anode.ne.BNode .and. Gneights(BNode,ANode).ne.0 ) then
              ncom = ncom + 1
              src(ncom) = BNode+1
              dst(ncom) = ANode+1
            endif
          enddo
        enddo

        comm%np = Nodes
C       reschedule the communications in order to minimize the time
        call scheduleComm( ncom, src, dst, comm )

        nullify( sereBUF )
        call re_alloc( sereBUF, 1, comm%ncol*3, name='sereBUF',
     &                 routine='dscfComm' )

C       Send the communication scheduling to the other processes
        do PP=2, Nodes
          com = 0
          do io= 1, comm%Ncol
           if (comm%ind(io,PP).ne.0) com = com + 1
          enddo
          Ncom = com
          
          com = 0
          do io= 1, comm%Ncol
            pos = comm%ind(io,PP)
            if (pos.ne.0) then
              com = com + 1
              ANode = src(pos)-1
              BNode = dst(pos)-1
              sereBUF(com       ) = ANode
              sereBUF(com+  Ncom) = BNode
              sereBUF(com+2*Ncom) = Gneights(ANode,BNode)
            endif
          enddo
#ifdef MPI
          call MPI_Send( Ncom, 1, MPI_integer,
     &                   PP-1, 1, MPI_Comm_World, MPIerror )
          call MPI_Send( sereBUF, Ncom*3, MPI_integer,
     &                   PP-1, 1, MPI_Comm_World, MPIerror )
#endif
        enddo
        call de_alloc( sereBUF, name='sereBUF', routine='dscfComm' )

C       Compute the master communication scheduling
        PP  = 1
        com = 0
        do io= 1, comm%Ncol
         if (comm%ind(io,PP).ne.0) com = com + 1
        enddo
        DCsize = com

        call re_alloc( DCBuffer, 1, DCsize*3, name='DCBuffer',
     &                 routine='dscfComm' )

        DCsrc => DCBuffer(         1:  DCsize)
        DCdst => DCBuffer(  DCsize+1:2*DCsize)
        DCsic => DCBuffer(2*DCsize+1:3*DCsize)

        com = 0
        do io= 1, comm%Ncol
          pos = comm%ind(io,PP)
          if (pos.ne.0) then
            com = com + 1
            ANode = src(pos)-1
            BNode = dst(pos)-1
            DCsrc(com) = ANode
            DCdst(com) = BNode
            DCsic(com) = Gneights(ANode,BNode)
          endif
        enddo

        call de_alloc( comm%ind, name='comm%ind', routine='dscfComm' )
        call de_alloc( src, name='src', routine='dscfComm' )
        call de_alloc( dst, name='dst', routine='dscfComm' )
      else
C       Receive the communication scheduling computed by the master
#ifdef MPI
        call mpi_recv( DCsize, 1, MPI_integer, 0,  1,
     &                 MPI_Comm_world, Status, MPIerror )
#endif

        nullify( DCBuffer )
        call re_alloc( DCBuffer, 1, DCsize*3, name='DCBuffer',
     &                 routine='dscfComm' )
 
#ifdef MPI
        call mpi_recv( DCBuffer, DCsize*3, MPI_integer, 0,  1,
     &                 MPI_Comm_world, Status, MPIerror )
#endif
        DCsrc => DCBuffer(         1:  DCsize)
        DCdst => DCBuffer(  DCsize+1:2*DCsize)
        DCsic => DCBuffer(2*DCsize+1:3*DCsize)
      endif
      call de_alloc( Gneights, name='Gneights', routine='dscfComm' )

C     Compute the number of elements that we should copy from DSCF to
C     DSCFL, including the data that is in the current process, the data
C     that we should send to other processes and the data that we should
C     recieve from other processes.
      DCself  = xneights(Node+1)-xneights(Node)
      DCtotal = DCself
      do ii= 1, DCsize
        DCtotal = DCtotal + DCsic(ii)
      enddo

C     Create a permutation vector that allows us to know the original ID
C     of every orbital that we should send, receive or copy.
      nullify( DCinvp )
      call re_alloc( DCinvp, 1, DCtotal, name='DCinvp',
     &               routine='dscfComm' )
      io               = xneights(Node)
      DCinvp(1:DCself) = invp(io:io+DCself-1)
      nn               = DCself+1
#ifdef MPI
      do ii= 1, DCsize
        if (DCsrc(ii).eq.Node) then
          call mpi_recv( DCinvp(nn), DCsic(ii), MPI_integer,
     &                   DCdst(ii),  1, MPI_Comm_world, Status,
     &                   MPIerror )
        else
          io = xneights(DCsrc(ii))
          DCinvp(nn:nn+DCsic(ii)-1) = invp(io:io+DCsic(ii)-1)
          call MPI_Send( invp(io), DCsic(ii), MPI_integer,
     &                   DCsrc(ii), 1, MPI_Comm_World, MPIerror )
        endif
        nn = nn + DCsic(ii)
      enddo
#endif

      call de_alloc( xneights, name='xneights', routine='dscfComm' )
      call de_alloc( invp, name='invp', routine='dscfComm' )

      DCmax = 0
      do ii= 1, DCsize
        DCmax = max(DCmax,DCsic(ii))
      enddo
!--------------------------------------------------------------------------- END
      end subroutine dscfComm
      
      subroutine resetdscfComm( )
      use m_dscfComm, only : DCBuffer, DCinvp
      use alloc,      only : de_alloc
      implicit none
!------------------------------------------------------------------------- BEGIN
      if (associated(DCBuffer))
     &  call de_alloc( DCBuffer, name='DCBuffer', routine='dscfComm' )
      if (associated(DCinvp))
     &  call de_alloc( DCinvp, name='DCinvp', routine='dscfComm' )
!--------------------------------------------------------------------------- END
      end subroutine resetdscfComm
