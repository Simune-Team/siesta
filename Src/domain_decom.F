      module domain_decom
      use alloc,         only : re_alloc, de_alloc
      use parallel,      only : Node, Nodes
      use printMat
      use scheComm
      implicit none

      PUBLIC :: preSetOrbitLimits, domainDecom, gatherGlobalMatrix,
     &          use_dd, ulimit, llimit, use_dd_perm,
     &          dd_perm, dd_invp, dd_cperm, dd_nuo, dd_nnode,
     &          dd_ncolum, dd_ncom, lni, lnb, dd_bsiz, dd_comm

      PRIVATE
      logical          :: use_dd
      logical          :: use_dd_perm
      integer          :: llimit       ! Lower limit
      integer          :: ulimit       ! Upper limit
      integer          :: dd_nuo
      integer          :: dd_ncom
      integer          :: dd_ncolum
      integer, pointer :: dd_perm(:)
      integer, pointer :: dd_invp(:)
      integer, pointer :: dd_cperm(:)
      integer, pointer :: dd_comm(:)
      integer, pointer :: dd_bsiz(:)
      integer, pointer :: dd_nnode(:)

      integer          :: gni          ! Global Node internal
      integer          :: gnb          ! Global Node boundary
      integer          :: lni          ! Local Node internal
      integer          :: lnb          ! Local Node boundary
      integer, pointer :: listhPtrTot(:)
      integer, pointer :: listhTot(:)
      integer, pointer :: partition(:)
      integer, pointer :: ginde(:,:)
      integer, pointer :: gperm(:)
      integer, pointer :: ginvp(:)
      integer, pointer :: nnzInt(:)
      integer, pointer :: nnzBou(:)
      integer, pointer :: nNeigh(:)
      logical, pointer :: DomNeigh(:)
      type(COMM_T)     :: Gcomm

      contains
      subroutine preSetOrbitLimits( no )
      implicit none
C     Input variables
      integer      :: no
C     Local variables
      integer      :: blocks, remain
      
      blocks = no/Nodes
      remain = no - blocks*Nodes
      
      llimit = blocks*Node + 1 + min(remain,Node)
      ulimit = blocks*(Node+1) + 1 + min(remain,Node+1)
      dd_nuo = ulimit - llimit
      end subroutine preSetOrbitLimits


      subroutine domainDecom( nuotot, nuo, nnz )
#ifdef MPI
      use mpi_siesta
#endif
      implicit none
C     Input variables
      integer          :: nuotot, nuo, nnz

      call reduceGlobalMatrix( nuotot, nuo, nnz )

      if (node.eq.0) then
C       Set orbitals to partitions
        call makePartitions( nuotot, nodes )

        call boundaryNodes( nuotot )

        call checkPartitionSize( nuotot )

        call createPerm( nuotot )

        call compuNeigh( )

        call compuComm( )
      endif

      call distriGlobalMatrix( nuotot, nuo, nnz )

      if (node.eq.0) then
!       Free local memory
        call de_alloc( Gcomm%ind, 'comm%ind', 'domainDecom' )
        call de_alloc( partition, 'partition', 'domainDecom' )
        call de_alloc( gperm, 'gperm', 'domainDecom' )
!        call de_alloc( ginvp, 'ginvp', 'domainDecom' )
!        call de_alloc( nnzBou, 'nnzBou', 'domainDecom' )
!        call de_alloc( nnzInt, 'nnzInt', 'domainDecom' )
!        call de_alloc( ginde, 'ginde', 'domainDecom' )
        call de_alloc( listhTot, 'listhTot', 'domainDecom' )
        call de_alloc( listhPtrTot, 'listhPtrTot', 'domainDecom' )
        call de_alloc( DomNeigh, 'DomNeigh', 'domainDecom' )
      endif
      write(23,*) 'domainDecom END'
      call pxfflush(23)
      end subroutine domainDecom


      subroutine reduceGlobalMatrix( nuotot, nuo, nnz )
      use sparse_matrices, only : numh, listh
#ifdef MPI
      use mpi_siesta
#endif
      implicit none
C     Input variables
      integer,           intent(in) :: nuotot, nuo, nnz
C     Local variables
      integer                       :: blocks, remain, PP, io, ind,
     &                                 nnzG, nnzL
      integer,              pointer :: nuos(:)
#ifdef MPI
      integer                       :: MPIerror, Status(MPI_Status_Size)
#endif

      if (node.eq.0) then
        nullify(nuos)
        call re_alloc( nuos, 0, Nodes-1, 'nuos', 'reduceGlobalMatrix' )
        blocks = nuotot/Nodes
        remain = nuotot - blocks*Nodes
        do PP= 0, Nodes-1
          if (PP.lt.remain) then
            nuos(PP) = blocks + 1
          else
            nuos(PP) = blocks
          endif
        enddo

        nullify(listhPtrTot)
        call re_alloc( listhPtrTot, 0, nuotot, 'listhPtrTot',
     &                 'reduceGlobalMatrix' )
        listhPtrTot(1:nuo) = numh
        ind = nuo + 1
#ifdef MPI
        do PP= 1, Nodes-1
          call MPI_recv( listhPtrTot(ind), nuos(PP), MPI_INTEGER,
     &                   PP, 1, MPI_Comm_world, Status, MPIerror )
          ind = ind + nuos(PP)
        enddo
#endif
        listhPtrTot(0) = 1
        do io= 1, nuotot
          listhPtrTot(io) = listhPtrTot(io) + listhPtrTot(io-1)
        enddo
        nnzG = listhPtrTot(nuotot) - 1

        nullify(listhTot)
        call re_alloc( listhTot, 1, nnzG, 'listhTot',
     &                 'reduceGlobalMatrix' )

        listhTot(1:nnz) = listh
        ind = nuo
#ifdef MPI
        do PP= 1, Nodes-1
          nnzL = listhPtrTot(ind+nuos(PP)) - listhPtrTot(ind)
          call MPI_recv( listhTot(listhPtrTot(ind)), nnzL, MPI_INTEGER,
     &                   PP, 1, MPI_Comm_world, Status, MPIerror )
          ind = ind + nuos(PP)
        enddo
#endif
        call printMatrix( nuotot, nuotot, listhTot, listhPtrTot,
     &                    'globmat1.ps', 'Global_matrix' )

        call de_alloc( nuos, 'nuos', 'domainDecom' )
      else
#ifdef MPI
        call MPI_send( numh, nuo, MPI_INTEGER, 0, 1,
     &                 MPI_Comm_world, MPIerror )
        call MPI_send( listh, nnz, MPI_INTEGER, 0, 1,
     &                 MPI_Comm_world, MPIerror )
#endif
      endif
      end subroutine reduceGlobalMatrix



      subroutine makePartitions( no, nparts )
      use sys,           only : die
      implicit none
!     Input/Output variables
      integer,           intent(in) :: no, nparts
!     Local variables
      integer                       :: ii, jj, ind, nnz, dummy(1),
     &                                 metis_options(1), edgecut
      integer,              pointer :: iia(:), jja(:)
      logical                       :: notfound

      nullify(partition)
      call re_alloc( partition, 1, no, 'partition', 'makePartitions' )

!     Drop diagonal nodes (autoconnection)
      nnz = listhPtrTot(no) - listhPtrTot(0) - no
      nullify(jja,iia)
      call re_alloc( jja, 0, no, 'jja', 'makePartitions' )
      call re_alloc( iia, 1, nnz, 'iia', 'makePartitions' )

      jja(0) = 1
      ind    = 1
      do jj= 1, no
        notfound = .true.
        DO ii= listhPtrTot(jj-1), listhPtrTot(jj)-1
          if (listhTot(ii).eq.jj) then
            notfound = .false.
          else
            iia(ind) = listhTot(ii)
            ind = ind + 1
          endif
        ENDDO
        if (notfound) call die( 'not diagonal element in this row' )
        jja(jj) = ind
      enddo

      call printmatrix( no, no, iia, jja,
     &                  'globmat2.ps', 'Global_matrix' )

      metis_options(1) = 0
      if (nparts .gt. 8 ) then
        call METIS_PartGraphKway( no, jja, iia, dummy, dummy, 0,
     &                            1, nparts, metis_options,
     &                            edgecut, partition )
      else
        call METIS_PartGraphRecursive( no, jja, iia, dummy, dummy, 0,
     &                                 1, nparts, metis_options,
     &                                 edgecut, partition )
      endif
      call de_alloc( jja, 'jja', 'makePartitions' )
      call de_alloc( iia, 'iia', 'makePartitions' )
      end subroutine makePartitions


      subroutine checkPartitionSize( no )
!     Input/Output variables
      integer,           intent(in) :: no
!     Local variables
      integer                       :: io, PP

      nullify(nnzInt,nnzBou)
      call re_alloc( nnzInt, 1, Nodes, 'nnzInt', 'checkPartitionSize' )
      call re_alloc( nnzBou, 1, Nodes, 'nnzBou', 'checkPartitionSize' )
      nnzInt = 0
      nnzBou = 0
      do io= 1, no
        PP = partition(io)
        if (PP.gt.0) then
          nnzInt(PP) = nnzInt(PP)+listhPtrTot(io)-listhPtrTot(io-1)
        else
          PP = -PP
          nnzBou(PP) = nnzBou(PP)+listhPtrTot(io)-listhPtrTot(io-1)
        endif
      enddo
      end subroutine checkPartitionSize


      subroutine boundaryNodes( no )
      implicit none
!     Input/Output variables
      integer,           intent(in) :: no
!     Local variables
      logical                       :: notfound
      integer                       :: ii, jj

      nullify(ginde)
      call re_alloc( ginde, 1, 4, 1, Nodes, 'ginde', 'boundaryNodes' )
      ginde = 0

      do jj = 1, no
!       look if all niegh. are in the same partition
        notfound = .true.
        ii = listhPtrTot(jj-1)
        do while (ii.lt.listhPtrTot(jj) .and. notfound)
          if (partition(jj).ne.abs(partition(listhTot(ii))))
     &      notfound = .false.
          ii = ii + 1
        enddo
        if ( .not.notfound ) then
          ginde(4,partition(jj)) = ginde(4,partition(jj)) + 1
          partition(jj) = -partition(jj)
        else
          ginde(2,partition(jj)) = ginde(2,partition(jj)) + 1
        endif
      enddo

      ginde(1,1) = 1
      do jj= 2, Nodes
        ginde(1,jj) = ginde(2,jj-1) + ginde(1,jj-1)
      enddo
      gni = ginde(1,Nodes) + ginde(2,Nodes) - 1

      ginde(3,1) = gni + 1
      do jj= 2, Nodes
        ginde(3,jj) = ginde(4,jj-1) + ginde(3,jj-1)
      enddo
      gnb = no - gni
      end subroutine boundaryNodes 

      subroutine createPerm( no )
      implicit none
!     Input/Output variables
      integer,           intent(in) :: no
!     Local variables
      integer               :: dom, maxni, maxnb, maxnz, noloc, boloc,
     &                         metisOpt(1), ii, jj, kk, offsetI,
     &                         offsetB
      integer,      pointer :: adj(:), xadj(:), perI(:), invI(:),
     &                         perB(:), invB(:), perR(:), invR(:)

      nullify(gperm,ginvp)
      call re_alloc( gperm, 1, no, 'gperm', 'createPerm' )
      call re_alloc( ginvp, 1, no, 'ginvp', 'createPerm' )

      maxni = 1
      do dom= 1, Nodes
        maxni = max(maxni,ginde(2,dom))
      enddo
      maxnb = 1
      do dom= 1, Nodes
        maxnb = max(maxnb,ginde(4,dom))
      enddo
      maxnz = 1
      do dom= 1, Nodes
        maxnz = max(maxnz,nnzInt(dom))
      enddo
      nullify(adj,xadj,perI,invI,perR,invR,perB,invB)
      call re_alloc(  adj, 1, maxnz,  'adj', 'createPerm' )
      call re_alloc( xadj, 0, maxni, 'xadj', 'createPerm' )
      call re_alloc( perI, 1,    no, 'perI', 'createPerm' )
      call re_alloc( invI, 1, maxni, 'invI', 'createPerm' )
      call re_alloc( perR, 1, maxni, 'perR', 'createPerm' )
      call re_alloc( invR, 1, maxni, 'invR', 'createPerm' )
      call re_alloc( perB, 1,    no, 'perB', 'createPerm' )
      call re_alloc( invB, 1, maxnb, 'invB', 'createPerm' )

      metisOpt = 0
      offsetI  = 0
      offsetB  = gni
      do dom= 1, Nodes
!       Loop over each subdomain: separate interior and boundary nodes
        call createSubgraph( no, listhTot, listhPtrTot, dom,
     &                       adj, xadj, perI, invI, perB, invB )
!       Number interior nodes
        noloc = ginde(2,dom)
        if (noloc.gt.0) then
          call METIS_NodeND( noloc, xadj, adj, 1, metisOpt, perR, invR )

          do ii = 1, noloc
            kk        = invI(ii)
            jj        = perR(ii) + offsetI
            gperm(kk) = jj
            ginvp(jj) = kk
          end do
          offsetI = offsetI + noloc 
        endif
!       Number boundary nodes
        boloc = ginde(4,dom)
        do ii = 1, boloc
          kk       = invB(ii)
          jj       = ii + offsetB
          gperm(kk) = jj
          ginvp(jj) = kk
        end do
        offsetB = offsetB + boloc
      enddo
      call de_alloc(  adj,  'adj', 'createPerm' )
      call de_alloc( xadj, 'xadj', 'createPerm' )
      call de_alloc( perI, 'perI', 'createPerm' )
      call de_alloc( invI, 'invI', 'createPerm' )
      call de_alloc( perR, 'perR', 'createPerm' )
      call de_alloc( invR, 'invR', 'createPerm' )
      call de_alloc( perB, 'perB', 'createPerm' )
      call de_alloc( invB, 'invB', 'createPerm' )


      call printPermMatrix( no, listhTot, listhPtrTot, gperm, ginvp,
     &                      'globmat3.ps', 'Global matrix' )

      end subroutine createPerm

      subroutine compuNeigh( )
      implicit none
!     Local variables
      integer                       :: tsize, io, jo, ko, lo, PP, QQ,
     &                                 ini, fin, ind

      tsize = (Nodes-1)*Nodes/2 
      nullify(DomNeigh)
      call re_alloc( DomNeigh, 1, tsize, 'DomNeigh', 'compuNeigh' )

      DomNeigh = .false.

      do PP=1, Nodes
        ini = ginde(3,PP)
        fin = ginde(3,PP) + ginde(4,PP) - 1

        do io= ini, fin
          jo = ginvp(io)
          do ko= listhPtrTot(jo-1), listhPtrTot(jo)-1
            lo = listhTot(ko)
            QQ = abs(partition(lo))
            if (PP.lt.QQ) then
              ind = PP + ((QQ-1)*(QQ-2))/2
              DomNeigh(ind) = .true.
            else if (PP.gt.QQ) then
              ind = QQ + ((PP-1)*(PP-2))/2
              DomNeigh(ind) = .true.
            endif
          enddo
        enddo
      enddo
      end subroutine compuNeigh

      subroutine compuComm( )
      implicit none
!     Local variables
      integer          :: ind, ncom, PP, QQ, KK
      integer, pointer :: src(:), dst(:)

      ncom   = 0
      ind    = 0
      do PP=2, Nodes
        do QQ=1, PP-1
          ind = ind + 1
          if (DomNeigh(ind)) ncom = ncom + 1
        enddo
      enddo

      nullify(src,dst)
      call re_alloc( src, 1, ncom, 'src', 'compuComm' )
      call re_alloc( dst, 1, ncom, 'dst', 'compuComm' )
      ncom = 0
      ind  = 0
      do PP=2, Nodes
        do QQ=1, PP-1
          ind = ind + 1
          if (DomNeigh(ind)) then
            ncom = ncom + 1
            src(ncom) = QQ
            dst(ncom) = PP
          endif
        enddo
      enddo

      Gcomm%np = Nodes
C     reschedule the communications in order to minimize the time
      call scheduleComm( ncom, src, dst, Gcomm )

      nullify(nNeigh)
      call re_alloc( nNeigh, 1, Nodes, 'nNeigh', 'compuComm' )
      do PP=1, Nodes
        KK = 0
        do QQ= 1, Gcomm%ncol
          ind = Gcomm%ind(QQ,PP)
          if (ind.ne.0) then
            KK = KK + 1
            if (src(ind).eq.PP) then
              Gcomm%ind(KK,PP) = dst(ind)
            else
              Gcomm%ind(KK,PP) = src(ind)
            endif
          endif
        enddo
        nNeigh(PP) = KK
      enddo

      call de_alloc( dst, 'src', 'compuComm' )
      call de_alloc( src, 'dst', 'compuComm' )
      end subroutine compuComm


      subroutine createSubgraph( no, ia, ja, dom, adj, xadj,
     &                           permI, invI, permB, invB )
      implicit none
!     Input Variables
      integer, intent(in)  :: no, dom
      integer, intent(in)  :: ia(*), ja(0:no)
!     Output Variables
      integer, intent(out) :: adj(*), xadj(0:*), permI(*), invI(*),
     &                        permB(*), invB(*)
!     Local Variables
      integer              :: next, vv, ww, ii, nodeI, nodeB

      nodeI   = 0
      nodeB   = 0
      xadj(0) = 1
      next    = 1
      do vv = 1, no
        if (partition(vv) == dom) then
          nodeI       = nodeI + 1
          permI(vv)   = nodeI
          invI(nodeI) = vv
          do ii = ja(vv-1), ja(vv)-1
            ww = ia(ii)
            if (ww/=vv) then
              if (partition(ww) == dom) then
                adj(next) = ww
                next = next + 1
              endif
            endif
          enddo
          xadj(nodeI) = next

        else if (partition(vv) == -dom) then
          nodeB       = nodeB + 1
          permB(vv)   = nodeB
          invB(nodeB) = vv
        endif
      enddo

      do ii = 1, next-1
        adj(ii) = permI(adj(ii))
      enddo

      end subroutine createSubgraph


      subroutine distriGlobalMatrix( notot, no, nnz )
      use sparse_matrices, only : numh, listhptr, listh
#ifdef MPI
      use mpi_siesta
#endif
      implicit none
!     Input Variables
      integer,  intent(in) :: notot
      integer, intent(out) :: no, nnz
!     Local Variables
      integer              :: PP, QQ, RR, l_no, l_lni, l_lnb, l_nnz,
     &                        l_nei, l_bou, io, jo, ll, ind, BUFF(5)
      integer,     pointer :: linvp(:), lbinvp(:), lbsize(:), iia(:),
     &                        jja(:)
#ifdef MPI
      integer              :: MPIerror, Status(MPI_Status_Size)
#endif

      nullify(dd_nnode)
      call re_alloc( dd_nnode, 1, notot, 'dd_nnode',
     &               'distriGlobalMatrix' )
      if (Node.eq.0) dd_nnode = abs(partition) - 1
#ifdef MPI
      call MPI_BCast( dd_nnode, notot, MPI_INTEGER, 0,
     &                MPI_Comm_world, MPIerror )
#endif
      if (Node.eq.0) then
        do PP=2, Nodes
          l_lni = ginde(2,PP)
          l_lnb = ginde(4,PP)
          l_no  = l_lni + l_lnb
          l_nnz = nnzInt(PP) + nnzBou(PP)
          l_nei = nNeigh(PP)
          l_bou = 0
          do QQ= 1, l_nei
            l_bou = l_bou + ginde(4,Gcomm%ind(QQ,PP))
          enddo

          nullify(linvp)
          call re_alloc( linvp, 1, l_no, 'linvp', 'distriGlobalMatrix' )
          linvp(1:ginde(2,PP)) =
     &      ginvp(ginde(1,PP):ginde(1,PP)+ginde(2,PP)-1)
          linvp(ginde(2,PP)+1:l_no) =
     &      ginvp(ginde(3,PP):ginde(3,PP)+ginde(4,PP)-1)

          nullify(lbinvp)
          call re_alloc( lbinvp, 1, l_bou, 'lbinvp',
     &                   'distriGlobalMatrix' )
          ind = 1
          do QQ= 1, l_nei
            RR = Gcomm%ind(QQ,PP)
            ll = ginde(4,RR)
            lbinvp(ind:ind+ll-1) =
     &        ginvp(ginde(3,RR):ginde(3,RR)+ll-1)
            ind = ind + ll
          enddo

          nullify(lbsize)
          call re_alloc( lbsize, 1, l_nei, 'lbsize',
     &                   'distriGlobalMatrix' )
          do QQ= 1, l_nei
            lbsize(QQ) = ginde(4,Gcomm%ind(QQ,PP))
          enddo

          nullify(jja,iia)
          call re_alloc( jja, 1, l_no, 'jja', 'distriGlobalMatrix' )
          call re_alloc( iia, 1, l_nnz, 'iia', 'distriGlobalMatrix' )

          ind = 1
          DO io= 1, l_no
            jo                = linvp(io)
            ll                = listhPtrTot(jo) - listhPtrTot(jo-1)
            jja(io)           = ll
            iia(ind:ind+ll-1) = listhTot(listhPtrTot(jo-1):
     &                                   listhPtrTot(jo)-1)
            ind               = ind + ll
          ENDDO

          BUFF(1) = l_lni
          BUFF(2) = l_lnb
          BUFF(3) = l_nnz
          BUFF(4) = l_nei
          BUFF(5) = l_bou

#ifdef MPI
          call MPI_Send( BUFF, 5, MPI_INTEGER, PP-1, 0,
     &                   MPI_Comm_world, MPIerror )

          call MPI_Send( linvp, l_no, MPI_INTEGER, PP-1, 0,
     &                   MPI_Comm_world, MPIerror )

          call MPI_Send( lbinvp, l_bou, MPI_INTEGER, PP-1, 0,
     &                   MPI_Comm_world, MPIerror )

          call MPI_Send( Gcomm%ind(1,PP), l_nei, MPI_INTEGER, PP-1, 0,
     &                   MPI_Comm_world, MPIerror )

          call MPI_Send( lbsize, l_nei, MPI_INTEGER, PP-1, 0,
     &                   MPI_Comm_world, MPIerror )

          call MPI_Send( jja, l_no, MPI_INTEGER, PP-1, 0,
     &                   MPI_Comm_world, MPIerror )

          call MPI_Send( iia, l_nnz, MPI_INTEGER, PP-1, 0,
     &                   MPI_Comm_world, MPIerror )
#endif

          call de_alloc( jja, 'jja', 'distriGlobalMatrix' )
          call de_alloc( iia, 'iia', 'distriGlobalMatrix' )
          call de_alloc( lbsize, 'lbsize', 'distriGlobalMatrix' )
          call de_alloc( lbinvp, 'lbinvp', 'distriGlobalMatrix' )
          call de_alloc( linvp, 'linvp', 'distriGlobalMatrix' )
        enddo
C       Compute data for local Process
        l_lni = ginde(2,1)
        l_lnb = ginde(4,1)
        l_no  = l_lni + l_lnb
        l_nnz = nnzInt(1) + nnzBou(1)
        l_nei = nNeigh(1)
        l_bou = 0
        do QQ= 1, l_nei
          l_bou = l_bou + ginde(4,Gcomm%ind(QQ,1))
        enddo

        nullify(dd_invp)
        call re_alloc( dd_invp, 1, l_no, 'dd_invp',
     &                 'distriGlobalMatrix' )
        dd_invp(1:ginde(2,1)) =
     &    ginvp(ginde(1,1):ginde(1,1)+ginde(2,1)-1)
        dd_invp(ginde(2,1)+1:l_no) =
     &    ginvp(ginde(3,1):ginde(3,1)+ginde(4,1)-1)

        nullify(lbinvp)
        call re_alloc( lbinvp, 1, l_bou, 'lbinvp',
     &                 'distriGlobalMatrix' )
        ind = 1
        do QQ= 1, l_nei
          RR = Gcomm%ind(QQ,1)
          ll = ginde(4,RR)
          lbinvp(ind:ind+ll-1) =
     &      ginvp(ginde(3,RR):ginde(3,RR)+ll-1)
          ind = ind + ll
        enddo

        nullify(dd_comm)
        call re_alloc( dd_comm, 1, l_nei, 'dd_comm',
     &                 'distriGlobalMatrix' )
        dd_comm = Gcomm%ind(1:l_nei,1)

        nullify(dd_bsiz)
        call re_alloc( dd_bsiz, 1, l_nei, 'dd_bsiz',
     &                 'distriGlobalMatrix' )
        do QQ= 1, l_nei
          dd_bsiz(QQ) = ginde(4,Gcomm%ind(QQ,1))
        enddo

C        call de_alloc( numh, 'numh', 'distriGlobalMatrix' )
        call re_alloc( numh, 1, l_no, 'numh', 'distriGlobalMatrix' )
C        call de_alloc( listh, 'listh', 'distriGlobalMatrix' )
        call re_alloc( listh, 1, l_nnz, 'listh', 'distriGlobalMatrix' )

        ind = 1
        DO io= 1, l_no
          jo                  = dd_invp(io)
          ll                  = listhPtrTot(jo) - listhPtrTot(jo-1)
          numh(io)            = ll
          listh(ind:ind+ll-1) = listhTot(listhPtrTot(jo-1):
     &                                 listhPtrTot(jo)-1)
          ind                 = ind + ll
        ENDDO
      else
C       Receive Data from process 0
#ifdef MPI
        call MPI_recv( BUFF, 5, MPI_INTEGER, 0, 0, MPI_Comm_world,
     &                 Status, MPIerror )
#endif

        l_lni = BUFF(1)
        l_lnb = BUFF(2)
        l_nnz = BUFF(3)
        l_nei = BUFF(4)
        l_bou = BUFF(5)
        l_no  = l_lni + l_lnb

        nullify(dd_invp)
        call re_alloc( dd_invp, 1, l_no, 'dd_invp',
     &                 'distriGlobalMatrix' )
#ifdef MPI
        call MPI_Recv( dd_invp, l_no, MPI_INTEGER, 0, 0,
     &                 MPI_Comm_world, Status, MPIerror )
#endif

        nullify(lbinvp)
        call re_alloc( lbinvp, 1, l_bou, 'lbinvp',
     &                 'distriGlobalMatrix' )
#ifdef MPI
        call MPI_Recv( lbinvp, l_bou, MPI_INTEGER, 0, 0,
     &                 MPI_Comm_world, Status, MPIerror )
#endif

        nullify(dd_comm)
        call re_alloc( dd_comm, 1, l_nei, 'dd_comm',
     &                 'distriGlobalMatrix' )
#ifdef MPI
        call MPI_Recv( dd_comm, l_nei, MPI_INTEGER, 0, 0,
     &                 MPI_Comm_world, Status, MPIerror )
#endif

        nullify(dd_bsiz)
        call re_alloc( dd_bsiz, 1, l_nei, 'dd_bsiz',
     &                 'distriGlobalMatrix' )
#ifdef MPI
        call MPI_Recv( dd_bsiz, l_nei, MPI_INTEGER, 0, 0,
     &                 MPI_Comm_world, Status, MPIerror )
#endif

C        call de_alloc( numh, 'numh', 'distriGlobalMatrix' )
        call re_alloc( numh, 1, l_no, 'numh', 'distriGlobalMatrix' )
C        call de_alloc( listh, 'listh', 'distriGlobalMatrix' )
        call re_alloc( listh, 1, l_nnz, 'listh', 'distriGlobalMatrix' )

#ifdef MPI
        call MPI_Recv( numh, l_no, MPI_INTEGER, 0, 0,
     &                 MPI_Comm_world, Status, MPIerror )

        call MPI_Recv( listh, l_nnz, MPI_INTEGER, 0, 0,
     &                 MPI_Comm_world, Status, MPIerror )
#endif
      endif

      no        = l_no
      nnz       = l_nnz
      lni       = l_lni
      lnb       = l_lnb
      dd_ncom   = l_nei
      dd_ncolum = l_no + l_bou
      write(23,*) 'l_no     ', l_no
      write(23,*) 'l_bou    ', l_bou
      write(23,*) 'dd_ncolum', dd_ncolum
      call pxfflush(23)
C      call de_alloc( listhptr, 'listhptr', 'distriGlobalMatrix' )
      call re_alloc( listhptr, 1, no, 'listhptr', 'distriGlobalMatrix' )
      ind = 0
      do io= 1, no
        listhptr(io) = ind
        ind          = ind + numh(io)
      enddo

      nullify(dd_perm)
      call re_alloc( dd_perm, 1, notot, 'dd_perm','distriGlobalMatrix' )
      dd_perm = 0
      do io= 1, no
        ind          = dd_invp(io)
        dd_perm(ind) = io
      enddo


      nullify(dd_cperm)
      call re_alloc( dd_cperm, 1, notot, 'dd_cperm',
     &               'distriGlobalMatrix' )
      dd_cperm = dd_perm
      do io= 1, l_bou
        ind           = lbinvp(io)
        dd_cperm(ind) = io + l_no
      enddo

      call de_alloc( lbinvp, 'lbinvp', 'distriGlobalMatrix' )


      use_dd_perm = .true.
      dd_nuo      = no
      write(23,*) 'distriGlobalMatrix END'
      call pxfflush(23)
      end subroutine distriGlobalMatrix

      subroutine gatherGlobalMatrix( nuotot, nuo, nnz, ML, filename )
#ifdef MPI
      use mpi_siesta
#endif
      use precision
      use sparse_matrices, only : numh, listh
      implicit none
!     Input Variables
      integer,  intent(in) :: nuotot, nuo, nnz
      real(dp),     target :: ML(:)
      character(*)         :: filename
!     Local Variables
      integer              :: PP, nuo_l, io, jo, ko, nnzG, nnz_l,
     &                        ll, ind, fd
      integer,     pointer :: ja(:), ia(:)
      real(dp),    pointer :: M(:), an(:)
#ifdef MPI
      integer              :: MPIerror, Status(MPI_Status_Size)
#endif
!
!     Need arrays: ginde, ginvp
!
      if (node.eq.0) then
        nullify(listhPtrTot)
        call re_alloc( listhPtrTot, 0, nuotot, 'listhPtrTot',
     &                 'gatherGlobalMatrix' )

        do PP=1, Nodes
          nuo_l = ginde(2,PP) + ginde(4,PP)
          if (PP.eq.1) then
            ja => numh
          else
            nullify(ja)
            call re_alloc( ja, 1, nuo_l, 'ja', 'gatherGlobalMatrix' )
#ifdef MPI
            call MPI_recv( ja, nuo_l, MPI_INTEGER,
     &                     PP-1, 1, MPI_Comm_world, Status, MPIerror )
#endif
          endif
          jo = ginde(1,PP)
          do io= 1, ginde(2,PP)
            ko              = ginvp(jo)
            listhPtrTot(ko) = ja(io)
            jo              = jo + 1
          enddo
          jo = ginde(3,PP)
          do io= io, nuo_l
            ko              = ginvp(jo)
            listhPtrTot(ko) = ja(io)
            jo              = jo + 1
          enddo
          if (PP.eq.1) then
            ja => null()
          else
            call de_alloc( ja, 'ja', 'gatherGlobalMatrix' )
          endif
        enddo
        listhPtrTot(0) = 1
        do io= 1, nuotot
          listhPtrTot(io) = listhPtrTot(io) + listhPtrTot(io-1)
        enddo
        nnzG = listhPtrTot(nuotot)-1

        nullify(listhTot,M)
        call re_alloc( listhTot, 1, nnzG, 'listhTot',
     &                 'gatherGlobalMatrix' )
        call re_alloc( M, 1, nnzG, 'M', 'gatherGlobalMatrix' )
        do PP=1, Nodes
          nuo_l = ginde(2,PP) + ginde(4,PP)
          nnz_l = nnzInt(PP) + nnzBou(PP)
          if (PP.eq.1) then
            ia => listh
            an => ML
          else
            nullify(ia,an)
            call re_alloc( ia, 1, nnz_l, 'ia', 'gatherGlobalMatrix' )
            call re_alloc( an, 1, nnz_l, 'an', 'gatherGlobalMatrix' )
#ifdef MPI
            call MPI_recv( ia, nnz_l, MPI_INTEGER,
     &                     PP-1, 1, MPI_Comm_world, Status, MPIerror )
            call MPI_recv( an, nnz_l, MPI_DOUBLE_PRECISION,
     &                     PP-1, 1, MPI_Comm_world, Status, MPIerror )
#endif
          endif
          jo   = ginde(1,PP)
          ind = 1
          do io= 1, ginde(2,PP)
            ko  = ginvp(jo)
            ll  = listhPtrTot(ko)-listhPtrTot(ko-1)
            listhTot(listhPtrTot(ko-1):listhPtrTot(ko)-1) =
     &        ia(ind:ind+ll-1)
            M(listhPtrTot(ko-1):listhPtrTot(ko)-1) =
     &        an(ind:ind+ll-1)
            jo  = jo + 1
            ind = ind + ll
          enddo
          jo = ginde(3,PP)
          do io= io, nuo_l
            ko  = ginvp(jo)
            ll  = listhPtrTot(ko)-listhPtrTot(ko-1)
            listhTot(listhPtrTot(ko-1):listhPtrTot(ko)-1) =
     &        ia(ind:ind+ll-1)
            M(listhPtrTot(ko-1):listhPtrTot(ko)-1) =
     &        an(ind:ind+ll-1)
            jo  = jo + 1
            ind = ind + ll
          enddo
          if (PP.eq.1) then
            ia => null()
            an => null()
          else
            call de_alloc( ia, 'ia', 'gatherGlobalMatrix' )
            call de_alloc( an, 'an', 'gatherGlobalMatrix' )
          endif
        enddo

        call io_assign( fd )
        open( unit=fd, file=filename, status='unknown' )
        write(fd,*) nuotot, nuotot, nnzG
        write(fd,*) listhPtrTot
        write(fd,*) listhTot
        write(fd,*) M

        do io= 1, nuotot
          do jo= listhPtrTot(io-1), listhPtrTot(io)-1
            write(fd,*) io, listhTot(jo), M(jo)
          enddo
        enddo
        call pxfflush(fd)
        call io_close( fd )

        call de_alloc( ja, 'ja', 'gatherGlobalMatrix' )
        call de_alloc( an, 'an', 'gatherGlobalMatrix' )
        call de_alloc( M, 'M', 'gatherGlobalMatrix' )
        call de_alloc( listhPtrTot, 'listhPtrTot','gatherGlobalMatrix' )
        call de_alloc( listhTot, 'listhTot', 'gatherGlobalMatrix' )
      else
#ifdef MPI
        call MPI_send( numh, nuo, MPI_INTEGER, 0, 1,
     &                 MPI_Comm_world, MPIerror )
        call MPI_send( listh, nnz, MPI_INTEGER, 0, 1,
     &                 MPI_Comm_world, MPIerror )
        call MPI_send( ML, nnz, MPI_DOUBLE_PRECISION, 0, 1,
     &                 MPI_Comm_world, MPIerror )
#endif
      endif
      end subroutine gatherGlobalMatrix

      end module domain_decom
