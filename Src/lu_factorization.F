      MODULE LU_FACTORIZATION
      use alloc,         only : re_alloc, de_alloc
      use superlu
      use precision
      implicit none

      public :: compCommDD2LU, distCommDD2LU, use_lu_fact, factOverlap,
     &          LU, reduce_array_lu

      private
      integer           :: lu_nn
      integer           :: lu_nnz
      integer           :: lu_ncom
      integer,  pointer :: lu_comm(:,:)
      integer,  pointer :: LUnNeigh(:)
      integer,  pointer :: lu_ia(:)
      integer,  pointer :: lu_ja(:)
      real(dp), pointer :: lu_an(:)
      logical           :: use_lu_fact
      type(superlu_t)   :: LU

      contains


      subroutine compCommDD2LU( NODES, N, GINDE )
!     In order to compute the LU factorization of the Overlap matrix
!     we need a new data distribution. This function computes
!     the communications needed to move data from one data distribution
!     to the other (Domain Decomposition/LU factorization).
!     This is only computed by process 0.
!
!     LUnNeigh(1:NODES)   Number of communications of every process
!     lu_comm(1,:)        Source process (DD data distribution)
!     lu_comm(2,:)        Destination process (LU data distribution)
!     lu_comm(3,:)        Size of the communication
!     lu_comm(4,:)        Index of the first element in Source process
!     lu_comm(5,:)        Index of the first element in Destination process
      use scheComm
      implicit none
!     Input Variables
      integer,        intent(in) :: NODES, N, GINDE(4,NODES)
!     Local Variables
      integer                    :: ndiv, nmod, PP, P1, P2, ncom, ind,
     &                              gni
      logical                    :: notend
      integer,           pointer :: dd_limits(:), lu_limits(:), src(:),
     &                              dst(:), ini(:), fin(:)
      type(COMM_T)               :: LUcomm

      nullify(lu_limits,dd_limits)
      call re_alloc( dd_limits, 1, 2*Nodes+1, 'dd_limits', 'LU_FAC' )
      call re_alloc( lu_limits, 1,   Nodes+1, 'lu_limits', 'LU_FAC' )

      ndiv = n/nodes
      nmod = mod(n,nodes)

!     Compute domain limits for "Domain Decomposition"
      dd_limits(1) = 1
      do PP= 1, nodes
        dd_limits(PP+1) = dd_limits(PP) + ginde(2,PP)
      enddo
      P1 = nodes+1
      do PP= 1, nodes
        dd_limits(P1+1) = dd_limits(P1) + ginde(4,PP)
        P1 = P1 + 1
      enddo

!     Compute domain limits for "LU factorization"
      lu_limits(1) = 1
      do PP= 1, nodes
        lu_limits(PP+1) = lu_limits(PP) + ndiv
        if (PP.le.nmod) lu_limits(PP+1) = lu_limits(PP+1) + 1
      enddo

!     Count the number of intersections between LU data distribution
!     and Domain Decomposition data distribution
      P1     = 1
      P2     = 1
      ncom   = 0
      notend = .true.
      do while( notend )
        do while( dd_limits(P1).eq.dd_limits(P1+1) .and.
     &            dd_limits(P1+1).lt.n+1 )
          P1   = P1+1
        enddo
        do while( lu_limits(P2).eq.lu_limits(P2+1) .and.
     &            lu_limits(P2+1).lt.n+1 )
          P2   = P2+1
        enddo
        ncom = ncom + 1
        if (dd_limits(P1+1).eq.n+1 .and. lu_limits(P2+1).eq.n+1) then
          notend = .false.
        else
          if (dd_limits(P1+1).lt.lu_limits(P2+1)) then
            P1 = P1+1
          else if (dd_limits(P1+1).eq.lu_limits(P2+1)) then
            P1 = P1+1
            P2 = P2+1
          else
            P2 = P2+1
          endif
        endif          
      enddo

!     Make a list of intersections between LU data distribution
!     and Domain Decomposition data distribution
      nullify(src,dst,ini,fin)
      call re_alloc( src, 1, ncom, 'src', 'LU_FAC' )
      call re_alloc( dst, 1, ncom, 'dst', 'LU_FAC' )
      call re_alloc( ini, 1, ncom, 'ini', 'LU_FAC' )
      call re_alloc( fin, 1, ncom, 'fin', 'LU_FAC' )

      P1     = 1
      P2     = 1
      ncom   = 0
      notend = .true.
      do while( notend )
        do while( dd_limits(P1).eq.dd_limits(P1+1) .and.
     &            dd_limits(P1+1).lt.n+1 )
          P1   = P1+1
        enddo
        do while( lu_limits(P2).eq.lu_limits(P2+1) .and.
     &            lu_limits(P2+1).lt.n+1 )
          P2   = P2+1
        enddo

        ncom      = ncom + 1
        src(ncom) = mod(P1-1,nodes)+1
        dst(ncom) = P2
        ini(ncom) = max(dd_limits(P1  ),lu_limits(P2  ))
        fin(ncom) = min(dd_limits(P1+1),lu_limits(P2+1))

        if (dd_limits(P1+1).eq.n+1 .and. lu_limits(P2+1).eq.n+1) then
          notend = .false.
        else
          if (dd_limits(P1+1).lt.lu_limits(P2+1)) then
            P1 = P1+1
          else if (dd_limits(P1+1).eq.lu_limits(P2+1)) then
            P1 = P1+1
            P2 = P2+1
          else
            P2 = P2+1
          endif
        endif
      enddo

      LUcomm%np = Nodes
C     reschedule the communications in order to minimize the time
      call scheduleComm( ncom, src, dst, LUcomm )

      nullify(LUnNeigh)
      call re_alloc( LUnNeigh, 1,Nodes, 'LUnNeigh', 'LU_FAC' )
      ncom = 0
      do P1=1, Nodes
        PP = 0
        do P2= 1, LUcomm%ncol
          ind = LUcomm%ind(P2,P1)
          if (ind.ne.0) PP = PP + 1
        enddo
        LUnNeigh(P1) = PP
        ncom         = ncom + PP
      enddo

      nullify(lu_comm)
      call re_alloc( lu_comm, 1, 5, 1, ncom, 'lu_comm', 'LU_FAC' )
      ncom = 0
      gni  = dd_limits(NODES+1)-1
      do P1=1, Nodes
        do P2= 1, LUcomm%ncol
          ind = LUcomm%ind(P2,P1)
          if (ind.ne.0) then
            ncom = ncom + 1
            lu_comm(1,ncom) = src(ind) - 1
            lu_comm(2,ncom) = dst(ind) - 1
            lu_comm(3,ncom) = fin(ind) - ini(ind)
            if (src(ind).eq.P1) then
              if (ini(ind).lt.gni) then
                lu_comm(4,ncom) = ini(ind) - dd_limits(P1) + 1
              else
                lu_comm(4,ncom) = ginde(2,P1) +
     &                            ini(ind) - dd_limits(P1+Nodes) + 1
              endif
            else
              lu_comm(4,ncom) = -1
            endif

            if (dst(ind).eq.P1) then
              lu_comm(5,ncom) = ini(ind) - lu_limits(P1) + 1
            else
              lu_comm(5,ncom) = -1
            endif
          endif
        enddo
      enddo

      call de_alloc( LUcomm%ind, 'comm%ind', 'scheduleComm' )

      call de_alloc( fin, 'fin', 'LU_FAC' )
      call de_alloc( ini, 'ini', 'LU_FAC' )
      call de_alloc( dst, 'dst', 'LU_FAC' )
      call de_alloc( src, 'src', 'LU_FAC' )

      call de_alloc( lu_limits, 'lu_limits', 'LU_FAC' )
      call de_alloc( dd_limits, 'dd_limits', 'LU_FAC' )

      end subroutine compCommDD2LU

      subroutine distCommDD2LU( NODE, NODES, N )
!     Process 0 distributes the information computed in function
!     compCommDD2LU to all the other processes.
!
!     lu_ncom    Number of communications of current process
!     lu_comm    Information needed of every communication
!     lu_nn      Number of orbitals in LU Data distribution
      use scheComm
#ifdef MPI
      use mpi
#endif
      implicit none
!     Input Variables
      integer,        intent(in) :: NODE, NODES, N
!     Local Variables
      integer                    :: pp, ierr
      integer,            target :: tmp(1,1)
      integer,           pointer :: displs(:), LU2DD(:,:)

      if (NODE.eq.0) then
        nullify(displs)
        call re_alloc( displs, 1, NODES, 'displs', 'LU_FAC' )
        displs(1) = 0
        do PP= 1, NODES-1
          displs(PP+1) = displs(PP) + LUnNeigh(PP)*5
        enddo
      else
        LUnNeigh => tmp(:,1)
        lu_comm  => tmp
        displs   => tmp(:,1)
      endif

      call MPI_Scatter( LUnNeigh, 1, MPI_INTEGER, lu_ncom, 1,
     &                  MPI_INTEGER, 0, MPI_Comm_world, ierr )

      nullify(LU2DD)
      call re_alloc( LU2DD, 1, 5, 1, lu_ncom, 'lu_comm', 'LU_FAC' )
      if (NODE.eq.0) LUnNeigh = LUnNeigh*5

      call MPI_Scatterv( lu_comm(1,1), LUnNeigh, displs, MPI_INTEGER,
     &                   LU2DD(1,1), lu_ncom*5, MPI_INTEGER, 0,
     &                   MPI_Comm_world, ierr )

      if (node.eq.0) then
        call de_alloc( lu_comm, 'lu_comm', 'LU_FAC' )
        call de_alloc( displs, 'displs', 'LU_FAC' )
      endif

      lu_comm => LU2DD
      end subroutine distCommDD2LU

      subroutine factOverlap( node, nodes, nn, dd_nn, dd_nnz, numh,
     &                        listhptr, listh, S, gperm )
!     Compute LU factorization of Overlap matrix.
!     First we need to redistribute the matrix and change arrays
!     to C-style mode.
#ifdef MPI
      use mpi
#endif
      implicit none
!     Input Variables
      integer,      intent(in) :: node, nodes, nn, dd_nn, dd_nnz,
     &                            numh(dd_nn), listhptr(dd_nn),
     &                            listh(dd_nnz), gperm(*)
      real(dp),     intent(in) :: S(dd_nnz)
!     Local Variables
      integer                  :: PP, P1, P2, ierr, I1, I2, LE, NNZ,
     &                            DI, DF, OI, OF, stype, dtype, mtype,
     &                            fsr
      real(dp)                 :: berr
      real(dp),        pointer :: bb(:)
#ifdef MPI
      integer                  :: status(MPI_STATUS_SIZE)
#endif
      lu_nn = nn/nodes
      pp    = mod(nn,nodes)
      fsr   = lu_nn*node + min(node,pp)
      if (Node.lt.pp) lu_nn = lu_nn + 1

      nullify(lu_ia)
      call re_alloc( lu_ia, 1, lu_nn+1, 'lu_ia', 'LU_FAC' )
      do pp= 1, lu_ncom
        P1 = lu_comm(1,PP)
        P2 = lu_comm(2,PP)
        LE = lu_comm(3,PP)
        I1 = lu_comm(4,PP)
        I2 = lu_comm(5,PP)
        if (Node.eq.P1) then
          if (Node.eq.P2) then
            lu_ia(I2:I2+LE-1) = numh(I1:I1+LE-1)
          else
#ifdef MPI
            call MPI_Send( numh(I1), LE, MPI_INTEGER, P2, 0,
     &                     MPI_COMM_WORLD, ierr )
#endif
          endif
        else
#ifdef MPI
          call MPI_Recv( lu_ia(I2), LE, MPI_INTEGER, P1, 0,
     &                   MPI_COMM_WORLD, status, ierr )
#endif
        endif
      enddo

      P1       = lu_ia(1)
      lu_ia(1) = 1
      do pp= 1, lu_nn
        p2          = lu_ia(pp+1)
        lu_ia(pp+1) = lu_ia(pp) + p1
        p1          = p2
      enddo
      lu_nnz = lu_ia(lu_nn+1)-1

      nullify(lu_ja)
      call re_alloc( lu_ja, 1, lu_nnz, 'lu_ja', 'LU_FAC' )
      call re_alloc( lu_an, 1, lu_nnz, 'lu_an', 'LU_FAC' )

      do pp= 1, lu_ncom
        P1 = lu_comm(1,PP)
        P2 = lu_comm(2,PP)
        LE = lu_comm(3,PP)
        I1 = lu_comm(4,PP)
        I2 = lu_comm(5,PP)
        if (node.eq.P2) then
          DI = lu_ia(I2)
          DF = lu_ia(I2+LE)-1
          if (node.eq.P1) then
            OI = listhptr(I1)+1
            OF = listhptr(I1+LE-1)+numh(I1+LE-1)
            lu_ja(DI:DF) = listh(OI:OF)
            lu_an(DI:DF) = S(OI:OF)
          else
            nnz = DF - DI + 1
#ifdef MPI
            call MPI_Recv( lu_ja(DI), nnz, MPI_INTEGER,
     &                     P1, 0, MPI_COMM_WORLD, status, ierr )
            call MPI_Recv( lu_an(DI), nnz, MPI_DOUBLE_PRECISION,
     &                     P1, 0, MPI_COMM_WORLD, status, ierr )
#endif
          endif
        else
          OI = listhptr(I1)+1
          OF = listhptr(I1+LE-1)+numh(I1+LE-1)
          nnz = OF - OI + 1
#ifdef MPI
          call MPI_Send( listh(OI), nnz, MPI_INTEGER, P2, 0,
     &                   MPI_COMM_WORLD, ierr )
          call MPI_Send( S(OI), nnz, MPI_DOUBLE_PRECISION, P2, 0,
     &                   MPI_COMM_WORLD, ierr )
#endif
        endif
      enddo
      do i1= 1, lu_nnz
        lu_ja(i1) = gperm(lu_ja(i1))
      enddo

!     Change array to C-style mode
      lu_ia = lu_ia - 1
      lu_ja = lu_ja - 1

      call create_superlu_t( LU, NODES )

      stype = SLU_NR_loc ! distributed compressed row format
      dtype = SLU_D      ! Double precision
      mtype = SLU_GE     ! General matrix
      call f_dCreate_CompRowLoc_Mat_dist( LU%A, nn, nn, lu_nnz, lu_nn,
     &                                    fsr, lu_an, lu_ja, lu_ia,
     &                                    stype, dtype, mtype )
      call f_ScalePermstructInit( nn, nn, LU%ScalePermstruct )
      call f_LUstructInit( nn, nn, LU%LUstruct )
      call f_PStatInit( LU%stat )

!     Call the linear equation solver
      nullify(bb)
      call re_alloc( bb, 1, lu_nn, 'bb', 'LU_FAC' )
      bb = 1.0

      call f_pdgssvx( LU%options, LU%A, LU%ScalePermstruct, bb, lu_nn,
     &                1, LU%grid, LU%LUstruct, LU%SOLVEstruct, berr,
     &                LU%stat, ierr )

      if (ierr.eq.0) then
        write(23,*) 'Backward error: ', berr
        flush(23)
      else
        write(23,*) 'INFO from f_pdgssvx = ', ierr
        flush(23)
      endif
      call de_alloc( bb, 'bb', 'LU_FAC' )

      call set_superlu_options( LU%options, Fact=FACTORED,
     &                          IterRefine=NOREFINE, PrintStat=NO )

      end subroutine factOverlap

      subroutine freeFactData( )
      implicit none

!     Free Memory from SuperLU structure 
      call destroy_superlu_t( LU )

!     Free Memory from factOverlap 
      call de_alloc( lu_an, 'lu_an', 'LU_FAC' )
      call de_alloc( lu_ja, 'lu_ja', 'LU_FAC' )
      call de_alloc( lu_ia, 'lu_ia', 'LU_FAC' )

!     Free Memory from compCommDD2LU 
      call de_alloc( lu_comm, 'lu_comm', 'LU_FAC' )

!     Free Memory from compCommDD2LU 
      call de_alloc( LUnNeigh, 'LUnNeigh', 'LU_FAC' )
      end subroutine freeFactData


      subroutine reduce_array_lu( dir, n1, x, n2, y )
      use parallel,     only : node
      implicit none
C     Input variables
      integer,          intent(in) :: dir, n1, n2
      real*8,   target, intent(in) :: x(n1)
C     Output variables
      real*8,          intent(out) :: y(n2)
C     Local variables
      integer                     :: icom, P1, P2, N, II, I1, I2,
     &                               ierror, status(MPI_STATUS_SIZE)

      do icom= 1, lu_ncom
        if (dir.eq.1) then
          P1 = lu_comm(1,icom)
          P2 = lu_comm(2,icom)
          N  = lu_comm(3,icom)
          I1 = lu_comm(4,icom)
          I2 = lu_comm(5,icom)
        else
          P2 = lu_comm(1,icom)
          P1 = lu_comm(2,icom)
          N  = lu_comm(3,icom)
          I2 = lu_comm(4,icom)
          I1 = lu_comm(5,icom)
        endif
        if (NODE.eq.P1) then
          if (NODE.eq.P2) then
            do II= 1, N
              y(I2) = x(I1)
              I1 = I1 + 1
              I2 = I2 + 1
            enddo
          else
            call MPI_Send ( x(I1), N, MPI_DOUBLE_PRECISION, P2, 0,
     &                      MPI_COMM_WORLD, ierror )
          endif
        else
          call MPI_Recv ( y(I2), N, MPI_DOUBLE_PRECISION, P1, 0,
     &                    MPI_COMM_WORLD, status, ierror )
        endif
      enddo
      end subroutine reduce_array_lu


      END MODULE LU_FACTORIZATION
