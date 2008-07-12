#ifdef GRID_DP
#define MPI_GRID_PRECISION MPI_DOUBLE_PRECISION
#else
#define MPI_GRID_PRECISION MPI_REAL
#endif
      subroutine reordMesh( fsrc, fdst, nm, nsm, itr )
      use precision, only: grid_p
      use MeshComm
      use parallel
#ifdef MPI
      use mpi_siesta
#endif
      implicit none
C     Input variables
      integer, intent(in)         :: nm(3), nsm, itr 
      real(grid_p), intent(inout) :: fsrc(*), fdst(*)
C     Local variables
      integer                     :: icom, Xsize, Ysize, Zsize,
     &                               MaxSize, XYsrc, XYdst, NSP,
     &                               I1, I2, I3, J1, J2, J3, K1, K2, K3,
     &                               N1, N2, N3, NN, JJ, KK
      integer                     :: NDST(3), NSRC(3)
      real(grid_p), allocatable   :: TBUF(:)
      integer, allocatable        :: JS(:)
#ifdef MPI
      integer                     :: MPIerror, Status(MPI_Status_Size)
#endif

!------------------------------------------------------------------------- BEGIN
      NSP = NSM*NSM*NSM
      allocate(JS(NSP))
      call memory('A','I',NSP,'reordMesh')

C     Compute the maximum size of the buffer needed to transfer data
C     among the several processes
      maxSize = 0
      Xsize   = nm(1)
      do icom= 1, MCsize
        if (MCsrc(icom).ne.MCdst(icom)) then
          Ysize = MCfiY(icom) - MCinY(icom) + 1
          Zsize = MCfiZ(icom) - MCinZ(icom) + 1
          MaxSize = max(MaxSize,Ysize*Zsize)
        endif
      enddo

      MaxSize = Xsize*MaxSize*nsp
      if (MaxSize.gt.0) then
        allocate(TBUF(MaxSize))
        call memory('A','X',MaxSize,'reordMesh')
      endif


      if (itr.eq.1) then
        NDST(1) = (udMeshLim(2,1) - udMeshLim(1,1) + 1)*nsm
        NDST(2) = (udMeshLim(2,2) - udMeshLim(1,2) + 1)*nsm
        NDST(3) = (udMeshLim(2,3) - udMeshLim(1,3) + 1)*nsm

        NSRC(1) = (nddMeshLim(2,1) - nddMeshLim(1,1) + 1)*nsm
        NSRC(2) = (nddMeshLim(2,2) - nddMeshLim(1,2) + 1)*nsm
        NSRC(3) = (nddMeshLim(2,3) - nddMeshLim(1,3) + 1)*nsm

        NN = 1
        I3 = 0
        DO N3=0, NSM-1
          I2 = 0
          DO N2= 0, NSM-1
            I1 = I2 + I3
            DO N1= 0, NSM-1
              JS(NN) = I1
              NN     = NN + 1
              I1     = I1 + 1
            ENDDO
            I2 = I2 + NDST(1)
          ENDDO
          I3 = I3 + NDST(1)*NDST(2)
        ENDDO
        do icom= 1, MCsize
          if (MCsrc(icom).eq.Node) then
            if (MCdst(icom).eq.Node) then
              J3 = (MCinZ(icom) - nddMeshLim(1,3))*NSRC(1)*NSRC(2)*NSM
              K3 = (MCinZ(icom) - udMeshLim(1,3))*NDST(1)*NDST(2)*NSM
              do I3 = MCinZ(icom), MCfiZ(icom)
                J2 = (MCinY(icom) - nddMeshLim(1,2))*NSRC(1)*NSM*NSM
                K2 = (MCinY(icom) - udMeshLim(1,2))*NDST(1)*NSM
                do I2 = MCinY(icom), MCfiY(icom)
                  J1 = 1 + J2 + J3
                  K1 = 1 + K2 + K3
                  do I1 = 1, Xsize
                    DO NN= 1, NSP
                      fdst(K1+JS(NN)) = fsrc(J1)
                      J1 = J1 + 1
                    ENDDO
                    K1 = K1 + NSM
                  enddo
                  J2 = J2 + NSRC(1)*NSM*NSM
                  K2 = K2 + NDST(1)*NSM
                enddo
                J3 = J3 + NSRC(1)*NSRC(2)*NSM
                K3 = K3 + NDST(1)*NDST(2)*NSM
              enddo
            else
              Ysize = MCfiY(icom) - MCinY(icom) + 1
              Zsize = MCfiZ(icom) - MCinZ(icom) + 1
              call mpi_recv( TBUF, Xsize*Ysize*Zsize*nsp,
     &                       MPI_GRID_PRECISION, MCdst(icom), 1,
     &                       MPI_Comm_world, Status, MPIerror )

              J1 = 1
              K3 = (MCinZ(icom) - udMeshLim(1,3))*NDST(1)*NDST(2)*NSM
              do I3 = MCinZ(icom), MCfiZ(icom)
                K2 = (MCinY(icom) - udMeshLim(1,2))*NDST(1)*NSM
                do I2 = MCinY(icom), MCfiY(icom)
                  K1 = 1 + K2 + K3
                  do I1 = 1, Xsize
                    DO NN= 1, NSP
                      fdst(K1+JS(NN)) = TBUF(J1)
                      J1 = J1 + 1
                    ENDDO
                    K1 = K1 + NSM
                  enddo
                  K2 = K2 + NDST(1)*NSM
                enddo
                K3 = K3 + NDST(1)*NDST(2)*NSM
              enddo
            endif
          else
            J3 = (MCinZ(icom) - nddMeshLim(1,3))*NSRC(1)*NSRC(2)*NSM
            K1 = 1
            do I3 = MCinZ(icom), MCfiZ(icom)
              J2 = (MCinY(icom) - nddMeshLim(1,2))*NSRC(1)*NSM*NSM
              do I2 = MCinY(icom), MCfiY(icom)
                J1 = 1 + J2 + J3
                do I1 = 1, Xsize
                  DO NN= 1, NSP
                    TBUF(K1) = fsrc(J1)
                    J1 = J1 + 1
                    K1 = K1 + 1
                  ENDDO
                enddo
                J2 = J2 + NSRC(1)*NSM*NSM
              enddo
              J3 = J3 + NSRC(1)*NSRC(2)*NSM
            enddo

            Ysize = MCfiY(icom) - MCinY(icom) + 1
            Zsize = MCfiZ(icom) - MCinZ(icom) + 1
            call MPI_Send( TBUF, Xsize*Ysize*Zsize*nsp,
     &                     MPI_GRID_PRECISION, MCsrc(icom), 1,
     &                     MPI_Comm_world, MPIerror )
          endif
        enddo
      else
        NSRC(1) = (udMeshLim(2,1) - udMeshLim(1,1) + 1)*nsm
        NSRC(2) = (udMeshLim(2,2) - udMeshLim(1,2) + 1)*nsm
        NSRC(3) = (udMeshLim(2,3) - udMeshLim(1,3) + 1)*nsm

        NDST(1) = (nddMeshLim(2,1) - nddMeshLim(1,1) + 1)*nsm
        NDST(2) = (nddMeshLim(2,2) - nddMeshLim(1,2) + 1)*nsm
        NDST(3) = (nddMeshLim(2,3) - nddMeshLim(1,3) + 1)*nsm

        NN = 1
        I3 = 0
        DO N3=0, NSM-1
          I2 = 0
          DO N2= 0, NSM-1
            I1 = I2 + I3
            DO N1= 0, NSM-1
              JS(NN) = I1
              NN     = NN + 1
              I1     = I1 + 1
            ENDDO
            I2 = I2 + NSRC(1)
          ENDDO
          I3 = I3 + NSRC(1)*NSRC(2)
        ENDDO

        do icom= 1, MCsize
          if (MCsrc(icom).eq.Node) then
            if (MCdst(icom).eq.Node) then
              J3 = (MCinZ(icom) - udMeshLim(1,3))*NSRC(1)*NSRC(2)*NSM
              K3 = (MCinZ(icom) - nddMeshLim(1,3))*NDST(1)*NDST(2)*NSM
              do I3 = MCinZ(icom), MCfiZ(icom)
                J2 = (MCinY(icom) - udMeshLim(1,2))*NSRC(1)*NSM
                K2 = (MCinY(icom) - nddMeshLim(1,2))*NDST(1)*NSM*NSM
                do I2 = MCinY(icom), MCfiY(icom)
                  J1 = 1 + J2 + J3
                  K1 = 1 + K2 + K3
                  do I1 = 1, Xsize
                    DO NN= 1, NSP
                      fdst(K1) = fsrc(J1+JS(NN))
                      K1 = K1 + 1
                    ENDDO
                    J1 = J1 + NSM
                  enddo
                  J2 = J2 + NSRC(1)*NSM
                  K2 = K2 + NDST(1)*NSM*NSM
                enddo
                J3 = J3 + NSRC(1)*NSRC(2)*NSM
                K3 = K3 + NDST(1)*NDST(2)*NSM
              enddo
            else
              J3 = (MCinZ(icom) - udMeshLim(1,3))*NSRC(1)*NSRC(2)*NSM
              K1 = 1
              do I3 = MCinZ(icom), MCfiZ(icom)
                J2 = (MCinY(icom) - udMeshLim(1,2))*NSRC(1)*NSM
                do I2 = MCinY(icom), MCfiY(icom)
                  J1 = 1 + J2 + J3
                  do I1 = 1, Xsize
                    DO NN= 1, NSP
                      TBUF(K1) = fsrc(J1+JS(NN))
                      K1       = K1 + 1
                    ENDDO
                    J1 = J1 + NSM
                  enddo
                  J2 = J2 + NSRC(1)*NSM
                enddo
                J3 = J3 + NSRC(1)*NSRC(2)*NSM
              enddo

              Ysize = MCfiY(icom) - MCinY(icom) + 1
              Zsize = MCfiZ(icom) - MCinZ(icom) + 1
              call MPI_Send( TBUF, Xsize*Ysize*Zsize*nsp,
     &                       MPI_GRID_PRECISION, MCdst(icom), 1,
     &                       MPI_Comm_world, MPIerror )
            endif
          else
            Ysize = MCfiY(icom) - MCinY(icom) + 1
            Zsize = MCfiZ(icom) - MCinZ(icom) + 1
            call mpi_recv( TBUF, Xsize*Ysize*Zsize*nsp,
     &                     MPI_GRID_PRECISION, MCsrc(icom), 1,
     &                     MPI_Comm_world, Status, MPIerror )
            J1 = 1
            K3 = (MCinZ(icom) - nddMeshLim(1,3))*NDST(1)*NDST(2)*NSM
            do I3 = MCinZ(icom), MCfiZ(icom)
              K2 = (MCinY(icom) - nddMeshLim(1,2))*NDST(1)*NSM*NSM
              do I2 = MCinY(icom), MCfiY(icom)
                K1 = 1 + K2 + K3
                do I1 = 1, Xsize
                  DO NN= 1, NSP
                    fdst(K1) = TBUF(J1)
                    K1 = K1 + 1
                    J1 = J1 + 1
                  ENDDO
                enddo
                K2 = K2 + NDST(1)*NSM*NSM
              enddo
              K3 = K3 + NDST(1)*NDST(2)*NSM
            enddo
          endif
        enddo
      endif

      if (MaxSize.gt.0) then
        deallocate(TBUF)
        call memory('D','X',MaxSize,'reordMesh')
      endif
      deallocate(JS)
      call memory('D','I',NSP,'reordMesh')

C     Syncronize processes in order to avoid problems
      call MPI_Barrier(MPI_Comm_World,MPIerror)

!--------------------------------------------------------------------------- END
      end subroutine reordMesh
