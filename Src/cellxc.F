      SUBROUTINE CELLXC( FUNCTL, AUTHOR, IREL, IDER, 
     .                   CELL, NMESH, NSPAN, MAXP, MTYPE, XMESH,
     .                   NSPIN, DENS,
     .                   EX, EC, DX, DC, VXC, IAUXDER, DVXCDN,
     .                   STRESSL, MAXAUX, AUX,
     .                   NMESHG, NSM, Node, Nodes )

C *******************************************************************
C Finds total exchange-correlation energy and potential in a
C   periodic cell.
C This version implements the Local (spin) Density Approximation and
C   the Generalized-Gradient-Aproximation with the 'explicit mesh 
C   functional' method of White & Bird, PRB 50, 4954 (1994).
C Gradients are 'defined' by numerical derivatives, using 2*NN+1 mesh
C   points, where NN is a parameter defined below
C Coded by L.C.Balbas and J.M.Soler. December 1996.
C Stress added by J.L.Martins and J.M.Soler. August 1997.
C ************************* INPUT ***********************************
C CHARACTER*(*) FUNCTL : Functional to be used:
C              'LDA' or 'LSD' => Local (spin) Density Approximation
C                       'GGA' => Generalized Gradient Corrections
C                                Uppercase is optional
C CHARACTER*(*) AUTHOR : Parametrization desired:
C     'CA' or 'PZ' => LSD Perdew & Zunger, PRB 23, 5075 (1981)
C           'PW92' => LSD Perdew & Wang, PRB, 45, 13244 (1992). This is
C                     the local density limit of the next:
C            'PBE' => GGA Perdew, Burke & Ernzerhof, PRL 77, 3865 (1996)
C                     Uppercase is optional
C INTEGER IREL         : Relativistic exchange? (0=>no, 1=>yes)
C INTEGER IDER         : Store derivatives of exchange-correlation potential?
C                        (0=>no, 1=>yes)
C REAL*8  CELL(3,3)    : Unit cell vectors CELL(ixyz,ivector)
C INTEGER NMESH(3)     : Number of mesh divisions of each vector
C INTEGER NSPAN(3)     : Physical dimensions of arrays XMESH, DENS and
C                        VXC (or memory span between array elements)
C                        See usage section for more information
C INTEGER MAXP         : Physical dimension of XMESH, DENS, and VXC
C INTEGER MTYPE        : Mesh type:
C                        0 => Uniform mesh
C                        1 => Adaptive mesh, given in cartesian coord
C                        2 => Adaptive mesh, given in cell-vector coord
C REAL    XMESH(3,MAXP): Mesh point coordinates (not used if MTYPE=0)
C                        When MTYPE=2, cartesian coordinates are
C                        Xcart(ix,im) = Sum_iv(CELL(ix,iv)*XMESH(iv,ip))
C                        Notice single precision in this version
C INTEGER NSPIN        : NSPIN=1 => unpolarized; NSPIN=2 => polarized;
C                        NSPIN=4 => non-collinear polarization
C REAL    DENS(MAXP,NSPIN) : Total (NSPIN=1) or spin (NSPIN=2) electron
C                        density at mesh points, ordered as
C                        IP = I1+NSPAN(1)*((I2-1)+NSPAN(2)*(I3-1)),
C                        with I1=1,...,NMESH(1), etc 
C                        For non-collinear polarization, the density
C                        matrix is given by: DENS(1)=D11, DENS(2)=D22,
C                        DENS(3)=Real(D12), DENS(4)=Im(D12)
C                        Notice single precision in this version
C INTEGER MAXAUX       : Physical dimension of array AUX.
C                        If FUNCTL.EQ.'GGA' and MTYPE.NE.0, must be
C                        MAXAUX .GE. NSPAN(1)*NSPAN(2)*NSPAN(3)
C                        Not used otherwise
C INTEGER NMESHG(3)    : Global number of mesh divisions of each vector
C INTEGER NSM          : Number of sub-mesh points along each axis
C INTEGER IAUXDER      : First dimension of DVXCDN,
C                        1 if Khon-Sham forces or ntpl if Harris forces.
C ************************* OUTPUT **********************************
C REAL*8  EX              : Total exchange energy
C REAL*8  EC              : Total correlation energy
C REAL*8  DX              : IntegralOf( rho * (eps_x - v_x) )
C REAL*8  DC              : IntegralOf( rho * (eps_c - v_c) )
C REAL    VXC(MAXP,NSPIN) : (Spin) exch-corr potential
C                           Notice single precision in this version
C REAL    DVXCDN(MAXP,NSPIN,NSPIN) : Derivatives of exchange-correlation
C                           potential respect the charge density
C REAL*8  STRESSL(3,3)    : xc contribution to the stress tensor,
C                           assuming constant density (not charge),
C                           i.e. r->r' => rho'(r') = rho(r)
C                           For plane-wave and grid (finite diff) basis
C                           sets, density rescaling gives an extra term
C                           (not included) (DX+DC-EX-EC)/cell_volume for
C                           the diagonal elements of stress. For other
C                           basis sets, the extra term is, in general:
C                           IntegralOf(v_xc * d_rho/d_strain) / cell_vol
C ************************* AUXILIARY *******************************
C REAL    AUX(MAXAUX)     : Not used unless 
C                           FUNCTL.EQ.'GGA' .AND. MTYPE.NE.0
C ************************ UNITS ************************************
C Distances in atomic units (Bohr).
C Densities in atomic units (electrons/Bohr**3)
C Energy unit depending of parameter EUNIT below
C Stress in EUNIT/Bohr**3
C ************************ USAGE ************************************
C Typical calls for different array dimensions:
C     PARAMETER ( MAXP = 1000000 )
C     DIMENSION NMESH(3), DENS(MAXP,2), VXC(MAXP,2)
C     Find CELL vectors
C     Find density at N1*N2*N3 mesh points (less than MAXP) and place 
C       them consecutively in array DENS
C     NMESH(1) = N1
C     NMESH(2) = N2
C     NMESH(3) = N3
C     CALL CELLXC( 'GGA', 'PBE', 0,
C    .              CELL, NMESH, NMESH, MAXP, 0, 0,
C    .              2, DENS,
C    .              EX, EC, DX, DC, VXC, STRESS, 0, 0 )
C Or alternatively:
C     PARAMETER ( M1=100, M2=100, M3=100 )
C     DIMENSION NMESH(3), NSPAN(3), DENS(M1,M2,M3,2), VXC(M1,M2,M3,2)
C     DATA NSPAN / M1, M2, M3 /
C     Find CELL vectors
C     Find DENS at N1*N2*N3 mesh points
C     NMESH(1) = N1
C     NMESH(2) = N2
C     NMESH(3) = N3
C     CALL CELLXC( 'GGA', 'PBE', 0,
C    .              CELL, NMESH, NSPAN, M1*M2*M3, 0, 0,
C    .              2, DENS,
C    .              EX, EC, DX, DC, VXC, STRESS, 0, 0 )
C ********* BEHAVIOUR ***********************************************
C - Notice that XMESH and AUX are not used if MTYPE=0. This means that
C   you do not even need to declare them in the calling program. See
C   usage section for calling examples.
C - If FNCTNL='LDA', DENS and VXC may be the same physical array.
C - Stops and prints a warning if arguments FUNCTL or AUTHOR are not
C   one of the allowed possibilities.
C - Since the exchange and correlation part is usually a small fraction
C   of a typical electronic structure calculation, this routine has
C   been coded with emphasis on simplicity and functionality, not in
C   eficiency. It is not designed to vectorize or parallelize well.
C ********* ROUTINES CALLED *****************************************
C GGAXC, LDAXC, RECLAT, TIMER, VOLCEL
C *******************************************************************

C
C  Modules
C
      use precision
#ifdef MPI
      use mpi
      use parallel
#endif

      IMPLICIT NONE

C Argument types and dimensions
      CHARACTER*(*)     FUNCTL, AUTHOR
      INTEGER           IREL, IDER, MAXAUX, MAXP, MTYPE, NMESHG(3), NSM,
     .                  NMESH(3), NSPAN(3), NSPIN, Node, Nodes, IAUXDER
      DOUBLE PRECISION  CELL(3,3), DC, DX, EC, EX, STRESSL(3,3)
C     If you change next line, please change also the argument
C     explanations in the header comments
*     DOUBLE PRECISION
      REAL
     .                  AUX(*), DENS(MAXP,NSPIN),
     .                  VXC(MAXP,NSPIN), XMESH(3,*)
      REAL, DIMENSION(IAUXDER,NSPIN,NSPIN) :: DVXCDN
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::
     .                  DVXDN, DVCDN
      INTEGER 
     .                  IS2

C Fix the order of the numerical derivatives
C NN is the number of points used in each coordinate and direction,
C i.e. a total of 6*NN neighbour points is used to find the gradients
      INTEGER NN
      PARAMETER ( NN = 3 )

C Fix energy unit:  EUNIT=1.0 => Hartrees,
C                   EUNIT=0.5 => Rydbergs,
C                   EUNIT=0.03674903 => eV
      DOUBLE PRECISION EUNIT
      PARAMETER ( EUNIT = 0.5D0 )

C Fix switch to skip points with zero density
      LOGICAL SKIP0
      DOUBLE PRECISION DENMIN
      PARAMETER ( SKIP0  = .TRUE. )
      PARAMETER ( DENMIN = 1.D-15 )

C MPI related variables
#ifdef MPI
      integer 
     .  MPIerror, ProcessorZ, BlockSizeY, BlockSizeZ, MeshNsm(3),
     .  BlockSizeYMax, BlockSizeZMax, NRemY, NRemZ, Py, Pz
      integer 
     .  NSwap, NDummy, IDen1, IDen2, IOut, INN, MeshDum(3),
     .  Requests(4), Request, Status(MPI_Status_Size), NReq,
     .  Statuss(MPI_Status_Size,4), Tag, Pleft2, Pright2,
     .  Pleft3, Pright3, PzMod, PyMod, NMeshPN, JPNN(3,-NN:NN)
      double precision 
     .  Buffer(4),Buffer2(4)
      logical 
     .  NoLocalPoints

      real, dimension(:,:,:), allocatable, save :: 
     .  bdens, bvxc
      integer, dimension(:,:), allocatable, save :: 
     .  PMesh, Statm
      integer, dimension(:), allocatable, save :: 
     .  Reqm
#endif

C Local variables and arrays
      LOGICAL           GGA
      INTEGER           I1, I2, I3, IC, IN, IP, IS, IX,
     .                  J1, J2, J3, JN, JP(3,-NN:NN), JX, NPG
      DOUBLE PRECISION  DENTOT, DGDM(-NN:NN), DGIDFJ(3,3,-NN:NN),
     .                  DMDX(3,3), DVOL, DXDM(3,3),
     .                  EPSC, EPSX, F1, F2, 
     .                  VOLCEL, VOLUME, STRESS(3,3)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE ::
     .                  D, DECDD, DEXDD
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE ::
     .                  DECDGD, DEXDGD, GD

      EXTERNAL          GGAXC, LDAXC, RECLAT, VOLCEL, memory

C Start time counter (intended only for debugging and development)
      CALL TIMER( 'CELLXC', 1 )


C Allocate local memory
      allocate(D(nspin))
      call memory('A','D',nspin,'cellxc')
      allocate(DECDD(nspin))
      call memory('A','D',nspin,'cellxc')
      allocate(DEXDD(nspin))
      call memory('A','D',nspin,'cellxc')
      allocate(DECDGD(3,nspin))
      call memory('A','D',3*nspin,'cellxc')
      allocate(DEXDGD(3,nspin))
      call memory('A','D',3*nspin,'cellxc')
      allocate(GD(3,nspin))
      call memory('A','D',3*nspin,'cellxc')
      allocate(DVXDN(NSPIN,NSPIN))
      call memory('A','D',nspin*nspin,'cellxc')
      allocate(DVCDN(NSPIN,NSPIN))
      call memory('A','D',nspin*nspin,'cellxc')

C Set GGA switch
      IF ( FUNCTL.EQ.'LDA' .OR. FUNCTL.EQ.'lda' .OR.
     .     FUNCTL.EQ.'LSD' .OR. FUNCTL.EQ.'lsd' ) THEN
        GGA = .FALSE.
      ELSEIF ( FUNCTL.EQ.'GGA' .OR. FUNCTL.EQ.'gga') THEN
        GGA = .TRUE.
      ELSE
        WRITE(6,*) 'CELLXC: Unknown functional ', FUNCTL
        STOP
      ENDIF

#ifdef MPI
C If GGA and the number of processors is greater than 1 then we need
C to exchange border densities so that the density derivatives can
C be calculated.
      if (GGA.and.Nodes.gt.1) then

C Allocate request and status storage
        allocate(Reqm(4*NMESH(3)))
        call memory('A','I',4*NMESH(3),'cellxc')
        allocate(Statm(MPI_Status_Size,4*NMESH(3)))
        call memory('A','I',4*NMESH(3)*MPI_Status_Size,'cellxc')

C Work out processor grid dimensions
        if (mod(Nodes,ProcessorY).gt.0) then
          write(6,'(''ERROR: ProcessorY must be a factor of the'',
     .      '' number of processors!'')')
          stop
        endif
        ProcessorZ = Nodes/ProcessorY

C Work out grid coordinates of local node in processor grid
        Py = (Node/ProcessorZ) + 1
        Pz = Node - (Py - 1)*ProcessorZ + 1

C Work out maximum blocksizes
        MeshNSM(1) = NMESHG(1)/NSM
        MeshNSM(2) = NMESHG(2)/NSM
        MeshNSM(3) = NMESHG(3)/NSM
        BlockSizeY = (MeshNsm(2)/ProcessorY)
        BlockSizeZ = (MeshNsm(3)/ProcessorZ)
        NRemY = MeshNsm(2) - BlockSizeY*ProcessorY
        NRemZ = MeshNsm(3) - BlockSizeZ*ProcessorZ
        if (NRemY.gt.0) then
          BlockSizeYMax = BlockSizeY + 1
        else
          BlockSizeYMax = BlockSizeY
        endif
        if (NRemZ.gt.0) then
          BlockSizeZMax = BlockSizeZ + 1
        else
          BlockSizeZMax = BlockSizeZ
        endif
        BlockSizeY = BlockSizeY*NSM
        BlockSizeZ = BlockSizeZ*NSM
        BlockSizeYMax = BlockSizeYMax*NSM
        BlockSizeZMax = BlockSizeZMax*NSM

C Check that there are enough points on each Node to cover the finite
C difference interval - if not then exit!
        if (NN.gt.BlockSizeY.or.NN.gt.BlockSizeZ) then
          if (Node.eq.0) then
            write(6,'(''  ERROR - number of fine mesh points per '',
     .        ''Node must be greater than finite difference order '')')
          endif
          stop
        endif

C Check that this Node was some points 
        call HowManyMeshPerNode(MeshNSM, Node, Nodes, NDummy, MeshDum)
        NoLocalPoints = (NDummy.eq.0)
        if (NoLocalPoints) then
          EX = 0.0d0
          EC = 0.0d0
          DX = 0.0d0
          DC = 0.0d0
          if (GGA) then
            do IS = 1,NSPIN
              do I3 = 0,NMESH(3)-1
                do I2 = 0,NMESH(2)-1
                  do I1 = 0,NMESH(1)-1
                    IP = 1+I1+NSPAN(1)*I2+NSPAN(1)*NSPAN(2)*I3
                    VXC(IP,IS) = 0.0
                  enddo
                enddo
              enddo
            enddo
          endif
          goto 991
        endif

C Allocate memory for border densities and potential
        allocate(bdens(NSPAN(1)*max(BlockSizeYMax,BlockSizeZMax),
     .    4*NN,NSPIN))
        call memory('A','S',NSPAN(1)*max(BlockSizeYMax,BlockSizeZMax)*
     .    4*NN*NSPIN,'cellxc')
        allocate(bvxc(NSPAN(1)*max(BlockSizeYMax,BlockSizeZMax),
     .    4*NN,NSPIN))
        call memory('A','S',NSPAN(1)*max(BlockSizeYMax,BlockSizeZMax)*
     .    4*NN*NSPIN,'cellxc')
        allocate(PMesh(3,4))
        call memory('A','I',12,'cellxc')

C Exchange density information in the Z-direction
        if (NMESH(3).ne.NMESHG(3)) then

C Find Nodes to exchange with
          PzMod = Pz - 1
          if (PzMod.lt.1) PzMod = PzMod + ProcessorZ
          Pleft3 = (Py - 1)*ProcessorZ + PzMod - 1
          PzMod = Pz + 1
          if (PzMod.gt.ProcessorZ) PzMod = PzMod - ProcessorZ
          Pright3 = (Py - 1)*ProcessorZ + PzMod - 1

C Find size of mesh at nodes to be exchanged with
          call HowManyMeshPerNode(MeshNsm,Pleft3,Nodes,NMeshPN,
     .      PMesh(1,3))
          call HowManyMeshPerNode(MeshNsm,Pright3,Nodes,NMeshPN,
     .      PMesh(1,4))
          do ix = 1,3
            PMesh(ix,3) = PMesh(ix,3)*NSM
            PMesh(ix,4) = PMesh(ix,4)*NSM
          enddo

C Loop over points of buffer region
          do in = 1,nn

C Loop over spins
            do is = 1,nspin

C Post density receive
              NSwap = PMesh(1,3)*PMesh(2,3)
              call MPI_IRecv(bdens(1,2*nn+in,is),NSwap,DAT_single,
     $             Pleft3,1,MPI_Comm_World,Requests(1),MPIerror)
              NSwap = PMesh(1,4)*PMesh(2,4)
              call MPI_IRecv(bdens(1,3*nn+in,is),NSwap,DAT_single,
     $             Pright3,2,MPI_Comm_World,Requests(2),MPIerror)

C Post density send
              NSwap = NSPAN(1)*NSPAN(2)
              IDen1 = NSwap*(NMESH(3)-in)+1
              call MPI_ISend(dens(IDen1,is),NSwap,DAT_single,Pright3,
     .          1,MPI_Comm_World,Requests(3),MPIerror)
              IDen1 = NSwap*(in-1)+1
              call MPI_ISend(dens(IDen1,is),NSwap,DAT_single,Pleft3,
     .          2,MPI_Comm_World,Requests(4),MPIerror)

C Wait for transfers to complete
              call MPI_WaitAll(4,Requests,Statuss,MPIerror)

            enddo

          enddo

        endif

C Exchange density information in the Y-direction
        if (NMESH(2).ne.NMESHG(2)) then

C Find Nodes to exchange with
          PyMod = Py - 1
          if (PyMod.lt.1) PyMod = PyMod + ProcessorY
          Pleft2 = (PyMod - 1)*ProcessorZ + Pz - 1
          PyMod = Py + 1
          if (PyMod.gt.ProcessorY) PyMod = PyMod - ProcessorY
          Pright2 = (PyMod - 1)*ProcessorZ + Pz - 1

C Find size of mesh at nodes to be exchanged with
          call HowManyMeshPerNode(MeshNsm,Pleft2,Nodes,NMeshPN,
     .      PMesh(1,1))
          call HowManyMeshPerNode(MeshNsm,Pright2,Nodes,NMeshPN,
     .      PMesh(1,2))
          do ix = 1,3
            PMesh(ix,1) = PMesh(ix,1)*NSM
            PMesh(ix,2) = PMesh(ix,2)*NSM
          enddo

C Loop over points of buffer region
          do in = 1,nn

C Loop over spins
            do is = 1,nspin

C Post density receives
              do ix = 1,PMesh(3,1)
                NSwap = PMesh(1,1)
                IDen1 = (ix-1)*PMesh(1,1)+1
                call MPI_IRecv(bdens(IDen1,in,is),NSwap,DAT_single,
     .            Pleft2,ix,MPI_Comm_World,Reqm(ix),MPIerror)
              enddo
              NReq = PMesh(3,1)
              do ix = 1,PMesh(3,2)
                NSwap = PMesh(1,2)
                IDen1 = (ix-1)*PMesh(1,2)+1
                call MPI_IRecv(bdens(IDen1,nn+in,is),NSwap,DAT_single,
     .            Pright2,PMesh(3,1)+ix,MPI_Comm_World,Reqm(NReq+ix),
     .            MPIerror)
              enddo
              NReq = NReq + PMesh(3,2)

C Post density sends
              NSwap = NSPAN(1)
              IDen1 = (NSPAN(2)-in)*NSPAN(1) + 1
              do ix = 1,NMESH(3)
                call MPI_ISend(dens(IDen1,is),NSwap,DAT_single,Pright2,
     .            ix,MPI_Comm_World,Reqm(NReq+ix),MPIerror)
                IDen1 = IDen1 + NSPAN(1)*NSPAN(2)
              enddo
              NReq = NReq + NMESH(3)
              IDen1 = (in-1)*NSPAN(1) + 1
              do ix = 1,NMESH(3)
                call MPI_ISend(dens(IDen1,is),NSwap,DAT_single,Pleft2,
     .            NMESH(3)+ix,MPI_Comm_World,Reqm(NReq+ix),MPIerror)
                IDen1 = IDen1 + NSPAN(1)*NSPAN(2)
              enddo
              NReq = NReq + NMESH(3)

C Wait for transfers to complete
              call MPI_WaitAll(NReq,Reqm,Statm,MPIerror)

            enddo

          enddo

        endif

      else

        NoLocalPoints = .false.

      endif

C Barrier to ensure that all densities are in place
  991 call MPI_Barrier(MPI_Comm_World,MPIerror)
      if (NoLocalPoints) goto 992
#endif

C Find weights of numerical derivation from Lagrange interp. formula
      DO 15 IN = -NN,NN
        F1 = 1.0d0
        F2 = 1.0d0
        DO 10 JN = -NN,NN
          IF (JN.NE.IN .AND. JN.NE.0) F1 = F1 * (0  - JN)
          IF (JN.NE.IN)               F2 = F2 * (IN - JN)
   10   CONTINUE
        DGDM(IN) = F1 / F2
   15 CONTINUE
      DGDM(0) = 0.0d0

C Find total number of mesh points
      NPG = NMESHG(1) * NMESHG(2) * NMESHG(3)

C Find Jacobian matrix dx/dmesh for uniform mesh
      IF (MTYPE.EQ.0) THEN

C       Find mesh cell volume 
        DVOL = VOLCEL( CELL ) / DBLE(NPG)

        IF (GGA) THEN

C         Find mesh unit vectors and reciprocal mesh vectors
          DO 25 IC = 1,3
            DO 20 IX = 1,3
              DXDM(IX,IC) = CELL(IX,IC) / NMESHG(IC)
   20       CONTINUE
   25     CONTINUE
          CALL RECLAT( DXDM, DMDX, 0 )

C         Find the weights for the derivative d(gradF(i))/d(F(j)) of
C         the gradient at point i with respect to the value at point j
          DO 50 IN = -NN,NN
            DO 40 IC = 1,3
              DO 30 IX = 1,3
                DGIDFJ(IX,IC,IN) = DMDX(IX,IC) * DGDM(IN)
   30         CONTINUE
   40       CONTINUE
   50     CONTINUE

        ENDIF

      ENDIF

C Initialize output
      EX = 0.0d0
      EC = 0.0d0
      DX = 0.0d0
      DC = 0.0d0
      IF (GGA) THEN
        DO 70 IS = 1,NSPIN
          DO 66 I3 = 0,NMESH(3)-1
            DO 64 I2 = 0,NMESH(2)-1
            DO 62 I1 = 0,NMESH(1)-1
              IP = 1 + I1 + NSPAN(1) * I2 + NSPAN(1) * NSPAN(2) * I3
              VXC(IP,IS) = 0.0
   62       CONTINUE
   64       CONTINUE
   66     CONTINUE
   70   CONTINUE
#ifdef MPI
C Initialise buffer regions of Vxc
        if (Nodes.gt.1) then
          do is = 1,nspin
            if (NMESH(2).ne.NMESHG(2)) then
              NSwap = NSPAN(1)*BlockSizeZMax
              do in = 1,2*nn
                do ix = 1,NSwap
                  bvxc(ix,in,is) = 0.0d0
                enddo
              enddo
            endif
            if (NMESH(3).ne.NMESHG(3)) then
              NSwap = NSPAN(1)*BlockSizeYMax
              do in = 1,2*nn
                do ix = 1,NSwap
                  bvxc(ix,2*nn+in,is) = 0.0d0
                enddo
              enddo
            endif
          enddo
        endif
#endif
      ENDIF
      DO 74 JX = 1,3
        DO 72 IX = 1,3
          STRESS(IX,JX) = 0.0D0
   72   CONTINUE
   74 CONTINUE

C Loop on mesh points
      DO 230 I3 = 0,NMESH(3)-1
      DO 220 I2 = 0,NMESH(2)-1
      DO 210 I1 = 0,NMESH(1)-1

C       Find mesh index of this point
        IP = 1 + I1 + NSPAN(1) * I2 + NSPAN(1) * NSPAN(2) * I3

C       Skip point if density=0
        IF (SKIP0) THEN
          DENTOT = 0.D0
          DO IS = 1,MIN(NSPIN,2)
            DENTOT = DENTOT + MAX(0.0,DENS(IP,IS))
          ENDDO
          IF (DENTOT .LT. DENMIN) THEN
            DO IS = 1,NSPIN
              VXC(IP,IS) = 0.0
            ENDDO
            GOTO 210
          ENDIF
        ENDIF

C Find mesh indexes of neighbour points
C Note : a negative index indicates a point from the buffer region
        IF (GGA .OR. MTYPE.NE.0) THEN

C X-direction
          DO IN = -NN,NN
            J1 = MOD( I1+IN+100*NMESH(1), NMESH(1) )
            JP(1,IN) = 1+J1+NSPAN(1)*I2+NSPAN(1)*NSPAN(2)*I3
          ENDDO

C Y-direction
#ifdef MPI
          if (NMESH(2).eq.NMESHG(2)) then
#endif
            DO IN = -NN,NN
              J2 = MOD( I2+IN+100*NMESH(2), NMESH(2) )
              JP(2,IN) = 1+I1+NSPAN(1)*J2+NSPAN(1)*NSPAN(2)*I3
            ENDDO
#ifdef MPI
          else
            DO IN = -NN,NN
              J2 = I2+IN
              if (J2.lt.0) then
C Out of block to the left - negative index
                IOut = -J2
                JP(2,IN) = -(1+I1+PMesh(1,1)*I3)
                JPNN(2,IN) = J2
              elseif (J2.gt.(NMESH(2)-1)) then
C Out of block to the right - negative index
                IOut = J2 - NMESH(2) + 1
                JP(2,IN) = -(1+I1+PMesh(1,2)*I3)
                JPNN(2,IN) = IOut
              else
C In block - positive index
                JP(2,IN) = 1+I1+NSPAN(1)*J2+NSPAN(1)*NSPAN(2)*I3
              endif
            ENDDO
          endif
#endif

C       Z-direction
#ifdef MPI
          if (NMESH(3).eq.NMESHG(3)) then
#endif
            DO IN = -NN,NN
              J3 = MOD( I3+IN+100*NMESH(3), NMESH(3) )
              JP(3,IN) = 1+I1+NSPAN(1)*I2+NSPAN(1)*NSPAN(2)*J3
            ENDDO
#ifdef MPI
          else
            DO IN = -NN,NN
              J3 = I3+IN
              if (J3.lt.0) then
C Out of block to the left - negative index
                IOut = -J3
                JP(3,IN) = -(1+I1+PMesh(1,3)*I2)
                JPNN(3,IN) = J3
              elseif (J3.gt.(NMESH(3)-1)) then
C Out of block to the right - negative index
                IOut = J3 - NMESH(3) + 1
                JP(3,IN) = -(1+I1+PMesh(1,4)*I2)
                JPNN(3,IN) = IOut
              else
C In block - positive index
                JP(3,IN) = 1+I1+NSPAN(1)*I2+NSPAN(1)*NSPAN(2)*J3
              endif
            ENDDO
          endif
#endif

        ENDIF

C       Find Jacobian matrix dx/dmesh for adaptative mesh
        IF (MTYPE .NE. 0) THEN

C         Find dx/dmesh
          DO 110 IC = 1,3
            DO 100 IX = 1,3
              DXDM(IX,IC) = 0.0d0
              DO 90 IN = -NN,NN
                IF (MTYPE .EQ. 1) THEN
                  DXDM(IX,IC) = DXDM(IX,IC) +
     .                          XMESH(IX,JP(IC,IN)) * DGDM(IN)
                ELSE
                  DXDM(IX,IC) = DXDM(IX,IC) +
     .                   ( CELL(IX,1) * XMESH(1,JP(IC,IN)) +
     .                     CELL(IX,2) * XMESH(2,JP(IC,IN)) +
     .                     CELL(IX,3) * XMESH(3,JP(IC,IN)) ) * DGDM(IN)
                ENDIF
   90         CONTINUE
  100       CONTINUE
  110     CONTINUE

C         Find inverse of matrix dx/dmesh
          CALL RECLAT( DXDM, DMDX, 0 )

C         Find differential of volume = determinant of Jacobian matrix
          DVOL = VOLCEL( DXDM )
          IF (GGA) AUX(IP) = DVOL

C         Find the weights for the derivative d(gradF(i))/d(F(j)), of
C         the gradient at point i with respect to the value at point j
          IF (GGA) THEN
            DO 140 IN = -NN,NN
              DO 130 IC = 1,3
                DO 120 IX = 1,3
                  DGIDFJ(IX,IC,IN) = DMDX(IX,IC) * DGDM(IN)
  120           CONTINUE
  130         CONTINUE
  140       CONTINUE
          ENDIF

        ENDIF

C  Find density and gradient of density at this point
        DO IS = 1,NSPIN
          D(IS) = DENS(IP,IS)
        ENDDO
C       Test to ensure that densities are always > 0 added to 
C       avoid floating point errors in ggaxc. JDG & JMS
        DO IS = 1,MIN(NSPIN,2)
          D(IS) = MAX(0.D0,D(IS))
*         D(IS) = max(DENMIN,D(IS))
        ENDDO
        IF (GGA) THEN
#ifdef MPI
          if (Nodes.eq.1) then
#endif
            DO IS = 1,NSPIN
              DO IX = 1,3
                GD(IX,IS) = 0.0d0
                DO IN = -NN,NN
                  GD(IX,IS) = GD(IX,IS) +
     .                      DGIDFJ(IX,1,IN) * DENS(JP(1,IN),IS) +
     .                      DGIDFJ(IX,2,IN) * DENS(JP(2,IN),IS) +
     .                      DGIDFJ(IX,3,IN) * DENS(JP(3,IN),IS)
                ENDDO
              ENDDO
            ENDDO
#ifdef MPI
          else
            DO IS = 1,NSPIN
              DO IX = 1,3
                GD(IX,IS) = 0.0d0
                DO IN = -NN,NN
                  GD(IX,IS) = GD(IX,IS) +
     .              DGIDFJ(IX,1,IN) * DENS(JP(1,IN),IS)
                ENDDO
                DO IN = -NN,NN
                  if (JP(2,IN).gt.0) then
                    GD(IX,IS) = GD(IX,IS) +
     .                DGIDFJ(IX,2,IN) * DENS(JP(2,IN),IS)
                  else
                    INN = JPNN(2,IN)
                    if (INN.lt.0) then
                      GD(IX,IS) = GD(IX,IS) +
     .                  DGIDFJ(IX,2,IN) * bdens(-JP(2,IN),-INN,IS)
                    else
                      GD(IX,IS) = GD(IX,IS) +
     .                  DGIDFJ(IX,2,IN) * bdens(-JP(2,IN),NN+INN,IS)
                    endif
                  endif
                ENDDO
                DO IN = -NN,NN
                  if (JP(3,IN).gt.0) then
                    GD(IX,IS) = GD(IX,IS) +
     .                DGIDFJ(IX,3,IN) * DENS(JP(3,IN),IS)
                  else
                    INN = JPNN(3,IN)
                    if (INN.lt.0) then
                      GD(IX,IS) = GD(IX,IS) +
     .                  DGIDFJ(IX,3,IN) * bdens(-JP(3,IN),2*NN-INN,IS)
                    else
                      GD(IX,IS) = GD(IX,IS) +
     .                  DGIDFJ(IX,3,IN) * bdens(-JP(3,IN),3*NN+INN,IS)
                    endif
                  endif
                ENDDO
              ENDDO
            ENDDO
          endif
#endif
        ENDIF

C       Find exchange and correlation energy densities and their 
C       derivatives with respect to density and density gradient
        IF (GGA) THEN
          CALL GGAXC( AUTHOR, IREL, NSPIN, D, GD,
     .                EPSX, EPSC, DEXDD, DECDD, DEXDGD, DECDGD )
        ELSE
          CALL LDAXC( AUTHOR, IREL, NSPIN, D, EPSX, EPSC, DEXDD, DECDD,
     .                DVXDN, DVCDN )
        ENDIF

C       Add contributions to exchange-correlation energy and its
C       derivatives with respect to density at all points
        DO 170 IS = 1,MIN(NSPIN,2)
          EX = EX + DVOL * D(IS) * EPSX
          EC = EC + DVOL * D(IS) * EPSC
          DX = DX + DVOL * D(IS) * EPSX
          DC = DC + DVOL * D(IS) * EPSC
  170   CONTINUE
        DO 200 IS = 1,NSPIN
          DX = DX - DVOL * D(IS) * DEXDD(IS)
          DC = DC - DVOL * D(IS) * DECDD(IS)
          IF (GGA) THEN
            VXC(IP,IS) = VXC(IP,IS) + DVOL * ( DEXDD(IS) + DECDD(IS) )
#ifdef MPI
            if (Nodes.eq.1) then
#endif
              DO IN = -NN,NN
                DO IC = 1,3
                  DO IX = 1,3
                    DX = DX - DVOL * DENS(JP(IC,IN),IS) *
     .                      DEXDGD(IX,IS) * DGIDFJ(IX,IC,IN)
                    DC = DC - DVOL * DENS(JP(IC,IN),IS) *
     .                      DECDGD(IX,IS) * DGIDFJ(IX,IC,IN)
                    VXC(JP(IC,IN),IS) = VXC(JP(IC,IN),IS) + DVOL *
     .                (DEXDGD(IX,IS)+DECDGD(IX,IS))*DGIDFJ(IX,IC,IN)
                  ENDDO
                ENDDO
              ENDDO
#ifdef MPI
            else
              DO IN = -NN,NN

C X-direction
                DO IX = 1,3
                  DX = DX - DVOL * DENS(JP(1,IN),IS) *
     .                    DEXDGD(IX,IS) * DGIDFJ(IX,1,IN)
                  DC = DC - DVOL * DENS(JP(1,IN),IS) *
     .                    DECDGD(IX,IS) * DGIDFJ(IX,1,IN)
                  VXC(JP(1,IN),IS) = VXC(JP(1,IN),IS) + DVOL *
     .              (DEXDGD(IX,IS)+DECDGD(IX,IS))*DGIDFJ(IX,1,IN)
                ENDDO

C Y-direction
                if (JP(2,IN).gt.0) then
                  DO IX = 1,3
                    DX = DX - DVOL * DENS(JP(2,IN),IS) *
     .                    DEXDGD(IX,IS) * DGIDFJ(IX,2,IN)
                    DC = DC - DVOL * DENS(JP(2,IN),IS) *
     .                    DECDGD(IX,IS) * DGIDFJ(IX,2,IN)
                    VXC(JP(2,IN),IS) = VXC(JP(2,IN),IS) + DVOL *
     .                (DEXDGD(IX,IS)+DECDGD(IX,IS))*DGIDFJ(IX,2,IN)
                  ENDDO
                else
                  INN = JPNN(2,IN)
                  if (INN.lt.0) then
                    DO IX = 1,3
                      DX = DX - DVOL * bdens(-JP(2,IN),-INN,IS) *
     .                      DEXDGD(IX,IS) * DGIDFJ(IX,2,IN)
                      DC = DC - DVOL * bdens(-JP(2,IN),-INN,IS) *
     .                      DECDGD(IX,IS) * DGIDFJ(IX,2,IN)
                      bvxc(-JP(2,IN),-INN,IS) = DVOL*(DEXDGD(IX,IS)+ 
     .                      DECDGD(IX,IS))*DGIDFJ(IX,2,IN) +
     .                      bvxc(-JP(2,IN),-INN,IS)
                    ENDDO
                  else
                    DO IX = 1,3
                      DX = DX - DVOL * bdens(-JP(2,IN),NN+INN,IS) *
     .                      DEXDGD(IX,IS) * DGIDFJ(IX,2,IN)
                      DC = DC - DVOL * bdens(-JP(2,IN),NN+INN,IS) *
     .                      DECDGD(IX,IS) * DGIDFJ(IX,2,IN)
                      bvxc(-JP(2,IN),NN+INN,IS) = DVOL*(DEXDGD(IX,IS)+ 
     .                      DECDGD(IX,IS))*DGIDFJ(IX,2,IN) +
     .                bvxc(-JP(2,IN),NN+INN,IS)
                    ENDDO
                  endif
                endif

C Z-direction
                if (JP(3,IN).gt.0) then
                  DO IX = 1,3
                    DX = DX - DVOL * DENS(JP(3,IN),IS) *
     .                    DEXDGD(IX,IS) * DGIDFJ(IX,3,IN)
                    DC = DC - DVOL * DENS(JP(3,IN),IS) *
     .                    DECDGD(IX,IS) * DGIDFJ(IX,3,IN)
                    VXC(JP(3,IN),IS) = VXC(JP(3,IN),IS) + DVOL *
     .                (DEXDGD(IX,IS)+DECDGD(IX,IS))*DGIDFJ(IX,3,IN)
                  ENDDO
                else
                  INN = JPNN(3,IN)
                  if (INN.lt.0) then
                    DO IX = 1,3
                      DX = DX - DVOL * bdens(-JP(3,IN),2*NN-INN,IS) *
     .                      DEXDGD(IX,IS) * DGIDFJ(IX,3,IN)
                      DC = DC - DVOL * bdens(-JP(3,IN),2*NN-INN,IS) *
     .                      DECDGD(IX,IS) * DGIDFJ(IX,3,IN)
                      bvxc(-JP(3,IN),2*NN-INN,IS) = DVOL*(DEXDGD(IX,IS)+ 
     .                      DECDGD(IX,IS)) * DGIDFJ(IX,3,IN) + 
     .                      bvxc(-JP(3,IN),2*NN-INN,IS)
                    ENDDO
                  else
                    DO IX = 1,3
                      DX = DX - DVOL * bdens(-JP(3,IN),3*NN+INN,IS) *
     .                      DEXDGD(IX,IS) * DGIDFJ(IX,3,IN)
                      DC = DC - DVOL * bdens(-JP(3,IN),3*NN+INN,IS) *
     .                      DECDGD(IX,IS) * DGIDFJ(IX,3,IN)
                      bvxc(-JP(3,IN),3*NN+INN,IS) = DVOL*(DEXDGD(IX,IS)+ 
     .                      DECDGD(IX,IS)) * DGIDFJ(IX,3,IN) +
     .                bvxc(-JP(3,IN),3*NN+INN,IS)
                    ENDDO
                  endif
                endif

              ENDDO
            endif

#endif
          ELSE
            VXC(IP,IS) = DEXDD(IS) + DECDD(IS)
            IF (IDER .EQ. 1) THEN
              DO IS2 = 1, NSPIN
                DVXCDN(IP,IS2,IS) = DVXDN(IS2,IS) + DVCDN(IS2,IS)
              ENDDO
            ENDIF
          ENDIF
  200   CONTINUE

C       Add contribution to stress due to change in gradient of density
C       originated by the deformation of the mesh with strain
        IF (GGA) THEN
          DO 206 JX = 1,3
            DO 204 IX = 1,3
              DO 202 IS = 1,NSPIN
                STRESS(IX,JX) = STRESS(IX,JX) - DVOL * GD(IX,IS) *
     .                           ( DEXDGD(JX,IS) + DECDGD(JX,IS) )
  202         CONTINUE
  204       CONTINUE
  206     CONTINUE
        ENDIF

  210 CONTINUE
  220 CONTINUE
  230 CONTINUE


#ifdef MPI
C Return buffer region contributions to VXC to their correct nodes
  992 if (GGA.and.Nodes.gt.1) then
        if (NoLocalPoints) goto 993

C Exchange density information in the Z-direction
        if (NMESH(3).ne.NMESHG(3)) then

C Loop over points of buffer region
          do in = 1,nn

C Loop over spins
            do is = 1,nspin

C Post VXC receive
              NSwap = NSPAN(1)*NSPAN(2)
              Tag = 8*in + 2*is
              call MPI_IRecv(bdens(1,2*nn+in,is),NSwap,DAT_single,
     $             Pleft3,Tag+1,MPI_Comm_World,Requests(1),MPIerror)
              call MPI_IRecv(bdens(1,3*nn+in,is),NSwap,DAT_single,
     $             Pright3,Tag+2,MPI_Comm_World,Requests(2),MPIerror)

C Post VXC send
              NSwap = PMesh(1,4)*PMesh(2,4)
              call MPI_ISend(bvxc(1,3*nn+in,is),NSwap,DAT_single,
     .          Pright3,Tag+1,MPI_Comm_World,Requests(3),MPIerror)
              NSwap = PMesh(1,3)*PMesh(2,3)
              call MPI_ISend(bvxc(1,2*nn+in,is),NSwap,DAT_single,
     .          Pleft3,Tag+2,MPI_Comm_World,Requests(4),MPIerror)

C Wait for transfers to complete
              call MPI_WaitAll(4,Requests,Statuss,MPIerror)

            enddo

          enddo

        endif

C Exchange density information in the Y-direction
        if (NMESH(2).ne.NMESHG(2)) then

C Loop over points of buffer region
          do in = 1,nn

C Loop over spins
            do is = 1,nspin

C Post VXC receives
              NSwap = NSPAN(1)*NSPAN(3)
              Tag = 8*in + 2*is
              call MPI_IRecv(bdens(1,in,is),NSwap,DAT_single,
     .          Pleft2,Tag+1,MPI_Comm_World,Requests(1),MPIerror)
              call MPI_IRecv(bdens(1,nn+in,is),NSwap,DAT_single,
     .          Pright2,Tag+2,MPI_Comm_World,Requests(2),MPIerror)

C Post VXC sends
              NSwap = PMesh(1,2)*PMesh(3,2)
              call MPI_ISend(bvxc(1,nn+in,is),NSwap,DAT_single,
     .          Pright2,Tag+1,MPI_Comm_World,Requests(3),MPIerror)
              NSwap = PMesh(1,1)*PMesh(3,1)
              call MPI_ISend(bvxc(1,in,is),NSwap,DAT_single,
     .          Pleft2,Tag+2,MPI_Comm_World,Requests(4),MPIerror)

C Wait for transfers to complete
              call MPI_WaitAll(4,Requests,Statuss,MPIerror)

            enddo

          enddo

        endif

C Wait for end of transfers
  993   call MPI_Barrier(MPI_Comm_World,MPIerror)
        if (NoLocalPoints) goto 999

C Return buffer VXC values to main array
        if (NMESH(2).ne.NMESHG(2)) then
          do is = 1,nspin
            do in = 1,nn
              IDen1 = (in-1)*NSPAN(1) 
              IDen2 = (NMESH(2)-in)*NSPAN(1) 
              do ix = 1,NMESH(3)
                do ic = 1,NMESH(1)
                  VXC(IDen1+ic,is) = VXC(IDen1+ic,is) +
     .              bdens((ix-1)*NSPAN(1)+ic,in,is)
                  VXC(IDen2+ic,is) = VXC(IDen2+ic,is) +
     .              bdens((ix-1)*NSPAN(1)+ic,nn+in,is)
                enddo
                IDen1 = IDen1 + NSPAN(1)*NSPAN(2)
                IDen2 = IDen2 + NSPAN(1)*NSPAN(2)
              enddo
            enddo
          enddo
        endif

        if (NMESH(3).ne.NMESHG(3)) then
          do is = 1,nspin
            do in = 1,nn
              IDen1 = NSPAN(1)*NSPAN(2)*(in-1)
              IDen2 = NSPAN(1)*NSPAN(2)*(NMESH(3)-in)
              do ix = 1,NSPAN(1)*NMESH(2)
                VXC(IDen1+ix,is) = VXC(IDen1+ix,is) +
     .            bdens(ix,2*nn+in,is)
                VXC(IDen2+ix,is) = VXC(IDen2+ix,is) +
     .            bdens(ix,3*nn+in,is)
              enddo
            enddo
          enddo
        endif

      endif
#endif

C Divide by volume element to obtain the potential (per electron)
      IF (GGA) THEN
        DO 270 IS = 1,NSPIN
          DO 260 I3 = 0,NMESH(3)-1
          DO 250 I2 = 0,NMESH(2)-1
          DO 240 I1 = 0,NMESH(1)-1
            IP = 1 + I1 + NSPAN(1) * I2 + NSPAN(1) * NSPAN(2) * I3
            VXC(IP,IS) = VXC(IP,IS) / DVOL
  240     CONTINUE
  250     CONTINUE
  260     CONTINUE
  270   CONTINUE
      ENDIF

C Add contribution to stress from the change of volume with strain and
C divide by volume to get correct stress definition (dE/dStrain)/Vol
      VOLUME = VOLCEL( CELL )
      DO 274 JX = 1,3
        STRESS(JX,JX) = STRESS(JX,JX) + EX + EC
        DO 272 IX = 1,3
          STRESS(IX,JX) = STRESS(IX,JX) / VOLUME
  272   CONTINUE
  274 CONTINUE
      
C Divide by energy unit
      EX = EX / EUNIT
      EC = EC / EUNIT
      DX = DX / EUNIT
      DC = DC / EUNIT
      DO 310 IS = 1,NSPIN
        DO 300 I3 = 0,NMESH(3)-1
        DO 290 I2 = 0,NMESH(2)-1
        DO 280 I1 = 0,NMESH(1)-1
          IP = 1 + I1 + NSPAN(1) * I2 + NSPAN(1) * NSPAN(2) * I3
          VXC(IP,IS) = VXC(IP,IS) / EUNIT
  280   CONTINUE
  290   CONTINUE
  300   CONTINUE
  310 CONTINUE
      DO 330 JX = 1,3
        DO 320 IX = 1,3
          STRESSL(IX,JX) = STRESSL(IX,JX) + (STRESS(IX,JX) / EUNIT)
  320   CONTINUE
  330 CONTINUE
      
      DO IS = 1, NSPIN
        DO IS2 = 1, NSPIN
          DO IP = 1, IAUXDER
            DVXCDN(IP,IS2,IS) = DVXCDN(IP,IS2,IS) / EUNIT
          ENDDO
        ENDDO
      ENDDO

C Re-entry point for nodes with no mesh points stored locally
  999 continue

C Collect together quantities
#ifdef MPI
      Buffer(1) = EX
      Buffer(2) = EC
      Buffer(3) = DX
      Buffer(4) = DC
      call MPI_AllReduce(Buffer,Buffer2,4,DAT_double,MPI_sum,
     .  MPI_Comm_World,MPIerror)
      EX = Buffer2(1)
      EC = Buffer2(2)
      DX = Buffer2(3)
      DC = Buffer2(4)

C Free memory
      if (GGA.and.Nodes.gt.1.and..not.NoLocalPoints) then
        call memory('D','S',size(bdens),'cellxc')
        deallocate(bdens)
        call memory('D','S',size(bvxc),'cellxc')
        deallocate(bvxc)
        call memory('D','I',size(PMesh),'cellxc')
        deallocate(PMesh)
      endif
      if (GGA.and.Nodes.gt.1) then
        call memory('D','I',size(Reqm),'cellxc')
        deallocate(Reqm)
        call memory('D','I',size(Statm),'cellxc')
        deallocate(Statm)
      endif
#endif

C Deallocate local memory
      call memory('D','D',size(D),'cellxc')
      deallocate(D)
      call memory('D','D',size(DECDD),'cellxc')
      deallocate(DECDD)
      call memory('D','D',size(DEXDD),'cellxc')
      deallocate(DEXDD)
      call memory('D','D',size(DECDGD),'cellxc')
      deallocate(DECDGD)
      call memory('D','D',size(DEXDGD),'cellxc')
      deallocate(DEXDGD)
      call memory('D','D',size(GD),'cellxc')
      deallocate(GD)
      call memory('D','D',size(DVXDN),'cellxc')
      deallocate(DVXDN)
      call memory('D','D',size(DVCDN),'cellxc')
      deallocate(DVCDN)

C Stop time counter
      CALL TIMER( 'CELLXC', 2 )

      END
