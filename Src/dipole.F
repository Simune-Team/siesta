      SUBROUTINE DIPOLE( CELL, NTM, M1, M2, M3, NSM, RHO, X0, DIPOL )
C ********************************************************************
C Finds the electric dipole
C Written by J.M.Soler. July 1997.
C Modified for distributed rho matrix using a 2D grid of processors.
C Routine now is based on intrinsic structure of grid distribution
C in order to calculate position of grid points in space for local
C matrix. J.D.Gale March 1999.
C *********** INPUT **************************************************
C REAL*8  CELL(3,3)     : Unit cell vectors
C INTEGER NTM(3)        : Global number of divisions of each lattice vector
C INTEGER M1,M2,M3      : Local number of divisions of each lattice vector
C INTEGER NSM           : Number of sub-points for each mesh point
C REAL    RHO(N1,N2,N3) : Minus neutral charge density at mesh points
C                         Notice single precision in this version
C REAL*8  X0(3)         : Origin in cartesian coordinates
C                         (center of molecule)
C *********** OUTPUT *************************************************
C REAL*8 DIPOL(3)   : Electric dipole
C *********** UNITS **************************************************
C CELL  in atomic units (Bohr)
C RHO   in atomic units (electrons/Bohr**3)
C X0    in atomic units (Bohr)
C DIPOL in atomic units (electrons*Bohr)
C ********************************************************************
C
C  Modules
C
      use precision
      use parallel
#ifdef MPI
      use mpi
#endif

      IMPLICIT          NONE

      INTEGER           M1, M2, M3, NTM(3), NSM
      REAL              RHO(M1,M2,M3)
      DOUBLE PRECISION  CELL(3,3), DIPOL(3), VOLCEL, X0(3)
      EXTERNAL          RECLAT, VOLCEL

C Internal variables and arrays
      INTEGER           I, I1, I2, I3, IX, Node, Nodes,
     .                  ProcessorZ, Py, Pz, I20, I30,
     .                  MG1, MG2, MG3, Yoffset, Zoffset,
     .                  BlockSizeY, BlockSizeZ, NRemY, NRemZ,
     .                  MeshNsm(3)
      DOUBLE PRECISION  D(3), DVOL, DX, RCELL(3,3), X0L(3)
#ifdef MPI
      integer           MPIerror
#endif

C Get the Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
#else
      Node = 0
      Nodes = 1
#endif

C Assign local variables
      MG1 = NTM(1)
      MG2 = NTM(2)
      MG3 = NTM(3)

C Find volume element
      DVOL = VOLCEL( CELL ) / (MG1*MG2*MG3)

C Find reciprocal cell vectors (without the factor 2*pi)
      CALL RECLAT( CELL, RCELL, 0 )

C Find origin in lattice coordinates
      DO 10 I = 1,3
        X0L(I)= X0(1)*RCELL(1,I) + X0(2)*RCELL(2,I) + X0(3)*RCELL(3,I)
   10 CONTINUE

C Initialize dipole
      DIPOL(1) = 0.D0
      DIPOL(2) = 0.D0
      DIPOL(3) = 0.D0

C Check that ProcessorY is a factor of the number of processors
      if (mod(Nodes,ProcessorY).gt.0) then
        write(6,'(''ERROR: ProcessorY must be a factor of the'',
     .    '' number of processors!'')')
        stop
      endif
      ProcessorZ = Nodes/ProcessorY

C Calculate coarse mesh sizes
      MeshNsm(1) = MG1/NSM
      MeshNsm(2) = MG2/NSM
      MeshNsm(3) = MG3/NSM

C Calculate blocking sizes
      BlockSizeY = (MeshNsm(2)/ProcessorY)*NSM
      NRemY = (MG2 - ProcessorY*BlockSizeY)/NSM
      BlockSizeZ = (MeshNsm(3)/ProcessorZ)*NSM
      NRemZ = (MG3 - ProcessorZ*BlockSizeZ)/NSM

C Calculate coordinates of current node in processor grid
      Py = (Node/ProcessorZ)+1
      Pz = Node - (Py - 1)*ProcessorZ + 1

C Calculate starting point for grid
      Yoffset = (Py-1)*BlockSizeY + NSM*min(Py-1,NRemY)
      Zoffset = (Pz-1)*BlockSizeZ + NSM*min(Pz-1,NRemZ)

C Find dipole by direct integration allowing for block distributed
C structure of rho
      I30 = Zoffset
      DO 80 I3 = 1,M3
        I30 = I30 + 1
        D(3) = DBLE(I30-1) / DBLE(MG3) - X0L(3)
        IF (D(3) .LT. -0.5D0) D(3) = D(3) + 1.D0
        IF (D(3) .GT. +0.5D0) D(3) = D(3) - 1.D0
        I20 = Yoffset
        DO 70 I2 = 1,M2
          I20 = I20 + 1
          D(2) = DBLE(I20-1) / DBLE(MG2) - X0L(2)
          IF (D(2) .LT. -0.5D0) D(2) = D(2) + 1.D0
          IF (D(2) .GT. +0.5D0) D(2) = D(2) - 1.D0
          DO 60 I1 = 1,M1
            D(1) = DBLE(I1-1) / DBLE(MG1) - X0L(1)
            IF (D(1) .LT. -0.5D0) D(1) = D(1) + 1.D0
            IF (D(1) .GT. +0.5D0) D(1) = D(1) - 1.D0
            DO 40 IX = 1,3
              DX = CELL(IX,1)*D(1) + CELL(IX,2)*D(2) + CELL(IX,3)*D(3)
              DIPOL(IX) = DIPOL(IX) - DX * RHO(I1,I2,I3) * DVOL
   40       CONTINUE
   60     CONTINUE
   70   CONTINUE
   80 CONTINUE

      END

