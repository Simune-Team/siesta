C $Id: dipole.f,v 1.3 1999/01/31 10:53:53 emilio Exp $

      SUBROUTINE DIPOLE( CELL, M1, M2, M3, RHO, X0, DIPOL )
C ********************************************************************
C Finds the electric dipole
C Written by J.M.Soler. July 1997.
C *********** INPUT **************************************************
C REAL*8  CELL(3,3)     : Unit cell vectors
C INTEGER M1,M2,M3      : Number of divisions of each lattice vector
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

      IMPLICIT          NONE
      INTEGER           M1, M2, M3
      REAL              RHO(M1,M2,M3)
      DOUBLE PRECISION  CELL(3,3), DIPOL(3), VOLCEL, X0(3)
      EXTERNAL          RECLAT, VOLCEL

C Internal variables and arrays
      INTEGER           I, I1, I2, I3, IX
      DOUBLE PRECISION  D(3), DVOL, DX, RCELL(3,3), X0L(3)

C Find volume element
      DVOL = VOLCEL( CELL ) / (M1*M2*M3)

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

C Find dipole by direct integration
      DO 80 I3 = 1,M3
      DO 70 I2 = 1,M2
      DO 60 I1 = 1,M1
        D(1) = DBLE(I1-1) / DBLE(M1) - X0L(1)
        D(2) = DBLE(I2-1) / DBLE(M2) - X0L(2)
        D(3) = DBLE(I3-1) / DBLE(M3) - X0L(3)
        DO 20 I = 1,3
          IF (D(I) .LT. -0.5D0) D(I) = D(I) + 1.D0
          IF (D(I) .GT. +0.5D0) D(I) = D(I) - 1.D0
   20   CONTINUE
        DO 40 IX = 1,3
          DX = CELL(IX,1)*D(1) + CELL(IX,2)*D(2) + CELL(IX,3)*D(3)
          DIPOL(IX) = DIPOL(IX) - DX * RHO(I1,I2,I3) * DVOL
   40   CONTINUE
   60 CONTINUE
   70 CONTINUE
   80 CONTINUE

      END

