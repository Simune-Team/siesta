C $Id: shaper.f,v 1.8 1999/05/05 17:25:34 emilio Exp $

      SUBROUTINE SHAPER( CELL, NA, ISA, XA, SHAPE, NV, VECS )

C **********************************************************************
C Finds the topology of a system (molecule, chain, slab or bulk)
C Written by J.M.Soler. July 1997, Feb. 1998.
C ************ INPUT ***************************************************
C REAL*8  CELL(3,3) : Lattice vectors CELL(Ixyz,Ivector)
C INTEGER NA        : Number of atoms
C INTEGER ISA(NA)   : Species index of each atom
C REAL*8  XA(3,NA)  : Cartesian atomic coordinates
C ************ OUTPUT **************************************************
C CHARACTER SHAPE*8 : atom|molecule|chain|slab|bulk
C INTEGER NV        : Number of linearly independent vectors
C                     defining the periodic bulk directions:
C                         NV=0 => Isolated atom or molecule
C                            1 => Chain
C                            2 => Slab
C                            3 => Bulk condensed system
C REAL*8  VECS(3,NV) : vectors defining the periodic bulk directions
C ************ ROUTINES REQUIRED ***************************************
C REAL*8 RCUT(IS,IO) : User defined function which will be called with
C                      IO=0 and must return the 'interaction radius'
C                      of an atom of species IS. Two atoms interact
C                      if their interaction spheres overlap, ie if
C                      RIJ < RCUT(ISA(I),0) + RCUT(ISA(J),0)
C ************ UNITS ***************************************************
C CELL and XA must be in the same units
C ************ BEHAVIOR ************************************************
C Assumes that all the atoms in one unit cell are connected, i.e. it
C will fail if one atom evaporates from one slab and reaches another
C image of the slab, reporting the system as bulk.
C **********************************************************************

      IMPLICIT          NONE
      INTEGER           NA, ISA(NA), NV
      DOUBLE PRECISION  CELL(3,3), XA(3,NA), RCUT, VECS(3,3)
      EXTERNAL          CHKDIM, NEIGHB, LIVEC, RCUT
      CHARACTER         SHAPE*(*)

C Internal variables and arrays
      INTEGER MAXNA
      PARAMETER (MAXNA = 1000)
      INTEGER          IA, IN, IS, JA, JAN(MAXNA), JS, NNA
      DOUBLE PRECISION RI, RIJ, RJ, RMAX, R2IJ(MAXNA),
     .                 XIJ(3,MAXNA), XXJ(3)

C Find maximum interaction range
      RMAX = 0.0D0
      DO 10 IA = 1,NA
        IS = ISA(IA)
        RMAX = MAX( RMAX, RCUT(IS,0) )
   10 CONTINUE

C Initialize neighbour-locater routine
      NNA = MAXNA
      CALL NEIGHB( CELL, 2*RMAX, NA, XA, 0, 0, NNA, JAN, XIJ, R2IJ )

C Main loop
      NV = 0
      DO 90 IA = 1,NA
        IS = ISA(IA)
        RI = RCUT(IS,0)

C       Find neighbours of atom IA
        NNA = MAXNA
        CALL NEIGHB( CELL, RI+RMAX, NA, XA, IA, 0,
     .               NNA, JAN, XIJ, R2IJ )
        CALL CHKDIM( 'SHAPER', 'MAXNA', MAXNA, NNA, 1 )

        DO 50 IN = 1,NNA
          JA = JAN(IN)
          JS = ISA(JA)
          RJ = RCUT(JS,0)
          RIJ = SQRT(R2IJ(IN))

C         Check if IA and JA interact
          IF (RIJ .LT. RI+RJ) THEN

C           Find vector between two images of atom JA
C           (we assume that all atoms in one unit cell are connected)
            XXJ(1) = XA(1,IA) + XIJ(1,IN) - XA(1,JA)
            XXJ(2) = XA(2,IA) + XIJ(2,IN) - XA(2,JA)
            XXJ(3) = XA(3,IA) + XIJ(3,IN) - XA(3,JA)

C           Add to set of linearly independent vectors
            CALL LIVEC( XXJ, NV, VECS )

            IF (NV .EQ. 3) GOTO 999
          ENDIF
   50   CONTINUE
   90 CONTINUE

C Exit point
  999 CONTINUE
      IF (NV.EQ.0 .AND. NA.EQ.1) SHAPE = 'atom'
      IF (NV.EQ.0 .AND. NA.GT.1) SHAPE = 'molecule'
      IF (NV.EQ.1) SHAPE = 'chain'
      IF (NV.EQ.2) SHAPE = 'slab'
      IF (NV.EQ.3) SHAPE = 'bulk'
      END



      SUBROUTINE LIVEC( X, NV, V )

C Adds vector X to set of NV linearly-independent vectors V (if true)
C Written by J.M.Soler. July 1997.

      IMPLICIT          NONE
      INTEGER           LV, NV
      DOUBLE PRECISION  TOL, V(3,3), VV(3), VV2, VVX, X(3), X2

      PARAMETER  ( TOL = 1.D-6 ) 

      LV = NV
      IF (NV .EQ. 0) THEN
        X2 = X(1)**2 + X(2)**2 + X(3)**2
        IF (X2 .GT. TOL) NV = 1
      ELSEIF (NV .EQ. 1) THEN
        VV(1) = V(2,1)*X(3) - V(3,1)*X(2)
        VV(2) = V(3,1)*X(1) - V(1,1)*X(3)
        VV(3) = V(1,1)*X(2) - V(2,1)*X(1)
        VV2 = VV(1)**2 + VV(2)**2 + VV(3)**2
        IF (VV2 .GT. TOL) NV = 2
      ELSEIF (NV .EQ. 2) THEN
        VV(1) = V(2,1)*V(3,2) - V(3,1)*V(2,2)
        VV(2) = V(3,1)*V(1,2) - V(1,1)*V(3,2)
        VV(3) = V(1,1)*V(2,2) - V(2,1)*V(1,2)
        VVX = VV(1)*X(1) + VV(2)*X(2) + VV(3)*X(3)
        IF (VVX .GT. TOL) NV = 3
      ENDIF

      IF (NV .GT. LV) THEN
        V(1,NV) = X(1)
        V(2,NV) = X(2)
        V(3,NV) = X(3)
      ENDIF

      END








