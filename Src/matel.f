C $Id: matel.f,v 1.13 1999/05/05 17:25:33 emilio Exp $

      SUBROUTINE MATEL( OPERAT, IS1, IS2, IO1, IO2, R12, S12, DSDR )
C *******************************************************************
C Finds the overlap or laplacian matrix elements between two atomic
C orbitals. Written by J.M.Soler. April 1995.
C ************************* INPUT ***********************************
C CHARACTER OPERAT : Operator to be used: 'S' => Unity (overlap)
C                                         'T' => -Laplacian
C INTEGER   IS1    : 		Species index of 1st orbital atom
C INTEGER   IS2    : 		Species index of 2nd orbital atom
C                    Allowed range of IS1, IS2: (1, MAXS), where
C                    MAXS and MAXL (below) are internal parameters.
C INTEGER   IO1    : Orbital index of 1st orbital (within atom)
C INTEGER   IO2    : Orbital index of 2nd orbital (within atom)
C                    Allowed range of IO1, IO2: (-MAXL**2, MAXL**2)
C                    Indexes IS1,IS2,IO1,IO2 are used to call routines
C                    RCUT and PHIATM (see below).
C REAL*8    R12(3) : Vector from first to second atom
C ************************* OUTPUT **********************************
C REAL*8 S12      : Matrix element between orbitals.
C REAL*8 DSDR(3)  : Derivative (gradient) of S12 with respect to R12.
C ************************* ROUTINES CALLED *************************
C The following functions must exist:
C
C REAL*8 FUNCTION RCUT(IS,IO)
C   Returns cutoff radius of orbitals and KB projectors.
C Input:
C     INTEGER IS : Species index
C     INTEGER IO : Orbital index
C
C SUBROUTINE PHIATM(IS,IO,R,PHI,GRPHI)
C    Finds values and gradients of:
C    a) basis orbitals (IO > 0)
C    b) KB proyectors  (IO < 0)
C    c) Local part of screened pseudopotential (IO = 0) ( b) and c) are
C       not required if MATEL is called only with IO > 0 ) 
C Input:
C   INTEGER IS   : Species index
C   INTEGER IO   : Orbital index
C   REAL*8  R(3) : Position with respect to atom
C Output:
C   REAL*8  PHI      : Value of orbital or KB projector at point R
C   REAL*8  GRPHI(3) : Gradient of PHI at point R  
C
C INTEGER FUNCTION LOMAXFIS(IS)
C    Returns the maximum angular momentum of orbitals 
C Input:
C     INTEGER IS : Species index
C 
C INTEGER FUNCTION LMXKBFIS(IS)
C    Returns the maximum angular momentum of KB projectors
C Input:
C     INTEGER IS : Species index
C
C INTEGER FUNCTION NZTFL(IS,L)
C    Returns the number of different basis functions with the
C    same angular momentum L.
C Input:
C     INTEGER IS : Species index
C     INTEGER L  : Angular momentum
C
C ************************* UNITS ***********************************
C Length units are arbitrary, but must be consistent in MATEL, RCUT
C   and PHIATM. The laplacian unit is (length unit)**(-2).
C ************************* BEHAVIOUR *******************************
C 1) Returns exactly zero if |R12| > RCUT(IS1,IO1) + RCUT(IS2,IO2)
C 2) If (IS1.LE.0 .OR. IS2.LE.0) all internal tables are erased for
C    reinitialization.
C *******************************************************************

C Next line is non-standard and may be suppressed -------------------
      IMPLICIT NONE
C -------------------------------------------------------------------

C Argument types and dimensions -------------------------------------
      CHARACTER         OPERAT
      INTEGER           IO1, IO2, IS1, IS2
      DOUBLE PRECISION  DSDR(3), RCUT, R12(3), S12
      EXTERNAL          PHIATM, RCUT
C -------------------------------------------------------------------

C Argument types and dimensions for entry MATEL0 --------------------
       INTEGER  NS, LOMAXFIS, LMXKBFIS, NZTFL
       EXTERNAL LOMAXFIS, LMXKBFIS, NZTFL        
C -------------------------------------------------------------------

C Internal precision parameters  ------------------------------------
C  NR is the number of radial points for matrix-element tables.
C  NQ is the number of radial points in reciprocal space.
C  Q2CUT is the required planewave cutoff for orbital expansion
C    (in Ry if lengths are in Bohr).
C  FFTOL is the tolerance for considering equal the radial part of
C    two orbitals.
C  TINY is a small number to add to a zero denominator
      INTEGER           NR, NQ
      DOUBLE PRECISION  Q2CUT, FFTOL, TINY
      PARAMETER ( NR     =  128    )
      PARAMETER ( NQ     =  1024   )
      PARAMETER ( Q2CUT  =  2.5D3  )
      PARAMETER ( FFTOL  =  1.D-8  )
      PARAMETER ( TINY   =  1.D-12 )
C -------------------------------------------------------------------

C Basic dimension parameters for internal variables -----------------
C  MAXS  : Maximum number of atomic species.
C  MAXL  : Maximum angular momentum of atomic basis orbitals and 
C          Kleinman-Bylander projectors.
C  MAXZ  : Maximum number of atomic basis orbitals per angular momentum
C          symmetry (number of 'zetas')
C  MAXY  : Maximum number of spherical-harmonic components of one
C          basis orbital or KB projector 
C  MAXR  : Maximum number of radial points.
C  MAXQ  : Maximum number of radial points in k-space.
*     INTEGER  MAXS, MAXL, MAXZ, MAXY, MAXR, MAXQ
*     PARAMETER ( MAXS =   5 )
*     PARAMETER ( MAXL =   3 )
*     PARAMETER ( MAXZ =   2 )
*     PARAMETER ( MAXY =  14 )
*     PARAMETER ( MAXR =  NR )
*     PARAMETER ( MAXQ =  NQ )
      INCLUDE 'matel.h'
C -------------------------------------------------------------------

C Derived dimension parameters --------------------------------------
C Ref: J.M.Soler notes of 4/4/96.
      INTEGER
     .  MAXF, MAXFF, MAXFFR, MAXFFY, MAXF2,
     .  MAXKB, MAXLO, MAXLKB, MAXLM, MAXLP1,
     .  MAXO, MAXSS, MAXZKB, MAXZO, MAXZZ
      PARAMETER ( MAXLO  = MAXL )
      PARAMETER ( MAXLKB = MAXL )
      PARAMETER ( MAXZO  = MAXZ )
      PARAMETER ( MAXZKB = 1 )
      PARAMETER ( MAXLP1 = MAXL + 1 )
      PARAMETER ( MAXLM  = (MAXL+1)*(MAXL+1) )
      PARAMETER ( MAXO   = MAXZO  * (MAXLO+1)  * (MAXLO+1)  )
      PARAMETER ( MAXKB  = MAXZKB * (MAXLKB+1) * (MAXLKB+1) )
      PARAMETER ( MAXSS  = MAXS * (MAXS+1) / 2 )
      PARAMETER ( MAXZZ  = MAXZO * (2*MAXZO+MAXZKB) + 2 )
      PARAMETER ( MAXF   = MAXS * (MAXO+MAXKB+1) )
      PARAMETER ( MAXFF  = MAXSS * MAXZZ * MAXLM * MAXLM )
      PARAMETER ( MAXFFR = MAXSS * MAXZZ *
     .                     (MAXL+1) * (MAXL+2) * (2*MAXL+3) / 6 )
      PARAMETER ( MAXFFY = MAXSS * MAXZZ * MAXLP1 *
     .                    (MAXLP1*MAXLP1*MAXLP1*(16*MAXLP1+15)-1) / 30 )
      PARAMETER ( MAXF2  = MAXF * MAXF * 2 )
C -------------------------------------------------------------------

C Internal variable types and dimensions ----------------------------
      INTEGER
     .  I, I1, I2, IF1, IF2, IFF, IFFR(0:2*MAXL), IFFY,
     .  ILM(MAXF+MAXY), ILMFF(MAXFFY),
     .  INDF(MAXS,-MAXKB:MAXO), INDFF(MAXF,MAXF,2), INDFFR(MAXFFY),
     .  INDFFY(0:MAXFF), IO, IOPER, IQ, IR, IS, IX,
     .  JF1, JF2, JFF, JFFR, JFFY, JLM,
     .  JO, JO1, JO2, JR, JS, 
     .  L, L1, L2, L3, LMAX, LOFILM,
     .  NF, NFF, NFFR, NFFY, NLM, NQQQ, NY, NZMAX

      DOUBLE PRECISION
     .  C, CFFR(0:2*MAXL), CPROP,
     .  DFFR0, DFFRMX, DSRDR, DYDR(3,4*MAXLM),
     .  F(0:MAXQ,MAXF+MAXY), FFQ(0:MAXQ), FFL(0:MAXQ),
     .  FFR(0:MAXR,2,MAXFFR), FFY(MAXFFY),
     .  PI, Q, QMAX, R, RI(0:MAXR), RMAX, 
     .  SR, X12(3), Y(4*MAXLM)

      LOGICAL
     .  FOUND, FRSTME, PROPOR

      EXTERNAL
     .  CHKDIM, LOFILM, PRMEM, PROPOR, RADFFT, RLYLM,
     .  SPLIN, SPLINU, TIMER, YLMYLM, YLMEXP

      SAVE
     .  F, FFR, FFY, FRSTME,
     .  I1, I2, ILM, ILMFF, INDF, INDFF, INDFFR, INDFFY,
     .  NF, NFF, NFFR, NFFY, QMAX, RMAX

      DATA
     .  I1, I2    /1, 2/,
     .  INDF      /MAXF * 0/,
     .  INDFF     /MAXF2 * 0/,
     .  INDFFY(0) /0/,
     .  NF        /0/,
     .  NFF       /0/,
     .  NFFR      /0/,
     .  NFFY      /0/,
     .  FRSTME    /.TRUE./
C -------------------------------------------------------------------

      ENTRY MATEL0( NS )
C *******************************************************************
C Checks if the array dimensions in MATEL are sufficient. If not,
C writes a new file matel.h with the correct ones, prints an error
C message and stops. 
C Modified by DSP, Aug. 1998
C *********** INPUT *************************************************
C INTEGER NS              : Number of species
C *******************************************************************

      LMAX = 0
      NZMAX = 0
      DO 320 IS = 1,NS
        LMAX = MAX( LMAX, LOMAXFIS(IS), LMXKBFIS(IS) )
        DO 310 L = 0,LOMAXFIS(IS)
          NZMAX = MAX( NZMAX, NZTFL(IS,L) )
  310   CONTINUE
  320 CONTINUE
      NY = (LMAX+1)**2

      IF (NS    .GT. MAXS .OR.
     .    LMAX  .GT. MAXL .OR.
     .    NZMAX .GT. MAXZ .OR.
     .    NY    .GT. MAXY .OR.
     .    NR    .GT. MAXR .OR.
     .    NQ    .GT. MAXQ) THEN

        OPEN( 1, FILE='matel.h', STATUS='UNKNOWN' )
        WRITE(1,'(A)') 'C Dimension parameters for MATEL'
        WRITE(1,'(6X,A)') 'INTEGER  MAXS,MAXL,MAXZ,MAXY,MAXR,MAXQ'
        WRITE(1,'(6X,A,I12,A)') 'PARAMETER ( MAXS   =', NS,   ' )'
        WRITE(1,'(6X,A,I12,A)') 'PARAMETER ( MAXL   =', LMAX, ' )'
        WRITE(1,'(6X,A,I12,A)') 'PARAMETER ( MAXZ   =', NZMAX,' )'
        WRITE(1,'(6X,A,I12,A)') 'PARAMETER ( MAXY   =', NY,   ' )'
        WRITE(1,'(6X,A,I12,A)') 'PARAMETER ( MAXR   =', NR,   ' )'
        WRITE(1,'(6X,A,I12,A)') 'PARAMETER ( MAXQ   =', NQ,   ' )'
        CLOSE( 1 )
        
        WRITE(6,'(a)') 'MATEL: Dimensions too small. RECOMPILE.'
        STOP 'MATEL: Dimensions too small. RECOMPILE.'
      ENDIF


      RETURN
c
c
c

      ENTRY MATEL( OPERAT, IS1, IS2, IO1, IO2, R12, S12, DSDR )
C *******************************************************************
C Finds the overlap or laplacian matrix elements between two atomic
C orbitals. Written by J.M.Soler. April 1995.
C ************************* INPUT ***********************************
C CHARACTER OPERAT : Operator to be used: 'S' => Unity (overlap)
C                                         'T' => -Laplacian
C INTEGER   IS1    : 		Species index of 1st orbital atom
C INTEGER   IS2    : 		Species index of 2nd orbital atom
C                    Allowed range of IS1, IS2: (1, MAXS), where
C                    MAXS and MAXL (below) are internal parameters.
C INTEGER   IO1    : Orbital index of 1st orbital (within atom)
C INTEGER   IO2    : Orbital index of 2nd orbital (within atom)
C                    Allowed range of IO1, IO2: (-MAXL**2, MAXL**2)
C                    Indexes IS1,IS2,IO1,IO2 are used to call routines
C                    RCUT and PHIATM (see below).
C REAL*8    R12(3) : Vector from first to second atom
C ************************* OUTPUT **********************************
C REAL*8 S12      : Matrix element between orbitals.
C REAL*8 DSDR(3)  : Derivative (gradient) of S12 with respect to R12.
C ************************* ROUTINES CALLED *************************
C The following functions must exist:
C
C REAL*8 FUNCTION RCUT(IS,IO)
C   Returns cutoff radius of orbitals and KB projectors.
C Input:
C     INTEGER IS : Species index
C     INTEGER IO : Orbital index
C
C SUBROUTINE PHIATM(IS,IO,R,PHI,GRPHI)
C    Finds values and gradients of:
C    a) basis orbitals (IO > 0)
C    b) KB proyectors  (IO < 0)
C    c) Local part of screened pseudopotential (IO = 0) ( b) and c) are
C       not required if MATEL is called only with IO > 0 ) 
C Input:
C   INTEGER IS   : Species index
C   INTEGER IO   : Orbital index
C   REAL*8  R(3) : Position with respect to atom
C Output:
C   REAL*8  PHI      : Value of orbital or KB projector at point R
C   REAL*8  GRPHI(3) : Gradient of PHI at point R  
C
C INTEGER FUNCTION LOMAXFIS(IS)
C    Returns the maximum angular momentum of orbitals 
C Input:
C     INTEGER IS : Species index
C 
C INTEGER FUNCTION LMXKBFIS(IS)
C    Returns the maximum angular momentum of KB projectors
C Input:
C     INTEGER IS : Species index
C
C INTEGER FUNCTION NZTFL(IS,L)
C    Returns the number of different basis functions with the
C    same angular momentum L.
C Input:
C     INTEGER IS : Species index
C     INTEGER L  : Angular momentum
C
C ************************* UNITS ***********************************
C Length units are arbitrary, but must be consistent in MATEL, RCUT
C   and PHIATM. The laplacian unit is (length unit)**(-2).
C ************************* BEHAVIOUR *******************************
C 1) Returns exactly zero if |R12| > RCUT(IS1,IO1) + RCUT(IS2,IO2)
C 2) If (IS1.LE.0 .OR. IS2.LE.0) all internal tables are erased for
C    reinitialization.
C *******************************************************************


C Start time counter ------------------------------------------------
*     CALL TIMER( 'MATEL', 1 )
C -------------------------------------------------------------------

      PI = 4.D0 * ATAN(1.D0)

C Print size of arrays ----------------------------------------------
      IF (FRSTME) THEN
        FRSTME = .FALSE.
        CALL PRMEM( 0, 'matel', 'ILM',    'I', MAXF+MAXY            )
        CALL PRMEM( 0, 'matel', 'ILMFF',  'I', MAXFFY               )
        CALL PRMEM( 0, 'matel', 'INDF',   'I', MAXF                 )
        CALL PRMEM( 0, 'matel', 'INDFF',  'I', MAXF2                )
        CALL PRMEM( 0, 'matel', 'INDFFR', 'I', MAXFFY               )
        CALL PRMEM( 0, 'matel', 'INDFFY', 'I', MAXFF+1              )
        CALL PRMEM( 0, 'matel', 'IFFR',   'I', MAXL+1               )
        CALL PRMEM( 0, 'matel', 'DYDR',   'D', 12*MAXLM             )
        CALL PRMEM( 0, 'matel', 'F',      'D', (MAXQ+1)*(MAXF+MAXY) )
        CALL PRMEM( 0, 'matel', 'FFR',    'D', (MAXR+1)*2*MAXFFR    )
        CALL PRMEM( 0, 'matel', 'FFQ',    'D', MAXQ+1               )
        CALL PRMEM( 0, 'matel', 'FFY',    'D', MAXFFY               )
        CALL PRMEM( 0, 'matel', 'RI',     'D', MAXR+1               )
        CALL PRMEM( 0, 'matel', 'Y',      'D', 4*MAXLM              )
        CALL PRMEM( 0, 'matel', ' ',      ' ', 0                    )
      ENDIF
C -------------------------------------------------------------------

C Check if tables must be re-initialized ----------------------------
      IF ( IS1.LE.0 .OR. IS2.LE.0 ) THEN
        DO 20 JS = 1,MAXS
          DO 10 JO = -MAXKB,MAXO
            INDF(JS,JO) = 0
   10     CONTINUE
   20   CONTINUE
        DO 40 JF2 = 1,MAXF
          DO 30 JF1 = 1,MAXF
            INDFF(JF1,JF2,1) = 0
            INDFF(JF1,JF2,2) = 0
   30     CONTINUE
   40   CONTINUE
        INDFFY(0) = 0
        NF   = 0
        NFF  = 0
        NFFR = 0
        NFFY = 0
        GOTO 999
      ELSEIF ( IS1.GT.MAXS .OR. IS2.GT.MAXS ) THEN
        STOP 'MATEL: parameter MAXS too small'
      ENDIF
C -------------------------------------------------------------------

C Check argument OPERAT ---------------------------------------------
      IF ( OPERAT .EQ. 'S' ) THEN
        IOPER = 1
      ELSEIF ( OPERAT .EQ. 'T' ) THEN
        IOPER = 2
      ELSE
        STOP 'MATEL: Invalid value of argument OPERAT'
      ENDIF
C -------------------------------------------------------------------

C Find radial expansion of each orbital -----------------------------
      DO 80 I = 1,2
        IF (I .EQ. 1) THEN
          IS = IS1
          IO = IO1
        ELSE
          IS = IS2
          IO = IO2
        ENDIF
        CALL CHKDIM( 'MATEL', 'MAXS',  MAXS,   IS, 1 )
        CALL CHKDIM( 'MATEL', 'MAXO',  MAXO,   IO, 1 )
        CALL CHKDIM( 'MATEL', 'MAXKB', MAXKB, -IO, 1 )
        IF (INDF(IS,IO) .EQ. 0) THEN
*         CALL TIMER( 'MATEL1', 1 )
          NF = NF + 1
          CALL CHKDIM( 'MATEL', 'MAXF', MAXF, NF, 1 )
          PI = 4.D0 * ATAN(1.D0)
          LMAX = MAX( MAXLO, MAXLKB )
C         Factor 2 below is because we will expand the product of
C         two orbitals
          QMAX = 2.D0 * SQRT( Q2CUT )
          RMAX = PI * NQ / QMAX
          IF ( RCUT(IS,IO) .GT. RMAX )
     .      STOP 'MATEL: NQ too small for required cutoff.'

C         Expand orbital in spherical harmonics
          CALL YLMEXP( LMAX, RLYLM, PHIATM, IS, IO, 0, NQ, RMAX,
     .                 NLM, ILM(NF), F(0,NF) )
          IF (NLM.GT.1) THEN
            WRITE(6,*) 'MATEL: PHIATM MUST BE A RADIAL FUNCTION ',
     .                         'TIMES A REAL SPHERICAL HARMONIC'
            STOP
          ENDIF

C         Find orbital in k-space
          L = LOFILM( ILM(NF) )
          CALL RADFFT( L, NQ, RMAX, F(0,NF), F(0,NF) )

          INDF(IS,IO) = NF
*         CALL TIMER( 'MATEL1', 2 )
        ENDIF
   80 CONTINUE
C -------------------------------------------------------------------

C Find radial expansion of overlap ----------------------------------
      IF1 = INDF(IS1,IO1)
      IF2 = INDF(IS2,IO2)
      IF ( INDFF(IF1,IF2,IOPER) .EQ. 0 ) THEN
*       CALL TIMER( 'MATEL2', 1 )

        NFF = NFF + 1
        CALL CHKDIM( 'MATEL', 'MAXFF', MAXFF, NFF, 1 )

C       Find orbitals convolution by multiplication in k-space
        C = ( 2.D0 * PI )**1.5D0
        DO 90 IQ = 0,NQ
          FFQ(IQ) = C * F(IQ,IF1) * F(IQ,IF2)
          IF ( OPERAT .EQ. 'T' ) THEN
            Q = IQ * QMAX / NQ
            FFQ(IQ) = FFQ(IQ) * Q * Q
          ENDIF
   90   CONTINUE

C       Loop on possible values of l quantum number of product
        L1 = LOFILM( ILM(IF1) )
        L2 = LOFILM( ILM(IF2) )
        DO 170 L3 = ABS(L1-L2), L1+L2, 2

C         Return to real space
          CALL RADFFT( L3, NQ, NQ*PI/RMAX, FFQ, FFL )
          IF (MOD(ABS(L1-L2-L3)/2,2) .NE. 0) THEN
            DO 100 IR = 0,NQ
              FFL(IR) = - FFL(IR)
  100       CONTINUE
          ENDIF

C         Divide by R**L
          IF (L3 .NE. 0) THEN
            DO 110 IR = 1,NQ
              R = IR * RMAX / NQ
              FFL(IR) = FFL(IR) / R**L3
  110       CONTINUE
C           Parabolic extrapolation to R=0. (I1=1,I2=2)
            FFL(0) = ( 4.D0 * FFL(I1) - FFL(I2) ) / 3.D0
          ENDIF

C         Select NR out of NQ points
          NQQQ = NQ
          IF (MOD(NQQQ,NR) .NE. 0)
     .      STOP 'MATEL: NQ must be multiple of NR'
          DO 120 IR = 0,NR
            JR = IR * NQ / NR
            FFL(IR) = FFL(JR)
  120     CONTINUE

C         Find if radial function is already in table
          FOUND = .FALSE.
          DO 150 JO1 = -MAXKB,MAXO
          DO 140 JO2 = -MAXKB,MAXO
            JF1 = INDF(IS1,JO1)
            JF2 = INDF(IS2,JO2)
            IF (JF1.NE.0 .AND. JF2.NE.0) THEN
              JFF = INDFF(JF1,JF2,IOPER)
              IF (JFF .NE. 0) THEN
                DO 130 JFFY = INDFFY(JFF-1)+1, INDFFY(JFF)
                  JFFR = INDFFR(JFFY)
                  IF ( PROPOR(NR,FFL(1),FFR(1,1,JFFR),
     .                        FFTOL,CPROP)           ) THEN
                    FOUND = .TRUE.
                    IFFR(L3) = JFFR
                    CFFR(L3) = CPROP
                    GOTO 160
                  ENDIF
  130           CONTINUE
              ENDIF
            ENDIF
  140     CONTINUE
  150     CONTINUE
  160     CONTINUE

C         Store new radial function and setup spline interpolation
          IF (.NOT.FOUND) THEN
            NFFR = NFFR + 1
            CALL CHKDIM( 'MATEL', 'MAXFFR', MAXFFR, NFFR, 1 )
            IFFR(L3) = NFFR
            CFFR(L3) = 1.D0
            DO 165 IR = 0,NR
              RI(IR) = IR * RMAX / NR
              FFR(IR,1,NFFR) = FFL(IR)
  165       CONTINUE
            DFFR0 = 0.D0
            DFFRMX = 0.D0
            CALL SPLIN( RI, FFR(0,1,NFFR), NR+1, DFFR0, DFFRMX,
     .                  FFR(0,2,NFFR) )
          ENDIF

  170   CONTINUE

C       Expand the product of two spherical harmonics (SH) also in SH
        CALL YLMEXP( L1+L2, RLYLM, YLMYLM, ILM(IF1), ILM(IF2),
     .                1, 1, 1.D0, NLM, ILMFF(NFFY+1), FFY(NFFY+1) )

C       Loop on possible values of lm quantum numbers of orbital product
        DO 180 I = 1,NLM
          NFFY = NFFY + 1
          CALL CHKDIM( 'MATEL', 'MAXFFY', MAXFFY, NFFY, 1 )
          JLM = ILMFF(NFFY)
          L3 = LOFILM( JLM )
          INDFFR(NFFY) = IFFR(L3)
          FFY(NFFY) = FFY(NFFY) * CFFR(L3)

  180   CONTINUE
        INDFFY(NFF) = NFFY
        INDFF(IF1,IF2,IOPER) = NFF
*       CALL TIMER( 'MATEL2', 2 )
      ENDIF
C End of initialization section -------------------------------------

C Find value of matrix element and its gradient ---------------------
*     CALL TIMER( 'MATEL3', 1 )

C     Initialize output
      S12 = 0.D0
      DSDR(1) = 0.D0
      DSDR(2) = 0.D0
      DSDR(3) = 0.D0

C     Exit if orbitals are out of range and avoid R12=0
      X12(1) = R12(1)
      X12(2) = R12(2)
      X12(3) = R12(3)
      R = SQRT( X12(1)*X12(1) + X12(2)*X12(2) + X12(3)*X12(3) )
      IF (R .GT. RCUT(IS1,IO1)+RCUT(IS2,IO2) ) THEN
C       Go to exit point
        GOTO 999
      ELSEIF (R .GT. RMAX) THEN
        STOP 'MATEL: NQ too small for required cutoff.'
      ELSEIF (R .LT. TINY) THEN
        X12(3) = TINY
        R = SQRT( X12(1)*X12(1) + X12(2)*X12(2) + X12(3)*X12(3) )
      ENDIF

C     Find spherical harmonics times R**L
      IFF = INDFF(IF1,IF2,IOPER)
      LMAX = 0
      DO 190 IFFY = INDFFY(IFF-1)+1, INDFFY(IFF)
        JLM = ILMFF(IFFY)
        LMAX = MAX( LMAX, LOFILM(JLM) )
  190 CONTINUE
      CALL RLYLM( LMAX, X12, Y, DYDR )

C     Interpolate radial functions and obtain SH expansion
      DO 210 IFFY = INDFFY(IFF-1)+1, INDFFY(IFF)
        JFFR = INDFFR(IFFY)
        CALL SPLINU( 0.D0, RMAX, FFR(0,1,JFFR), FFR(0,2,JFFR),
     .               NR+1, R, SR, DSRDR )
        JLM = ILMFF(IFFY)
        S12 = S12 + SR * FFY(IFFY) * Y(JLM)
        DO 200 IX = 1,3
          DSDR(IX) = DSDR(IX) +
     .               DSRDR * FFY(IFFY) * Y(JLM) * X12(IX) / R +
     .               SR * FFY(IFFY) * DYDR(IX,JLM)
*         write(6,'(a,2i4,3f12.6)')
*    .       'matel: ilm, ix, dsrdr, dydr, ddsdr =',
*    .       ilm, ix, dsrdr, dydr(ix,ilm),
*    .       DSRDR * FFY(IFFY) * Y(JLM) * X12(IX) / R +
*    .       SR * FFY(IFFY) * DYDR(IX,JLM)
  200   CONTINUE
  210 CONTINUE
*     CALL TIMER( 'MATEL3', 2 )
C -------------------------------------------------------------------

C Stop time counter -------------------------------------------------
  999 CONTINUE
*     CALL TIMER( 'MATEL', 2 )
C -------------------------------------------------------------------
      RETURN


      END





