C $Id: ylmexp.f,v 1.5 1999/01/31 11:45:18 emilio Exp $

      SUBROUTINE YLMEXP( LMAX, RLYLM, FUNC, IS, IO, IR1, NR, RMAX,
     .                   NY, ILM, FLM )
C *******************************************************************
C Makes a radial times spherical-harmonic expansion of a function.
C Written by J.M.Soler. September 1995.
C ************************* INPUT ***********************************
C INTEGER  LMAX                     : Max. ang. momentum quantum number
C EXTERNAL RLYLM(Lmax,Rvec,Y,GradY) : Returns spherical harmonics,
C                                     times r**l
C EXTERNAL FUNC(IS,IO,Rvec,F,GradF) : Function to be expanded.
C INTEGER  IS, IO                   : Indexes to call FUNC.
C INTEGER  IR1                      : First radial index.
C INTEGER  NR                       : Last radial index.
C INTEGER  RMAX                     : Maximum radius: R(IR)=RMAX*IR/NR
C ************************* OUTPUT **********************************
C INTEGER  NY            : Number of spherical harmonics required.
C INTEGER  ILM(NY)       : Spherical harmonics indexes: ILM=L*L+L+M+1.
C REAL*8   FLM(IR1:NR,NY): Radial expansion for each spherical harmonic:
C   F(Rvec) = Sum_iy( FLM(Rmod,iy) * Y(Rdir,ILM(iy)) ),
C   with   Rmod = RMAX*IR/NR
C ************************* UNITS ***********************************
C RMAX must be in the same units of the argument Rvec in FUNC.
C FLM is in the same units of the argument F in FUNC.
C ************************* BEHAVIOUR *******************************
C If function FUNC contains angular components with L higher than LMAX,
C   they will corrupt those computed with lower L.
C *******************************************************************

C Next line is non-standard and may be suppressed -------------------
      IMPLICIT NONE
C -------------------------------------------------------------------

C Declare argument types and dimensions -----------------------------
      INTEGER           ILM(*), IR1, IS, IO, LMAX, NY, NR
      DOUBLE PRECISION  FLM(IR1:NR,*), RMAX
      EXTERNAL          FUNC, RLYLM
C -------------------------------------------------------------------

C Tolerance for FLM -------------------------------------------------
      DOUBLE PRECISION FTOL
      PARAMETER ( FTOL = 1.D-12 )
C -------------------------------------------------------------------

C Dimension parameters for internal variables -----------------------
      INTEGER MAXL, MAXLM, MAXSP
      PARAMETER ( MAXL  = 8 )
      PARAMETER ( MAXLM = (MAXL+1) * (MAXL+1) )
      PARAMETER ( MAXSP = (MAXL+1) * (2*MAXL+1) )
C -------------------------------------------------------------------

C Declare internal variables ----------------------------------------
      INTEGER
     .  IM, IR, ISP, IZ, JLM, JR, JY, NLM, NSP
      DOUBLE PRECISION
     .  DFDX(3), DOT, DYDX(3,MAXLM),
     .  F(MAXSP), FY, PHI, PI, R, RX(3), THETA, 
     .  W(MAXL+1), WSP,
     .  X(3,MAXSP), Y(MAXLM), YSP(MAXSP,MAXLM),
     .  Z(MAXL+1)
      EXTERNAL
     .  CHKDIM, DOT, GAULEG, TIMER
C -------------------------------------------------------------------

C Start time counter ------------------------------------------------
*     CALL TIMER( 'YLMEXP', 1 )
C -------------------------------------------------------------------

C Check maximum angular momentum ------------------------------------
      CALL CHKDIM( 'YLMEXP', 'MAXL', MAXL, LMAX, 1 )
C -------------------------------------------------------------------

C Find special points and weights for gaussian quadrature -----------
      CALL GAULEG( -1.D0, 1.D0, Z, W, LMAX+1 )
C -------------------------------------------------------------------

C Find weighted spherical harmonics at special points ---------------
      PI = 4.D0 * ATAN(1.D0)
      NLM = (LMAX+1)**2
      NSP = 0
      DO 30 IZ = 1,LMAX+1
        THETA = ACOS( Z(IZ) )
        DO 20 IM = 0,2*LMAX
          NSP = NSP + 1
          PHI = IM * 2.D0 * PI / (2*LMAX+1)
          X(1,NSP) = SIN(THETA) * COS(PHI)
          X(2,NSP) = SIN(THETA) * SIN(PHI)
          X(3,NSP) = COS(THETA)
          CALL RLYLM( LMAX, X(1,NSP), Y, DYDX )
          WSP = W(IZ) * (2.D0*PI) / (2*LMAX+1)
          DO 10 JLM = 1,NLM
            YSP(NSP,JLM) = WSP * Y(JLM)
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
C -------------------------------------------------------------------

C Expand FUNC in spherical harmonics at each radius -----------------
      NY = 0

C     Loop on radial points
      DO 80 IR = IR1,NR
        R = IR * RMAX / NR

C       Find function at points on a sphere of radius R
        DO 40 ISP = 1,NSP
          RX(1) = R * X(1,ISP)
          RX(2) = R * X(2,ISP)
          RX(3) = R * X(3,ISP)
          CALL FUNC( IS, IO, RX, F(ISP), DFDX )
   40   CONTINUE

C       Expand F(R) in spherical harmonics
        DO 70 JLM = 1,NLM
          FY = DOT( F, YSP(1,JLM), NSP )
          IF ( ABS(FY) .GT. FTOL ) THEN

C           Find JY corresponding to JLM
            DO 50 JY = 1,NY
              IF ( ILM(JY) .EQ. JLM ) THEN
                FLM(IR,JY) = FY
                GOTO 70
              ENDIF
   50       CONTINUE

C           New spherical harmonic required.
            NY = NY + 1
            ILM(NY) = JLM
            DO 60 JR = IR1,NR
              FLM(JR,NY) = 0.D0
   60       CONTINUE
            FLM(IR,NY) = FY
          ENDIF
   70   CONTINUE

   80 CONTINUE
C -------------------------------------------------------------------

C Special case for zero function ------------------------------------
      IF (NY .EQ. 0) THEN
        NY = 1
        ILM(1) = 1
        DO 90 IR = IR1,NR
          FLM(IR,1) = 0.D0
   90   CONTINUE
      ENDIF
C -------------------------------------------------------------------

C Stop time counter -------------------------------------------------
*     CALL TIMER( 'YLMEXP', 2 )
C -------------------------------------------------------------------

      END





