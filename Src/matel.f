      MODULE MATEL_MODULE
c
c     This module contains two 'entries': matel0 and matel.
c
c     The declaration section should come first.
c 
c     This version of matel now uses dynamic memory and
c     consequently matel.h has been removed. JDG Sept 99
c

C
C  Modules
C
      use precision
      use atmfuncs, only: lmxkbfis, lomaxfis, nztfl, nkbl_func,
     $                    rcut, phiatm, xphiatm, yphiatm, zphiatm

      PUBLIC :: MATEL_INIT, MATEL

      PRIVATE

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

C Derived dimension parameters --------------------------------------
C Ref: J.M.Soler notes of 4/4/96.
      INTEGER
     .  MAXF, MAXFF, MAXFFR, MAXFFY, MAXF2,
     .  MAXKB, MAXLO, MAXLKB, MAXLM, MAXLP1,
     .  MAXO, MAXSS, MAXZKB, MAXZO, MAXZZ,
     .  NOPERAT, MAXIND, MAXLP2, MAXF1,
     .  MAXLM2
C -------------------------------------------------------------------

C Internal variable types and dimensions ----------------------------

      INTEGER
     .  I, I1, I2, IF1, IF2, IFF, IFFY,
     .  IO, IOPER, IQ, IR, IS, IX,
     .  JF1, JF2, JFF, JFFR, JFFY, JLM,
     .  JO, JO1, JO2, JR, JS,
     .  L, L1, L2, L3, LMAX, LOFILM,
     .  NF, NFF, NFFR, NFFY, MAXS,
     .  NQQQ, NY, NZMAX, NKBMAX_MATEL, 
     .  IOP, NNLM, IL,
     .  NNLM1, NNLM2, IL1, IL2, NF0, IIF1, IIF2,
     .  IAUXM(3,3), NGRPOL, JJF1, JJF2, JJLM1,
     .  JJLM2, JL1, JL2, J, K

 
      integer, dimension(:), allocatable, save ::
     .  IFFR

      integer, dimension(:), allocatable, save ::
     .  ILM, ILMFF, INDFFY, INDFFR

      integer, dimension(:,:,:), allocatable, save ::
     .  INDFF, INDF, NLM

      DOUBLE PRECISION
     .  C, CPROP, DFFR0, DFFRMX, DSRDR, 
     .  PI, Q, QMAX, R, RMAX,
     .  SR, X12(3), AUXV(3)

      double precision, dimension(:), allocatable, save ::
     .  CFFR, FFQ, FFL, RI, Y

      double precision, dimension(:,:), allocatable, save ::
     .  DYDR

      double precision, dimension(:), allocatable, save ::
     .  FFY

      double precision, dimension(:,:), allocatable, save ::
     .  F

      double precision, dimension(:,:,:), allocatable, save ::
     .  FFR

      LOGICAL
     .  FOUND, FRSTME, PROPOR

      EXTERNAL
     .  CHKDIM, LOFILM, PROPOR, RADFFT, RLYLM,
     .  SPLIN, SPLINU, TIMER, YLMYLM, YLMEXP, memory

      SAVE
     .  FRSTME, I1, I2, NF, NFF, NFFR, NFFY, QMAX, RMAX

C  Save parameters that control dimensions
      SAVE
     .  LMAX, NZMAX, NKBMAX_MATEL, NY, IOPER, MAXS, 
     .  MAXLO, MAXLKB, MAXZO, MAXZKB, MAXLP1, MAXLP2,
     .  MAXLM, MAXLM2, MAXO, MAXKB, MAXSS, MAXZZ,
     .  MAXF1, MAXF, MAXIND, MAXFF, MAXFFR, MAXFFY,
     .  NOPERAT, MAXF2

      DATA
     .  I1, I2    /1, 2/,
     .  NF        /0/,
     .  NFF       /0/,
     .  NFFR      /0/,
     .  NFFY      /0/,
     .  FRSTME    /.TRUE./
C -------------------------------------------------------------------

      CONTAINS

      SUBROUTINE MATEL_INIT( NS )
C *******************************************************************
C Calculates dimension parameters for matel
C Modified by DSP, Aug. 1998, and June 1999
C *********** INPUT *************************************************
C INTEGER NS              : Number of species
C *******************************************************************

      IMPLICIT NONE
      INTEGER NS

      LMAX = 0
      NZMAX = 0
      NKBMAX_MATEL = 0
      MAXS = NS
      DO 320 IS = 1,NS
        LMAX = MAX( LMAX, LOMAXFIS(IS), LMXKBFIS(IS) )
        DO 310 L = 0,LOMAXFIS(IS)
          NZMAX = MAX( NZMAX, NZTFL(IS,L) )
  310   CONTINUE
        DO 315 L = 0,LMXKBFIS(IS)
          NKBMAX_MATEL = MAX( NKBMAX_MATEL, NKBL_FUNC(IS,L))
  315 CONTINUE
  320 CONTINUE
      NY = (LMAX+1)**2

      CALL REPOL(IAUXM,AUXV)

      NGRPOL=0
      DO IS=1,3
        NGRPOL=NGRPOL+ IAUXM(1,IS)*IAUXM(2,IS)*IAUXM(3,IS)
      ENDDO
      IF (NGRPOL.EQ.0) THEN
        IOPER=2
      ELSE
        IOPER=5
      ENDIF

C  Calculate derived array dimensions needed for matel
      MAXLO  = LMAX
      MAXLKB = LMAX
      MAXZO  = NZMAX
      MAXZKB = NKBMAX_MATEL
      NOPERAT = IOPER
      MAXLP1 = LMAX + 1
      MAXLP2 = MAXLP1 + 1
      MAXLM  = (LMAX+1)*(LMAX+1)
      MAXLM2 = (LMAX+2)*(LMAX+2)
      MAXO   = MAXZO  * (MAXLO+1)  * (MAXLO+1)
      MAXKB  = MAXZKB * (MAXLKB+1) * (MAXLKB+1)
      MAXSS  = MAXS * (MAXS+1) / 2
      MAXZZ  = MAXZO * (2*MAXZO+MAXZKB)  + 2
      MAXF1  = MAXS * (MAXO+MAXKB+1)
      MAXF   = MAXS * (MAXO+MAXKB+1) *
     .         (1+4*(IOPER-2))
      MAXIND = MAXS * (MAXO+MAXKB+1) * IOPER
      MAXFF  = MAXSS * MAXZZ * MAXLM * MAXLM
     .             + 4*(IOPER-2)*MAXSS*MAXZZ*MAXLM2*MAXLM2
      MAXFFR = MAXSS * MAXZZ *
     .                     (LMAX+1) * (LMAX+2) * (2*LMAX+3) / 6
     .             + 4*(IOPER-2)* MAXSS * MAXZZ *
     .               (MAXLP1+1) * (MAXLP1+2) * (2*MAXLP1+3) / 6
      MAXFFY =  MAXSS * MAXZZ * MAXLP1 *
     .            (MAXLP1*MAXLP1*MAXLP1*(16*MAXLP1+15)-1) / 30
     .             + 4*(IOPER-2)*MAXSS * MAXZZ * MAXLP2 *
     .            (MAXLP2*MAXLP2*MAXLP2*(16*MAXLP2+15)-1) / 30
      MAXF2  = MAXF1 * MAXF * IOPER

C  Allocate local arrays that must be preserved between calls
      if (.not.allocated(ILM)) then
        allocate(ILM(MAXF+NY))
        call memory('A','I',MAXF+NY,'matel')
      endif
      if (.not.allocated(ILMFF)) then
        allocate(ILMFF(MAXFFY))
        call memory('A','I',MAXFFY,'matel')
      endif
      if (.not.allocated(INDFFY)) then
        allocate(INDFFY(0:MAXFF))
        call memory('A','I',MAXFF+1,'matel')
      endif
      if (.not.allocated(INDFFR)) then
        allocate(INDFFR(MAXFFY))
        call memory('A','I',MAXFFY,'matel')
      endif
      if (.not.allocated(INDFF)) then
        allocate(INDFF(MAXF1,MAXF,NOPERAT))
        call memory('A','I',MAXF1*MAXF*NOPERAT,'matel')
      endif
      if (.not.allocated(INDF)) then
        allocate(INDF(MAXS,-MAXKB:MAXO,NOPERAT))
        call memory('A','I',MAXS*(MAXKB+MAXO+1)*NOPERAT,'matel')
      endif
      if (.not.allocated(NLM)) then
        allocate(NLM(MAXS,-MAXKB:MAXO,NOPERAT))
        call memory('A','I',MAXS*(MAXKB+MAXO+1)*NOPERAT,'matel')
      endif
      if (.not.allocated(FFY)) then
        allocate(FFY(MAXFFY))
        call memory('A','D',MAXFFY,'matel')
      endif
      if (.not.allocated(FFR)) then
        allocate(FFR(0:NR,2,MAXFFR))
        call memory('A','D',(NR+1)*2*MAXFFR,'matel')
      endif
      if (.not.allocated(F)) then
        allocate(F(0:NQ,MAXF+NY))
        call memory('A','D',(NQ+1)*(MAXF+NY),'matel')
      endif

C  Initialise local arrays that have just been allocated
      DO I = 1,NOPERAT
        DO J = -MAXKB,MAXO
          DO K = 1,MAXS
            INDF(K,J,I) = 0
            NLM(K,J,I) = 0
          ENDDO
        ENDDO
      ENDDO
      DO I = 1,NOPERAT
        DO J = 1,MAXF
          DO K = 1,MAXF1
            INDFF(K,J,I) = 0
          ENDDO
        ENDDO
      ENDDO
      INDFFY(0) = 0

      END SUBROUTINE MATEL_INIT
c
c
c

      SUBROUTINE MATEL( OPERAT, IS1, IS2, IO1, IO2, R12, S12, DSDR )
C *******************************************************************
C Finds the overlap or laplacian matrix elements between two atomic
C orbitals. Written by J.M.Soler. April 1995.
C Modified to calculate the matrix elements of the position operator
C by DSP. June, 1999
C ************************* INPUT ***********************************
C CHARACTER OPERAT : Operator to be used: 'S' => Unity (overlap)
C                                         'T' => -Laplacian
C                                         'T' => -Laplacian
C                                         'X' =>  x
C                                         'Y' =>  y
C                                         'Z' =>  z
C What is actually computed is the matrix elements:
C           'X' =>       < phi1(r-R12)|x|phi2(r)>
C           'Y' =>       < phi1(r-R12)|y|phi2(r)>
C           'Z' =>       < phi1(r-R12)|z|phi2(r)>
C i.e the origin is taken at the position of the second atom.
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
C SUBROUTINE XPHIATM(IS,IO,R,PHI,GRPHI)
C   The same as PHIATM but multiply by x.
C
C SUBROUTINE YPHIATM(IS,IO,R,PHI,GRPHI)
C   The same as PHIATM but multiply by y.
C
C SUBROUTINE ZPHIATM(IS,IO,R,PHI,GRPHI)
C   The same as PHIATM but multiply by z.
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


C Argument types and dimensions -------------------------------------
      IMPLICIT          NONE
      CHARACTER         OPERAT
      INTEGER           IO1, IO2, IS1, IS2
      DOUBLE PRECISION  DSDR(3), R12(3), S12
C -------------------------------------------------------------------

C Start time counter ------------------------------------------------
*     CALL TIMER( 'MATEL', 1 )
C -------------------------------------------------------------------

C Allocate local memory
      allocate(IFFR(0:2*MAXLP1))
      call memory('A','I',2*MAXLP1+1,'matel')
      allocate(CFFR(0:2*MAXLP1))
      call memory('A','D',2*MAXLP1+1,'matel')
      allocate(DYDR(3,4*MAXLM2))
      call memory('A','D',3*(4*MAXLM2),'matel')
      allocate(FFQ(0:NQ))
      call memory('A','D',NQ+1,'matel')
      allocate(FFL(0:NQ))
      call memory('A','D',NQ+1,'matel')
      allocate(RI(0:NR))
      call memory('A','D',NR+1,'matel')
      allocate(Y(4*MAXLM2))
      call memory('A','D',4*MAXLM2,'matel')

      PI = 4.D0 * ATAN(1.D0)

C Print size of arrays ----------------------------------------------
      IF (FRSTME) THEN
        FRSTME = .FALSE.
      ENDIF
C -------------------------------------------------------------------

C Check if tables must be re-initialized ----------------------------
      IF ( IS1.LE.0 .OR. IS2.LE.0 ) THEN
       DO 45  IOP = 1, NOPERAT
          DO 20 JS = 1,MAXS
            DO 10 JO = -MAXKB,MAXO
               INDF(JS,JO,IOP) = 0
               NLM(JS,JO,IOP) = 0
   10       CONTINUE
   20     CONTINUE
          DO 40 JF2 = 1,MAXF
           DO 30 JF1 = 1,MAXF
            INDFF(JF1,JF2,IOP) = 0
            INDFF(JF1,JF2,IOP) = 0
   30     CONTINUE
   40    CONTINUE
   45   CONTINUE
        INDFFY(0) = 0
        NF   = 0
        NFF  = 0
        NFFR = 0
        NFFY = 0
        GOTO 999
      ELSEIF ( IS1.GT.MAXS .OR. IS2.GT.MAXS ) THEN
         call die('MATEL: parameter MAXS too small')
      ENDIF
C -------------------------------------------------------------------

C Check argument OPERAT ---------------------------------------------
      IF ( OPERAT .EQ. 'S' ) THEN
        IOPER = 1
      ELSEIF ( OPERAT .EQ. 'T' ) THEN
        IOPER = 2
      ELSEIF ( OPERAT .EQ. 'X' ) THEN
        IOPER = 3
      ELSEIF ( OPERAT .EQ. 'Y' ) THEN
        IOPER = 4
      ELSEIF ( OPERAT .EQ. 'Z' ) THEN
        IOPER = 5
      ELSE
        call die('MATEL: Invalid value of argument OPERAT')
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
        IF (((INDF(IS,IO,IOPER) .EQ. 0).AND.(I.EQ.2)).OR.
     .      ((INDF(IS,IO,1).EQ. 0).AND.(I.EQ.1))) THEN
*         CALL TIMER( 'MATEL1', 1 )
          NF = NF + 1
          CALL CHKDIM( 'MATEL', 'MAXF', MAXF, NF, 1 )
          PI = 4.D0 * ATAN(1.D0)
          LMAX = MAX( MAXLO, MAXLKB )
C         Factor 2 below is because we will expand the product of
C         two orbitals
          QMAX = 2.D0 * SQRT( Q2CUT )
          RMAX = PI * NQ / QMAX
          IF ( RCUT(IS,IO) .GT. RMAX ) THEN
            call die('MATEL: NQ too small for required cutoff.')
          ENDIF

C         Expand orbital in spherical harmonics
          IF ((I.EQ.1).OR.
     .        ((I.EQ.2).AND.(IOPER.LE.2))) THEN
            CALL YLMEXP( LMAX, RLYLM, PHIATM, IS, IO, 0, NQ, RMAX,
     .                 NNLM, ILM(NF), F(0,NF) )
            IF (NNLM.GT.1) THEN
               call die('MATEL: PHIATM MUST BE A RADIAL FUNCTION '
     $                 // 'TIMES A REAL SPHERICAL HARMONIC')
            ENDIF
             INDF(IS,IO,1) = NF
             INDF(IS,IO,2) = NF
             NLM(IS,IO,1) =  1
             NLM(IS,IO,2) =  1
          ELSE
            IF(IOPER.EQ.3) THEN
               CALL YLMEXP(LMAX+1, RLYLM, XPHIATM, IS, IO, 0, NQ, RMAX,
     .                NNLM, ILM(NF), F(0,NF) )
            ELSEIF(IOPER.EQ.4) THEN
               CALL YLMEXP(LMAX+1, RLYLM, YPHIATM, IS, IO, 0, NQ, RMAX,
     .                NNLM, ILM(NF), F(0,NF) )
            ELSE
               CALL YLMEXP(LMAX+1, RLYLM, ZPHIATM, IS, IO, 0, NQ, RMAX,
     .                NNLM, ILM(NF), F(0,NF) )
            ENDIF
            INDF(IS,IO,IOPER) = NF
            NLM(IS,IO,IOPER) = NNLM
          ENDIF

          NF0=NF
          DO IL=1,NNLM
             NF=NF0+IL-1
C         Find orbital in k-space
             L = LOFILM( ILM(NF) )
             CALL RADFFT( L, NQ, RMAX, F(0,NF), F(0,NF) )
          ENDDO
*         CALL TIMER( 'MATEL1', 2 )
        ENDIF
   80 CONTINUE
C -------------------------------------------------------------------

C Find radial expansion of overlap ----------------------------------
      IF1 = INDF(IS1,IO1,1)
      IF2 = INDF(IS2,IO2,IOPER)
      NNLM1=NLM(IS1,IO1,1)
      NNLM2=NLM(IS2,IO2,IOPER)

      DO 186  IL1=1, NNLM1
         DO 185 IL2= 1, NNLM2
            IIF1=IF1+(IL1-1)
            IIF2=IF2+(IL2-1)

       IF ( INDFF(IIF1,IIF2,IOPER) .EQ. 0 ) THEN
*       CALL TIMER( 'MATEL2', 1 )

        NFF = NFF + 1
        CALL CHKDIM( 'MATEL', 'MAXFF', MAXFF, NFF, 1 )

C       Find orbitals convolution by multiplication in k-space
        C = ( 2.D0 * PI )**1.5D0
        DO 90 IQ = 0,NQ
          FFQ(IQ) = C * F(IQ,IIF1) * F(IQ,IIF2)
          IF ( OPERAT .EQ. 'T' ) THEN
            Q = IQ * QMAX / NQ
            FFQ(IQ) = FFQ(IQ) * Q * Q
          ENDIF
   90   CONTINUE

C       Loop on possible values of l quantum number of product
        L1 = LOFILM( ILM(IIF1) )
        L2 = LOFILM( ILM(IIF2) )
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
     $         call die('MATEL: NQ must be multiple of NR')
          DO 120 IR = 0,NR
            JR = IR * NQ / NR
            FFL(IR) = FFL(JR)
  120     CONTINUE

C         Find if radial function is already in table
          FOUND = .FALSE.
          DO 150 JO1 = -MAXKB,MAXO
          DO 140 JO2 = -MAXKB,MAXO
            JF1 = INDF(IS1,JO1,1)
            JF2 = INDF(IS2,JO2,IOPER)
            IF (JF1.NE.0 .AND. JF2.NE.0) THEN
             JJLM1=NLM(IS1,JO1,1)
             JJLM2=NLM(IS2,JO2,IOPER)
             DO 143  JL1=1, JJLM1
               DO 142 JL2=1, JJLM2
                  JJF1=JF1+(JL1-1)
                  JJF2=JF2+(JL2-1)
                  JFF = INDFF(JJF1,JJF2,IOPER)
                  IF (JFF .NE. 0) THEN
                    DO 130 JFFY = INDFFY(JFF-1)+1, INDFFY(JFF)
                       JFFR = INDFFR(JFFY)
                       IF ( PROPOR(NR,FFL(1),FFR(1,1,JFFR),
     .                       FFTOL,CPROP)           ) THEN
                          FOUND = .TRUE.
                          IFFR(L3) = JFFR
                          CFFR(L3) = CPROP
                          GOTO 160
                       ENDIF
  130               CONTINUE
                   ENDIF
  142          CONTINUE
  143        CONTINUE
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
            DFFR0 =  huge(1.d0)
            DFFRMX = 0.D0
            CALL SPLIN( RI, FFR(0,1,NFFR), NR+1, DFFR0, DFFRMX,
     .                  FFR(0,2,NFFR) )
          ENDIF
  170   CONTINUE


C       Expand the product of two spherical harmonics (SH) also in SH
        CALL YLMEXP( L1+L2, RLYLM, YLMYLM, ILM(IIF1), ILM(IIF2),
     .                1, 1, 1.D0, NNLM, ILMFF(NFFY+1), FFY(NFFY+1) )

C    Loop on possible values of lm quantum numbers of orbital product
        DO 180 I = 1,NNLM
          NFFY = NFFY + 1
          CALL CHKDIM( 'MATEL', 'MAXFFY', MAXFFY, NFFY, 1 )
          JLM = ILMFF(NFFY)
          L3 = LOFILM( JLM )
          INDFFR(NFFY) = IFFR(L3)
          FFY(NFFY) = FFY(NFFY) * CFFR(L3)
  180   CONTINUE

        INDFFY(NFF) = NFFY
        INDFF(IIF1,IIF2,IOPER) = NFF
*       CALL TIMER( 'MATEL2', 2 )
       
        ENDIF
  185   CONTINUE
  186   CONTINUE

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
          call die('MATEL: NQ too small for required cutoff.')
      ELSEIF (R .LT. TINY) THEN
        X12(3) = TINY
        R = SQRT( X12(1)*X12(1) + X12(2)*X12(2) + X12(3)*X12(3) )
      ENDIF

C     Find spherical harmonics times R**L
      IF1 = INDF(IS1,IO1,1)
      IF2 = INDF(IS2,IO2,IOPER)
      NNLM1=NLM(IS1,IO1,1)
      NNLM2=NLM(IS2,IO2,IOPER)
      DO 250  IL1=1, NNLM1
         DO 240 IL2= 1, NNLM2
            IIF1=IF1+(IL1-1)
            IIF2=IF2+(IL2-1)
            IFF = INDFF(IIF1,IIF2,IOPER)
            LMAX = 0
            DO 190 IFFY = INDFFY(IFF-1)+1, INDFFY(IFF)
               JLM = ILMFF(IFFY)
               LMAX = MAX( LMAX, LOFILM(JLM) )
  190       CONTINUE
            CALL RLYLM( LMAX, X12, Y, DYDR )

C     Interpolate radial functions and obtain SH expansion
            DO 210 IFFY = INDFFY(IFF-1)+1, INDFFY(IFF)
                JFFR = INDFFR(IFFY)
                 CALL SPLINU( 0.D0, RMAX, FFR(0,1,JFFR),
     .                 FFR(0,2,JFFR), NR+1, R, SR, DSRDR )
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
  200            CONTINUE
  210      CONTINUE
  240    CONTINUE
  250   CONTINUE
*     CALL TIMER( 'MATEL3', 2 )
C -------------------------------------------------------------------

C Stop time counter -------------------------------------------------
  999 CONTINUE

C Deallocate local memory
      call memory('D','I',size(IFFR),'matel')
      deallocate(IFFR)
      call memory('D','D',size(CFFR),'matel')
      deallocate(CFFR)
      call memory('D','D',size(DYDR),'matel')
      deallocate(DYDR)
      call memory('D','D',size(FFQ),'matel')
      deallocate(FFQ)
      call memory('D','D',size(FFL),'matel')
      deallocate(FFL)
      call memory('D','D',size(RI),'matel')
      deallocate(RI)
      call memory('D','D',size(Y),'matel')
      deallocate(Y)

*     CALL TIMER( 'MATEL', 2 )
C -------------------------------------------------------------------
      RETURN

      END SUBROUTINE MATEL

      END MODULE MATEL_MODULE

