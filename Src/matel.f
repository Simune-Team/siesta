      SUBROUTINE MATEL( OPERAT, IS1, IS2, IO1, IO2, R12, S12, DSDR )
C *******************************************************************
C Finds the overlap or laplacian matrix elements between two atomic
C orbitals. Written by J.M.Soler. April 1995.
C Matrix elements of the position operator by DSP. June, 1999
C Last modifications by JMS. January 2001.
C ************************* INPUT ***********************************
C CHARACTER OPERAT : Operator to be used: 'S' => Unity (overlap)
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
C INTEGER FUNCTION LOFIO(IS,IO)
C   Returns the total angular momentum of orbitals and KB proyectors.
C Input:
C     INTEGER IS : Species index
C     INTEGER IO : Orbital index
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
C ************************* UNITS ***********************************
C Length units are arbitrary, but must be consistent in MATEL, RCUT
C   and PHIATM. The laplacian unit is (length unit)**(-2).
C ************************* BEHAVIOUR *******************************
C 1) Returns exactly zero if |R12| > RCUT(IS1,IO1) + RCUT(IS2,IO2)
C 2) If (IS1.LE.0 .OR. IS2.LE.0) all internal tables are erased for
C    reinitialization.
C *******************************************************************
C    6  10        20        30        40        50        60        7072

C Modules -----------------------------------------------------------
      use precision, only : dp
      USE ALLOC
      USE ATMFUNCS, ONLY: 
     .  LOFIO, PHIATM, RCUT, XPHIATM, YPHIATM, ZPHIATM
      use spher_harm
      use m_radfft
C -------------------------------------------------------------------

C Argument types and dimensions -------------------------------------
      IMPLICIT          NONE
      CHARACTER         OPERAT
      INTEGER           IO1, IO2, IS1, IS2
      real(dp)          DSDR(3), R12(3), S12
C -------------------------------------------------------------------

C Internal precision parameters  ------------------------------------
C  NR is the number of radial points for matrix-element tables.
C  NQ is the number of radial points in reciprocal space.
C  EXPAND is a factor to expand some array sizes
C  Q2CUT is the required planewave cutoff for orbital expansion
C    (in Ry if lengths are in Bohr).
C  FFTOL is the tolerance for considering equal the radial part of
C    two orbitals.
C  TINY is a small number to add to a zero denominator
      INTEGER,          PARAMETER :: NR     =  128
      INTEGER,          PARAMETER :: NQ     =  1024
      real(dp),         PARAMETER :: EXPAND =  1.20D0
      real(dp),         PARAMETER :: Q2CUT  =  2.50D3
      real(dp),         PARAMETER :: FFTOL  =  1.D-8
      real(dp),         PARAMETER :: TINY   =  1.D-12
      CHARACTER(LEN=*), PARAMETER :: MYNAME =  'MATEL '
C -------------------------------------------------------------------

C Internal variable types and dimensions ----------------------------
      INTEGER ::
     .  I, IF1, IF2, IFF, IFFY, IFLM1, IFLM2, 
     .  IO, IOPER, IQ, IR, IS, IX,
     .  JF1, JF2, JFF, JFFR, JFFY, JFLM1, JFLM2, JLM, 
     .  JO1, JO2, JR,
     .  L, L1, L2, L3, LMAX,
     .  N, NILM, NILM1, NILM2, NJLM1, NJLM2
      INTEGER, SAVE ::
     .  MF=0, MFF=0, MFFR=0, MFFY=0, 
     .  NF=0, NFF=0, NFFR=0, NFFY=0

      INTEGER, POINTER, SAVE ::
     .  IFFR(:), ILM(:), ILMFF(:), INDF(:,:,:), INDFF(:,:,:),
     .  INDFFR(:), INDFFY(:), NLM(:,:,:)

      real(dp) ::
     .  C, CPROP, DFFR0, DFFRMX, DSRDR, FFL(0:NQ), FFQ(0:NQ),
     .  Q, R, SR, X12(3)

      real(dp), SAVE ::
     .  PI, QMAX, RMAX

      real(dp), POINTER, SAVE ::
     .  CFFR(:), DYDR(:,:), F(:,:), FFR(:,:,:), FFY(:), Y(:)

      LOGICAL ::
     .  FOUND, PROPOR

      LOGICAL, SAVE ::
     .  NULLIFIED=.FALSE.

!      logical, save :: opened=.false.   ! JMS debug

      TYPE(allocDefaults) ::
     .  OLDEFS

      EXTERNAL
     .  PROPOR, SPLINE, SPLINT, TIMER

C -------------------------------------------------------------------

C Start time counter 
*     CALL TIMER( MYNAME, 1 )

C Nullify pointers 
      IF (.NOT.NULLIFIED) THEN
        NULLIFY( IFFR, ILM, ILMFF, INDF, INDFF, INDFFR, INDFFY, NLM,
     .           CFFR, DYDR, F, FFR, FFY, Y )
        ALLOCATE( INDF(0,0,0) )
        CALL RE_ALLOC( INDFFY, 0,MFF, MYNAME//'INDFFY' )
        INDFFY(0) = 0
        NULLIFIED = .TRUE.
      ENDIF

C Set allocation defaults 
      CALL ALLOC_DEFAULT( OLD=OLDEFS, ROUTINE=MYNAME, 
     .                    COPY=.TRUE., SHRINK=.FALSE. )

C Check if tables must be re-initialized 
      IF ( IS1.LE.0 .OR. IS2.LE.0 ) THEN
        CALL DE_ALLOC( IFFR,   MYNAME//'IFFR'   )
        CALL DE_ALLOC( ILM,    MYNAME//'ILM'    )
        CALL DE_ALLOC( ILMFF,  MYNAME//'ILMFF'  )
        CALL DE_ALLOC( INDF,   MYNAME//'INDF'   )
        CALL DE_ALLOC( INDFF,  MYNAME//'INDFF'  )
        CALL DE_ALLOC( INDFFR, MYNAME//'INDFFR' )
        CALL DE_ALLOC( INDFFY, MYNAME//'INDFFY' )
        CALL DE_ALLOC( NLM,    MYNAME//'NLM'    )
        CALL DE_ALLOC( CFFR,   MYNAME//'CFFR'   )
        CALL DE_ALLOC( DYDR,   MYNAME//'DYDR'   )
        CALL DE_ALLOC( F,      MYNAME//'F'      )
        CALL DE_ALLOC( FFR,    MYNAME//'FFR'    )
        CALL DE_ALLOC( FFY,    MYNAME//'FFY'    )
        CALL DE_ALLOC( Y,      MYNAME//'Y'      )
        MF   = 0
        MFF  = 0
        MFFR = 0
        MFFY = 0
        NF   = 0
        NFF  = 0
        NFFR = 0
        NFFY = 0
        CALL RE_ALLOC( INDFFY, 0,MFF, MYNAME//'INDFFY' )
        INDFFY(0) = 0
      ENDIF

C Check argument OPERAT 
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

C Check size of orbital index table
      IF ( MAX(IS1,IS2).GT.SIZE(INDF,1)   .OR.
     .     MIN(IO1,IO2).LT.LBOUND(INDF,2) .OR.
     .     MAX(IO1,IO2).GT.UBOUND(INDF,2) .OR.
     .     MAX(IOPER,2).GT.SIZE(INDF,3) ) THEN
        CALL RE_ALLOC( INDF, 1,MAX(SIZE(INDF,1),IS1,IS2),
     .                MIN(LBOUND(INDF,2),IO1,IO2),MAX(UBOUND(INDF,2),
     .                IO1,IO2),1,MAX(SIZE(INDF,3),IOPER,2),
     .                MYNAME//'INDF' )
        CALL RE_ALLOC( NLM,  1,MAX(SIZE(INDF,1),IS1,IS2),
     .                MIN(LBOUND(INDF,2),IO1,IO2),MAX(UBOUND(INDF,2),
     .                IO1,IO2),1,MAX(SIZE(INDF,3),IOPER,2),
     .                MYNAME//'NLM'  )
      ENDIF

C JMS debug
c      if (.not.opened) then
c        open(41,file='matel.out',form='formatted',status='unknown')
c        opened = .true.
c      endif

C Find radial expansion of each orbital -----------------------------
      DO I = 1,2
        IF (I .EQ. 1) THEN
          IS = IS1
          IO = IO1
          FOUND = (INDF(IS,IO,1) .NE. 0)
        ELSE
          IS = IS2
          IO = IO2
          FOUND = (INDF(IS,IO,IOPER) .NE. 0)
        ENDIF
        IF (.NOT.FOUND) THEN
*         CALL TIMER( 'MATEL1', 1 )
          PI = 4.D0 * ATAN(1.D0)
C         Factor 2 below is because we will expand the product of
C         two orbitals
          QMAX = 2.D0 * SQRT( Q2CUT )
          RMAX = PI * NQ / QMAX
          IF ( RCUT(IS,IO) .GT. RMAX ) THEN
            call die('MATEL: NQ too small for required cutoff.')
          ENDIF

C         Reallocate arrays
          L = LOFIO( IS, IO )
          NILM = (L+2)**2
          IF (NF+NILM .GT. MF) MF = EXPAND * (NF+NILM)
          CALL RE_ALLOC( F, 0,NQ, 1,MF, MYNAME//'F'   )
          CALL RE_ALLOC( ILM,     1,MF, MYNAME//'ILM' )
          CALL RE_ALLOC( INDFF,   1,MF, 1,MF, 1,MAX(IOPER,2),
     .                  MYNAME//'INDFF' )

C         Expand orbital in spherical harmonics
          IF ((I.EQ.1) .OR. (IOPER.LE.2)) THEN
            CALL YLMEXP( L, RLYLM, PHIATM, IS, IO, 0, NQ, RMAX,
     .                   NILM, ILM(NF+1:), F(0:,NF+1:) )
            INDF(IS,IO,1) = NF+1
            INDF(IS,IO,2) = NF+1
            NLM(IS,IO,1) = NILM
            NLM(IS,IO,2) = NILM
          ELSE
            IF(IOPER.EQ.3) THEN
              CALL YLMEXP( L+1, RLYLM, XPHIATM, IS, IO, 0, NQ, RMAX,
     .                     NILM, ILM(NF+1:), F(0:,NF+1:) )
            ELSEIF(IOPER.EQ.4) THEN
              CALL YLMEXP( L+1, RLYLM, YPHIATM, IS, IO, 0, NQ, RMAX,
     .                     NILM, ILM(NF+1:), F(0:,NF+1:) )
            ELSE
              CALL YLMEXP( L+1, RLYLM, ZPHIATM, IS, IO, 0, NQ, RMAX,
     .                     NILM, ILM(NF+1:), F(0:,NF+1:) )
            ENDIF
            INDF(IS,IO,IOPER) = NF+1
            NLM(IS,IO,IOPER) = NILM
          ENDIF

C         Store orbital in k-space
          DO JLM = 1,NILM
            NF = NF + 1
            L = LOFILM( ILM(NF) )
            CALL RADFFT( L, NQ, RMAX, F(0:NQ,NF), F(0:NQ,NF) )
*           F(NQ,NF) = 0.D0
          ENDDO

C JMS debug
c          write(41,'(/,a,4i6)') 'is,io,l,nf =', is, io, l, nf
c          write(41,'(i6,e25.12)') (iq,f(iq,nf),iq=0,nq)

*         CALL TIMER( 'MATEL1', 2 )
        ENDIF
      ENDDO
C -------------------------------------------------------------------

C Find radial expansion of overlap ----------------------------------
      IF1 = INDF(IS1,IO1,1)
      IF2 = INDF(IS2,IO2,IOPER)

      IF ( INDFF(IF1,IF2,IOPER) .EQ. 0 ) THEN
*       CALL TIMER( 'MATEL2', 1 )

        NILM1 = NLM(IS1,IO1,1)
        NILM2 = NLM(IS2,IO2,IOPER)
        DO IFLM1 = IF1,IF1+NILM1-1
        DO IFLM2 = IF2,IF2+NILM2-1

C         Check interaction range
          IF (RCUT(IS1,IO1)+RCUT(IS2,IO2) .GT. RMAX) THEN
            call die('MATEL: NQ too small for required cutoff.')
          ENDIF

C         Find orbitals convolution by multiplication in k-space
          C = ( 2.D0 * PI )**1.5D0
          DO IQ = 0,NQ
            FFQ(IQ) = C * F(IQ,IFLM1) * F(IQ,IFLM2)
            IF ( OPERAT .EQ. 'T' ) THEN
              Q = IQ * QMAX / NQ
              FFQ(IQ) = FFQ(IQ) * Q * Q
            ENDIF
          ENDDO

C         Loop on possible values of l quantum number of product
          L1 = LOFILM( ILM(IFLM1) )
          L2 = LOFILM( ILM(IFLM2) )
          CALL RE_ALLOC( IFFR, 0,L1+L2, MYNAME//'IFFR' )
          CALL RE_ALLOC( CFFR, 0,L1+L2, MYNAME//'CFFR' )
          DO L3 = ABS(L1-L2), L1+L2, 2

C           Return to real space
            CALL RADFFT( L3, NQ, NQ*PI/RMAX, FFQ, FFL )
*           FFL(NQ) = 0.D0
            IF (MOD(ABS(L1-L2-L3)/2,2) .NE. 0) THEN
              DO IR = 0,NQ
                FFL(IR) = - FFL(IR)
              ENDDO
            ENDIF

C           Divide by R**L
            IF (L3 .NE. 0) THEN
              DO IR = 1,NQ
                R = IR * RMAX / NQ
                FFL(IR) = FFL(IR) / R**L3
              ENDDO
C             Parabolic extrapolation to R=0
              FFL(0) = ( 4.D0 * FFL(1) - FFL(2) ) / 3.D0
            ENDIF

C           Select NR out of NQ points
C           Copy NQ to a variable, to fool the compiler
            N = NQ
            IF (MOD(N,NR) .NE. 0)
     .           call die('MATEL: NQ must be multiple of NR')
            DO IR = 0,NR
              JR = IR * NQ / NR
              FFL(IR) = FFL(JR)
            ENDDO

C           Find if radial function is already in table
            FOUND = .FALSE.
            SEARCH: DO JO1 = LBOUND(INDF,2),UBOUND(INDF,2)
                    DO JO2 = LBOUND(INDF,2),UBOUND(INDF,2)
              JF1 = INDF(IS1,JO1,1)
              JF2 = INDF(IS2,JO2,IOPER)
              IF (JF1.NE.0 .AND. JF2.NE.0) THEN
                NJLM1 = NLM(IS1,JO1,1)
                NJLM2 = NLM(IS2,JO2,IOPER)
                DO JFLM1 = JF1,JF1+NJLM1-1
                DO JFLM2 = JF2,JF2+NJLM2-1
                  JFF = INDFF(JFLM1,JFLM2,IOPER)
                  IF (JFF .NE. 0) THEN
                    DO JFFY = INDFFY(JFF-1)+1, INDFFY(JFF)
                      JFFR = INDFFR(JFFY)
                      IF ( PROPOR(NR,FFL(1),FFR(1,1,JFFR),
     .                            FFTOL,CPROP)           ) THEN
                        FOUND = .TRUE.
                        IFFR(L3) = JFFR
                        CFFR(L3) = CPROP
                        EXIT SEARCH
                      ENDIF
                    ENDDO
                  ENDIF
                ENDDO
                ENDDO
              ENDIF
            ENDDO
            ENDDO SEARCH

C           Store new radial function and setup spline interpolation
            IF (.NOT.FOUND) THEN
              NFFR = NFFR + 1
              IF (NFFR .GT. MFFR) THEN
                MFFR = EXPAND * NFFR
                CALL RE_ALLOC( FFR, 0,NR, 1,2, 1,MFFR, MYNAME//'FFR' )
              ENDIF
              IFFR(L3) = NFFR
              CFFR(L3) = 1.D0
              DO IR = 0,NR
                FFR(IR,1,NFFR) = FFL(IR)
              ENDDO
              DFFR0 = HUGE(1.D0)
              DFFRMX = 0.D0
              CALL SPLINE( RMAX/NR, FFR(0,1,NFFR), NR+1, DFFR0, DFFRMX,
     .                     FFR(0,2,NFFR) )

C JMS debug
c              write(41,'(/,a,4i6)') 'if1,if2,nffr =', if1,if2,nffr
c              write(41,'(i6,2e25.12)') 
c     .          (ir,ffr(ir,1,nffr),ffr(ir,2,nffr),ir=0,nr)

            ENDIF
          ENDDO

C         Reallocate some arrays
          NILM = (L1+L2+1)**2
          IF (NFFY+NILM .GT. MFFY) THEN
            MFFY = EXPAND * (NFFY+NILM)
            CALL RE_ALLOC( FFY,    1,MFFY, MYNAME//'FFY'    )
            CALL RE_ALLOC( ILMFF,  1,MFFY, MYNAME//'ILMFF'  )
            CALL RE_ALLOC( INDFFR, 1,MFFY, MYNAME//'INDFFR' )
          ENDIF
          CALL RE_ALLOC( Y,         1,NILM, MYNAME//'Y'    )
          CALL RE_ALLOC( DYDR, 1,3, 1,NILM, MYNAME//'DYDR' )

C         Expand the product of two spherical harmonics (SH) also in SH
          CALL YLMEXP( L1+L2, RLYLM, YLMYLM, ILM(IFLM1), ILM(IFLM2),
     .                 1, 1, 1.D0, NILM, ILMFF(NFFY+1:MFFY), 
     .                 FFY(NFFY+1:MFFY))

C         Loop on possible lm values of orbital product
          DO I = 1,NILM
            NFFY = NFFY + 1
            JLM = ILMFF(NFFY)
            L3 = LOFILM( JLM )
            INDFFR(NFFY) = IFFR(L3)
            FFY(NFFY) = FFY(NFFY) * CFFR(L3)
          ENDDO

          NFF = NFF + 1
          IF (NFF .GT. MFF) THEN
            MFF = EXPAND * NFF
            CALL RE_ALLOC( INDFFY, 0,MFF, MYNAME//'INDFFY' )
          ENDIF
          INDFFY(NFF) = NFFY
          INDFF(IFLM1,IFLM2,IOPER) = NFF
        ENDDO
        ENDDO

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

C     Find if orbitals are out of range and avoid R12=0
      X12(1) = R12(1)
      X12(2) = R12(2)
      X12(3) = R12(3)
      R = SQRT( X12(1)*X12(1) + X12(2)*X12(2) + X12(3)*X12(3) )
      FOUND = .FALSE.
      IF (R .GT. RCUT(IS1,IO1)+RCUT(IS2,IO2) ) THEN
        FOUND = .TRUE.
      ELSEIF (R .LT. TINY) THEN
        X12(3) = TINY
        R = SQRT( X12(1)*X12(1) + X12(2)*X12(2) + X12(3)*X12(3) )
      ENDIF

C     Find spherical harmonics times R**L
      IF (.NOT.FOUND) THEN
        IF1 = INDF(IS1,IO1,1)
        IF2 = INDF(IS2,IO2,IOPER)
        NILM1 = NLM(IS1,IO1,1)
        NILM2 = NLM(IS2,IO2,IOPER)
        DO IFLM1 = IF1,IF1+NILM1-1
        DO IFLM2 = IF2,IF2+NILM2-1
          IFF = INDFF(IFLM1,IFLM2,IOPER)
          LMAX = 0
          DO IFFY = INDFFY(IFF-1)+1, INDFFY(IFF)
            JLM = ILMFF(IFFY)
            LMAX = MAX( LMAX, LOFILM(JLM) )
          ENDDO
          CALL RLYLM( LMAX, X12, Y, DYDR )

C         Interpolate radial functions and obtain SH expansion
          DO IFFY = INDFFY(IFF-1)+1, INDFFY(IFF)
            JFFR = INDFFR(IFFY)
            CALL SPLINT( RMAX/NR, FFR(0,1,JFFR),
     .                   FFR(0,2,JFFR), NR+1, R, SR, DSRDR )
            JLM = ILMFF(IFFY)
            S12 = S12 + SR * FFY(IFFY) * Y(JLM)
            DO IX = 1,3
              DSDR(IX) = DSDR(IX) +
     .                   DSRDR * FFY(IFFY) * Y(JLM) * X12(IX) / R +
     .                   SR * FFY(IFFY) * DYDR(IX,JLM)
            ENDDO
          ENDDO
        ENDDO
        ENDDO
      ENDIF

*     CALL TIMER( 'MATEL3', 2 )
C -------------------------------------------------------------------

C Restore allocation defaults 
      CALL ALLOC_DEFAULT( RESTORE=OLDEFS )

C Stop time counter
*     CALL TIMER( MYNAME, 2 )

      END SUBROUTINE MATEL
