C $Id: fermid.f,v 1.5 1999/01/31 11:14:53 emilio Exp $

      SUBROUTINE FERMID( NSPIN, MAXSPN, NK, WK, MAXE, NE, E, 
     .                   TEMP, QTOT, WKE, EF )

C *********************************************************************
C Finds the Fermi energy and the occupation weights of states.
C Written by J.M.Soler. August'96.
C ********** INPUT ****************************************************
C INTEGER NSPIN    : Number of different spin polarizations (1 or 2)
C INTEGER MAXSPN   : Maximum number of different spin polarizations (1 or 2)
C                    for E and WKE matrices dimensions
C INTEGER NK       : Number of K-points
C REAL*8  WK(NK)   : Sampling weights of k-points (must sum 1)
C INTEGER MAXE     : First dimension of E and WKE
C INTEGER NE       : Number of bands
C REAL*8  E(MAXE,MAXSPN,NK) : State eigenvalues
C REAL*8  TEMP     : Temperature (in the same units of E)
C REAL*8  QTOT     : Total valence charge (number of electrons)
C ********** OUTPUT ***************************************************
C REAL*8  WKE(MAXE,MAXSPN,NK) : Occupations multiplied by k-point weights
C                               (sum QTOT)
C REAL*8  EF                 : Fermi energy
C *********************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION E(MAXE,MAXSPN,NK),WKE(MAXE,MAXSPN,NK),WK(NK)
      PARAMETER (TOL=1.D-10,NITMAX=50)

      SUMQ=0.D0
      EMIN=E(1,1,1)
      EMAX=E(1,1,1)
      DO 20 IK=1,NK
        DO 15 ISPIN=1,NSPIN
          DO 10 IE=1,NE
            WKE(IE,ISPIN,IK)=WK(IK)*2.D0/NSPIN
            SUMQ=SUMQ+WKE(IE,ISPIN,IK)
            EMIN=MIN(EMIN,E(IE,ISPIN,IK))
            EMAX=MAX(EMAX,E(IE,ISPIN,IK))
  10      CONTINUE
  15    CONTINUE
  20  CONTINUE
      EF=EMAX
      IF (ABS(SUMQ-QTOT).LT.TOL) RETURN
      IF (SUMQ.LT.QTOT) THEN
         WRITE (6,*) 'FERMID: NOT ENOUGH STATES'
         WRITE (6,*) 'FERMID: QTOT,SUMQ=',QTOT,SUMQ
         STOP
      ENDIF
      T=MAX(TEMP,1.D-6)
      DRANGE=T*SQRT(-LOG(TOL*.01D0))
      EMIN=EMIN-DRANGE
      EMAX=EMAX+DRANGE
      DO 50 ITER=1,NITMAX
         EF=0.5D0*(EMIN+EMAX)
         SUMQ=0.D0
         SUME=0.D0
         DO 40 IK=1,NK
           DO 35 ISPIN=1,NSPIN
             DO 30 IE=1,NE
               WKE(IE,ISPIN,IK)=WK(IK)*
     .                          STEPF((E(IE,ISPIN,IK)-EF)/T)/NSPIN
               SUMQ=SUMQ+WKE(IE,ISPIN,IK)
  30         CONTINUE
  35       CONTINUE
  40     CONTINUE
         IF (ABS(SUMQ-QTOT).LT.TOL) RETURN
         IF (SUMQ.LE.QTOT) EMIN=EF
         IF (SUMQ.GE.QTOT) EMAX=EF
  50  CONTINUE
      WRITE (6,*) 'FERMID: ITERATION HAS NOT CONVERGED.'
      WRITE (6,*) 'FERMID: QTOT,SUMQ=',QTOT,SUMQ
      STOP 'FERMID: ITERATION HAS NOT CONVERGED.'
      END



      DOUBLE PRECISION FUNCTION STEPF(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C     Complementary error function. Ref: Fu & Ho, PRB 28, 5480 (1983)
*     STEPF=DERFC(X)

C     Improved step function. Ref: Methfessel & Paxton PRB40 (15/Aug/89)
*     PARAMETER (C=0.5641895835D0)
*     STEPF=DERFC(X)-C*X*DEXP(-X*X)

C     Fermi-Dirac distribution
      IF (X.GT.100.D0) THEN
        STEPF = 0.D0
      ELSEIF (X.LT.-100.D0) THEN
        STEPF = 2.D0
      ELSE
        STEPF = 2.D0 / ( 1.D0 + EXP(X) )
      ENDIF

      END



      DOUBLE PRECISION FUNCTION DERFC (X)

C  COMPLEMENTARY ERROR FUNCTION FROM "NUMERICAL RECIPES"
C  NOTE: SINGLE PRECISION ACCURACY

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      Z=ABS(X)
      T=1.D0/(1.D0+0.5D0*Z)
      DERFC=T*EXP(-(Z*Z)-1.26551223D0+T*(1.00002368D0+T*(0.37409196D0+
     .      T*(0.09678418D0+T*(-0.18628806D0+
     .      T*(0.27886807D0+T*(-1.13520398D0+
     .      T*(1.48851587D0+T*(-0.82215223D0+T*.17087277D0)))))))))
      IF (X.LT.0.D0) DERFC=2.D0-DERFC
      END



      DOUBLE PRECISION FUNCTION DERF (X)

C  ERROR FUNCTION FROM "NUMERICAL RECIPES"
C  NOTE: SINGLE PRECISION ACCURACY

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      Z=ABS(X)
      T=1.D0/(1.D0+0.5D0*Z)
      DERF= T*EXP(-(Z*Z)-1.26551223D0+T*(1.00002368D0+T*(0.37409196D0+
     .      T*(0.09678418D0+T*(-0.18628806D0+
     .      T*(0.27886807D0+T*(-1.13520398D0+
     .      T*(1.48851587D0+T*(-0.82215223D0+T*.17087277D0)))))))))
      IF (X.LT.0.D0) DERF=2.D0-DERF

      DERF = 1.D0 - DERF
      END
