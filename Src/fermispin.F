      SUBROUTINE FERMISPIN( NSPIN, MAXSPN, NK, WK, MAXE, NE, E, 
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
C REAL*8  QTOT(MAXSPN) : Total valence charge (number of electrons)
C                         for each spin component
C ********** OUTPUT ***************************************************
C REAL*8  WKE(MAXE,MAXSPN,NK) : Occupations multiplied by k-point weights
C                               (sum QTOT)
C REAL*8  EF(NSPIN)           : Fermi energy (for each spin, if QTOT
C                              is different for each spin component.
C *********************************************************************

C
C  Modules
C
      use precision
#ifdef MPI
      use mpi_siesta
#endif

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION E(MAXE,MAXSPN,NK),EMIN(4),EMAX(4),
     .                 EF(NSPIN),QTOT(MAXSPN),SUMQ(4),
     .                 WKE(MAXE,MAXSPN,NK),WK(NK)
      INTEGER          Node
      LOGICAL          CONV
      PARAMETER (TOL=1.0D-10,NITMAX=50)

C Local variables
#ifdef MPI
      integer MPIerror
#endif

C Get the Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
#else
      Node = 0
#endif

      CONV = .FALSE.

      DO ISPIN = 1,NSPIN
       SUMQ(ISPIN)=0.0D0
      ENDDO
      DO ISPIN = 1,NSPIN
        EMIN(ISPIN)=E(1,ISPIN,1)
        EMAX(ISPIN)=E(1,ISPIN,1)
      ENDDO
      DO 20 IK=1,NK
        DO 15 ISPIN=1,NSPIN
          DO 10 IE=1,NE
            WKE(IE,ISPIN,IK)=WK(IK)*2.D0/NSPIN
            SUMQ(ISPIN)=SUMQ(ISPIN)+WKE(IE,ISPIN,IK)
            EMIN(ISPIN)=MIN(EMIN(ISPIN),E(IE,ISPIN,IK))
            EMAX(ISPIN)=MAX(EMAX(ISPIN),E(IE,ISPIN,IK))
  10      CONTINUE
  15    CONTINUE
  20  CONTINUE
      DO ISPIN=1,NSPIN
        EF(ISPIN)=EMAX(ISPIN)
      ENDDO
      CONV = .TRUE.
      DO ISPIN = 1,NSPIN
        IF (ABS(SUMQ(ISPIN)-QTOT(ISPIN)).GT.TOL) CONV = .FALSE.
      ENDDO
      IF (CONV) RETURN
      DO ISPIN = 1,NSPIN
        IF (SUMQ(ISPIN).LT.QTOT(ISPIN)) THEN
          if (Node.eq.0) then
            WRITE (6,*) 'FERMID: NOT ENOUGH STATES'
            WRITE (6,*) 'FERMID: ISPIN,QTOT,SUMQ=',
     .                 ISPIN,QTOT(ISPIN),SUMQ(ISPIN)
          endif
          STOP
        ENDIF
      ENDDO
      T=MAX(TEMP,1.D-6)
      DRANGE=T*SQRT(-LOG(TOL*.01D0))
      DO ISPIN = 1,NSPIN
        EMIN(ISPIN)=EMIN(ISPIN)-DRANGE
        EMAX(ISPIN)=EMAX(ISPIN)+DRANGE
      ENDDO
      DO 50 ITER=1,NITMAX
        DO ISPIN = 1,NSPIN
          EF(ISPIN)=0.5D0*(EMIN(ISPIN)+EMAX(ISPIN))
          SUMQ(ISPIN)=0.D0
        ENDDO
        DO 40 IK=1,NK
          DO 35 ISPIN=1,NSPIN
            DO 30 IE=1,NE
              WKE(IE,ISPIN,IK)=WK(IK)*
     .             STEPF((E(IE,ISPIN,IK)-EF(ISPIN))/T)*2.0/NSPIN
              SUMQ(ISPIN)=SUMQ(ISPIN)+WKE(IE,ISPIN,IK)
  30        CONTINUE
  35      CONTINUE
  40    CONTINUE
        CONV = .TRUE.
        DO ISPIN = 1,NSPIN
          IF (ABS(SUMQ(ISPIN)-QTOT(ISPIN)).GT.TOL) CONV = .FALSE.
        ENDDO
        IF (CONV) RETURN
        DO ISPIN = 1,NSPIN
          IF (SUMQ(ISPIN).LE.QTOT(ISPIN)) EMIN(ISPIN)=EF(ISPIN)
          IF (SUMQ(ISPIN).GE.QTOT(ISPIN)) EMAX(ISPIN)=EF(ISPIN)
        ENDDO
  50  CONTINUE
      if (Node.eq.0) then
        WRITE (6,*) 'FERMID: ITERATION HAS NOT CONVERGED.'
        DO ISPIN = 1,NSPIN
          WRITE (6,*) 'FERMID: ISPIN,QTOT,SUMQ=',
     .               ISPIN,QTOT(ISPIN),SUMQ(ISPIN)
        ENDDO 
        STOP 'FERMID: ITERATION HAS NOT CONVERGED.'
      else
        STOP
      endif
      END

