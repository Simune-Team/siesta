      SUBROUTINE FERMID( NSPIN, MAXSPN, NK, WK, MAXE, NE, E, 
     .                   TEMP, QTOT, WKE, EF, ENTROPY )

C *********************************************************************
C Finds the Fermi energy and the occupation weights of states.
C Written by J.M.Soler. August'96.
C Simple single excitation introduced by E. Artacho August 2002.
C Alternative occupation functions introduced by P. Ordejon, May'03.
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
C REAL*8  ENTROPY            : Entropy contribution to the electronic
C                              Free Energy
C *********************************************************************

C
C  Modules
C
      use precision
      use fdf
#ifdef MPI
      use mpi_siesta
#endif

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      PARAMETER (TOL=1.D-10,NITMAX=150)
      DOUBLE PRECISION E(MAXE,MAXSPN,NK),WKE(MAXE,MAXSPN,NK),WK(NK)
      INTEGER Node, ief

C Local variables
      logical blread, excitd
#ifdef MPI
      integer MPIerror
#endif
      
      save blread, excitd
      data blread, excitd /.false., .false./


C Get the Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
#else
      Node = 0
#endif

c Reading whether excited state -------

      if ( .not. blread) then
         if ( Node .eq. 0 ) then
            excitd = fdf_boolean('SingleExcitation', .false.)
            if ( excitd ) write(6,'(/a)') 
     .          'fermid: Calculating for lowest-exciton excited state'
         endif
#ifdef MPI
         call MPI_Bcast(excitd,1,MPI_logical,0,
     .                  MPI_Comm_World,MPIerror)
#endif
         blread = .true.
      endif   
c--------------------------------------


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
      IF (ABS(SUMQ-QTOT).LT.TOL) THEN
        if (excitd) then
           if (Node.eq.0) then
              write (6,'(/a)') 
     .               'FERMID: BANDS FULL, NO EXCITATION POSSIBLE'
           endif
           stop
        else
           return
        endif
      ENDIF
      IF (SUMQ.LT.QTOT) THEN
        if (Node.eq.0) then
          WRITE (6,*) 'FERMID: NOT ENOUGH STATES'
          WRITE (6,*) 'FERMID: QTOT,SUMQ=',QTOT,SUMQ
        endif
        STOP
      ENDIF
      T=MAX(TEMP,1.D-6)
      DRANGE=T*SQRT(-LOG(TOL*.01D0))
      EMIN=EMIN-DRANGE
      EMAX=EMAX+DRANGE
      DO 50 ITER=1,NITMAX
         EF=0.5D0*(EMIN+EMAX)
         SUMQ=0.D0
         DO 40 IK=1,NK
           DO 35 ISPIN=1,NSPIN
             DO 30 IE=1,NE
               WKE(IE,ISPIN,IK)=WK(IK)*
     .             STEPF((E(IE,ISPIN,IK)-EF)/T)*2.0/NSPIN
               SUMQ=SUMQ+WKE(IE,ISPIN,IK)
  30         CONTINUE
  35       CONTINUE
  40     CONTINUE

c If the Fermi level was found..................... 
         IF (ABS(SUMQ-QTOT).LT.TOL) THEN

c If excited state is to be calculated, find the level above Ef for
c k=1 and spin=1
            if (excitd) then
               do ie = 1, ne
                  if ( e(ie,1,1) .gt. EF ) then
                     ief = ie
                     go to 314
                  endif
               enddo
  314          continue
c and swap populations (meaningful only for T close to 0):
c if nspin =1 populations of homo and lumo are just swapped for is=1
c if nspin =2 populations of homo and lumo become equal.
               wkebuf = ( wke(ief-1,1,1) - wke(ief,1,1) )/(3-nspin)
               wke(ief,1,1) = wke(ief,1,1) + wkebuf
               wke(ief-1,1,1) = wke(ief-1,1,1) - wkebuf
            endif


C Obtain the electronic entropy

100         ENTROPY=0.0D0

            DO 140 IK=1,NK
              DO 135 IE=1,NE
                DO 130 ISPIN=1,NSPIN

	          W = (NSPIN / 2.0D0) * WKE(IE,ISPIN,IK) / WK(IK)
	          EIK = (E(IE,ISPIN,IK)-EF) / T

                  ENTROPY = ENTROPY + ( 2.0D0 * WK(IK) / NSPIN ) *
     .                      ENPY(EIK,W)

130             CONTINUE
135           CONTINUE
140         CONTINUE

            RETURN

         ENDIF
c...................................................
         IF (SUMQ.LE.QTOT) EMIN=EF
         IF (SUMQ.GE.QTOT) EMAX=EF
  50  CONTINUE
      if (Node.eq.0) then
        WRITE (6,*) 'FERMID: ITERATION HAS NOT CONVERGED.'
        WRITE (6,*) 'FERMID: QTOT,SUMQ=',QTOT,SUMQ
        STOP 'FERMID: ITERATION HAS NOT CONVERGED.'
      else
        STOP
      endif
      END


      DOUBLE PRECISION FUNCTION STEPF(X)
      use fdf
#ifdef MPI
      use mpi_siesta
#endif

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C Local variables
      character(len=22), save :: ocf = 'FD'
      integer,           save :: nh = 1
      integer                 :: Node
      integer,           save :: ocupfnct
      logical                 :: leqi
      logical,           save :: ocfread = .false.
#ifdef MPI
      integer MPIerror
#endif
      
      PARAMETER (PI = 3.14159265358979D0)

c Reading which electronic occupation function to use -------

      if ( .not. ocfread) then

C Get the Node number
#ifdef MPI
        call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
#else
        Node = 0
#endif

        if ( Node .eq. 0 ) then
          ocf = fdf_string('OccupationFunction','FD')

          if (leqi(ocf,'FD')) then
            ocupfnct=1
            write(6,'(/a)') 
     .     'stepf: Fermi-Dirac step function'
          else if (leqi(ocf,'MP')) then
            ocupfnct=2
            nh = fdf_integer('OccupationMPOrder',1)
            write(6,'(/a,i2)') 
     .     'stepf: Methfessel-Paxton step function'
            write(6,'(a,i2)') 
     .     '       Using Hermite-Gauss polynomials of order ',nh
          else
            stop 'fermid: Error: Allowed values for OccupationFunction 
     . are FD and MP'
          endif
        endif
#ifdef MPI
        call MPI_Bcast(ocupfnct,1,MPI_integer,0,
     .                 MPI_Comm_World,MPIerror)
        call MPI_Bcast(nh,1,MPI_integer,0,
     .                 MPI_Comm_World,MPIerror)
#endif
        ocfread = .true.
      endif   
c--------------------------------------



C     Complementary error function. Ref: Fu & Ho, PRB 28, 5480 (1983)
*     STEPF=DERFC(X)  -  not available

      IF (OCUPFNCT .EQ. 1) THEN

C     Fermi-Dirac distribution
        IF (X.GT.100.D0) THEN
          STEPF = 0.D0
        ELSEIF (X.LT.-100.D0) THEN
          STEPF = 1.D0
        ELSE
          STEPF = 1.D0 / ( 1.D0 + EXP(X) )
        ENDIF


      ELSE IF (OCUPFNCT .EQ. 2) THEN

C     Improved step function. Ref: Methfessel & Paxton PRB40 (15/Aug/89)
C     NH is the order of the Hemite polynomial expansion.



        STEPF =  0.5D0 * DERFC(X)
        A = 1.0D0/SQRT(PI)

        DO I=1,NH

C Get coefficients in Hermite-Gauss expansion
          A = -A / (I * 4.0D0)

C Get contribution to step function at order I
          GAUSS = DEXP(-X*X)
          J = 2*I -1 
          IF (GAUSS .GT. 1.D-20) STEPF = STEPF + A * HP(X,J) * GAUSS

        ENDDO

      ELSE

        STOP 'STEPF: INCORRECT STEP FUNCTION'

      ENDIF

      END



      DOUBLE PRECISION FUNCTION ENPY(E,W)
      use fdf
#ifdef MPI
      use mpi_siesta
#endif

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C Computes the contribution of a given state with energy E
C (refered to the Fermi energy, in units of the smearing 
C temperature) and occupation W to the electronic entropy
C P. Ordejon, June 2003

C Local variables
      character(len=22), save :: ocf = 'FD'
      integer,           save :: nh = 1
      integer                 :: Node 
      integer,           save :: ocupfnct
      logical                 :: leqi
      logical,           save :: ocfread = .false.
#ifdef MPI
      integer MPIerror
#endif
      
      DATA tiny /1.d-15/

c Reading which electronic occupation function to use -------

      if ( .not. ocfread) then

C Get the Node number
#ifdef MPI
        call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
#else
        Node = 0
#endif

        if ( Node .eq. 0 ) then
          ocf = fdf_string('OccupationFunction', 'FD')

          if (leqi(ocf,'FD')) then
            ocupfnct=1
          else if (leqi(ocf,'MP')) then
            ocupfnct=2
            nh = fdf_integer('OccupationMPOrder',1)
          else
            stop 'fermid: Error: Allowed values for OccupationFunction 
     . are FD and MP'
          endif
        endif
#ifdef MPI
        call MPI_Bcast(ocupfnct,1,MPI_integer,0,
     .                  MPI_Comm_World,MPIerror)
        call MPI_Bcast(nh,1,MPI_integer,0,
     .                  MPI_Comm_World,MPIerror)
#endif
        ocfread = .true.
      endif   
c--------------------------------------



      IF (OCUPFNCT .EQ. 1) THEN

C     Mermin entropy for the Fermi-Dirac distribution
        wo = max( w, tiny )
        we = 1.0d0 - wo
        we = max( we, tiny )

        ENPY = - 1.0 * ( wo*log(wo) + we*log(we) )
     .         

      ELSE IF (OCUPFNCT .EQ. 2) THEN

C     Entropy for the Improved step function. 
C     Ref: Methfessel & Paxton PRB40 (15/Aug/89)

        ENPY = WHG(E,NH)

      ELSE

        STOP 'STEPF: INCORRECT STEP FUNCTION'

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




      DOUBLE PRECISION FUNCTION WHG(X,N)
C
C  Computes the factors to get the entropy term 
C  for the Methfessel-Paxton smearing with Hermite
C  polynomials of order N
C
C  P. Ordejon, June '03
C

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      PI = 3.14159265358979D0
      X2=X**2.0D0


C get coefficients

      A = 1.0D0/SQRT(PI)
      DO I=1,N
        A = -A / (I * 4.0D0)
      ENDDO

      GAUSS = DEXP(-X2)
      WHG = 0.D0
      IF (GAUSS .GT. 1.D-20) WHG = 0.5D0 * A * HP(X,2*N) * GAUSS

      RETURN
      END


      DOUBLE PRECISION FUNCTION HP(X,N)

C 
C  Returns the value of the Hermite polynomial of degree N
C  evaluated at X.
C
C  H_0  (x) = 1
C  H_1  (x) = 2x
C  ...
C  H_n+1(x) = 2 x H_n(x) - 2 n H_n-1(x)
C
C  P. Ordejon, June 2003
C

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N

      IF (N .GT. 1000) 
     .     stop 'fermid: Order of Hermite polynomial too large'


      HP = 1.0D0
      HM2 = 0.0D0
      HM1 = 1.0D0

      DO I=1,N

        HP = 2.0 * (X * HM1 - (I-1) * HM2)
        HM2 = HM1
        HM1 = HP

      ENDDO
      
      END

      
      
