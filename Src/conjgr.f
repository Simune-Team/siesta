C $Id: conjgr.f,v 1.5 2004/07/15 17:16:55 wdpgaara Exp $

      SUBROUTINE CONJGR (N,X,G,DXMAX,GTOL,CNTROL,H)

      use precision

C DIRECTS A CONJUGATE-GRADIENT MINIMIZATION OF A FUNCTION WHICH
C IS EVALUATED BY THE CALLING PROGRAM.
C  N     : INPUT SPACE DIMENSIONALITY
C  X     : INPUT POSITION AT WHICH GRADIENT HAS BEEN EVALUATED AND
C          OUTPUT NEW POSITION AT WHICH GRADIENT MUST BE EVALUATED NEXT
C  G     : INPUT GRADIENT (WITH A MINUS SIGN, OR ACTIVATE A LINE BELOW)
C  DXMAX : INPUT MAXIMUM ALLOWED DISPLACEMENT IN EACH COORDINATE
C  GTOL  : INPUT MAXIMUM FINAL VALUE OF EACH GRADIENT COMPONENT
C  CNTROL: CONTROL ARRAY. FIRST ELEMENT MUST BE MADE ZERO BEFORE
C          FIRST CALL. IF IT IS ZERO ON OUTPUT, MINIMIZATION IS
C          CONVERGED. OTHERWISE, CALCULATE GRADIENT AT THE NEW
C          POSITION AND CALL AGAIN THIS ROUTINE. DO NOT MODIFY ANY
C          ARGUMENT OTHER THAN G BETWEEN CALLS.
C  H     : AUXILIARY ARRAY WHICH MUST NOT BE MODIFIED BETWEEN CALLS
C  IOPT  : PARAMETER BELOW WHICH DETERMINES METHOD USED AND AUXILIARY
C          STORAGE REQUIRED: IOP=1 => FLETCHER-REEVES. IOPT=2 =>
C          POLAK-RIBIERE. DETAILS IN SECT. 10.6 OF 'NUMERICAL RECIPES'
C WRITTEN BY J.SOLER. JAN/91. BASED ON ROUTINES IN 'NUMERICAL RECIPES'

      integer, intent(in)   :: n
      integer,  PARAMETER   :: IOPT=2

      real(dp) ::  X(N),G(N),H(N,IOPT),CNTROL(0:19)
      real(dp), intent(in)  :: dxmax, gtol

      real(dp) :: gmax, gg, gamma
      integer  :: i, j

      real(dp) :: dot
      external  dot

C     IF GRADIENT IS SMALLER THAN TOLERENCE, RETURN
      GMAX=ABS(G(1))
      DO 10 J=1,N
*       G(J)=-G(J)
        GMAX=MAX(GMAX,ABS(G(J)))
   10 CONTINUE
      IF (GMAX.LE.GTOL) THEN
        CNTROL(0)=0
        GOTO 60
      ENDIF

C     FIRST-CALL INITIALIZATIONS
      IF (NINT(CNTROL(0)).EQ.0) THEN
        DO 30 I=1,IOPT
          DO 20 J=1,N
            H(J,I)=G(J)
   20     CONTINUE
   30   CONTINUE
        CNTROL(0)=1
        CNTROL(1)=1
        CNTROL(2)=DOT(G,G,N)
        CNTROL(10)=0
        CNTROL(18)=DXMAX
      ENDIF

C     LINE MINIMIZATION IS ALWAYS CALLED
   40 CALL LINMIN (N,X,H,G,DXMAX,CNTROL(10))

C     IF LINE MINIMIZATION IS FINISHED, FIND NEW LINE DIRECTION
      IF (NINT(CNTROL(10)).EQ.0) THEN
        GG=DOT(G,G,N)
        IF (IOPT.EQ.2) GG=GG-DOT(G,H(1,2),N)
        GAMMA=GG/CNTROL(2)
        DO 50 J=1,N
          H(J,1)=G(J)+GAMMA*H(J,1)
          IF (IOPT.EQ.2) H(J,2)=G(J)
   50   CONTINUE
        CNTROL(1)=CNTROL(1)+1
        CNTROL(2)=DOT(G,G,N)
*       WRITE(6,'(A,I4,F15.6)')
*    .     ' CONJGR: NEW LINE DIRECTION. N,DX=',N,CNTROL(18)
        GOTO 40
      ENDIF

   60 CONTINUE
      END

      SUBROUTINE LINMIN (N,XVEC,HVEC,GVEC,DXMAX,CNTROL)

      use precision

      integer, intent(in) :: n
      real(dp) ::  XVEC(N),HVEC(N),GVEC(N),CNTROL(0:9)
      real(dp), intent(in) :: dxmax

      real(dp), PARAMETER :: FACTOR=1.6D0

      integer  :: i, icntrl
      real(dp) :: x1, x2, y1, y2, x0, hmod, hmax, dx, x, y

      real(dp) :: dot
      external  dot

C     TRANSLATE CONTROL PARAMETERS
      ICNTRL=NINT(CNTROL(0))
      X1=CNTROL(1)
      X2=CNTROL(2)
      Y1=CNTROL(3)
      Y2=CNTROL(4)
      X0=CNTROL(5)
      HMOD=CNTROL(6)
      HMAX=CNTROL(7)
      DX=CNTROL(8)
      X=X0
      Y=DOT(GVEC,HVEC,N)
*     WRITE(6,'(A,I4,2F12.6)') ' LINMIN: ICNTRL,X,Y=',ICNTRL,X,Y

      IF (ICNTRL.EQ.0) THEN
C       INITIALIZE X1,Y1 ON FIRST CALL
        X1=0.D0
        Y1=Y
C       PREPARE SECOND POINT
        ICNTRL=1
        X0=0.D0
        IF (DX.EQ.0.D0) DX=DXMAX
        HMOD=SQRT(DOT(HVEC,HVEC,N))
        HMAX=0.D0
        DO 10 I=1,N
          HMAX=MAX(HMAX,ABS(HVEC(I)))
   10   CONTINUE
        X=MIN(DX/HMOD,DXMAX/HMAX)
        GOTO 20
      ELSEIF (ICNTRL.EQ.1) THEN
C       INITIALIZE X2,Y2 ON SECOND CALL
        X2=X
        Y2=Y
        ICNTRL=2
      ELSEIF (ICNTRL.EQ.2) THEN
C       SHIFT INTERVAL USING NEW POINT
        X1=X2
        Y1=Y2
        X2=X
        Y2=Y
      ELSEIF (ICNTRL.EQ.3) THEN
C       IF ROOT WAS FOUND IN LAST CALL, ALL IS DONE NOW
        ICNTRL=0
        GOTO 20
      ENDIF

      IF (Y2.GT.0.D0) THEN
C       ROOT NOT BRACKETED YET. TRY NEW RIGHT BRACKET
        X=X2+MIN(FACTOR*(X2-X1),DXMAX/HMAX)
      ELSE
C       INTERPOLATE FOR ROOT AND RETURN TO CALCULATE LAST GRADIENT
        X=(X1*Y2-X2*Y1)/(Y2-Y1)
        ICNTRL=3
      ENDIF

C     STORE CONTROL PARAMETERS AND SET NEW POINT
   20 CNTROL(0)=ICNTRL
      CNTROL(1)=X1
      CNTROL(2)=X2
      CNTROL(3)=Y1
      CNTROL(4)=Y2
      CNTROL(5)=X
      CNTROL(6)=HMOD
      CNTROL(7)=HMAX
      CNTROL(8)=ABS(X)*HMOD
      DO 30 I=1,N
        XVEC(I)=XVEC(I)+HVEC(I)*(X-X0)
   30 CONTINUE
      END











