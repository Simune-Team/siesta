
      SUBROUTINE COLINEAR( COORPO, COLIN )
C **********************************************************************
C Checks if three points lyes in the same straight line (if they are 
C colinear).
C Coded by J. Junquera November'98
C **********************************************************************

      IMPLICIT NONE

      DOUBLE PRECISION
     .  COORPO(3,3)
 
      LOGICAL
     .  COLIN

C **********************************************************************
C REAL*8 COORPO(3,3)   : Coordinates of three points. COORPO(POINT,IX)
C LOGICAL COLIN        : True => Three points are colinear
C                        False=> Three points NOT colinear
C **********************************************************************

      INTEGER 
     .  IX
 
      DOUBLE PRECISION
     .  VEC1(3), VEC2(3), NORMAL(3), EPS

      EXTERNAL
     .  CROSSV

      DATA EPS /1.D-12/

      DO IX = 1,3
        VEC1(IX) = COORPO(2,IX) - COORPO(1,IX)
        VEC2(IX) = COORPO(3,IX) - COORPO(1,IX)
      ENDDO

      CALL CROSSV( VEC1, VEC2, NORMAL )
     
      IF( (ABS(NORMAL(1)) .LT. EPS) .AND. 
     .    (ABS(NORMAL(2)) .LT. EPS) .AND.
     .    (ABS(NORMAL(3)) .LT. EPS) ) COLIN = .TRUE.

      END


      
      
