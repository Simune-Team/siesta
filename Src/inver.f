C $Id: inver.f,v 1.3 2000/01/11 10:19:11 jgale Exp $

        SUBROUTINE INVER(A,B,N)
        implicit double precision (a-h,o-z)
        double precision A(N,N),B(N,N),X
        DO I=1,N
          DO J=1,N
            B(I,J)=A(I,J)
          ENDDO
        ENDDO
        DO I=1,N
          X=B(I,I)
          B(I,I)=1.0d0
          DO J=1,N
            B(J,I)=B(J,I)/X
          ENDDO
          DO K=1,N
            IF ((K-I).NE.0) THEN 
              X=B(I,K)
              B(I,K)=0.0d0
              DO J=1,N
                B(J,K)=B(J,K)-B(J,I)*X
              ENDDO
            ENDIF
          ENDDO
        ENDDO
        RETURN
        END
