*************************************************************
        SUBROUTINE INVER(A,B,N)
        implicit double precision (a-h,o-z)
        double precision A(N,N),B(N,N),X
        DO 20 I=1,N
        DO 20 J=1,N
20      B(I,J)=A(I,J)
        DO 4 I=1,N
        X=B(I,I)
        B(I,I)=1.
        DO 1 J=1,N
1       B(J,I)=B(J,I)/X
        DO 4 K=1,N
        IF(K-I) 2,4,2
2       X=B(I,K)
        B(I,K)=0.0d0
        DO 3 J=1,N
3       B(J,K)=B(J,K)-B(J,I)*X
4       CONTINUE
        RETURN
        END
