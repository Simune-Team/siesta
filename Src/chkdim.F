      SUBROUTINE CHKDIM (SUB,VAR,ND,N,IOPT)
C **********************************************************************
C Checks a dimension parameter, warning and stopping if it is too small
C Written by J.M.Soler. 1987.
C ************ INPUT ***************************************************
C CHARACTER SUB*(*) : Name of calling subroutine
C CHARACTER VAR*(*) : Name of dimension parameter
C INTEGER   ND      : Actual value of dimension parameter
C INTEGER   N       : Required value of dimension parameter
C INTEGER   IOPT    : Option switch: IOPT=0 => Require that ND.EQ.N
C                                    IOPT=1 => Require that ND.GE.N
C **********************************************************************
C
C  Modules
C
      use precision
#ifdef MPI
      use mpi
#endif
C
#ifdef MPI
      integer 
     .  MPIerror,MPIstop
#endif
      character SUB*(*),VAR*(*)

      IF ( IOPT.EQ.0 ) THEN
        IF ( ND.EQ.N ) RETURN
        WRITE (6,'(/5A,I8,A,I8)') 'chkdim: ERROR: In ', SUB,
     .        ', dimension ',VAR,' =',ND,'. It must be exactly ',N
#ifdef MPI
        MPIstop = 1
        call MPI_Abort(MPI_Comm_World,MPIstop,MPIerror)
#else
        STOP
#endif
      ELSE
        IF ( ND.GE.N ) RETURN
        WRITE (6,'(/5A,I8,A,I8)') 'chkdim: ERROR: In ', SUB,
     .        ', dimension ',VAR,' =',ND,'. It must be at least ',N
#ifdef MPI
        MPIstop = 1
        call MPI_Abort(MPI_Comm_World,MPIstop,MPIerror)
#else
        STOP
#endif
      ENDIF
      END



      SUBROUTINE CHKDIME (ND,N,OVERFLOW,NM)
      LOGICAL OVERFLOW
      NM = MAX(N,NM)
      IF ( ND.GE.N ) RETURN
      OVERFLOW = .TRUE.
      RETURN
      END

