      INTEGER FUNCTION LOFILM( ILM )

C *********************************************************************
C Finds the total angular momentum quantum number L which corresponds
C to the combined index ILM == (L,M) implicit in the nested loop
C    ILM = 0
C    DO L=0,LMAX
C      DO M=-L,L
C        ILM = ILM + 1
C with ILM = 1,2,...,L**2
C Written by J.M.Soler. April 1996.
C *********************************************************************

C
C  Modules
C
      use precision
#ifdef MPI
      use mpi
#endif

      IMPLICIT NONE

#ifdef MPI
      integer 
     .  MPIerror
#endif

      INTEGER ILM, JLM, MAXL, Node
      PARAMETER ( MAXL = 100 )

C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
#else
      Node = 0
#endif

      if ( ILM .LE. 0 ) then
        if (Node.eq.0) then
          STOP 'LOFILM: ILM not allowed'
        else
          STOP
        endif
      endif
      JLM = 0
      DO 10 LOFILM = 0,MAXL
        JLM = JLM + 2*LOFILM + 1
        IF ( JLM .GE. ILM ) RETURN
   10 CONTINUE
      if (Node.eq.0) then
        STOP 'LOFILM: ILM too large'
      else
        STOP
      endif
      END

