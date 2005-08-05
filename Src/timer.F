      SUBROUTINE TIMER (PROG,IOPT)
      character(len=*) prog
      integer iopt
c
c     Two-in-one: CPU plus elapsed times 
c     (elapsed times go to file WALLTIME)
c
      call timer1(prog,iopt)
      call elap1(prog,iopt)

      end


      SUBROUTINE TIMER1 (PROG,IOPT)


C  FINDS AND PRINTS THE CPU TIME SPENT IN DIFFERENT ROUTINES AND/OR
C   PROGRAM SECTIONS. IT MUST BE CALLED WITH IOPT=0 AT THE BEGINNING
C   OF EACH ROUTINE AND WITH IOPT=1 AT THE END OF IT.
C  ARGUMENTS:
C    PROG: INPUT ARBITRARY NAME FOR THE ROUTINE AND/OR PROGRAM SECTION
C    IOPT: INPUT OPTION PARAMETER:
C      IOPT = 0  => SET ALL TIMES TO ZERO (OPTIONAL)
C      IOPT = 1  => BEGIN COUNTING TIME FOR A ROUTINE
C      IOPT = 2  => STOP  COUNTING TIME FOR A ROUTINE
C      IOPT = 3  => PRINT TIME FOR A ROUTINE OR FOR ALL (IF PROG='ALL')
C  ROUTINE TIMES INCLUDE THAT SPENT IN THE ROUTINES THEY CALL
C  WRITTEN BY J.SOLER (JSOLER AT EMDUAM11) DEC/90

C
C  Modules
C
      use precision
      use parallel,   only : Node
#ifdef MPI
      use mpi_siesta
#endif

      implicit real(dp)(A-H,O-Z)

      real :: treal

      integer :: nmax, iopt, nprogs, iprog
      PARAMETER (NMAX=500,ZERO=0.D0,HUNDRD=100.D0,TIMMIN=1.D-6)
      DIMENSION TIME1(NMAX),TIMET(NMAX)
#ifdef MPI
      dimension buffer1(3),buffer2(3)
      integer mpierror
#endif
      INTEGER NCALLS(NMAX)
      CHARACTER*10 PROGS(NMAX),PROG*(*)
      SAVE PROGS,NPROGS,TIME0,TIME1,TIMET,NCALLS
      DATA NPROGS,TIME0 /0,ZERO/

!!      CALL CPUTIM (TIME)
      call cpu_time( treal)         ! Standard Fortran95
      TIME = treal

      IF (IOPT.EQ.0) THEN
         NPROGS=0
         TIME0=TIME
      ELSEIF (IOPT.EQ.1 .OR. IOPT.EQ.2) THEN
         DO 10 IPROG=1,NPROGS
            IF (PROGS(IPROG).EQ.PROG) GO TO 20
   10    CONTINUE
            NPROGS=NPROGS+1
            IF (NPROGS.GT.NMAX) THEN
               if (Node.eq.0) then
                 WRITE (6,*) 'timer: NMAX IS SATURATED. PROG = ',PROG
               endif
               RETURN
            ENDIF
            IPROG=NPROGS
            PROGS(IPROG)=PROG
            NCALLS(IPROG)=0
            TIMET(IPROG)=ZERO
   20    CONTINUE
         IF (IOPT.EQ.1) THEN
            NCALLS(IPROG)=NCALLS(IPROG)+1
            TIME1(IPROG)=TIME
         ELSE
            TIMET(IPROG)=TIMET(IPROG)+TIME-TIME1(IPROG)
         ENDIF
      ELSEIF (IOPT.EQ.3) THEN
         TIMTOT=TIME-TIME0

C Sum TIMTOT across all Nodes and ensure that all Nodes have same value here
#ifdef MPI
         buffer1(1)=timtot
         call MPI_AllReduce(buffer1(1),timtot,1,MPI_double_precision,
     .     MPI_sum,MPI_Comm_World,MPIerror)
#endif

         IF (TIMTOT.LT.TIMMIN) RETURN

         IF (PROG.EQ.'ALL' .OR. PROG.EQ.'all') THEN
           if (Node.eq.0) then
             WRITE (6,'(/,A)') 'timer: CPU execution times:'
             WRITE (6,'(A,2X,A10,A9,2A12,A9)') 'timer:',
     .         'Routine   ', 'Calls', 'Time/call', 'Tot.time', '%'
           endif
           DO 40 IPROG=1,NPROGS
             TIMETL=TIMET(IPROG)
             AVGTME=TIMET(IPROG)/NCALLS(IPROG)
             FRACTN=TIMET(IPROG)/TIMTOT*HUNDRD

C Sum values across Nodes
#ifdef MPI
             buffer1(1)=TIMETL
             buffer1(2)=AVGTME
             buffer1(3)=FRACTN
             call MPI_Reduce(buffer1,buffer2,3,MPI_double_precision,
     .         MPI_sum,0,MPI_Comm_World,MPIerror)
             TIMETL=buffer2(1)
             AVGTME=buffer2(2)
             FRACTN=buffer2(3)
#endif

             if (Node.eq.0) then
               WRITE(6,'(A,2X,A10,I9,2F12.3,F9.2)') 'timer:',
     .           PROGS(IPROG),NCALLS(IPROG),AVGTME,TIMETL,FRACTN
             endif
   40      CONTINUE
           if (Node.eq.0) then
             WRITE(6,*) ' '
           endif
         ELSE
           DO 50 IPROG=1,NPROGS
             IF (PROGS(IPROG).NE.PROG) GOTO 50
             TIMETL=TIMET(IPROG)
             FRACTN=TIMET(IPROG)/TIMTOT*HUNDRD

C Sum values across Nodes
#ifdef MPI
             buffer1(1)=TIMETL
             buffer1(2)=FRACTN
             call MPI_Reduce(buffer1,buffer2,2,MPI_double_precision,
     .         MPI_sum,0,MPI_Comm_World,MPIerror)
             TIMETL=buffer2(1)
             FRACTN=buffer2(2)
#endif

             if (Node.eq.0) then
               WRITE(6,'(A,A10,I6,F12.3,F7.2)')
     .          'timer: Routine,Calls,Time,% = ',
     .           PROGS(IPROG),NCALLS(IPROG),TIMETL,FRACTN
             endif
   50      CONTINUE
         ENDIF
      ELSE
         if (Node.eq.0) then
           WRITE(6,*) 'timer: INVALID OPTION IOPT =',IOPT
         endif
      ENDIF
      END
C
c-----------------------------------------------------------------
C
      SUBROUTINE elap1 (PROG,IOPT)


C  FINDS AND PRINTS THE WALL CLOCK SPENT IN DIFFERENT ROUTINES AND/OR
C   PROGRAM SECTIONS. IT MUST BE CALLED WITH IOPT=0 AT THE BEGINNING
C   OF EACH ROUTINE AND WITH IOPT=1 AT THE END OF IT.
C  ARGUMENTS:
C    PROG: INPUT ARBITRARY NAME FOR THE ROUTINE AND/OR PROGRAM SECTION
C    IOPT: INPUT OPTION PARAMETER:
C      IOPT = 0  => SET ALL TIMES TO ZERO (OPTIONAL)
C      IOPT = 1  => BEGIN COUNTING TIME FOR A ROUTINE
C      IOPT = 2  => STOP  COUNTING TIME FOR A ROUTINE
C      IOPT = 3  => PRINT TIME FOR A ROUTINE OR FOR ALL (IF PROG='ALL')
C  ROUTINE TIMES INCLUDE THAT SPENT IN THE ROUTINES THEY CALL
C  WRITTEN BY Alberto Garcia, Feb 2000, stealing code from
C  J.SOLER (JSOLER AT EMDUAM11) DEC/90

C
C  Modules
C
      use precision
      use parallel, only : Node
#ifdef MPI
      use mpi_siesta
#endif

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      integer :: nmax, iopt, nprogs, iprog
      PARAMETER (NMAX=500,ZERO=0.D0,HUNDRD=100.D0,TIMMIN=1.D-6)
      DIMENSION TIME1(NMAX),TIMET(NMAX)
      integer time_int, count_rate_int, wt
      double precision count_rate
      logical first
#ifdef MPI
      dimension buffer1(3),buffer2(3)
      integer mpierror
#endif
      INTEGER NCALLS(NMAX)
      CHARACTER*10 PROGS(NMAX),PROG*(*)
      SAVE PROGS,NPROGS,TIME0,TIME1,TIMET,NCALLS,count_rate,wt
      DATA NPROGS,TIME0 /0,ZERO/
      data first /.true./

      if (first) then
         CALL system_clock (count_rate=count_rate_int)
         count_rate = dble(count_rate_int)
         if (node.eq.0) then
          call io_assign(wt)
          open(wt,file='WALLTIME',form='formatted',status='unknown')
         endif
         first = .false.
      endif

      CALL system_clock (TIME_int)
      time = time_int / count_rate

      IF (IOPT.EQ.0) THEN
         NPROGS=0
         TIME0=TIME
      ELSEIF (IOPT.EQ.1 .OR. IOPT.EQ.2) THEN
         DO 10 IPROG=1,NPROGS
            IF (PROGS(IPROG).EQ.PROG) GO TO 20
   10    CONTINUE
            NPROGS=NPROGS+1
            IF (NPROGS.GT.NMAX) THEN
               if (Node.eq.0) then
                 WRITE (wt,*) 'elap: NMAX IS SATURATED. PROG = ',PROG
               endif
               RETURN
            ENDIF
            IPROG=NPROGS
            PROGS(IPROG)=PROG
            NCALLS(IPROG)=0
            TIMET(IPROG)=ZERO
   20    CONTINUE
         IF (IOPT.EQ.1) THEN
            NCALLS(IPROG)=NCALLS(IPROG)+1
            TIME1(IPROG)=TIME
         ELSE
            TIMET(IPROG)=TIMET(IPROG)+TIME-TIME1(IPROG)
         ENDIF
      ELSEIF (IOPT.EQ.3) THEN
         TIMTOT=TIME-TIME0

C Sum TIMTOT across all Nodes and ensure that all Nodes have same value here
#ifdef MPI
         buffer1(1)=timtot
         call MPI_AllReduce(buffer1(1),timtot,1,MPI_double_precision,
     .     MPI_sum,MPI_Comm_World,MPIerror)
#endif

         IF (TIMTOT.LT.TIMMIN) RETURN

         IF (PROG.EQ.'ALL' .OR. PROG.EQ.'all') THEN
           if (Node.eq.0) then
             WRITE (wt,'(/,A)') 'el: ELAPSED times:'
             WRITE (wt,'(A,2X,A10,A9,2A12,A9)') 'el:',
     .         'Routine   ', 'Calls', 'Time/call', 'Tot.time', '%'
           endif
           DO 40 IPROG=1,NPROGS
             TIMETL=TIMET(IPROG)
             AVGTME=TIMET(IPROG)/NCALLS(IPROG)
             FRACTN=TIMET(IPROG)/TIMTOT*HUNDRD

C Sum values across Nodes
#ifdef MPI
             buffer1(1)=TIMETL
             buffer1(2)=AVGTME
             buffer1(3)=FRACTN
             call MPI_Reduce(buffer1,buffer2,3,MPI_double_precision,
     .         MPI_sum,0,MPI_Comm_World,MPIerror)
             TIMETL=buffer2(1)
             AVGTME=buffer2(2)
             FRACTN=buffer2(3)
#endif

             if (Node.eq.0) then
               WRITE(wt,'(A,2X,A10,I9,2F12.3,F9.2)') 'elap:',
     .           PROGS(IPROG),NCALLS(IPROG),AVGTME,TIMETL,FRACTN
             endif
   40      CONTINUE
           if (Node.eq.0) then
             WRITE(wt,*) ' '
           endif
         ELSE
           DO 50 IPROG=1,NPROGS
             IF (PROGS(IPROG).NE.PROG) GOTO 50
             TIMETL=TIMET(IPROG)
             FRACTN=TIMET(IPROG)/TIMTOT*HUNDRD

C Sum values across Nodes
#ifdef MPI
             buffer1(1)=TIMETL
             buffer1(2)=FRACTN
             call MPI_Reduce(buffer1,buffer2,2,MPI_double_precision,
     .         MPI_sum,0,MPI_Comm_World,MPIerror)
             TIMETL=buffer2(1)
             FRACTN=buffer2(2)
#endif

             if (Node.eq.0) then
               WRITE(wt,'(A,A10,I6,F12.3,F7.2)')
     .          'el: Routine,Calls,Time Elapsed,% = ',
     .           PROGS(IPROG),NCALLS(IPROG),TIMETL,FRACTN
             endif
   50      CONTINUE
         ENDIF
      ELSE
         if (Node.eq.0) then
           WRITE(wt,*) 'elap: INVALID OPTION IOPT =',IOPT
         endif
      ENDIF
      END


