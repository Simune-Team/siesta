C $Id: prmem.f,v 1.6 1999/11/26 18:28:27 wdpgaara Exp $

      SUBROUTINE PRMEM( IUNIT, PROG, VAR, TYPE, DIM )
C ********************************************************************
C Prints memory used by arrays
C Written by J.M.Soler. June 1997.
C ********** Input **************************************************
C INTEGER   IUNIT    : I/O unit to be printed on
C CHARACTER PROG*(*) : Name of calling routine 
C CHARACTER VAR*(*)  : Name of array, or ' ' (one blank) for Total
C CHARACTER TYPE     : Array type : 'C' = CHARACTER
C                                   'I' = INTEGER
C                                   'R' = REAL
C                                   'D' = DOUBLE PRECISION
C                                   'L' = LOGICAL
C INTEGER   DIM      : Size of array
C ********** Behaviour ***********************************************
C If IUNIT=0, it first opens a file called PROG.size and then prints
C   sizes of subsequent calls (with IUNIT=0) into it. To print the 
C   total of all routines CALL PRMEM(0,' ',' ',' ',0)
C ********************************************************************


      IMPLICIT          NONE
      INTEGER           MAXPR, NC, NTY
      PARAMETER       ( MAXPR = 90, NC = 2, NTY = 7 )
      CHARACTER         PROG*(*), VAR*(*)
      CHARACTER         CTY(NC,NTY), TYPE, FNAME*14
      CHARACTER*10      PROGS(MAXPR)
C     CHARACTER*10      NAMETY(NTY)
      INTEGER           DIM, IC, IPR, ITY, IU, IUNIT, NPR, MEMTY(NTY),
     .                  MYUNIT
      DOUBLE PRECISION  MEMPR(MAXPR), MEMTOT, MBYTE
      SAVE              CTY, IU, MEMPR, MEMTOT, MYUNIT, MEMTY,
C    .                  NAMETY, 
     .                  NPR, PROGS

C     DATA ( NAMETY(ITY), (CTY(IC,ITY),IC=1,NC),
C    .       MEMTY(ITY), ITY=1,NTY) /
C    .  'Character', 'c', 'C', 1,
C    .  'Short',     's', 'S', 2,
C    .  'Integer',   'i', 'I', 4,
C    .  'Logical',   'l', 'L', 4,
C    .  'Real',      'r', 'R', 4,
C    .  'Float',     'f', 'F', 4,
C    .  'Double',    'd', 'D', 8 /
      DATA ( (CTY(IC,ITY),IC=1,NC),
     .       MEMTY(ITY), ITY=1,NTY) /
     .  'c', 'C', 1,
     .  's', 'S', 2,
     .  'i', 'I', 4,
     .  'l', 'L', 4,
     .  'r', 'R', 4,
     .  'f', 'F', 4,
     .  'd', 'D', 8 /
      DATA MBYTE / 1.D6 /
      DATA NPR / 0 /
      DATA MYUNIT / 0 /

      DO 10 IPR = 1,NPR
        IF ( PROG .EQ. PROGS(IPR) ) GOTO 20
   10 CONTINUE
        NPR = NPR + 1
        CALL CHKDIM( 'PRMEM', 'MAXPR', MAXPR, NPR, 1 )
        PROGS(NPR) = PROG
        MEMPR(NPR) = 0.D0
        IPR = NPR
        IF (IUNIT .LE. 0) THEN
          IF (MYUNIT .EQ. 0) THEN
            CALL IO_ASSIGN(MYUNIT)
            FNAME = PROG//'.size'
            OPEN( UNIT=MYUNIT, FILE=FNAME, STATUS='UNKNOWN' )
            WRITE(MYUNIT,'(/,2A)')
     .        'C PRMEM: Memory for internal arrays in ', PROG
            IF (VAR.EQ.' ') RETURN
          ENDIF
          IU = MYUNIT
        ELSE
          IU = IUNIT
          WRITE(IU,'(/,2A)')
     .      'C PRMEM: Memory for internal arrays in ', PROG
        ENDIF
   20 CONTINUE

      IF (PROG.EQ.' ') THEN
        WRITE(IU,99) 'Total', 'Total', MEMTOT / MBYTE
        WRITE(IU,*) ' '
      ELSEIF (VAR.EQ.' ') THEN
        WRITE(IU,99) PROG, 'Total', MEMPR(IPR) / MBYTE
        WRITE(IU,*) ' '
        PROGS(IPR) = PROGS(NPR)
        MEMPR(IPR) = MEMPR(NPR)
        NPR = NPR - 1
      ELSE
        DO 40 ITY = 1,NTY
          DO 30 IC = 1,NC
            IF ( CTY(IC,ITY) .EQ. TYPE ) THEN
              WRITE(IU,99) PROG, VAR, MEMTY(ITY) * DIM / MBYTE
              MEMPR(IPR) = MEMPR(IPR) + MEMTY(ITY) * DIM
              MEMTOT = MEMTOT + MEMTY(ITY) * DIM
              RETURN
            ENDIF
   30     CONTINUE
   40   CONTINUE
        WRITE(6,*) 'PRMEM: TYPE  ', TYPE, '  NOT FOUND'
      ENDIF

   99 FORMAT('C PRMEM:', 2A12, F12.3, ' MBytes')
      END

