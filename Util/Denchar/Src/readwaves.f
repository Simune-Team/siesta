      SUBROUTINE READWAVES(NSPIN,NORB,IFLAG,NWF,PSI,E,IND)

C Reads the wavefunctions and energies from a file written by Siesta
C
C P. Ordejon, July 2003
C **************** INPUT ********************************************
C INTEGER NSPIN     : Number of spin components
C INTEGER NORB      : Number of basis orbitals
C INTEGER IFLAG     : 0=only read and return number of wavefunctions
C                     1=actually read wavefunctions
C **************** INPUT OR OUTPUT **********************************
C INTEGER NWF       : Number of wavefunctions to read 
C                     input(output) if IFLAG=0(1)
C **************** OUTPUT *******************************************
C REAL*8 PSI(NORB,NWF,NSPIN): Wavefunctions
C REAL*8 E(NWF,NSPIN)       : Eigenvalues
C INTEGER IND(NWF)          : List of indexes of wavefunctions
C *******************************************************************

C Modules

      use fdf


      IMPLICIT NONE

      INTEGER IFLAG, NSPIN, NORB, NWF
      INTEGER IND(NWF)
      DOUBLE PRECISION PSI(NORB,NWF,NSPIN), E(NWF,NSPIN)

C INTERNAL VARIABLES .............
      INTEGER UNIT, NK, NSP, NUO, ISPIN, IISPIN, IWF, IIWF, IORB
      INTEGER IDUMB
      DOUBLE PRECISION REPSI,IMPSI

      CHARACTER PASTE*33
      CHARACTER, SAVE :: SNAME*30, FNAME*33
      CHARACTER CHDUMB*20

      SAVE UNIT
      EXTERNAL PASTE
C ..................


c      write(6,*) NORB,NWF,NSPIN

c      write(6,*) 'In readwaves with iflag=',iflag

      IF (IFLAG .EQ. 0) THEN
        SNAME = FDF_STRING('SystemLabel','siesta')
        FNAME = PASTE(SNAME,'.WFS')

        CALL IO_ASSIGN(UNIT)
        OPEN (UNIT, FILE=FNAME, FORM='unformatted', STATUS='unknown')

c        write(6,*) 'opening unit=', unit

        READ(UNIT) NK
        IF (NK .NE. 1) THEN
          WRITE(6,*) 'Wavefunctions file contains more then 1 k-point'
          WRITE(6,*) 'DENCHAR can only handle the Gamma point!!'
          STOP
        ENDIF
        READ(UNIT) NSP
        IF (NSP .NE. NSPIN) THEN
          WRITE(6,*) 'NSPIN is not consistent between data files!'
          STOP
        ENDIF
        READ(UNIT) NUO
        IF (NUO .NE. NORB) THEN
          WRITE(6,*) 'Nr. of orbs is not consistent between data files!'
          STOP
        ENDIF

        READ(UNIT) 
        READ(UNIT) IISPIN
        READ(UNIT)NWF

        REWIND(UNIT)
        
c        write(6,*) 'Exiting iflag= ',iflag
         
        
        RETURN

      ELSE IF (IFLAG .EQ. 1) THEN

c        write(6,*) 'trying to read in unit=', unit
        READ(UNIT) NK
        IF (NK .NE. 1) THEN
          WRITE(6,*) 'Wavefunctions file contains more then 1 k-point'
          WRITE(6,*) 'DENCHAR can only handle the Gamma point!!'
          STOP
        ENDIF
        READ(UNIT) NSP
        IF (NSP .NE. NSPIN) THEN
          WRITE(6,*) 'NSPIN is not consistent between data files!'
          STOP
        ENDIF
        READ(UNIT) NUO
        IF (NUO .NE. NORB) THEN
          WRITE(6,*) 'Nr. of orbs is not consistent between data files!'
          STOP
        ENDIF

        DO ISPIN = 1,NSPIN
          READ(UNIT) 
          READ(UNIT) IISPIN
c          write(6,*) 'spin=',iispin
          IF (IISPIN .NE. ISPIN) THEN
            WRITE(6,*) 'Inconsistent order of spins in wavefuncs. file!'
            STOP
          ENDIF
          READ(UNIT)NWF
c          write(6,*) 'nwfs=',nwf
c          write(6,*) 'number of orbitals=',norb
 

          DO IWF=1,NWF
            READ(UNIT) IND(IWF)
c            write(6,*) 'orbital index ',iwf,' = ',ind(iwf)
            READ(UNIT) E(IWF,ISPIN)
c            write(6,*) 'energy = ',e(iwf,ispin)
            DO IORB = 1, NORB
              READ(UNIT) IDUMB,CHDUMB,IDUMB,IDUMB,IDUMB,CHDUMB,
     .                   REPSI,IMPSI
c              write(6,*) repsi,impsi
              IF (DABS(IMPSI) .GT. 1.0D-10) 
     .          WRITE(6,*) 'Warning: complex wavefunctions in file!'
              PSI(IORB,IWF,ISPIN)=REPSI
            ENDDO
          ENDDO
        ENDDO

        CLOSE (UNIT)
        CALL IO_CLOSE(UNIT)
c        write(6,*) 'Exiting iflag= ',iflag
      ELSE
        WRITE(6,*) 'IFLAG must be either 0 or 1 in READWAVE!!'
        STOP
      ENDIF
           

      RETURN

      END
        

