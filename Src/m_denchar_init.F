      module m_denchar_init

      integer, dimension(3), save  :: num_version = (/1,1,0/)
      character(len=80), parameter :: version_str =
     . "DENCHAR 1.1.0 -- [Translation into F90] (Feb 04, 2002)"

      external :: io_assign, io_close

      CONTAINS

      SUBROUTINE REINIT( MAXO, MAXA, MAXUO, MAXNH, MAXNA, NSPIN, IDIMEN,
     .                   CHARGE, WAVES )

C **********************************************************************
C Read some variables from SIESTA in order to define
C the dimensions of some arrays in DENCHAR 
C
C Coded by J. Junquera 07/01
C Modifications for 3D and wavefunctions by P. Ordejon, June 2003
C **********************************************************************

      USE FDF

      IMPLICIT NONE

      INTEGER, INTENT(OUT) ::
     .  MAXO, MAXA, MAXUO, NSPIN, MAXNH, MAXNA, IDIMEN
      LOGICAL, INTENT(OUT) ::
     .  CHARGE, WAVES

C **********************************************************************
C INTEGER MAXO           : Maximum number of atomic orbitals in supercell
C INTEGER MAXA           : Maximum number of atoms in supercell
C INTEGER MAXUO          : Number of atomic orbitals in unit cell.
C INTEGER MAXNH          : Maximum number
C                          of basis orbitals interacting, either directly
C                          or through a KB projector, with any orbital
C INTEGER MAXNA          : Maximum number of neighbour of any atom
C INTEGER NSPIN          : Number of different spin polarizations
C                          Nspin = 1 => Non polarized. Nspin = 2 => Polarized
C INTEGER IDIMEN         : Type of calculation: 2D or 3D
C LOGICAL CHARGE         : Should charge density be computed?
C LOGICAL WAVES          : Should wave functions be computed?
C **********************************************************************

C Internal variables --------------------------------------------------

      CHARACTER*33 PASTE

      CHARACTER*30
     .  SNAME, FNAME1

      CHARACTER 
     .  TYPRUN*2, TYPRUN_DEFECT*2


      INTEGER
     .  UNIT1 

      LOGICAL
     .  LEQI

      EXTERNAL PASTE, LEQI

C Assign the name of the output file -----------------------------------
      SNAME = FDF_STRING('SystemLabel','siesta')
      FNAME1 = PASTE(sname,'.DIM')

      TYPRUN_DEFECT = '2D'
      TYPRUN = FDF_STRING('Denchar.TypeOfRun',TYPRUN_DEFECT)
      IF (LEQI(TYPRUN,'2D')) THEN
        IDIMEN  = 2
      ELSE IF (LEQI(TYPRUN,'3D')) THEN
        IDIMEN  = 3
      ELSE
        WRITE(6,'(A)')' readpla:  Wrong type of run; must be 2D or 3D  '
        STOP
      ENDIF

      CHARGE = FDF_BOOLEAN('Denchar.PlotCharge',.FALSE.)
      WAVES  = FDF_BOOLEAN('Denchar.PlotWaveFunctions',.FALSE.)


      IF (.NOT. CHARGE .AND. .NOT. WAVES) THEN
        WRITE(6,*)'Denchar.PlotCharge and Denchar.PlotWaveFunctions'
        WRITE(6,*)'are both .FALSE.'
        WRITE(6,*)'At least one of them should be .TRUE.'
        STOP
      ENDIF

      CALL IO_ASSIGN(UNIT1)
        OPEN ( UNIT = UNIT1, FILE = FNAME1, FORM = 'UNFORMATTED',
     .         STATUS = 'UNKNOWN' )

          READ(UNIT1)MAXA
          READ(UNIT1)MAXO
          READ(UNIT1)MAXUO 
          READ(UNIT1)NSPIN
          READ(UNIT1)MAXNH
          READ(UNIT1)MAXNA

      CALL IO_CLOSE(UNIT1)

      RETURN

      END subroutine reinit

!
      SUBROUTINE REDATA( MAXO, MAXA, 
     .                   ISA, IPHORB, INDXUO, LASTO,
     .                   CELL, NSC, XA, RMAXO, DATM )

C **********************************************************************
C Read the data files to plot charge density at the points of a plane 
C or at a 3D grid in real space.
C The information is written by the subroutine plcharge in SIESTA,
C when WriteDenchar option is set up to .true. in the fdf input data 
C file.
C
C Coded by J. Junquera 11/98
C Modified by DSP, July 1999
C Modified by J. Junquera, 7/01
C Modified by P. Ordejon to include 3D capabilities, June 2003
C **********************************************************************

      USE FDF

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::  MAXO, MAXA

      INTEGER, INTENT(OUT) ::
     .  LASTO(0:MAXA), ISA(MAXA), IPHORB(MAXO), INDXUO(MAXO), NSC(3)

      DOUBLE PRECISION, INTENT(OUT) ::
     .  CELL(3,3), XA(3,MAXA), RMAXO, DATM(MAXO)


C **** INPUT ***********************************************************
C INTEGER MAXO           : Maximum number of atomic orbitals in supercell
C INTEGER MAXA           : Maximum number of atoms in supercell
C **** OUTPUT **********************************************************
C INTEGER ISA(MAXA)      : Species index of each atom in the supercell
C INTEGER IPHORB(MAXO)   : Orbital index (within atom) of each orbital
C INTEGER INDXUO(MAXO)   : Equivalent orbital in unit cell
C INTEGER LASTO(0:MAXA)  : Last orbital of each atom in array iphorb
C REAL*8  CELL(3,3)      : Supercell vectors CELL(IXYZ,IVECT)
C                          (in bohrs)
C INTEGER NSC(3)         : Num. of unit cells in each supercell direction
C REAL*8  XA(3,MAXA)     : Atomic positions in cartesian coordinates
C                          (in bohrs)
C REAL*8  RMAXO          : Maximum range of basis orbitals
C REAL*8  DATM(MAXO)     : Occupations of basis orbitals in free atom
C **********************************************************************

C Internal variables ---------------------------------------------------

      CHARACTER*33 PASTE

      CHARACTER*30
     .  SNAME, FNAME

      INTEGER
     .  UNIT1, IL, IA, J

      EXTERNAL  PASTE


C Assign the name of the output file -----------------------------------
      SNAME = FDF_STRING('SystemLabel','siesta')
      FNAME = PASTE(sname,'.PLD')

      CALL IO_ASSIGN(UNIT1)

      OPEN ( UNIT = UNIT1, FILE = FNAME, FORM = 'UNFORMATTED',
     .       STATUS = 'UNKNOWN' )

        READ(UNIT1)RMAXO

        DO IL = 1, MAXO
          READ(UNIT1)IPHORB(IL), INDXUO(IL), DATM(IL)
        ENDDO

        DO IA = 1, MAXA
          READ(UNIT1)ISA(IA)
        ENDDO

        DO IA = 0, MAXA
          READ(UNIT1)LASTO(IA)
        ENDDO

        DO IA = 1,3
          READ(UNIT1)(CELL(J,IA),J=1,3)
        ENDDO

        READ(UNIT1)(NSC(IA),IA=1,3)

        DO IA = 1, MAXA
          READ(UNIT1)(XA(J,IA),J=1,3)
        ENDDO

      CALL IO_CLOSE(UNIT1)

      END SUBROUTINE REDATA

!---------------------------------------------------------------------
      subroutine prversion
c
c     Simple routine to print the version string. Could be extended to
c     provide more information, if needed.
c
      write(6,'(a)') trim(version_str)
#ifdef MPI
      write(6,'(a)') 'PARALLEL version'
#else
      write(6,'(a)') 'SERIAL version'
#endif

#ifdef CDF
      write(6,'(a)') 'NetCDF-capable'
#endif

      end subroutine prversion
!------------------------------------------------------------------------
      subroutine get_version(v)
      integer, intent(out)  :: v(3)
      v = num_version
      end subroutine get_version

      END MODULE m_denchar_init








