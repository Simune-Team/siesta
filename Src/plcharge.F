      SUBROUTINE PLCHARGE( MAXO, MAXA, MAXUO, MAXNH, MAXNA, NSPIN, 
     .                     ISA, IPHORB, INDXUO, LASTO, 
     .                     CELL, NSC, XA, RMAXO, DATM )

C **********************************************************************
C Prepare the data files to plot charge density at the points of a plane 
C in real space.
C The information is to be read by the external DENCHAR
C program, to plot charge density contours in 2D
C
C Coded by J. Junquera 11/98
C Modified by DSP, July 1999
C Modified by J. Junquera 7/01
C Modified by J. Junquera 2/02
C **********************************************************************

      use precision
      use old_atmfuncs
      use fdf
      use parallel,      only : Node, Nodes
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

      INTEGER, INTENT(IN) ::
     .  MAXO, MAXA, MAXUO, MAXNH, MAXNA, NSPIN

      INTEGER, INTENT(IN) ::
     .  LASTO(0:MAXA), ISA(MAXA), IPHORB(MAXO), INDXUO(MAXO), NSC(3)

      real(dp), INTENT(IN) ::
     .  CELL(3,3), XA(3,MAXA), RMAXO, DATM(MAXO)

C ****** INPUT *********************************************************
C INTEGER MAXO           : Maximum number of atomic orbitals in supercell
C INTEGER MAXA           : Maximum number of atoms in supercell
C INTEGER MAXUO          : Maximum number of atomic orbitals in unit cell.
C INTEGER MAXNH          : Maximum number
C                          of basis orbitals interacting, either directly
C                          or through a KB projector, with any orbital
C INTEGER MAXNA          : Maximum numbers of neighbour for any atom 
C INTEGER NSPIN          : Number of different spin polarizations
C                          Nspin = 1 => Non polarized. Nspin = 2 => Polarized
C INTEGER LASTO(0:MAXA)  : Last orbital of each atom in array iphorb
C INTEGER ISA(MAXA)      : Species index of each atom in the supercell
C INTEGER IPHORB(MAXO)   : Orbital index (within atom) of each orbital
C INTEGER INDXUO(MAXO)   : Equivalent orbital in unit cell
C INTEGER NSC(3)         : Num. of unit cells in each supercell direction
C REAL*8  CELL(3,3)      : Supercell vectors CELL(IXYZ,IVECT)
C                          (in bohrs)
C REAL*8  XA(3,MAXA)     : Atomic positions in cartesian coordinates
C                          (in bohrs)
C REAL*8  RMAXO          : Maximum range of basis orbitals
C REAL*8  DATM(MAXO)     : Occupations of basis orbitals in free atom
C **********************************************************************

C Internal variables ---------------------------------------------------

      CHARACTER*33 PASTE

      CHARACTER*30
     .  SNAME, FNAME1, FNAME2

      INTEGER
     .  UNIT1, UNIT2, IL, IA, J, MAXNHG

#ifdef MPI
      INTEGER MPIerror
#endif

      EXTERNAL
     .  IO_ASSIGN, IO_CLOSE, PASTE

#ifdef MPI
      call MPI_AllReduce(MAXNH,MAXNHG,1,MPI_integer,MPI_sum,
     .                   MPI_Comm_World,MPIerror)
#endif

C Write only if we are in the node 0 -----------------------------------
      IF (Node .EQ. 0) THEN
C Assign the name of the output file -----------------------------------
        SNAME = FDF_STRING('SystemLabel','siesta')
        FNAME1 = PASTE(sname,'.DIM')
        FNAME2 = PASTE(sname,'.PLD')

        CALL IO_ASSIGN(UNIT1)
          OPEN ( UNIT = UNIT1, FILE = FNAME1, FORM = 'UNFORMATTED',
     .           STATUS = 'UNKNOWN' )

            WRITE(UNIT1)MAXA
            WRITE(UNIT1)MAXO
            WRITE(UNIT1)MAXUO 
            WRITE(UNIT1)NSPIN
#ifdef MPI
            WRITE(UNIT1)MAXNHG
#else
            WRITE(UNIT1)MAXNH
#endif
            WRITE(UNIT1)MAXNA

        CALL IO_CLOSE(UNIT1)

        CALL IO_ASSIGN(UNIT2)

        OPEN ( UNIT = UNIT2, FILE = FNAME2, FORM = 'UNFORMATTED',
     .         STATUS = 'UNKNOWN' )
C Dump the tables into a file ------------------------------------------

        WRITE(UNIT2) RMAXO
C       WRITE(6,*) RMAXO

        DO IL = 1, MAXO
          WRITE(UNIT2)IPHORB(IL), INDXUO(IL), DATM(IL)
        ENDDO

        DO IA = 1, MAXA
          WRITE(UNIT2)ISA(IA)
        ENDDO

        DO IA = 0, MAXA
          WRITE(UNIT2)LASTO(IA)
        ENDDO
        
        DO IA = 1,3
          WRITE(UNIT2)(CELL(J,IA),J=1,3)
        ENDDO

        WRITE(UNIT2)(NSC(IA),IA=1,3)

        DO IA = 1, MAXA
          WRITE(UNIT2)(XA(J,IA),J=1,3)
        ENDDO

        CALL IO_CLOSE(UNIT2)
 
      ENDIF

      RETURN

      END
