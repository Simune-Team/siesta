      SUBROUTINE PLCHARGE( NSPIN, NUOTOT, NO, NA, MAXO, MAXA, MAXNH, 
     .                     CELL, RMAXO, XA, 
     .                     LASTO, ISA, IPHORB, DATM,
     .                     NUMH, LISTHPTR, LISTH, INDXUO, Nodes)
C **********************************************************************
C Prepare the data files to plot charge density at the points of a plane 
C in real space.
C The information is to be read by the external DENCHAR
C program, to plot charge density contours in 2D
C
C Coded by J. Junquera 11/98
C Modified by DSP, July 1999
C **********************************************************************

      use old_atmfuncs
      use fdf
      use parallel
#ifdef MPI
      use mpi
#endif

      IMPLICIT NONE


      CHARACTER*33 PASTE

      CHARACTER*30
     .  SNAME,FNAME

      INTEGER
     .  UNIT1, MAXO, MAXA, MAXNH, NUOTOT, NO, NA, NSPIN, 
     .  LASTO(0:MAXA), ISA(MAXA), IPHORB(MAXO),
     .  IL, IS, IB, IA, IIO, IO, J, NUMH(*), LISTH(MAXNH),
     .  LISTHPTR(*), INDXUO(MAXO), Node, Nodes

      DOUBLE PRECISION
     .  CELL(3,3), RMAXO, XA(3,MAXA), DATM(MAXO)


#ifdef MPI
      integer MPIerror
#endif

      EXTERNAL
     .  IO_ASSIGN, IO_CLOSE, PASTE


C **********************************************************************
C INTEGER MAXO           : Max. total number of basis orbitals
C INTEGER MAXA           : Max. number of atoms
C INTEGER MAXNH          : Dimension of listh
C INTEGER NSPIN          : Number of different spin polarizations
C                          Nspin = 1 => Non polarized. Nspin = 2 => Polarized
C INTEGER NUOTOT         : Number of atomic orbitals in unit cell
C INTEGER NO             : Total number of orbitals in supercell
C INTEGER NA             : Number of atoms in supercell
C REAL*8  CELL(3,3)      : Supercell vectors CELL(IXYZ,IVECT)
C REAL*8  XA(3,MAXA)     : Atomic positions in cartesian coordinates
C REAL*8  RMAXO          : Maximum cutoff of atomic orbitals
C INTEGER LASTO(0:MAXA)  : Position of last orbital of each atom
C INTEGER ISA(MAXA)      : Species index of each atom
C INTEGER IPHORB(MAXO)   : Orbital index of each atom in its atom
C REAL*8  DATM(MAXO)     : Neutral atom charge of each orbital
C INTEGER NUMH(NUOTOT)   : Control vector of Density Matrix
C                          (number of nonzero elements of each row)
C INTEGER LISTHPTR(NUOTOT) : Control vector of Density Matrix
C INTEGER LISTH(MAXNH)   : Control vector of Density Matrix
C                          (list of nonzero elements of each row)
C INTEGER INDXUO(MAXO)   : Index of equivalent orbital in unit cell
C **********************************************************************

#ifdef MPI
      call MPI_Comm_Rank( MPI_Comm_World, Node, MPIerror )
      call MPI_Comm_Size( MPI_Comm_World, Nodes, MPIerror )
#else
      Node =  0
      Nodes = 1
#endif

C Assign the name of the output file -----------------------------------
      SNAME = FDF_STRING('SystemLabel','siesta')
      FNAME = PASTE(sname,'.PLD')

      CALL IO_ASSIGN(UNIT1)
C     REWIND(UNIT1)

      OPEN ( UNIT = UNIT1, FILE = FNAME, FORM = 'UNFORMATTED',
     .       STATUS = 'UNKNOWN' )
C Dump the tables into a file ------------------------------------------

CIn the previous version this loop only was executed from
CIL= -LMAXD+1 instead of -(LMAXD+1).
C I belive that now is right! DSP.
      DO IB = 1, NTBMAX+2
        DO IL = -(LMAXD+1), NSEMX*NZETMX*(LMAXD+1)
          DO IS = 1, NSMAX
            WRITE(UNIT1)TABLE(IB,IL,IS)
          ENDDO
        ENDDO
      ENDDO
 
      DO IB = 1, NTBMAX+2
        DO IL = 1, NSEMX*NZETMX*(LMAXD+1)
          DO IS = 1, NSMAX
            WRITE(UNIT1)TABPOL(IB,IL,IS)
          ENDDO
        ENDDO
      ENDDO
 
CIn the previous version this loop only was executed from
CIL= -LMAXD+1 instead of -(LMAXD+1).
C I belive that now is right! DSP.
      DO IB = 1, NTBMAX
        DO IL = -(LMAXD+1), NSEMX*NZETMX*(LMAXD+1)
          DO IS = 1,NSMAX
            WRITE(UNIT1)TAB2(IB,IL,IS)
          ENDDO
        ENDDO
      ENDDO
 

      DO IB = 1, NTBMAX
        DO IL = 1, NSEMX*NZETMX*(LMAXD+1)
          DO IS = 1,NSMAX
            WRITE(UNIT1)TAB2POL(IB,IL,IS)
          ENDDO
        ENDDO
      ENDDO

      DO IS = 1,NSMAX
        WRITE(UNIT1)LMXOSAVE(IS)
      ENDDO 

      DO IS = 1, NSMAX
        WRITE(UNIT1)LMXKBSAVE(IS)
      ENDDO 


      DO IL = 0,LMAXD
        DO IB= 1, NSEMX
          DO IS = 1,NSMAX
            WRITE(UNIT1)NZETASAVE(IL,IB,IS)
          ENDDO
        ENDDO 
      ENDDO
      
      DO IS = 1, NSMAX 
        DO IB = 1, NSEMX
          DO IL =0 , LMAXD
            WRITE(UNIT1)NPOLORBSAVE(IL,IB,IS)
          ENDDO 
        ENDDO 
      ENDDO 
      
      DO IS=1, NSMAX
        DO IL= 0, LMAXD
          WRITE(UNIT1)NSEMICSAVE(IL,IS)
        ENDDO 
      ENDDO

      DO IS=1, NSMAX
        DO IL= 0, LMAXD
          WRITE(UNIT1)NKBLSAVE(IL,IS)
        ENDDO 
      ENDDO

      WRITE(UNIT1)ISMAX
               
      DO IS = 1,NSMAX 
        WRITE(UNIT1)NOMAX(IS), NKBMAX(IS)
      ENDDO

      DO IL = 1, MAXO
        WRITE(UNIT1)DATM(IL)
      ENDDO

      WRITE(UNIT1)NUOTOT
      WRITE(UNIT1)NO
      WRITE(UNIT1)NA
      WRITE(UNIT1)NSPIN

      DO IA = 1,3
        WRITE(UNIT1)(CELL(J,IA),J=1,3)
      ENDDO

      DO IA = 1,MAXA
        WRITE(UNIT1)(XA(J,IA),J=1,3),ISA(IA)
      ENDDO

      DO IA = 0,MAXA
        WRITE(UNIT1)LASTO(IA)
      ENDDO

      WRITE(UNIT1)RMAXO

      DO IO = 1,MAXO
        WRITE(UNIT1)IPHORB(IO),INDXUO(IO)
      ENDDO

      DO IO = 1,NUOTOT
        call GlobalToLocalOrb(IO,Node,Nodes,IIO)
        if (IIO.gt.0) then
          WRITE(UNIT1)NUMH(IIO)
          DO J = 1,NUMH(IIO)
            WRITE(UNIT1)LISTH(LISTHPTR(IIO)+J)
          ENDDO
        endif
#ifdef MPI
        call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif
      ENDDO

      CALL IO_CLOSE(UNIT1)

      RETURN

      END
