      SUBROUTINE PDOS( NO, NSPIN, MAXSPN, MAXUO, MAXNH, 
     .                 MAXO, NUMH, LISTHPTR, LISTH, H, S, 
     .                 E1, E2, SIGMA, NHIST,
     .                 GAMMA, XIJ, INDXUO, NK, KPOINT, WK, EO,
     .                 NBK, BK, EBK, NUOTOT )
C **********************************************************************
C Subroutine to calculate the proyected density of states on the
C atomic orbitals for a given eigenvalue spectra
C Written by J. Junquera and E. Artacho, November 1999.
C ***********  INPUT  **************************************************
C INTEGER NO                  : Number of basis orbitals in the supercell
C INTEGER NSPIN               : Spin polarizations (1 or 2)
C INTEGER MAXSPN              : Second dimension of eo and qo
C                               (Max number of different spin polarizations)
C INTEGER MAXUO               : Maximum number of atomic orbitals in the unit 
C                               cell. First dimension of eo, qo, last of xij
C                               Must be at least max(indxuo)
C INTEGER MAXNH               : Maximum number of orbitals interacting 
C                               with any orbital
C INTEGER MAXO                : First dimension of eo 
C INTEGER NUMH(NUO)           : Number of nonzero elements of each row
C                               of hamiltonian matrix
C INTEGER LISTH(MAXNH)        : Nonzero hamiltonian-matrix element
C                               column indexes for each matrix row
C INTEGER LISTHPTR(NUO)       : Pointer to each row (-1) of the
C                               density matrix
C REAL*8  H(MAXNH,NSPIN)      : Hamiltonian in sparse format
C REAL*8  S(MAXNH)            : Overlap in sparse format
C REAL*8  E1, E2              : Energy range for density-matrix states
C                               (to find local density of states)
C                               Not used if e1 > e2
C REAL*8  SIGMA               : Width of the gaussian to expand the eigenvalues
C INTEGER NHIST               : Number of subdivisions of the histogram
C LOGICAL GAMMA               : Only gamma point?
C REAL*8  XIJ(3,MAXNH)        : Vectors between orbital centers (sparse)
C                               (not used if only gamma point)
C INTEGER INDXUO(NO)          : Index of equivalent orbital in unit cell
C                               Unit cell orbitals must be the first in
C                               orbital lists, i.e. indxuo.le.nuo, with
C                               nuo the nuber of orbitals in the unit cell
C INTEGER NK                  : Number of k points
C REAL*8  KPOINT(3,NK)        : k point vectors
C REAL*8  WK(NK)              : k point weights (must sum one)
C REAL*8  EO(MAXUO,MAXSPN,NK) : Eigenvalues
C INTEGER NBK                 : Number of bands k-points
C REAL*8  BK(3,NBK)           : k-point vectors to compute the band structure
C REAL*8  EBK(MAXUO,MAXSPN,NBK): Eigenvalues for the band structure
C INTEGER NUOTOT              : Total number of orbitals in unit cell
C **********************************************************************

      USE FDF
      use atmfuncs
      use atomlist, dummy => indxuo
      use xml

      IMPLICIT NONE

      INTEGER
     .  NO, NSPIN, MAXSPN, MAXUO, MAXNH, NK, NHIST, NBK, 
     .  MAXO, NUOTOT

      INTEGER 
     .  NUMH(*), LISTH(MAXNH), LISTHPTR(*), INDXUO(NO)

      DOUBLE PRECISION
     .  H(MAXNH,NSPIN), S(MAXNH), E1, E2, SIGMA,  
     .  XIJ(3,MAXNH), KPOINT(3,NK), WK(NK), BK(3,NBK),
     .  EO(MAXO,MAXSPN,NK), EBK(MAXO,MAXSPN,NBK)

      LOGICAL
     .  GAMMA

C Dynamic arrays -------------------------------------------------------
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE :: HAUX, SAUX
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE :: PSI, AUX
      INTEGER         , DIMENSION(:), ALLOCATABLE, SAVE :: MUO
      DOUBLE PRECISION, DIMENSION(:,:)  , ALLOCATABLE :: DTOT
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: DPR

C Internal variables ---------------------------------------------------
      INTEGER
     .  NUO, IO, NHS, NPSI, NAUX, IUO, IHIST, ISPIN, 
     .  NTY, ITY, IA, IS, IPHI, IN, IUNIT1, IUNIT2, I

      integer iat, spec, ii, iorb

      CHARACTER*7  IL
      character*40 pos_string

      CHARACTER
     .  SNAME*30, FNAMETOT*37, FNAMEPRO*37, FNAMEPAR*53, 
     .  PARTIAL1*11, PARTIAL2*18,
     .  PARTIAL3*31, PARTIAL4*22, PASTE*53, ICHARA*10, 
     .  INTCHAR*10, POINT

      DOUBLE PRECISION
     .  DELTA, ENER, EV

      EXTERNAL
     .  INTCHAR, IO_ASSIGN, IO_CLOSE, PASTE, 
     .  PDOSK, TIMER


      CALL TIMER( 'pdos', 1)

C Find the intervals between the subdivisions in the energy scale ------
      DELTA = (E2 - E1) / NHIST
      EV = 13.6058D0
      SNAME = FDF_STRING('SystemLabel','siesta')
      POINT = '.'

C Find number of orbitals per unit cell and check argument sizes -------
      NUO = 0
      DO IO = 1, NO
        NUO = MAX( NUO, INDXUO(IO) )
      ENDDO

C Check internal dimensions --------------------------------------------
      IF( NSPIN.LE.2 .AND. GAMMA) THEN
         NHS  = NUOTOT * NUO
         NPSI = NUOTOT * MAXUO * NSPIN
         NAUX = NUOTOT * 5
      ELSE IF( NSPIN.LE.2 .AND. .NOT.GAMMA) THEN
         NHS  = 2 * NUOTOT * NUO
         NPSI = 2 * NUOTOT * NUO
         NAUX = 2 * NUOTOT * 5    !!! Following AG in diagon.F
      ENDIF

C Allocate local arrays ------------------------------------------------
      ALLOCATE(HAUX(NHS))
      CALL MEMORY('A','D',NHS,'pdos')
      ALLOCATE(SAUX(NHS))
      CALL MEMORY('A','D',NHS,'pdos')
      ALLOCATE(PSI(NPSI))
      CALL MEMORY('A','D',NPSI,'pdos')
      ALLOCATE(AUX(NAUX))
      CALL MEMORY('A','D',NAUX,'pdos')
      ALLOCATE(MUO(NUO))
      CALL MEMORY('A','I',NUO,'pdos')
      ALLOCATE(DTOT(NHIST,2))
      CALL MEMORY('A','D',2*NHIST,'pdos')
      ALLOCATE(DPR(NHIST,NUO,2))
      CALL MEMORY('A','D',2*NHIST*NUO,'pdos')

C Check indxuo ---------------------------------------------------------
      DO IUO = 1, NUO
        MUO(IUO) = 0
      ENDDO

      DO IO = 1, NO
        IUO = INDXUO(IO)
        IF( IUO.LE.0 .OR. IUO.GT.NUOTOT ) THEN
          WRITE(6,*)'pdos: ERROR: invalid index: io, indxuo = ',
     .              IO, INDXUO(IO)
          STOP 'pdos: ERROR: invalid indxuo'
        ENDIF
        MUO(IUO) = MUO(IUO) + 1
      ENDDO

      DO IO = 1, NUO
        IF( MUO(IUO) . NE. MUO(1)) THEN
          WRITE(6,'(/,2a,3i6)')' pdos: ERROR: inconsistent indxuo.',
     .            ' IUO, MUO(IUO), MUO(1) =', IUO, MUO(IUO), MUO(1)
          STOP 'pdos: ERROR: inconsistent indxuo.'
        ENDIF
      ENDDO

C Initialize the projected density of states ---------------------------
      DO ISPIN = 1, 2
        DO IHIST = 1, NHIST
          DTOT(IHIST,ISPIN) = 0.D0
          DO IUO = 1, NUO
            DPR(IHIST,IUO,ISPIN) = 0.D0
          ENDDO
        ENDDO
      ENDDO

C Call appropiate routine ----------------------------------------------
      IF( NSPIN.LE.2 .AND. GAMMA) THEN
        WRITE(6,*)'Pdos only at Gamma point not yet implemented'
      ELSE IF( NSPIN.LE.2 .AND. .NOT.GAMMA) THEN
        CALL PDOSK( NSPIN, NUO, NO, MAXSPN, MAXUO, MAXNH,
     .              MAXO, NUMH, LISTHPTR, LISTH, H, S,
     .              E1, E2, NHIST, SIGMA, 
     .              XIJ, INDXUO, NK, KPOINT, EO,
     .              HAUX, SAUX, PSI, AUX, DTOT, DPR, NUOTOT )
c        IF( NBK .GT. 0) THEN
c         CALL PDOSK( NSPIN, NUO, NO, MAXSPN, MAXUO, MAXNH, 
c     .               MAXO,  NUMH, LISTHPTR, LISTH, H, S,
c     .               E1, E2, NHIST, SIGMA, 
c     .               XIJ, INDXUO, NBK, BK, EBK,
c     .               HAUX, SAUX, PSI, AUX, DTOT, DPR, NUOTOT )
c        ENDIF

        FNAMETOT = PASTE(SNAME,'.DOS')

        CALL IO_ASSIGN(IUNIT1)
          OPEN(UNIT=IUNIT1, FILE=FNAMETOT, FORM='formatted', 
     .         STATUS='unknown') 
          DO IHIST = 1, NHIST
            ENER = E1 + (IHIST-1) * DELTA
            WRITE(IUNIT1,'(3f20.5)') ENER*EV,DTOT(IHIST,1)/EV,
     .           DTOT(IHIST,2)/EV
          ENDDO
        CALL IO_CLOSE(IUNIT1)

CCC New writing
        FNAMEPRO = PASTE(SNAME,'.PDOS')
        CALL IO_ASSIGN(IUNIT2)
        open(iunit2,file=FNAMEPRO,form='formatted',status='unknown')
        write(iunit2,'(a)') '<pdos>'
        write(iunit2,'(a,i1,a)') '<nspin>', nspin, '</nspin>'
        write(iunit2,'(a,i4,a)') '<norbitals>', nuo, '</norbitals>'
        write(iunit2,'(a)') '<energy_values units="eV">'
        do ihist=1,nhist
           ENER = E1 + (IHIST-1) * DELTA
           write(iunit2,'(f20.5)') ener*eV
        enddo
        write(iunit2,'(a)') '</energy_values>'

        do i = 1, nuo
           iat = iaorb(i)
           iorb = iphorb(i)
           spec = isa(iat)

           write(iunit2,'(a)') '<orbital '
           call xml_dump_attribute(iunit2,"index",str(i))
           call xml_dump_attribute(iunit2,"species",
     $                           trim(labelfis(spec)))
           write(pos_string,'(3f11.6)') (xa(ii,iat),ii=1,3)
           call xml_dump_attribute(iunit2,"position",pos_string)
           call xml_dump_attribute(iunit2,"n",str(cnfigfio(spec,iorb)))
           call xml_dump_attribute(iunit2,"l",str(lofio(spec,iorb)))
           call xml_dump_attribute(iunit2,"m",str(mofio(spec,iorb)))
           call xml_dump_attribute(iunit2,"z",str(zetafio(spec,iorb)))
           write(iunit2,'(a)') '> '

           write(iunit2,'(a)') '<data>'
           do ihist=1,nhist
              if (nspin.eq.1) then
                 write(iunit2,'(f20.5)') dpr(ihist,i,1)/eV
              else if (nspin .eq. 2) then
                 write(iunit2,'(2f20.5)') ener*eV, dpr(ihist,i,1)/eV,
     $                                 dpr(ihist,i,2)/eV
              endif
           enddo
           write(iunit2,'(a)') '</data>'
           write(iunit2,'(a)') '</orbital>'
        enddo
        write(iunit2,'(a)') '</pdos>'
        call io_close(iunit2)
      ENDIF
            

C Free local arrays ----------------------------------------------------
      CALL MEMORY('D','D',SIZE(HAUX),'pdos')
      DEALLOCATE(HAUX)
      CALL MEMORY('D','D',SIZE(SAUX),'pdos')
      DEALLOCATE(SAUX)
      CALL MEMORY('D','D',SIZE(PSI),'pdos')
      DEALLOCATE(PSI)
      CALL MEMORY('D','D',SIZE(AUX),'pdos')
      DEALLOCATE(AUX)
      CALL MEMORY('D','I',SIZE(MUO),'pdos')
      DEALLOCATE(MUO)
      CALL MEMORY('D','D',SIZE(DTOT),'pdos')
      DEALLOCATE(DTOT)
      CALL MEMORY('D','D',SIZE(DPR),'pdos')
      DEALLOCATE(DPR)


      CALL TIMER( 'pdos', 2)

      RETURN
      END




C ----------------------------------------------------------------------

      CHARACTER*10 FUNCTION INTCHAR(IA)

      IMPLICIT NONE

      INTEGER IA
      CHARACTER*1 MYCHAR, ICHAR(10)

      INTEGER I, IN, J, IFIN

      DO I = 1, 10
        ICHAR(I) = ' '
      ENDDO
      
      INTCHAR = ' '

      I = IA
 
      DO IN = 1, 10
        IF(I.GT.0) THEN
          J = MOD(I,10)
          I = I/10
          ICHAR(IN) = MYCHAR(J)
          IFIN = IN
        ENDIF
      ENDDO  
        
      DO IN = 1, IFIN
        INTCHAR = ICHAR(IN)//INTCHAR
      ENDDO

      END

      CHARACTER FUNCTION MYCHAR(I)
      CHARACTER CHARI(0:9)
      DATA CHARI/'0','1','2','3','4','5','6','7','8','9'/
      MYCHAR = CHARI(I)
      END

