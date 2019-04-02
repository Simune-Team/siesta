! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

      SUBROUTINE WAVOFR( NA, NO, no_u, MAXNA, NSPIN, nspin_blocks,
     .                   non_coll,
     .                   ISA, IPHORB, INDXUO, LASTO, XA, CELL,
     .                   wf_unit, NK, gamma_wfsx,
     .                   IDIMEN, IOPTION, XMIN, XMAX, YMIN, YMAX, 
     .                   ZMIN, ZMAX, NPX, NPY, NPZ, COORPO, NORMAL, 
     .                   DIRVER1, DIRVER2, 
     .                   ARMUNI, IUNITCD, ISCALE, RMAXO )
C **********************************************************************
C Compute the wave functions at the points of a plane or a 3D grid
C in real space
C Coded by P. Ordejon, from Junquera's rhoofr. July 2003
C **********************************************************************

      use precision
      USE FDF
      USE ATMFUNCS
      USE CHEMICAL
      USE LISTSC_MODULE, ONLY: LISTSC
      use planed, only: plane

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::
     .  NA, NO, NO_U, IOPTION, NPX, NPY, NPZ, ISCALE, IUNITCD,
     .  IDIMEN, NSPIN, nspin_blocks, MAXNA, NK,
     .  ISA(NA), IPHORB(NO), INDXUO(NO), LASTO(0:NA)
      real(dp), INTENT(IN) ::
     .  XA(3,NA), XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX,
     .  COORPO(3,3), NORMAL(3), DIRVER1(3), DIRVER2(3), ARMUNI,
     .  RMAXO
      logical, intent(in) :: non_coll, gamma_wfsx
      integer, intent(in) :: wf_unit
      real(dp), INTENT(IN) :: CELL(3,3)


C ****** INPUT *********************************************************
C INTEGER NA               : Total number of atoms in Supercell
C INTEGER NO               : Total number of orbitals in Supercell
C INTEGER NO_U             : Total number of orbitals in Unit Cell
C INTEGER MAXNA            : Maximum number of neighbours of any atom
C INTEGER NSPIN            : Number of different spin polarizations: 1, 2, 4
! integer nspin_blocks     : blocks in WFSX file       
! logical non_coll         : NC/SOC wfs data?
C INTEGER ISA(NA)          : Species index of each atom
C INTEGER IPHORB(NO)       : Orital index of each orbital in its atom
C INTEGER INDXUO(NO)       : Equivalent orbital in unit cell
C INTEGER LASTO(0:NA)      : Last orbital of each atom in array iphorb
C REAL*8  XA(3,NA)         : Atomic positions in cartesian coordinates
C                            (in bohr)
C REAL*8  CELL(3,3)        : Supercell vectors CELL(IXYZ,IVECT)
C                            (in bohr)
C INTEGER NK               : Number of k-points
C INTEGER IDIMEN           : Specify if the run is to plot quantities
C                            in a plane or in a 3D grid (2 or 3, respect)
C INTEGER IOPTION          : Option to generate the plane
C                            1 = Normal vector
C                            2 = Two vectors contained in the plane
C                            3 = Three points in the plane
C INTEGER NPX,NPY,NPZ      : Number of points along x and y and z
C REAL*8  XMIN, XMAX       : Limits of the plane in the PRF for x-direction
C REAL*8  YMIN, YMAX       : Limits of the plane in the PRF for y-direction
C REAL*8  ZMIN, ZMAX       : Limits of the z-direction
C REAL*8  COORPO(3,3)      : Coordinates of the three points used to define
C                            the plane (Only used if ioption = 3)
C REAL*8  NORMAL(3)        : Components of the normal vector used to define 
C                            the plane
C REAL*8  DIRVER1(3)       : Components of the first vector contained 
C                            in the plane
C                            (Only used if ioption = 2)
C REAL*8  DIRVER2(3)       : Components of the second vector contained 
C                            in the plane
C                            (Only used if ioption = 2)
C REAL*8  ARMUNI           : Conversion factor for the charge density
C INTEGER IUNITCD          : Unit of the charge density
C INTEGER ISCALE           : Unit if the points of the plane
C REAL*8  RMAXO            : Maximum range of basis orbitals
C **********************************************************************

      INTEGER  NPLAMAX, NAPLA, NAINCELL
      INTEGER, DIMENSION(:), ALLOCATABLE ::  INDICES, JNA

      REAL(SP), ALLOCATABLE :: wf_single(:,:)
      COMPLEX(DP), ALLOCATABLE :: wf(:,:)
      complex(dp), allocatable :: CWAVE(:)

      real(dp) :: k(3)
      REAL, DIMENSION(:,:), allocatable :: RWF, IMWF, MWF, PWF

      real(dp), DIMENSION(:), ALLOCATABLE ::  R2IJ
      real(dp), DIMENSION(:,:), ALLOCATABLE ::
     .   PLAPO, POINRE, XAPLA, XIJ, XAINCELL

      INTEGER
     .  NPO, IA, ISEL, NNA, UNITRE1, UNITRE2, UNITIM1, UNITIM2,
     .  UNITPH1, UNITPH2, UNITMO1, UNITMO2
 
      INTEGER
     .  I, J, IN, IAT1, IAT2, JO, IO, IUO, IAVEC, IAVEC1, IAVEC2, 
     .  IS1, IS2, IPHI1, IPHI2, IND, IX, IY, IZ, NX, NY, NZ, IWF, 
     .  INDWF, IZA(NA), IK

      integer :: spinor_comps, ispin, iwf_orig
      integer :: idummy, number_of_wfns
      real(dp) :: ener
      complex(dp) :: expphi
      
      real(dp)
     .  RMAX, XPO(3), RMAX2, XVEC1(3),
     .  XVEC2(3), PHIMU, GRPHIMU(3),
     .  OCELL(3,3), PHASE, SI, CO, PI

      real(dp)
     .  RWAVE, RWAVEUP, RWAVEDN,
     .  IWAVE, IWAVEUP, IWAVEDN,
     .  MWAVE, MWAVEUP, MWAVEDN,
     .  PWAVE, PWAVEUP, PWAVEDN
 
      LOGICAL FIRST

      CHARACTER
     .  SNAME*40, FNAMEWFRE*60, FNAMEWFIM*60, 
     .  FNAMEWFURE*60, FNAMEWFUIM*60, FNAMEWFDRE*60, FNAMEWFDIM*60, 
     .  FNAMEWFMO*60, FNAMEWFPH*60,
     .  FNAMEWFUMO*60, FNAMEWFUPH*60, FNAMEWFDMO*60, FNAMEWFDPH*60,
     .  CHAR1*10, CHAR2*10, ITOCHAR*10, 
     .  EXT*20, EXT2*25

      EXTERNAL
     .  IO_ASSIGN, IO_CLOSE,
     .  NEIGHB, WROUT, ITOCHAR

C **********************************************************************
C INTEGER NPLAMAX          : Maximum number of points in the plane
C REAL*8  PLAPO(NPLAMAX,3) : Coordinates of the points of the plane in PRF
C REAL*8  POINRE(NPLAMAX,3): Coordinates of the points of the plane in Lattice
C                            Reference Frame
C INTEGER NAPLA            : Number of atoms whose coordinates will be rotated
C INTEGER INDICES(NA)      : Indices of the atoms whose coordinates will 
C                            be rotated from the lattice reference frame 
C                            to the in-plane reference frame
C REAL*8  XAPLA(3,NA)      : Atomic coordinates in plane reference frame
C INTEGER IA               : Atom whose neighbours are needed.
C                            A routine initialization must be done by
C                            a first call with IA = 0
C                            If IA0=0, point X0 is used as origin instead
C INTEGER ISEL             : Single-counting switch (0=No, 1=Yes). If ISEL=1,
C                            only neighbours with JA.LE.IA are included in JNA
C INTEGER NNA              : Number of non-zero orbitals at a point in 
C                            real space
C INTEGER JNA(MAXNA)       : Atom index of neighbours. The neighbours
C                            atoms might be in the supercell
C REAL*8  XIJ(3,MAXNA)     : Vectors from point in real space to orbitals
C REAL*8  R2IJ(MAXNA)      : Squared distance to atomic orbitals
C REAL*8  XPO(3)           : Coordinates of the point of the plane respect
C                            we are going to calculate the neighbours orbitals
C REAL RWF(NPO,NSPIN)      : Wave fnctn at each point of the grid (real part)
C REAL IMWF(NPO,NSPIN)     : Wave fnctn at each point of the grid (imag part)
C INTEGER IZA(NA)          : Atomic number of each atom
C **********************************************************************

      ! The first dimension of wf_single is the number of real numbers per orbital
      ! to be read from the WFSX file:
      ! 1 for real wfs, 2 for complex, and four for the two spinor components
      ! wf is a complex array which holds either a wfn or a two-component spinor.

      if (non_coll) then
        allocate(wf_single(4,1:no_u))
        allocate(wf(1:no_u,2))
        spinor_comps = 2
      else
        spinor_comps = 1
        if (gamma_wfsx) then
           allocate(wf_single(1,1:no_u))
           allocate(wf(1:no_u,1))
        else
           allocate(wf_single(2,1:no_u))
           allocate(wf(1:no_u,1))
        endif
      endif
      allocate(CWAVE(spinor_comps))

C     Allocate some variables ---------------------------------------------

      PI = 4.0D0 * ATAN(1.0D0)

      NPLAMAX = NPX * NPY * NPZ

      ALLOCATE(PLAPO(NPLAMAX,3))
      ALLOCATE(POINRE(NPLAMAX,3))
      ALLOCATE(INDICES(NA))
      ALLOCATE(XAPLA(3,NA))

      ALLOCATE(XAINCELL(3,NA))

C Build the plane ------------------------------------------------------
          CALL PLANE( NA, NPLAMAX, IDIMEN, IOPTION, 
     .                XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, 
     .                NPX, NPY, NPZ, COORPO, NORMAL, 
     .                DIRVER1, DIRVER2, 
     .                XA, NAPLA, INDICES, ISCALE,
     .                POINRE, PLAPO, XAPLA )   

C Initialize neighbour subroutine --------------------------------------
      IA = 0
      ISEL = 0
      RMAX = RMAXO
      NNA  = MAXNA
      IF (ALLOCATED(JNA)) THEN
        DEALLOCATE(JNA)
      ENDIF
      IF (ALLOCATED(R2IJ)) THEN
        DEALLOCATE(R2IJ)
      ENDIF
      IF (ALLOCATED(XIJ)) THEN
        DEALLOCATE(XIJ)
      ENDIF
      ALLOCATE(JNA(MAXNA))
      ALLOCATE(R2IJ(MAXNA))
      ALLOCATE(XIJ(3,MAXNA))

      FIRST = .TRUE.
      DO I = 1,3
        XPO(I) = 0.D0
      ENDDO

      CALL NEIGHB( CELL, RMAX, NA, XA, XPO, IA, ISEL, 
     .             NNA, JNA, XIJ, R2IJ, FIRST )
 
      FIRST = .FALSE.
      RMAX2 =  RMAXO**2

! Stream over wavefunctions in file

      DO IK  = 1, NK
        do ispin = 1, nspin_blocks
         read(wf_unit) idummy, k(1:3)
            if (idummy /= ik) then
               write(6,*) "ik index mismatch in WFS file"
               WRITE(6,*) "ik in file, ik: ", idummy, ik
            endif
         read(wf_unit) idummy
            if (idummy /= ispin) then
               write(6,*) "ispin index mismatch in WFS file"
               WRITE(6,*) "ispin in file, ispin: ", idummy, ispin
            endif
         read(wf_unit) number_of_wfns

         WRITE(6,*) 'stm:  Processing kpoint ',IK
         WRITE(6,*) 'stm:  nwf: ', number_of_wfns
         WRITE(6,*) '     --------------------------------'
         DO IWF = 1, number_of_wfns
            read(wf_unit) iwf_orig
            if (iwf_orig /= iwf) then
               ! The file holds a subset of wfs, with the original indexes...
               WRITE(6,*) 'Original wf index: ', iwf_orig
            endif
            read(wf_unit) ener

            read(wf_unit) (wf_single(:,io), io=1,no_u)
            ! Use a double precision complex form in what follows
            if ( non_coll) then
               wf(:,1) = cmplx(wf_single(1,:), wf_single(2,:), kind=dp)
               wf(:,2) = cmplx(wf_single(3,:), wf_single(4,:), kind=dp)
            else
               if (gamma_wfsx) then
                  wf(:,1) = cmplx(wf_single(1,:), 0.0_sp, kind=dp)
               else
                  wf(:,1) = cmplx(wf_single(1,:),wf_single(2,:),kind=dp)
               endif
            endif


            CHAR1 = ITOCHAR(Iwf_Orig)
            CHAR2 = ITOCHAR(IK)
!     Open files to store wave functions -----------------
            SNAME = FDF_STRING('SystemLabel','siesta')

         IF (NSPIN .EQ. 1) THEN

          IF (IDIMEN .EQ. 2) THEN
            FNAMEWFRE = TRIM(SNAME)//'.CON.K' //
     $            trim(char2) // '.WF.' // trim(CHAR1)
            FNAMEWFPH = TRIM(FNAMEWFRE)//'.PHASE'
            FNAMEWFMO = TRIM(FNAMEWFRE)//'.MOD'
            FNAMEWFIM = TRIM(FNAMEWFRE)//'.IMAG'
            FNAMEWFRE = TRIM(FNAMEWFRE)//'.REAL'
          ELSEIF (IDIMEN .EQ. 3) THEN
            FNAMEWFRE = TRIM(SNAME)//'.K' //
     $            trim(char2) // '.WF.' // trim(CHAR1)
            FNAMEWFPH = TRIM(FNAMEWFRE)//'.PHASE.cube'
            FNAMEWFMO = TRIM(FNAMEWFRE)//'.MOD.cube'
            FNAMEWFIM = TRIM(FNAMEWFRE)//'.IMAG.cube'
            FNAMEWFRE = TRIM(FNAMEWFRE)//'.REAL.cube'
          ENDIF

          CALL IO_ASSIGN(UNITRE1)
          OPEN(UNIT = UNITRE1, FILE = FNAMEWFRE, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITRE1)
          CALL IO_ASSIGN(UNITIM1)
          OPEN(UNIT = UNITIM1, FILE = FNAMEWFIM, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITIM1)
          CALL IO_ASSIGN(UNITMO1)
          OPEN(UNIT = UNITMO1, FILE = FNAMEWFMO, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITMO1)
          CALL IO_ASSIGN(UNITPH1)
          OPEN(UNIT = UNITPH1, FILE = FNAMEWFPH, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITPH1)

       ELSEIF ((NSPIN .EQ. 2) .or. (NSPIN == 4))  THEN
          ! We will reuse 'up' and 'down' for the spinor components
           IF (IDIMEN .EQ. 2) THEN
            FNAMEWFURE = TRIM(SNAME)//'.CON.K' //
     $            trim(char2) // '.WF.' // trim(CHAR1)
            FNAMEWFDRE = TRIM(FNAMEWFURE)//'.DOWN'
            FNAMEWFDPH = TRIM(FNAMEWFDRE)//'.PHASE'
            FNAMEWFDMO = TRIM(FNAMEWFDRE)//'.MOD'
            FNAMEWFDIM = TRIM(FNAMEWFDRE)//'.IMAG'
            FNAMEWFDRE = TRIM(FNAMEWFDRE)//'.REAL'

            FNAMEWFURE = TRIM(FNAMEWFURE)//'.UP'
            FNAMEWFUPH = TRIM(FNAMEWFURE)//'.PHASE'
            FNAMEWFUMO = TRIM(FNAMEWFURE)//'.MOD'
            FNAMEWFUIM = TRIM(FNAMEWFURE)//'.IMAG'
            FNAMEWFURE = TRIM(FNAMEWFURE)//'.REAL'

         ELSE IF (IDIMEN .EQ. 3) THEN
            FNAMEWFURE = TRIM(SNAME)//'.K' //
     $            trim(char2) // '.WF.' // trim(CHAR1)
            FNAMEWFDPH = TRIM(FNAMEWFURE)//'.DOWN.PHASE.cube'
            FNAMEWFDMO = TRIM(FNAMEWFURE)//'.DOWN.MOD.cube'
            FNAMEWFDIM = TRIM(FNAMEWFURE)//'.DOWN.IMAG.cube'
            FNAMEWFDRE = TRIM(FNAMEWFURE)//'.DOWN.REAL.cube'
            FNAMEWFUPH = TRIM(FNAMEWFURE)//'.UP.PHASE.cube'
            FNAMEWFUMO = TRIM(FNAMEWFURE)//'.UP.MOD.cube'
            FNAMEWFUIM = TRIM(FNAMEWFURE)//'.UP.IMAG.cube'
            FNAMEWFURE = TRIM(FNAMEWFURE)//'.UP.REAL.cube'
          ENDIF

          CALL IO_ASSIGN(UNITRE1)
          OPEN(UNIT = UNITRE1, FILE = FNAMEWFURE, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITRE1)
          CALL IO_ASSIGN(UNITRE2)
          OPEN(UNIT = UNITRE2, FILE = FNAMEWFDRE, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITRE2)
          CALL IO_ASSIGN(UNITIM1)
          OPEN(UNIT = UNITIM1, FILE = FNAMEWFUIM, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITIM1)
          CALL IO_ASSIGN(UNITIM2)
          OPEN(UNIT = UNITIM2, FILE = FNAMEWFDIM, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITIM2)
          CALL IO_ASSIGN(UNITMO1)
          OPEN(UNIT = UNITMO1, FILE = FNAMEWFUMO, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITMO1)
          CALL IO_ASSIGN(UNITMO2)
          OPEN(UNIT = UNITMO2, FILE = FNAMEWFDMO, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITMO2)
          CALL IO_ASSIGN(UNITPH1)
          OPEN(UNIT = UNITPH1, FILE = FNAMEWFUPH, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITPH1)
          CALL IO_ASSIGN(UNITPH2)
          OPEN(UNIT = UNITPH2, FILE = FNAMEWFDPH, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITPH2)
        ELSE
          WRITE(6,*)'BAD NUMBER NSPIN IN WAVOFR.F'
          WRITE(6,*)'NSPIN = ',NSPIN
          WRITE(6,*)'IT MUST BE 1, 2, or 4'
          STOP
        ENDIF

        IF (IDIMEN .EQ. 2) THEN
C Select all atoms in list to be printed out
          NAINCELL=NAPLA
          DO IA=1,NAPLA
            DO IX=1,3
              XAINCELL(IX,IA)=XAPLA(IX,IA)
            ENDDO
          ENDDO

        ELSE IF (IDIMEN .EQ.3) THEN
          DO IX = 1,3
            DO IY = 1,3
              OCELL(IX,IY)=0.D0
            ENDDO
          ENDDO
C   Determine cell size
          OCELL(1,1) = DABS(XMAX-XMIN)
          OCELL(2,2) = DABS(YMAX-YMIN)
          OCELL(3,3) = DABS(ZMAX-ZMIN)
C   Determine atoms which are within the plotting box
          NAINCELL=0
          DO IA=1,NA
            IF ((XAPLA(1,IA).LT.XMIN*1.1).OR.(XAPLA(1,IA).GT.XMAX*1.1)
     .     .OR. (XAPLA(2,IA).LT.YMIN*1.1).OR.(XAPLA(2,IA).GT.YMAX*1.1)
     .     .OR. (XAPLA(3,IA).LT.ZMIN*1.1).OR.(XAPLA(3,IA).GT.ZMAX*1.1))
     .      GOTO 90
            NAINCELL=NAINCELL+1
            IZA(NAINCELL) = ATOMIC_NUMBER(ISA(IA))
            DO IX=1,3
              XAINCELL(IX,NAINCELL)=XAPLA(IX,IA)
            ENDDO
90          CONTINUE
          ENDDO

          IF (NSPIN .EQ. 1) THEN
            call write_cube_header(unitre1,fnamewfre)
            call write_cube_header(unitim1,fnamewfim)
            call write_cube_header(unitmo1,fnamewfmo)
            call write_cube_header(unitph1,fnamewfph)
      
          ELSE IF (NSPIN .EQ. 2) THEN
             if (ispin == 1 ) then
                call write_cube_header(unitre1,fnamewfure)
                call write_cube_header(unitim1,fnamewfuim)
                call write_cube_header(unitmo1,fnamewfumo)
                call write_cube_header(unitph1,fnamewfuph)
             else
                call write_cube_header(unitre2,fnamewfdre)
                call write_cube_header(unitim2,fnamewfdim)
                call write_cube_header(unitmo2,fnamewfdmo)
                call write_cube_header(unitph2,fnamewfdph)
             endif
          ELSE IF (NSPIN .EQ. 4) THEN
                call write_cube_header(unitre1,fnamewfure)
                call write_cube_header(unitim1,fnamewfuim)
                call write_cube_header(unitmo1,fnamewfumo)
                call write_cube_header(unitph1,fnamewfuph)

                call write_cube_header(unitre2,fnamewfdre)
                call write_cube_header(unitim2,fnamewfdim)
                call write_cube_header(unitmo2,fnamewfdmo)
                call write_cube_header(unitph2,fnamewfdph)

          ENDIF
        ENDIF


        CALL WROUT(IDIMEN, .FALSE., .TRUE., IOPTION, NORMAL, COORPO,
     .             DIRVER1, DIRVER2,
     .             NPX, NPY, NPZ, XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX,
     .             IUNITCD,  NA, NAINCELL, INDICES, XAINCELL )

        WRITE(6,'(A)')
        WRITE(6,'(A,i0,a,i0)')
     .    '   Generating grid values for wf... k-point:',IK,
     $       ' wf #:', iwf_orig


C Allocate space for wave functions in 3D-grid

        IF (.NOT.ALLOCATED(RWF)) THEN
          ALLOCATE(RWF(NPX*NPY*NPZ,NSPIN))
          ALLOCATE(IMWF(NPX*NPY*NPZ,NSPIN))
          ALLOCATE(MWF(NPX*NPY*NPZ,NSPIN))
          ALLOCATE(PWF(NPX*NPY*NPZ,NSPIN))
        ENDIF

      
C Loop over all points in real space -----------------------------------

        NPO = 0
        DO 102 NZ = 1,NPZ
        DO 101 NY = 1,NPY
        DO 100 NX = 1,NPX
          NPO = NPO + 1

          ! Initialize the wave function at each point -----------------------
          CWAVE(:) = 0.0_dp

C Localize non-zero orbitals at each point in real space ---------------
          DO IX = 1,3
            XPO(IX) = POINRE(NPO,IX)
          ENDDO
     
          IA   = 0
          ISEL = 0
          NNA  = MAXNA

          CALL NEIGHB( CELL, RMAX, NA, XA, XPO, IA, ISEL, 
     .                 NNA, JNA, XIJ, R2IJ, FIRST )

C Loop over Non-zero orbitals ------------------------------------------ 
           DO IAT1 = 1, NNA
             IF( R2IJ(IAT1) .GT. RMAX2 ) EXIT

             IAVEC1   = JNA(IAT1)
             IS1      = ISA(IAVEC1)
             XVEC1(1) = -XIJ(1,IAT1)
             XVEC1(2) = -XIJ(2,IAT1)
             XVEC1(3) = -XIJ(3,IAT1)

             PHASE = K(1)*(XPO(1)+XIJ(1,IAT1))+
     .               K(2)*(XPO(2)+XIJ(2,IAT1))+
     .               K(3)*(XPO(3)+XIJ(3,IAT1))

             SI=SIN(PHASE)
             CO=COS(PHASE)
             EXPPHI=CMPLX(CO,SI,dp)

             DO IO = LASTO(IAVEC1-1) + 1, LASTO(IAVEC1)
               IPHI1 = IPHORB(IO)
               IUO   = INDXUO(IO)
               CALL PHIATM( IS1, IPHI1, XVEC1, PHIMU, GRPHIMU )
               ! Note implicit loop over spinor components
               CWAVE(:) = CWAVE(:) + PHIMU * WF(iuo,:) * EXPPHI 

            ENDDO       ! orbitals
         ENDDO ! atoms
           
           IF ( NSPIN .EQ. 1 ) THEN

              rwave = real(cwave(1), dp)
              iwave = aimag(cwave(1))

              call mod_and_phase(rwave,iwave,mwave,pwave)

      
           ELSEIF (NSPIN .EQ. 2) THEN

              if (ispin == 1) then
                 rwaveup = real(cwave(1), dp)
                 iwaveup = aimag(cwave(1))
                 call mod_and_phase(rwaveup,iwaveup,mwaveup,pwaveup)
              else
                 rwavedn = real(cwave(1), dp)
                 iwavedn = aimag(cwave(1))
                 call mod_and_phase(rwavedn,iwavedn,mwavedn,pwavedn)
              endif

           ELSEIF (NSPIN .EQ. 4) THEN

              rwaveup = real(cwave(1), dp)
              iwaveup = aimag(cwave(1))
              rwavedn = real(cwave(2), dp)
              iwavedn = aimag(cwave(2))
              call mod_and_phase(rwaveup,iwaveup,mwaveup,pwaveup)
              call mod_and_phase(rwavedn,iwavedn,mwavedn,pwavedn)

           ENDIF

           IF (IDIMEN .EQ. 2) THEN
             IF ( NSPIN .EQ. 1 ) THEN
                WRITE(UNITRE1,'(3F12.5)')
     .               PLAPO(NPO,1),PLAPO(NPO,2),RWAVE*SQRT(ARMUNI)
                WRITE(UNITIM1,'(3F12.5)')
     .               PLAPO(NPO,1),PLAPO(NPO,2),IWAVE*SQRT(ARMUNI)
                WRITE(UNITMO1,'(3F12.5)')
     .               PLAPO(NPO,1),PLAPO(NPO,2),MWAVE*SQRT(ARMUNI)
                WRITE(UNITPH1,'(3F12.5)')
     .               PLAPO(NPO,1),PLAPO(NPO,2),PWAVE
             ELSEIF ( NSPIN .EQ. 2 ) THEN
                if (ispin == 1) then              
                   WRITE(UNITRE1,'(3F12.5)')
     .                  PLAPO(NPO,1),PLAPO(NPO,2),RWAVEUP*SQRT(ARMUNI)
                   WRITE(UNITIM1,'(3F12.5)')
     .                  PLAPO(NPO,1),PLAPO(NPO,2),IWAVEUP*SQRT(ARMUNI)
                   WRITE(UNITMO1,'(3F12.5)')
     .                  PLAPO(NPO,1),PLAPO(NPO,2),MWAVEUP*SQRT(ARMUNI)
                   WRITE(UNITPH1,'(3F12.5)')
     .                  PLAPO(NPO,1),PLAPO(NPO,2),PWAVEUP
                else
                   WRITE(UNITRE2,'(3F12.5)')
     .                  PLAPO(NPO,1),PLAPO(NPO,2),RWAVEDN*SQRT(ARMUNI)
                   WRITE(UNITIM2,'(3F12.5)')
     .                  PLAPO(NPO,1),PLAPO(NPO,2),IWAVEDN*SQRT(ARMUNI)
                   WRITE(UNITMO2,'(3F12.5)')
     .                  PLAPO(NPO,1),PLAPO(NPO,2),MWAVEDN*SQRT(ARMUNI)
                   WRITE(UNITPH2,'(3F12.5)')
     .                  PLAPO(NPO,1),PLAPO(NPO,2),PWAVEDN
                endif
             ELSEIF ( NSPIN .EQ. 4 ) THEN
                WRITE(UNITRE1,'(3F12.5)')
     .               PLAPO(NPO,1),PLAPO(NPO,2),RWAVEUP*SQRT(ARMUNI)
                WRITE(UNITIM1,'(3F12.5)')
     .               PLAPO(NPO,1),PLAPO(NPO,2),IWAVEUP*SQRT(ARMUNI)
                WRITE(UNITMO1,'(3F12.5)')
     .               PLAPO(NPO,1),PLAPO(NPO,2),MWAVEUP*SQRT(ARMUNI)
                WRITE(UNITPH1,'(3F12.5)')
     .               PLAPO(NPO,1),PLAPO(NPO,2),PWAVEUP

                WRITE(UNITRE2,'(3F12.5)')
     .               PLAPO(NPO,1),PLAPO(NPO,2),RWAVEDN*SQRT(ARMUNI)
                WRITE(UNITIM2,'(3F12.5)')
     .               PLAPO(NPO,1),PLAPO(NPO,2),IWAVEDN*SQRT(ARMUNI)
                WRITE(UNITMO2,'(3F12.5)')
     .               PLAPO(NPO,1),PLAPO(NPO,2),MWAVEDN*SQRT(ARMUNI)
                WRITE(UNITPH2,'(3F12.5)')
     .               PLAPO(NPO,1),PLAPO(NPO,2),PWAVEDN
             ENDIF
          ELSE IF (IDIMEN .EQ. 3) THEN
             IF (NSPIN .EQ. 1) THEN
                RWF(NPO,1) = RWAVE*SQRT(ARMUNI)
                IMWF(NPO,1) = IWAVE*SQRT(ARMUNI)
                MWF(NPO,1) = MWAVE*SQRT(ARMUNI)
                PWF(NPO,1) = PWAVE
             ELSE IF (NSPIN .EQ.2) THEN
               if (ispin == 1) then
                  RWF(NPO,1) = RWAVEUP*SQRT(ARMUNI)
                  IMWF(NPO,1) = IWAVEUP*SQRT(ARMUNI)
                  MWF(NPO,1) = MWAVEUP*SQRT(ARMUNI)
                  PWF(NPO,1) = PWAVEUP
               else
                  RWF(NPO,2) = RWAVEDN*SQRT(ARMUNI)
                  IMWF(NPO,2) = IWAVEDN*SQRT(ARMUNI)
                  MWF(NPO,2) = MWAVEDN*SQRT(ARMUNI)
                  PWF(NPO,2) = PWAVEDN
               endif
             ELSE IF (NSPIN .EQ. 4) THEN

                RWF(NPO,1) = RWAVEUP*SQRT(ARMUNI)
                IMWF(NPO,1) = IWAVEUP*SQRT(ARMUNI)
                MWF(NPO,1) = MWAVEUP*SQRT(ARMUNI)
                PWF(NPO,1) = PWAVEUP

                RWF(NPO,2) = RWAVEDN*SQRT(ARMUNI)
                IMWF(NPO,2) = IWAVEDN*SQRT(ARMUNI)
                MWF(NPO,2) = MWAVEDN*SQRT(ARMUNI)
                PWF(NPO,2) = PWAVEDN

            ENDIF
         ENDIF

           IF (IDIMEN .EQ. 2) THEN
             IF ( MOD(NPO,NPX) .EQ. 0 ) THEN
               WRITE(UNITRE1,*)
               WRITE(UNITIM1,*)
               IF ( NSPIN .EQ. 2 ) THEN
                 WRITE(UNITRE2,*)
                 WRITE(UNITIM2,*)
               ENDIF
             ENDIF
           ENDIF


C End x loop
 100    ENDDO  
C End y and z loops
 101    ENDDO  
 102    ENDDO  

        IF (IDIMEN .EQ. 3) THEN
          IF (NSPIN .EQ. 1) THEN
            DO NX=1,NPX
              DO NY=1,NPY
                WRITE(UNITRE1,'(6e13.5)')
     .            (RWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,1),NZ=1,NPZ)
                WRITE(UNITIM1,'(6e13.5)')
     .            (IMWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,1),NZ=1,NPZ)
                WRITE(UNITMO1,'(6e13.5)')
     .            (MWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,1),NZ=1,NPZ)
                WRITE(UNITPH1,'(6e13.5)')
     .            (PWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,1),NZ=1,NPZ)
              ENDDO
            ENDDO
         ELSE IF (NSPIN .EQ. 2) THEN
            if (ispin == 1) then
               DO NX=1,NPX
                  DO NY=1,NPY
                     WRITE(UNITRE1,'(6e13.5)')
     .                    (RWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,1),NZ=1,NPZ)
                     WRITE(UNITIM1,'(6e13.5)')
     .                   (IMWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,1),NZ=1,NPZ)
                     WRITE(UNITMO1,'(6e13.5)')
     .                    (MWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,1),NZ=1,NPZ)
                     WRITE(UNITPH1,'(6e13.5)')
     .                    (PWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,1),NZ=1,NPZ)
                  ENDDO
               ENDDO
            else
               DO NX=1,NPX
                  DO NY=1,NPY
                     WRITE(UNITRE2,'(6e13.5)')
     .                    (RWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,2),NZ=1,NPZ)
                     WRITE(UNITIM2,'(6e13.5)')
     .                   (IMWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,2),NZ=1,NPZ)
                     WRITE(UNITMO2,'(6e13.5)')
     .                    (MWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,2),NZ=1,NPZ)
                     WRITE(UNITPH2,'(6e13.5)')
     .                    (PWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,2),NZ=1,NPZ)
                  ENDDO
               ENDDO
            endif
         ELSE IF (NSPIN .EQ. 4) THEN

               DO NX=1,NPX
                  DO NY=1,NPY
                     WRITE(UNITRE1,'(6e13.5)')
     .                    (RWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,1),NZ=1,NPZ)
                     WRITE(UNITIM1,'(6e13.5)')
     .                   (IMWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,1),NZ=1,NPZ)
                     WRITE(UNITMO1,'(6e13.5)')
     .                    (MWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,1),NZ=1,NPZ)
                     WRITE(UNITPH1,'(6e13.5)')
     .                    (PWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,1),NZ=1,NPZ)
                  ENDDO
               ENDDO

               DO NX=1,NPX
                  DO NY=1,NPY
                     WRITE(UNITRE2,'(6e13.5)')
     .                    (RWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,2),NZ=1,NPZ)
                     WRITE(UNITIM2,'(6e13.5)')
     .                   (IMWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,2),NZ=1,NPZ)
                     WRITE(UNITMO2,'(6e13.5)')
     .                    (MWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,2),NZ=1,NPZ)
                     WRITE(UNITPH2,'(6e13.5)')
     .                    (PWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,2),NZ=1,NPZ)
                  ENDDO
               ENDDO

         ENDIF

          WRITE(6,'(A)')
          WRITE(6,'(A)')
     .      '   Your output files are:'
          IF (NSPIN .EQ. 1) THEN
            WRITE(6,'(A,A)') '   ',FNAMEWFRE
            WRITE(6,'(A,A)') '   ',FNAMEWFIM
            WRITE(6,'(A,A)') '   ',FNAMEWFMO
            WRITE(6,'(A,A)') '   ',FNAMEWFPH
          ELSE IF ((NSPIN .EQ. 2) .or. (NSPIN .EQ. 4)) THEN
            WRITE(6,'(A,A)') '   ',FNAMEWFURE
            WRITE(6,'(A,A)') '   ',FNAMEWFUIM
            WRITE(6,'(A,A)') '   ',FNAMEWFDRE
            WRITE(6,'(A,A)') '   ',FNAMEWFDIM
            WRITE(6,'(A,A)') '   ',FNAMEWFUMO
            WRITE(6,'(A,A)') '   ',FNAMEWFUPH
            WRITE(6,'(A,A)') '   ',FNAMEWFDMO
            WRITE(6,'(A,A)') '   ',FNAMEWFDPH
            if (nspin == 4) write(6,"(a)")
     $           "... up and down for spinor components"
          ENDIF
        ENDIF ! 3D


        CALL IO_CLOSE(UNITRE1)
        CALL IO_CLOSE(UNITIM1)
        IF (NSPIN >= 2) CALL IO_CLOSE(UNITRE2)
        IF (NSPIN >= 2) CALL IO_CLOSE(UNITIM2)
        CALL IO_CLOSE(UNITMO1)
        CALL IO_CLOSE(UNITPH1)
        IF (NSPIN >= 2) CALL IO_CLOSE(UNITMO2)
        IF (NSPIN >= 2) CALL IO_CLOSE(UNITPH2)
     
          
      ENDDO  ! wfn
      ENDDO  ! spin_block
      ENDDO  ! k-point

      DEALLOCATE(RWF)
      DEALLOCATE(IMWF)
      DEALLOCATE(MWF)
      DEALLOCATE(PWF)

        CONTAINS
      
       subroutine write_cube_header(lun,fname)
       integer, intent(in) :: lun
       character(len=*), intent(in) :: fname
       ! rest by host association
      
            WRITE(LUN,*) FNAME
            WRITE(LUN,*) FNAME
            WRITE(LUN,'(i5,4f12.6)') NAINCELL, XMIN, YMIN, ZMIN
            WRITE(LUN,'(i5,4f12.6)') NPX,(OCELL(1,J)/(NPX-1),J=1,3)
            WRITE(LUN,'(i5,4f12.6)') NPY,(OCELL(2,J)/(NPY-1),J=1,3)
            WRITE(LUN,'(i5,4f12.6)') NPZ,(OCELL(3,J)/(NPZ-1),J=1,3)
            DO IA = 1,NAINCELL
              WRITE(LUN,'(i5,4f12.6)') IZA(IA),0.0,
     .                                   (XAINCELL(IX,IA),IX=1,3)
            ENDDO
       end subroutine write_cube_header

      END

      subroutine mod_and_phase(rwave,iwave,mwave,pwave)
      integer, parameter :: dp = selected_real_kind(10,100)
      real(dp), intent(in) :: rwave
      real(dp), intent(inout) :: iwave
      real(dp), intent(out) :: mwave, pwave

      real(dp), parameter :: pi = 3.1141592653589_dp
      
      MWAVE = DSQRT(RWAVE**2 + IWAVE**2)
      IF (DABS(IWAVE) .LT. 1.D-6) IWAVE = DABS(IWAVE)
      IF (DABS(RWAVE) .LT. 1.D-12) THEN
         IF (IWAVE .GT. 0.0D0) PWAVE = PI/2.0D0
         IF (IWAVE .LT. 0.0D0) PWAVE = -PI/2.0D0
      ELSE
         PWAVE = DATAN(IWAVE/RWAVE)
      ENDIF
      IF (RWAVE .LT. 0.0D0 .AND. IWAVE .GE. 0.0d0)
     .     PWAVE = PWAVE + PI
      IF (RWAVE .LT. 0.0D0 .AND. IWAVE .LT. 0.0D0)
     .     PWAVE = PWAVE - PI
      end
