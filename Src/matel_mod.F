!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
!
      module matel_mod
      use precision, only : dp
      use alloc,     only : re_alloc, de_alloc, alloc_default,
     &                      allocDefaults
      use sys,       only: die
      use parallel,  only : Node, Nodes
      use m_radfft,  only : radfft, reset_radfft
      use m_matel_registry, ONLY : EVALUATE, LCUT, RCUT
      use spher_harm, only : RLYLM, YLMEXP, ylmylm, lofilm
      use spher_harm, only : reset_spher_harm
      use mpi_siesta,  only : MPI_INTEGER, MPI_COMM_WORLD
      use mpi_siesta,  only : MPI_DOUBLE_PRECISION, MPI_IN_PLACE
      private
      public :: init_matel
      public :: get_matel_s
      public :: get_matel_t

C     Module Variables
      integer           :: no     ! orbitals + KB projectors
      integer           :: num_kb ! Number of different KB projectors
      integer           :: num_or ! Number of different orbitals

C     parameters
      integer,          parameter :: NQ        =  1024
      integer,          parameter :: NR        =  NQ
      integer,          parameter :: NRTAB     =  1024
      real(dp),         parameter :: Q2CUT     =  2.50e3_dp
      real(dp),         parameter :: EXPAND    =  1.20_dp
      real(dp),         parameter :: FFTOL     =  1.e-8_dp
      integer,          parameter :: MINEXPAND =  32
      CHARACTER(LEN=*), parameter :: MYNAME =  'MATEL_MOD'

      real(dp)          :: RMAX
      real(dp)          :: PI
      real(dp)          :: QMAX
      real(dp)          :: DQ
      real(dp)          :: DR
      real(dp)          :: DRTAB

      TYPE MATEL
      integer           :: MODE
      logical           :: initiallized = .false.

      integer           :: NC                ! Number of columns of matrix

      integer           :: MF, NF, NF_MINOR
      real(dp), POINTER :: F(:,:)
      integer,  POINTER :: ILM(:), INDF(:)

      integer           :: MFFR, NFFR
      real(dp), POINTER :: FFR(:,:,:)
      integer,  POINTER :: INDFFR(:)

      integer           :: MFFY, NFFY
      real(dp), POINTER :: FFY(:,:)
      integer,  POINTER :: INDFFY(:), ILMFF(:)

      integer           :: MILM
      real(dp), POINTER :: Y(:), DYDR(:,:)

      integer           :: L1L2
      integer,  pointer :: IFFR(:)
      real(dp), pointer :: FFQ(:), FFL(:), CFFR(:)
      
      CONTAINS
      procedure :: INIT
      procedure :: ORBITALTOHARMONIC
      procedure :: COMPUTE_SPHERICAL_HARMONICS
      procedure :: SWAP_SPHERICAL_HARMONICS
      procedure :: FIND_RADIAL_EXPANSION
      procedure :: SWAP_RADIAL_EXP
      procedure :: GET_MATEL

      END TYPE MATEL

      TYPE(MATEL) :: S
      TYPE(MATEL) :: T

      contains

      subroutine INIT( this, IOPER )
C *******************************************************************
C ************************* INPUT ***********************************
C ************************* OUTPUT **********************************
      use precision, only : dp
      implicit none
      class(MATEL)        :: this
      integer             :: IOPER
C     Local variables
      integer             :: IG, IG1, IG2, ILOCAL, IFLM1, IFLM2,
     &                       fi, la, lorb
      real(dp), external  :: MPI_WTIME
      real(dp)            :: lo_time

      call MPI_Barrier( MPI_Comm_World, IG )
      lo_time = MPI_WTIME( )

      ! Set the MODE
      this%MODE = IOPER
      if (this%MODE==1) then
        this%NC = NO
      else
        this%NC = num_or
      endif

      ! Initiallize only once
      if (this%initiallized) return
      this%initiallized = .true.

      ! Distribute Orbitals between processes
      call GET_LOOP_LIMITS( this%NC, NODE, fi, la )
      lorb = la - fi + 1
      ! Allocate memory for spherical harmonics
      this%MF = MAX(INT(this%NC*EXPAND),this%NC+MINEXPAND)
      CALL RE_ALLOC( this%F, 0, NQ, 1, this%MF, 'F', MYNAME )
      CALL RE_ALLOC( this%ILM, 1, this%MF, 'ILM', MYNAME )
      CALL RE_ALLOC( this%INDF, 1, this%NC+1, 'INDF', MYNAME )
      this%NF = 0
      do IG = fi, la
        this%INDF(IG) = this%NF
        call this%compute_spherical_harmonics( IG )
        this%INDF(IG) = this%NF - this%INDF(IG)
      enddo
      if (NODE==0) write(*,*) "TIME(sph_Ha)", MPI_Wtime() - lo_time
      lo_time = MPI_Wtime()

      call this%swap_spherical_harmonics( )
      if (NODE==0) write(*,*) "TIME(swap_H)", MPI_Wtime() - lo_time
      lo_time = MPI_Wtime()

      ! Find radial expansion of overlap ----------------------------------
      ! Get loop indeces
      this%NF_MINOR = this%INDF(this%NC+1)-1
      call GET_LOOP_LIMITS( this%NC*NUM_OR, NODE, fi, la )
      lorb =  la-fi+1

      nullify(this%FFR,this%INDFFR)
      this%MFFR = MAX( INT(lorb*EXPAND), lorb+MINEXPAND )
      CALL RE_ALLOC( this%FFR, 0, NRTAB, 1, 2, 1,
     &  this%MFFR, 'FFR', MYNAME )
      this%NFFR = 0

      nullify(this%FFY,this%ILMFF,this%INDFFY)
      this%MFFY = this%MFFR
      CALL RE_ALLOC( this%FFY, 1, 1, 1, this%MFFY, 'FFY', MYNAME )
      CALL RE_ALLOC( this%INDFFR, 1, this%MFFY, 'INDFFR', MYNAME )
      CALL RE_ALLOC( this%ILMFF, 1, this%MFFY, 'ILMFF', MYNAME )
      CALL RE_ALLOC( this%INDFFY, 1, this%NF*this%NF_MINOR+1,
     &  'INDFFY', MYNAME )
      this%NFFY = 0

      nullify(this%Y,this%DYDR)
      this%MILM = (MINEXPAND*MINEXPAND+1)**2
      CALL RE_ALLOC( this%Y, 1, this%MILM, 'Y', MYNAME )
      CALL RE_ALLOC( this%DYDR, 1, 3, 1, this%MILM, 'DYDR', MYNAME )

      nullify(this%IFFR,this%CFFR,this%FFQ,this%FFL)
      this%L1L2 = MINEXPAND
      CALL RE_ALLOC( this%IFFR, 0, this%L1L2, 'IFFR', MYNAME )
      CALL RE_ALLOC( this%CFFR, 0, this%L1L2, 'CFFR', MYNAME )
      CALL RE_ALLOC( this%FFQ, 0, NQ, 'FFQ', MYNAME )
      CALL RE_ALLOC( this%FFL, 0, NQ, 'FFL', MYNAME )

      DO IG = fi, la
        IG1 = MOD(IG-1,this%NC)+1
        IG2 = (IG-1)/this%NC + 1
        ILOCAL = this%OrbitalToHarmonic( IG1, IG2 )
        DO IFLM1 = this%INDF(IG1),this%INDF(IG1+1)-1
          DO IFLM2 = this%INDF(IG2),this%INDF(IG2+1)-1
            this%INDFFY(ILOCAL) = this%NFFY
            call this%Find_radial_expansion( IG1, IG2, IFLM1, IFLM2 )
            this%INDFFY(ILOCAL) = this%NFFY - this%INDFFY(ILOCAL)
            ILOCAL = ILOCAL + 1
          ENDDO
        ENDDO
      ENDDO

      CALL DE_ALLOC( this%FFL, 'FFL', MYNAME )
      CALL DE_ALLOC( this%FFQ, 'FFQ', MYNAME )
      CALL DE_ALLOC( this%CFFR, 'CFFR', MYNAME )
      CALL DE_ALLOC( this%IFFR, 'IFFR', MYNAME )

      if (NODE==0) write(*,*)"TIME(rad_ex)",MPI_Wtime()-lo_time
      lo_time = MPI_Wtime()

      call this%swap_radial_exp( )

      if (NODE==0) write(*,*) "TIME(swap_r)", MPI_Wtime() - lo_time
      !call MPI_Finalize( IG )
      !stop

      end subroutine INIT

      subroutine compute_spherical_harmonics( this, IG )
C *******************************************************************
C ************************* INPUT ***********************************
C ************************* OUTPUT **********************************
      implicit none
      class(MATEL) :: this
      integer      :: IG
      ! Local variables
      integer :: L, NILM, JLM

      if ( rcut(IG) .GT. RMAX )
     &   call die('MATEL: NQ too small for required cutoff.')

      L = LCUT( IG )
      NILM = (L+2)**2
      IF (this%NF+NILM>this%MF) THEN
        this%MF=(this%NF+NILM)*EXPAND
        CALL RE_ALLOC( this%F, 0, NQ, 1, this%MF, 'F', MYNAME )
        CALL RE_ALLOC( this%ILM, 1, this%MF, 'ILM', MYNAME )
      ENDIF

      CALL YLMEXP( L, RLYLM, EVALUATE, IG, 0, NQ, RMAX, NILM,
     &             this%ILM(this%NF+1:), this%F(:,this%NF+1:) )
      ! Store orbital in k-space
      DO JLM = 1, NILM
        this%NF = this%NF + 1
        L = LOFILM( this%ILM(this%NF) )
        CALL RADFFT( L, NQ, RMAX, this%F(0:NQ,this%NF),
     &    this%F(0:NQ,this%NF) )
        ! F(NQ,NF) = 0._dp
      ENDDO
      end subroutine compute_spherical_harmonics

      subroutine Find_radial_expansion( this, IG1, IG2, IFLM1, IFLM2 )
C *******************************************************************
C ************************* INPUT ***********************************
C ************************* OUTPUT **********************************
      use interpolation, only: spline
      implicit none
      class(MATEL)      :: this
      integer           :: IG1, IG2, IFLM1, IFLM2
      ! Local variables
      integer           :: IQ, NILM, I, L1, L2, L3, IR, JR, JG, JLM
      real(dp)          :: C, R, DFFR0, DFFRMX, CPROP, Q
      real(dp), pointer :: FQ1(:), FQ2(:)
      logical, external :: propor

      ! Check interaction range
      IF (RCUT(IG1)+RCUT(IG2) .GT. RMAX) THEN
        call die('MATEL: NQ too small for required cutoff.')
      ENDIF

      FQ1(0:) => this%F(:,IFLM1)
      FQ2(0:) => this%F(:,IFLM2)

      ! Find orbitals convolution by multiplication in k-space
      C = ( 2.0_dp * PI )**1.5_dp
      DO IQ = 0, NQ
        this%FFQ(IQ) = C * FQ1(IQ) * FQ2(IQ)
        IF (this%MODE==2) THEN
          Q = IQ * DQ
          this%FFQ(IQ) = this%FFQ(IQ) * Q*Q
        ENDIF
      ENDDO

      ! Loop on possible values of l quantum number of product
      L1 = LOFILM( this%ILM(IFLM1) )
      L2 = LOFILM( this%ILM(IFLM2) )
      if (L1+L2 .GT. this%L1L2) then
        this%L1L2 = L1+L2
        CALL RE_ALLOC( this%IFFR, 0, this%L1L2, 'IFFR', MYNAME )
        CALL RE_ALLOC( this%CFFR, 0, this%L1L2, 'CFFR', MYNAME )
      endif
      DO L3 = ABS(L1-L2), L1+L2, 2
        ! Return to real space
        CALL RADFFT( L3, NQ, NQ*PI/RMAX, this%FFQ, this%FFL )
        ! FFL(NQ) = 0._dp
        IF (MOD(ABS(L1-L2-L3)/2,2) .NE. 0) THEN
          DO IR = 0,NR
            this%FFL(IR) = - this%FFL(IR)
          ENDDO
        ENDIF
        ! Divide by R**L
        IF (L3 .NE. 0) THEN
          DO IR = 1, NR
            R = IR * DR
            this%FFL(IR) = this%FFL(IR) / R**L3
          ENDDO
          ! Parabolic extrapolation to R=0
          this%FFL(0) = ( 4.0_dp * this%FFL(1) - this%FFL(2) ) / 3.0_dp
        ENDIF
        ! Select NRTAB out of NR points
        IF (MOD(NR,NRTAB) .NE. 0)
     &    CALL DIE('matel ERROR: NQ must be multiple of NRTAB')
        DO IR = 0, NRTAB
          JR = IR * NR / NRTAB
          this%FFL(IR) = this%FFL(JR)
        ENDDO
        this%IFFR(L3) = 0
        DO JG= 1, this%NFFR
          IF ( PROPOR( NRTAB, this%FFL(1), this%FFR(1,1,JG),
     &                 FFTOL, CPROP)) THEN
            this%IFFR(L3) = JG
            this%CFFR(L3) = CPROP
            EXIT
          ENDIF
        ENDDO
        ! If not found, store new radial function
        if (this%IFFR(L3)==0) then
          this%NFFR = this%NFFR + 1
          IF (this%NFFR .GT. this%MFFR) THEN
            this%MFFR = EXPAND * this%NFFR
            CALL RE_ALLOC( this%FFR, 0, NRTAB, 1, 2, 1, this%MFFR,
     &                     'FFR', MYNAME )
          ENDIF
          this%IFFR(L3) = this%NFFR
          this%CFFR(L3) = 1._dp
          DO IR = 0,NRTAB
            this%FFR(IR,1,this%NFFR) = this%FFL(IR)
          ENDDO
          ! Setup spline interpolation
          ! Force derivative, rather than second derivative, to zero
          ! DFFR0 = HUGE(1.0_dp)
          DFFR0 = 0.0_dp
          DFFRMX = 0.0_dp
          CALL SPLINE( RMAX/NRTAB, this%FFR(0:NRTAB,1,this%NFFR),
     &                 NRTAB+1, DFFR0, DFFRMX,
     &                 this%FFR(0:NRTAB,2,this%NFFR) )
        ENDIF
      ENDDO

      ! Reallocate some arrays
      NILM = (L1+L2+1)**2
      IF (NILM.GT.this%MILM) THEN
        CALL RE_ALLOC( this%Y, 1, NILM, 'Y', MYNAME, .FALSE. )
        CALL RE_ALLOC( this%DYDR, 1, 3, 1, NILM, 'DYDR', MYNAME,.FALSE.)
        this%MILM = NILM
      ENDIF
      IF (this%NFFY+NILM .GT. this%MFFY) THEN
        this%MFFY = EXPAND * (this%NFFY+NILM)
        CALL RE_ALLOC( this%FFY, 1, 1, 1, this%MFFY, 'FFY', MYNAME )
        CALL RE_ALLOC( this%ILMFF, 1, this%MFFY, 'ILMFF', MYNAME )
        CALL RE_ALLOC( this%INDFFR, 1, this%MFFY, 'INDFFR', MYNAME )
      ENDIF
      ! Expand the product of two spherical harmonics (SH) also in SH
      CALL YLMEXP( L1+L2, RLYLM, YLMYLM, this%ILM(IFLM1),
     &             this%ILM(IFLM2), 1, 1, 1.0_dp, NILM,
     &             this%ILMFF(this%NFFY+1:), this%FFY(:,this%NFFY+1:))
C     Loop on possible lm values of orbital product
      DO I = 1,NILM
        this%NFFY = this%NFFY + 1
        JLM = this%ILMFF(this%NFFY)
        L3 = LOFILM( JLM )
        this%INDFFR(this%NFFY) = this%IFFR(L3)
        this%FFY(1,this%NFFY) = this%FFY(1,this%NFFY) * this%CFFR(L3)
      ENDDO
      end subroutine Find_radial_expansion

      subroutine swap_spherical_harmonics( this )
C *******************************************************************
C ************************* INPUT ***********************************
C ************************* OUTPUT **********************************
      implicit none
      class(MATEL)     :: this
      integer          :: i, prev, next
#ifdef MPI
      integer           :: fi, la, MPIERR
      integer,  pointer :: count(:), displ(:), COUN2(:), DISP2(:)
      real(dp), pointer :: G_F(:,:)
      integer,  pointer :: G_ILM(:)

      nullify(count,displ,coun2,disp2)
      CALL RE_ALLOC( COUNT, 1, NODES, 'COUNT,', MYNAME )
      CALL RE_ALLOC( DISPL, 1, NODES+1, 'DISPL,', MYNAME )
      CALL RE_ALLOC( COUN2, 1, NODES, 'COUN2,', MYNAME )
      CALL RE_ALLOC( DISP2, 1, NODES+1, 'DISP2,', MYNAME )

      DISPL(1) = 0
      DO I= 1, NODES
        call GET_LOOP_LIMITS( this%NC, I-1, fi, la )
        DISPL(I+1) = DISPL(I) + la - fi +1
        COUNT(I) = DISPL(I+1) - DISPL(I)
      ENDDO
      ! Reduce INDF
      call MPI_Allgatherv( MPI_IN_PLACE, COUNT(NODE+1), MPI_INTEGER,
     &  this%INDF, COUNT, DISPL, MPI_INTEGER, MPI_COMM_WORLD, MPIERR )
#endif
      prev = this%INDF(1)
      this%INDF(1) = 1
      DO I= 1, this%NC
        next = this%INDF(I+1)
        this%INDF(I+1) = this%INDF(I) + prev
        prev = next
      ENDDO

#ifdef MPI
      DISP2(1) = 0
      this%NF = 0
      DO I= 1, NODES
        this%NF = this%NF+this%INDF(DISPL(I+1)+1)-this%INDF(DISPL(I)+1)
        DISP2(I+1) = this%NF
        COUN2(I) = DISP2(I+1) - DISP2(I)
      ENDDO

      ! Allocate arrays
      nullify(G_ILM)
      CALL RE_ALLOC( G_ILM, 1, this%NF, 'ILM', MYNAME )
      call MPI_Allgatherv( this%ILM, COUN2(NODE+1), MPI_INTEGER,
     &  G_ILM, COUN2, DISP2, MPI_INTEGER, MPI_COMM_WORLD, MPIERR )
      CALL DE_ALLOC( this%ILM, 'ILM', MYNAME )
      this%ILM => G_ILM

      nullify(G_F)
      CALL RE_ALLOC( G_F, 0, NQ, 1, this%NF, 'F', MYNAME )
      COUN2 = COUN2*(NQ+1)
      DISP2 = DISP2*(NQ+1)
      call MPI_Allgatherv( this%F, COUN2(NODE+1), MPI_DOUBLE_PRECISION,
     &  G_F, COUN2, DISP2, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD,
     &  MPIERR )
      CALL DE_ALLOC( this%F, 'F', MYNAME )
      this%F => G_F

      CALL DE_ALLOC( COUNT, 'COUNT', MYNAME )
      CALL DE_ALLOC( DISPL, 'DISPL', MYNAME )
      CALL DE_ALLOC( COUNT, 'COUN2', MYNAME )
      CALL DE_ALLOC( DISP2, 'DISP2', MYNAME )
#endif
      end subroutine swap_spherical_harmonics

      subroutine swap_radial_exp( this )
C *******************************************************************
C ************************* INPUT ***********************************
C ************************* OUTPUT **********************************
      implicit none
      class(MATEL) :: this
      ! Local variables
      integer           :: I,PACUM, NACUM
#ifdef MPI
      integer           :: J, FI, LA,MPIERR, IG1, IG2
      integer,  pointer :: COUNT(:), DISPL(:), COUN2(:), DISP2(:),
     &                     G_ILMFF(:), G_INDFFR(:)
      real(dp), pointer :: G_FFR(:,:,:), G_FFY(:,:)
      ! Allocate space for temporal arrays
      nullify(count,displ,coun2,disp2)
      CALL RE_ALLOC( COUNT, 1, NODES, 'COUNT,', MYNAME )
      CALL RE_ALLOC( DISPL, 1, NODES+1, 'DISPL,', MYNAME )
      CALL RE_ALLOC( COUN2, 1, NODES, 'COUN2,', MYNAME )
      CALL RE_ALLOC( DISP2, 1, NODES+1, 'DISP2,', MYNAME )
      DISPL(1) = 0
      DO I= 1, NODES
        call GET_LOOP_LIMITS( this%NC*NUM_OR, I-1, fi, la )
        DISPL(I+1) = this%OrbitalToHarmonic(MOD(LA-1,this%NC)+1,
     &               (LA-1)/this%NC+1)
        COUNT(I) = DISPL(I+1) - DISPL(I)
      ENDDO
      ! Reduce INDFFY
      call MPI_Allgatherv( MPI_IN_PLACE, COUNT(NODE+1), MPI_INTEGER,
     &  this%INDFFY, COUNT, DISPL, MPI_INTEGER, MPI_COMM_WORLD, MPIERR )
#endif
      ! Create the acumulated version of INDFFY
      PACUM = 1
      DO I= 1, this%NF*this%NF_MINOR
        NACUM = PACUM + this%INDFFY(I)
        this%INDFFY(I) = PACUM
        PACUM = NACUM
      ENDDO
      this%INDFFY(this%NF*this%NF_MINOR+1) = PACUM
#ifdef MPI
      DISP2(1) = 0
      DO I= 1, NODES
        COUN2(I) = this%INDFFY(DISPL(I+1)+1)-this%INDFFY(DISPL(I)+1)
        DISP2(I+1) = DISP2(I) + COUN2(I)
      ENDDO

      ! Get the accumulated NFFR
      call MPI_AllGather( this%NFFR, 1, MPI_INTEGER, COUNT,
     &                    1, MPI_INTEGER, MPI_COMM_WORLD, MPIERR )
      DISPL(1) = 0
      DO I= 1, NODES
        DISPL(I+1) = DISPL(I) + COUNT(I)
      ENDDO

      ! Update INDFFR to local indices
      this%INDFFR(1:this%NFFY) = this%INDFFR(1:this%NFFY)+DISPL(NODE+1)
      ! Get the global NFFY
      this%NFFY = DISP2(NODES+1)

      ! Reduce FFY
      nullify(G_FFY)
      CALL RE_ALLOC( G_FFY, 1, 1, 1, this%NFFY, 'FFY', MYNAME )
      call MPI_Allgatherv( this%FFY, COUN2(NODE+1),MPI_DOUBLE_PRECISION,
     &  G_FFY, COUN2, DISP2, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD,
     &  MPIERR )
      CALL DE_ALLOC( this%FFY, 'FFY', MYNAME )
      this%FFY => G_FFY

      ! Reduce ILMFF
      nullify(G_ILMFF)
      CALL RE_ALLOC( G_ILMFF, 1, this%NFFY, 'ILMFF', MYNAME )
      call MPI_Allgatherv( this%ILMFF, COUN2(NODE+1), MPI_INTEGER,
     &  G_ILMFF, COUN2, DISP2, MPI_INTEGER, MPI_COMM_WORLD,
     &  MPIERR )
      CALL DE_ALLOC( this%ILMFF, 'ILMFF', MYNAME )
      this%ILMFF => G_ILMFF


      ! Reduce ILMFF
      nullify(G_INDFFR)
      CALL RE_ALLOC( G_INDFFR, 1, this%NFFY, 'INDFFR', MYNAME )
      call MPI_Allgatherv( this%INDFFR, COUN2(NODE+1), MPI_INTEGER,
     &  G_INDFFR, COUN2, DISP2, MPI_INTEGER, MPI_COMM_WORLD,
     &  MPIERR )
      CALL DE_ALLOC( this%INDFFR, 'INDFFR', MYNAME )
      this%INDFFR => G_INDFFR


      this%NFFR = DISPL(NODES+1)
      COUN2 = COUNT*2*(NRTAB+1)
      DISP2 = DISPL*2*(NRTAB+1)
      nullify(G_FFR)
      CALL RE_ALLOC( G_FFR, 0, NRTAB, 1, 2, 1, this%NFFR,'FFR',MYNAME )
      call MPI_Allgatherv( this%FFR, COUN2(NODE+1),MPI_DOUBLE_PRECISION,
     &  G_FFR, COUN2, DISP2, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD,
     &  MPIERR )
      CALL DE_ALLOC( this%FFR, 'FFR', MYNAME )
      this%FFR => G_FFR

      CALL DE_ALLOC( COUNT, 'COUNT', MYNAME )
      CALL DE_ALLOC( DISPL, 'DISPL', MYNAME )
      CALL DE_ALLOC( COUNT, 'COUN2', MYNAME )
      CALL DE_ALLOC( DISP2, 'DISP2', MYNAME )
#endif
      end subroutine swap_radial_exp

      function OrbitalToHarmonic( this, IG1, IG2 ) result(IG)
C *******************************************************************
C ************************* INPUT ***********************************
C ************************* OUTPUT **********************************
      implicit none
      class(MATEL)        :: this
      integer :: IG1, IG2
      integer :: IG
      IG = this%NF_MINOR*(this%INDF(IG2)-1) +
     &         (this%INDF(IG2+1)-this%INDF(IG2))*(this%INDF(IG1)-1)+1
      end function OrbitalToHarmonic

      subroutine get_matel( THIS, IG1, IG2, R12, S12, DSDR )
C *******************************************************************
C ************************* INPUT ***********************************
C ************************* OUTPUT **********************************
      use interpolation, only: splint
      implicit none
      class(MATEL) :: THIS
      integer      :: IG1, IG2
      real(dp)     :: R12(3), DSDR(3), S12
C     Local Parameter
      real(dp), parameter :: TINY = 1.e-12_dp
C     Local variables
      integer  :: IX, IH, LMAX, IFFY, JLM, JFFR, IFLM1, IFLM2
      real(dp) ::  X12(3), R, SR, DSRDR

      if (IG2.gt.NUM_OR .or. IG1.GT.this%NC) then
        write(*,*) '=== IG1', IG1, 'IG2', IG2
        call die('MATEL: check orbital matrix.')
      endif
C     Initialize output
      S12 = 0.0_dp
      DSDR(1) = 0.0_dp
      DSDR(2) = 0.0_dp
      DSDR(3) = 0.0_dp

C     Avoid R12=0
      X12(1) = R12(1)
      X12(2) = R12(2)
      X12(3) = R12(3)
      R = SQRT( X12(1)*X12(1) + X12(2)*X12(2) + X12(3)*X12(3) )
      IF (R .LT. TINY) THEN
        X12(3) = TINY
        R = SQRT( X12(1)*X12(1) + X12(2)*X12(2) + X12(3)*X12(3) )
      ENDIF

C     Find if orbitals are far (out of range)
      IF (R .LE. RCUT(IG1)+RCUT(IG2) ) THEN
C       Find spherical harmonics times R**L
        IH = this%OrbitalToHarmonic( IG1, IG2 )
        DO IFLM1 = this%INDF(IG1),this%INDF(IG1+1)-1
          DO IFLM2 = this%INDF(IG2), this%INDF(IG2+1)-1
            LMAX = 0
            DO IFFY = this%INDFFY(IH), this%INDFFY(IH+1)-1
              JLM = this%ILMFF(IFFY)
              LMAX = MAX( LMAX, LOFILM(JLM) )
            ENDDO
            CALL RLYLM( LMAX, X12, this%Y, this%DYDR )
C           Interpolate radial functions and obtain SH expansion
            DO IFFY = this%INDFFY(IH), this%INDFFY(IH+1)-1
              JFFR = this%INDFFR(IFFY)
              CALL SPLINT( RMAX/NRTAB, this%FFR(0:NRTAB,1,JFFR),
     &          this%FFR(0:NRTAB,2,JFFR), NRTAB+1, R, SR, DSRDR )
              JLM = this%ILMFF(IFFY)
              S12 = S12 + SR * this%FFY(1,IFFY) * this%Y(JLM)
              DO IX = 1,3
                DSDR(IX) = DSDR(IX) +
     &            DSRDR * this%FFY(1,IFFY) * this%Y(JLM) * X12(IX) / R +
     &            SR * this%FFY(1,IFFY) * this%DYDR(IX,JLM)
              ENDDO
            ENDDO
            IH = IH + 1
          ENDDO
        ENDDO
      ENDIF
      end subroutine get_matel

      subroutine GET_LOOP_LIMITS( IND, ID,  FIRST, LAST )
C *******************************************************************
C ************************* INPUT ***********************************
C ************************* OUTPUT **********************************
      implicit none
      integer :: IND, ID, FIRST, LAST
      integer :: di, mo
      if ( NODES .eq. 1 ) THEN
        FIRST = 1
        LAST  = IND
      else
        di   = IND / NODES
        mo   = MOD(IND,NODES)
        FIRST = 1 + ID*di + MIN(ID,mo)
        LAST  = (ID+1)*di + MIN(ID+1,mo)
      endif
      end subroutine GET_LOOP_LIMITS

      subroutine init_matel( )
C *******************************************************************
C ************************* INPUT ***********************************
C ************************* OUTPUT **********************************
      use atm_types, only : nspecies
      use atmfuncs,  only : nofis, nkbfis, orb_gindex, kbproj_gindex
      implicit none
      ! Local Variables
      integer :: IO, KO
      TYPE(allocDefaults) :: OLDEFS
      ! Set allocation defaults
      call alloc_default( old=oldefs, copy=.true., shrink=.false. )

      CALL RESET_SPHER_HARM( )
      CALL RESET_RADFFT( )

      ! Get the number of orbitals and KB orbitals
      IO = NOFIS( nspecies )
      KO = NKBFIS( nspecies )
      NUM_OR = orb_gindex( nspecies, IO )
      NO     = kbproj_gindex( nspecies, -KO )
      num_kb = NO - NUM_OR
      if (NODE==0) write(*,*) '=== init_matel ===', NUM_OR, num_kb, NO
      ! Init constant variables
      PI = 4._dp * ATAN(1._dp)
      QMAX = 2._dp * SQRT( Q2CUT )
      DQ = QMAX / NQ
      DR = PI / QMAX
      RMAX = NR * DR
      DRTAB = RMAX / NRTAB

      ! Allocate S Matel
      call S%init( 1 )
      call T%init( 2 )

      CALL alloc_default( restore=oldefs )

      end subroutine init_matel

      subroutine get_matel_s( IG1, IG2, R12, S12, DSDR )
C *******************************************************************
C ************************* INPUT ***********************************
C ************************* OUTPUT **********************************
      use interpolation, only: splint
      implicit none
      integer  :: IG1, IG2
      real(dp) ::  R12(3), DSDR(3), S12
      ! Local Variables
      call S%GET_MATEL( IG1, IG2, R12, S12, DSDR )
      end subroutine get_matel_s

      subroutine get_matel_t( IG1, IG2, R12, S12, DSDR )
C *******************************************************************
C ************************* INPUT ***********************************
C ************************* OUTPUT **********************************
      use interpolation, only: splint
      implicit none
      integer  :: IG1, IG2
      real(dp) ::  R12(3), DSDR(3), S12
      ! Local Variables
      call T%GET_MATEL( IG1, IG2, R12, S12, DSDR )
      end subroutine get_matel_t


      end module matel_mod
