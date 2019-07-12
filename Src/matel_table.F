!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
!
      !> Re-designed module to allow a pre-computation of the needed matrix elements
      !> in parallel, followed by a globalization of the data among all the MPI processes.
      !> Once the interpolation tables are setup, further calls to the matrix-element
      !> evaluator (here renamed 'get_matel') are cheap. This has a dramatic effect in
      !> some routines (such as nlefsm) that had to perform the table-building operations
      !> under conditions that did not scale in parallel.

      !> Concept: Rogeli Grima (BSC) and Alberto Garcia (ICMAB)
      !> Initial implementation:  Rogeli Grima (BSC)

      module matel_table

      use matel_params
      use matel_ylm, only: spher_harm_t
      
      use precision, only : dp
      use alloc,     only : re_alloc, de_alloc, alloc_default,
     &                      allocDefaults

      use parallel,  only : Node, Nodes
      use m_radfft,  only : radfft
      use m_matel_registry, ONLY : RCUT
      use spher_harm, only : RLYLM, YLMEXP, ylmylm, lofilm
#ifdef MPI
      use mpi_siesta,  only : MPI_INTEGER, MPI_COMM_WORLD
      use mpi_siesta,  only : MPI_DOUBLE_PRECISION, MPI_IN_PLACE
#endif      
      private

      integer, parameter, public :: MODE_S = 1, MODE_T = 2, MODE_XYZ = 3

      !> The main type that holds the interpolation tables and all the indexes needed to find
      !> the right radial tables. Note that it contains pointers to "spherical-harmonic decompositions",
      !> with two versions (1 and 2), but '2' is not always needed.
      TYPE, public :: MATEL
        integer           :: MODE             ! One of MODE_S, MODE_T or MODE_XYZ

        integer           :: min_ig1              ! Registry index of first "F_1" function
        integer           :: max_ig1              ! Registry index of last "F_1" function
        integer           :: min_ig2              ! Registry index of first "F_2" function
        integer           :: max_ig2              ! Registry index of last "F_2" function

        integer           :: MAT_ROWS         ! Leading size of matrix (in num of harmonics)
        integer           :: MAT_COLS         ! Columns of matrix (in num of harmonics)

        TYPE(SPHER_HARM_t), pointer :: SpHa1  ! Orbitals & projectors
        TYPE(SPHER_HARM_t), pointer :: SpHa2
        integer           :: BASE_SPHA1        ! Base index for block of interest in SpHa1
        integer           :: BASE_SPHA2        ! Base index for block of interest in SpHa2

        real(dp), POINTER :: FFR(:,:,:)       ! Radial functions. To save space we only store
        integer           :: MFFR, NFFR       ! linear independent functions.
        integer,  POINTER :: INDFFR(:)        !

        real(dp), POINTER :: FFY(:,:)         ! Expand of the product of two spherical harmonics
        integer,  POINTER :: ILMFF(:)         ! Spherical harmonics indexes
        integer           :: MFFY, NFFY       !
        integer,  POINTER :: INDFFY(:)        ! Number of radial functions (after COMPUTE_RADEXP)
                                              ! Index of radial functions (after REDUCE_RADEXP)

        integer           :: MILM
        !> Work space for computation of rlylm
        real(dp), POINTER :: Y(:), DYDR(:,:)

        CONTAINS
        procedure :: INIT
        procedure :: GET_MATEL
        ! called only internally
        procedure :: COMPUTE_RADEXP
        procedure :: REDUCE_RADEXP

      END TYPE MATEL


      real(dp),         parameter :: EXPAND    =  1.20_dp
      integer,          parameter :: MINEXPAND =  32
      real(dp),         parameter :: FFTOL     =  1.e-8_dp
      CHARACTER(LEN=*), parameter :: MYNAME =  'MATEL_TABLE'

      external :: die

      CONTAINS
      
      subroutine INIT( this, IOPER, SpHa1, lo1, hi1,
     $                              SpHa2, lo2, hi2)

!     Initialize a MATEL table structure using the previously calculated expansions
!     in spherical harmonics of the fourier transforms of the functions.

!      M(R) =  < F_1(r) | OP | F_2(r-R) >
!      M(R) =  < F_1(r) | G_2(r-R) >

!           G_2(r) is F_2(r) when OP is the identity (simple overlap)
!           G_2(r) can be  {x,y,z}*F_2  when OP = {x,y,z}
!           When OP is the laplacian, this routine will just insert a factor of Q^2
!           in the product of F_1(Q) and F_2(Q).

      implicit none

      class(MATEL)        :: this

      !> One of MODE_S, MODE_T or MODE_XYZ
      integer, intent(in) :: IOPER

      !> Structures with the expansions of F_1 and F_2 in spherical harmonics
      type(SPHER_HARM_t), intent(in), target :: SpHa1, SpHa2

      !> First and last entries to consider in SpHa1
      integer, intent(in) :: lo1, hi1
      !> First and last entries to consider in SpHa2
      integer, intent(in) :: lo2, hi2

      integer :: n_harms_hi, n_harms_lo
      
      this%MODE  = IOPER

      this%min_ig1 = lo1 + spha1%base
      this%max_ig1 = hi1 + spha1%base
      this%min_ig2 = lo2 + spha2%base
      this%max_ig2 = hi2 + spha2%base

      this%SpHa1 => SpHa1
      this%SpHa2 => SpHa2

      ! Find total number of harmonics in the block of SpHa1 we are
      ! interested in

      n_harms_hi = this%SpHa1%get_Dim( hi1 )
      n_harms_lo = 0
      if (lo1 > 1 ) then
         n_harms_lo = this%SpHa1%get_Dim( lo1 - 1 )
      endif
      this%MAT_ROWS = n_harms_hi - n_harms_lo
      this%BASE_SPHA1 = n_harms_lo

      ! Find total number of harmonics in the block of SpHa2 we are
      ! interested in

      n_harms_hi = this%SpHa2%get_Dim( hi2 )
      n_harms_lo = 0
      if (lo2 > 1 ) then
         n_harms_lo = this%SpHa2%get_Dim( lo2 - 1 )
      endif
      this%MAT_COLS = n_harms_hi - n_harms_lo
      this%BASE_SPHA2 = n_harms_lo

      ! Compute and reduce Radial Expansion
      call this%compute_RadExp( )
      call this%reduce_RadExp( )

      end subroutine INIT

      subroutine compute_RadExp( this )
C *******************************************************************
C  Fill a matrix with the radial expansion functions
C ************************* DESCRI **********************************
C  Compute the radial expansion between all harmonics in parallel.
C  For very harmonic permutation, we compute its Radial functions and
C  the expansion of the product of spherical harmonics
C *******************************************************************
      use interpolation, only: spline
      implicit none
      class(MATEL)      :: this
      integer           :: dim1, dim2, nlocal, first, last, i, j,
     &                     IR, JR, JG, ih1, ih2, ig1, ig2, IQ,
     &                     L1, L2, L3, L1L2, JLM, NILM
      real(dp)          :: C, Q, R, DFFR0, DFFRMX, CPROP
      integer,  pointer :: IFFR(:)
      real(dp), pointer :: FQ1(:), FQ2(:), FFQ(:), FFL(:), CFFR(:)
      logical, external :: propor

      ! Initialize temporal arrays
      nullify(this%Y,this%DYDR)
      this%MILM = (MINEXPAND*MINEXPAND+1)**2
      CALL RE_ALLOC( this%Y, 1, this%MILM, 'Y', MYNAME )
      CALL RE_ALLOC( this%DYDR, 1, 3, 1, this%MILM, 'DYDR', MYNAME )

      dim1 = this%MAT_ROWS
      dim2 = this%MAT_COLS
      nullify( this%INDFFY )
      call RE_ALLOC( this%INDFFY, 1, dim1*dim2+1, 'INDFFY', MYNAME )

      nlocal = GET_LOOP_LIMITS( dim1*dim2, Node, first, last )

      nullify(this%FFR,this%INDFFR)
      this%MFFR = MAX( INT(dim1*EXPAND), dim1+MINEXPAND )
      CALL RE_ALLOC( this%FFR, 0, NRTAB, 1, 2, 1,
     &               this%MFFR, 'FFR', MYNAME )
      this%NFFR = 0

      nullify(this%FFY,this%ILMFF,this%INDFFR)
      this%MFFY = MAX( INT(nlocal*EXPAND), nlocal+MINEXPAND )
      CALL RE_ALLOC( this%FFY, 1, 1, 1, this%MFFY, 'FFY', MYNAME )
      CALL RE_ALLOC( this%INDFFR, 1, this%MFFY, 'INDFFR', MYNAME )
      CALL RE_ALLOC( this%ILMFF, 1, this%MFFY, 'ILMFF', MYNAME )
      this%NFFY = 0

      nullify(IFFR,CFFR,FFQ,FFL)
      CALL RE_ALLOC( FFQ, 0, NQ, 'FFQ', MYNAME )
      CALL RE_ALLOC( FFL, 0, NQ, 'FFL', MYNAME )
      L1L2 = MINEXPAND
      CALL RE_ALLOC( IFFR, 0, L1L2, 'IFFR', MYNAME )
      CALL RE_ALLOC( CFFR, 0, L1L2, 'CFFR', MYNAME )

      ! Iterate all combinations of harmonics
      do i= first, last
        ! Get the harmonic terms IDS
        ! Note accounting for block of interest
        ih1 = MOD(i-1,dim1) + 1 + this%BASE_SPHA1
        ih2 = (i-1)/dim1 + 1 + this%BASE_SPHA2
        ! Get the matel-registry IDS of the functions
        ! from which the harmonic terms originate
        ig1 = this%SpHa1%Harm2Orb( ih1 )
        ig2 = this%SpHa2%Harm2Orb( ih2 )

        ! Check interaction range
        IF (RCUT(IG1)+RCUT(IG2) .GT. RMAX) THEN
          call die('MATEL: NQ too small for required cutoff.')
        ENDIF
        FQ1(0:) => this%SpHa1%F(:,ih1)
        FQ2(0:) => this%SpHa2%F(:,ih2)

        ! Find orbitals convolution by multiplication in k-space
        C = ( 2.0_dp * PI )**1.5_dp
        DO IQ = 0, NQ
          FFQ(IQ) = C * FQ1(IQ) * FQ2(IQ)
          IF (this%MODE==MODE_T) THEN
            Q = IQ * DQ
            FFQ(IQ) = FFQ(IQ) * Q*Q
          ENDIF
        ENDDO
        ! Loop on possible values of l quantum number of product
        L1 = LOFILM( this%SpHa1%ILM(IH1) )
        L2 = LOFILM( this%SpHa2%ILM(IH2) )
        if (L1+L2 .GT. L1L2) then
          ! Reallocate more memory if necessary
          L1L2 = L1+L2
          CALL RE_ALLOC( IFFR, 0, L1L2, 'IFFR', MYNAME )
          CALL RE_ALLOC( CFFR, 0, L1L2, 'CFFR', MYNAME )
        endif
        DO L3 = ABS(L1-L2), L1+L2, 2
          ! Return to real space
          CALL RADFFT( L3, NQ, NQ*PI/RMAX, FFQ, FFL )
          ! FFL(NQ) = 0._dp
          IF (MOD(ABS(L1-L2-L3)/2,2) .NE. 0) THEN
            DO IR = 0, NR
              FFL(IR) = - FFL(IR)
            ENDDO
          ENDIF
          ! Divide by R**L
          IF (L3 .NE. 0) THEN
            DO IR = 1, NR
              R = IR * DR
              FFL(IR) = FFL(IR) / R**L3
            ENDDO
            ! Parabolic extrapolation to R=0
            FFL(0) = ( 4.0_dp * FFL(1) - FFL(2) ) / 3.0_dp
          ENDIF
          ! Select NRTAB out of NR points
          IF (MOD(NR,NRTAB) .NE. 0)
     &      CALL DIE('matel ERROR: NQ must be multiple of NRTAB')
          DO IR = 0, NRTAB
            JR = IR * NR / NRTAB
            FFL(IR) = FFL(JR)
          ENDDO

          ! To save space, we check if the new radial function has
          ! been seen before (save a global factor of proportionality)
          IFFR(L3) = 0
          DO JG= 1, this%NFFR
            IF ( PROPOR( NRTAB, FFL(1), this%FFR(1,1,JG),
     &                   FFTOL, CPROP)) THEN
              ! If found, save the index and the coefficient of proportionality
              IFFR(L3) = JG
              CFFR(L3) = CPROP
              EXIT
            ENDIF
          ENDDO
          ! If not found, store new radial function
          if (IFFR(L3)==0) then
            this%NFFR = this%NFFR + 1
            IF (this%NFFR .GT. this%MFFR) THEN
              this%MFFR = EXPAND * this%NFFR
              CALL RE_ALLOC( this%FFR, 0, NRTAB, 1, 2, 1, this%MFFR,
     &                       'FFR', MYNAME )
            ENDIF
            IFFR(L3) = this%NFFR
            CFFR(L3) = 1._dp
            DO IR = 0,NRTAB
              this%FFR(IR,1,this%NFFR) = FFL(IR)
            ENDDO
            ! Setup spline interpolation
            ! Force derivative, rather than second derivative, to zero
            ! DFFR0 = HUGE(1.0_dp)
            DFFR0 = 0.0_dp
            DFFRMX = 0.0_dp
            CALL SPLINE( RMAX/NRTAB, this%FFR(0:NRTAB,1,this%NFFR),
     &                   NRTAB+1, DFFR0, DFFRMX,
     &                   this%FFR(0:NRTAB,2,this%NFFR) )
          ENDIF
        ENDDO

        ! Reallocate some arrays for the angular expansion
        NILM = (L1+L2+1)**2
        IF (NILM.GT.this%MILM) THEN
          CALL RE_ALLOC( this%Y, 1, NILM, 'Y', MYNAME, .FALSE. )
          CALL RE_ALLOC( this%DYDR, 1, 3, 1, NILM, 'DYDR',
     &                   MYNAME, .FALSE. )
          this%MILM = NILM
        ENDIF
        IF (this%NFFY+NILM .GT. this%MFFY) THEN
          this%MFFY = EXPAND * (this%NFFY+NILM)
          CALL RE_ALLOC( this%FFY, 1, 1, 1, this%MFFY, 'FFY', MYNAME )
          CALL RE_ALLOC( this%ILMFF, 1, this%MFFY, 'ILMFF', MYNAME )
          CALL RE_ALLOC( this%INDFFR, 1, this%MFFY, 'INDFFR', MYNAME )
        ENDIF
        ! Expand the product of two spherical harmonics (SH) also in SH
        CALL YLMEXP( L1+L2, RLYLM, YLMYLM, this%SpHa1%ILM(IH1),
     &               this%SpHa2%ILM(IH2), 1, 1, 1.0_dp, NILM,
     &               this%ILMFF(this%NFFY+1:), this%FFY(:,this%NFFY+1:))
C       Loop on possible lm values of orbital product
        DO J = 1, NILM
          this%NFFY = this%NFFY + 1
          JLM = this%ILMFF(this%NFFY)
          L3 = LOFILM( JLM )
          this%INDFFR(this%NFFY) = IFFR(L3)
          this%FFY(1,this%NFFY) = this%FFY(1,this%NFFY) * CFFR(L3)
        ENDDO
        this%INDFFY(I) = NILM
      enddo

      CALL DE_ALLOC( CFFR, 'CFFR', MYNAME )
      CALL DE_ALLOC( IFFR, 'IFFR', MYNAME )
      CALL DE_ALLOC( FFQ, 'FFQ', MYNAME )
      CALL DE_ALLOC( FFL, 'FFL', MYNAME )
      end subroutine compute_RadExp

      subroutine reduce_RadExp( this )
C *******************************************************************
C  Reduce Radial functions and expansion of spherical harmonics
C ************************* DESCRI **********************************
C  Reduce the values of the different arrays of MATEL
C  At input, INDFFY contains the number of harmonics of every
C  interaction. At exit, it can be used as and index.
C  NOTE: TO DO!!! The reduction of FFR can be improved. We are saving
C  functions that are not linearly independent (since we did the
C  work on different processors). We should choose between saving time
C  or saving space.
C *******************************************************************
      implicit none
      class(MATEL) :: this
      ! Local variables
      integer           :: I, PREV, NACUM, DIM1, DIM2
#ifdef MPI
      integer           :: J, FIRST, LAST, MPIERR, IG1, IG2
      integer,  pointer :: COUNT(:), DISPL(:), COUN2(:), DISP2(:),
     &                     G_ILMFF(:), G_INDFFR(:)
      real(dp), pointer :: G_FFR(:,:,:), G_FFY(:,:)
#endif
      dim1 = this%MAT_ROWS
      dim2 = this%MAT_COLS
#ifdef MPI
      if (NODES.GT.1) then
        ! Allocate space for temporal arrays
        nullify(count,displ,coun2,disp2)
        CALL RE_ALLOC( COUNT, 0, NODES-1, 'COUNT,', MYNAME )
        CALL RE_ALLOC( DISPL, 0, NODES, 'DISPL,', MYNAME )
        CALL RE_ALLOC( COUN2, 0, NODES-1, 'COUN2,', MYNAME )
        CALL RE_ALLOC( DISP2, 0, NODES, 'DISP2,', MYNAME )
        DISPL(1) = 0
        DO I= 0, NODES-1
          COUNT(I) = GET_LOOP_LIMITS( dim1*dim2, I, first, last )
          DISPL(I+1) = DISPL(I) + COUNT(I)
        ENDDO
        ! Reduce INDFFY
        call MPI_Allgatherv( MPI_IN_PLACE, COUNT(NODE), MPI_INTEGER,
     &                       this%INDFFY, COUNT, DISPL, MPI_INTEGER,
     &                       MPI_COMM_WORLD, MPIERR )
      endif
#endif
      ! Create the acumulated version of INDFFY
      PREV           = this%INDFFY(1)
      this%INDFFY(1) = 1
      DO I= 1, dim1*dim2
        NACUM            = PREV + this%INDFFY(I)
        PREV             = this%INDFFY(I+1)
        this%INDFFY(I+1) = NACUM
      ENDDO
#ifdef MPI
      if (NODES.GT.1) then
        DISP2(0) = 0
        DO I= 0, NODES-1
          COUN2(I) = this%INDFFY(DISPL(I+1)+1)-this%INDFFY(DISPL(I)+1)
          DISP2(I+1) = DISP2(I) + COUN2(I)
        ENDDO

        ! Get the accumulated NFFR
        call MPI_AllGather( this%NFFR, 1, MPI_INTEGER, COUNT,
     &                      1, MPI_INTEGER, MPI_COMM_WORLD, MPIERR )
        DISPL(0) = 0
        DO I= 0, NODES-1
          DISPL(I+1) = DISPL(I) + COUNT(I)
        ENDDO

        ! Update INDFFR to local indices
        this%INDFFR(1:this%NFFY) = this%INDFFR(1:this%NFFY)+DISPL(NODE)
        ! Get the global NFFY
        this%NFFY = DISP2(NODES)

        ! Reduce FFY
        nullify(G_FFY)
        CALL RE_ALLOC( G_FFY, 1, 1, 1, this%NFFY, 'FFY', MYNAME )
        call MPI_Allgatherv( this%FFY, COUN2(NODE),MPI_DOUBLE_PRECISION,
     &    G_FFY, COUN2, DISP2, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD,
     &    MPIERR )
        CALL DE_ALLOC( this%FFY, 'FFY', MYNAME )
        this%FFY => G_FFY

        ! Reduce ILMFF
        nullify(G_ILMFF)
        CALL RE_ALLOC( G_ILMFF, 1, this%NFFY, 'ILMFF', MYNAME )
        call MPI_Allgatherv( this%ILMFF, COUN2(NODE), MPI_INTEGER,
     &    G_ILMFF, COUN2, DISP2, MPI_INTEGER, MPI_COMM_WORLD,
     &    MPIERR )
        CALL DE_ALLOC( this%ILMFF, 'ILMFF', MYNAME )
        this%ILMFF => G_ILMFF


        ! Reduce ILMFF
        nullify(G_INDFFR)
        CALL RE_ALLOC( G_INDFFR, 1, this%NFFY, 'INDFFR', MYNAME )
        call MPI_Allgatherv( this%INDFFR, COUN2(NODE), MPI_INTEGER,
     &    G_INDFFR, COUN2, DISP2, MPI_INTEGER, MPI_COMM_WORLD,
     &    MPIERR )
        CALL DE_ALLOC( this%INDFFR, 'INDFFR', MYNAME )
        this%INDFFR => G_INDFFR


        this%NFFR = DISPL(NODES)
        COUN2 = COUNT*2*(NRTAB+1)
        DISP2 = DISPL*2*(NRTAB+1)
        nullify(G_FFR)
        CALL RE_ALLOC( G_FFR, 0, NRTAB, 1, 2, 1, this%NFFR,'FFR',MYNAME)
        call MPI_Allgatherv( this%FFR, COUN2(NODE),MPI_DOUBLE_PRECISION,
     &    G_FFR, COUN2, DISP2, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD,
     &    MPIERR )
        CALL DE_ALLOC( this%FFR, 'FFR', MYNAME )
        this%FFR => G_FFR

        CALL DE_ALLOC( COUNT, 'COUNT', MYNAME )
        CALL DE_ALLOC( DISPL, 'DISPL', MYNAME )
        CALL DE_ALLOC( COUNT, 'COUN2', MYNAME )
        CALL DE_ALLOC( DISP2, 'DISP2', MYNAME )
      endif
#endif
      end subroutine reduce_RadExp


      !> Utility function      
      function GET_LOOP_LIMITS( N, ID,  FIRST, LAST ) RESULT(SIZE)
C *******************************************************************
C  Get the loop limits adequate to split N elements among Nodes
C ************************* INPUT ***********************************
C  integer  N : Number of elements to distribute
C  integer ID : Id of the Node
C ************************* OUTPUT **********************************
C  integer  FIRST : first element of the loop
C  integer   LAST : last element of the loop
C  integer   SIZE : Number of elements of the loop
C *******************************************************************
      implicit none
      integer,  intent(in) :: N, ID
      integer, intent(out) :: FIRST, LAST
      integer              :: size, di, mo
      if ( NODES .eq. 1 ) THEN
        FIRST = 1
        LAST  = N
      else
        di   = N / NODES
        mo   = MOD(N,NODES)
        FIRST = 1 + ID*di + MIN(ID,mo)
        LAST  = (ID+1)*di + MIN(ID+1,mo)
      endif
      SIZE = LAST - FIRST + 1
      end function GET_LOOP_LIMITS

      !>  Get the matrix element and its derivative
      subroutine get_matel( THIS, IG1, IG2, R12, S12, DSDR )

      use interpolation, only: splint

      implicit none

      !> Matel_table object
      class(MATEL) :: THIS

      !> Global index of 1st function
      integer, intent(in)      :: IG1
      !> Global index of 2nd function
      integer, intent(in)      :: IG2

      !> Vector from first to second atom
      real(dp), intent(in)    :: R12(3)

      !> Matrix element
      real(dp), intent(out)    :: S12
      !> Derivative (gradient) of S12 with respect to R12.
      real(dp), intent(out)    :: DSDR(3)

      real(dp), parameter :: TINY = 1.e-12_dp

      integer  :: IX, IH, LMAX, IFFY, JLM, JFFR, IFLM1, IFLM2
      real(dp) ::  X12(3), R, SR, DSRDR

      if ( (ig1 < this%min_ig1) .or. (ig1 > this%max_ig1) ) then
         call die("IG1 out of range in get_matel")
      endif
      if ( (ig2 < this%min_ig2) .or. (ig2 > this%max_ig2) ) then
         call die("IG2 out of range in get_matel")
      endif

!     Initialize output
      S12 = 0.0_dp
      DSDR(1) = 0.0_dp
      DSDR(2) = 0.0_dp
      DSDR(3) = 0.0_dp

!     Avoid R12=0
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
        DO IFLM1= this%SpHa1%INDF(IG1),this%SpHa1%INDF(IG1+1)-1
          DO IFLM2= this%SpHa2%INDF(IG2), this%SpHa2%INDF(IG2+1)-1
            ! Note
            IH = (IFLM1-this%BASE_SPHA1)  +
     $           (IFLM2-this%BASE_SPHA2-1)*this%MAT_ROWS
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
          ENDDO
        ENDDO
      ENDIF
      end subroutine get_matel

      end module matel_table
