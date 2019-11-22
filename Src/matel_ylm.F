!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
!
!> This module implements a data structure and associated procedures to
!> handle the expansion in spherical harmonics of the functions in the
!> 'matel registry', and their fourier transforms, as a first step in
!> the computation of two-center integrals (See Siesta paper). In this
!> context, several "functions" might be associated to the same object
!> in the registry: for example, a given PAO with index 'ig' might be
!> just 'evaluate'd, or 'x_evaluate'd, etc.  Note that *most* of the
!> functions we deal with are "pure" (PAOs, KBs, and LDAU projectors are
!> of the form f(r)Y_lm(theta,phi) already, so most of the dereferencing
!> done here is superfluous.  We have, on the other hand,
!> {x,y,z}*{PAO,KB} (but they might be simply Y_(l+1),m for a given m?)
!> and the wannier projectors, which might not be "pure".

!> The computations are done in parallel and the results replicated
!> ("all-reduced") in in all processors.
      
!> Initial implementation: Rogeli Grima (BSC)
!> Clarifications and some enhancements: Alberto Garcia (ICMAB)

      module matel_ylm

      use matel_params, only: NQ, RMAX
      
      use precision, only : dp
      use alloc,     only : re_alloc, de_alloc, alloc_default,
     &                      allocDefaults
      use sys,       only: die
      use parallel,  only : Node, Nodes
      use m_radfft,  only : radfft, reset_radfft

      use m_matel_registry, only: LCUT
      
      use spher_harm, only : RLYLM, YLMEXP, lofilm
      private

C     parameters
      real(dp),         parameter :: EXPAND    =  1.20_dp
      integer,          parameter :: MINEXPAND =  32
      CHARACTER(LEN=*), parameter :: MYNAME =  'MATEL_YLM'


      !> This really holds the expansion of a set of functions in spherical-harmonics 
      !> Indexing arrays could be avoided if we kept the information of a given function
      !> in a simpler derived type instead, and making this an _array_ of the simpler elements.
   
      TYPE, public :: SPHER_HARM_t
        integer           :: N           ! Number of orbitals/projectors
        integer           :: BASE        ! Base for registry index count
        real(dp), POINTER :: F(:,:)      ! Radial expansion for each spherical harmonic
        integer,  POINTER :: ILM(:)      ! Spherical harmonics indexes (packed l,m)
        integer,  POINTER :: IND(:)      ! Number of spherical harmonics (after COMPUTE_SPHA)
                                         ! Index of harmonics (after REDUCE_SPHA)
        integer,  POINTER :: INVP(:)     ! Map from harmonic components to the associated function

        CONTAINS
          procedure :: COMPUTE_SPHA
          procedure :: GET_DIM
          procedure :: HARM2ORB
          procedure :: INDF
          ! now called internally by 'compute_spha'
          procedure, private :: REDUCE_SPHA
      END TYPE SPHER_HARM_t

      contains

      subroutine COMPUTE_SPHA( this, first_gindex, last_gindex,
     $                         func, l_inc )

!     Find the spherical-harmonic decomposition of a set of functions
!     from the 'matel-registry'      
!     The operations are performed in parallel

      implicit none
      
      class(SPHER_HARM_t) :: this

      !> First and last global registry indexes of the set of functions
      integer, intent(in) :: first_gindex, last_gindex

      !> Function to be decomposed. Typically it is "evaluate", as in
      !> the 'matel-registry' interface.  It can include an implicit
      !> matrix element function (for example, 'evaluate_x' includes the
      !> product by 'x' of radial functions).

      !> The evaluator selects the function on the basis of a single index 'ig'
      !> that is de-referenced in the 'matel_registry'. It returns the value and the gradient.
      INTERFACE
        subroutine func(ig,rvec,yy,grady)
          use precision, only : dp
          integer, intent(in)   :: ig
          real(dp), intent(in)  :: rvec(3)
          real(dp), intent(out) :: yy, grady(3)
        end subroutine func
      END INTERFACE

      !> When the implied matrix element function increases the maximum l
      !> e.g. for x*\Phi, the maximum l is increased by 1, so l_inc=1
      integer, intent(in) :: l_inc

!     Local variables
      integer             :: i, ig, il, first, last, nlocal, nf, mf,
     &                       l, nilm, jlm

      this%N      = last_gindex - first_gindex + 1
      this%BASE   = first_gindex - 1  ! this is the base for counting

      ! Split the work on the N functions among processors
      nlocal      = GET_LOOP_LIMITS( this%N, Node, first, last )
      nf          = 0

      ! Allocate space
      nullify(this%IND,this%F,this%ILM,this%INVP)
      CALL RE_ALLOC( this%IND, 1, this%n+1, 'F', MYNAME )
      mf = MAX(INT(nlocal*EXPAND),nlocal+MINEXPAND)
      CALL RE_ALLOC( this%INVP, 1, MF, 'INVP', MYNAME )
      CALL RE_ALLOC( this%F, 0, NQ, 1, MF, 'F', MYNAME )
      CALL RE_ALLOC( this%ILM, 1, MF, 'ILM', MYNAME )

      ! Iterate local orbitals
      do IL = first, last
        ! Get the Global Index of the orbital
        IG = this%BASE+IL
        l = LCUT( IG )            ! Maximum l
        nilm = (l+2)**2           ! Total maximum number of l,m indexes
        ! Check if we have enough memory
        if (nf+nilm > mf ) THEN
          mf = (nf+NILM)*EXPAND
          CALL RE_ALLOC( this%INVP, 1, mf, 'INVP', MYNAME )
          CALL RE_ALLOC( this%F, 0, NQ, 1, mf, 'F', MYNAME )
          CALL RE_ALLOC( this%ILM, 1, mf, 'ILM', MYNAME )
        endif
        ! Compute Spherical harmonic expansion of orbital IG
        CALL YLMEXP( L+L_INC, RLYLM, func, IG, 0, NQ, RMAX, NILM,
     &             this%ILM(nf+1:), this%F(:,nf+1:) )
        ! Store orbital in Fourier space
        DO JLM = 1, NILM
          nf = nf + 1
          ! Save the orbital ID of current harmonic
          this%INVP(nf) = IL
          l = LOFILM( this%ILM(nf) )
          CALL RADFFT( L, NQ, RMAX, this%F(0:NQ,NF),
     &      this%F(0:NQ,NF) )
          ! F(NQ,NF) = 0._dp
        ENDDO
        this%IND(IL) = NILM
      enddo

      ! Replicate data to all processors
      call REDUCE_SPHA( this)
      
      end subroutine COMPUTE_SPHA

      subroutine REDUCE_SPHA( this )
C *******************************************************************
C  Reduce the Spherical Harmonics
C ************************* DESCRI **********************************
C  Reduce the values of the different arrays of SPHER_HARM_t.
C  At input, IND contains the number of harmonics of every orbital.
C  At exit, it can be used to index the harmonics of every orbital
C *******************************************************************
#ifdef MPI
      use mpi_siesta,  only : MPI_INTEGER, MPI_COMM_WORLD
      use mpi_siesta,  only : MPI_DOUBLE_PRECISION, MPI_IN_PLACE
#endif      
      implicit none

      class(SPHER_HARM_t) :: this

!     Local variables
      integer             :: i, ig, first, last, nf, prev
#ifdef MPI
      integer             :: MPIERR
      integer,    pointer :: count(:), displ(:), coun2(:), disp2(:)
      real(dp),   pointer :: G_F(:,:)
      integer,    pointer :: G_ILM(:), G_INVP(:)
#endif
#ifdef MPI
      ! reduce IND array in place
      if (NODES.GT.1) then
        nullify(count,displ)
        CALL RE_ALLOC( COUNT, 0, NODES-1, 'COUNT,', MYNAME )
        CALL RE_ALLOC( DISPL, 0, NODES, 'DISPL,', MYNAME )
        DISPL(0) = 0
        DO I= 0, NODES-1
          COUNT(I) = GET_LOOP_LIMITS( this%n, I, first, last )
          DISPL(I+1) = DISPL(I) + COUNT(I)
        ENDDO
        call MPI_Allgatherv( MPI_IN_PLACE, COUNT(NODE), MPI_INTEGER,
     &   this%IND, COUNT, DISPL, MPI_INTEGER, MPI_COMM_WORLD, MPIERR )
      endif
#endif
      ! IND array should be acumulative
      prev        = this%IND(1)
      this%IND(1) = 1
      DO I= 1, this%N
        IG = prev + this%IND(I)
        prev = this%IND(I+1)
        this%IND(I+1) = IG
      ENDDO
#ifdef MPI
      if (NODES.GT.1) then
        ! Reduce INVP, ILM and F between processors
          nullify(coun2,disp2)
          CALL RE_ALLOC( COUN2, 0, NODES-1, 'COUN2,', MYNAME )
          CALL RE_ALLOC( DISP2, 0, NODES, 'DISP2,', MYNAME )
          DISP2(0) = 0
          nf       = 0
          DO I= 0, NODES-1
            DISP2(I+1) = this%IND(DISPL(I+1)+1)-1
            COUN2(I)   = DISP2(I+1) - DISP2(I)
          ENDDO
          nf = DISP2(NODES)

          ! Allocate arrays and reduce INVP
          nullify(G_INVP)
          CALL RE_ALLOC( G_INVP, 1, nf, 'INVP', MYNAME )
          call MPI_Allgatherv( this%INVP, COUN2(NODE), MPI_INTEGER,
     &      G_INVP, COUN2, DISP2, MPI_INTEGER, MPI_COMM_WORLD, MPIERR)
          CALL DE_ALLOC( this%INVP, 'INVP', MYNAME )
          this%INVP => G_INVP

          ! Allocate arrays and reduce ILM
          nullify(G_ILM)
          CALL RE_ALLOC( G_ILM, 1, nf, 'ILM', MYNAME )
          call MPI_Allgatherv( this%ILM, COUN2(NODE), MPI_INTEGER,
     &      G_ILM, COUN2, DISP2, MPI_INTEGER, MPI_COMM_WORLD, MPIERR )
          CALL DE_ALLOC( this%ILM, 'ILM', MYNAME )
          this%ILM => G_ILM

          ! Allocate arrays and reduce F
          COUN2 = COUN2*(NQ+1)
          DISP2 = DISP2*(NQ+1)
          nullify(G_F)
          CALL RE_ALLOC( G_F, 0, NQ, 1, nf, 'F', MYNAME )
          call MPI_Allgatherv( this%F, COUN2(NODE),
     &      MPI_DOUBLE_PRECISION, G_F, COUN2, DISP2,
     &      MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, MPIERR )
          CALL DE_ALLOC( this%F, 'F', MYNAME )
          this%F => G_F

          ! Dealloc temporal arrays
          CALL DE_ALLOC( COUN2, 'COUN2,', MYNAME )
          CALL DE_ALLOC( DISP2, 'DISP2,', MYNAME )
          CALL DE_ALLOC( COUNT, 'COUNT,', MYNAME )
          CALL DE_ALLOC( DISPL, 'DISPL,', MYNAME )
      endif
#endif
      end subroutine REDUCE_SPHA

      function GET_DIM( this, n ) result(d)
C *******************************************************************
!  Returns the accumulated number of spherical-harmonic terms at position n
!  of the 'Spha' data structure, as tracked by the pointer array IND.
!  The total number of terms is get_dim(this%N)
C ************************* INPUT ***********************************
C  integer n      : position to check. (1 <= n <= this%N)
C ************************* OUTPUT **********************************
C  integer d      : Number of terms
C *******************************************************************
      implicit none
      class(SPHER_HARM_t) :: this
      integer, intent(in) :: n
      integer             :: d
      d = this%IND(n+1)-1
      end function GET_DIM

      function Harm2Orb( this, id ) result(d)
C *******************************************************************
! Returns the global registry Id of the function from which the current harmonic
! originates. 
C ************************* INPUT ***********************************
C  integer id      : harmonic Id
C ************************* OUTPUT **********************************
C  integer n      : Orbital Id
C *******************************************************************
      implicit none
      class(SPHER_HARM_t) :: this
      integer, intent(in) :: id
      integer             :: d
      d = this%Invp(id)
      end function Harm2Orb

      function INDF( this, id ) result(d)
C *******************************************************************
C Return the index of the first harmonic of an orbital
C ************************* INPUT ***********************************
C  integer id      : Orbital Id
C ************************* OUTPUT **********************************
C  integer n      : Harmonic Id
C *******************************************************************
      implicit none
      class(SPHER_HARM_t) :: this
      integer, intent(in) :: id
      integer             :: d
      d = this%Ind(id-this%BASE)
      end function INDF

      function GET_LOOP_LIMITS( N, ID,  FIRST, LAST ) RESULT(SIZE)
C *******************************************************************
C  Get the loop limits resulting to split N elements between Nodes
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

      end module matel_ylm
