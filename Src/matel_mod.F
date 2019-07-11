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
      !> evaluator (here renamed 'get_matel_?') are cheap. This has a dramatic effect in
      !> some routines (such as nlefsm) that had to perform the table-building operations
      !> under conditions that did not scale in parallel.

      !> Concept: Rogeli Grima (BSC) and Alberto Garcia (ICMAB)
      !> Initial implementation:  Rogeli Grima (BSC)

      module matel_mod

!      use matel_params
      use matel_ylm, only: spher_harm_t
      use matel_table, only: matel
      use matel_table, only: MODE_S, MODE_T, MODE_XYZ
      
      use precision, only : dp
      use alloc,     only : alloc_default, allocDefaults

      use spher_harm, only : reset_spher_harm
      use m_radfft,  only : reset_radfft
      use m_matel_registry, ONLY : EVALUATE,
     &                             EVALUATE_X, EVALUATE_Y, EVALUATE_Z


      PRIVATE
      
      public :: init_matel     ! Initialize S, T, TA, X, Y, Z, S_op, X_op, Y_op, Z_op
      public :: init_matel_wannier  ! Initialize S_wann
      public :: get_matel_s    ! Get a value from S, S_wann, or S_op
      public :: get_matel_t    ! Get a value from T or TA
      public :: get_matel_x    ! Get a value from X or X_op
      public :: get_matel_y    ! Get a value from Y or Y_op
      public :: get_matel_z    ! Get a value from Z or Z_op

      !> These are the functions recorded in the 'matel registry'; the number of each kind
      !> serve as markers for the different sections needed in the tables.
      !> 
      integer           :: num_orb  ! Number of different orbitals
      integer           :: num_kb  ! Number of different KB projectors
      integer           :: num_ldau ! Number of different LDA+U functions
      integer           :: num_va  ! Number of different Vna functions

      !> k-space spherical-harmonic decomposition of
      !> registered functions
      type(spher_harm_t), target :: ylmk_val_all
      type(spher_harm_t), target :: ylmk_x_orbs
      type(spher_harm_t), target :: ylmk_y_orbs
      type(spher_harm_t), target :: ylmk_z_orbs
      type(spher_harm_t), target :: ylmk_x_kbs
      type(spher_harm_t), target :: ylmk_y_kbs
      type(spher_harm_t), target :: ylmk_z_kbs
      type(spher_harm_t), target :: ylmk_val_wannier_projs


      !> There are different tables, each appropriate to a given operation and kind of function.
      !> The indexing is a bit cumbersome due to the one-dimensional nature of the matel_registry,
      !> in which all functions are stored in the same section: first orbitals, then KB projectors,
      !> then Vna, and finally LDA+U projectors.
      !> In addition, when using the Wannier interface, "trial orbitals" (called also "projectors"
      !> in the code (numproj of them) need to be dealt with.

      !> Overlaps of PAOs, KBs, or LDAUprojs with PAOs
      TYPE(MATEL) :: S           ! Unity (overlap). 
      !> Laplacian among PAOs
      TYPE(MATEL) :: T           ! -Laplacian.      
      !> Laplacian among Vna functions
      TYPE(MATEL) :: TA         ! -Laplacian.
      
      !> X, Y, Z among orbitals (for dielectric polarization calculations)
      TYPE(MATEL) :: X           ! X projection.   
      TYPE(MATEL) :: Y           ! Y projection.   
      TYPE(MATEL) :: Z           ! Z projection.   

      !> Overlaps of PAOs with wannier projectors
      TYPE(MATEL) :: S_wann      ! Unity (overlap).

      !> Overlaps of PAOs with KBs for optical calculations (different from table "S" because
      !> the KBs are "function 2")
      TYPE(MATEL) :: S_op        ! Unity (overlap)

      !> X,Y,Z  PAOs with KBs for optical calculations
      TYPE(MATEL) :: X_op        ! X projection
      TYPE(MATEL) :: Y_op        ! Y projection
      TYPE(MATEL) :: Z_op        ! Z projection

      CONTAINS

      !>  Initialize the MATEL tables S, T, TA, X, Y, Z, S_op,
      !>  X_op, Y_op, Z_op
      subroutine init_matel( )

      use fdf,       only : fdf_boolean
      use atm_types, only : nspecies, species
      use atmfuncs,  only : nofis, nkbfis, orb_gindex, kbproj_gindex
      use atmfuncs,  only : vna_gindex, ldau_gindex
      use siesta_options, only : bornz, w90_processing
      use m_siesta2wannier90, only : numproj
      use ldau_specs, only : switch_ldau
      use m_matel_registry, only : get_nfuncs

      implicit none
      ! Local Variables
      integer             :: is, IO, KO, top, top_vna
      logical             :: optical
      logical,       save :: init = .false.
      TYPE(allocDefaults) :: OLDEFS

      
      ! Initialize only once
      if (init) return
      init = .true.

       ! Set allocation defaults
      call alloc_default( old=oldefs, copy=.true., shrink=.false. )

      CALL RESET_SPHER_HARM( )
      CALL RESET_RADFFT( )

      ! Get the number of orbitals and KB orbitals
      IO = NOFIS( nspecies )
      KO = NKBFIS( nspecies )
      NUM_ORB = orb_gindex( nspecies, IO )
      NUM_KB = kbproj_gindex( nspecies, -KO ) - NUM_ORB
      NUM_VA = nspecies
      NUM_LDAU = 0
      do is= 1, nspecies
        NUM_LDAU = NUM_LDAU + species(is)%nprojsldau
      enddo

      ! Find the spherical harmonic decomposition and
      ! Fourier transform of all the registry objects
      ! except the wannier trial functions
      top = num_orb + num_kb + num_ldau 
      top_vna = top + num_va
      call ylmk_val_all%compute_spha(1, top_vna, EVALUATE, 0)

      ! Overlap table from PAOs, KBs, and LDAU projs to PAOs.
      call S%init(MODE_S, ylmk_val_all, 1, top,
     $                    ylmk_val_all, 1, num_orb)

      ! Laplacian table from PAOs to PAOs.
      call T%init(MODE_T, ylmk_val_all, 1, num_orb,
     $                    ylmk_val_all, 1, num_orb)
      
      ! Laplacian table from Vnas to Vnas
      call TA%init(MODE_T, ylmk_val_all, top+1, top_vna,
     $                     ylmk_val_all, top+1, top_vna)


      ! This section executed only if polarization calculations
      ! are requested (they need matrix elements of the position operator)
      if (bornz) then

      ! Find the spherical harmonic decomposition and
      ! Fourier transform of {x,y,z}*\phi, where \phi is a PAO

      ! The X, Y, and Z tables are then overlaps of PAOs and xPAOs

         call ylmk_x_orbs%compute_spha(1, num_orb, EVALUATE_X, 1)
         call X%init( MODE_XYZ, ylmk_val_all, 1, num_orb,
     $                         ylmk_x_orbs, 1, num_orb)

         call ylmk_y_orbs%compute_spha(1, num_orb, EVALUATE_Y, 1)
         call Y%init( MODE_XYZ, ylmk_val_all, 1, num_orb,
     $                         ylmk_y_orbs, 1, num_orb)

         call ylmk_z_orbs%compute_spha(1, num_orb, EVALUATE_Z, 1)
         call X%init( MODE_XYZ, ylmk_val_all, 1, num_orb,
     $                         ylmk_z_orbs, 1, num_orb)

      endif

      optical = fdf_boolean('OpticalCalculation',.false.)

      ! This section executed only if an 'optical' calculation
      ! is requested 
      if (optical) then
         
      ! Overlap table from PAOs to KBs.
         call S_op%init(MODE_S, ylmk_val_all, 1, num_orb,
     $                    ylmk_val_all, num_orb+1, num_kb)

      ! Find the spherical harmonic decomposition and
      ! Fourier transform of {x,y,z}*\chi, where \chi is a KB proj.

      ! The X_op, Y_op, and Z_op tables are then overlaps of PAOs and xKBs

         call ylmk_x_kbs%compute_spha(num_orb+1, num_orb+num_kb,
     $                                EVALUATE_X, 1)

         call X_op%init( MODE_XYZ, ylmk_val_all, 1, num_orb,
     $                         ylmk_x_kbs, 1, num_kb)

         call ylmk_y_kbs%compute_spha(num_orb+1, num_orb+num_kb,
     $                                EVALUATE_Y, 1)

         call Y_op%init( MODE_XYZ, ylmk_val_all, 1, num_orb,
     $                         ylmk_y_kbs, 1, num_kb)

         call ylmk_z_kbs%compute_spha(num_orb+1, num_orb+num_kb,
     $                                EVALUATE_Z, 1)

         call Z_op%init( MODE_XYZ, ylmk_val_all, 1, num_orb,
     $                         ylmk_z_kbs, 1, num_kb)

         
      endif

      ! Restore allocation defaults
      CALL alloc_default( restore=oldefs )
      
      end subroutine init_matel
!--------------------------------------------------------------------
      subroutine init_matel_wannier( numproj )
!
!     Initialize the table of overlaps between orbitals and
!     Wannier trial functions, once the number of functions is known

      implicit none
      integer, intent(in) :: numproj

      ! Local Variables
      logical, save :: init = .false.
      TYPE(allocDefaults) :: OLDEFS
      integer :: top_vna
      
      ! Initialize only once
      if (init) return
      init = .true.

      ! Set allocation defaults
      call alloc_default( old=oldefs, copy=.true., shrink=.false. )

      ! This table refers to overlaps between orbitals and
      ! Wannier trial functions, which are registered in the
      ! pool after every other function

      top_vna = num_orb + num_kb + num_ldau + num_va

      call ylmk_val_wannier_projs%compute_spha(top_vna+1,
     $                           top_vna+numproj, EVALUATE, 0)

      call S_wann%init(MODE_S, ylmk_val_all, 1, num_orb,
     $                     ylmk_val_wannier_projs, 1, numproj)

      ! Restore allocation defaults
      CALL alloc_default( restore=oldefs )
      end subroutine init_matel_wannier

      subroutine get_matel_s( IG1, IG2, R12, S12, DSDR )
C *******************************************************************
C  Get the Matrix element between orbitals and its derivative from
C  one of the S matrices
C ************************* INPUT ***********************************
C INTEGER IG1    : Global index of 1st function
C INTEGER IG2    : Gloabal index of 2nd function
C REAL*8  R12(3) : Vector from first to second atom
C ************************* OUTPUT **********************************
C REAL*8 S12      : Matrix element between orbitals.
C REAL*8 DSDR(3)  : Derivative (gradient) of S12 with respect to R12.
C *******************************************************************
      implicit none
      integer  :: IG1, IG2
      real(dp) ::  R12(3), DSDR(3), S12
      ! Local Variables
      IF (IG2.LE.num_orb) then
        call S%GET_MATEL( IG1, IG2, R12, S12, DSDR )
      ELSE IF (IG2.LE.num_orb+num_kb) then
        call S_op%GET_MATEL( IG1, IG2, R12, S12, DSDR )
      ELSE
        call S_wann%GET_MATEL( IG1, IG2, R12, S12, DSDR )
      ENDIF
      end subroutine get_matel_s

      subroutine get_matel_t( IG1, IG2, R12, S12, DSDR )
C *******************************************************************
C  Get the Matrix element between orbitals and its derivative from
C  one of the T matrices
C ************************* INPUT ***********************************
C INTEGER IG1    : Global index of 1st function
C INTEGER IG2    : Gloabal index of 2nd function
C REAL*8  R12(3) : Vector from first to second atom
C ************************* OUTPUT **********************************
C REAL*8 S12      : Matrix element between orbitals.
C REAL*8 DSDR(3)  : Derivative (gradient) of S12 with respect to R12.
C *******************************************************************
      implicit none
      integer  :: IG1, IG2
      real(dp) ::  R12(3), DSDR(3), S12
      if (IG1.le.NUM_ORB .and. IG2.le.NUM_ORB) then
        call T%GET_MATEL( IG1, IG2, R12, S12, DSDR )
      else
        call TA%GET_MATEL( IG1, IG2, R12, S12, DSDR )
      endif
      end subroutine get_matel_t

      subroutine get_matel_x( IG1, IG2, R12, S12, DSDR )
C *******************************************************************
C  Get the Matrix element between orbitals and its derivative from
C  one of the X matrices
C ************************* INPUT ***********************************
C INTEGER IG1    : Global index of 1st function
C INTEGER IG2    : Gloabal index of 2nd function
C REAL*8  R12(3) : Vector from first to second atom
C ************************* OUTPUT **********************************
C REAL*8 S12      : Matrix element between orbitals.
C REAL*8 DSDR(3)  : Derivative (gradient) of S12 with respect to R12.
C *******************************************************************
      implicit none
      integer  :: IG1, IG2
      real(dp) ::  R12(3), DSDR(3), S12
      IF (IG2.LE.num_orb) then
        call X%GET_MATEL( IG1, IG2, R12, S12, DSDR )
      ELSE
        call X_op%GET_MATEL( IG1, IG2, R12, S12, DSDR )
      ENDIF
      end subroutine get_matel_x

      subroutine get_matel_y( IG1, IG2, R12, S12, DSDR )
C *******************************************************************
C  Get the Matrix element between orbitals and its derivative from
C  one of the Y matrices
C ************************* INPUT ***********************************
C INTEGER IG1    : Global index of 1st function
C INTEGER IG2    : Gloabal index of 2nd function
C REAL*8  R12(3) : Vector from first to second atom
C ************************* OUTPUT **********************************
C REAL*8 S12      : Matrix element between orbitals.
C REAL*8 DSDR(3)  : Derivative (gradient) of S12 with respect to R12.
C *******************************************************************
      implicit none
      integer  :: IG1, IG2
      real(dp) ::  R12(3), DSDR(3), S12
      IF (IG2.LE.num_orb) then
        call Y%GET_MATEL( IG1, IG2, R12, S12, DSDR )
      ELSE
        call Y_op%GET_MATEL( IG1, IG2, R12, S12, DSDR )
      ENDIF
      end subroutine get_matel_y

      subroutine get_matel_z( IG1, IG2, R12, S12, DSDR )
C *******************************************************************
C  Get the Matrix element between orbitals and its derivative from
C  one of the Z matrices
C ************************* INPUT ***********************************
C INTEGER IG1    : Global index of 1st function
C INTEGER IG2    : Gloabal index of 2nd function
C REAL*8  R12(3) : Vector from first to second atom
C ************************* OUTPUT **********************************
C REAL*8 S12      : Matrix element between orbitals.
C REAL*8 DSDR(3)  : Derivative (gradient) of S12 with respect to R12.
C *******************************************************************
      implicit none
      integer  :: IG1, IG2
      real(dp) ::  R12(3), DSDR(3), S12
      IF (IG2.LE.num_orb) then
        call Z%GET_MATEL( IG1, IG2, R12, S12, DSDR )
      ELSE
        call Z_op%GET_MATEL( IG1, IG2, R12, S12, DSDR )
      ENDIF
      end subroutine get_matel_z

      end module matel_mod
