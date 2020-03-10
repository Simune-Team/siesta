!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
!
      !> Re-designed module to allow a pre-computation of the needed
      !> matrix elements in parallel, followed by a globalization of the
      !> data among all the MPI processes.  Once the interpolation
      !> tables are setup, further calls to the matrix-element evaluator
      !> are cheap. This has a dramatic effect in some routines (such as
      !> nlefsm) that had to perform the table-building operations under
      !> conditions that did not scale in parallel.

      !> Concept: Rogeli Grima (BSC) and Alberto Garcia (ICMAB)
      !> Initial implementation:  Rogeli Grima (BSC)

      module matel_mod

      use matel_ylm, only: spher_harm_t
      use matel_table, only: matel
      use matel_table, only: MODE_S, MODE_T, MODE_XYZ
      
      use precision, only : dp
      use alloc,     only : alloc_default, allocDefaults

      use m_radfft,  only : reset_radfft
      use m_matel_registry, ONLY : EVALUATE,
     &                             EVALUATE_X, EVALUATE_Y, EVALUATE_Z

      PRIVATE
      
      public :: init_matel     ! Initialize tables S, T, TA, X, Y, Z, S_opt, X_opt, Y_opt, Z_opt
      public :: init_matel_wannier ! Initialize table S_wann
      
      public :: new_matel   ! Evaluates matrix elements with the old interface
      
      !> These relate to the functions recorded in the 'matel registry'; the number of each kind
      !> serve as markers for the different sections needed in the tables.
      !> 
      integer           :: num_orb  ! Number of different orbitals
      integer           :: num_kb  ! Number of different KB projectors
      integer           :: num_ldau ! Number of different LDA+U functions
      integer           :: num_va  ! Number of different Vna functions
      integer           :: num_wannier_projs  ! Number of different Wannier trial orbitals

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
      TYPE(MATEL) :: tab_S           ! Unity (overlap). 
      !> Laplacian among PAOs
      TYPE(MATEL) :: tab_T           ! -Laplacian.      
      !> Laplacian among Vna functions
      TYPE(MATEL) :: tab_TA         ! -Laplacian.
      
      !> X, Y, Z among orbitals (for dielectric polarization calculations)
      TYPE(MATEL) :: tab_X           ! X projection.   
      TYPE(MATEL) :: tab_Y           ! Y projection.   
      TYPE(MATEL) :: tab_Z           ! Z projection.   

      !> Overlaps of PAOs with wannier projectors
      TYPE(MATEL) :: tab_S_wann      ! Unity (overlap).

      !> Overlaps of PAOs with KBs for optical calculations (different from table "S" because
      !> the KBs are "function 2")
      TYPE(MATEL) :: tab_S_opt        ! Unity (overlap)

      !> X,Y,Z  PAOs with KBs for optical calculations
      TYPE(MATEL) :: tab_X_opt        ! X projection
      TYPE(MATEL) :: tab_Y_opt        ! Y projection
      TYPE(MATEL) :: tab_Z_opt        ! Z projection

      CONTAINS

      !>  Initialize the MATEL tables S, T, TA, X, Y, Z, S_opt,
      !>  X_opt, Y_opt, Z_opt

      !> *NOTE: Optical and Polarization tables are initialized
      !> on the basis of "static" information in the fdf file. If
      !> and when those operations can be launched "after-the-fact"
      !> this might have to be changed.
      !> Same if the 'LDAU' features are turned on after the first reading.

      subroutine init_matel( )

      use fdf,       only : fdf_boolean
      use atm_types, only : nspecies, species
      use siesta_options, only : bornz
      use ldau_specs, only : switch_ldau

      implicit none
      ! Local Variables
      integer             :: is, top, top_vna
      logical             :: optical
      logical,       save :: init = .false.
      TYPE(allocDefaults) :: OLDEFS

      
      ! Initialize only once
      if (init) return
      init = .true.

       ! Set allocation defaults
      call alloc_default( old=oldefs, copy=.true., shrink=.false. )

      CALL RESET_RADFFT( )

      ! Get the number of orbitals, KB projs, LDAU projs and Neutral-atom potentials

      NUM_ORB = 0
      NUM_KB = 0
      NUM_LDAU = 0
      do is= 1, nspecies
        NUM_ORB  = NUM_ORB  + species(is)%norbs
        NUM_KB   = NUM_KB   + species(is)%nprojs
        NUM_LDAU = NUM_LDAU + species(is)%nprojsldau
      enddo
      NUM_VA = nspecies

      ! Find the spherical harmonic decomposition and
      ! Fourier transform of all the registry objects
      ! except the wannier trial functions
      top = num_orb + num_kb + num_ldau 
      top_vna = top + num_va
      call ylmk_val_all%compute_spha(1, top_vna, EVALUATE, 0)

      ! Overlap table from PAOs, KBs, and LDAU projs to PAOs.
      call tab_S%init(MODE_S, ylmk_val_all, 1, top,
     $                    ylmk_val_all, 1, num_orb)

      ! Laplacian table from PAOs to PAOs.
      call tab_T%init(MODE_T, ylmk_val_all, 1, num_orb,
     $                    ylmk_val_all, 1, num_orb)
      
      ! Laplacian table from Vnas to Vnas
      call tab_TA%init(MODE_T, ylmk_val_all, top+1, top_vna,
     $                     ylmk_val_all, top+1, top_vna)


      ! This section executed only if polarization calculations
      ! are requested (they need matrix elements of the position operator)
      if (bornz) then

      ! Find the spherical harmonic decomposition and
      ! Fourier transform of {x,y,z}*\phi, where \phi is a PAO

      ! The X, Y, and Z tables are then overlaps of PAOs and xPAOs

         call ylmk_x_orbs%compute_spha(1, num_orb, EVALUATE_X, 1)
         call tab_X%init( MODE_XYZ, ylmk_val_all, 1, num_orb,
     $                         ylmk_x_orbs, 1, num_orb)

         call ylmk_y_orbs%compute_spha(1, num_orb, EVALUATE_Y, 1)
         call tab_Y%init( MODE_XYZ, ylmk_val_all, 1, num_orb,
     $                         ylmk_y_orbs, 1, num_orb)

         call ylmk_z_orbs%compute_spha(1, num_orb, EVALUATE_Z, 1)
         call tab_X%init( MODE_XYZ, ylmk_val_all, 1, num_orb,
     $                         ylmk_z_orbs, 1, num_orb)

      endif

      optical = fdf_boolean('OpticalCalculation',.false.)

      ! This section executed only if an 'optical' calculation
      ! is requested 
      if (optical) then
         
      ! Overlap table from PAOs to KBs.
         call tab_S_opt%init(MODE_S, ylmk_val_all, 1, num_orb,
     $                    ylmk_val_all, num_orb+1, num_kb)

      ! Find the spherical harmonic decomposition and
      ! Fourier transform of {x,y,z}*\chi, where \chi is a KB proj.

      ! The X_opt, Y_opt, and Z_opt tables are then overlaps of PAOs and xKBs

         call ylmk_x_kbs%compute_spha(num_orb+1, num_orb+num_kb,
     $                                EVALUATE_X, 1)

         call tab_X_opt%init( MODE_XYZ, ylmk_val_all, 1, num_orb,
     $                         ylmk_x_kbs, 1, num_kb)

         call ylmk_y_kbs%compute_spha(num_orb+1, num_orb+num_kb,
     $                                EVALUATE_Y, 1)

         call tab_Y_opt%init( MODE_XYZ, ylmk_val_all, 1, num_orb,
     $                         ylmk_y_kbs, 1, num_kb)

         call ylmk_z_kbs%compute_spha(num_orb+1, num_orb+num_kb,
     $                                EVALUATE_Z, 1)

         call tab_Z_opt%init( MODE_XYZ, ylmk_val_all, 1, num_orb,
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

      num_wannier_projs = numproj ! set module variable
      top_vna = num_orb + num_kb + num_ldau + num_va

      call ylmk_val_wannier_projs%compute_spha(top_vna+1,
     $                           top_vna+numproj, EVALUATE, 0)

      call tab_S_wann%init(MODE_S, ylmk_val_all, 1, num_orb,
     $                     ylmk_val_wannier_projs, 1, numproj)

      ! Restore allocation defaults
      CALL alloc_default( restore=oldefs )
      end subroutine init_matel_wannier

      SUBROUTINE new_MATEL( OPERAT, IG1, IG2, R12, S12, DSDR )

! Finds two-center matrix elements between 'atomic orbitals' 
! with finite radial and angular momentum cutoffs.

! The interface is the same as "new_matel" in the m_new_matel
! module, but the routine itself only dispatches execution
! to the new table objects (Alberto Garcia, July 2019, after Rogeli Grima)      

! The tables must be initialized previously by calls to init_matel and
! (if needed) init_matel_wannier.
      
! INPUT 

C CHARACTER OPERAT : Operator to be used. The valid options are:
C   'S' => Unity (overlap). Uppercase required for all values.
C   'T' => -Laplacian
C   'U' => 1/|r'-r| (with evaluate returning charge distributions) *NOT IMPLEMENTED
C   'X' => x, returning <phi1(r-R12)|x|phi2(r)> (origin on second atom)
C   'Y' => y, returning <phi1(r-R12)|y|phi2(r)>
C   'Z' => z, returning <phi1(r-R12)|z|phi2(r)>
C INTEGER IG1    : Global index of 1st function (must be positive)
C INTEGER IG2    : Global index of 2nd function
C                    Indexes IG1, IG2 are used only to call 
C                    routines LCUT, RCUT and EVALUATE (in matel_registry), and 
C                    may have other meanings within those routines
C REAL*8  R12(3) : Vector from first to second atom
C ************************* OUTPUT **********************************
C REAL*8 S12      : Matrix element between orbitals.
C REAL*8 DSDR(3)  : Derivative (gradient) of S12 with respect to R12.

C Length units are arbitrary, but must be consistent in MATEL, RCUT
C   and EVALUATE. The laplacian unit is (length unit)**(-2).
C ************************* BEHAVIOUR *******************************
C If |R12| > RCUT(IS1,IO1) + RCUT(IS2,IO2), returns exactly zero.

      !> Operation code
      character(len=1), intent(in) :: OPERAT
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

      integer :: top
      
      top = num_orb + num_kb + num_ldau 

      select case (operat)
      case ( 'S' )

         if ( (ig1 <= top) .and. (is_orb(ig2)) ) then
            call tab_S%get_matel(IG1, IG2, R12, S12, DSDR)
         else if (is_orb(ig1) .and. (is_kb(ig2))) then
            call tab_S_opt%get_matel(IG1, IG2, R12, S12, DSDR)
         else if (is_orb(ig1) .and. (is_wannier_proj(ig2))) then
            call tab_S_wann%get_matel(IG1, IG2, R12, S12, DSDR)
         else
            call die("Cannot process 'S' in matel")
         endif

      case ( 'T' )

         if ( is_orb(ig1) .and. (is_orb(ig2))) then
            call tab_T%get_matel(IG1, IG2, R12, S12, DSDR)
         else if  (is_vna(ig1) .and. (is_vna(ig2))) then
            call tab_TA%get_matel(IG1, IG2, R12, S12, DSDR)
         else
            call die("Cannot process 'T' in matel")
         endif
         
      case ( 'X' )

         if  (is_orb(ig1) .and. (is_orb(ig2))) then
            call tab_X%get_matel(IG1, IG2, R12, S12, DSDR)
         else if (is_orb(ig1) .and. (is_kb(ig2))) then
            call tab_X_opt%get_matel(IG1, IG2, R12, S12, DSDR)
         else
            call die("Cannot process 'X' in matel")
         endif
         
      case ( 'Y' )

         if ( is_orb(ig1) .and. (is_orb(ig2))) then
            call tab_Y%get_matel(IG1, IG2, R12, S12, DSDR)
         else if (is_orb(ig1) .and. (is_kb(ig2))) then
            call tab_Y_opt%get_matel(IG1, IG2, R12, S12, DSDR)
         else
            call die("Cannot process 'Y' in matel")
         endif
         
      case ( 'Z' )

         if  (is_orb(ig1) .and. (is_orb(ig2))) then
            call tab_Z%get_matel(IG1, IG2, R12, S12, DSDR)
         else if (is_orb(ig1) .and. (is_kb(ig2))) then
            call tab_Z_opt%get_matel(IG1, IG2, R12, S12, DSDR)
         else
            call die("Cannot process 'Z' in matel")
         endif
         
         case ( 'U' )
            call die("Cannot process 'U' in matel")

         case default
            call die("Unrecognized 'operat' in matel")

         end select

      end subroutine new_matel

      function is_orb(ig) result(res)
      integer, intent(in) :: ig
      logical             :: res

      res = ( ig > 0 .and. ig <= num_orb)
      end function is_orb

      function is_kb(ig) result(res)
      integer, intent(in) :: ig
      logical             :: res

      res = ( ig > num_orb .and. ig <= (num_orb+num_kb))
      end function is_kb

      function is_vna(ig) result(res)
      integer, intent(in) :: ig
      logical             :: res

      integer :: top

      top = num_orb + num_kb + num_ldau 
      res = ( ig > top .and. ig <= (top + num_va))
      end function is_vna

      function is_wannier_proj(ig) result(res)
      integer, intent(in) :: ig
      logical             :: res

      integer :: top_vna

      top_vna = num_orb + num_kb + num_ldau + num_va
      
      res = ( ig > top_vna .and. ig <= (top_vna + num_wannier_projs))
      end function is_wannier_proj
      
      end module matel_mod
