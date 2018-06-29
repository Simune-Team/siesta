! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module atm_types

      use precision, only: dp
      use radial, only: rad_func
!
!     Derived types for orbitals,  KB projectors, and LDA+U projectors
!
      implicit none
!
!     Storage of orbital and projector real-space tables and other
!     characteristics
!
!     These parameters are over-dimensioned, but there is no storage
!     penalty, as the real information is packed and indexed.
!
      integer, parameter, public  :: maxnorbs = 100
!       Maximum number of nlm orbitals
!
      integer, parameter, public  :: maxn_pjnl = 20
!       Maximum number of projectors (not counting different "m" copies)
      integer, parameter, public  :: maxn_orbnl = 200
!       Maximum number of nl orbitals (not counting different "m" copies)
!       Now very large to accommodate filteret basis sets
      integer, parameter, public  :: maxnprojs = 100
!       Maximum number of nlm projectors
!

!
!     Species_info: Consolidate all the pieces of information in one place
!
      type, public :: species_info
         character(len=2)                ::  symbol
         character(len=20)               ::  label
         integer                         ::  z          ! Atomic number
         real(dp)                        ::  mass
         real(dp)                        ::  zval       ! Valence charge
         real(dp)                        ::  self_energy !Electrostatic
                                                         !self-energy
!
!        Orbitals
!             We keep track of just one orbital for each
!             "nl" family
!
         integer                         ::  n_orbnl=0  ! num of nl orbs
         integer                         ::  lmax_basis=0 ! basis l cutoff
         integer, dimension(maxn_orbnl)  ::  orbnl_l    ! l of each nl orb
         integer, dimension(maxn_orbnl)  ::  orbnl_n    ! n of each nl orb
         integer, dimension(maxn_orbnl)  ::  orbnl_z    ! z of each nl orb
         logical, dimension(maxn_orbnl)  ::  orbnl_ispol! is it a pol. orb?

         real(dp),
     $            dimension(maxn_orbnl)  ::  orbnl_pop  ! pop. of nl orb
                                                        ! (total of 2l+1
                                                        ! components)
!
!        KB Projectors
!             For each l, there can be several projectors. Formally, we 
!             can can use the "nl" terminology for them. n will run from
!             1 to the total number of projectors at that l.
!             
!
         logical                         ::  lj_projs = .false.
         integer                         ::  n_pjnl=0   ! num of "nl" projs
         integer                         ::  lmax_projs=0 ! l cutoff for projs
         integer, dimension(maxn_pjnl)   ::  pjnl_l     ! l of each nl proj
         real(dp), dimension(maxn_pjnl)  ::  pjnl_j     ! j of each nl proj
         integer, dimension(maxn_pjnl)   ::  pjnl_n     ! n of each nl proj
         real(dp), dimension(maxn_pjnl)
     $                                   ::  pjnl_ekb   ! energy of

!
!        Aggregate numbers of orbitals and projectors (including 2l+1
!        copies for each "nl"), and index arrays keeping track of
!        which "nl" family they belong to, and their n, l, and m (to avoid
!        a further dereference)
!
         integer                         ::  norbs=0
         integer, dimension(maxnorbs)    ::  orb_index
         integer, dimension(maxnorbs)    ::  orb_n
         integer, dimension(maxnorbs)    ::  orb_l
         integer, dimension(maxnorbs)    ::  orb_m
         integer, dimension(maxnorbs)    ::  orb_gindex
         real(dp),
     $            dimension(maxnorbs)    ::  orb_pop   ! pop. of nl orb

         integer                         ::  nprojs=0
         integer, dimension(maxnprojs)   ::  pj_index
         integer, dimension(maxnprojs)   ::  pj_n
         integer, dimension(maxnprojs)   ::  pj_l
         real(dp), dimension(maxnprojs)  ::  pj_j
         integer, dimension(maxnprojs)   ::  pj_m
         integer, dimension(maxnprojs)   ::  pj_gindex
!----------------------------
!        LDA+U Projectors
!        Here we follow the scheme used for the KB projectors
!        
         integer                         ::  n_pjldaunl = 0
                                             ! num of "nl" projs
                                             ! not counting the "m copies"
         integer                         ::  lmax_ldau_projs = 0
                                             ! l cutoff for LDA+U proj
         integer, dimension(maxn_pjnl)   ::  pjldaunl_l ! l of each nl proj
         integer, dimension(maxn_pjnl)   ::  pjldaunl_n ! n of each nl proj
                                             ! Here, n is not the principal
                                             ! quantum number, but a sequential
                                             ! index from 1 to the total 
                                             ! number of projectors for that l.
                                             ! In the case of LDA+U projectors,
                                             ! It is always equal to 1.
         real(dp), dimension(maxn_pjnl)  ::  pjldaunl_U ! U of each nl projector
         real(dp), dimension(maxn_pjnl)  ::  pjldaunl_J ! J of each nl projector

         integer                         ::  nprojsldau = 0
                                             ! Total number of LDA+U proj.
                                             ! counting the "m copies"
                                             ! (including the (2l + 1) factor))
         integer, dimension(maxnprojs)   ::  pjldau_index
         integer, dimension(maxnprojs)   ::  pjldau_n
         integer, dimension(maxnprojs)   ::  pjldau_l
         integer, dimension(maxnprojs)   ::  pjldau_m
         integer, dimension(maxnprojs)   ::  pjldau_gindex
!
         type(rad_func), dimension(:), pointer       ::  orbnl
         type(rad_func), dimension(:), pointer       ::  pjnl
         type(rad_func), dimension(:), pointer       ::  pjldau
         type(rad_func)                              ::  vna
         integer                                     ::  vna_gindex=0
         type(rad_func)                              ::  chlocal
         type(rad_func)                              ::  reduced_vlocal
         logical                                     ::  there_is_core
         type(rad_func)                              ::  core

         logical                        :: read_from_file
      end type species_info

!
      integer, save, public             :: nspecies
      integer, save, public             :: npairs

      type(species_info), target, allocatable,
     $                            save, public   ::  species(:)
      type(rad_func), allocatable, target,
     $                            save, public   ::  elec_corr(:)
      

      private

      end module atm_types
