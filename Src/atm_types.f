      module atm_types

      use radial
!
!     Derived types for orbitals and KB projectors
!
      implicit none

!
!     Storage of orbital and projector real-space tables and other
!     characteristics
!
         
!
!     These parameters are over-dimensioned, but there is no storage
!     penalty, as the real information is packed and indexed.
!
!
      integer, parameter  :: maxn_pjnl = 10
!       Maximum number of projectors (not counting different "m" copies)
      integer, parameter  :: maxn_orbnl = 20
!       Maximum number of nl orbitals (not counting different "m" copies)

      integer, parameter  :: maxnorbs = 100
!       Maximum number of nlm orbitals
      integer, parameter  :: maxnprojs = 50
!       Maximum number of nlm projectors
!

!
!     Species_info: Consolidate all the pieces of information in one place
!
      type species_info
         character(len=2)                ::  symbol
         character(len=20)               ::  label
         integer                         ::  z          ! Atomic number
         double precision                ::  mass
         integer                         ::  zval       ! Valence charge
         double precision                ::  self_energy !Electrostatic
                                                         !self-energy
!
!        Orbitals
!             We keep track of just one orbital for each
!             "nl" family
!
         integer                         ::  n_orbnl    ! num of nl orbs
         integer                         ::  lmax_basis ! basis l cutoff
         integer, dimension(maxn_orbnl)  ::  orbnl_l    ! l of each nl orb
         integer, dimension(maxn_orbnl)  ::  orbnl_n    ! n of each nl orb
         integer, dimension(maxn_orbnl)  ::  orbnl_z    ! z of each nl orb
         logical, dimension(maxn_orbnl)  ::  orbnl_ispol! is it a pol. orb?

         double precision,
     $            dimension(maxn_orbnl)  ::  orbnl_pop  ! pop. of nl orb
                                                        ! (total of 2l+1
                                                        ! components)
!
!        Projectors
!             For each l, there can be several projectors. Formally, we 
!             can can use the "nl" terminology for them. n will run from
!             1 to the total number of projectors at that l.
!             
!
         integer                         ::  n_pjnl     ! num of "nl" projs
         integer                         ::  lmax_projs ! l cutoff for projs
         integer, dimension(maxn_pjnl)   ::  pjnl_l     ! l of each nl proj
         integer, dimension(maxn_pjnl)   ::  pjnl_n     ! n of each nl proj
         double precision, dimension(maxn_pjnl)
     $                                   ::  pjnl_ekb   ! energy of
                                                         ! each nl proj
!
!                        ------------------------------
!
!        Aggregate numbers of orbitals and projectors (including 2l+1
!        copies for each "nl"), and index arrays keeping track of
!        which "nl" family they belong to, and their n, l, and m (to avoid
!        a further dereference)
!
         integer                         ::  norbs
         integer, dimension(maxnorbs)    ::  orb_index
         integer, dimension(maxnorbs)    ::  orb_n
         integer, dimension(maxnorbs)    ::  orb_l
         integer, dimension(maxnorbs)    ::  orb_m
         double precision,
     $            dimension(maxnorbs)    ::  orb_pop   ! pop. of nl orb

         integer                         ::  nprojs
         integer, dimension(maxnprojs)   ::  pj_index
         integer, dimension(maxnprojs)   ::  pj_n
         integer, dimension(maxnprojs)   ::  pj_l
         integer, dimension(maxnprojs)   ::  pj_m
!
         type(rad_func), dimension(:), pointer         ::  orbnl
         type(rad_func), dimension(:), pointer         ::  pjnl
         type(rad_func)                                ::  vlocal
         type(rad_func)                                ::  chlocal
         logical                                       ::  there_is_core
         type(rad_func)                                ::  core


      end type species_info

!
      integer, save             :: nspecies
      integer, save             :: npairs

      type(species_info), target, allocatable, save   ::  species(:)
      type(rad_func), allocatable, target, save     ::  elec_corr(:)
!

      end module atm_types
