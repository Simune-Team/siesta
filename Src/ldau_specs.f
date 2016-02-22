! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      module ldau_specs
! 
! Javier Junquera, January 2016, based on previous reldauproj
! 
! Processes the information in an fdf file 
! to generate the projectors for the LDA+U simulations, 
! and populates the "projectors specifications" data structures.
! 
! Here is a guide to the behavior of the main routine "read_ldau_specs":
! 
! * Find the generation method of the projectors:
!     The LDA+U projectors are the localized functions used
!     to calculate the local populations used in a Hubbard-like term
!     that modifies the LDA Hamiltonian and  energy.
!     It is important to recall that LDA+U projectors should be
!     quite localized functions.
!     Otherwise the calculated populations lose their atomic character
!     and physical meaning. Even more importantly,
!     the interaction range can increase so much that jeopardizes
!     the efficiency of the calculation.
!
!     Two methods are currently implemented (accepted values are 1 and 2):
!        - If method_gen_ldau_proj = 1 
!          Projectors are slightly-excited numerical atomic orbitals
!          similar to those used as an automatic basis set by  SIESTA.
!          The radii of these orbitals are controlled using
!          the parameter LDAU.EnergyShift and/or the data
!          in block LDAU.proj (quite similar to the data block PAO.Basis used
!          to specify the basis set, see below).
!        - If method_gen_ldau_proj = 2 
!          Projectors are exact solutions of the pseudoatomic
!          problem (and, in principle, are not strictly localized) which are
!          cutoff using a Fermi function $1/\{1+\exp[(r-r_c)\omega]\}$.
!          The values of $r_c$ and $\omega$ are controlled using
!          the parameter LDAU.CutoffNorm and/or the  data
!          block LDAU.proj.
!     The default value is method_gen_ldau_proj = 2
!
! * Find the energy shift to generate the LDA+U projectors
!     Energy increased used to define the localization radious
!     of the LDAU projectors (similar to the parameter PAO.EnergyShift).
!
! * Allocate storage for the data structures
!   in particular the projector pointer that will be used later
!   in 
! * Determine any "global" basis specification parameters:
!   - basis_size   (sz, dz, dzp, etc)
!   - basis_type   (split, nodes, etc) ("USER" is no longer valid)
!   LDAU.proj - This is the most complex block, very flexible but in
!               need  of spelling-out the specific restrictions.
!               It follows the same spirit as the PAO.Basis block.
!               Line by line, the specific info is:
! 
!   1st:   Species_label number_of_l_shells [basis_type] [ionic_charge] 
! 
!   For each l_shell:
!     [n= n] l nzeta [P [ nzeta_pol]] [E vcte rinn] [Q qcoe [qyuk [qwid]]] [F cutoff]
!   where 'n= n' is *mandatory* if the species has any semicore states,
!   and the 'P' (polarization), E (soft confinement potential), 
!   'Q' (charge confinement), and 'F' (filteret) sections are optional and can appear
!   in any order after nzeta. 
! 
!          rc_1  rc_2 .... rc_nzeta  
! 
!   are the cutoff radii in bohrs. This line is mandatory and there must
!   be at least nzeta values (the extra ones are discarded)
! 
!   A line containing contraction (scale) factors is optional, but if it
!   appears, the values *must be* real numbers, and there must be at
!   least nzeta of them.
!   --------------------------------------------------------------------
! 
!   After processing LDAU.proj block, whatever PAO information
!   which is not already available is determined in routine 'autobasis', 
!   using the following defaults:
! 
!   rc(1:nzeta) is set to 0.0
!   lambda(1:nzeta) is set to 1.0  (this is a change from old practice)
! 
!  ----------------------------------
!  
! =======================================================================
!
      use precision
      use fdf

      use sys,         only : die               ! Termination routine
      use basis_specs, only : label2species     ! Function that reads the
                                                !   label of a species and
                                                !   converts it to the 
                                                !   corresponding index 
                                                !   according to the 
                                                !   Chemical_Species_Label block
      use basis_types, only : basis_def_t       ! Derived type where all the
                                                !   information relative to the
                                                !   definition of the basis set
                                                !   is defined
      use basis_types, only : ldaushell_t       ! Derived type where all the
                                                !   information relative to the
                                                !   atomic orbitals where the U
                                                !   correction will be applied 
                                                !   is defined
      use basis_types, only : basis_parameters  ! Derived type where all the
                                                !   information about the 
                                                !   - basis set
                                                !   - Kleinman-Bylander proj.
                                                !   - LDA+U proj. 
                                                !   ...
                                                !   for all the species 
                                                !   are defined
      use basis_types, only: initialize         ! Subroutine to initialize
                                                !   the values of some derived
                                                !   types
      use basis_types, only: print_ldaushell    ! Subroutine to print 
                                                !   the values of the projectors
                                                !   for LDA+U calculations
      use basis_types, only : nsp               ! Number of different 
                                                !   chemical species
      use basis_types, only : charge            ! Ionic charge to generate the 
                                                !   the basis set
      use atmparams,   only : nrmax             ! Maximum number of points 
                                                !   in the logarithmic grid
      use atmparams,   only : lmaxd             ! Maximum angular momentum 
                                                !   for both orbitals and
                                                !   projectors.
      use atmparams,   only : NTBMAX            ! Maximum number of points 
                                                !   in the tables defining
                                                !   orbitals, projectors and
                                                !   local neutral-atom pseudo
      use pseudopotential, only: pseudopotential_t ! Derived type where all
                                                !   the information about
                                                !   the pseudopotential 
                                                !   is stored
      use atom,        only : schro_eq          ! Subroutine to solve the
                                                !   radial part of the
                                                !   Schrodinger equation
      use atom,        only : rc_vs_e           ! Subroutine to determine
                                                !   the cutoff radius from the 
                                                !   energy shift
      use atom,        only : build_vsoft       ! Subroutine to construct 
                                                !   the soft-confinement potent.
      use atom_options,only: write_ion_plot_files ! Subroutine to plot the 
                                                !   basis functions and other
                                                !   atomic functions
      use atm_types,   only : species_info      ! Derived type with all the info
                                                !   about the radial functions
                                                !   (PAOs, KB projectors, 
                                                !   LDA+U proj,
                                                !   VNA potentials, etc)
                                                !   for a given atomic specie
      use atm_types,   only : species           ! Actual array where the  
                                                !   previous information is 
                                                !   stored
      use atm_types,   only : nspecies          ! Total number of different  
                                                !   atomic species
      use units,       only : pi                ! Value of pi
      use alloc,       only : re_alloc          ! Allocation routines
      use radial                                ! Derived type for the radial
                                                !   functions
      use interpolation, only: spline           ! set spline interpolation
      use interpolation, only: polint           ! polynomial interpolation


      implicit none

      integer :: method_gen_ldau_proj ! Method used to generate the 
                                      !   LDA+U projectors
                                      !   Default value: exact solution 
                                      !   of the pseudoatomic problem 
                                      !   cutted with a Fermi function
      real(dp) :: energy_shift_ldau   ! Energy increase used to define 
                                      !   the localization radious of the LDA+U 
                                      !   projectors (similar to the parameter
                                      !   PAO.EnergyShift)
                                      !   Default value: 0.05 Ry
      real(dp) :: dnrm_rc             ! Parameter used to define the cutoff 
                                      !   radius that enters the 
                                      !   Fermi distribution to cut the 
                                      !   LDA+U projectors. 
                                      !   Only used if method_gen_ldau_proj = 2
                                      !   It is the norm of the original 
                                      !   pseudoatomic orbital contained in 
                                      !   a sphere of radius r_c.
                                      !   Default value: 0.90
      real(dp) :: width_fermi_ldau    ! Parameter used to define the width of 
                                      !   Fermi distribution to cut the 
                                      !   LDA+U projectors. 
                                      !   Only used if method_gen_ldau_proj = 2
                                      !   Default value: 0.05
      real(dp), pointer :: projector(:,:,:) ! Radial part of the LDA+U projector
      integer,  save, public, pointer  ::  nprojsldau(:)
                                      ! Total number of LDA+U projectors
                                      !   (including the different angular 
                                      !   dependencies): i.e. a radial projector
                                      !   with d-character counts as 5 different
                                      !   LDA+U projectors
      integer      :: nrval           ! Actual number of points in the
                                      !   logarithmic grid


      type(basis_def_t), pointer :: basp
      type(ldaushell_t), pointer :: ldau
      type(ldaushell_t), pointer :: lsldau

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

      public :: read_ldau_specs
      public :: ldau_proj_gen
      public :: populate_species_info_ldau

      private

      CONTAINS

! subroutine read_ldau_specs           : Subroutine that reads all the 
!                                        info in the fdf file related with the
!                                        LDA+U projectors and 
!                                        allocate some space for 
!                                        the projector pointer
! subroutine ldau_proj_gen             : Subroutine that solves the 
!                                        Schrodinger equation for the 
!                                        isolated atom and generates the
!                                        LDA+U projectors 
! subroutine fermicutoff               : Subroutine that computes the Fermi
!                                        function used to cut the long 
!                                        atomic wave functions and produce 
!                                        the LDA+U projectors.
!                                        only used if 
!                                        method_gen_ldau_proj = 2
! subroutine populate_species_info_ldau: Subroutine that populates the 
!                                        data structures related with the LDA+U
!                                        projectors in the species 
!                                        derived types.
!                                        Called from the atm_transfer subroutine

!---
      subroutine read_ldau_specs()

      integer :: isp                ! Dummy parameter to account for the 
                                    !   species label
      integer :: ish, jsh           ! Dummy parameters to account for the 
                                    !   loop on shells
      integer :: indexp             ! Dummy parameters to account for the 
                                    !   reading of lines in LDAU.proj block
      integer :: l                  ! Angular quantum number
      integer :: maxnumberproj      ! Maximum number of projectors 
                                    !   considered in a given species

!     Default generation method for the LDA+U projectors
      integer, parameter          :: method_gen_default= 2

!     Default value of the energy-shift to define the cut-off of the LDA+U proj.
      real(dp), parameter         :: energy_shift_ldau_default=0.05_dp

!     Default value of parameter used to define the cutoff radius that enter 
!     the Fermi distribution to produce the LDA+U projectors
!     (only used if method_gen_default= 2)
      real(dp), parameter         :: dnrm_rc_default = 0.90_dp

!     Default value of parameter used to define the width of
!     the Fermi distribution to produce the LDA+U projectors
!     (only used if method_gen_default= 2)
      real(dp), parameter         :: width_fermi_ldau_default = 0.05_dp

!     Default Soft-confinement parameters set by the user
      logical,  save  :: lsoft
      real(dp), save  :: softRc, softPt

!------------------------------------------------------------------------
!     Read the generation method for the LDA+U projectors
      method_gen_ldau_proj = 
     .  fdf_integer('LDAU.ProjectorGenerationMethod',method_gen_default)

!     Read the energy-shift to define the cut-off radius of the LDA+U projectors
      energy_shift_ldau = 
     .  fdf_physical('LDAU.EnergyShift',energy_shift_ldau_default,'Ry')

!     Read the parameter used to define the cutoff radius used in the Fermi
!     distribution 
      dnrm_rc = fdf_double('LDAU.CutoffNorm',dnrm_rc_default)

!     Read information about defaults for soft confinement
      lsoft  = fdf_boolean('PAO.SoftDefault',     .false. )
      softRc = fdf_double( 'PAO.SoftInnerRadius', 0.9d0   )
      softPt = fdf_double( 'PAO.SoftPotential',   40.0d0  )
!     Sanity checks on values
      softRc = max(softRc,0.00d0)
      softRc = min(softRc,0.99d0)
      softPt = abs(softPt)

!     Allocate and initialize the array with the number of projectors per 
!     atomic specie
      nullify( nprojsldau )
      call re_alloc( nprojsldau, 1, nsp, 'nprojsldau', 
     .    'read_ldau_specs' )
      nprojsldau(:) = 0

!     Read the LDAU.proj block
      if (.not. fdf_block('LDAU.proj',bfdf)) RETURN

      do while(fdf_bline(bfdf, pline))     !! over species
        if (.not. fdf_bmatch(pline,'ni'))
     .      call die('Wrong format in LDAU.proj')
        isp = label2species(fdf_bnames(pline,1))
        if (isp .eq. 0) then
          write(6,'(a,1x,a)')
     .      'WRONG species symbol in LDAU.proj:',
     .      trim(fdf_bnames(pline,1))
          call die()
        endif

        basp => basis_parameters(isp)

!       Read on how many orbitals of a given atomic species
!       we are going to apply the U correction
        basp%nldaushells_tmp = fdf_bintegers(pline,1)

!       Allocate space in the derived type basis_parameters
!       to host the information on the atomic orbitals where 
!       the U will be applied
        allocate(basp%tmp_ldaushell(basp%nldaushells_tmp))

!       Loop on all the different orbitals where the U will be applied
        shells: do ish = 1, basp%nldaushells_tmp
          ldau => basp%tmp_ldaushell(ish)
          call initialize(ldau)

          if (.not. fdf_bline(bfdf, pline)) 
     .         call die('Not enough information on the AO in LDAU.proj')

!         Read the principal and angular quantum numbers of the atomic orbital
!         where the U will be applied
!         Also we check what is the maximum value of the angular quantum
!         number between all of them that are read

!         In the LDAU.proj block, the information about the projectors
!         can be given as:
!         n=3    2            # n, l
!         i.e. with a string "n=" and then two integers...
          if (fdf_bmatch(pline,'nii')) then
            ldau%n = fdf_bintegers(pline,1)
            ldau%l = fdf_bintegers(pline,2)
            basp%lmxldaupj = max(basp%lmxldaupj,ldau%l)

!         or deleting the string "n="
!           3    2            # n, l
!         i.e. only with two integers
          elseif (fdf_bmatch(pline,'ii')) then
            ldau%n = fdf_bintegers(pline,1)
            ldau%l = fdf_bintegers(pline,2)
            basp%lmxldaupj = max(basp%lmxldaupj,ldau%l)

!         or only with one integer. In this case, this is the 
!         angular quantum number.
!         If the semicore states is included in the valence, then 
!         this is not valid and the principal quantum number has to be given
!         explicitly
!           2            # l
          elseif (fdf_bmatch(pline,'i')) then
             if (basp%semic) 
     .         call die('Please specify n if there are semicore states')
             ldau%l = fdf_bintegers(pline,1)
             ldau%n = basp%ground_state%n(ldau%l)
             basp%lmxldaupj = max(basp%lmxldaupj,ldau%l)
          else
             call die('Bad format of (n), l line in LDAU.proj')
          endif

!         Check for consistency in the sequence of principal and
!         angular quantum numbers
          do jsh = 1, ish-1
             if( ldau%l .eq. basp%tmp_ldaushell(jsh)%l .and. 
     .           ldau%n .eq. basp%tmp_ldaushell(jsh)%n ) 
     .       call die(
     .        'LDAU projs. with the same l need different values of n')        
          enddo

!         Check whether soft-confinement will be used
          if (fdf_bsearch(pline,'E',indexp)) then
            if (fdf_bmatch(pline,'vv',after=indexp)) then
              ldau%vcte = fdf_bvalues(pline,ind=1,after=indexp)
              ldau%rinn = fdf_bvalues(pline,ind=2,after=indexp)
            else
              call die('Need vcte and rinn after E in LDAU.proj')
            endif
          elseif (lsoft) then
            ldau%vcte = softPt
            ldau%rinn = -softRc
          else
            ldau%vcte = 0.0_dp
            ldau%rinn = 0.0_dp
          endif

!         Read the U and J parameters for this atomic orbital
          if (.not. fdf_bline(bfdf, pline)) 
     .      call die('No information for the U and J parameters...')

          if ( fdf_bnvalues(pline) .ne. 2)
     .      call die('Insert values for the U and J parameters')

          if (fdf_bmatch(pline,'vv')) then
              ldau%u = fdf_bvalues(pline,1)
              ldau%j = fdf_bvalues(pline,2)
          endif

!         Read the cutoff radii (rc) to generate the projectors
!         and the contraction functions (lambda) 
          if ( .not. fdf_bline(bfdf, pline) ) 
     .      call die('No information for the rc for projectors...')

          if( method_gen_ldau_proj .eq. 1 ) then
            if ( fdf_bnvalues(pline) .ne. 1 )
     .        call die('Insert one value for the rc')
            ldau%rc    = fdf_bvalues(pline,1)

          elseif( method_gen_ldau_proj .eq. 2 ) then
            if ( fdf_bnvalues(pline) .ne. 2 )
     .        call die('Insert one value for the rc and width')
            ldau%rc      = fdf_bvalues(pline,1)
            ldau%dnrm_rc = dnrm_rc
            if ( fdf_bvalues(pline,2) .lt. 1.d-4 ) then
              ldau%width   = width_fermi_ldau_default
            else
              ldau%width   = fdf_bvalues(pline,2)
            endif

          endif

!         Optional: read the value for the contraction factor (lambda)
          if ( .not. fdf_bline(bfdf,pline) ) then
             if (ish.ne.basp%nldaushells_tmp)
     .         call die('Not enough shells')
!            Dafault values for the scale factors
          else
            if (.not. fdf_bmatch(pline,'r')) then
              ! New shell or species
              ! Default values for the scale factors
              if ( .not. fdf_bbackspace(bfdf) )
     .          call die('read_ldau_specs: ERROR in LDAU.proj block')
!              cycle shells
            else
              if ( fdf_bnreals(pline) .ne. 1 )
     .          call die('One optional value of lambda')
              ldau%lambda = fdf_breals(pline,1)
            endif
          endif

!         For debugging
          call print_ldaushell(ldau)
!         End debugging

        enddo shells     ! end of loop over shells for species isp


!       Count the total number of projectors
        nprojsldau(isp) = 0
        do ish = 1, basp%nldaushells_tmp
          ldau => basp%tmp_ldaushell(ish)
          l      = ldau%l
          nprojsldau(isp) = nprojsldau(isp) + (2*l + 1)
        enddo
        basp%nldauprojs_lm = nprojsldau(isp)
!!       For debugging
!        write(6,'(a,i5)')'read_ldau_specs: lmxkb     = ', 
!     .    basp%lmxkb
!        write(6,'(a,i5)')'read_ldau_specs: lmxldaupj = ', 
!     .    basp%lmxldaupj
!        write(6,'(a,i5)')'read_ldau_specs: nprojsldau(isp) = ', 
!     .    nprojsldau(isp)
!!       End debugging

      enddo  ! end loop over species

!!     Is this needed??
!!     tmp_ldaushell contains all the required information for the 
!!     LDA+U projectors.
!!     It does not matter if some of them do have the same value of l.
!!     They can be easily distinguished by the principal quantum number.
!!     I think this is not necessary at all.
!!     Even more, I do not think that this will be used ...
!!     I can not see a problem where a U will be required at the same time
!!     for the 3d and 4d. One of them will be completely full...
!!
!!     OK, now classify the states by l-shells
!!
!      do isp = 1, nsp
!        basp => basis_parameters(isp)
!        if ( basp%lmxldaupj .eq. -1 ) cycle !! Species not in block
!
!        allocate (basp%ldaushell(0:basp%lmxldaupj))
!
!        loop_l: do l = 0, basp%lmxldaupj
!          lsldau => basp%ldaushell(l)
!          call initialize(lsldau)
!          lsldau%l = l
!
!!         Search for tmp_shells with given l
!          nn = 0
!          do ish = 1, basp%nldaushells_tmp
!            ldau=>basp%tmp_ldaushell(ish)
!            if (ldau%l .eq. l) nn = nn + 1
!          enddo
!
!        enddo loop_l ! End loop on angular momentum
!
!      enddo  ! End loop on number of species

!     Allocate and initialize the array where the radial part of the 
!     projectors in the logarithmic grid will be stored 
      maxnumberproj = 0
      do isp = 1, nsp
        basp => basis_parameters(isp)
        maxnumberproj = max( maxnumberproj, basp%nldaushells_tmp ) 
      enddo
      nullify( projector )
      call re_alloc( projector, 
     .               1, nsp, 
     .               1, maxnumberproj,
     .               1, nrmax, 
     .               'projector', 'ldau_proj_gen' )
      projector = 0.0_dp


!!     For debugging
!      call die('Testing read_ldau_specs')
!!     End debugging

            
      end subroutine read_ldau_specs

! ---------------------------------------------------------------------
! ---------------------------------------------------------------------

      subroutine ldau_proj_gen( isp )
! ---------------------------------------------------------------------
!     Generation of LDA+U projectors
!     LDA+U projectors are, basically, pseudo-atomic-orbitals 
!     with artificially small radii.
!     Written by D. Sanchez-Portal, Aug. 2008 after module basis_gen
!     Rewritten by Javier Junquera to merge with the top of the trunk
!     (Siesta 4.0), Feb. 2016
! ---------------------------------------------------------------------
      use basis_specs, only : restricted_grid
      use basis_specs, only : rmax_radial_grid

      use siestaXC,    only : setXC
      use siestaXC,    only : atomXC

      integer, intent(in)   :: isp   ! Species index

!     Internal variables
      integer  :: n           ! Principal quantum number of the projector
      integer  :: l           ! Angular quantum number of the projector
      integer  :: lpseudo     ! Angular quantum number of the pseudopotential
      integer  :: iproj       ! Counter for the loop on projectors
      integer  :: ir          ! Counter for the loop on points in the log grid
      integer  :: ndown       ! Counter for the loop on l for the pseudos
      integer  :: nldaupj     ! Number of LDA+U projectors that will be computed
                              !    for a given specie (here we consider only
                              !    different radial parts)
      integer  :: nodd        ! Check whether we have and odd number of points
                              !    in the logarithmic grid
      integer  :: nnodes      ! Number of nodes in the radial part of the 
                              !    eigenfunctions of the Schrodinger equation
      integer  :: nprin       ! Principal quantum number within the pseudoatom
      real(dp) :: U           ! Value of the U parameter
      real(dp) :: J           ! Value of the J parameter
      real(dp) :: r2          ! Square of the distance to the nuclei 
      real(dp) :: rco         ! Cutoff radius
      real(dp) :: rc          ! Cutoff radius (auxiliary variable to fit in an
                              !   odd number of points in the log grid)
      real(dp) :: phi         ! Wave function times r at a given point in
                              !   the log grid
      real(dp) :: lambda      ! Contraction factor
      real(dp) :: el          ! Energy of the eigenvalue after adding the
                              !    energy shift
      real(dp) :: dnorm       ! Norm of the projector
      real(dp) :: rinn        ! Inner radius where the soft-confinement potent.
                              !   starts off
      real(dp) :: vcte        ! Prefactor of the soft-confinement potent.
      real(dp) :: ionic_charge! Ionic charge to generate the basis set.
      logical  :: switch      ! Logical variable that determines whether the 
                              !    the calculation of LDA+U projectors is 
                              !    required for this species
!     Variables used only in the call to atomxc
      real(dp) :: ex          ! Total exchange energy 
      real(dp) :: ec          ! Total correlation energy
      real(dp) :: dx          ! IntegralOf( rho * (eps_x - v_x) )
      real(dp) :: dc          ! IntegralOf( rho * (eps_c - v_c) )

      real(dp) :: eigen(0:lmaxd)      ! Eigenvalues  of the Schrodinger equation
      real(dp) :: rphi(nrmax,0:lmaxd) ! Eigenvectors of the Schrodinger equation
      real(dp) :: vsoft(nrmax)        ! Soft-confinement potential
      real(dp) :: fermi_func(nrmax)   ! Fermi function used to cut the 
                                      !    long pseudowave functions and 
                                      !    produce the LDA+U projectors

!
!     Derived types where some information on the different shells are stored
!

      type(basis_def_t),       pointer :: basp  ! Parameters that define the
                                                !   basis set, KB projectors,
                                                !   LDA+U projectors, pseudopot
                                                !   etc for a given species
      type(ldaushell_t),       pointer :: shell ! Information about 
                                                !   LDA+U projectors
      type(pseudopotential_t), pointer :: vps   ! Pseudopotential information

!
!     Variables related with the radial logarithmic grid
!
      integer      :: nr                     ! Number of points required to
                                             !   store the pseudopotential and
                                             !   the wave functions in the
                                             !   logarithmic grid
                                             !   (directly read from the
                                             !   pseudopotential file)
      integer      :: nrc                    ! Number of points required to 
                                             !   store the pseudowave functions
                                             !   in the logarithmic grid
                                             !   after being strictly confined.
      real(dp)     :: a                      ! Step parameter of log. grid
                                             !   (directly read from the
                                             !   pseudopotential file)
      real(dp)     :: b                      ! Scale parameter of log. grid
                                             !   (directly read from the
                                             !   pseudopotential file)
      real(dp)     :: rofi(nrmax)            ! Radial points of the
                                             !   logarithmic grid
                                             !   rofi(r)=b*[exp(a*(i-1)) - 1]
                                             !   (directly read from the
                                             !   pseudopotential file)
      real(dp)     :: drdi(nrmax)            ! Derivative of the radial
                                             !   distance respect the mesh index
                                             !   Computed after the radial mesh
                                             !    is read
      real(dp)     :: s(nrmax)               ! Metric array
                                             !   Computed after the radial mesh
                                             !    is read
      real(dp)     :: rpb, ea                ! Local variables used in the
                                             !   calculation of the log. grid

!
!     Variable used to store the semilocal component of the pseudopotential 
!
!
      character*4  ::  nicore                ! Flag that determines whether
                                             !   non-linear core corrections
                                             !   are included
      character*3  ::  irel                  ! Flag that determines whether
                                             !   the atomic calculation is
                                             !   relativistic or not
      real(dp)     :: vpseudo(nrmax,0:lmaxd) ! Semilocal components of the
                                             !   pseudopotentials
                                             !   (directly read from the
                                             !   pseudopotential file)
      real(dp)     :: zval                   ! Valence charge of the atom
                                             !   (directly read from the
                                             !   pseudopotential file)


!
!     Variable used to store the semilocal component of the pseudopotential 
!
      real(dp)                 :: chgvps     ! Valence charge of the pseudoion
                                             !   for which the pseudo was
                                             !   generated in the ATM code
                                             !   (it might not correspond with
                                             !   the nominal valence charge
                                             !   of the atom if the pseudo
                                             !   has been generated for an ionic
                                             !   configuration, for instance
                                             !   when semicore has been
                                             !   explicitly included in the
                                             !   valence).
                                             !   For instance, for Ba with
                                             !   the semicore in valence,
                                             !   (atomic reference configuration
                                             !   5s2 5p6 5d0 4f0),
                                             !   chgvps = 8  (two in the 5s
                                             !                and six in the 5p)
                                             !   zval   = 10 (two in the 5s,
                                             !                six in the 5p,
                                             !                and two in the 6s.
                                             !   These two last not included in
                                             !   reference atomic configuration)
      real(dp)     :: rho(nrmax)             ! Valence charge density
                                             !   As read from the pseudo file,
                                             !   it is angularly integrated
                                             !   (i.e. multiplied by 4*pi*r^2).
      real(dp)     :: rho_PAO(nrmax)         ! Valence charge density
      real(dp)     :: ve(nrmax)              ! Electrostatic potential
                                             !   generated by the valence charge
                                             !   density, readed from the
                                             !   pseudo file
      real(dp)     :: vePAO(nrmax)           ! Electrostatic potential 
                                             !   generated by the "scaled" 
                                             !   valence charge density
      real(dp)     :: vePAOsoft(nrmax)       ! vePAO + the soft-confinement pot.
      real(dp)     :: vxc(nrmax)             ! Exchange and correlation potentil
      real(dp)     :: chcore(nrmax)          ! Core charge density
                                             !   As read from the pseudo file,
                                             !   it is angularly integrated
                                             !   (i.e. multiplied by 4*pi*r^2).
      real(dp)     :: auxrho(nrmax)          !  Sum of the valence charge and
                                             !   core charge (if NLCC included)
                                             !   densities to compute the
                                             !   atomic exchange and correl.
                                             !   potential.
                                             !   auxrho is NOT angularly integr.
                                             !   (not multiplied by 4*pi*r^2)
      integer      :: irelt                  ! Flag that determines whether the
                                             !   atomic calculation to
                                             !   generate the pseudopotential
                                             !   was relativistic (irelt = 1)
                                             !   or no relativistic (irelt = 0)


!     Associate the pointer so it points to the variable where all the
!     parameters defining the basis sets of the given species are stored
      basp => basis_parameters(isp)

!     Determine if something has to be done regarding the 
!     generation of the LDA+U projectors.
!     If LDA+U is not required (number of LDA+U projectors equal to zero), 
!     then do nothing and return.

!     Compute how many LDA+U projector we are going to compute
!     for this species
      nldaupj = basp%nldaushells_tmp

!     Determine whether the calculation of LDA+U projectors is required or not
!     for this atomic species
      switch=.false.
      if( nldaupj .gt. 0 ) switch=.true.
 
      if( .not. switch ) return 

!     Associate the pointer so it points to the variable where all the
!     parameters defining the basis sets of the given species are stored
      vps => basp%pseudopotential 

!
!     Read all the required information from the pseudopotentials that
!     will be required to solve the Schrodinger equation for the isolated atoms
!
      nr     = vps%nr
      b      = vps%b
      a      = vps%a
      zval   = vps%zval
      nicore = vps%nicore
      irel   = vps%irel
      ionic_charge = charge(isp)

      nrval = nr + 1
      if (rmax_radial_grid /= 0.0_dp) then
         nrval = nint(log(rmax_radial_grid/b+1.0d0)/a)+1
         write(6,"(a,f10.5,i5)")
     .     'Maximum radius (at nrval) set to ',
     .     rmax_radial_grid, nrval
      endif

      if (restricted_grid) then
        nodd  = mod(nrval,2)
        nrval = nrval -1 + nodd ! Will be less than or equal to vp%nrval
      endif

      if ( nrval .gt. nrmax ) then
        write(6,'(a,i4)')
     .   'ldau_proj_gen: ERROR: Nrmax must be increased to at least',
     .    nrval
        call die
      endif

!     Read the radial logarithmic mesh
      rofi(1:nrval) = vps%r(1:nrval)

!     Calculate drdi and s
!     drdi is the derivative of the radial distance respect to the mesh index
!     i.e. rofi(ir)= b*[ exp( a*(i-1) ) - 1 ] and therefore
!     drdi=dr/di =a*b*exp(a*(i-1))= a*[rofi(ir)+b]

      rpb = b
      ea  = dexp(a)
      do ir = 1, nrval
        drdi(ir) = a * rpb
        s(ir)    = dsqrt( a * rpb )
        rpb      = rpb * ea
      enddo

!!     For debugging
!!     Differences with respect Daniel's grid implementation
!!     can appear in the number of points nrval.
!!     In this latest version, the option restricted_grid is activated
!!     by default.
!!     This was not yet implemented in the version where Daniel started
!!     Even more, the value of s printed by Daniel corresponds
!!     with the redefinition given below for the integration of the
!!     Schrodinger equation (s = drdi^2).
!!     The one defined here is required for the solution of the Poisson equation
!      do ir = 1, nrval
!        write(6,'(i5,3f20.12)') ir, rofi(ir), drdi(ir), s(ir)
!      enddo
!      call die()
!!     End debugging

!       
!     Read the ionic pseudopotentials (Only 'down' used)
!     These are required to solve the Schrodinger equation for the isolated
!     atoms.
!     Here we read all the semilocal components of the pseudopotential
!     independently of whether the LDA+U projector for a particular
!     angular momentum is required or not.
!
      do 20 ndown = 1, basp%lmxldaupj+1

        lpseudo = vps%ldown(ndown)

        if ( lpseudo .ne. ndown-1 ) then
           write(6,'(a)')
     . 'ldau_proj_gen: Unexpected angular momentum  for pseudopotential'      
           write(6,'(a)')
     . 'ldau_proj_gen: Pseudopot. should be ordered by increasing l'      
        endif

        vpseudo(1:nrval,lpseudo) = vps%vdown(ndown,1:nrval)

        do ir = 2, nrval
          vpseudo(ir,lpseudo) = vpseudo(ir,lpseudo)/rofi(ir)
        enddo
        vpseudo(1,lpseudo) = vpseudo(2,lpseudo)     ! AG

  20  enddo

!!     For debugging
!!     Up to this point, these are the same pseudos as read in
!!     Daniel's version of LDA+U
!!     The only difference might be at the number of points in
!!     the log grid
!      do lpseudo = 0, basp%lmxldaupj 
!        write(6,'(/a,i5)')
!     .    ' ldau_proj_gen: Reading pseudopotential for l = ',
!     .    lpseudo
!
!        do ir = 1, nrval
!          write(6,'(a,i5,2f20.12)')
!     .      ' ir, rofi, vpseudo = ', ir, rofi(ir), vpseudo(ir,lpseudo)
!        enddo
!      enddo
!!     End debugging

!     Read the valence charge density from the pseudo file
!     and scale it if the ionic charge of the reference configuration
!     is not the same as the nominal valence charge of the atom
      chgvps = vps%gen_zval
      do ir = 1, nrval
        rho(ir) = chgvps * vps%chval(ir)/zval
      enddo

!     Find the Hartree potential created by a radial electron density
!     using the Numerov's method to integrate the radial Poisson equation.
!     The input charge density at this point has to be angularly integrated.
      call vhrtre( rho, ve, rofi, drdi, s, nrval, a )

!!     For debugging
!      do ir = 1, nrval
!        write(6,'(a,i5,3f20.12)')
!     .    ' ir, rofi, rho, ve = ', 
!     .      ir, rofi(ir), rho(ir), ve(ir)
!      enddo
!      call die()
!!     End debugging

!      Set 'charge':
!      1. If 'charge' is not set in the fdf file
!         then set it to zval-chgvps.
!      2. If 'charge' is equal to zval-chgvps, set it to that.
!
       if( ( abs(ionic_charge) .eq. huge(1.0_dp) ) .or.
     .     ( abs( ionic_charge-(zval-chgvps) ) .lt. 1.0d-3) ) then
        ionic_charge = zval - chgvps
       endif

!      For LDA+U projector calculations
!      We use the "scaled" charge density of an ion of total charge "charge"
!      As seen above, this ion could be the one involved in ps generation,
!      or another specified at the time of basis generation.
!      Example: Ba: ps ionic charge: +2
!               basis gen specified charge: +0.7

       do ir = 2,nrval
         rho_PAO(ir) = (zval-ionic_charge) * rho(ir) / chgvps
       enddo
       call vhrtre( rho_PAO, vePAO, rofi, drdi, s, nrval, a )

!     Read the core charge density from the pseudo file
      chcore(1:nrval) = vps%chcore(1:nrval)

!     Compute the exchange and correlation potential in the atom
!     Note that atomxc expects true rho(r), not 4 * pi * r^2 * rho(r)
!     We use auxrho for that
!
      do ir = 2, nrval
        r2 = rofi(ir)**2
        r2 = 4.0_dp * pi * r2
        dc = rho(ir) / r2
        if( nicore .ne. 'nc  ')  dc = dc + chcore(ir) / r2
        auxrho(ir) = dc
      enddo

      r2        = rofi(2) / (rofi(3)-rofi(2))
      auxrho(1) = auxrho(2) - ( auxrho(3) - auxrho(2) ) * r2

!     Determine whether the atomic calculation to generate the pseudopotential
!     is relativistic or not
      if (irel.eq.'rel') irelt=1
      if (irel.ne.'rel') irelt=0

!     Set the exchange and correlation functional
      call setXC( 1, (/'GGA'/), (/'PBE'/), (/1._dp/), (/1._dp/) )

!     Compute the exchange and correlation potential
      call atomxc( irelt, nrval, nrmax, rofi,
     .             1, auxrho, ex, ec, dx, dc, vxc )

!!     For debugging
!      write(6,'(a,i5)') 'irelt = ', irelt
!      write(6,'(a,i5)') 'nrval = ', nrval
!      write(6,'(a,i5)') 'nrmax = ', nrmax
!      do ir = 1, nrval
!        write(6,'(a,i5,3f20.12)')
!     .    ' ir, rofi, auxrho, vxc = ', 
!     .      ir, rofi(ir), auxrho(ir), vxc(ir)
!      enddo
!      call die()
!!     End debugging

!     Add the exchange and correlation potential to the Hartree potential
      ve(1:nrval) = ve(1:nrval) + vxc(1:nrval) 

!!     For debugging
!      write(6,'(a,f20.12)')' chg = ', chgvps
!      write(6,'(a,f20.12)')' a   = ', a
!      write(6,'(a,f20.12)')' b   = ', b
!      do ir = 1, nrval
!        write(6,'(a,i5,3f20.12)') 
!     .      ' ir, rofi, vxc+vhr, vpseudo = ', 
!     .        ir, rofi(ir), ve(ir), vpseudo(ir,0)
!      enddo
!      call die()
!!     End debugging

      do ir = 2,nrval
        r2 = rofi(ir)**2
        r2 = 4.0d0*pi*r2
        dc = rho_PAO(ir)/r2
        if (nicore.ne.'nc  ') dc = dc + chcore(ir)/r2
        auxrho(ir) = dc
      enddo

      r2 = rofi(2)/(rofi(3)-rofi(2))
      auxrho(1) = auxrho(2) -(auxrho(3)-auxrho(2))*r2

      call atomxc( irelt, nrval, nrmax, rofi,
     .             1, auxrho, ex, ec, dx, dc, vxc )

      vePAO(1:nrval) = vePAO(1:nrval) + vxc(1:nrval)

!
!     Redefine the array s for the Schrodinger equation integration
!
      s(2:nrval) = drdi(2:nrval) * drdi(2:nrval)
      s(1) = s(2)

!     Loop over all the projectors that will be generated
      loop_projectors: do iproj = 1, nldaupj
         shell => basp%tmp_ldaushell(iproj)
         n      = shell%n
         l      = shell%l
         U      = shell%u
         J      = shell%j
         rco    = shell%rc
         rinn   = shell%rinn
         vcte   = shell%vcte

!        If the compression factor is negative or zero,
!        the orbitals are left untouched
         if( shell%lambda .le. 0.0d0 ) shell%lambda=1.0d0
         lambda = shell%lambda

!        Check whether the cutoff radius for the LDA+U projector
!        is explicitly determined in the input file or automatically
!        controlled by 
!        - the EnergyShift parameter            (method_gen_ldau_proj = 1)
!        - the cutoff of the Fermi distribution (method_gen_ldau_proj = 2)
!        If the latter is the case, then we need to solve the
!        Schrodinger equation for the isolated atom
         if ( rco .lt. 1.0d-5 ) then

!          Determine the number of nodes in the radial part of
!          the eigenfunction
!          THIS HAS TO BE UPDATED WITH THE SUBROUTINES
!          OF THE NEW PSEUDOS:
!          FROM THE KNOWLEDGE OF n AND l, IT SHOULD BE POSSIBLE
!          TO DETERMINE THE NUMBER OF NODES
           nnodes = 1

!          Determine the principal quantum number within the pseudoatom
           nprin  = l + 1

!!          For debugging
!           write(6,'(/a,i2)')  
!     .       'LDAUprojs with principal quantum number n = ', n
!           write(6,'(a,i2)')   'LDAUprojs with angular momentum l = ', l
!           write(6,'(a,f12.5)')'LDAUprojs with U        = ',U
!           write(6,'(a,f12.5)')'LDAUprojs with J        = ',J
!           write(6,'(a,f12.5)')'LDAUprojs with lambda   = ',shell%lambda
!           write(6,'(a,f12.5)')'LDAUprojs with rc       = ',shell%rc
!           write(6,'(a,i5)')   'LDAUprojs with nnodes   = ',nnodes
!           write(6,'(a,i5)')   'LDAUprojs with nprin    = ',nprin
!           write(6,'(a,i5)')   'LDAUprojs with nrval    = ',nrval
!           write(6,'(a,i5)')   'LDAUprojs with nrmax    = ',nrmax
!           write(6,'(a,f12.5)')'LDAUprojs with zval     = ',zval
!           write(6,'(a,f12.5)')'LDAUprojs with a        = ',a
!           write(6,'(a,f12.5)')'LDAUprojs with b        = ',b
!!          End debugging

!          Initialize the eigenfunctions
           rphi(:,l) = 0.0_dp

           call schro_eq( zval, rofi, vpseudo(1,l), ve, s, drdi,
     .                    nrval, l, a, b, nnodes, nprin,
     .                    eigen(l), rphi(1,l) )

!!          For debugging
!           write(6,'(/a,i5)')  '# l = '          , l 
!           write(6,'(a,f12.5)')'# Eigenvalue =  ', eigen(l)
!           write(6,'(a)')      '# Eigenfunction '
!           do ir = 1, nrval
!             write(6,'(2f20.12)')rofi(ir), rphi(ir,l)
!           enddo 
!!          End debugging

!           Cutoff controled by the energy shift parameter:
            if( method_gen_ldau_proj .eq. 1) then
!
!             Compute the cutoff radius of the LDA+U projectors 
!             as given by energy_shift_ldau
!
              if( eigen(l) .gt. 0.0_dp ) then
                write(6,'(/a,i2,a)')
     .          'ldau_proj_gen: ERROR Orbital with angular momentum L=',
     .          l, ' not bound in the atom'
                write(6,'(a)')
     .          'ldau_proj_gen: an rc  radius must be explicitly given'     
                call die()
              endif

              if( abs(energy_shift_ldau) .gt. 1.0d-5 ) then
                el = eigen(l) + energy_shift_ldau
                call rc_vs_e( a, b, rofi, vpseudo(1,l), ve, nrval, l, 
     .                        el, nnodes, rco )
              else
                rco = rofi(nrval-2)
              endif

!             Store the new variable for the cutoff radii 
!             automatically determined
              shell%rc = rco 

              write(6,'(/,a,/,a,f10.6,a)')
     .          'ldau_proj_gen: PAO cut-off radius determined from an',
     .          'ldau_proj_gen: energy shift =',energy_shift_ldau,' Ry'
              write(6,'(a,f10.6,a)')
     .          'ldau_proj_gen: rco =',rco,' Bohr'

!           Cutoff controled by the Fermi distribution
            else if( method_gen_ldau_proj .eq. 2) then
              call fermicutoff( nrmax, nrval, rofi, drdi, 
     .                          rphi(:,l), shell, fermi_func )
            endif 

         endif     ! End if automatic determination of the rc

!        At this point, independently of the method,
!        we should now the cutoff radius of the LDA+U projector.
!        Now, we compute it
         rco = shell%rc

!        Store the radial point of the logarithmic grid where the
!        LDA+U projector vanishes
         nrc = nint(log(rco/b+1.0_dp)/a)+1
         shell%nrc = nrc

         if( method_gen_ldau_proj .eq. 1) then
!          Build the soft confinement potential
           vsoft = 0.0_dp
!          Scale the orbitals with the contraction factor
           rc  = rco / lambda
           call build_vsoft( isp, l, 1, rinn, vcte,
     .                       0.0_dp, 0.0_dp, 0.0_dp,
     .                       a, b, rc, rofi, nrval,
     .                       vsoft, plot=write_ion_plot_files )
!!          For debugging
!           write(6,'(/a,i5)')  '# l = '          , l 
!           write(6,'(a,f12.5)')'# Eigenvalue =  ', eigen(l)
!           write(6,'(a)')      '# Soft-confinement      '
!           write(6,'(a,f12.5)')'# Inner radius  = ' , rinn
!           write(6,'(a,f12.5)')'# Prefactor     =  ', vcte
!           write(6,'(a,f12.5)')'# Cutoff radius =  ', rco
!           do ir = 1, nrval
!             write(6,'(2f20.12)')rofi(ir), vsoft(ir)
!           enddo 
!!          End debugging

           do ir = 1, nrval
             vePAOsoft(ir) = vePAO(ir) + vsoft(ir)
           enddo

!
!          If rc is negative, treat it as a fractional value

           if (rco .lt. 0.0_dp) then
             call die("rc < 0 for first-zeta orbital")
           endif


!          Find the required number of points in the logarithmic grid
!          to solve the Scrodingcer equation
           nrc = nint(log(rc/b+1.0_dp)/a)+1

!          Note that rco is redefined here, to make it fall on an odd-numbered
!          grid point.
!
           if (restricted_grid) then
             nodd = mod(nrc,2)
             if( nodd .eq. 0 ) then
               nrc = nrc + 1
             endif
           endif

           rc  = b*(exp(a*(nrc-1))-1.0d0)
           rco = rc * lambda

!          Solve the Schrodinger equation for the required cutoff
!          and with the Hartree potential from the scaled charge density

!          Determine the number of nodes
           nnodes = 1
!          Determine the principal quantum number within the pseudoatom
           nprin  = l + 1

!          Initialize the eigenfunctions
           rphi(:,l) = 0.0_dp
           call schro_eq( zval, rofi, vpseudo(1,l), vePAOsoft, s, drdi,
     .                    nrc, l, a, b, nnodes, nprin,
     .                    eigen(l), rphi(1,l) )

!          Normalize the eigenfunctions
!          and divide them by r^(l+1)
!          In the previous subroutine, we compute r * phi,
!          where phi is the radial part of the wave functions.
!          In Siesta, we store in the tables phi/r^l.
!          Therefore, we need to divide the previous solution by 
!          r^(l+1)
           dnorm = 0.0_dp
           do ir = 2, nrc
             phi  = rphi(ir,l)
             dnorm = dnorm + drdi(ir) * phi * phi
             projector(isp,iproj,ir)=rphi(ir,l)/(rofi(ir)**(l+1))
           enddo
           projector(isp,iproj,1)=projector(isp,iproj,2)

         else if( method_gen_ldau_proj .eq. 2) then
           do ir = 1, nrval
             projector(isp,iproj,ir) = fermi_func(ir) * rphi(ir,l)
           enddo 

!          Normalize the projector
           dnorm = 0.0_dp
           do ir = 1, nrval
!            Here we have computed  r*projector, where projector is the 
!            radial part of the LDA+U projector.
!            To compute the norm in spherical coordinates,
!            we have to integrate \int r^{2} R^{2} dr,
!            and this implies just to take projector**2
             dnorm = dnorm + drdi(ir) * projector(isp,iproj,ir)**2
           enddo
           dnorm = dsqrt(dnorm)
           projector(isp,iproj,:) = projector(isp,iproj,:) / dnorm

!          To store the projector in the radial tables, 
!          we need R/r^l, where R is the radial part of the projector.
!          Since up to this point we have r*R in the array projector,
!          we have to divide it by r^(l+1)
           do ir = 2, nrval
             projector(isp,iproj,ir) = projector(isp,iproj,ir) / 
     .                                 rofi(ir)**(l+1)
           enddo 
           projector(isp,iproj,1) = projector(isp,iproj,2)

         endif     
        
!        For debugging
         dnorm = 0.0_dp
         do ir = 1, nrval
!          The projector that has been stored in the array projector
!          is written in the same format as the atomic orbitals in the
!          inners of Siesta, i. e., in the format of R/r^l,
!          where R is the radial part of the projector.
!          To check if it is normalized,
!          we have to integrate \int r^{2} R^{2} dr,
!          and this implies just to take projector**2
!          and multiply by r^(2l+2) = r^(2*(l+1))
           dnorm = dnorm + drdi(ir)*(projector(isp,iproj,ir)**2)*
     .             rofi(ir)**(2*(l+1))
         enddo

         write(6,'(/a,i5)')  '# l = '          , l 
         write(6,'(a,f12.5)')'# Eigenvalue =  ', eigen(l)
         write(6,'(a)')      '# Projector      '
         write(6,'(a,f12.5)')'# Norm       =  ', dnorm
         do ir = 1, nrval
           write(6,'(2f20.12)')rofi(ir), projector(isp,iproj,ir)
         enddo 
!        End debugging

      enddo loop_projectors     ! End the loop on projectors

!!     For debugging
!      call die("Testing ldau_proj_gen")
!!     End debugging


      end subroutine ldau_proj_gen
! ---------------------------------------------------------------------

      subroutine fermicutoff( nrmax, nrval, rofi, drdi, rphi,
     .                        ldaushell, fermi_func )
!
! This subroutine defines the fermi function used to cut the long 
! atomic wave functions and produce the LDA+U projectors
! Only used if method_gen_ldau_proj = 2
!

      integer,          intent(in)     :: nrmax       ! Parameter of the
      integer,          intent(in)     :: nrval       ! Parameter of the
      real(dp),         intent(in)     :: rofi(nrmax) ! Parameter of the
      real(dp),         intent(in)     :: drdi(nrmax) ! Parameter of the
      real(dp),         intent(in)     :: rphi(nrmax) ! Eigenvectors 
                                                      !   of the Schrodinger
                                                      !   equation
      type(ldaushell_t),intent(inout)  :: ldaushell
      real(dp),         intent(out)    :: fermi_func(nrmax)  ! Fermi function

!     Internal vars
      integer               :: ir      ! Counter for the loops on real space
                                       !   grids
      integer               :: l       ! Angular momentum of the shell
      real(dp)              :: rc      ! "Fermi energy" of the Fermi function
      real(dp)              :: width   ! Width of the Fermi function
      real(dp)              :: a       ! Auxiliary function to compute the
                                       !   Fermi function
      real(dp)              :: dnorm   ! Norm of the original pseudoatomic
                                       !   wave function
      real(dp), parameter   :: gexp = 60.0_dp
      real(dp), parameter   :: eps  = 1.0e-4_dp  ! A small value (epsilon)
                                                 !    for comparison

!     Initialize the angular momentum quantum number.
      l   = ldaushell%l

!     If no cutoff distance is explicitly given in the input file 
!     (LDAU.proj block) then compute the cutoff distance for the Fermi function
!     For this, we have to check at which radial distance
!     the norm of the original pseudo atomic orbital equals 
!     the value introduced in LDAU.CutoffNorm

      if ( ldaushell%rc .lt. eps ) then

         dnorm = 0.0_dp
         do ir = 1, nrmax
!          In rphi we have r*phi, where phi is the radial part of the
!          wave function.
!          To compute the norm in spherical coordinates,
!          we have to integrate \int r^{2} R^{2} dr,
!          and this implies just to take rphi**2
           dnorm = dnorm + drdi(ir) * rphi(ir)**2
           if( dnorm .gt. dnrm_rc ) exit
         enddo
         ldaushell%rc = rofi(ir)

      endif 

!     Initialize Fermi function
      fermi_func = 0.0_dp

!     Determine the parameters of the Fermi distribution
      rc    = ldaushell%rc
      width = ldaushell%width

      do ir = 1, nrval
        a = ( rofi(ir) - rc ) / width
        if( a .lt. -gexp ) then
          fermi_func(ir) = 1.0_dp
        else if( a .gt. gexp ) then
          fermi_func(ir) = 0.0_dp
        else
          fermi_func(ir) = 1.0_dp / ( 1.0_dp+dexp(a) )
        endif
      enddo

!!     For debugging
!      write(6,'(a,f12.5)')'# Fermi function computed with rc = ', rc
!      write(6,'(a,f12.5)')'#  and width  = ', width
!      do ir = 1, nrval
!        write(6,'(2f20.12)')
!     .    rofi(ir), fermi_func(ir)
!      enddo
!      call die()
!!     End debugging

      end subroutine fermicutoff

! ----------------------------------------------------------------------
      subroutine populate_species_info_ldau
!
!     In this subroutine, we populate the variables in the species_info
!     derived type related with the LDA+U projectors.
!     It is called from atm_transfer. 
!
      type(species_info),      pointer :: spp
      type(basis_def_t),       pointer :: basp
      type(ldaushell_t),       pointer :: ldaushell
      type(rad_func),          pointer :: pp
      type(pseudopotential_t), pointer :: vps   

!     Internal variables
      integer  :: is      ! Counter for the loop on atomic species
      integer  :: iproj   ! Counter for the loop on projectors
      integer  :: ir      ! Counter for the loop on real space points
      integer  :: l       ! Quantum angular momentum of a given LDA+U proj.
      integer  :: im      ! Counter for the loop magnetic quantum number
      integer  :: imcount ! 
      integer  :: nr      ! Point in the log. grid closest to the linear grid
      integer  :: nn      ! Total number of points in the log grid considered
                          !   for the interpolation
      integer  :: nmin    ! nr - npoint (see below for the meaning of npoint)
      integer  :: nmax    ! nr + npoint (see below for the meaning of npoint)
      real(dp) :: rc      ! Cutoff radius of the different LDA+U proj.
      integer  :: nrc     ! Point in the log. grid where the LDA+U proj. vanish
      real(dp) :: delta   ! Interval between consecutive points in the grid
                          !   where the LDA+U projectors are stored
      real(dp) :: rpoint  ! Coordinate of the real space points
      real(dp) :: projint ! Interpolated value of the LDA+U projector at rpoint
      real(dp) :: dy      ! Function derivative at point rpoint
      real(dp) :: a       ! Parameters of the logarithmic grid
      real(dp) :: b       ! Parameters of the logarithmic grid
      real(dp) :: yp1     ! First derivative at the first point of the grid
      real(dp) :: ypn     ! First derivative at the last point of the grid
      real(dp) :: rofi(nrmax)          ! Radial points of the
                                       !   logarithmic grid
                                       !   rofi(r)=b*[exp(a*(i-1)) - 1]
                                       !   (directly read from the
                                       !   pseudopotential file)
      real(dp) :: projinputint(nrmax)  ! Radial part of the projector that
                                       !   enters the interpolation routines

      integer, parameter  :: npoint = 4  ! Number of points used by polint 
                                         !    for the interpolation

!     Loop on different atomic species
      loop_species: do is = 1, nspecies
        spp  => species(is)
        basp => basis_parameters(is)
        vps  => basp%pseudopotential 

!       Read the parameters for the logarithmic grid
        a = vps%a
        b = vps%b
!       Read the radial logarithmic mesh
        rofi(1:nrval) = vps%r(1:nrval)

!       Store the total number of LDA+U projectors 
!       counting the "m copies"
!       (including the (2l + 1) factor for each l).
        spp%nprojsldau = nprojsldau(is)

!       Number of LDA+U projectors 
!       not counting the "m copies"
        spp%n_pjldaunl = basp%nldaushells_tmp

!       Store the maximum angular momentum of the LDA+U projectors
!       for each atomic specie
        spp%lmax_ldau_projs = basp%lmxldaupj

!       Loop on all the projectors for a given specie
!       This loop is done only on the different radial shapes,
!       without considering the (2l + 1) possible angular dependencies
        imcount = 0
        loop_projectors: do iproj = 1, spp%n_pjldaunl
          ldaushell => basp%tmp_ldaushell(iproj)

          spp%pjldaunl_n(iproj) = 1
          spp%pjldaunl_l(iproj) = ldaushell%l
          l = spp%pjldaunl_l(iproj)

          do im = -l, l
            imcount = imcount + 1
            spp%pjldau_n(imcount)     = ldaushell%n
            spp%pjldau_l(imcount)     = ldaushell%l
            spp%pjldau_m(imcount)     = im
            spp%pjldau_index(imcount) = iproj
          enddo 
        enddo loop_projectors ! End loop on projectors for a given specie
        if( imcount .ne. spp%nprojsldau ) call die('LDA+U indexing...')

!       Allocate the derived types pjldau, of radial kind,
!       where the radial components of the LDA+U projectors will be stored
!       There will be as many radial functions of this kind
!       as different LDA+U projectors, without including the m copies.
        allocate ( spp%pjldau(spp%n_pjldaunl) )

        do iproj = 1, spp%n_pjldaunl
          ldaushell => basp%tmp_ldaushell(iproj)
          pp => spp%pjldau(iproj)
          call rad_alloc(pp,NTBMAX)
          rc        = ldaushell%rc
          nrc       = ldaushell%nrc
          delta     = rc/(dble(ntbmax-1)+1.0d-20)
          pp%cutoff = rc 
          pp%delta  = delta

          projinputint(:) = projector(is,iproj,:)

!         Interpolate the projectors from the logarithmic grid to the 
!         linear grid
          do ir = 1, ntbmax-1
            rpoint = delta * (ir-1)
            nr     = nint(log(rpoint/b+1.0d0)/a)+1
            nmin   = max( 1,   nr-npoint )
            nmax   = min( nrc, nr+npoint )
            nn     = nmax - nmin + 1
            call polint( rofi(nmin), projinputint(nmin), 
     .                   nn, rpoint, projint, dy )
            pp%f(ir) = projint
          enddo

!         Compute the second derivative of the projectors 
          call rad_setup_d2(pp,yp1=0.0_dp,ypn=huge(1.0_dp))

        enddo

      enddo loop_species ! End loop on atomic species

!     For debugging
      do is = 1, nspecies
        write(6,'(/a,i5)')
     .    '# populate_species_info_ldau: specie number              = ',
     .    is

        spp  => species(is)

        write(6,'(a,i5)')
     .    '#populate_species_info_ldau: specie, spp%lmax_ldau_projs = ',
     .    spp%lmax_ldau_projs
        write(6,'(a,i5)')
     .    '#populate_species_info_ldau: specie, spp%n_pjldaunl      = ',
     .    spp%n_pjldaunl
        write(6,'(a)')
     .    '#populate_species_info_ldau: Loop over different projectors'
        write(6,'(a)')
     .    '#populate_species_info_ldau: not considering m copies '
        write(6,'(a)')
     .    '#populate_species_info_ldau: iproj, pjldau_n, pjldaunl_l'

        do iproj = 1, spp%n_pjldaunl
          write(6,'(a,3i5)')
     .      '#populate_species_info_ldau:', 
     .       iproj, spp%pjldaunl_n(iproj), spp%pjldaunl_l(iproj)
          pp => spp%pjldau(iproj)
          write(6,'(a,f20.12)')
     .      '#populate_species_info_ldau: cutoff = ', pp%cutoff
          write(6,'(a,f20.12)')
     .      '#populate_species_info_ldau: delta  = ', pp%delta
          do ir = 1, ntbmax-1
            rpoint = pp%delta * (ir-1)
            write(6,'(3f20.12)') rpoint, pp%f(ir), pp%d2(ir)
          enddo
          write(6,*) 
        enddo

        write(6,'(a,i5)')
     .    '#populate_species_info_ldau: specie, spp%nprojsldau      = ',
     .    spp%nprojsldau 
        write(6,'(a)')
     .    '#populate_species_info_ldau: Loop over different projectors'
        write(6,'(a)')
     .    '#populate_species_info_ldau: considering m copies '
        write(6,'(a)')
     .    '#populate_species_info_ldau: iproj, pjldau_n, pjldaunl_l'
        write(6,'(a)')
     .    '#populate_species_info_ldau: pjldau_n, l , m, index'
        do iproj = 1, spp%nprojsldau
          write(6,'(4i5)')
     .     spp%pjldau_n(iproj), spp%pjldau_l(iproj), 
     .     spp%pjldau_m(iproj), spp%pjldau_index(iproj) 
        enddo 
      enddo 
!      call die('End testing populate_species_info_ldau')
!     End debugging
      

      end subroutine populate_species_info_ldau

      end module ldau_specs
