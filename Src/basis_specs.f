      module basis_specs
! 
! Alberto Garcia, August 2000, based on 'classic' redbasis.
! 
! Processes the basis information in an fdf file and populates
! the "basis specifications" data structures. This new version
! makes use of a new "parse" module for increased clarity. 
! 
! Here is a guide to the behavior of the main routine "read_basis_specs":
! 
! * Find out the number of species
! * Allocate storage for the data structures
! * Determine any "global" basis specification parameters:
!   - basis_size   (sz, dz, dzp, etc)
!   - basis_type   (split, nodes, etc) ("USER" is no longer valid)
! * Process the mandatory "Chemical-species-label" block to read the
!   species labels (not the same as chemical symbols) and atomic numbers.
!   (note that negative atomic numbers have a special meaning)
! * Assign default values to 'basis_size', 'basis_type'.
! * Assign default values to the 'mass' based on the atomic number.
! * Set up the information about the ground state of the atom.
! * Read information about the valence charge density in the
!   pseudopotential file and determine whether there are semicore
!   states. (Note that in this case  the user is explicitly asked 
!   (see below) to input the 'n' quantum numbers in the PAO.Basis block.)
! * Read the optional fdf blocks:
!   AtomicMass  (routine remass)
!   PAO.BasisSize - Overrides the default 'basis_size' on a per-species basis.
!   PAO.Basis - This is the most complex block, very flexible but in
!               need  of spelling-out the specific restrictions. Line by
!               line, these are:
! 
!   1st:   Species_label number_of_l_shells [basis_type] [ionic_charge] 
! 
!   For each l_shell:
!          [n= n] l nzeta [P [ nzeta_pol]] [E vcte rinn]
!   where 'n= n' is *mandatory* if the species has any semicore states,
!   and the 'P' (polarization) and 'E' (soft confinement potential)
!   sections are optional and can appear in any order after nzeta. 
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
!   After processing PAO.Basis (if it exists), whatever PAO information
!   which is not already available is determined in routine 'autobasis', 
!   using the following defaults:
! 
!   (There is a first check to make sure that species with semicore
!   states, or with Bessel floating orbitals,  have been included in 
!   the PAO.Basis block)
!   
!   The maximum angular momentum for the basis (lmxo) (excluding any
!   polarization orbitals) is that of the ground state, as returned by
!   routine 'lmxofz' from Z.
! 
!   Each l-shell contains just one shell of PAOS, with n set
!   to the appropriate ground state value.
! 
!   Nzeta is determined by the first two characters of the 'basis_size'
!   string: 1 for 'sz' and  2 for 'dz'.
! 
!   There are no 'per l-shell' polarization orbitals, except if the
!   third character of 'basis_size' is 'p' (as in 'dzp'), in which
!   case polarization orbitals are defined so they have the minimum      
!   angular momentum l such that there are not occupied orbitals 
!   with the same l in the valence shell of the ground-state 
!   atomic configuration. They polarize the corresponding l-1 shell.

! 
!   The Soft-Confinement parameters 'rinn' and 'vcte' are set to 0.0
! 
!   rc(1:nzeta) is set to 0.0
!   lambda(1:nzeta) is set to 1.0  (this is a change from old practice)
! 
!  ----------------------------------
!  
!   Next come the blocks associated to the KB projectors:
! 
!   Block PS.lmax (routine relmxkb) assigns a per-species maximum l for
!   the KB projectors. Internally, it sets the 'lmxb_requested'
!   component of the data structure, which by default is -1.
!  
!   Block PS.KBprojectors (optional) sets the number of KB projectors
!   per angular momentum. The same routine which reads this block
!   (readkb) also sets any appropriate defaults:
! 
! 
!   If the species does not appear in the  PS.KBprojectors block, then
!       Setting of lmxkb:
!       If lmxkb is set in PS.lmax, 
!          set it to that
!       else
!          Use PAO information: (routine set_default_lmxkb)
!          (Set it to a meaningless value for floating orbitals)
!           Set it to lmxo+1, or to lpol+1, where lpol is the angular
!           momentum of the highest-l polarization orbital.
!       endif
! 
!       The  number of KB projectors per l is set to the number of
!       n-shells in the corresponding PAO shell with the same l. For l
!       greater than lmxo, it is set to 1. The reference energies are
!       in all cases set to huge(1.d0) to mark the default.
!       
! =======================================================================
!
      use precision
      use basis_types, only: basis_def_t, shell_t, lshell_t, kbshell_t
      use basis_types, only: nsp, basis_parameters, ground_state_t
      use basis_types, only: destroy, copy_shell, initialize
      use pseudopotential
      use periodic_table
      use chemical
      use sys

      Use fdf
      use parse

      Implicit None

      interface
         function leqi(s1,s2)
         logical leqi
         character(len=*), intent(in)   :: s1, s2
         end function leqi
      end interface

      type(basis_def_t), pointer::   basp
      type(shell_t), pointer::  s
      type(lshell_t), pointer::  ls
      type(kbshell_t), pointer:: k

      type(block), pointer  :: bp
      type(parsed_line), pointer  :: p
      character(len=132) line

      character(len=*), parameter   :: defunit='Ry'
      character(len=1), parameter   ::
     $                           sym(0:4) = (/ 's','p','d','f','g' /)

      integer           isp         ! just an index dummy variable for the whole module

      public :: read_basis_specs

      private

      CONTAINS

!---
      subroutine read_basis_specs()

      character(len=15), parameter  :: basis_size_default='standard'
      character(len=10), parameter  :: basistype_default='split'

      character(len=15) :: basis_size
      character(len=10) :: basistype_generic


      basis_size=fdf_string('PAO.BasisSize',basis_size_default)
      call size_name(basis_size)
      basistype_generic=fdf_string('PAO.BasisType',basistype_default)
      call type_name(basistype_generic)

!
!     Use standard routine in chemical module to process the
!     chemical species
!
      call read_chemical_types
      nsp = number_of_species()

      allocate(basis_parameters(nsp))
      do isp=1,nsp
         call initialize(basis_parameters(isp))
      enddo

      do isp=1,nsp
         basp=>basis_parameters(isp)

         basp%label = species_label(isp)
         basp%z = atomic_number(isp)
         basp%floating = is_floating(isp)
         basp%bessel = is_bessel(isp)

         basp%basis_size = basis_size
         basp%basis_type = basistype_generic
         if (basp%floating) then
            basp%mass = 1.d40   ! big but not too big, as it is used
                                ! later in computations
         else
            basp%mass = atmass(abs(basp%z))
         endif
         if (.not. basp%bessel) then
            call ground_state(abs(basp%z),basp%ground_state)
            call pseudo_read(basp%label,basp%pseudopotential)
         endif
         call semicore_check(isp)
      enddo

      call remass
      call resizes
      call repaobasis
      call autobasis
      call relmxkb
      call readkb

!      do isp=1,nsp
!         call print_basis_def(basis_parameters(isp))
!      enddo

      end subroutine read_basis_specs
!-----------------------------------------------------------------------

      function label2species(str) result(index)
      integer index
      character(len=*), intent(in) ::  str

      integer i

      index = 0
      do i=1,nsp
         if(leqi(basis_parameters(i)%label,str)) index = i
      enddo
      end function label2species
!-----------------------------------------------------------------------

      subroutine ground_state(z,gs)
      integer, intent(in)               ::  z
      type(ground_state_t), intent(out) :: gs
!
!     Determines ground state valence configuration from Z
!
      integer l, latm

      gs%z_valence = 0.d0
      do l=0,3
         gs%occupation(l)=0.0d0
      enddo

      call lmxofz(z,gs%lmax_valence,latm)
      call qvlofz(z,gs%occupation(:))
      do l=0,gs%lmax_valence
         gs%z_valence = gs%z_valence + gs%occupation(l)
      enddo
      call cnfig(z,gs%n(0:3))

      write(6,'(a,i2)',advance='no')
     $     'Ground state valence configuration: '
      gs%occupied(4) = .false.         !! always
      do l=0,3
        gs%occupied(l) =  (gs%occupation(l).gt.0.1d0)
        if (gs%occupied(l))
     $   write(6,'(2x,i1,a1,i2.2)',advance='no')
     $       gs%n(l),sym(l),nint(gs%occupation(l))
      enddo
      write(6,'(a)') ''

      end subroutine ground_state

!---------------------------------------------------------------
      subroutine readkb
      integer lpol, isp, ish, i, l
      character(len=20) unitstr

      lpol = 0

      nullify(bp)
      if (.not. fdf_block('PS.KBprojectors',bp) ) goto 2000
 
! First pass to find out about lmxkb and set any defaults.

      loop: DO     !! over species
         if (.not. fdf_bline(bp,line)) exit loop
         p => digest(line)
         if (.not. match(p,"ni"))
     $        call die("Wrong format in PS.KBprojectors")
         isp = label2species(names(p,1))
         if (isp .eq. 0) then
            write(6,'(a,1x,a)')
     $        "WRONG species symbol in PS.KBprojectors:",
     $                          trim(names(p,1))
            call die
         endif
         basp=>basis_parameters(isp)
         basp%nkbshells = integers(p,1)
         call destroy(p)
         do ish=1,basp%nkbshells
            if (.not. fdf_bline(bp,line)) call die("No l nkbl")
            p => digest(line)
            if (.not. match(p,"ii")) call die("Wrong format l nkbl")
            l = integers(p,1)
            if (l .gt. basp%lmxkb) basp%lmxkb = l
            call destroy(p)
            if (.not. fdf_bline(bp,line)) then
                   if (ish .ne. basp%nkbshells)
     $              call die("Not enough shells for this species...")
                   ! There is no line with ref energies
            else 
               p => digest(line)
               if (match(p,"ni")) then
                  ! We are seeing the next species' section
                  if (ish .ne. basp%nkbshells)
     $              call die("Not enough shells for this species...")
                  call backspace(bp)
                  call destroy(p)
               else
                  ! would set erefs
                  call destroy(p)
               endif
            endif
         enddo       ! end of loop over shells for species isp

      enddo loop
      call destroy(bp)

 2000 CONTINUE

      do isp=1, nsp
!
!      Fix defaults and allocate kbshells
!
         basp=>basis_parameters(isp)
         if (basp%lmxkb .eq. -1) then ! not set in KBprojectors 
            if(basp%lmxkb_requested.eq.-1) then ! not set in PS.lmax
               basp%lmxkb = set_default_lmxkb(isp) ! Use PAO info
            else
               basp%lmxkb = basp%lmxkb_requested
            endif
            allocate(basp%kbshell(0:basp%lmxkb))
            do l=0,basp%lmxkb
               call initialize(basp%kbshell(l))
               k=>basp%kbshell(l)
               k%l = l
               if (l.gt.basp%lmxo) then
                  k%nkbl = 1
               else           ! Set equal to the number of PAO shells
                  k%nkbl = basp%lshell(l)%nn
                  if (k%nkbl.eq.0) then
                     write(6,*) 'Warning: Empty PAO shell. l =', l
                     write(6,*) 'Will have a KB projector anyway...'
                     k%nkbl = 1
                  endif
               endif
               allocate(k%erefkb(1:k%nkbl))
               k%erefkb(1:k%nkbl) = huge(1.d0)
            enddo
         else          ! Set in KBprojectors
            if(basp%lmxkb_requested.ne.-1) then ! set in PS.lmax
               if (basp%lmxkb.ne.basp%lmxkb_requested) then
                  call die("LmaxKB conflict")
               endif
            endif
            !! OK, we have a genuine lmxkb
            allocate(basp%kbshell(0:basp%lmxkb))
            do l=0,basp%lmxkb
               call initialize(basp%kbshell(l))
            enddo
         endif
         if (basp%z .le. 0) then
            if (basp%lmxkb .ne. -1)
     $        call die("Floating orbs cannot have KB projectors...")
         endif
      enddo

!
!     Now re-scan the block (if it exists) and fill in as instructed
!            
      if (.not. fdf_block('PS.KBprojectors',bp) ) return

      loopkb: DO     !! over species
         if (.not. fdf_bline(bp,line)) exit loopkb
         p => digest(line)
         if (.not. match(p,"ni"))
     $        call die("Wrong format in PS.KBprojectors")
         isp = label2species(names(p,1))
         if (isp .eq. 0) then
            write(6,'(a,1x,a)')
     $        "WRONG species symbol in PS.KBprojectors:",
     $                          trim(names(p,1))
            call die
         endif
         basp=>basis_parameters(isp)
         basp%nkbshells = integers(p,1)
         call destroy(p)
         do ish=1,basp%nkbshells
            if (.not. fdf_bline(bp,line)) call die("No l nkbl")
            p => digest(line)
            if (.not. match(p,"ii")) call die("Wrong format l nkbl")
            l = integers(p,1)
            k=>basp%kbshell(l)
            k%l = l
            k%nkbl = integers(p,2)
            call destroy(p)
            allocate(k%erefkb(k%nkbl))
            if (.not. fdf_bline(bp,line)) then
                   if (ish .ne. basp%nkbshells)
     $              call die("Not enough shells for this species...")
                   ! There is no line with ref energies
                   ! Use default values
                   k%erefKB(1:k%nkbl)=huge(1.d0)
            else 
               p => digest(line)
               if (match(p,"ni")) then
                  ! We are seeing the next species' section
                  if (ish .ne. basp%nkbshells)
     $              call die("Not enough shells for this species...")
                   ! Use default values for ref energies
                  k%erefKB(1:k%nkbl)=huge(1.d0)
                  call backspace(bp)
                  call destroy(p)
               else
                  if (nvalues(p) .ne. k%nkbl)
     $                 call die("Wrong number of energies")
                  unitstr = defunit
                  if (nnames(p) .eq. 1) unitstr = names(p,1)
                  ! Insert ref energies in erefkb
                  do i=1,k%nkbl
                     k%erefKB(i) =
     $                    values(p,i)*fdf_convfac(unitstr,defunit)
                  enddo
                  call destroy(p)
               endif
            endif
         enddo            ! end of loop over shells for species isp

         ! For those l's not specified in block, use default values
         do l=0, basp%lmxkb
            k=>basp%kbshell(l)
            if (k%l.eq.-1) then
               k%l = l
               k%nkbl = 1
               allocate(k%erefkb(1))
               k%erefkb(1) = huge(1.d0)
            endif
         enddo
      enddo loopkb   !! Over species

      end subroutine readkb
!---------------------------------------------------------------

      subroutine repaobasis

      integer isp, ish, nn, i, ind, l, indexp

      nullify(bp)
      If (.not.fdf_block('PAO.Basis',bp)) RETURN

      loop: DO     !! over species
         if (.not. fdf_bline(bp,line)) exit loop
         p => digest(line)
         if (.not. match(p,"ni"))
     $        call die("Wrong format in PAO.Basis")
         isp = label2species(names(p,1))
         if (isp .eq. 0) then
            write(6,'(a,1x,a)')
     $        "WRONG species symbol in PAO.Basis:",
     $                          trim(names(p,1))
            call die
         endif

         basp=>basis_parameters(isp)
         basp%label=names(p,1)
         basp%nshells_tmp = integers(p,1)
         basp%lmxo = 0
         !! Check whether there are optional type and ionic charge
         if (nnames(p).eq.2) basp%basis_type=names(p,2)
         if (nvalues(p).eq.2) basp%ionic_charge=values(p,2)
         call destroy(p)
         allocate(basp%tmp_shell(basp%nshells_tmp))

         shells: do ish=1,basp%nshells_tmp
            s=>basp%tmp_shell(ish)
            call initialize(s)
            if (.not. fdf_bline(bp,line)) call die("No l nzeta, etc")

            p => digest(line)
            if (match(p,"niii")) then
              s%n = integers(p,1)
              s%l = integers(p,2)
              basp%lmxo = max(basp%lmxo,s%l)
              s%nzeta = integers(p,3)
            else if (match(p,"ii")) then
              !    l, nzeta

              if (basp%semic)
     $        call die("Please specify n if there are semicore states")

              s%l = integers(p,1)
              s%n = basp%ground_state%n(s%l)
              s%nzeta = integers(p,2)
              basp%lmxo = max(basp%lmxo,s%l)
            else
               call die("Bad format of (n), l, nzeta line in PAO.Basis")
            endif
            ! Optional stuff: Polarization and Soft-confinement Potential
            if (search(p,"P",indexp)) then
               s%polarized = .true.
               if (match(p,"i",after=indexp)) then
                  s%nzeta_pol=integers(p,ind=1,after=indexp)
               else
                  s%nzeta_pol = 1
               endif
            endif
            if (search(p,"E",indexp)) then
               if (match(p,"vv",after=indexp)) then
                  s%vcte = values(p,ind=1,after=indexp)
                  s%rinn = values(p,ind=2,after=indexp)
               else
                  call die("Need vcte and rinn after E in PAO.Basis")
               endif
            endif
            call destroy(p)

            allocate(s%rc(s%nzeta),s%lambda(s%nzeta))
            s%rc(:) = 0.d0
            s%lambda(:) = 1.d0
            if (.not. fdf_bline(bp,line)) call die("No rc's")
            p => digest(line)
            if (nvalues(p).ne.s%nzeta) call die("Wrong number of rc's")
            do i=1,s%nzeta
               s%rc(i) = values(p,i)
            enddo
            call destroy(p)
            ! Optional scale factors. They MUST be reals, or else...
            if (.not. fdf_bline(bp,line)) then
               if (ish.ne.basp%nshells_tmp)
     $                    call die("Not enough shells")
               ! Default values for scale factors
            else
               p => digest(line)
               if (.not.match(p,"r")) then
                  ! New shell or species
                  ! Default values for the scale factors
                  call backspace(bp)
                  cycle shells
               else
                  if (nreals(p).ne.s%nzeta)
     $                 call die("Wrong number of lambda's")
                  do i=1,s%nzeta
                     s%lambda(i) = reals(p,i)
                  enddo
               endif
               call destroy(p)
            endif
         enddo shells
         ! Clean up for this species
      enddo loop

!
!        OK, now classify the states by l-shells
!
         do isp = 1, nsp
            basp=>basis_parameters(isp)
            if (basp%lmxo.eq.-1) cycle !! Species not in block
                                       !! 
            allocate (basp%lshell(0:basp%lmxo))
            loop_l: do l=0,basp%lmxo
               ls=>basp%lshell(l)
               call initialize(ls)
               ls%l = l
               ! Search for tmp_shells with given l
               nn = 0
               do ish=1, basp%nshells_tmp
                  s=>basp%tmp_shell(ish)
                  if (s%l .eq. l) nn=nn+1
               enddo
               ls%nn = nn
               if (nn.eq.0) then
                  !! What else do we do here?
                  cycle loop_l  
               endif
                                          !! 
               allocate(ls%shell(1:nn))
               ! Collect previously allocated shells
               ind = 0
               do ish=1, basp%nshells_tmp
                  s=>basp%tmp_shell(ish)
                  if (s%l .eq. l) then
                     ind = ind+1
                     call copy_shell(source=s,target=ls%shell(ind))
                  endif
               enddo
               if (nn.eq.1) then
                  ! If n was not specified, set it to ground state n
                  if (ls%shell(1)%n.eq.-1)
     $                 ls%shell(1)%n=basp%ground_state%n(l)
               endif
               !! Do we have to sort by n value????
               !!
            enddo loop_l
            !! Destroy temporary shells in basp
            !! Warning: This does seem to destroy information!!
            call destroy(basp%tmp_shell)
         enddo
         call destroy(bp)

         end subroutine repaobasis
!_______________________________________________________________________


         function set_default_lmxkb(is) result(lmxkb)
         integer lmxkb
         integer, intent(in) :: is

         integer lpol, l, i

         lmxkb = -1
         basp=>basis_parameters(is)
         if (basp%z .le. 0) return     ! Leave it at -1 for floating orbs.

         lmxkb = basp%lmxo + 1

         ! But watch out for polarization orbitals...
         ! We ASSUME that these do not count towards lshell%nn...

         lpol = 0
         do l = 0, basp%lmxo
            ls=>basp%lshell(l)
            do i = 1, ls%nn
               s=>ls%shell(i)
               if (s%polarized) lpol = s%l + 1
            enddo
         enddo
         if (lpol .gt. basp%lmxo) lmxkb = lpol + 1

         write(6,'(3a,i1,/,2a,/,a)') 'For ', trim(basp%label),
     $              ', standard SIESTA heuristics set lmxkb to ',
     $              lmxkb,
     $              ' (one more than the basis l,',
     $              ' including polarization orbitals).',
     $  'Use PS.lmax or PS.KBprojectors blocks to override.'
!
!        But there is an upper limit for sanity: f is the highest
!
         if (lmxkb.gt.3) then
            write(6,'(3a,i1)') 'Warning: For ', trim(basp%label),
     $           ' lmxkb would have been set to ', lmxkb
            write(6,'(a)')
     $           'Setting it to maximum value of 3 (f projector)'
            lmxkb = 3
         endif

         end function set_default_lmxkb

!-----------------------------------------------------------------------

      subroutine resizes

c Reading atomic basis sizes for different species.
c
c Reads fdf block. Not necessarily all species have to be given. The
c ones not given at input will be assumed to have the basis sizes
c given by the general input PAO.BasisSize, or its default value.

      type(block), pointer  :: bp
      type(parsed_line), pointer  :: p
      character(len=132) line

      integer isp

      nullify(bp)
      if (.not. fdf_block('PAO.BasisSizes',bp) ) return
      loop: DO
         if (.not.fdf_bline(bp,line)) exit loop
         p => digest(line)
         if (.not. match(p,"nn"))
     $        call die("Wrong format in PAO.BasisSizes")
         isp = label2species(names(p,1))
         if (isp .eq. 0) then
            write(6,'(a,1x,a)')
     $           "WRONG species symbol in PAO.BasisSizes:",
     $            trim(names(p,1))
            call die
         else
            basp=>basis_parameters(isp)
            basp%basis_size = names(p,2)
            call size_name(basp%basis_size)   !!! DEPRECATED
            write(6,'(4a)')
     .           'resizes: Read basis size for species ',
     .             trim(names(p,1)),' = ',basp%basis_size
         endif
         call destroy(p)
      enddo loop
      call destroy(bp)

      end subroutine resizes

!-----------------------------------------------------------------------

      subroutine relmxkb

c Reads the maximum angular momentum of the Kleinman-Bylander
c projectors for the different species.
c
c Reads fdf block. Not necessarily all species have to be given.

      type(block), pointer  :: bp
      type(parsed_line), pointer  :: p
      character(len=132) line

      integer isp

      nullify(bp)
      if (.not. fdf_block('PS.lmax',bp) ) return
      loop: DO
         if (.not.fdf_bline(bp,line)) exit loop
         p => digest(line)
         if (.not. match(p,"ni")) call die("Wrong format in PS.lmax")
         isp = label2species(names(p,1))
         if (isp .eq. 0) then
            write(6,'(a,1x,a)') "WRONG species symbol in PS.lmax:",
     $                          trim(names(p,1))
            call die
         else
            basp=>basis_parameters(isp)
            basp%lmxkb_requested = integers(p,1)
            write(6,"(a, i4, 2a)")
     .            'relmxkb: Read Max KB Ang. Momentum= ',
     $             basp%lmxkb_requested,
     .            ' for species ', trim(names(p,1))
         endif
         call destroy(p)
      enddo loop
      call destroy(bp)

      end subroutine relmxkb
!-----------------------------------------------------------------------

      subroutine remass


      type(block), pointer  :: bp
      type(parsed_line), pointer  :: p
      character(len=132) line

      integer isp

c Read atomic masses of different species.
c
c Reads fdf block. Not necessarily all species have to be given. The
c ones not given at input will be assumed to have their natural mass
c (according to atmass subroutine).

      nullify(bp)
      if (.not. fdf_block('AtomicMass',bp) ) return
      loop: DO
        if (.not. fdf_bline(bp,line)) exit loop
        p => digest(line)
        if (.not. match(p,"iv"))
     $       call die("Wrong format in AtomicMass")
        isp = integers(p,1)
        if (isp .gt. nsp .or. isp .lt. 1)
     $       call die("Wrong specnum in AtomicMass")
        basp=>basis_parameters(isp)
        basp%mass = values(p,2)
        write(6,"(a, i4, a, f12.5)")
     .       'remass: Read atomic mass for species ', isp,
     .       ' as ', basp%mass
        call destroy(p)
      enddo loop
      call destroy(bp)

      end subroutine remass
!-----------------------------------------------------------------------

      subroutine size_name(str)
      character(len=*), intent(inout)  ::  str


      if(leqi(str,'STANDARD')) str='dzp'
      if(leqi(str,'DZP'))  str='dzp'
      if(leqi(str,'DZ')) str='dz'
      if(leqi(str,'MINIMAL')) str='sz'
      if(leqi(str,'SZ'))  str='sz'
      if(leqi(str,'SZP')) str='szp'
 
      if( (str.ne.'szp').and.(str.ne.'sz').and.
     .    (str.ne.'dz') .and.(str.ne.'dzp') ) then

         write(6,'(/,2a,/,a)')
     .   'size_name: Incorrect basis-size option specified,',
     .   ' active options are:',
     .   'SZ, SZP, DZ, DZP, and STANDARD'

         call die
      endif

      end subroutine size_name
!-----------------------------------------------------------------------

        subroutine type_name(basistype)

         character basistype*(*)

         if(leqi(basistype,'NODES')) then
               basistype='nodes'
         elseif(leqi(basistype,'NONODES')) then
               basistype='nonodes'
         elseif(leqi(basistype,'SPLIT')) then
               basistype='split'
         elseif(leqi(basistype,'SPLITGAUSS')) then
               basistype='splitgauss'
         else
              write(6,'(/,2a,(/,5(3x,a)),(/,2(3x,a)))')
     .        'type_name: Incorrect basis-type option specified,',
     .        ' active options are:',
     .        'NODES','SPLIT','SPLITGAUSS','NONODES'
              call die
         endif

        end subroutine type_name
!-----------------------------------------------------------------------
!
!       Find out whether semicore states are implied by the valence
!       charge density in the pseudopotential file.
!
        subroutine semicore_check(is)
        integer, intent(in)  :: is

        real*8, parameter :: tiny = 1.d-5
        integer ndiff
        real*8 zval, zval_vps, charge_loc

        basp => basis_parameters(is)

        basp%semic = .false.
        if (basp%bessel) return

        zval_vps = basp%pseudopotential%zval
        zval = basp%ground_state%z_valence

        if(abs(Zval-zval_vps).lt.tiny) return

        ndiff=nint(abs(Zval-zval_vps))
        if(abs(ndiff-abs(Zval-zval_vps)).gt.tiny) then
              write(6,'(2a)')
     .   'ERROR expected valence charge for species ',
     .        basp%label
              write(6,'(a)')
     .  'ERROR and the value read from the vps file'
              write(6,'(a,f6.3,a,f6.3)')
     .  'ERROR differ:  Zval(expected)= ', Zval,' Zval(vps)= ',zval_vps
              call die
        endif

        basp%semic=.true.
        charge_loc=Zval_vps-Zval
        write(6,'(a,i2,a)')
     .       'Semicore shell(s) with ', nint(charge_loc),
     $       ' electrons included in the valence for', trim(basp%label)

        end subroutine semicore_check
!----------------------------------------------------------------------
      subroutine autobasis

!
!     It sets the defaults if a species has not been included
!     in the PAO.Basis block
!
      integer l, nzeta

      loop: do isp=1, nsp
         basp=>basis_parameters(isp)
         if (basp%lmxo .ne. -1) cycle loop   ! Species already set
                                             ! in PAO.Basis block
         if (basp%semic) then
            write(6,'(2a)') basp%label,
     $           ' must be in PAO.Basis (it has semicore states)'
            call die
         endif
         if (basp%bessel) then
            write(6,'(2a)') basp%label,
     $      ' must be in PAO.Basis (it is a floating Bessel function)'
            call die
         endif
         !
         ! Set the default max l 
         !
         basp%lmxo = basp%ground_state%lmax_valence

         allocate (basp%lshell(0:basp%lmxo))

         if (basp%basis_size(1:2) .eq. 'sz') nzeta = 1
         if (basp%basis_size(1:2) .eq. 'dz') nzeta = 2

         loop_l: do l=0, basp%lmxo
            ls=>basp%lshell(l)
            call initialize(ls)
            ls%l = l
            ls%nn = 1
            allocate(ls%shell(1:1))
            s => ls%shell(1)
            s%l = l
            s%n = basp%ground_state%n(l)
            if (basp%ground_state%occupied(l)) then
               s%nzeta = nzeta
            else
               s%nzeta = 0
            endif
            s%polarized = .false.
            s%nzeta_pol = 0
            s%rinn = 0.d0
            s%vcte = 0.d0
            if (s%nzeta .ne.0) then
               allocate(s%rc(1:s%nzeta))
               allocate(s%lambda(1:s%nzeta))
               s%rc(1:s%nzeta) = 0.d0
               s%lambda(1:s%nzeta) = 1.d0
            endif
         enddo loop_l

         if (basp%basis_size(3:3) .eq. 'p') then

         ! Polarization orbitals are defined so they have the minimum      
         ! angular momentum l such that there are not occupied orbitals 
         ! with the same l in the valence shell of the ground-state 
         ! atomic configuration. They polarize the corresponding l-1    
         ! shell.

         ! note that we go up to l=4 to make the loop simpler.
         ! (that is the reason why 'occupied' is dimensioned to 0:4)

            loop_angmom: do l=1,4
            if (.not. basp%ground_state%occupied(l)) then
                  ls=>basp%lshell(l-1)
                  ls%shell(1)%polarized = .true.
                  ls%shell(1)%nzeta_pol = 1
                  exit loop_angmom
               endif
            enddo loop_angmom

         endif

      enddo loop
            
      end subroutine autobasis

      End module basis_specs










