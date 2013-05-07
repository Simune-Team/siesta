      MODULE m_new_dm

!     Prepares a starting density matrix for a new geometry iteration
!     This DM can be:
!     1. Synthesized directly from atomic occupations (not idempotent)
!     2. Read from file
!     3. Re-used (with possible extrapolation) from previous geometry step(s).
!
!     In cases 2 and 3, the structure of the read or extrapolated DM 
!     is automatically adapted to the current sparsity pattern.
!
!     Special cases:
!            Harris: The matrix is always initialized
!            Force calculation: The DM should be written to disk
!                               at the time of the "no displacement"
!                               calculation and read from file at
!                               every subsequent step.
!            Variable-cell calculation:
!              If the auxiliary cell changes, the DM is forced to be
!              initialized (conceivably one could rescue some important
!              information from an old DM, but it is too much trouble
!              for now). NOTE that this is a change in policy with respect
!              to previous versions of the program, in which a (blind?)
!              re-use was allowed, except if 'ReInitialiseDM' was 'true'.
!              Now 'ReInitialiseDM' is 'true' by default. Setting it to
!              'false' is not recommended. (This fdf variable maps to the
!              'initdmaux' module variable)
!
!              In all other cases (including "server operation"), the
!              default is to allow DM re-use (with possible extrapolation)
!              from previous geometry steps.
!              The fdf variables 'DM.AllowReuse' and 'DM.AllowExtrapolation'
!              (mapped to 'allow_dm_reuse' and 'allow_dm_extrapolation', and
!              both 'true' by default) can be used to change this behavior.
!
!              There is no re-use of the DM for "Forces", and "Phonon"
!              dynamics types (i.e., the DM is re-initialized)
!
!              For "CG" calculations, the default is not to extrapolate the
!              DM (unless requested by setting 'DM.AllowExtrapolation' to
!              "true"). The previous step's DM is reused.
!
!     Alberto Garcia, September 2007, April 2012
!

      use sys, only: die
      use precision, only: dp
      use alloc, only: re_alloc, de_alloc
      use parallel,  only: IOnode
      use fdf,       only: fdf_get

      implicit none

      character(len=*),parameter:: modName = 'm_new_dm '

      private
      public :: new_dm
      public :: get_allowed_history_depth

      CONTAINS

!=====================================================================

      subroutine  new_dm( auxchanged, DM_history, DMnew)

      USE siesta_options
      use siesta_geom,      only: xa, na_u
      use sparse_matrices,  only: sparse_pattern, block_dist
      use atomlist,         only: datm, iaorb, lasto, no_u, no_l
      use m_steps,          only: istp
      use m_spin,   only: nspin

      use class_dSpData2D
      use class_Sparsity
      use class_Pair_Geometry_dSpData2D
      use class_Fstack_Pair_Geometry_dSpData2D

      implicit none

      logical, intent(in) :: auxchanged ! Has auxiliary supercell changed?
      type(Fstack_Pair_Geometry_dSpData2D), intent(inout)      :: DM_history
      type(dSpData2D), intent(inout)   :: DMnew

!     Local variables

      logical :: dminit     ! Initialize density matrix?
      logical :: try_to_read_from_file
      integer :: n_dms_in_history, n_depth

      if (IOnode) then
         write(6,"(a,i5)") "New_DM. Step: ", istp
      endif


!     In principle we allow the re-use of the DM (i.e, we do not initialize it)
!
      dminit = .false.
      try_to_read_from_file = usesavedm     ! As per defaults
!
!     Except if there are explicit instructions
!
      if (.not. allow_dm_reuse) then
         dminit = .true.
         try_to_read_from_file = .false.  ! In case the user has a fossil DM.UseSaveDM
         if (IOnode) then
            write(6,"(a)") "DM re-use not allowed. Resetting always"
            if (usesavedm) then
               write(6,"(a)") "DM.UseSaveDM  overriden !!"
            endif
         endif
      endif
!
!     or using Harris...
!
      if (harrisfun) dminit = .true.
!
!     or we are in the first step, or performing force-constant calculations

      n_dms_in_history = n_items(DM_history)

      if (n_dms_in_history== 0) then
         dminit = .true.
      else
         if ((idyn .eq. 6)             &   ! Force Constants
              .and. usesavedm .and. writedm)  dminit = .true.
         if ((idyn .eq. 7)             &   ! Phonon series (writedm??)
               .and. usesavedm)  dminit = .true.
      endif

!
!     ... or if the auxiliary cell has changed
!     (in this case we have to  avoid reading back saved copy from file)
!
      if (auxchanged) then
         if (initdmaux) then
            dminit = .true.
            try_to_read_from_file = .false.
            if (IOnode) then
               write(6,"(a)") "DM history reset as supercell changed."
            endif
	    call get_allowed_history_depth(n_depth)
            call new(DM_history,n_depth,"(reset DM history stack)")
         else
            if (IOnode) then
               write(6,"(a)") "** Warning: DM history NOT reset upon supercell change"
               write(6,"(a)") "** Warning: since 'ReinitialiseDM' is set to .false."
            endif
         endif
      endif


      if (dminit) then
         if (IOnode) then
            write(6,"(a)") "Initializing Density Matrix..."
         endif

         call initdm(Datm, DMnew, sparse_pattern, block_dist,   &
                     lasto, no_u,                       &
                     no_l, nspin, na_u,           &
                     iaorb, inspn,                             &
                     try_to_read_from_file)

      else    ! not initializing the DM

         if (IOnode) then
            write(6,"(a)") "Re-using DM from previous geometries..."
         endif

        ! Extrapolation or simple re-structuring

         if (ionode) print "(a,i0)", "N DMs in history: ", n_dms_in_history
         ! if (ionode) call print_type(DM_history)
         call extrapolate_dm_with_coords(DM_history,na_u,xa(:,1:na_u),sparse_pattern,DMnew)
         if (ionode)  print "(a)", "New DM after history re-use:"
         if (ionode)  call print_type(DMnew)

      endif

      END subroutine new_dm

!====================================================================

      subroutine initdm(Datm, DMnew, sparse_pattern, block_dist, &
                        lasto, no_u,                              &
                        no_l, nspin, na_u,     &
                        iaorb, inspn,                              &
                        try_dm_from_file)

! Density matrix initialization
!
!    If Try_Dm_From_File is true, it is read from file if present.
!    Otherwise it is generated assuming atomic charging
!      (filling up atomic orbitals).


! logical try_dm_from_file     : whether DM has to be read from files or not
! logical found         : whether DM was found in files
! logical inspn         : true : AF ordering according to atom ordering
!                                if no DM files, no DM.InitSpin, ispin=2
!                         false: Ferro ordering  (fdf DM.InitSpinAF)
! integer na_u           : Number of atoms in the unit cell
! integer no_l           : Number of orbitals in the unit cell (local)
! integer no_u           : Number of orbitals in the unit cell (global)
! integer nspin         : Number of spin components
! integer lasto(0:na_u) : List with last orbital of each atom
! integer iaorb(no_u)   : List saying to what atom an orbital belongs
! double Datm(no_u)       : Occupations of basis orbitals in free atom

      use precision, only : dp
      use files,     only : slabel
      use class_Sparsity
      use class_dSpData2D
      use class_OrbitalDistribution
      use class_dData2D
      use m_readSpData2D, only: readdSpData2D
#ifdef TRANSIESTA
      use sparse_matrices, only : EDM, Escf
      use m_ts_iodm, only       : ts_init_dm
      use m_energies, only: ef  ! Transiesta uses the EF obtained in a initial SIESTA run
                                ! to place the electrodes and scattering region energy
                                ! levels at the appropriate relative position, so it is
                                ! stored in the TSDE file.
      use m_ts_options,   only : TSmode, ImmediateTSmode
#endif /* TRANSIESTA */

      implicit          none

      logical           dm_found, inspn, try_dm_from_file
      integer           no_l, na_u, no_u, nspin
      integer           lasto(0:na_u), iaorb(no_u)
      real(dp)          Datm(no_u)

      type(dSpData2D), intent(inout)      :: DMnew
      type(Sparsity), intent(in) :: sparse_pattern
      type(OrbitalDistribution), intent(in) :: block_dist

! ---------------------------------------------------------------------

      character(len=*),parameter:: myName = 'initdm'

      integer :: nspin_read
      real(dp), pointer              :: Dscf(:,:)
      integer, pointer, dimension(:) :: numh, listhptr, listh
      type(dSpData2D)                 :: DMread
      type(dData2D)                 :: dm_a2d
#ifdef TRANSIESTA
      logical                        :: tsde_found
      type(dSpData2D)                 :: EDMread
#endif

! Try to read DM from disk if wanted (DM.UseSaveDM true) ---------------

#ifdef TRANSIESTA
      dm_found = .false.
      tsde_found = .false.
      if (try_dm_from_file) then
         if (TSmode) then
            if (ionode) print *, "Attempting to read DM,EDM from TSDE file..."
            call readdSpData2D(trim(slabel)//".TSDE",   &
                 DMread,tsde_found,block_dist,EDMread,ef)
            !
            call ts_init_dm(tsde_found)
            !
            if (.not. tsde_found) then
               if (ionode) print *, "Attempting to read DM from file (TSmode)"
               call readdSpData2D(trim(slabel)//".DM",   &
                    DMread,dm_found,block_dist)
            endif
            dm_found = (tsde_found .or. dm_found)

            ! if the user requests to start the transiesta SCF immediately.
            ! We will allow this if a DM file is found.
            !
            if ( dm_found .and. ImmediateTSmode .and. &
                 .not. tsde_found .and. .FALSE. ) then
               ! We need a way to ensure a correct Escf for the transiesta
               ! i.e. we should read in the E scf in some way as well
               ! Otherwise the energies are incorrect.
               ! So for now this will not be used
               ! We will not reset tsde_found as that will
               ! mess up the Escf array (see further down)
               call ts_init_dm(.true.)
            end if

         else  ! Not TSmode

            if (ionode) print *, "Attempting to read DM from file..."
            call readdSpData2D(trim(slabel)//".DM",   &
                 DMread,dm_found,block_dist)

         endif
      endif
#else
      dm_found = .false.
      if (try_dm_from_file) then
         if (ionode) print *, "Attempting to read DM from file..."
         call readdSpData2D(trim(slabel)//".DM",   &
                           DMread,dm_found,block_dist)
      endif
#endif

! If DM found, check and update, otherwise initialize with neutral atoms

      if (dm_found) then
        ! Various degrees of sanity checks

        nspin_read = size(val(DMread),dim=2)
        if (nspin_read /= Nspin) then
           if (IOnode) then
              write(6,"(a,i6,/,a)")                   &
              "WARNING: Wrong nspin in DM file: ",  nspin_read,  &
              "WARNING: Falling back to atomic initialization of DM."
           endif
           dm_found = .false.
        endif

        if (nrows_g(DMread) /= nrows_g(sparse_pattern)) then
           if (IONode) then
              write(6,"(a,/,a)")                             &
             "WARNING: Wrong number of orbs in DM file. ",     &
             "WARNING: Falling back to atomic initialization of DM."
           endif
           dm_found = .false.
        endif

      endif
      
      ! Density matrix
      if (dm_found) then

	call restructdSpData2D(DMread,sparse_pattern,DMnew)
        if (ionode) print *, "DMread after reading file:"
        if (ionode) call print_type(Dmread)

      else

	call newdData2D(dm_a2d,nnzs(sparse_pattern),nspin,"(DMatomic)")
	Dscf => val(dm_a2d)
        numh     => n_col(sparse_pattern)
        listhptr => list_ptr(sparse_pattern)
        listh    => list_col(sparse_pattern)
	
        call   fill_dscf_from_atom_info(Datm, Dscf,              &
                        numh, listhptr, listh, lasto,         &
                        no_u, na_u, no_l, nspin,     &
                        iaorb, inspn)

        call newdSpData2D(sparse_pattern,dm_a2d,block_dist,DMnew,  &
                         "(DM initialized from atoms)")
        call delete(dm_a2d)
        if (ionode) print *, "DMnew after filling with atomic data:"
        if (ionode) call print_type(DMnew)

       endif

       ! Energy-density matrix
#ifdef TRANSIESTA
       if (dm_found) then
          if (tsde_found) then
             call restructdSpData2D(EDMread,sparse_pattern,EDM)
             if (ionode) print *, "EDMread after reading file:"
             if (ionode) call print_type(EDMread)
             Escf => val(EDM)
          else
             ! Escf remains associated to old EDM
          endif
       else
          ! Escf remains associated to old EDM
       endif
#else
       ! Escf remains associated to old EDM
#endif           

       ! Put deletes here to avoid complicating the logic
       call delete(DMread)
#ifdef TRANSIESTA
       call delete(EDMread)  
#endif

      end subroutine initdm

!======================================================================
      subroutine fill_dscf_from_atom_info(Datm, Dscf,              &
                        numh, listhptr, listh, lasto,        &
                        no_u,  na_u, no_l, nspin,     &
                        iaorb, inspn)


!      The DM is generated assuming atomic charging
!      (filling up atomic orbitals). The DM originated that way is
!      not a good DM due to overlaps, but the SCF cycling corrects
!      that for the next cycle.

!    Spin polarized calculations starting from atoms:
!      Default: All atoms with maximum polarization compatible with
!               atomic configuration. In Ferromagnetic ordering (up).
!      If DM.InitSpinAF is true, as default but in Antiferro order:
!               even atoms have spin down, odd up.
!      If fdf %block DM.InitSpin is present it overwrites previous
!         schemes: magnetic moments are explicitly given for some atoms.
!         Atoms not mentioned in the block are initialized non polarized.

! Written by E. Artacho. December 1997. Taken from the original piece
! of siesta.f written by P. Ordejon.
! Non-collinear spin added by J.M.Soler, May 1998.
! ********* INPUT ***************************************************
! logical try_dm_from_file     : whether DM has to be read from files or not
! logical found         : whether DM was found in files
! logical inspn         : true : AF ordering according to atom ordering
!                                if no DM files, no DM.InitSpin, ispin=2
!                         false: Ferro ordering  (fdf DM.InitSpinAF)
! integer na_u           : Number of atoms in the unit cell
! integer no_l           : Number of orbitals in the unit cell
! integer nspin         : Number of spin components
! integer no_u          : Max. number of orbitals (globally)
! integer lasto(0:maxa) : List with last orbital of each atom
! integer iaorb(no_u)   : List saying to what atom an orbital belongs
! double Datm(no)       : Occupations of basis orbitals in free atom
! ********* OUTPUT **************************************************
! double Dscf(:,:) : Density matrix in sparse form
! *******************************************************************

!
!  Modules
!
      use precision
      use parallel,     only : Node, Nodes, IOnode
      use parallelsubs, only : LocalToGlobalOrb, GlobalToLocalOrb
      use fdf
      use parsing
      use sys,          only : die
      use alloc,        only : re_alloc, de_alloc

#ifdef MPI
      use mpi_siesta
#endif
      use units, only : pi

      implicit          none

      logical, intent(in)       :: inspn
      integer, intent(in)       :: no_l, na_u, nspin, no_u
      integer, intent(in)       :: lasto(0:na_u), iaorb(no_u)
      real(dp), intent(in)      :: Datm(no_u)
      integer, intent(in), dimension(:) :: numh, listhptr, listh

      real(dp), intent(out)     :: Dscf(:,:)

! ---------------------------------------------------------------------

! Internal variables and arrays
      character(len=*),parameter:: myName = 'fill_dscf_from_atom_info'
      character         updo*1, msg*80
      logical           noncol, peratm, badsyntax
      integer           nh, ni, nn, nr, nv, iat, nat, ia, iu,   &
                        i1, i2, in, ind, ispin, jo, io,         &
                        iio, maxatnew

      integer, save ::  maxat

      integer           integs(4), lastc, lc(0:3)

      integer, pointer, save ::  atom(:)

#ifdef MPI
      integer  MPIerror
      logical  lbuffer
#endif
      real(dp)          aspin, cosph, costh, epsilon,        &
                        qio, rate, reals(4),                 &
                        sinph, sinth, spinat, spio, values(4)

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

      real(dp), pointer, save :: phi(:), spin(:), theta(:)

      integer :: is

      data maxat / 1000 /
      data epsilon / 1.d-8 /

! See whether specific initial spins are given in a DM.InitSpin block
! and read them in a loop on atoms where lines are read and parsed
!   integer nat       : how many atoms to polarize
!   integer atom(nat) : which atoms
!   double  spin(nat) : what polarization -----------------------------

        noncol = .false.
        peratm = fdf_block('DM.InitSpin',bfdf)
        noncol = .false.
        if (Node.eq.0) then
          if (peratm .and. nspin.lt.2) write(6,'(/,a)')             &
          'initdm: WARNING: DM.InitSpin not used because nspin < 2'
        endif

        if (peratm .and. nspin.ge.2) then

! Allocate local memory
          nullify(atom,phi,spin,theta)
          call re_alloc( atom, 1, maxat, 'atom', 'initdm' )
          call re_alloc( phi, 1, maxat, 'phi', 'initdm' )
          call re_alloc( spin, 1, maxat, 'spin', 'initdm' )
          call re_alloc( theta, 1, maxat, 'theta', 'initdm' )

          nat = 0
          badsyntax = .FALSE.
          do while(fdf_bline(bfdf,pline) .and. (nat .lt. na_u) .and.  &
                 (.not. badsyntax))

            nn = fdf_bnnames(pline)
            ni = fdf_bnintegers(pline)
            nr = fdf_bnreals(pline)

            if (ni .eq. 1) then
              if (nat .eq. maxat) then
                maxatnew = nat + nint(0.1*nat)
!
                call re_alloc(atom, 1, maxatnew, 'atom', 'initdm',copy=.true.)
!
                call re_alloc( phi, 1, maxatnew, 'phi', 'initdm',    &
                               copy=.true. )
                call re_alloc( spin, 1, maxatnew, 'spin', 'initdm',  &
                               copy=.true. )
                call re_alloc( theta, 1, maxatnew, 'theta', 'initdm', &
                               copy=.true. )
!
                maxat = maxatnew
              endif
              nat = nat + 1
              atom(nat) = fdf_bintegers(pline,1)

              if (nn .eq. 0) then
! Read value of spin
                if (nr .eq. 3) then
! Read spin value and direction
                  spin(nat)  = fdf_breals(pline,1)
                  theta(nat) = fdf_breals(pline,2) * pi/180.0d0
                  phi(nat)   = fdf_breals(pline,3) * pi/180.0d0
                elseif (nr .eq. 1) then
! Read spin value. Default direction.
                  spin(nat)  = fdf_breals(pline,1)
                  theta(nat) = 0.d0
                  phi(nat)   = 0.d0
                else
! Print bad-syntax error and stop
                  badsyntax = .TRUE.
                endif
              elseif (nn .eq. 1) then
! Read spin as + or - (maximun value)
                updo = fdf_bnames(pline,1)
                if (updo .eq. '+') then
                  spin(nat) =  100.d0
                elseif (updo .eq. '-') then
                  spin(nat) = -100.d0
                else
! Print bad-syntax error and stop
                  badsyntax = .TRUE.
                endif
                if (nr .eq. 2) then
                  theta(nat) = fdf_breals(pline,1) * pi/180.0d0
                  phi(nat)   = fdf_breals(pline,2) * pi/180.0d0
                elseif (nr .eq. 0) then
                  theta(nat) = 0.d0
                  phi(nat)   = 0.d0
                else
! Print bad-syntax error and stop
                  badsyntax = .TRUE.
                endif
              else
! Print bad-syntax error and stop
                badsyntax = .TRUE.
              endif

              if ((atom(nat) .lt. 1) .or. (atom(nat) .gt. na_u)) then
                write(msg,'(a,a,i4)') 'intdm: ERROR: Bad atom ' //    &
                 'index in DM.InitSpin, line', nat+1
                call die(TRIM(msg))
              endif
              if (abs(theta(nat)) .gt. 1.d-12) noncol = .true.
            else
! Print bad-syntax error and stop
              badsyntax = .TRUE.
            endif
          enddo

          if (badsyntax) then
            write(msg,'(a,i4)')   &
             'initdm: ERROR: bad syntax in DM.InitSpin, line', nat+1
            call die(msg)
          endif

          if (noncol .and. nspin.lt.4) then
            if (Node.eq.0) then
            write(6,'(/,2a)') 'initdm: WARNING: noncolinear spins ',  &
                     'in DM.InitSpin not used because nspin < 4'
            endif
            noncol = .false.
          endif

! Initialize to 0

          Dscf(:,1:nspin) = 0.0_dp

! Initialize all paramagnetic

          do ia = 1, na_u
            do io = lasto(ia-1) + 1, lasto(ia)
              call GlobalToLocalOrb(io,Node,Nodes,iio)
              if (iio.gt.0) then
                do in = 1, numh(iio)
                  ind = listhptr(iio)+in
                  jo = listh(ind)
                  if (io .eq. jo) then
                    Dscf(ind,1) = 0.5d0 * Datm(io)
                    Dscf(ind,2) = Dscf(ind,1)
                  endif
                enddo
              endif
            enddo
          enddo

! Loop on atoms with spin

          do iat = 1, nat
            ia = atom(iat)

! Find maximum atomic moment that the atoms involved can carry

            spinat = 0.d0
            do io = lasto(ia-1) + 1, lasto(ia)
              spinat = spinat + min( Datm(io), 2.d0 - Datm(io) )
            enddo
            if (spinat.lt.epsilon .and. Node.eq.0) print'(a,i6,a)',  &
             'initdm: WARNING: atom ', atom(iat),                    &
             ' has a closed-shell and cannot be polarized'

! If given spin is larger than possible, make it to max atomic

            aspin = abs(spin(iat))
            if ((aspin .gt. spinat) .and. (aspin .gt. epsilon))    &
                spin(iat) = spinat*spin(iat)/aspin

! Initialize orbitals with same rate as atom

            rate = spin(iat) / (spinat+epsilon)
            do io = lasto(ia-1) + 1, lasto(ia)
              call GlobalToLocalOrb(io,Node,Nodes,iio)
              if (iio.gt.0) then
                qio = Datm(io)
                spio = rate * min( Datm(io), 2.d0 - Datm(io) )
                do in = 1, numh(iio)
                  ind = listhptr(iio)+in
                  jo = listh(ind)
                  if (io .eq. jo) then
                    if (noncol) then
! Store non-collinear-spin density matrix as
!   ispin=1 => D11, ispin=2 => D22;
!   ispin=3 => Real(D12); ispin=4 => Imag(D12)
                      costh = cos(theta(iat))
                      sinth = sin(theta(iat))
                      cosph = cos(phi(iat))
                      sinph = sin(phi(iat))
                      Dscf(ind,1) = (qio + spio * costh) / 2
                      Dscf(ind,2) = (qio - spio * costh) / 2
                      Dscf(ind,3) =   spio * sinth * cosph / 2
                      Dscf(ind,4) = - spio * sinth * sinph / 2
                    else
                      Dscf(ind,1) = (qio + spio) / 2
                      Dscf(ind,2) = (qio - spio) / 2
                    endif
                  endif
                enddo
              endif
            enddo

          enddo

! Deallocate local memory
          call de_alloc( atom, 'atom', 'initdm' )
          call de_alloc( phi, 'phi', 'initdm' )
          call de_alloc( spin, 'spin', 'initdm' )
          call de_alloc( theta, 'theta', 'initdm' )

! ---------------------------------------------------------------------

        else

! Initialize to 0
          Dscf(:,1:nspin) = 0.0d0

! Automatic, for non magnetic (nspin=1) or for Ferro or Antiferro -----
          do io = 1, no_l
            call LocalToGlobalOrb(io,Node,Nodes,iio)
            do in = 1,numh(io)
              ind = listhptr(io)+in
              jo = listh(ind)
              if (iio .eq. jo) then
                if (nspin .eq. 1) then

! No spin polarization

                  Dscf(ind,1) = Datm(iio)
                else

! Spin polarization

                  i1 = 1
                  i2 = 2

! Ferro or antiferro according to DM.InitSpinAF (inspn)

                  if (inspn) then
                    if (mod(iaorb(iio),2).eq.0) then
                      i1 = 2
                      i2 = 1
                    endif
                  endif
                  Dscf(ind,i1) = min( Datm(iio), 1.d0 )
                  Dscf(ind,i2) = Datm(iio) - Dscf(ind,i1)
                endif
              endif
            enddo
          enddo

        endif

      call print_initial_spin()


      CONTAINS

      subroutine print_initial_spin()
      use m_mpi_utils, only: Globalize_sum
      use sparse_matrices, only: S

      real(dp) :: qspin(nspin)
      integer  :: io, j, ispin, ind
#ifdef MPI
      real(dp) :: qtmp(nspin)
#endif

! Print spin polarization
      if (nspin .ge. 2) then
        do ispin = 1,nspin
          qspin(ispin) = 0.0_dp
          do io = 1,no_l
            do j = 1,numh(io)
              ind = listhptr(io)+j
              qspin(ispin) = qspin(ispin) + Dscf(ind,ispin) * S(ind)
            enddo
          enddo
        enddo

#ifdef MPI
! Global reduction of spin components
        call globalize_sum(qspin(1:nspin),qtmp(1:nspin))
        qspin(1:nspin) = qtmp(1:nspin)
#endif
        if (nspin .eq. 2) then
           if (IOnode) then
              write(6,'(/,a,f12.6)')   &
                  'initdm: Initial spin polarization (Qup-Qdown) =',  &
                  qspin(1) - qspin(2)
           endif
        endif
      endif
      end subroutine print_initial_spin

      end subroutine fill_dscf_from_atom_info

      subroutine extrapolate_dm_with_coords(DM_history,na_u,xa,sparse_pattern,DMnew)
        use class_Sparsity
        use class_dData2D
        use class_OrbitalDistribution
        use class_dSpData2D
        use class_Geometry
        use class_Pair_Geometry_dSpData2D
        use class_Fstack_Pair_Geometry_dSpData2D

        use fdf, only: fdf_get

        type(Fstack_Pair_Geometry_dSpData2D), intent(in) :: DM_history
	integer, intent(in)                             :: na_u 
        real(dp), intent(in)                            :: xa(:,:)
        type(Sparsity), intent(in)                      :: sparse_pattern
        type(dSpData2D), intent(inout)                   :: DMnew

        integer :: n, i, nspin, nnzs_out
        real(dp), allocatable   :: c(:)
        real(dp), allocatable   :: xan(:,:,:), dummy_cell(:,:,:)
        type(Geometry), pointer :: geom 
        type(dSpData2D), pointer :: dm
        type(OrbitalDistribution), pointer    :: orb_dist
        type(Pair_Geometry_dSpData2D), pointer :: pair

        type(dSpData2D)       :: DMtmp
        type(dData2D)       :: a_out

        real(dp), dimension(:,:) , pointer  :: a, ai, xp

        n = n_items(DM_history)
        allocate(c(n))

        allocate(xan(3,na_u,n),dummy_cell(3,3,n))

        do i = 1, n
           pair => get_pointer(DM_history,i)
           call firstp(pair,geom)
           xp => coords(geom)
           xan(:,:,i) = xp(:,:)
           dummy_cell(:,:,i) = 1.0_dp
        enddo
        if (fdf_get("UseDIISforDMExtrapolation",.true.)) then
           ! Cast Jose Soler's idea into a DIIS framework

           if (fdf_get("UseSVDExperimental",.false.)) then
              ! Attempt to use the "alternate" KF method with
              ! first differences.  It does not work well yet
              call extrapolate_diis_svd_new(na_u,n,dummy_cell,  &
                                        xan,dummy_cell(:,:,1),xa,c)
           else if (fdf_get("UseSVD",.true.)) then
              ! Straightforward SVD
              call extrapolate_diis_svd(na_u,n,dummy_cell,  &
                                        xan,dummy_cell(:,:,1),xa,c)
           else
              call extrapolate_diis(na_u,n,dummy_cell, &
                                    xan,dummy_cell(:,:,1),xa,c)
           endif
        else  
           ! Use Jose Soler's original method
           call extrapolate(na_u,n,dummy_cell,xan,dummy_cell(:,:,1),xa,c)
        endif
        if (ionode) then
           print *, "DM extrapolation coefficients: "
           do i = 1, n
              print "(i0,f10.5)", i, c(i)
           enddo
        endif

        pair => get_pointer(DM_history,1)
        call secondp(pair,dm)

        ! We assume that all DMs in the history stack have the same orbital distribution...
        orb_dist => dist(dm)
        a => val(dm)
        nspin = size(a,dim=2)
        nnzs_out = nnzs(sparse_pattern)

        ! Scratch array to accumulate the elements
        call newdData2D(a_out,nnzs_out, nspin,name="(temp array for extrapolation)")
        a => val(a_out)
        a(:,:) = 0.0_dp

        do i = 1, n
           pair => get_pointer(DM_history,i)
           call secondp(pair,dm)
!           if (.not. associated(orb_dist,dist(dm))) then
!              call die("Different orbital distributions in DM history stack")
!           endif
           call restructdSpData2D(dm,sparse_pattern,DMtmp)
           ai => val(DMtmp)
           a = a + c(i) * ai
        enddo

        call newdSpData2D(sparse_pattern,a_out,orb_dist, &
                         DMnew,name="SpM extrapolated using coords")
        call delete(a_out)
        call delete(DMtmp)
        deallocate(xan,c,dummy_cell)

      end subroutine extrapolate_dm_with_coords

!
!
      subroutine get_allowed_history_depth(n)
      ! 
      ! Encapsulates the logic of DM extrapolation and re-use
      ! (work in progress)

        use siesta_options, only: DM_history_depth, harrisfun, idyn

      integer, intent(out)         :: n

      n = DM_history_depth     ! As set by default or by the user

      if (harrisfun) then
         n = 0
         if (ionode) print "(a)", &
          "DM_history_depth set to zero for 'Harris' run"
         return
      else if (.not. fdf_get("DM.AllowReuse",.true.)) then
         n = 0
         if (ionode) print "(a)",   &
            "DM_history_depth set to zero since no re-use is allowed"
         return
      else if (.not. fdf_get("DM.AllowExtrapolation",.true.)) then
         n = 1
         if (ionode) print "(a)", &
         "DM_history_depth set to one since no extrapolation is allowed"
         return
      else if (idyn .eq. 0)  then   ! Geometry relaxation
	 if (fdf_get("DM.AllowExtrapolation",.false.)) then
           if (ionode) print "(a,i0)", "Requested Extrapolation for geometry relaxation. DM_history_depth: ", n
         else	
           n = 1
           if (ionode) print "(a)", &
            "DM_history_depth set to one: no extrapolation allowed by default for geometry relaxation"
         endif
         return
      else if ((idyn .eq. 6) .OR. (idyn .eq. 7))  then   ! Forces or Phonon series
         n = 0
         if (ionode) print "(a)", &
         "DM_history_depth set to zero for 'Forces' run"
         return
      endif
      end subroutine get_allowed_history_depth


! *******************************************************************
! SUBROUTINE extrapolate( na, n, cell, xa, cell0, x0, c )
! *******************************************************************
! Finds optimal coefficients for an approximate expasion of the form
! D_0 = sum_i c_i D_i, where D_i is the density matrix in the i'th
! previous iteration, or any other function of the atomic positions. 
! In practice, given points x_0 and x_i, i=1,...,n, it finds the
! coefficients c_i such that sum_i c_i*x_i, is closest to x_0, with
! sum_i c_i = 1. Unit cell vectors are included for completeness of
! the geometry specification only. Though not used in this version,
! this routine can be used if cell vectors change (see algorithms).
! Written by J.M.Soler. Feb.2010.
! ************************* INPUT ***********************************
!  integer  na          ! Number of atoms
!  integer  n           ! Number of previous atomic positions
!  real(dp) cell(3,3,n) ! n previous unit cell vectors
!  real(dp) xa(3,na,n)  ! n previous atomic positions (most recent first)
!  real(dp) cell0(3,3)  ! Present unit cell vectors
!  real(dp) x0(3,na)    ! Present atomic positions
! ************************* OUTPUT **********************************
!  real(dp) c(n)        ! Expansion coefficients
! ************************ UNITS ************************************
! Unit of distance is arbitrary, but must be the same in all arguments
! ********* BEHAVIOUR ***********************************************
! - Returns without any action if n<1 or na<1
! - Stops with an error message if matrix inversion fails
! ********* DEPENDENCIES ********************************************
! Routines called: 
!   inver     : Matrix inversion
! Modules used:
!   precision : defines parameter 'dp' (double precision real kind)
!   sys       : provides the stopping subroutine 'die'
! ********* ALGORITHMS **********************************************
! Using a linear approximation for D(x), imposing that sum(c)=1, and
! assuming that <dD/dx_i*dD/dx_j>=0, <(dD/dx_i)**2>=<(dD/dx_j)**2>
! (i.e. assuming no knowledge on the parcial derivatives), we find
! (D(x0)-D(sum_i c_i*x_i))**2 = const * (x0-sum_i c_i*x_i)**2
! Therefore, given points x_0 and x_i, i=1,...,n, we look for the
! coefficients c_i such that sum_i c_i*x_i, is closest to x_0, with
! sum_i c_i = 1. It is straighforward to show that this reduces to 
! solving S*c=s0, where S_ij=x_i*x_j and s0_i=x_i*x0, under the
! constraint sum(c)=1. To impose it, we rewrite the expansion as
! xmean + sum_i c_i*(x_i-xmean), where xmean=(1/n)*sum_i x_i.
! Since the vectors (x_i-xmean) are linearly dependent, the matrix
! S_ij=(x_i-xmean)*(x_j-xmean) is singular. Therefore, we substitute
! the last row of S and s0 by 1, thus imposing sum(c)=1.
! Unit cell vectors are not used in this version, but this routine
! can be safely used even if cell vectors change, since D depends
! on the interatomic distances, and it can be shown that
! sum_ia,ja (xa(:,ia,i)-xa(:,ja,i))*(xa(:,ia,j)-xa(:,ja,j))
!   = 2*na*sum_ia (xa(:,ia,i)-xmean(:,ia))*(xa(:,ia,j)-xmean(:,ia))
!   = 2*na*S_ij
! This implies that approximating the present ineratomic distances,
! as an expansion of previous ones, is equivalent to approximating
! the atomic positions.
! *******************************************************************

SUBROUTINE extrapolate( na, n, cell, xa, cell0, x0, c )

! Used procedures and parameters
  USE sys,       only: message           ! Termination routine
  USE precision, only: dp            ! Double precision real kind
  use parallel,  only: Node

! Passed arguments
  implicit none
  integer, intent(in) :: na          ! Number of atoms
  integer, intent(in) :: n           ! Number of previous atomic positions
  real(dp),intent(in) :: cell(3,3,n) ! n previous unit cell vectors
  real(dp),intent(in) :: xa(3,na,n)  ! n previous atomic positions
  real(dp),intent(in) :: cell0(3,3)  ! Present unit cell vectors
  real(dp),intent(in) :: x0(3,na)    ! Present atomic positions
  real(dp),intent(out):: c(n)        ! Expansion coefficients

! Internal variables and arrays
  integer :: i, ierr, j, m, ix, ia
  real(dp):: s(n,n), s0(n), si(n,n), xmean(3,na)

! Trap special cases
  if (na<1 .or. n<1) then
    return
  else if (n==1) then
    c(1) = 1
    return
  end if

! Find average of previous positions
  do ia=1,na
     do ix=1,3
        xmean(ix,ia) = sum(xa(ix,ia,1:n)) / n
     enddo
  enddo
 
! The above is equivalent to
!   xmean = sum(xa,dim=3)/n

! Find matrix s of dot products. Subtract xmean to place origin within the
! hyperplane of x vectors

  do j = 1,n
    s0(j) = sum( (x0-xmean) * (xa(:,:,j)-xmean) )
    do i = j,n
      s(i,j) = sum( (xa(:,:,i)-xmean) * (xa(:,:,j)-xmean) )
      s(j,i) = s(i,j)
    end do
  end do

! Find the largest number of (first) m linearly independent vectors xa-xmean.
! Notice that m<n because we have subtracted xmean. Optimally, we should not
! restrict ourselves to the first (most recent) m vectors, but that would
! complicate the code too much.
  do m = n-1,0,-1
    if (m==0) exit                   ! Do not try to invert a 0x0 'matrix'
    call inver( s, si, m, n, ierr )
    if (ierr==0) exit                ! Not too elegant, but it will do.
  end do

  if (node==0) print *, " # of linearly independent vectors: ", m

! Trap the case in which the first two xa vectors are equal
  if (m==0) then  ! Just use most recent point only
    c(n) = 1
    c(1:n-1) = 0
    return
  end if

! Set one more row for equation sum(c)=1.
  m = m+1
  s0(m) = 1
  s(m,1:m) = 1

! Invert the full equations matrix
  call inver( s, si, m, n, ierr )
  if (ierr/=0) then
     c(:) = 0.0_dp
     c(n) = 1.0_dp
     call message('extrapolate: matrix inversion failed')
     call message('extrapolate: using last item in history')
     return
  endif

! Find expansion coefficients
! This is wrong regarding the order...
  c(1:m) = matmul( si(1:m,1:m), s0(1:m) )
  c(m+1:n) = 0

END SUBROUTINE extrapolate

! *******************************************************************
! SUBROUTINE extrapolate_diis( na, n, cell, xa, cell0, x0, c )
! *******************************************************************
! Finds optimal coefficients for an approximate expasion of the form
! D_0 = sum_i c_i D_i, where D_i is the density matrix in the i'th
! previous iteration, or any other function of the atomic positions. 
! In practice, given points x_0 and x_i, i=1,...,n, it finds the
! coefficients c_i such that sum_i c_i*x_i, is closest to x_0, with
! sum_i c_i = 1. Unit cell vectors are included for completeness of
! the geometry specification only. Though not used in this version,
! this routine can be used if cell vectors change (see algorithms).
! Original version Written by J.M.Soler. Feb.2010.
! Couched in DIIS language by A. Garcia, Nov. 2012.
! ************************* INPUT ***********************************
!  integer  na          ! Number of atoms
!  integer  n           ! Number of previous atomic positions
!  real(dp) cell(3,3,n) ! n previous unit cell vectors
!  real(dp) xa(3,na,n)  ! n previous atomic positions (most recent first)
!  real(dp) cell0(3,3)  ! Present unit cell vectors
!  real(dp) x0(3,na)    ! Present atomic positions
! ************************* OUTPUT **********************************
!  real(dp) c(n)        ! Expansion coefficients
! ************************ UNITS ************************************
! Unit of distance is arbitrary, but must be the same in all arguments
! ********* BEHAVIOUR ***********************************************
! - Returns without any action if n<1 or na<1
! - Stops with an error message if matrix inversion fails
! ********* DEPENDENCIES ********************************************
! Routines called: 
!   inver     : Matrix inversion
! Modules used:
!   precision : defines parameter 'dp' (double precision real kind)
!   sys       : provides the stopping subroutine 'die'
! ********* ALGORITHMS **********************************************

SUBROUTINE extrapolate_diis( na, n, cell, xa, cell0, x0, c )

! Used procedures and parameters
  USE sys,       only: message, die
  USE precision, only: dp            ! Double precision real kind
  use parallel,  only: Node

! Passed arguments
  implicit none
  integer, intent(in) :: na          ! Number of atoms
  integer, intent(in) :: n           ! Number of previous atomic positions
  real(dp),intent(in) :: cell(3,3,n) ! n previous unit cell vectors
  real(dp),intent(inout) :: xa(3,na,n)  ! n previous atomic positions
  real(dp),intent(in) :: cell0(3,3)  ! Present unit cell vectors
  real(dp),intent(in) :: x0(3,na)    ! Present atomic positions
  real(dp),intent(out):: c(n)        ! Expansion coefficients

! Internal variables and arrays
  integer :: i, info, j, m, ix, ia
  real(dp), allocatable, dimension(:,:) :: s, si

  allocate (s(n+1,n+1), si(n+1,n+1))

! Trap special cases
  if (na<1 .or. n<1) then
    return
  else if (n==1) then
    c(1) = 1
    return
  end if

! Find residuals with respect to x0 
  do i = 1, n
  do ia=1,na
     do ix=1,3
        xa(ix,ia,i) = xa(ix,ia,i) - x0(ix,ia)
     enddo
  enddo
  enddo
 
! Find matrix s of dot products.

  do j = 1,n
    do i = j,n
      s(i,j) = sum( xa(:,:,i) * xa(:,:,j) )
      s(j,i) = s(i,j)
    end do
    ! Now extend the matrix with ones in an extra column
    ! and row ...
    s(j,n+1)=1.0_dp   ! This value is really arbitrary
    s(n+1,j)=1.0_dp   ! This represents the Sum_{Ci} = 1 constraint
  end do
  s(n+1,n+1) = 0.0_dp

  if (Node == 0) then
     call print_mat(s,n+1)
  endif

  ! Find rank of the xi*xj matrix
  do m = n,0,-1
    if (m==0) exit                   ! Do not try to invert a 0x0 'matrix'
    call inver( s, si, m, n+1, info )
    if (info==0) exit                ! Not too elegant, but it will do.
  end do
  if (node==0) print *, " Rank of DIIS matrix: ", m

  if (m<n) then
     info = -1
  else
     ! Invert the full matrix
     call inver( s, si, n+1, n+1, info)
  endif
  c(:) = 0.0_dp

    ! If inver was successful, get coefficients for DIIS/Pulay mixing
    ! (Last column of the inverse matrix, corresponding to solving a
    ! linear system with (0,0,0,...,0,1) in the right-hand side)
    if (info .eq. 0) then
       do i=1,n
          c(i)=si(i,n+1)
       enddo
    else
       ! Otherwise, use only last step
       if (Node == 0) then
         write(6,"(a,i5)")  &
         "Warning: unstable inversion in DIIS - use last item only"
       endif
       c(n) = 1.0_dp
    endif

    deallocate(s,si)
end SUBROUTINE extrapolate_diis

SUBROUTINE extrapolate_diis_svd( na, n, cell, xa, cell0, x0, c )

! Used procedures and parameters
  USE sys,       only: message, die
  USE precision, only: dp            ! Double precision real kind
  use parallel,  only: Node
  use m_svd,     only: solve_with_svd

! Passed arguments
  implicit none
  integer, intent(in) :: na          ! Number of atoms
  integer, intent(in) :: n           ! Number of previous atomic positions
  real(dp),intent(in) :: cell(3,3,n) ! n previous unit cell vectors
  real(dp),intent(inout) :: xa(3,na,n)  ! n previous atomic positions
  real(dp),intent(in) :: cell0(3,3)  ! Present unit cell vectors
  real(dp),intent(in) :: x0(3,na)    ! Present atomic positions
  real(dp),intent(out):: c(n)        ! Expansion coefficients

! Internal variables and arrays
  integer :: i, info, j, m, ix, ia, rank
  real(dp), allocatable, dimension(:,:) :: s, si
  real(dp), allocatable, dimension(:)   :: rhs, sigma, beta

  allocate (s(n+1,n+1), si(n+1,n+1), sigma(n+1), rhs(n+1))
  allocate (beta(n+1))  ! full set of coefficients

! Trap special cases
  if (na<1 .or. n<1) then
    return
  else if (n==1) then
    c(1) = 1
    return
  end if

! Find residuals with respect to x0 
  do i = 1, n
  do ia=1,na
     do ix=1,3
        xa(ix,ia,i) = xa(ix,ia,i) - x0(ix,ia)
     enddo
  enddo
  enddo
 
! Find matrix s of dot products.

  do j = 1,n
    do i = j,n
      s(i,j) = sum( xa(:,:,i) * xa(:,:,j) )
      s(j,i) = s(i,j)
    end do
    ! Now extend the matrix with ones in an extra column
    ! and row ...
    s(j,n+1)=1.0_dp   ! This value is really arbitrary
    s(n+1,j)=1.0_dp   ! This represents the Sum_{Ci} = 1 constraint
  end do
  s(n+1,n+1) = 0.0_dp

  if (Node == 0) then
    ! call print_mat(s,n+1)
  endif

  ! Find rank of the xi*xj matrix
  do m = n,0,-1
    if (m==0) exit                   ! Do not try to invert a 0x0 'matrix'
    call inver( s, si, m, n+1, info )
    if (info==0) exit                ! Not too elegant, but it will do.
  end do
  if (node==0) print *, " Estimated Rank of xi*xj matrix: ", m

  rhs(1:n) = 0.0_dp
  rhs(n+1) = 1.0_dp
  !
  call solve_with_svd(s,rhs,beta,info,sigma=sigma,rank_out=rank)
  !
  if (node==0) then
     print *, 'matrix rank: ', rank
     print "(a,/,6g12.5)", 'Singular values: ', sigma(:)
  endif
  if (info == 0) then
     c(1:n) = beta(1:n)
  else
     if (node==0) then
        print *, 'The SVD algorithm failed to converge'
        print *, 'Re-using last item only...'
        print *, "info: ", info
     endif
     c(:) = 0.0_dp
     c(n) = 1.0_dp
  endif

  deallocate(s,si,sigma,rhs,beta)
end SUBROUTINE extrapolate_diis_svd

SUBROUTINE extrapolate_diis_svd_new( na, n, cell, xa, cell0, x0, c )

! Used procedures and parameters
  USE sys,       only: message, die
  USE precision, only: dp            ! Double precision real kind
  use parallel,  only: Node
  use m_svd,     only: solve_with_svd

! Passed arguments
  implicit none
  integer, intent(in) :: na          ! Number of atoms
  integer, intent(in) :: n           ! Number of previous atomic positions
  real(dp),intent(in) :: cell(3,3,n) ! n previous unit cell vectors
  real(dp),intent(inout) :: xa(3,na,n)  ! n previous atomic positions
  real(dp),intent(in) :: cell0(3,3)  ! Present unit cell vectors
  real(dp),intent(in) :: x0(3,na)    ! Present atomic positions
  real(dp),intent(out):: c(n)        ! Expansion coefficients

! Internal variables and arrays
  integer :: i, info, j, m, ix, ia, rank
  real(dp), allocatable, dimension(:,:) :: s
  real(dp), allocatable, dimension(:)   :: rhs, sigma
  real(dp), allocatable, dimension(:)   :: beta ! intermediate coeffs

  allocate (s(n-1,n-1), rhs(n-1), beta(n-1), sigma(n-1))

! Trap special cases
  if (na<1 .or. n<1) then
    return
  else if (n==1) then
    c(1) = 1
    return
  end if

! Find residual of Nth term with respect to x0 
  do ia=1,na
     do ix=1,3
        xa(ix,ia,n) = xa(ix,ia,n) - x0(ix,ia)
     enddo
  enddo

! Find first differences
! As in alternate method in KF
  do j=1,n-1
     do ia=1,na
        do ix=1,3
           xa(ix,ia,j) = xa(ix,ia,j+1) - xa(ix,ia,j)
        enddo
     enddo
  enddo
 
! Find matrix s of dot products.

  do j = 1,n-1
    do i = j,n-1
      s(i,j) = sum( xa(:,:,i) * xa(:,:,j) )
      s(j,i) = s(i,j)
    end do
    ! And right-hand side
    rhs(j) = -sum( xa(:,:,j) * xa(:,:,n) )
  end do

  if (Node == 0) then
    ! call print_mat(s,n+1)
  endif

  call solve_with_svd(s,rhs,beta,info,rank_out=rank,sigma=sigma)
  if (node==0) then
     print *, 'matrix rank: ', rank
     print "(a,/,6g12.5)", 'Singular values: ', sigma(:)
     print "(a,/,6g12.5)", 'Beta coeffs: ', beta(:)
  endif
  if (info /= 0) then
     if (node==0) then
        print *, 'The SVD algorithm failed to converge'
        print *, 'Re-using last item only...'
        print *, "info: ", info
     endif
     c(:) = 0.0_dp
     c(n) = 1.0_dp
  endif

  ! Re-shuffle coefficients
  c(n) = 1.0_dp
  c(1) = -beta(1)
  do j = 2, n-1
     c(j) = beta(j-1)-beta(j)
  enddo

  deallocate(s,beta,rhs,sigma)
end SUBROUTINE extrapolate_diis_svd_new

    subroutine print_mat(a,n)
      integer, intent(in)  :: n
      real(dp), intent(in) :: a(n,n)
      integer :: i, j

      print *, "mat:"
      do i = 1, n
         print "(6g15.7)", (a(i,j),j=1,n)
      enddo
      print *, "------"
    end subroutine print_mat

    SUBROUTINE inverse(A,B,N,NDIM,INFO,debug_inverse)
      use parallel, only: node

      IMPLICIT NONE
      INTEGER, intent(in) ::  N,NDIM
      real(dp), intent(in)  ::  A(NDIM,NDIM)
      real(dp), intent(out) ::  B(NDIM,NDIM)
      integer, intent(out) :: info
      logical, intent(in)  :: debug_inverse

      real(dp)  ::  X

!!$C Routine to obtain the inverse of a general, nonsymmetric
!!$C square matrix, by calling the appropriate LAPACK routines
!!$C The matrix A is of order N, but is defined with a 
!!$C size NDIM >= N.
!!$C If the LAPACK routines fail, try the good old Numerical Recipes
!!$C algorithm
!!$C
!!$C P. Ordejon, June 2003
!!$C
!!$C **** INPUT ****
!!$C A(NDIM,NDIM)   real*8      Matrix to be inverted
!!$C N              integer     Size of the matrix
!!$C NDIM           integer     Defined size of matrix
!!$C **** OUTPUT ****
!!$C B(NDIM,NDIM)   real*8      Inverse of A
!!$C INFO           integer     If inversion was unsucessful, 
!!$C                            INFO <> 0 is returned
!!$C                            Is successfull, INFO = 0 returned
!!$C ***************


      real(dp) ::  C,ERROR,DELTA,TOL
      real(dp), allocatable:: work(:), pm(:,:)
      integer, dimension(:), allocatable  ::  IWORK, ipiv
      INTEGER I,J,K
      real(dp) :: ANORM, RCOND

      logical :: debug
      real(dp), external :: dlange

      allocate (iwork(n), ipiv(n))
      allocate (work(n), pm(n,n))

      debug = debug_inverse .and. (node == 0)

      TOL=1.0D-4
      INFO = 0

      DO I=1,N
      DO J=1,N
        B(I,J)=A(I,J)
      ENDDO
      ENDDO
      ANORM = DLANGE( "1", N, N, B, NDIM, WORK )
      CALL DGETRF(N,N,B,NDIM,IPIV,INFO)

      IF (INFO .NE. 0) THEN
       if (debug) print *,   &
           'inver:  ERROR: DGETRF exited with error message', INFO
        GOTO 100
      ENDIF
      CALL DGECON( "1", N, B, NDIM, ANORM, RCOND, WORK, IWORK, INFO )
      IF (INFO .NE. 0) THEN
       if (debug) print *,   &
           'inver:  ERROR: DGECON exited with error message', INFO
        GOTO 100
      ENDIF
      if (debug) print "(a,g20.8)", "Reciprocal condition number: ", rcond

      CALL DGETRI(N,B,NDIM,IPIV,WORK,N,INFO)

      IF (INFO .NE. 0) THEN
       if (debug) print *,   &
          'inver:  ERROR: DGETRI exited with error message', INFO
        GOTO 100
      ENDIF

! CHECK THAT THE INVERSE WAS CORRECTLY CALCULATED

      pm = matmul(a(1:n,1:n),b(1:n,1:n))

      ERROR=0.0D0
      DO I=1,N
      DO J=1,N
        C=0.0D0
        DO K=1,N
          C=C+A(I,K)*B(K,J)
        ENDDO
        DELTA=0.0D0
        IF (I.EQ.J)  DELTA=1.0D0
        ERROR=ERROR+DABS(C-DELTA)
        C=0.0D0
        DO K=1,N
          C=C+B(I,K)*A(K,J)
        ENDDO
        DELTA=0.0D0
        IF (I.EQ.J)  DELTA=1.0D0
        ERROR=ERROR+DABS(C-DELTA)
      ENDDO
      ENDDO

      IF (ERROR/N .GT. TOL) THEN
       if (debug) then
          print *,   &
          'inver:  ERROR in lapack inverse. Error: ', error
          call print_mat(a,n)
          call print_mat(b,n)
          call print_mat(pm,n)
       endif
        GOTO 100
      ENDIF

      INFO = 0
      deallocate(work,iwork,ipiv,pm)
      RETURN

100   CONTINUE

! Try simple, old algorithm:

      DO I=1,N
        DO J=1,N
          B(I,J)=A(I,J)
        ENDDO
      ENDDO
      DO I=1,N
        IF (B(I,I) .EQ. 0.0D0) THEN
           if (debug) print *,   &
            'inver:  zero pivot in fallback algorithm'
          INFO = -777
          deallocate(work,iwork,ipiv,pm)
          RETURN
        ENDIF
        X=B(I,I)
        B(I,I)=1.0d0
        DO J=1,N
          B(J,I)=B(J,I)/X
        ENDDO
        DO K=1,N
          IF ((K-I).NE.0) THEN 
            X=B(I,K)
            B(I,K)=0.0d0
            DO J=1,N
              B(J,K)=B(J,K)-B(J,I)*X
            ENDDO
          ENDIF
        ENDDO
      ENDDO

! CHECK THAT THE INVERSE WAS CORRECTLY CALCULATED

      pm = matmul(a(1:n,1:n),b(1:n,1:n))
      ERROR=0.0D0
      DO I=1,N
      DO J=1,N
        C=0.0D0
        DO K=1,N
          C=C+A(I,K)*B(K,J)
        ENDDO
        DELTA=0.0D0
        IF (I.EQ.J)  DELTA=1.0D0
        ERROR=ERROR+DABS(C-DELTA)
        C=0.0D0
        DO K=1,N
          C=C+B(I,K)*A(K,J)
        ENDDO
        DELTA=0.0D0
        IF (I.EQ.J)  DELTA=1.0D0
        ERROR=ERROR+DABS(C-DELTA)
      ENDDO
      ENDDO

      IF (ERROR/N .GT. TOL) THEN
       if (debug) then
          print *,   &
            'inver:  INVER unsuccessful, ERROR = ', ERROR
          call print_mat(a,n)
          call print_mat(b,n)
          call print_mat(pm,n)
       endif
        INFO = -888
        deallocate(work,iwork,ipiv,pm)
        RETURN
      ENDIF

      INFO = 0
      deallocate(work,iwork,ipiv,pm)
      RETURN

    end SUBROUTINE inverse


      end module m_new_dm
