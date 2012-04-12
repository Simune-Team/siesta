      MODULE m_new_dm

!     Prepares a starting density matrix for a new geometry iteration
!     This DM can be:
!     1. Synthesized directly from atomic occupations (not idempotent)
!     2. Read from file
!     3. Re-used (with possible extrapolation) from previous geometry step(s).
!
!     In cases 2 and 3, a check is done to guarantee that the structure
!     of the read or extrapolated DM conforms to the current sparsity.
!     If it does not, the information is re-arranged.
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

      CONTAINS

!=====================================================================

      subroutine  new_dm( auxchanged, DM_history, DMnew)

      USE siesta_options
      use siesta_geom,      only: xa, na_u
      use sparse_matrices,  only: sparse_pattern, block_dist
      use atomlist,         only: datm, iaorb, lasto, no_u, no_l
      use m_steps,          only: istp
      use m_spin,   only: nspin

      use class_SpMatrix
      use class_Sparsity
      use class_Pair_Geometry_SpMatrix
      use class_Fstack_Pair_Geometry_SpMatrix

      implicit none

      logical, intent(in) :: auxchanged ! Has auxiliary supercell changed?
      type(Fstack_Pair_Geometry_SpMatrix), intent(inout)      :: DM_history
      type(SpMatrix), intent(inout)   :: DMnew

!     Local variables

      logical :: dminit     ! Initialize density matrix?
      logical :: try_to_read_from_file
      integer :: n_dms_in_history

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
      if (initdmaux.and.auxchanged) then
        dminit = .true.
        try_to_read_from_file = .false.
        if (IOnode) then
           write(6,"(a)") "DM history reset as supercell changed."
        endif
        call new(DM_history,DM_history_depth,"(DM history stack)")
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

         if (ionode) print *, "N DMs in history: ", n_dms_in_history
         if (ionode) call print_type(DM_history)
         call extrapolate_dm_with_coords(DM_history,xa,sparse_pattern,DMnew)
         if (ionode)  print *, "DMnew after DM reuse:"
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
      use class_SpMatrix
      use class_OrbitalDistribution
      use class_Array2D

      implicit          none

      logical           found, inspn, try_dm_from_file
      integer           no_l, na_u, no_u, nspin
      integer           lasto(0:na_u), iaorb(no_u)
      real(dp)          Datm(no_u)

      type(SpMatrix), intent(inout)      :: DMnew
      type(Sparsity), intent(in) :: sparse_pattern
      type(OrbitalDistribution), intent(in) :: block_dist

! ---------------------------------------------------------------------

      character(len=*),parameter:: myName = 'initdm'

      integer :: nspin_read
      real(dp), pointer              :: Dscf(:,:)
      integer, pointer, dimension(:) :: numh, listhptr, listh
      type(SpMatrix)                 :: DMread
      type(Array2D)                  :: dm_a2d

! Try to read DM from disk if wanted (DM.UseSaveDM true) ---------------

      found = .false.
      if (try_dm_from_file) then
         if (ionode) print *, "Attempting to read DM from file..."
         call readSpMatrix(trim(slabel)//".DM",   &
                           DMread,found,block_dist)
      endif

! If found, check and update, otherwise initialize with neutral atoms

      if (found) then
        ! Various degrees of sanity checks

        nspin_read = size(val(DMread),dim=2)
        if (nspin_read /= Nspin) then
           if (IOnode) then
              write(6,"(a,i6,/,a)")                   &
              "WARNING: Wrong nspin in DM file: ",  nspin_read,  &
              "WARNING: Falling back to atomic initialization of DM."
           endif
           found = .false.
	   call delete(DMread)
        endif

        if (nrows_g(DMread) /= nrows_g(sparse_pattern)) then
           if (IONode) then
              write(6,"(a,/,a)")                             &
             "WARNING: Wrong number of orbs in DM file. ",     &
             "WARNING: Falling back to atomic initialization of DM."
           endif
           found = .false.
	   call delete(DMread)
        endif

      endif
      
      if (found) then

	call restructSpMatrix(DMread,sparse_pattern,DMnew)
        if (ionode) print *, "DMread after reading file:"
        if (ionode) call print_type(Dmread)
        call delete(DMread)

      else

	call newArray2D(dm_a2d,nnzs(sparse_pattern),nspin,"(DMatomic)")
	Dscf => val(dm_a2d)
        numh     => n_col(sparse_pattern)
        listhptr => list_ptr(sparse_pattern)
        listh    => list_col(sparse_pattern)
	
        call   fill_dscf_from_atom_info(Datm, Dscf,              &
                        numh, listhptr, listh, lasto,         &
                        no_u, na_u, no_l, nspin,     &
                        iaorb, inspn)

        call newSpMatrix(sparse_pattern,dm_a2d,block_dist,DMnew,  &
                         "(DM initialized from atoms)")
        call delete(dm_a2d)
        if (ionode) print *, "DMnew after filling with atomic data:"
        if (ionode) call print_type(DMnew)

       endif
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

      end subroutine fill_dscf_from_atom_info

      subroutine extrapolate_dm_with_coords(DM_history,xa,sparse_pattern,DMnew)
        use class_Sparsity
        use class_Array2D
        use class_OrbitalDistribution
        use class_SpMatrix
        use class_Geometry
        use class_Pair_Geometry_SpMatrix
        use class_Fstack_Pair_Geometry_SpMatrix

        type(Fstack_Pair_Geometry_SpMatrix), intent(in) :: DM_history
        real(dp), intent(in)                            :: xa(:,:)
        type(Sparsity), intent(in)                      :: sparse_pattern
        type(SpMatrix), intent(inout)                   :: DMnew

        integer :: n, na, i, nspin, nnzs_out
        real(dp), allocatable   :: c(:)
        real(dp), allocatable   :: xan(:,:,:), dummy_cell(:,:,:)
        type(Geometry), pointer :: geom 
        type(SpMatrix), pointer :: dm
        type(OrbitalDistribution), pointer    :: orb_dist
        type(Pair_Geometry_SpMatrix), pointer :: pair

        type(SpMatrix)       :: DMtmp
        type(Array2D)        :: a_out

        real(dp), dimension(:,:) , pointer  :: a, ai, xp

        n = n_items(DM_history)
        allocate(c(n))

        na = size(xa,dim=2)
        allocate(xan(3,na,n),dummy_cell(3,3,n))

        do i = 1, n
           pair => get_pointer(DM_history,i)
           call firstp(pair,geom)
           xp => coords(geom)
           xan(:,:,i) = xp(:,:)
           dummy_cell(:,:,i) = 1.0_dp
        enddo
        call extrapolate(na,n,dummy_cell,xan,dummy_cell(:,:,1),xa,c)
        print *, "Extrapolation coefficients: "
        do i = 1, n
           print "(i0,f10.6)", i, c(i)
        enddo

        pair => get_pointer(DM_history,1)
        call secondp(pair,dm)

        ! We assume that all DMs in the history stack have the same orbital distribution...
        orb_dist => dist(dm)
        a => val(dm)
        nspin = size(a,dim=2)
        nnzs_out = nnzs(sparse_pattern)

        ! Scratch array to accumulate the elements
        call newArray2D(a_out,nnzs_out, nspin,name="(temp array for extrapolation)")
        a => val(a_out)
        a(:,:) = 0.0_dp

        do i = 1, n
           pair => get_pointer(DM_history,i)
           call secondp(pair,dm)
!           if (.not. associated(orb_dist,dist(dm))) then
!              call die("Different orbital distributions in DM history stack")
!           endif
           call restructSpMatrix(dm,sparse_pattern,DMtmp)
           ai => val(DMtmp)
           a = a + c(i) * ai
        enddo

        call newSpMatrix(sparse_pattern,a_out,orb_dist, &
                         DMnew,name="SpM extrapolated using coords")
        call delete(a_out)
        call delete(DMtmp)

      end subroutine extrapolate_dm_with_coords

      end module m_new_dm
