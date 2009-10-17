      MODULE m_new_dm

!     Prepares a starting density matrix for a new geometry iteration
!     This DM can be:
!     1. Synthesized directly from atomic occupations (not idempotent)
!     2. Read from file
!     3. Extrapolated from (two) previous steps
!     3.a The DM of the previous iteration
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
!
!      For the purposes of extrapolation, this module keeps DMsaved, which
!      holds the DM from a previous geometry step (used, together with the
!      current DM to extrapolate).
!      The module also keeps the sparsity at the time of initialization,
!      so that it can be compared to the current one. If restructuring is
!      needed, the saved sparsity is updated.
!
!     Alberto Garcia, September 2007
!
      use m_sparse,  only : same_sparsity, change_sparsity_inplace
      use sys,       only : die
      use precision, only : dp
      use alloc,     only : re_alloc, de_alloc
      use parallel,  only : IOnode
      use m_norm
#ifdef MPI
      use mpi_siesta
#endif
      implicit none
      real(dp), pointer, private :: DMsaved(:,:)
      integer, pointer, private  ::
     &                   listold(:), listptrold(:), numold(:)

      private
      public :: new_dm
!      real(dp) :: anorm

      CONTAINS

      subroutine new_dm( auxchanged )
      ! We should probably pass it Dscf and the rest....
      use sparse_matrices
      USE siesta_options
      use siesta_geom
      use atomlist, only : datm, no_s, iaorb, lasto, no_u, no_l
      use m_spin,   only : nspin
      use m_steps,  only : istp
#ifdef TRANSIESTA
!      TSS Begin
      use m_energies,   only : ef
      use sys,          only : die
      use m_ts_options, only : VIn, VOut, TSinit, TSrun, foundts, TSmode
      use m_ts_io,      only : TSiodm 
!     TSS End 
#endif /* TRANSIESTA */
      implicit none
      logical, intent(in) :: auxchanged ! Has auxiliary supercell changed?
      logical             :: dminit     ! Initialize density matrix?
      logical             :: found, try_to_read_from_file

#ifdef TRANSIESTA
!     TSS Begin
      integer         :: ts_alloc_stat
      logical, save   :: first_time = .true.
!     TSS End
#endif /* TRANSIESTA */

      if (IOnode) then
         write(6,"(a,i5)") "New_DM. Step: ", istp
      endif

!     Initialize density matrix
      if (istp == 1) then
         nullify(numold,listold,listptrold)
         nullify(DMsaved)
      endif
!
!     In principle we allow the re-use of the DM (i.e, we do not initialize it)
!
      dminit = .false. 
      try_to_read_from_file = usesavedm     ! As per defaults
!
!     Except if there are explicit instructions
!
      if (.not. allow_dm_reuse) then
         dminit = .true.
         try_to_read_from_file = .false.  ! Just in case the user has a fossil DM.UseSaveDM
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
!
      if (istp .eq. 1) then
         dminit = .true.
      else
         if ((idyn .eq. 6) ! Force Constants
     $         .and. usesavedm .and. writedm)  dminit = .true.
         if ((idyn .eq. 7) ! Phonon series (writedm??)
     $         .and. usesavedm)  dminit = .true.
      endif

!
!     ... or if the auxiliary cell has changed
!     (in this case we have to  avoid reading back saved copy from file)
!
      if (initdmaux.and.auxchanged) then
        dminit = .true.
        try_to_read_from_file = .false.
        if (IOnode) then
           write(6,"(a)") "DM reset as supercell changed."
        endif
      endif

#ifdef TRANSIESTA
!     TSS Begin: Mixing arrays
      if (first_time) then
        nullify(VIn,VOut)
        first_time = .false.
      endif
      call re_alloc( VIn,  1, maxnh, 1, nspin, 'Vin',  'new_dm' )
      call re_alloc( Vout, 1, maxnh, 1, nspin, 'Vout', 'new_dm' )
!     TSS End
#endif /* TRANSIESTA */

      if (dminit) then
         if (IOnode) then
            write(6,"(a)") "Initializing Density Matrix..."
         endif

         ! Reset DMsaved to null
         call de_alloc(DMsaved)
         nullify(DMsaved)
         !
         call initdm(Datm, Dscf, lasto, na_s,
     .               maxnh, no_s, no_l, nspin, na_u, no_l, nspin,
     .               numh, listhptr, listh, iaorb, found, inspn, 
     .               try_to_read_from_file)
         !
         ! Store sparsity
         call re_alloc(numold,1,no_l)
         call re_alloc(listptrold,1,no_l)
         ! NOTE possible shrinking, as the test for changed
         ! sparsity might involve the length of the listh array
         call re_alloc(listold,1,size(listh),shrink=.true.)  
         numold = numh
         listold = listh
         listptrold = listhptr

!         anorm = compute_norm(Dscf,maxnh,nspin)
!         if (IONode) print *, "DM norm after initdm: ", anorm

#ifdef TRANSIESTA
!     TSS Begin
!     TSS check of TSDM stuff
      if (istp .eq. 1 .and. TSmode) then
!       try to read TSDM
        call TSiodm( 'read', maxnh, no_l, nspin, numh, listhptr,
     &              listh, Dscf, Escf, ef, found=foundts )

        if(.not. foundts )  then ! not a TS continuation run
          TSinit = .true.     ! start to converge a diagon run
          TSrun = .false.

          if(IONode) then
             write(6,'(a)')
     &            'TRANSIESTA: No TS-DensityMatrix file found'
             write(6,'(a)')
     &            'TRANSIESTA: Initialization runs using diagon'
          endif
        else
          TSinit = .false.
          TSrun = .true.

          if(IONode) then
             write(6,'(a)') 'TRANSIESTA: Continuation run'
             write(6,'(/a)')
     &            '                     ************************  '
             write(6,'(a)')
     &            '                     *   TRANSIESTA BEGIN   *  '
             write(6,'(a)')
     &            '                     ************************  '
          endif

!   copy Dsfc -> Dold
          Dold(1:nh,1:nspin) = Dscf(1:nh,1:nspin)
          Eold(1:nh,1:nspin) = Escf(1:nh,1:nspin)
        end if ! found TSDM-file
      end if ! istp == 1  and TSmode
!     TSS End
#endif /* TRANSIESTA */
      else

         if (IOnode) then
            write(6,"(a)") "Re-using DM from previous geometry..."
         endif

        ! Extrapolate density matrix between steps

!         anorm =  compute_norm(Dscf,maxnh,nspin)
!         if (IONode) print *, "DM norm of previous geom: ", anorm

         call new_extrapol(no_l, numh, listhptr, listh, Dscf)

!         anorm =  compute_norm(Dscf,maxnh,nspin)
!         if (IONode) print *, "DM norm after extrapolation: ", anorm

      endif

      ! Initialize energy-density matrix to zero for first call to overfsm
      Escf(:,:) = 0.0_dp

      END subroutine new_dm

!---------------------------------------------------

      subroutine new_extrapol(nrows,num, listptr,list, DM)

      USE siesta_options

      ! This routine extrapolates DM
      ! If it is called for the first time after initialization
      ! it just passes DM through, possibly re-structured.

      integer, intent(in) :: nrows
      integer, intent(in) :: num(nrows), listptr(nrows)
      integer, intent(in) :: list(:)
      real(dp), pointer   :: DM(:,:)
      

      logical   ::  sparsity_has_changed
      real(dp)  ::  msave
      integer    ::  i, j, ispin, nspin, nnz

#ifdef MPI
      integer  MPIerror
      logical  lbuffer
#endif

      sparsity_has_changed = 
     $  (.not.  same_sparsity (nrows,numold,listold, num,list))

#ifdef MPI
      call MPI_AllReduce(sparsity_has_changed,lbuffer,1,MPI_logical,
     $                   MPI_lor,
     $                   MPI_Comm_World,MPIerror)
      sparsity_has_changed = lbuffer
#endif

      ! Here is an interesting thought:
      ! Can we do the re-structuring node by node, only on those
      ! in which it is necessary? I very much think so, so the 
      ! above block is superfluous.

      ! Another thought: Is it possible that no_l changes during
      ! the calculation (i.e., the distribution of orbitals among
      ! nodes)? If it does, we are toast. ! ***** AG
      ! In that case we might simply request an initialization of
      ! the density matrix, or write a super extrapolator.

      if (sparsity_has_changed) then
   
         ! Fix structure of DM

         call change_sparsity_inplace(nrows,
     $                                numold,listptrold,listold,
     $                                num, listptr, list,
     $                                DM )

      endif

      nspin = size(DM,dim=2)
      nnz = size(DM,dim=1)

      if (associated(DMsaved)) then

         if (IOnode) then
            write(6,"(a)") "Extrapolating Density Matrix..."
         endif

         if (sparsity_has_changed) then
            !Fix structure of DMsaved
            call change_sparsity_inplace(nrows,
     $                                numold,listptrold,listold,
     $                                num, listptr, list,
     $                                DMsaved )
         endif

         ! Extrapolate DM

         do ispin = 1, nspin
            do i = 1, nnz
               msave = DM(i,ispin)
               DM(i,ispin) = 2.0_dp * DM(i,ispin) - DMsaved(i,ispin)
               DMsaved(i,ispin) = msave
            enddo
         enddo

      else

         if (IOnode) then
            write(6,"(a)") "Re-using DM without extrapolation..."
         endif

         ! Copy DM to DMsaved
         ! Note that we use explicit bounds
         ! If we do not request a shrink, DMsaved might end up
         ! being larger than DM 

         if (allow_dm_extrapolation) then
            call re_alloc(DMsaved,1,nnz,1,nspin)
            DMsaved(1:nnz,1:nspin) = DM(1:nnz,1:nspin)
         else
            if (IOnode) then
               write(6,"(a)") "(DM extrapolation not allowed)"
            endif
         endif
      endif

      if (sparsity_has_changed) then

         if (IOnode) then
            write(6,"(a)") "Density Matrix sparsity pattern changed."
         endif

         ! Store sparsity
         call re_alloc(numold,1,size(num),shrink=.true.)
         call re_alloc(listptrold,1,size(listptr),shrink=.true.)
         call re_alloc(listold,1,size(list),shrink=.true.)
         numold = num
         listold = list
         listptrold = listptr

      endif
         
      end subroutine new_extrapol

!
      subroutine initdm(Datm, Dscf, lasto, maxa,
     .                  maxnh, maxo, maxuo, maxspn, nua, no_l, nspin, 
     .                  numh, listhptr, listh, iaorb, found, inspn, 
     .                  try_dm_from_file)

c *******************************************************************
c Density matrix initialization
c
c    If Try_Dm_From_File is true, it is read from file if present.
c    Otherwise it is generated assuming atomic charging 
c      (filling up atomic orbitals). The DM originated that way is
c      not a good DM due to overlaps, but the SCF cycling corrects
c      that for the next cycle.
c    Spin polarized calculations starting from atoms:
c      Default: All atoms with maximum polarization compatible with 
c               atomic configuration. In Ferromagnetic ordering (up).
c      If DM.InitSpinAF is true, as default but in Antiferro order:
c               even atoms have spin down, odd up.
c      If fdf %block DM.InitSpin is present it overwrites previous 
c         schemes: magnetic moments are explicitly given for some atoms.
c         Atoms not mentioned in the block are initialized non polarized.
c      
c Written by E. Artacho. December 1997. Taken from the original piece
c of siesta.f written by P. Ordejon.
c Non-collinear spin added by J.M.Soler, May 1998.
c ********* INPUT ***************************************************
c logical try_dm_from_file     : whether DM has to be read from files or not
c logical found         : whether DM was found in files
c logical inspn         : true : AF ordering according to atom ordering
c                                if no DM files, no DM.InitSpin, ispin=2
c                         false: Ferro ordering  (fdf DM.InitSpinAF)
c integer nua           : Number of atoms in the unit cell
c integer no_l           : Number of orbitals in the unit cell
c integer nspin         : Number of spin components
c integer maxa          : Max num. atoms for dimension
c integer maxo          : Max. number of orbitals (globally)
c integer maxuo         : Max. number of orbitals (locally)
c integer maxnh         : Max number of nonzero interactions
c integer maxspn        : Max number of spin components in dimensions
c integer lasto(0:maxa) : List with last orbital of each atom
c integer numh(*)       : Dscf matrix sparse information
c integer listhptr(*)   :  "
c integer listh(maxnh)  :  "
c integer iaorb(maxo)   : List saying to what atom an orbital belongs 
c double Datm(no)       : Occupations of basis orbitals in free atom
c ********* OUTPUT **************************************************
c double Dscf(maxnh,maxspn) : Density matrix in sparse form
c *******************************************************************

C
C  Modules
C
      use precision
      use parallel,     only : Node, Nodes, IOnode
      use parallelsubs, only : LocalToGlobalOrb, GlobalToLocalOrb
      use fdf
      use parsing
      use sys,          only : die
      use m_iodm,       only : read_dynamic_dm
      use m_sparse,     only : same_sparsity, change_sparsity
      use alloc,        only : re_alloc, de_alloc
#ifdef MPI
      use mpi_siesta
#endif
      use units, only : pi

      implicit          none

      logical           found, inspn, try_dm_from_file
      integer           no_l, nua, maxo, maxuo, maxnh, nspin, maxa, 
     .                  maxspn
      integer           lasto(0:maxa), numh(maxuo), iaorb(maxo)
      integer, pointer ::  listhptr(:), listh(:)
      real(dp), pointer :: Dscf(:,:)
      real(dp)          Datm(maxo)
      external          memory

c ---------------------------------------------------------------------

C Internal variables and arrays
 
      character         updo*1, msg*80
      logical           noncol, peratm, badsyntax
      integer           nh, ni, nn, nr, nv, iat, nat, ia, iu,
     .                  i1, i2, in, ind, ispin, jo, io, 
     .                  iio, maxatnew

      integer, save ::  maxat

      integer           integs(4), lastc, lc(0:3)

      integer, pointer, save ::  atom(:)
#ifdef MPI
      integer  MPIerror
      logical  lbuffer
#endif
      real(dp)          aspin, cosph, costh, epsilon, 
     .                  qio, rate, reals(4),
     .                  sinph, sinth, spinat, spio, values(4)

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

      real(dp), pointer, save :: phi(:), spin(:), theta(:)

      integer, pointer :: list_tmp(:), num_tmp(:), listptr_tmp(:)
      real(dp), pointer :: D_tmp(:,:), work(:)

      integer :: nh_tmp, nspin_read, is
      logical :: same_sp

      data maxat / 1000 /
      data epsilon / 1.d-8 /

C Try to read DM from disk if wanted (DM.UseSaveDM true) ---------------

      if (try_dm_from_file) then

         nullify(list_tmp,D_tmp,num_tmp,listptr_tmp)
         call re_alloc(num_tmp,1,no_l)
         call re_alloc(listptr_tmp,1,no_l)

         call read_dynamic_dm(nh_tmp, no_l, nspin_read, num_tmp, 
     .             listptr_tmp, list_tmp, D_tmp, found )
      else
        found = .false.
      endif

C If found, check and update, otherwise initialize with neutral atoms

      if (found) then
        ! Various degrees of sanity checks
        if (nspin_read /= nspin) then
           if (Node.eq.0) then
              write(6,"(a,i6,/,a)")
     $        "WARNING: Wrong nspin in DM file: ",
     $             nspin_read,
     $        "WARNING: Falling back to atomic initialization of DM."
              write(0,"(a,i6,/,a)")
     $        "WARNING: Wrong nspin in DM file: ",
     $             nspin_read,
     $        "WARNING: Falling back to atomic initialization of DM."
           endif
           found = .false.
        endif
      endif

      if (found) then

!         anorm = compute_norm(D_tmp,nh_tmp,nspin_read)
!         if (IONode) print *, "DM norm of file contents: ", anorm

         same_sp = same_sparsity(no_l,num_tmp,list_tmp,
     $                                 numh,listh)
#ifdef MPI
        call MPI_AllReduce(same_sp,lbuffer,1,MPI_logical,MPI_land,
     .    MPI_Comm_World,MPIerror)
        same_sp = lbuffer
#endif
        if (same_sp) then
           Dscf(:,:) = D_tmp(:,:)
        else
           if (IOnode) then
              write(6,*) "Read DM has different structure. Fixing..."
           endif
           nullify(work)
           do is= 1, nspin
              call change_sparsity (no_l,num_tmp,listptr_tmp,list_tmp,
     $                                  numh,listhptr,listh,
     $                                  D_tmp(:,is), Dscf(:,is),
     $                            work)          
           enddo
           call de_alloc(num_tmp)
           call de_alloc(listptr_tmp)
           call de_alloc(list_tmp)
           call de_alloc(D_tmp)
           call de_alloc(work)  ! Allocated by change_sparsity
        endif


      else

C See whether specific initial spins are given in a DM.InitSpin block
C and read them in a loop on atoms where lines are read and parsed
C   integer nat       : how many atoms to polarize
C   integer atom(nat) : which atoms
C   double  spin(nat) : what polarization -----------------------------

        noncol = .false.
        peratm = fdf_block('DM.InitSpin',bfdf)
        noncol = .false.
        if (Node.eq.0) then
          if (peratm .and. nspin.lt.2) write(6,'(/,a)')
     .    'initdm: WARNING: DM.InitSpin not used because nspin < 2'
        endif

        if (peratm .and. nspin.ge.2) then

C Allocate local memory
          nullify(atom,phi,spin,theta)
          call re_alloc( atom, 1, maxat, 'atom', 'initdm' )
          call re_alloc( phi, 1, maxat, 'phi', 'initdm' )
          call re_alloc( spin, 1, maxat, 'spin', 'initdm' )
          call re_alloc( theta, 1, maxat, 'theta', 'initdm' )

          nat = 0
          badsyntax = .FALSE.
          do while(fdf_bline(bfdf,pline) .and. (nat .lt. nua) .and.
     .            (.not. badsyntax))

            nn = fdf_bnnames(pline)
            ni = fdf_bnintegers(pline)
            nr = fdf_bnreals(pline)

            if (ni .eq. 1) then
              if (nat .eq. maxat) then
                maxatnew = nat + nint(0.1*nat)
C
                call re_alloc( atom, 1, maxatnew, 'atom', 'initdm',
     &                         copy=.true. )
C
                call re_alloc( phi, 1, maxatnew, 'phi', 'initdm',
     &                         copy=.true. )
                call re_alloc( spin, 1, maxatnew, 'spin', 'initdm',
     &                         copy=.true. )
                call re_alloc( theta, 1, maxatnew, 'theta', 'initdm',
     &                         copy=.true. )
C
                maxat = maxatnew
              endif
              nat = nat + 1
              atom(nat) = fdf_bintegers(pline,1)

              if (nn .eq. 0) then
C Read value of spin
                if (nr .eq. 3) then
C Read spin value and direction
                  spin(nat)  = fdf_breals(pline,1)
                  theta(nat) = fdf_breals(pline,2) * pi/180.0d0
                  phi(nat)   = fdf_breals(pline,3) * pi/180.0d0
                elseif (nr .eq. 1) then
C Read spin value. Default direction.
                  spin(nat)  = fdf_breals(pline,1)
                  theta(nat) = 0.d0
                  phi(nat)   = 0.d0
                else
C Print bad-syntax error and stop
                  badsyntax = .TRUE.
                endif
              elseif (nn .eq. 1) then
C Read spin as + or - (maximun value)
                updo = fdf_bnames(pline,1)
                if (updo .eq. '+') then
                  spin(nat) =  100.d0
                elseif (updo .eq. '-') then
                  spin(nat) = -100.d0
                else
C Print bad-syntax error and stop
                  badsyntax = .TRUE.
                endif
                if (nr .eq. 2) then
                  theta(nat) = fdf_breals(pline,1) * pi/180.0d0
                  phi(nat)   = fdf_breals(pline,2) * pi/180.0d0
                elseif (nr .eq. 0) then
                  theta(nat) = 0.d0
                  phi(nat)   = 0.d0
                else
C Print bad-syntax error and stop
                  badsyntax = .TRUE.
                endif
              else
C Print bad-syntax error and stop
                badsyntax = .TRUE.
              endif

              if ((atom(nat) .lt. 1) .or. (atom(nat) .gt. nua)) then
                write(msg,'(a,a,i4)') 'intdm: ERROR: Bad atom ' //
     .            'index in DM.InitSpin, line', nat+1
                call die(TRIM(msg))
              endif
              if (abs(theta(nat)) .gt. 1.d-12) noncol = .true.
            else
C Print bad-syntax error and stop
              badsyntax = .TRUE.
            endif
          enddo

          if (badsyntax) then
            write(msg,'(a,i4)')
     .        'initdm: ERROR: bad syntax in DM.InitSpin, line', nat+1
            call die(msg)
          endif

          if (noncol .and. nspin.lt.4) then
            if (Node.eq.0) then
            write(6,'(/,2a)') 'initdm: WARNING: noncolinear spins ',
     .                 'in DM.InitSpin not used because nspin < 4'
            endif
            noncol = .false.
          endif

C Initialize to 0

          Dscf(1:maxnh,1:nspin) = 0.0d0

C Initialize all paramagnetic 

          do ia = 1, nua
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

C Loop on atoms with spin

          do iat = 1, nat
            ia = atom(iat)

C Find maximum atomic moment that the atoms involved can carry
          
            spinat = 0.d0
            do io = lasto(ia-1) + 1, lasto(ia)
              spinat = spinat + min( Datm(io), 2.d0 - Datm(io) )
            enddo
            if (spinat.lt.epsilon .and. Node.eq.0) print'(a,i6,a)', 
     .        'initdm: WARNING: atom ', atom(iat),
     .        ' has a closed-shell and cannot be polarized'

C If given spin is larger than possible, make it to max atomic

            aspin = abs(spin(iat))
            if ((aspin .gt. spinat) .and. (aspin .gt. epsilon)) 
     .         spin(iat) = spinat*spin(iat)/aspin 

C Initialize orbitals with same rate as atom

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
C Store non-collinear-spin density matrix as
C   ispin=1 => D11, ispin=2 => D22;
C   ispin=3 => Real(D12); ispin=4 => Imag(D12)
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

C Deallocate local memory
          call de_alloc( atom, 'atom', 'initdm' )
          call de_alloc( phi, 'phi', 'initdm' )
          call de_alloc( spin, 'spin', 'initdm' )
          call de_alloc( theta, 'theta', 'initdm' )

c ---------------------------------------------------------------------

        else

C Initialize to 0
          Dscf(1:maxnh,1:nspin) = 0.0d0

C Automatic, for non magnetic (nspin=1) or for Ferro or Antiferro -----
          do io = 1, no_l
            call LocalToGlobalOrb(io,Node,Nodes,iio)
            do in = 1,numh(io)
              ind = listhptr(io)+in
              jo = listh(ind)
              if (iio .eq. jo) then
                if (nspin .eq. 1) then

C No spin polarization

                  Dscf(ind,1) = Datm(iio)
                else

C Spin polarization

                  i1 = 1
                  i2 = 2

C Ferro or antiferro according to DM.InitSpinAF (inspn)

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

      endif

!      anorm =  compute_norm(Dscf,maxnh,nspin)
!      if (IONode) print *, "DM norm at end of initdm: ", anorm

      end subroutine initdm

      END MODULE m_new_dm




