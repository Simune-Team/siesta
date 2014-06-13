module m_ts_electype

  use precision, only : dp

  use class_Sparsity
  use class_dSpData1D
  use class_dSpData2D

  use m_ts_chem_pot, only : ts_mu, name

  implicit none

  private

  ! Name is part of the class-system, we need
  ! it as an interface
  interface name
     module procedure name_
  end interface name

  interface delete
     module procedure delete_
  end interface delete

  interface operator(.eq.)
     module procedure equal_el_el
     module procedure equal_el_str
     module procedure equal_str_el
  end interface

  public :: Elec
  public :: Name, HSfile, GFFile, GFTitle
  public :: TotUsedAtoms, TotUsedOrbs
  public :: AtomInElec, OrbInElec
  public :: Rep
  public :: q_exp

  public :: fdf_nElec, fdf_elec

  public :: assign, read_Elec
  public :: create_sp2sp01
  public :: print_Elec
  public :: check_Elec, check_connectivity
  public :: delete

  public :: operator(.eq.)

  ! 300 chars for a full path should be fine
  integer, parameter, public :: FILE_LEN = 300
  integer, parameter, public :: NAME_LEN = 100

  integer, parameter, public :: INF_NEGATIVE = 0 ! old 'left'
  integer, parameter, public :: INF_POSITIVE = 1 ! old 'right'

  type :: Elec
     character(len=FILE_LEN) :: HSfile = ' ', GFfile  = ' '
     character(len=NAME_LEN) :: Name   = ' ', GFtitle = ' '
     ! These variables are relative to the big system
     integer :: idx_na = 0, idx_no = 0
     ! atoms used
     integer :: na_used = 0
     ! orbitals used
     integer :: no_used = 0
     ! repetitions
     integer :: RepA1 = 1, RepA2 = 1, RepA3 = 1
     ! chemical potential of the electrode
     type(ts_mu), pointer :: mu => null()
     ! infinity direction
     integer :: inf_dir = INF_NEGATIVE
     ! transport direction (determines H01)
     integer :: t_dir = 3 
     ! whether the electrode should be bulk
     logical :: Bulk = .true.
     logical :: DM_CrossTerms = .false.
     logical :: BandBottom = .false.
     ! whether to re-calculate the GF-file
     logical :: ReUseGF = .false. 
     ! Create a GF-file or re-calculate the self-energies everytime
     logical :: out_of_core = .true. 
     ! In case of 'out_of_core == .false.' we can reduce the number of operations
     ! by skipping the copying of H00, S00. Hence we need to compare when the 
     ! k-point changes...
     real(dp), pointer :: bkpt_cur(:) => null()
     ! If the user requests to assign different "spill-in fermi-level" 
     ! we allow that
     real(dp) :: Ef_frac_CT = 0._dp
     ! Used xa and lasto
     real(dp), pointer :: xa_used(:,:) => null()
     integer,  pointer :: lasto_used(:) => null()

     ! Advanced
     logical :: kcell_check = .true.

     ! ---v--- Below we have the content of the TSHS file
     integer  :: nspin = 0, na_u = 0, no_u = 0, no_s = 0
     real(dp) :: ucell(3,3) = 0._dp, Ef = 0._dp, Qtot = 0._dp
     real(dp), pointer :: xa(:,:) => null()
     integer,  pointer :: lasto(:) => null()
     type(Sparsity)  :: sp
     type(dSpData2D) :: H, xij
     type(dSpData1D) :: S
     ! --- --- completed the content of the TSHS file
     ! Below we create the content for the self-energy creation
     ! Notice that we can save some elements simply by extracting the 0-1 connections
     ! for large systems this is a non-negligeble part of the memory...
     type(Sparsity)  :: sp00
     type(dSpData2D) :: H00, xij00, xijo00
     type(dSpData1D) :: S00
     type(Sparsity)  :: sp01
     type(dSpData2D) :: H01, xij01, xijo01
     type(dSpData1D) :: S01

     ! These arrays are used to construct the full Hamiltonian and overlap and Green's function
     complex(dp), pointer :: HA(:,:,:), SA(:,:,:), GA(:,:)
     ! Arrays needed to partition the scattering matrix and self-energies
     ! Notice that Gamma should "ALWAYS" contain the transposed
     complex(dp), pointer :: Gamma(:,:), Sigma(:)

  end type Elec

contains

  function fdf_nElec(prefix,this_n) result(n)
    use fdf

    character(len=*), intent(in) :: prefix
    type(Elec), allocatable :: this_n(:)
    integer :: n

    ! prepare to read in the data...
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline => null()
    integer :: i
    
    logical :: found

    n = 0
    found = fdf_block(trim(prefix)//'.Elecs',bfdf)
    if ( .not. found ) return

    ! first count the number of electrodes
    n = 0
    do while ( fdf_bline(bfdf,pline) )
       if ( fdf_bnnames(pline) == 0 ) cycle
       n = n + 1 
    end do

    allocate(this_n(n))

    ! rewind to read again
    call fdf_brewind(bfdf)

    n = 0
    do while ( fdf_bline(bfdf,pline) )
       if ( fdf_bnnames(pline) == 0 ) cycle
       n = n + 1 
       this_n(n)%Name = trim(fdf_bnames(pline,1))
       if ( n > 1 ) then
          ! Check that no name is the same
          do i = 1 , n - 1 
             if ( leqi(this_n(i)%name,this_n(n)%name) ) then
                call die('Electrode names must not be the same')
             end if
          end do
       end if
    end do

  end function fdf_nElec

  function fdf_Elec(prefix,slabel,this,N_mu,mus) result(found)
    use fdf
    use m_ts_io, only : ts_read_TSHS_opt

    character(len=*), intent(in) :: prefix,slabel
    type(Elec), intent(inout) :: this
    integer, intent(in) :: N_mu
    type(ts_mu), intent(in), target :: mus(N_mu)
    logical :: found

    ! prepare to read in the data...
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline => null()
    logical :: info(4)
    integer :: i, j
    integer :: idx_na 

    character(len=200) :: name, ln, tmp

    name = trim(this%name)
    found = fdf_block(trim(prefix)//'.Elec.'//trim(name),bfdf)
    if ( .not. found ) return

    info(:) = .false.
    idx_na = 0

    ! We default a lot of the options
    this%GFtitle = 'Surface-Greens function for '//trim(name)
    this%GFfile = trim(slabel)//'.'//trim(prefix)//'GF'//trim(name)
    this%na_used = -1
    
    do while ( fdf_bline(bfdf,pline) )
       if ( fdf_bnnames(pline) == 0 ) cycle
       
       ln = fdf_bnames(pline,1) 
       
       ! We select the input
       if ( leqi(ln,'TSHS') .or. &
            leqi(ln,'TSHS-file') ) then
          if ( fdf_bnnames(pline) < 2 ) call die('TSHS name not supplied')
          this%HSfile = trim(fdf_bnames(pline,2))
          info(1) = .true.

       else if ( leqi(ln,'semi-inf-direction') .or. &
            leqi(ln,'semi-inf-dir') .or. leqi(ln,'semi-inf') ) then

          tmp = 'Semi-infinite direction not understood correctly, &
                  &allowed format: [-+][a-c|a[1-3]]'

          ! This possibility exists
          !  semi-inf [-+][ ][a-c|a[1-3]] -> [direction] [vector]
          if ( fdf_bnnames(pline) < 2 ) then
             call die(trim(tmp))
          end if
          
          ln = fdf_bnames(pline,2)
          if ( fdf_bnnames(pline) > 2 ) then
             if ( len_trim(ln) /= 1 ) then
                call die(trim(tmp))
             end if

             ln = trim(ln) // fdf_bnames(pline,3)

          end if
          
          ! now for testing
          if ( ln(1:1) == '+' ) then
             this%inf_dir = INF_POSITIVE
          else if ( ln(1:1) == '-' ) then
             this%inf_dir = INF_NEGATIVE
          else
             call die(trim(tmp))
          end if

          ! copy over remaining part...
          ln = ln(2:)
          if ( leqi(ln,'a') .or. leqi(ln,'a1') ) then
             this%t_dir = 1
          else if ( leqi(ln,'b') .or. leqi(ln,'a2') ) then
             this%t_dir = 2
          else if ( leqi(ln,'c') .or. leqi(ln,'a3') ) then
             this%t_dir = 3
          else
             call die(trim(tmp))
          end if
          info(2) = .true.
          
       else if ( leqi(ln,'chemical-potential') .or. &
            leqi(ln,'chem-pot') .or. leqi(ln,'mu') ) then
          if ( fdf_bnnames(pline) < 2 ) &
               call die('Name of chemical-potential not supplied')
          ln = fdf_bnames(pline,2)
          nullify(this%mu)
          do i = 1 , N_mu
             if ( leqi(ln,mus(i)%name) ) then
                this%mu => mus(i)
                exit
             end if
          end do
          if ( .not. associated(this%mu) ) then
             call die('Could not find the chemical potential "'//trim(ln)//&
                  '" for the electrode: "'//trim(name)//'". '//&
                  'Please supply an existing name.')
          end if
          info(3) = .true.
          
       else if ( leqi(ln,'electrode-position') .or. &
            leqi(ln,'elec-pos') ) then
          idx_na      = 0
          this%idx_na = 0
          if ( fdf_bnnames(pline) > 1 ) then
             ! the user is requesting on a string basis
             ln = fdf_bnames(pline,2)
             if ( leqi(ln,'start') .or. leqi(ln,'begin') ) then
                if ( fdf_bnintegers(pline) > 0 ) then
                   this%idx_na = fdf_bintegers(pline,1)
                else
                   this%idx_na = 1 ! default starting position
                end if
             else if ( leqi(ln,'end') ) then
                idx_na      = -1
                this%idx_na =  0
                if ( fdf_bnintegers(pline) > 0 ) then
                   this%idx_na = fdf_bintegers(pline,1)
                end if
             end if
          else
             if ( fdf_bnintegers(pline) < 1 ) &
                  call die('Atomic position not found in input line: '//trim(ln))
             this%idx_na = fdf_bintegers(pline,1)
          end if
          info(4) = .true.
          
       else if ( leqi(ln,'bulk') ) then
          this%Bulk = fdf_bboolean(pline,1,after=1)

       else if ( leqi(ln,'calculate-band-bottom') ) then
          this%BandBottom = fdf_bboolean(pline,1,after=1)

       else if ( leqi(ln,'update-cross-terms') ) then
          this%DM_CrossTerms = fdf_bboolean(pline,1,after=1)

       else if ( leqi(ln,'Ef-fraction') ) then

          ! highly experimental feature,
          ! it determines the fraction of the electrode fermi-level
          ! that shifts the energy, of the H_{E-C} region.
          ! instead of the energy at Ef
          if ( fdf_bnvalues(pline) < 1 ) &
               call die('Fraction specification missing.')
          
          this%Ef_frac_CT = fdf_bvalues(pline,1,after=1)
          if ( this%Ef_frac_CT < 0._dp .or. &
               1._dp < this%Ef_frac_CT ) then
             call die('Fraction for fermi-level must be in [0;1] range')
          end if

       else if ( leqi(ln,'GF-title') ) then
          if ( fdf_bnnames(pline) < 2 ) &
               call die('GF-title not supplied')
          this%GFtitle = trim(fdf_bnames(pline,2))

       else if ( leqi(ln,'GF') .or. &
            leqi(ln,'GF-file') ) then
          if ( fdf_bnnames(pline) < 2 ) call die('GF-file not supplied')
          this%GFfile = trim(fdf_bnames(pline,2))

       else if ( leqi(ln,'GF-ReUse') ) then

          this%ReUseGF = fdf_bboolean(pline,1,after=1)

       else if ( leqi(ln,'used-atoms') ) then
          if ( fdf_bnintegers(pline) < 1 ) &
               call die('Number of atoms used not supplied')
          this%na_used = fdf_bintegers(pline,1)

       else if ( leqi(ln,'replicate-a') .or. leqi(ln,'rep-a') .or. &
            leqi(ln,'replicate-a1') .or. leqi(ln,'rep-a1') ) then
          if ( fdf_bnintegers(pline) < 1 ) &
               call die('Repetition A1 is not supplied')
          this%RepA1 = fdf_bintegers(pline,1)

       else if ( leqi(ln,'replicate-b') .or. leqi(ln,'rep-b') .or. &
            leqi(ln,'replicate-a2') .or. leqi(ln,'rep-a2') ) then
          if ( fdf_bnintegers(pline) < 1 ) &
               call die('Repetition A2 is not supplied')
          this%RepA2 = fdf_bintegers(pline,1)

       else if ( leqi(ln,'replicate-c') .or. leqi(ln,'rep-c') .or. &
            leqi(ln,'replicate-a3') .or. leqi(ln,'rep-a3') ) then
          if ( fdf_bnintegers(pline) < 1 ) &
               call die('Repetition A3 is not supplied')
          this%RepA3 = fdf_bintegers(pline,1)

       else if ( leqi(ln,'replicate') .or. leqi(ln,'rep') ) then
          if ( fdf_bnintegers(pline) < 3 ) &
               call die('Repetition for all directions are not supplied <A1> <A2> <A3>')
          this%RepA1 = fdf_bintegers(pline,1)
          this%RepA2 = fdf_bintegers(pline,2)
          this%RepA3 = fdf_bintegers(pline,3)

       else if ( leqi(ln,'out-of-core') ) then

          this%out_of_core = fdf_bboolean(pline,1,after=1)

       else if ( leqi(ln,'check-kgrid') ) then

          ! Advanced (NOT recommended)
          ! if false, it will not check the k-point sampling
          ! in the electrode
          ! I.e. it will expect the system to be converged as
          ! the system
          ! However, it only checks the k-cell
          ! and still suspects a non-Gamma calculation
          this%kcell_check = fdf_bboolean(pline,1,after=1)

       else

          ! we should always die in case something non-understandable 
          ! is passed, if that is the case, the chances are that the
          ! user has made a typo is high
          call die('Unrecognized option "'//trim(ln)//'" &
               &for electrode: '//trim(name))

       end if

    end do
    
    if ( this%RepA1 < 1 .or. this%RepA2 < 1 .or. this%RepA3 < 1 ) &
         call die("Repetition in "//trim(name)//" electrode must be >= 1.")

    if ( .not. all(info) ) then
       write(*,*)'You need to supply at least:'
       write(*,*)' - TSHS'
       write(*,*)' - semi-inf-direction'
       write(*,*)' - chemical-potential'
       write(*,*)' - electrode-position'
       call die('You have not supplied all electrode information')
    end if

    inquire(file=trim(HSFile(this)), exist=info(1))
    if ( .not. info(1) ) then
       call die("Electrode file does not exist. &
            &Please create electrode '"//trim(HSFile(this))//"' first.")
    end if

    ! Read in the number of atoms in the HSfile
    call ts_read_TSHS_opt(HSFile(this),no_u=this%no_u,na_u=this%na_u, &
         nspin=this%nspin, Ef=this%Ef, ucell=this%ucell, &
         Bcast=.true.)

    allocate(this%xa(3,this%na_u),this%lasto(0:this%na_u))
    call ts_read_TSHS_opt(HSFile(this),xa=this%xa,lasto=this%lasto, &
         Bcast=.true.)

    ! in case the number of used atoms has not been set
    if ( this%na_used <= 0 ) this%na_used = this%na_u

    if ( this%na_used <= 0 ) &
         call die("None atoms requested for electrode calculation.")

    if ( this%na_u < this%na_used ) then
       write(*,*) "# of requested atoms is larger than available."
       write(*,*) "Requested: ",this%na_used
       write(*,*) "Available: ",this%na_u
       call die("Error on requested atoms, please correct input.")
    end if

    allocate(this%lasto_used(0:this%na_used),this%xa_used(3,this%na_used))
    this%lasto_used(0) = 0
    this%no_used = 0
    if ( this%inf_dir == INF_NEGATIVE ) then ! same as old 'left'
       ! We use the last atoms
       j = 0
       do i = this%na_u - this%na_used + 1 , this%na_u
          j = j + 1
          this%lasto_used(j) = this%lasto_used(j-1) + this%lasto(i)-this%lasto(i-1)
          this%xa_used(:,j)  = this%xa(:,i)
       end do

    else if ( this%inf_dir == INF_POSITIVE ) then ! same as old 'right'
       ! We use the first atoms
       do i = 1 , this%na_used
          this%lasto_used(i) = this%lasto_used(i-1) + this%lasto(i)-this%lasto(i-1)
          this%xa_used(:,i)  = this%xa(:,i)
       end do

    else
       call die('Unknown direction for the semi-infinite lead')

    end if
    this%no_used = this%lasto_used(this%na_used)

    ! We deallocate xa and lasto as they are not needed
    deallocate(this%xa,this%lasto)

    ! Check that the repetition is not in the transport-direction
    select case ( this%t_dir )
    case ( 1 )
       if ( this%RepA1 /= 1 ) then
          call die('Repetition in the transport direction &
               &is not allowed.')
       end if
    case ( 2 )
       if ( this%RepA2 /= 1 ) then
          call die('Repetition in the transport direction &
               &is not allowed.')
       end if
    case ( 3 )
       if ( this%RepA3 /= 1 ) then
          call die('Repetition in the transport direction &
               &is not allowed.')
       end if
    end select

    ! if the electrode does not use a bulk electrode we need to update
    ! the cross-terms
    if ( .not. this%Bulk ) then
       this%DM_CrossTerms = .true.
    end if

    if ( this%DM_CrossTerms .or. this%mu%mu == 0._dp ) then
       ! when updating the cross-terms
       ! there is no needing to also shift the energy
       ! I think this will battle each other out...
       this%Ef_frac_CT = 0._dp
    end if

    ! if the user has specified text for the electrode position
    if ( idx_na == -1 ) then
       this%idx_na = this%idx_na + 1 - TotUsedAtoms(this)
    end if

  end function fdf_Elec

  function equal_el_el(this1,this2) result(equal)
    type(Elec), intent(in) :: this1, this2
    logical :: equal
    equal = this1%name == this2%name
  end function equal_el_el

  function equal_el_str(this,str) result(equal)
    type(Elec), intent(in) :: this
    character(len=*), intent(in) :: str
    logical :: equal
    equal = this%name == trim(str)
  end function equal_el_str

  function equal_str_el(str,this) result(equal)
    character(len=*), intent(in) :: str
    type(Elec), intent(in) :: this
    logical :: equal
    equal = this%name == trim(str)
  end function equal_str_el

  subroutine assign(this,D,HSfile,GFfile,GFtitle, &
       na_u,na_used,no_u,no_s,no_used, &
       RepA1,RepA2,RepA3)
    type(Elec), intent(inout) :: this
    character, optional, intent(in) :: D
    character(len=*), intent(in), optional :: HSfile, GFfile, Gftitle
    integer, intent(in), optional :: na_u, na_used, no_u, no_s, no_used
    integer, intent(in), optional :: RepA1, RepA2, RepA3

    if (present(D)) call die('Wrong usage of Elec-type assign')

    if (present(HSfile)) this%HSfile = HSfile
    if (present(GFfile)) this%GFfile = GFfile
    if (present(GFtitle)) this%GFtitle = GFtitle
    if (present(na_u)) this%na_u = na_u
    if (present(na_used)) this%na_used = na_used
    if (present(no_u)) this%no_u = no_u
    if (present(no_s)) this%no_s = no_s
    if (present(no_used)) this%no_used = no_used
    if (present(RepA1)) this%RepA1 = RepA1
    if (present(RepA2)) this%RepA2 = RepA2
    if (present(RepA3)) this%RepA3 = RepA3

  end subroutine assign
  
  elemental function Name_(this)
    type(Elec), intent(in) :: this
    character(len=NAME_LEN) :: Name_
    Name_ = this%Name
  end function Name_

  elemental function HSfile(this)
    type(Elec), intent(in) :: this
    character(len=FILE_LEN) :: HSfile
    HSfile = this%HSfile
  end function HSfile

  elemental function GFfile(this)
    type(Elec), intent(in) :: this
    character(len=FILE_LEN) :: GFfile
    GFfile = this%GFfile
  end function GFfile

  elemental function GFtitle(this)
    type(Elec), intent(in) :: this
    character(len=NAME_LEN) :: GFtitle
    GFtitle = this%GFtitle
  end function GFtitle

  elemental function TotUsedAtoms(this) result(val)
    type(Elec), intent(in) :: this
    integer :: val
    val = this%na_used * Rep(this)
  end function TotUsedAtoms

  elemental function Rep(this) result(val)
    type(Elec), intent(in) :: this
    integer :: val
    val = this%RepA1 * this%RepA2 * this%RepA3
  end function Rep

  pure function q_exp_all(this,i,j,k) result(q)
    type(Elec), intent(in) :: this
    integer, intent(in) :: i,j,k
    real(dp) :: q(3)
    q(1) = 1._dp*(i-1) / real(this%RepA1,dp)
    q(2) = 1._dp*(j-1) / real(this%RepA2,dp)
    q(3) = 1._dp*(k-1) / real(this%RepA3,dp)
  end function q_exp_all

  pure function q_exp(this,idx) result(q)
    type(Elec), intent(in) :: this
    integer, intent(in) :: idx
    real(dp) :: q(3)
    integer :: i,j,k,ii
    i =     this%RepA1
    j = i * this%RepA2
    k = j * this%RepA3
    if ( idx <= i ) then
       q = q_exp_all(this,idx,1,1)
    else if ( idx <= j ) then
       j = idx / i
       if ( MOD(idx,i) /= 0 ) j = j + 1
       i = idx - (j-1) * i
       q = q_exp_all(this,i,j,1)
    else if ( idx <= k ) then
       ! this seems to work well!
       k = idx / j
       if ( MOD(idx,j) /= 0 ) k = k + 1
       ii = idx - (k-1) * j
       j = ii / i
       if ( MOD(ii,i) /= 0 ) j = j + 1
       i = ii - (j-1) * i
       q = q_exp_all(this,i,j,k)
    else
       q = 0._dp
    end if
  end function q_exp

  elemental function TotUsedOrbs(this) result(val)
    type(Elec), intent(in) :: this
    integer :: val
    val = this%no_used * Rep(this)
  end function TotUsedOrbs

  elemental function OrbInElec(this,io) result(in)
    type(Elec), intent(in) :: this
    integer, intent(in) :: io
    logical :: in
    in = this%idx_no <= io .and. io < (this%idx_no + TotUsedOrbs(this))
  end function OrbInElec

  elemental function AtomInElec(this,ia) result(in)
    type(Elec), intent(in) :: this
    integer, intent(in) :: ia
    logical :: in
    in = this%idx_na <= ia .and. ia < (this%idx_na + TotUsedAtoms(this))
  end function AtomInElec

  subroutine read_Elec(this,Bcast,io)
    use fdf
    use parallel
    use class_OrbitalDistribution

    use m_ts_io
#ifdef MPI
    use mpi_siesta
#endif

    type(Elec), intent(inout) :: this
    logical, intent(in), optional :: Bcast
    logical, intent(in), optional :: io

    character(len=200) :: fN
    integer :: fL, kscell(3,3), istep, ia1
    logical :: onlyS, Gamma_file, TSGamma, lio
    real(dp) :: temp, kdispl(3)
    
    ! Sparsity pattern
    integer, pointer :: iza(:), numh(:), listhptr(:), listh(:), indxuo(:)
    real(dp), pointer :: H(:,:), S(:), xij(:,:), t2D(:,:), t1D(:)
    integer :: n_nzs

    type(OrbitalDistribution) :: fdist

    lio = .true.
    if ( present(io) ) lio = io

    fN = trim(HSfile(this))
    ! We read in the information
    fL = len_trim(fN)
    if ( leqi(fN(fL-4:fL),'.TSHS') ) then
       call ts_read_tshs(fN, &
            onlyS, Gamma_file, TSGamma, &
            this%ucell, this%na_u, this%no_u, &
            this%no_u, this%no_s, n_nzs, this%nspin,  &
            kscell, kdispl, &
            this%xa, iza, this%lasto, &
            numh, listhptr, listh, xij, indxuo, &
            H, S, this%Ef, &
            this%Qtot, Temp, & ! Qtot, Temp
            istep, ia1, &
            Bcast=Bcast)
    else
       call die('Could not infer the file type of the &
            &electrode file: '//trim(fN))
    end if
    deallocate(indxuo,iza)
    call memory('D','I',this%na_u+this%no_s,'iohs')

    call newSparsity(this%sp,this%no_u,this%no_u, &
         n_nzs, numh, listhptr, listh, "Electrode "//trim(fN))

    ! Create the Fake distribution
    ! The Block-size is the number of orbitals, i.e. all on the first processor
    ! Notice that we DO need it to be the SIESTA size.
#ifdef MPI
    call newDistribution(this%no_u,MPI_Comm_Self,fdist,name='TS-fake dist')
#else
    call newDistribution(this%no_u,-1           ,fdist,name='TS-fake dist')
#endif

    ! Create containers
    call newdSpData2D(this%sp,this%nspin,fdist,this%H,name='E spH')
    call newdSpData1D(this%sp,fdist,this%S,name='E spS')

    ! Copy data
    t2D => val(this%H)
    t2D(:,:) = H(:,:)
    t1D => val(this%S)
    t1D(:) = S(:)

    ! clean-up
    deallocate(H,S)
    call memory('D','D',n_nzs*(this%nspin+1),'iohs')

    ! copy xij
    call newdSpData2D(this%sp,3,fdist,this%xij,name='E spxij', &
         sparsity_dim=2)
    t2D => val(this%xij)
    t2D = xij
    deallocate(xij)
    call memory('D','D',n_nzs*3,'iohs')

    deallocate(numh,listhptr,listh)
    call memory('D','I',this%no_u*2+n_nzs,'iohs')

    if ( IONode .and. lio ) call print_type(this%sp)

  end subroutine read_Elec

  subroutine create_sp2sp01(this,calc_xijo,IO)

    use parallel, only : IONode, Node

    use class_OrbitalDistribution

    use create_Sparsity_SC
    use geom_helper, only : iaorb,cell_c
#ifdef MPI
    use mpi_siesta
#endif

    type(Elec), intent(inout) :: this
    logical, intent(in), optional :: calc_xijo
    logical, intent(in), optional :: io

    logical :: lcalc_xijo, lio

    real(dp), pointer :: xij(:,:), xij00(:,:), xij01(:,:)
    real(dp), pointer :: xijo00(:,:), xijo01(:,:)
    real(dp), pointer :: H(:,:), H00(:,:), H01(:,:)
    real(dp), pointer :: S(:), S00(:), S01(:)
    type(OrbitalDistribution), pointer :: fdist

    integer, pointer :: l_ncol(:) => null()
    integer, pointer :: l_ptr(:)  => null()
    integer, pointer :: l_col(:)  => null()
    integer, pointer :: ncol00(:) => null()
    integer, pointer :: ptr00(:)  => null()
    integer, pointer :: col00(:)  => null()
    integer, pointer :: ncol01(:) => null()
    integer, pointer :: ptr01(:)  => null()
    integer, pointer :: col01(:)  => null()

    integer :: no_l, i, iio, j, ind, ind00, ind01, ia, ja
    integer :: tm(3)

    ! Retrieve distribution
    fdist => dist(this%H)

    lcalc_xijo = .false.
    if ( present(calc_xijo) ) lcalc_xijo = calc_xijo
    lio = .true.
    if ( present(io) ) lio = io

    H   => val(this%H)
    S   => val(this%S)
    xij => val(this%xij)
    tm(:)          = TM_ALL
    tm(this%t_dir) = 0
    call crtSparsity_SC(this%sp,this%sp00, &
         TM=tm, ucell=this%ucell, &
         lasto=this%lasto, xa=this%xa, xij=xij)

    ! Notice that we create the correct electrode transfer hamiltonian...
    if ( this%inf_dir == INF_NEGATIVE ) then
       tm(this%t_dir) = -1
    else if ( this%inf_dir == INF_POSITIVE ) then
       tm(this%t_dir) =  1
    else
       call die('Electrode direction not recognized')
    end if
    call crtSparsity_SC(this%sp,this%sp01, &
         TM=tm, ucell=this%ucell, &
         lasto=this%lasto, xa=this%xa, xij=xij)
    
    ! create data
    call newdSpData2D(this%sp00,this%nspin,fdist,this%H00,name='E spH00')
    H00 => val(this%H00)
    call newdSpData2D(this%sp01,this%nspin,fdist,this%H01,name='E spH01')
    H01 => val(this%H01)
    call newdSpData1D(this%sp00,fdist,this%S00,name='E spS00')
    S00 => val(this%S00)
    call newdSpData1D(this%sp01,fdist,this%S01,name='E spS01')
    S01 => val(this%S01)
    call newdSpData2D(this%sp00,3,fdist,this%xij00,name='E spxij00', &
         sparsity_dim=2)
    xij00 => val(this%xij00)
    call newdSpData2D(this%sp01,3,fdist,this%xij01,name='E spxij01', &
         sparsity_dim=2)
    xij01 => val(this%xij01)

    if ( lcalc_xijo ) then
       call newdSpData2D(this%sp00,3,fdist,this%xijo00,name='E spxijo00', &
            sparsity_dim=2)
       xijo00 => val(this%xijo00)
       call newdSpData2D(this%sp01,3,fdist,this%xijo01,name='E spxijo01', &
            sparsity_dim=2)
       xijo01 => val(this%xijo01)
    end if

    call attach(this%sp,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=no_l)
    call attach(this%sp00,n_col=ncol00,list_ptr=ptr00,list_col=col00, &
         nrows=iio)
    if ( iio /= no_l ) call die('Could not do index matching due to &
         &inconsistent sparsity patterns')
    call attach(this%sp01,n_col=ncol01,list_ptr=ptr01,list_col=col01, &
         nrows=iio)
    if ( iio /= no_l ) call die('Could not do index matching due to &
         &inconsistent sparsity patterns')

    ! loop and assign data elements
    do i = 1 , no_l

       ! Shift out of the buffer region
       iio = index_local_to_global(fdist,i,Node)
       ia = iaorb(iio,this%lasto)

       ! Loop number of entries in the row...
       do j = 1 , ncol00(i)

          ! The index in the pointer array is retrieved
          ind00 = ptr00(i) + j

          ! Loop in the super-set sparsity pattern
          idx00: do ind = l_ptr(i) + 1 , l_ptr(i) + l_ncol(i)

             ! If we have the same column index it must be
             ! the same entry they represent
             if ( col00(ind00) == l_col(ind) ) then

                H00(ind00,:)   = H(ind,:)
                S00(ind00)     = S(ind)
                xij00(:,ind00) = xij(:,ind)
                
                if ( lcalc_xijo ) then
                   ja = iaorb(col00(ind00),this%lasto)
                   xijo00(:,ind00) = xij(:,ind) - &
                        (this%xa(:,ja) - this%xa(:,ia))
                end if
                exit idx00
             end if

          end do idx00

       end do

       ! Loop number of entries in the row...
       do j = 1 , ncol01(i)

          ! The index in the pointer array is retrieved
          ind01 = ptr01(i) + j

          ! Loop in the super-set sparsity pattern
          idx01: do ind = l_ptr(i) + 1 , l_ptr(i) + l_ncol(i)

             ! If we have the same column index it must be
             ! the same entry they represent
             if ( col01(ind01) == l_col(ind) ) then

                H01(ind01,:)   = H(ind,:)
                S01(ind01)     = S(ind)
                xij01(:,ind01) = xij(:,ind)

                if ( lcalc_xijo ) then
                   ja = iaorb(col01(ind01),this%lasto)
                   xijo01(:,ind01) = xij(:,ind) - &
                        (this%xa(:,ja) - this%xa(:,ia))
                end if
                exit idx01
             end if
             
          end do idx01

       end do
    end do

    if ( IONode .and. lio ) then
       call print_type(this%sp00)
       call print_type(this%sp01)
    end if

  end subroutine create_sp2sp01

  subroutine delete_(this)
    type(Elec), intent(inout) :: this

    call delete(this%H)
    call delete(this%S)
    call delete(this%xij)
    call delete(this%H00)
    call delete(this%H01)
    call delete(this%S00)
    call delete(this%S01)
    call delete(this%xij00)
    call delete(this%xij01)
    call delete(this%xijo00)
    call delete(this%xijo01)
    call delete(this%sp00)
    call delete(this%sp01)
    call delete(this%sp)
    if ( associated(this%xa) ) deallocate(this%xa)
    if ( associated(this%lasto) ) deallocate(this%lasto)
    if ( associated(this%bkpt_cur) ) deallocate(this%bkpt_cur)
    nullify(this%xa,this%lasto,this%bkpt_cur)
    !if ( associated(this%xa_used) ) deallocate(this%xa_used)
    !if ( associated(this%lasto_used) ) deallocate(this%lasto_used)
    !nullify(this%xa_used,this%lasto_used)

  end subroutine delete_

  ! Routine for checking the validity of the electrode against the 
  ! system setup in transiesta
  subroutine check_Elec(this,nspin,na_u,xa,lasto,xa_EPS, &
       kcell,kdispl)

    use intrinsic_missing, only : VNORM
    use parallel, only : IONode
    use units, only : Ang
    use m_ts_io, only : ts_read_TSHS_opt
#ifdef MPI
    use mpi_siesta, only : MPI_Bcast, MPI_Logical, MPI_Comm_World
#endif

    type(Elec), intent(in) :: this
    integer, intent(in) :: nspin,na_u
    real(dp), intent(in) :: xa(3,na_u)
    integer, intent(in) :: lasto(0:na_u)
    real(dp), intent(in) :: xa_EPS
    integer, intent(in), optional :: kcell(3,3) ! optional as then we can control when to check
    real(dp), intent(in), optional :: kdispl(3)

    ! Local variables
    integer :: i,j,k, ia, iaa, this_kcell(3,3)
    logical :: er, this_er, Gamma
    real(dp) :: xa_o(3), this_xa_o(3), ucell(3,3), this_kdispl(3)
    real(dp) :: max_xa(3), cur_xa(3)
    real(dp), pointer :: this_xa(:,:)

#ifdef MPI
    integer :: MPIerror
#endif

    er = .false.

    this_xa => this%xa_used
    xa_o(:) = xa(:,this%idx_na)
    this_xa_o(:) = this_xa(:,1)
    ucell = this%ucell

    max_xa = 0._dp
    iaa = this%idx_na
    do ia = 1 , this%na_used
       
       do k = 0 , this%RepA3 - 1
       do j = 0 , this%RepA2 - 1
       do i = 0 , this%RepA1 - 1
          ! Calculate repetition vector
          cur_xa(1) = sum(ucell(1,:)*(/i,j,k/))
          cur_xa(2) = sum(ucell(2,:)*(/i,j,k/))
          cur_xa(3) = sum(ucell(3,:)*(/i,j,k/))
          
          ! Add the electrode distance for the ELEC atom
          cur_xa(:) = cur_xa(:) + this_xa(:,ia)-this_xa_o(:)
          ! Subtract the SYSTEM position
          cur_xa(:) = cur_xa(:) - xa(:,iaa) + xa_o(:)
          if ( VNORM(cur_xa) > VNORM(max_xa) ) then
             max_xa = cur_xa
             er = er .or. any( abs(max_xa) > xa_EPS )
          end if

          iaa = iaa + 1
       end do
       end do
       end do
    end do

    if ( er ) then
       
       ! We need to write out the correct formatting of the electrodes.
       ! We shift to the first electrode coordinate, and add the 
       ! first system electrode atom. 
       ! If the user has the %block AtomicCoord... in:
       !    LatticeConstant         1. Ang
       !    AtomicCoordinatesFormat    Ang
       ! then it is direct copy paste ! :)
       xa_o(:) = -this_xa(:,1) + xa(:,this%idx_na)

       if ( IONode ) then
          write(*,'(a)') "Coordinates from the electrode repeated &
               &out to an FDF file"
          write(*,'(a,i0,a)') "NOTICE: that these coordinates are &
               &arranged with respect to atom ", this%idx_na," in your FDF file"
          write(*,'(a)') "NOTICE: that you need to add the species label again"
          write(*,'(a,3(tr1,g10.4))') "Maximal offset in position (Ang):",max_xa/Ang
          write(*,'(a)') "For the same species in the electrode you can do:"
          write(*,'(a,/)') "awk '{print $1,$2,$3,1}' <OUT-file>"
          write(*,'(t3,3a20)') "X (Ang)","Y (Ang)","Z (Ang)"
          iaa = this%idx_na
          do ia = 1 , this%na_used
             do k=0,this%RepA3-1
             do j=0,this%RepA2-1
             do i=0,this%RepA1-1
                write(*,'(t2,3(tr1,f20.10))') &
                     (this_xa(1,ia)+xa_o(1)+sum(ucell(1,:)*(/i,j,k/)))/Ang, &
                     (this_xa(2,ia)+xa_o(2)+sum(ucell(2,:)*(/i,j,k/)))/Ang, &
                     (this_xa(3,ia)+xa_o(3)+sum(ucell(3,:)*(/i,j,k/)))/Ang
             end do
             end do
             end do
          end do
          
       end if
     
    end if

    if ( nspin /= this%nspin ) then
       write(*,*)"ERROR: Electrode: "//trim(this%name)
       write(*,*) '  nspin=',nspin,' expected:', this%nspin
       er = .true.
    end if

    this_er = er

    if ( present(kcell) .and. this%kcell_check ) then

       er = .false.

       call ts_read_TSHS_opt(HSFile(this), &
            kscell=this_kcell,kdispl=this_kdispl, &
            Gamma=Gamma, &
            Bcast=.true.)

       ! If the system is not a Gamma calculation, then the file must
       ! not be either (the repetition will only increase the number of
       ! k-points, hence the above)
       do j = 1 , 3
          select case ( j ) 
          case ( 1 ) 
             k = this%RepA1
          case ( 2 )
             k = this%RepA2
          case ( 3 )
             k = this%RepA3
          end select
          if ( j == this%t_dir ) cycle
          do i = 1 , 3
             if ( i == this%t_dir ) cycle
             if ( j == i ) then
                er = er .or. ( this_kcell(i,j) /= kcell(i,j)*k )
             else 
                er = er .or. ( this_kcell(i,j) /= kcell(i,j) )
             end if
          end do
          er = er .or. ( abs(this_kdispl(j) - kdispl(j)) > 1.e-7_dp )
       end do
       
       ! We still require that the offset in the T-direction is the same
       ! is this even necessary?
       er = er .or. ( abs(this_kdispl(this%t_dir) - kdispl(this%t_dir)) > 1.e-7_dp )
       if ( er .and. IONode ) then
          write(*,'(a)') 'Incompatible k-grids...'
          write(*,'(a)') 'Electrode file k-grid:'
          do j = 1 , 3
             write(*,'(3(i4,tr1),f8.4)') (this_kcell(i,j),i=1,3),this_kdispl(j)
          end do
          write(*,'(a)') 'System k-grid:'
          do j = 1 , 3
             write(*,'(3(i4,tr1),f8.4)') (kcell(i,j),i=1,3),kdispl(j)
          end do
          write(*,'(a)') 'Electrode file k-grid should be:'
          this_kcell(:,1) = kcell(:,1) * this%RepA1
          this_kcell(:,2) = kcell(:,2) * this%RepA2
          this_kcell(:,3) = kcell(:,3) * this%RepA3
          do j = 1 , 3
             write(*,'(3(i4,tr1),f8.4)') (this_kcell(i,j),i=1,3),kdispl(j)
          end do
       end if

    else
       
       call ts_read_TSHS_opt(HSFile(this), &
            Gamma=Gamma, &
            Bcast=.true.)

    end if

    er = this_er .or. er

#ifdef MPI
    call MPI_Bcast(er,1,MPI_Logical,0,MPI_Comm_World,MPIerror)
#endif

    if ( Gamma ) then
       write(*,*) 'Electrode : '//trim(this%name)//' is a Gamma-only calculation &
            &this is not feasible.'
       er = .true.
    end if

    if ( er ) then
       call die("The electrode does not conform with the system settings. &
            &Please correct accordingly.")
    end if

  end subroutine check_Elec

  subroutine check_connectivity(this)

    use parallel, only : IONode, Node
    use units, only : eV

    use class_OrbitalDistribution
    use class_Sparsity

    use create_Sparsity_SC
    use geom_helper, only : iaorb, ucorb
#ifdef MPI
    use mpi_siesta
#endif

    type(Elec), intent(inout) :: this

    real(dp), pointer :: H(:,:)
    real(dp), pointer :: S(:)
    real(dp), pointer :: xij(:,:)
    type(OrbitalDistribution), pointer :: fdist
    type(Sparsity) :: sp02

    integer, pointer :: l_ncol(:) => null()
    integer, pointer :: l_ptr(:)  => null()
    integer, pointer :: l_col(:)  => null()
    integer, pointer :: ncol02(:) => null()
    integer, pointer :: ptr02(:)  => null()
    integer, pointer :: col02(:)  => null()

    integer :: no_l, no_u, i, io, j, ind, ind02, ia
    integer :: tm(3)
    real(dp) :: maxH, maxS
    integer :: maxi, maxj, maxia, maxja

    ! Retrieve distribution
    fdist => dist(this%H)

    if ( .not. initialized(this%H) ) then
       call die('check_connectivity: Error in code')
    end if
    H   => val(this%H)
    S   => val(this%S)
    xij => val(this%xij)

    tm(:) = TM_ALL
    if ( this%inf_dir == INF_POSITIVE ) then
       tm(this%t_dir) =  2
    else
       tm(this%t_dir) = -2
    end if
    call crtSparsity_SC(this%sp,sp02, &
         TM=tm, ucell=this%ucell, lasto=this%lasto, &
         xa=this%xa, xij=xij)

    call attach(this%sp,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=no_l,nrows_g=no_u)

    call attach(sp02,n_col=ncol02,list_ptr=ptr02,list_col=col02)

    ! initialize
    maxH = 0._dp
    maxS = 0._dp
    maxi = 0
    maxj = 0

    ! loop and assign data elements
    do i = 1 , no_l

       ! Shift out of the buffer region
       io = index_local_to_global(fdist,i,Node)
       ia = iaorb(io,this%lasto)

       ! Loop number of entries in the row...
       do j = 1 , ncol02(i)

          ! The index in the pointer array is retrieved
          ind02 = ptr02(i) + j

          ! Loop in the super-set sparsity pattern
          idx02: do ind = l_ptr(i) + 1 , l_ptr(i) + l_ncol(i)

             ! If we have the same column index it must be
             ! the same entry they represent
             if ( col02(ind02) == l_col(ind) ) then

                if ( any(abs(H(ind,:)) > maxH) ) then
                   maxH  = maxval(abs(H(ind,:)))
                   maxS  = S(ind)
                   maxi  = io
                   maxia = ia
                   maxj  = col02(ind02)
                   maxja = iaorb(col02(ind02),this%lasto)
                end if
                exit idx02
             end if
             
          end do idx02
          
       end do
       
    end do

    ind = nnzs(sp02)
    ! clean-up
    call delete(sp02)

    if ( .not. IONode ) return

    if ( ind == 0 ) then
       write(*,'(t2,a)') trim(this%name)//' principal cell is perfect!'
    else if ( maxi == 0 ) then
       write(*,'(t2,a,i0,a)') trim(this%name)//' principal cell is extending out &
            &with all zeroes ',ind,' elements'
    else
       write(*,'(t2,a,i0,a)') trim(this%name)//' principal cell is extending out with ',ind,' elements:'
       write(*,'(t5,2(a,i0))') 'Atom ',maxia,' connects with atom ',maxja
       write(*,'(t5,2(a,i0))') 'Orbital ',UCORB(maxi,no_u) , &
            ' connects with orbital ',UCORB(maxj,no_u)
       write(*,'(t5,3(a,i0),a,g10.3,a)') 'Hamiltonian value: |H(',&
            maxi,',',maxj,')|@R=',tm(this%t_dir),' = ',maxH/eV,' eV'
       write(*,'(t5,3(a,i0),a,g10.3)') 'Overlap          :  S(',&
            maxi,',',maxj,')|@R=',tm(this%t_dir),' = ',maxS
    end if

  end subroutine check_connectivity

  ! Routine for checking the validity of the electrode against the 
  ! system setup in transiesta
  subroutine print_Elec(this,na_u,xa)
    use parallel, only : IONode
    use units, only : Ang
    type(Elec), intent(in) :: this
    integer, intent(in) :: na_u
    real(dp), intent(in) :: xa(3,na_u)

    ! Local variables
    integer :: i,j,k, ia, iaa
    real(dp) :: xa_o(3), this_xa_o(3), ucell(3,3)
    real(dp), pointer :: this_xa(:,:)

    if ( .not. IONode ) return

    write(*,*) trim(this%name)//' unit cell (Ang):'
    write(*,'(2(3(f8.4),/),3(f8.4))') this%ucell/Ang
    
    write(*,'(a,t35,a)') &
         " Structure of "//trim(this%name)//" electrode","| System electrode:"
    write(*,'(t3,3a10,''  |'',3a10)') &
         "X (Ang)","Y (Ang)","Z (Ang)", "X (Ang)","Y (Ang)","Z (Ang)"

    this_xa => this%xa_used
    xa_o(:) = xa(:,this%idx_na)
    this_xa_o(:) = this_xa(:,1)
    ucell = this%ucell

    iaa = this%idx_na
    do ia = 1 , this%na_used
       do k = 0 , this%RepA3 - 1
       do j = 0 , this%RepA2 - 1
       do i = 0 , this%RepA1 - 1
          write(*,'(t3,3f10.5,''  |'',3f10.5)') &
               (this_xa(1,ia)-this_xa_o(1)+sum(ucell(1,:)*(/i,j,k/)))/Ang, &
               (this_xa(2,ia)-this_xa_o(2)+sum(ucell(2,:)*(/i,j,k/)))/Ang, &
               (this_xa(3,ia)-this_xa_o(3)+sum(ucell(3,:)*(/i,j,k/)))/Ang, &
               (xa(1,iaa)-xa_o(1))/Ang, &
               (xa(2,iaa)-xa_o(2))/Ang, &
               (xa(3,iaa)-xa_o(3))/Ang
          iaa = iaa + 1
       end do
       end do
       end do
    end do
    
    write(*,*) ! new-line

  end subroutine print_Elec

end module m_ts_electype
