module m_ts_electype

  use precision, only : dp

  use class_Sparsity
  use class_dSpData1D
  use class_dSpData2D
  use m_region

  use m_geom_plane, only: geo_plane_delta
  use m_geom_plane, only: in_basal_Elec => voxel_in_plane_delta
  use m_geom_box, only: geo_box_delta

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

  public :: Elec, Name, Elec_idx
  public :: TotUsedAtoms, TotUsedOrbs
  public :: AtomInElec, OrbInElec
  public :: q_exp, q_exp_all

  public :: fdf_nElec, fdf_elec

  public :: assign, read_Elec
  public :: create_sp2sp01
  public :: print_Elec
  public :: print_settings
  public :: init_Elec_sim
  public :: check_Elec, check_connectivity
  public :: delete

  public :: in_basal_Elec

  public :: operator(.eq.)

  public :: copy_DM

  ! 300 chars for a full path should be fine
  integer, parameter, public :: FILE_LEN = 300
  integer, parameter, public :: NAME_LEN = 50

  integer, parameter, public :: INF_NEGATIVE = 0 ! old 'left'
  integer, parameter, public :: INF_POSITIVE = 1 ! old 'right'

  type :: Elec
     integer :: ID = 0
     character(len=FILE_LEN) :: HSfile = ' ', DEfile = ' ', GFfile  = ' '
     character(len=NAME_LEN) :: Name   = ' '
     ! These variables are relative to the big system
     integer :: idx_a = 0, idx_o = 0
     ! atoms used
     integer :: na_used = 0
     ! orbitals used
     integer :: no_used = 0
     ! repetitions
     integer :: Rep(3) = 1
     ! Pre-expand before saving Gf (we default to all)
     ! In this way will the user 
     integer :: pre_expand = 2
     ! chemical potential of the electrode
     type(ts_mu), pointer :: mu => null()
     ! infinity direction
     integer :: inf_dir = INF_NEGATIVE
     ! transport direction (determines H01)
     ! And is considered with respect to the electrode direction...
     !  t_dir is with respect to the electrode unit-cell
     integer :: t_dir = 3
     ! This is the cell pivoting table from the
     ! electrode unit-cell to the simulation unit-cell
     ! So: cell(:,pvt(1)) ~= this%cell(:,1)
     integer :: pvt(3)
     ! whether the electrode should be bulk
     logical :: Bulk = .true.
     integer :: DM_update = 0 ! This determines the update scheme for the crossterms
                              ! == 0 means no update
                              ! == 1 means update cross-terms
                              ! == 2 means update everything (no matter Bulk)
     ! whether to re-calculate the GF-file
     logical :: ReUseGF = .true.
     ! Create a GF-file or re-calculate the self-energies everytime
     logical :: out_of_core = .true. 
     ! In case of 'out_of_core == .false.' we can reduce the number of operations
     ! by skipping the copying of H00, S00. Hence we need to compare when the 
     ! k-point changes...
     real(dp) :: bkpt_cur(3)
     ! Used xa and lasto
     real(dp), pointer :: xa_used(:,:) => null()
     integer,  pointer :: lasto_used(:) => null()

     ! Advanced stuff...
     logical :: kcell_check = .true.
     ! If the user requests to assign different "spill-in fermi-level" 
     ! we allow that
     real(dp) :: Ef_frac_CT = 0._dp

     ! ---v--- Below we have the content of the TSHS file
     integer  :: nspin = 0, na_u = 0, no_u = 0, no_s = 0
     real(dp) :: cell(3,3) = 0._dp, Ef = 0._dp, Qtot = 0._dp
     real(dp), pointer :: xa(:,:) => null()
     integer,  pointer :: lasto(:) => null()
     type(Sparsity)  :: sp
     type(dSpData2D) :: H
     type(dSpData1D) :: S
     ! Supercell offsets
     integer :: nsc(3)
     integer, pointer :: isc_off(:,:) => null()
     ! --- --- completed the content of the TSHS file
     ! Below we create the content for the self-energy creation
     ! Notice that we can save some elements simply by extracting the 0-1 connections
     ! for large systems this is a non-negligeble part of the memory...
     type(Sparsity)  :: sp00, sp01
     type(dSpData2D) :: H00, H01
     type(dSpData1D) :: S00, S01

     ! These arrays are used to construct the full Hamiltonian and overlap and Green function
     complex(dp), pointer :: HA(:,:,:), SA(:,:,:), GA(:)

     ! Arrays needed to partition the scattering matrix and self-energies

     ! Gamma stored is actually this: (Sigma - Sigma^\dagger) ^ T
     ! and NOT: i (Sigma - Sigma^\dagger)
     complex(dp), pointer :: Gamma(:), Sigma(:)

     ! The imaginary part in the electrode
     real(dp) :: Eta = 7.3498067e-7_dp ! corresponds to 0.00001 eV

     ! The region of the down-folded region
     type(tRgn) :: o_inD, inDpvt

     ! The basal plane of the electrode
     type(geo_plane_delta) :: p
     ! A box containing all atoms of the electrode in the
     ! simulation box
     type(geo_box_delta) :: box

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
       this_n(n)%ID = n
       if ( index(this_n(n)%name,'.') > 0 ) then
          call die('Electrodes cannot be named with .!')
       end if
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

  function fdf_Elec(prefix,slabel,this,N_mu,mus,idx_a,name_prefix) result(found)
    use fdf
    use m_os, only : file_exist
    use m_ts_io, only : ts_read_TSHS_opt
    use m_ts_io_ctype, only : pline_E_parse

    character(len=*), intent(in) :: prefix,slabel
    type(Elec), intent(inout) :: this
    integer, intent(in) :: N_mu
    type(ts_mu), intent(in), target :: mus(N_mu)
    integer, intent(in), optional :: idx_a
    character(len=*), intent(in), optional :: name_prefix

    logical :: found

    ! prepare to read in the data...
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline => null()
    logical :: info(5)
    integer :: i, j
    integer :: cidx_a 

    character(len=200) :: bName, name, ln, tmp

    info(:) = .false.

    bName = trim(prefix)//'.Elec.'//trim(this%name)
    found = fdf_block(trim(bName),bfdf)

    ! Allow the filename to be read in individually
    name = trim(bName)//'.TSHS'
    if ( fdf_defined(trim(name)) ) then
       this%HSfile = trim(fdf_get(name,''))
       info(1) = .true.
    end if
    name = trim(bName)//'.Rep.A1'
    if ( fdf_defined(trim(name)) ) this%Rep(1) = fdf_get(name,1)
    name = trim(bName)//'.Rep.A2'
    if ( fdf_defined(trim(name)) ) this%Rep(2) = fdf_get(name,1)
    name = trim(bName)//'.Rep.A3'
    if ( fdf_defined(trim(name)) ) this%Rep(3) = fdf_get(name,1)
    name = trim(bName)//'.GF'
    if ( fdf_defined(trim(name)) ) this%GFfile = trim(fdf_get(name,''))

    if ( .not. found ) return

    cidx_a = 0

    if ( present(idx_a) ) then
       if ( idx_a /= 0 ) then
          this%idx_a = idx_a
          if ( idx_a < 0 ) cidx_a = -1
          info(4) = .true.
       end if
    end if

    ! We default a lot of the options
    if ( len_trim(this%GFfile) == 0 ) then
       name = 'GF'//trim(this%name)
       if ( present(name_prefix) ) then
          this%GFfile = trim(slabel)//'.'//trim(name_prefix)//trim(name)
       else
          this%GFfile = trim(slabel)//'.'//trim(prefix)//trim(name)
       end if
    end if
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
          cidx_a     = 0
          this%idx_a = 0
          if ( fdf_bnnames(pline) > 1 ) then
             ! the user is requesting on a string basis
             ln = fdf_bnames(pline,2)
             if ( leqi(ln,'start') .or. leqi(ln,'begin') ) then
                if ( fdf_bnintegers(pline) > 0 ) then
                   this%idx_a = fdf_bintegers(pline,1)
                else
                   this%idx_a = 1 ! default starting position
                end if
             else if ( leqi(ln,'end') ) then
                cidx_a = -1
                if ( fdf_bnintegers(pline) > 0 ) then
                   this%idx_a = fdf_bintegers(pline,1)
                else
                   this%idx_a = -1
                end if
             end if
          else
             if ( fdf_bnintegers(pline) < 1 ) &
                  call die('Atomic position not found in input line: '//trim(ln))
             this%idx_a = fdf_bintegers(pline,1)
          end if
          info(4) = .true.
          
       else if ( leqi(ln,'bulk') ) then
          this%Bulk = fdf_bboolean(pline,1,after=1)

       else if ( leqi(ln,'pre-expand') ) then
          tmp = fdf_bnames(pline,2)
          if ( leqi(tmp,'all') .or. leqi(tmp,'everything') ) then
             ! We expand H, S and GS before writing...
             this%pre_expand = 2
          else if ( leqi(tmp,'Green') .or. &
               leqi(tmp,'surface') ) then
             ! We expand only the surface Green function
             this%pre_expand = 1
          else if ( leqi(tmp,'none') ) then
             ! We do not expand anything
             this%pre_expand = 0
          else
             call die('Error in option ''pre-expand'', please &
                  &correct. Must be: [all|everything,Green|Surface,none]!')
          end if

       else if ( leqi(ln,'DM-update') ) then
          if ( fdf_bnnames(pline) < 2 ) &
               call die('Update scheme not supplied')

          tmp = fdf_bnames(pline,2)
          if ( leqi(tmp,'none') ) then
             this%DM_update = 0
          else if ( leqi(tmp,'cross-terms') ) then
             this%DM_update = 1
          else if ( leqi(tmp,'all') ) then
             this%DM_update = 2
          else
             call die('DM-update: unrecognized option: '//trim(tmp))
          end if
          info(5) = .true.

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

       else if ( leqi(ln,'GF') .or. &
            leqi(ln,'GF-file') ) then
          if ( fdf_bnnames(pline) < 2 ) call die('GF-file not supplied')
          this%GFfile = trim(fdf_bnames(pline,2))

       else if ( leqi(ln,'GF.ReUse') ) then

          this%ReUseGF = fdf_bboolean(pline,1,after=1)

       else if ( leqi(ln,'used-atoms') ) then
          if ( fdf_bnintegers(pline) < 1 ) &
               call die('Number of atoms used not supplied')
          this%na_used = fdf_bintegers(pline,1)

       else if ( leqi(ln,'replicate-a') .or. leqi(ln,'rep-a') .or. &
            leqi(ln,'replicate-a1') .or. leqi(ln,'rep-a1') ) then
          if ( fdf_bnintegers(pline) < 1 ) &
               call die('Repetition A1 is not supplied')
          this%Rep(1) = fdf_bintegers(pline,1)

       else if ( leqi(ln,'replicate-b') .or. leqi(ln,'rep-b') .or. &
            leqi(ln,'replicate-a2') .or. leqi(ln,'rep-a2') ) then
          if ( fdf_bnintegers(pline) < 1 ) &
               call die('Repetition A2 is not supplied')
          this%Rep(2) = fdf_bintegers(pline,1)

       else if ( leqi(ln,'replicate-c') .or. leqi(ln,'rep-c') .or. &
            leqi(ln,'replicate-a3') .or. leqi(ln,'rep-a3') ) then
          if ( fdf_bnintegers(pline) < 1 ) &
               call die('Repetition A3 is not supplied')
          this%Rep(3) = fdf_bintegers(pline,1)

       else if ( leqi(ln,'replicate') .or. leqi(ln,'rep') ) then
          if ( fdf_bnintegers(pline) < 3 ) &
               call die('Repetition for all directions are not supplied <A1> <A2> <A3>')
          this%Rep(1) = fdf_bintegers(pline,1)
          this%Rep(2) = fdf_bintegers(pline,2)
          this%Rep(3) = fdf_bintegers(pline,3)

#ifdef TBTRANS
       else if ( leqi(ln,'tbt.out-of-core') .or. &
            leqi(ln,'out-of-core') ) then
#else
       else if ( leqi(ln,'out-of-core') ) then
#endif

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

       else if ( leqi(ln,'TSDE') .or. &
            leqi(ln,'TSDE-file') ) then
          if ( fdf_bnnames(pline) < 2 ) call die('TSDE name not supplied')
          this%DEfile = trim(fdf_bnames(pline,2))

#ifdef TBTRANS
       else if ( leqi(ln,'tbt.Eta') .or. leqi(ln,'Eta') ) then
#else
       else if ( leqi(ln,'Eta') ) then
#endif
          call pline_E_parse(pline,1,ln, &
               val = this%Eta, before=3)

#ifdef TBTRANS
       else if ( leqi(ln,'tbt.GF') .or. &
            leqi(ln,'tbt.GF-file') ) then
          if ( fdf_bnnames(pline) < 2 ) call die('tbt.GF-file not supplied')
          this%GFfile = trim(fdf_bnames(pline,2))

       else if ( leqi(ln,'tbt.GF.ReUse') ) then

          this%ReUseGF = fdf_bboolean(pline,1,after=1)

#else
       else if ( leqi(ln(1:3),'tbt') ) then
          ! by-pass
          ! All options that are meant for tbtrans are discarded :)
#endif

       else

          ! we should always die in case something non-understandable 
          ! is passed, if that is the case, the chances are that the
          ! user has made a typo is high
          call die('Unrecognized option "'//trim(ln)//'" &
               &for electrode: '//trim(name))

       end if

    end do
    
    if ( any(this%Rep(:) < 1) ) &
         call die("Repetition in "//trim(name)//" electrode must be >= 1.")

    if ( .not. all(info(1:4)) ) then
       write(*,*)'You need to supply at least:'
       write(*,*)' - TSHS'
       write(*,*)' - semi-inf-direction'
       write(*,*)' - chemical-potential'
       write(*,*)' - electrode-position'
       call die('You have not supplied all electrode information')
    end if

#ifndef TBTRANS
    if ( this%Eta <= 0._dp ) then
       call die('We do not allow the advanced Green function &
            &to be calculated. Please ensure a positive imaginary &
            &part of the energy (non-zero).')
    end if
#endif
    
    if ( .not. file_exist(this%HSfile, Bcast = .true.) ) then
       call die("Electrode file does not exist. &
            &Please create electrode '"//trim(this%HSfile)//"' first.")
    end if

    ! Read in the number of atoms in the HSfile
    call ts_read_TSHS_opt(this%HSfile,no_u=this%no_u,na_u=this%na_u, &
         nspin=this%nspin, Ef=this%Ef, ucell=this%cell, Qtot=this%Qtot, &
         nsc = this%nsc , &
         Bcast=.true.)

    allocate(this%xa(3,this%na_u),this%lasto(0:this%na_u))
    call ts_read_TSHS_opt(this%HSfile,xa=this%xa,lasto=this%lasto, &
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
    if ( this%Rep(this%t_dir) /= 1 ) then
       call die('Repetition in the transport direction &
            &is not allowed.')
    end if

    ! if the electrode does not use a bulk electrode we need to update
    ! the cross-terms
    if ( .not. this%Bulk ) then
       this%DM_update = 2
    end if

                                 ! Same criteria as IsVolt
    if ( this%DM_update > 0 .or. abs(this%mu%mu) < 0.000000735_dp ) then
       ! when updating the cross-terms
       ! there is no need to also shift the energy
       ! I think this will battle each other out...
       this%Ef_frac_CT = 0._dp
    end if

    ! if the user has specified text for the electrode position
    if ( cidx_a == -1 ) then
       this%idx_a = this%idx_a + 1 - TotUsedAtoms(this)
    end if

    ! For in-core calculations of the Self-energy,
    ! the pre-expansion does not make sense.
    ! Hence, we immediately set it to 0
    if ( .not. this%out_of_core ) then
       this%pre_expand = 0
    end if

    ! In case the user has not supplied a DM file for the
    ! electrode we might as well try and guess one... :)
    if ( len_trim(this%DEfile) == 0 ) then
       i = len_trim(this%HSfile)
#ifdef NCDF_4
       ! If the TSHS file is a netcdf file we re-use it
       if ( this%HSfile(i-1:i) == 'nc' ) then
          this%DEfile = this%HSfile
       else
          this%DEfile = this%HSfile(1:i-4)//'TSDE'
       end if
#else
       this%DEfile = this%HSfile(1:i-4)//'TSDE'
#endif
    end if

  end function fdf_Elec

  ! Initialize variables for the electrode according
  ! to the simulation variables
  subroutine init_Elec_sim(this,cell,na_u,xa)

    use parallel, only : IONode
    use intrinsic_missing, only : VNORM, SPC_PROJ, VEC_PROJ, IDX_SPC_PROJ

    ! The electrode that needs to be processed
    type(Elec), intent(inout) :: this
    ! The simulation parameters.
    ! Unit-cell of simulation cell
    real(dp), intent(in) :: cell(3,3)
    integer, intent(in) :: na_u
    real(dp), intent(in) :: xa(3,na_u)
    
    real(dp) :: p(3), max_xa(3), contrib, min_bond
    integer :: i, j, n_cell, ia, na

    ! First figure out the minimum bond-length
    ia = this%idx_a
    na = TotUsedAtoms(this)
    if ( na == 1 ) then
       call die('One atom electrodes are not allowed')
    end if

    min_bond = huge(1._dp)
    do i = 1 , na - 1
       p(:) = xa(:,ia-1+i)
       do j = i + 1 , na
          contrib = VNORM(xa(:,ia-1+j)-p)
          if ( contrib < min_bond ) then
             min_bond = contrib
          end if
       end do
    end do
    ! Take half the bond-length to get in between two atoms
    min_bond = min_bond * 0.5_dp

    ! Calculate the pivoting table
    do i = 1 , 3
       this%pvt(i) = IDX_SPC_PROJ( cell,this%cell(:,i) )
    end do
    if ( sum(this%pvt) /= 6 .or. count(this%pvt==2) /= 1 ) then
       print *, this%pvt
       call die('The pivoting table for the electrode unit-cell, &
            &onto the simulation unit-cell is not unique. &
            &Please check your simulation cells.')
    end if

    ! Print out a warning if the electrode uses several cell-vectors
    p = SPC_PROJ(cell,this%cell(:,this%t_dir))
    n_cell = 0
    do i = 1 , 3 

       ! project the unit-cell vector onto each cell component
       contrib = VNORM(VEC_PROJ(cell(:,i),p))

       ! If the contribution in this cell direction is too
       ! small we consider it not to be important.
       ! TODO this might in certain skewed examples be a bad choice.
       if ( contrib < 1.e-6_dp ) cycle
       n_cell = n_cell + 1
    end do
    if ( n_cell > 1 .and. IONode ) then
       ! The electrode uses more than one cell-vector
       ! to describe the transport direction.
       write(*,'(a)')' *** Electrode '//trim(this%name)//' has the transport direction'
       write(*,'(a,i0)')'     aligning with ',n_cell,' cell vectors.'
       write(*,'(a)')'     Responsibility has been relieved of transiesta/tbtrans!'
    end if

    ! The cell-vector along the transport direction.
    p = this%cell(:,this%t_dir)
    ! We add a vector with length of half the minimal bond length
    ! to the vector, to do the averaging 
    ! not on-top of an electrode atom.
    p = p / VNORM(p) * min_bond
    
    ! Create the basal plane of the electrode
    ! Decide which end of the electrode we use
    ! TODO, correct for systems not having the last electrode atom
    ! farthest from the device region (or correct intrinsically)
    if ( this%inf_dir == INF_POSITIVE ) then
       ! We need to utilize the last atom
       this%p%c = xa(:,ia+na-1)
       this%p%c = this%p%c + p ! add vector
    else
       this%p%c = xa(:,ia)
       this%p%c = this%p%c - p ! subtract vector
    end if
    
    ! Normal vector to electrode transport direction
    this%p%n = SPC_PROJ(cell,this%cell(:,this%t_dir))
    this%p%n = this%p%n / VNORM(this%p%n) ! normalize

    ! The distance parameter
    this%p%d = sum( this%p%n(:)*this%p%c(:) )

    ! Create the box cell for the electrode
    ! This enables the Hartree correction for the 
    ! electrode region
    ! The box-cell will be extended by the bond-length
    ! to make "as big as box as viable"
    min_bond = min_bond * 2._dp ! re-scale to actual bond-length
    do i = 1 , 3
       ! Get the lower left corner of the electrode
       p(i) = minval(xa(i,ia:ia+na-1))
       ! The box extends back one bond-length in each direction.
       p(i) = p(i) - min_bond
       ! Get the maximum in each direction
       ! In principle this should probably be the electrode
       ! unit-cell vectors, however, for non-periodic electrodes
       ! this might not be such a good choice.
       ! The best thing would be to check the periodic directions,
       ! and in those directions assign the vectors in those directions
       ! always having the third vector equalling the cell transport
       ! direction.
       ! NOTE TODO : this is for now a "stupid" implementation.
       max_xa(i) = maxval(xa(i,ia:ia+na-1))
       max_xa(i) = max_xa(i) + min_bond
       this%box%v(:,i) = 0._dp
       this%box%v(i,i) = max_xa(i) - p(i)

    end do
    do i = 1 , 3
       ! If there is no periodicity of the electrode in this
       ! direction we move the center of the box by the 
       ! unit-cell. This ensures that for non-periodic
       ! cells, we still have an equal weight of the
       ! lifting of the Hartree potential.
       if ( this%nsc(i) == 1 ) then
          p(:) = p(:) - 0.5_dp * this%cell(:,i)
       end if
    end do
    ! p now contains the origin of the box
    this%box%c(:) = p(:)
    ! For now the box vectors are defined using the unit-cell of
    ! the electrode.
    ! Note that we extend the unit-cell by the bond-length
    ! We will really try to extend the box to be
    ! as large as possible to retain the potential surrounding
    ! the electrode as much as possible.
    ! This is really important for periodic calculations
    ! where it is necessary to have the lifting potential
    ! in the entire box of the electrode
    do i = 1 , 3
       ! Remember that we start one bond-length "below"
       ! and the unit-cell is always "one bond length" too
       ! long (for periodic reasons)
       p(:) = min_bond
       ! Create a unit-vector along the unit-cell direction
       max_xa(:) = this%cell(:,i) / VNORM(this%cell(:,i))
       p(:) = VEC_PROJ(max_xa,p)
       p(:) = this%cell(:,i) + p(:)
       this%box%v(:,i) = SPC_PROJ(cell,p)
    end do

  end subroutine init_Elec_sim

  function equal_el_el(this1,this2) result(equal)
    type(Elec), intent(in) :: this1, this2
    logical :: equal
    equal = trim(this1%name) == trim(this2%name)
  end function equal_el_el

  function equal_el_str(this,str) result(equal)
    type(Elec), intent(in) :: this
    character(len=*), intent(in) :: str
    logical :: equal
    equal = trim(this%name) == trim(str)
  end function equal_el_str

  function equal_str_el(str,this) result(equal)
    character(len=*), intent(in) :: str
    type(Elec), intent(in) :: this
    logical :: equal
    equal = trim(this%name) == trim(str)
  end function equal_str_el

  subroutine assign(this,D,HSfile,GFfile, &
       na_u,na_used,no_u,no_s,no_used, &
       RepA1,RepA2,RepA3)
    type(Elec), intent(inout) :: this
    character, optional, intent(in) :: D
    character(len=*), intent(in), optional :: HSfile, GFfile
    integer, intent(in), optional :: na_u, na_used, no_u, no_s, no_used
    integer, intent(in), optional :: RepA1, RepA2, RepA3

    if (present(D)) call die('Wrong usage of Elec-type assign')

    if (present(HSfile)) this%HSfile = HSfile
    if (present(GFfile)) this%GFfile = GFfile
    if (present(na_u)) this%na_u = na_u
    if (present(na_used)) this%na_used = na_used
    if (present(no_u)) this%no_u = no_u
    if (present(no_s)) this%no_s = no_s
    if (present(no_used)) this%no_used = no_used
    if (present(RepA1)) this%Rep(1) = RepA1
    if (present(RepA2)) this%Rep(2) = RepA2
    if (present(RepA3)) this%Rep(3) = RepA3

  end subroutine assign
  
  elemental function Name_(this) result(name)
    type(Elec), intent(in) :: this
    character(len=NAME_LEN) :: Name
    name = this%name
  end function Name_

  elemental function TotUsedAtoms(this) result(val)
    type(Elec), intent(in) :: this
    integer :: val
    val = this%na_used * product(this%Rep)
  end function TotUsedAtoms

  pure function q_exp_all(this,i,j,k) result(q)
    type(Elec), intent(in) :: this
    integer, intent(in) :: i,j,k
    real(dp) :: q(3)
    
    ! TODO, the current implementation assumes k-symmetry!
    ! Hence, using repetition with non-symmetry will produce
    ! wrong results.
    ! Luckily this is not a problem currently.
    ! Perhaps one should consider this in tbtrans

    q(1) = 1._dp*(i-1) / real(this%Rep(1),dp)
    q(2) = 1._dp*(j-1) / real(this%Rep(2),dp)
    q(3) = 1._dp*(k-1) / real(this%Rep(3),dp)
  end function q_exp_all

  pure function q_exp(this,idx) result(q)
    type(Elec), intent(in) :: this
    integer, intent(in) :: idx
    real(dp) :: q(3)
    integer :: i,j,k,ii
    i =     this%Rep(1)
    j = i * this%Rep(2)
    k = j * this%Rep(3)
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
    val = this%no_used * product(this%Rep)
  end function TotUsedOrbs

  elemental function OrbInElec(this,io) result(in)
    type(Elec), intent(in) :: this
    integer, intent(in) :: io
    logical :: in
    in = this%idx_o <= io .and. io < (this%idx_o + TotUsedOrbs(this))
  end function OrbInElec

  elemental function AtomInElec(this,ia) result(in)
    type(Elec), intent(in) :: this
    integer, intent(in) :: ia
    logical :: in
    in = this%idx_a <= ia .and. ia < (this%idx_a + TotUsedAtoms(this))
  end function AtomInElec

  subroutine read_Elec(this,Bcast,io,ispin)
    use fdf
    use parallel
    use class_OrbitalDistribution

    use m_handle_sparse, only : reduce_spin_size
    use m_ts_io
#ifdef MPI
    use mpi_siesta
#endif

    type(Elec), intent(inout) :: this
    logical, intent(in), optional :: Bcast ! Bcast information
    logical, intent(in), optional :: IO ! Write to STD-out
    integer, intent(in), optional :: ispin ! select one spin-channel

    character(len=200) :: fN
    integer :: fL, kscell(3,3), istep, ia1
    logical :: onlyS, Gamma_file, TSGamma, lio
    real(dp) :: Temp, kdispl(3), Qtot, Ef
    
    ! Sparsity pattern
    integer :: nsc(3)

    lio = .true.
    if ( present(io) ) lio = io

    fN = trim(this%HSfile)
    ! We read in the information
    fL = len_trim(fN)
    call ts_read_tshs(fN, &
         onlyS, Gamma_file, TSGamma, &
         this%cell, nsc, this%na_u, this%no_u, this%nspin,  &
         kscell, kdispl, &
         this%xa, this%lasto, &
         this%sp, this%H, this%S, this%isc_off, &
         Ef, Qtot, Temp, &
         istep, ia1, &
         Bcast=Bcast)

    if ( present(ispin) ) then
       if ( ispin > 0 ) then
          call reduce_spin_size(ispin,this%H)
          this%nspin = 1
       end if
    end if

    if ( IONode .and. lio ) call print_type(this%sp)

  end subroutine read_Elec

  subroutine create_sp2sp01(this,IO)

    use parallel, only : IONode

    use class_OrbitalDistribution

    use create_Sparsity_SC
    use geom_helper, only : iaorb
#ifdef MPI
    use mpi_siesta
#endif

    type(Elec), intent(inout) :: this
    logical, intent(in), optional :: io

    logical :: lio

    real(dp), pointer :: H(:,:), H00(:,:), H01(:,:)
    real(dp), pointer :: S(:), S00(:), S01(:)
    type(OrbitalDistribution), pointer :: fdist

    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer, pointer :: ncol00(:), ptr00(:), col00(:)
    integer, pointer :: ncol01(:), ptr01(:), col01(:)

    integer :: no_l, i, iio, j, ind, ind00, ind01, ia
    integer :: tm(3)

    ! Retrieve distribution
    fdist => dist(this%H)

    lio = .true.
    if ( present(io) ) lio = io

    H => val(this%H)
    S => val(this%S)
    tm(:)          = TM_ALL
    tm(this%t_dir) = 0
    call crtSparsity_SC(this%sp,this%sp00, &
         TM=tm, ucell=this%cell, &
         isc_off=this%isc_off)

    ! Notice that we create the correct electrode transfer hamiltonian...
    if ( this%inf_dir == INF_NEGATIVE ) then
       tm(this%t_dir) = -1
    else if ( this%inf_dir == INF_POSITIVE ) then
       tm(this%t_dir) =  1
    else
       call die('Electrode direction not recognized')
    end if
    call crtSparsity_SC(this%sp,this%sp01, &
         TM=tm, ucell=this%cell, &
         isc_off=this%isc_off)
    
    ! create data
    call newdSpData2D(this%sp00,this%nspin,fdist,this%H00,name='E spH00')
    H00 => val(this%H00)
    call newdSpData2D(this%sp01,this%nspin,fdist,this%H01,name='E spH01')
    H01 => val(this%H01)
    call newdSpData1D(this%sp00,fdist,this%S00,name='E spS00')
    S00 => val(this%S00)
    call newdSpData1D(this%sp01,fdist,this%S01,name='E spS01')
    S01 => val(this%S01)

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
       iio = index_local_to_global(fdist,i)
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

                H00(ind00,:) = H(ind,:)
                S00(ind00)   = S(ind)

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

                H01(ind01,:) = H(ind,:)
                S01(ind01)   = S(ind)

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
    call delete(this%H00)
    call delete(this%H01)
    call delete(this%S00)
    call delete(this%S01)
    call delete(this%sp00)
    call delete(this%sp01)
    call delete(this%sp)
    if ( associated(this%xa) ) deallocate(this%xa)
    if ( associated(this%lasto) ) deallocate(this%lasto)
    nullify(this%xa,this%lasto)
    !if ( associated(this%xa_used) ) deallocate(this%xa_used)
    !if ( associated(this%lasto_used) ) deallocate(this%lasto_used)
    !nullify(this%xa_used,this%lasto_used)

  end subroutine delete_

  ! Routine for checking the validity of the electrode against the 
  ! system setup in transiesta
  subroutine check_Elec(this,nspin,ucell,na_u,xa,lasto,xa_EPS, &
       kcell,kdispl)

    use intrinsic_missing, only : VNORM, SPC_PROJ, VEC_PROJ
    use parallel, only : IONode
    use units, only : Ang
    use m_ts_io, only : ts_read_TSHS_opt
#ifdef MPI
    use mpi_siesta, only : MPI_Bcast, MPI_Logical, MPI_Comm_World
#endif

    type(Elec), intent(inout) :: this
    integer, intent(in) :: nspin,na_u
    real(dp), intent(in) :: ucell(3,3)
    real(dp), intent(in) :: xa(3,na_u)
    integer, intent(in) :: lasto(0:na_u)
    real(dp), intent(in) :: xa_EPS
    integer, intent(in), optional :: kcell(3,3) ! optional as then we can control when to check
    real(dp), intent(in), optional :: kdispl(3)

    ! Local variables
    integer :: i,j,k, ia, iaa, this_kcell(3,3)
    logical :: ldie, er, Gamma
    real(dp) :: xa_o(3), this_xa_o(3), cell(3,3), this_kdispl(3)
    real(dp) :: max_xa(3), cur_xa(3)
    real(dp), pointer :: this_xa(:,:)
    integer :: pvt(3)
#ifdef MPI
    integer :: MPIerror
#endif

    ldie = .false.

    this_xa => this%xa_used
    xa_o(:) = xa(:,this%idx_a)
    this_xa_o(:) = this_xa(:,1)
    cell = this%cell
    pvt = this%pvt

    max_xa = 0._dp
    iaa = this%idx_a
    er = .false.
    do ia = 1 , this%na_used
       
       do k = 0 , this%Rep(3) - 1
       do j = 0 , this%Rep(2) - 1
       do i = 0 , this%Rep(1) - 1
          ! Calculate repetition vector
          cur_xa(1) = sum(cell(1,:)*(/i,j,k/))
          cur_xa(2) = sum(cell(2,:)*(/i,j,k/))
          cur_xa(3) = sum(cell(3,:)*(/i,j,k/))
          
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

    ldie = ldie .or. er

    if ( er ) then
       
       ! We need to write out the correct formatting of the electrodes.
       ! We shift to the first electrode coordinate, and add the 
       ! first system electrode atom. 
       ! If the user has the %block AtomicCoord... in:
       !    LatticeConstant         1. Ang
       !    AtomicCoordinatesFormat    Ang
       ! then it is direct copy paste ! :)
       xa_o(:) = -this_xa(:,1) + xa(:,this%idx_a)

       if ( IONode ) then
          i = this%idx_a
          j = this%idx_a + TotUsedAtoms(this) - 1
          write(*,'(a,e10.5,a)') 'The electrode coordinates does not overlap within the &
               &required accuracy: ',xa_EPS/Ang,' Ang'
          write(*,'(a,3(tr1,g10.4))') 'The maximal offset vector is (Ang):',max_xa/Ang
          write(*,'(2(a,i0),2a)') 'The system coordinates of atoms ',i,' to ',j,  &
               ' does not coincide with the electrode coordinates found in: ',trim(this%HSfile)
          write(*,'(a)') 'This is a requirement to place the self-energy terms correctly.'
          write(*,'(a,/,2(a,i0),a)') 'To ensure the electrode coordinates conform with &
               &the system coordinates you can take the following list of coordinates','and &
               &replace atoms ',i,' to ',j,' in your system AtomicCoordinatesAndAtomicSpecies block'
          write(*,'(a,i0,a)') 'NOTICE: The listed coordinates are already arranged &
               &with respect to atom ',i,' in your system AtomicCoordinatesAndAtomicSpecies block'
          write(*,'(a)') 'NOTICE: You have to add the correct species label for the atoms'
          write(*,'(a)') 'NOTICE: You can possibly do this by using this awk-command on this output:'
          write(*,'(a,/)') "        awk '{print $1,$2,$3,1}' <OUT-file>"

          write(*,'(t3,3a20)') 'X (Ang)','Y (Ang)','Z (Ang)'
          iaa = this%idx_a
          do ia = 1 , this%na_used
             do k = 0 , this%Rep(3)-1
             do j = 0 , this%Rep(2)-1
             do i = 0 , this%Rep(1)-1
                write(*,'(t2,3(tr1,f20.10))') &
                     (this_xa(1,ia)+xa_o(1)+sum(cell(1,:)*(/i,j,k/)))/Ang, &
                     (this_xa(2,ia)+xa_o(2)+sum(cell(2,:)*(/i,j,k/)))/Ang, &
                     (this_xa(3,ia)+xa_o(3)+sum(cell(3,:)*(/i,j,k/)))/Ang
             end do
             end do
             end do
          end do
          
       end if
       
    end if

    iaa = this%idx_a
    er = .false.
    do ia = 1 , this%na_used
       do k = 0 , this%Rep(3) - 1
       do j = 0 , this%Rep(2) - 1
       do i = 0 , this%Rep(1) - 1

          ! Check number of orbitals for this electrode
          ! atom
          if ( lasto(iaa) - lasto(iaa-1) /= &
               this%lasto_used(ia) - this%lasto_used(ia-1) ) then
             er = .true.
          end if

          iaa = iaa + 1
       end do
       end do
       end do
    end do
    
    ldie = ldie .or. er
    if ( er ) then
       if ( IONode ) then
          write(*,'(a)') "Number of orbitals per atom in the electrode does not match the system electrode"
          write(*,'(a)') 'Have you changed your basis size?'
          write(*,'(t3,3a20)') "ia system","n_orb_el","n_orb_sys"
          iaa = this%idx_a
          do ia = 1 , this%na_used
             do k = 0 , this%Rep(3) - 1
             do j = 0 , this%Rep(2) - 1
             do i = 0 , this%Rep(1) - 1

              ! Check number of orbitals for this electrode
              ! atom
                write(*,'(t3,3(i20))')iaa, &
                     this%lasto_used(ia) - this%lasto_used(ia-1), &
                     lasto(iaa) - lasto(iaa-1)
                iaa = iaa + 1
             end do
             end do
             end do
          end do
       end if

    end if

    if ( nspin /= this%nspin ) then
       write(*,*)"ERROR: Electrode: "//trim(this%name)
       write(*,*) '  nspin=',nspin,' expected:', this%nspin
       ldie = .true.
    end if

    er = .false.
    if ( present(kcell) .and. this%kcell_check ) then

       call ts_read_TSHS_opt(this%HSfile, &
            kscell=this_kcell,kdispl=this_kdispl, &
            Gamma=Gamma, &
            Bcast=.true.)

       ! If the system is not a Gamma calculation, then the file must
       ! not be either (the repetition will only increase the number of
       ! k-points, hence the above)
       do j = 1 , 3
          k = this%Rep(j)
          if ( j == this%t_dir ) cycle
          do i = 1 , 3
             if ( i == this%t_dir ) cycle
             if ( j == i ) then
                er = er .or. ( this_kcell(i,j) /= kcell(pvt(i),pvt(j))*k )
             else 
                er = er .or. ( this_kcell(i,j) /= kcell(pvt(i),pvt(j)) )
             end if
          end do
          er = er .or. ( abs(this_kdispl(j) - kdispl(pvt(j))) > 1.e-7_dp )
       end do
       
       ! We still require that the offset in the T-direction is the same
       ! is this even necessary?
       er = er .or. ( abs(this_kdispl(this%t_dir) - kdispl(pvt(this%t_dir))) > 1.e-7_dp )
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
          this_kcell(:,1) = kcell(:,pvt(1)) * this%Rep(1)
          this_kcell(:,2) = kcell(:,pvt(2)) * this%Rep(2)
          this_kcell(:,3) = kcell(:,pvt(3)) * this%Rep(3)
          do j = 1 , 3
             write(*,'(3(i4,tr1),f8.4)') (this_kcell(i,j),i=1,3),kdispl(pvt(j))
          end do
       end if

    else
       
       call ts_read_TSHS_opt(this%HSfile, Gamma=Gamma, Bcast=.true.)

    end if

    ldie = ldie .or. er

#ifdef MPI
    call MPI_Bcast(ldie,1,MPI_Logical,0,MPI_Comm_World,MPIerror)
#endif

    if ( Gamma ) then
       write(*,*) 'Electrode : '//trim(this%name)//' is a Gamma-only calculation &
            &this is not feasible.'
       ldie = .true.
    end if

    if ( ldie ) then
       call die("The electrode does not conform with the system settings. &
            &Please correct accordingly.")
    end if

  end subroutine check_Elec

  subroutine check_connectivity(this)

    use parallel, only : IONode
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
    H => val(this%H)
    S => val(this%S)

    tm(:) = TM_ALL
    if ( this%inf_dir == INF_POSITIVE ) then
       tm(this%t_dir) =  2
    else
       tm(this%t_dir) = -2
    end if
    call crtSparsity_SC(this%sp,sp02, TM=tm, &
         ucell=this%cell, isc_off=this%isc_off)

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
       io = index_local_to_global(fdist,i)
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
    use intrinsic_missing, only : VNORM
    type(Elec), intent(in) :: this
    integer,  intent(in) :: na_u
    real(dp), intent(in) :: xa(3,na_u)

    ! Local variables
    logical :: do_print
    integer  :: i,j,k, ia, iaa
    real(dp), parameter :: check_xa = 0.0005_dp * Ang
    real(dp) :: xa_o(3), this_xa_o(3), ucell(3,3), tmp(3)
    real(dp), pointer :: this_xa(:,:)

    if ( .not. IONode ) return

    this_xa      => this%xa_used
    xa_o(:)      =  xa(:,this%idx_a)
    this_xa_o(:) =  this_xa(:,1)
    ucell        =  this%cell

    ! We only print out this structure if it does not fit the coordinates
    do_print = .false.
    
    iaa = this%idx_a
    do ia = 1 , this%na_used
       do k = 0 , this%Rep(3) - 1
       do j = 0 , this%Rep(2) - 1
       do i = 0 , this%Rep(1) - 1
          tmp(1) = this_xa(1,ia)-this_xa_o(1)+sum(ucell(1,:)*(/i,j,k/))
          tmp(2) = this_xa(2,ia)-this_xa_o(2)+sum(ucell(2,:)*(/i,j,k/))
          tmp(3) = this_xa(3,ia)-this_xa_o(3)+sum(ucell(3,:)*(/i,j,k/))
          do_print = do_print .or. VNORM(xa(:,iaa) - xa_o - tmp) > check_xa
          iaa = iaa + 1
       end do
       end do
       end do
    end do

    if ( .not. do_print ) return

    write(*,*) trim(this%name)//' unit cell (Ang):'
    write(*,'(2(3(tr1,f10.5),/),3(tr1,f10.5))') this%cell/Ang
    
    write(*,'(a,t35,a)') &
         " Structure of "//trim(this%name)//" electrode","| System electrode:"
    write(*,'(t3,3a10,''  |'',3a10,''  | '',a10)') &
         "X (Ang)","Y (Ang)","Z (Ang)", "X (Ang)","Y (Ang)","Z (Ang)","|r_S-r_E|"

    iaa = this%idx_a
    do ia = 1 , this%na_used
       do k = 0 , this%Rep(3) - 1
       do j = 0 , this%Rep(2) - 1
       do i = 0 , this%Rep(1) - 1
          tmp(1) = this_xa(1,ia)-this_xa_o(1)+sum(ucell(1,:)*(/i,j,k/))
          tmp(2) = this_xa(2,ia)-this_xa_o(2)+sum(ucell(2,:)*(/i,j,k/))
          tmp(3) = this_xa(3,ia)-this_xa_o(3)+sum(ucell(3,:)*(/i,j,k/))
          write(*,'(t3,3f10.5,''  |'',3f10.5,''  | '',e10.5)') &
               tmp / Ang, (xa(:,iaa) - xa_o(:))/Ang, &
               VNORM(xa(:,iaa) - xa_o - tmp) / Ang
          iaa = iaa + 1
       end do
       end do
       end do
    end do
    
    write(*,*) ! new-line

  end subroutine print_Elec

  subroutine copy_DM(this,na_u,xa,lasto,nsc,isc_off,cell,DM_2D, EDM_2D, &
       na_a, allowed)
    use m_handle_sparse
    use m_iodm
    use m_ts_iodm
    ! We will copy over the density matrix and "fix it in the leads
    type(Elec), intent(inout) :: this
    integer, intent(in) :: na_u, lasto(0:na_u), nsc(3), isc_off(3,product(nsc))
    real(dp), intent(in) :: xa(3,na_u), cell(3,3)
    type(dSpData2D), intent(inout) :: DM_2D, EDM_2D
    integer, intent(in) :: na_a, allowed(na_a)

    type(OrbitalDistribution) :: fake_dit
    type(Sparsity), pointer :: sp
    type(dSpData2D) :: f_DM_2D, f_EDM_2D
    real(dp), pointer :: DM(:,:), EDM(:,:)
    real(dp) :: tmp, Ef
    integer :: i
    logical :: found, alloc(3), is_TSDE

    ! Check and see if we should de-allocate
    alloc(1) = associated(this%xa)
    alloc(2) = associated(this%lasto)
    alloc(3) = associated(this%isc_off)

    ! We *must* read in the isc_off array
    call read_Elec(this, Bcast=.true., IO=.false.)

    i = len_trim(this%DEfile)
    is_TSDE = ( this%DEfile(i-3:i) == 'TSDE' )
    if ( is_TSDE ) then
       call read_ts_dm( this%DEfile, this%nspin, fake_dit, &
            this%no_u, f_DM_2D, f_EDM_2D, Ef, found, &
            Bcast = .true. )
    else
       call read_dm( this%DEfile, this%nspin, fake_dit, &
            this%no_u, f_DM_2D, found, &
            Bcast = .true. )
    end if
    if ( .not. found ) call die('Could not read file: '//trim(this%DEfile))
    sp => spar(f_DM_2D)
    if ( .not. equivalent(sp,this%sp) ) then
       call die('Bulk electrode expansion, read in sparsity pattern, &
            &does not match the TSHS sparsity pattern.')
    end if

    ! We must delete the additional arrays read in
    call delete(this%H)
    call delete(this%S)
    call delete(this%sp)

    ! TODO, there is some memory that could be leaking with the
    ! electrode arrays.

    if ( is_TSDE ) then
       ! Shift the energy matrix to the chemical potential :)
       DM  => val(f_DM_2D)
       EDM => val(f_EDM_2D)
       i = size(DM)
       tmp = -( Ef + this%mu%mu )
       call daxpy(i,tmp,DM(1,1),1,EDM(1,1),1)
    end if

    if ( this%inf_dir == INF_POSITIVE ) then
       i = 1
    else
       i = this%na_u - this%na_used + 1
    end if

    call expand_spd2spd_2D(i,this%na_used, &
         this%na_u,this%lasto,this%xa,f_DM_2D,&
         this%cell, (/1,1,1/), this%Rep, &
         product(this%nsc), this%isc_off, &
         na_u,xa,lasto,DM_2D,cell,product(nsc),isc_off, this%idx_a, &
         print = .true., allowed_a = allowed)

    if ( is_TSDE ) then
       call expand_spd2spd_2D(i,this%na_used, &
            this%na_u,this%lasto,this%xa,f_EDM_2D, &
            this%cell, (/1,1,1/), this%Rep, &
            product(this%nsc), this%isc_off, &
            na_u,xa,lasto,EDM_2D,cell,product(nsc),isc_off, this%idx_a, &
            allowed_a = allowed)
    end if
       
    if ( .not. alloc(1) ) deallocate(this%xa) ; nullify(this%xa)
    if ( .not. alloc(2) ) deallocate(this%lasto) ; nullify(this%lasto)
    if ( .not. alloc(3) ) deallocate(this%isc_off) ; nullify(this%isc_off)

    if ( is_TSDE ) then
       call delete(f_EDM_2D)
    end if
    call delete(f_DM_2D)

  end subroutine copy_DM

  subroutine print_settings(this,prefix,plane,box)
    use units, only : eV, Ang, Kelvin
    use parallel, only : Node
    type(Elec), intent(in) :: this
    character(len=*), intent(in) :: prefix
    logical, intent(in), optional :: plane, box

    character(len=100) :: chars
    character(len=60) :: f1, f5, f20, f6, f7, f8, f9, f10, f11, f15, f3, f16

    if ( Node /= 0 ) return
    
    ! Write out the settings
    ! First create the different out-put options
    f1  = '('''//trim(prefix)//': '',a,t53,''='',4x,l1)'
    f3  = '('''//trim(prefix)//': '',a,t53,''= { '',2(e12.5,'','',tr1),e12.5,''}'',a)'
    f5  = '('''//trim(prefix)//': '',a,t53,''='',i5,a)'
    f20 = '('''//trim(prefix)//': '',a,t53,''= '',i0,'' -- '',i0)'
    f6  = '('''//trim(prefix)//': '',a,t53,''='',f10.4,tr1,a)'
    f7  = '('''//trim(prefix)//': '',a,t53,''='',f12.6,tr1,a)'
    f8  = '('''//trim(prefix)//': '',a,t53,''='',f10.4)'
    f9  = '('''//trim(prefix)//': '',a,t53,''='',e12.4,tr1,a)'
    f10 = '('''//trim(prefix)//': '',a,t53,''='',4x,a)'
    f11 = '('''//trim(prefix)//': '',a)'
    f15 = '('''//trim(prefix)//': '',a,t53,''= '',2(i0,'' x ''),i0)'
    f16 = '('''//trim(prefix)//': '',a,t53,''= A'',2(i0,'', A''),i0)'


    write(*,f11) '>> '//trim(name(this))
    write(*,f16) '  Electrode cell pivoting: A1, A2, A3', this%pvt
    if ( this%out_of_core ) then
       write(*,f10) '  GF file', trim(this%GFfile)
       write(*,f1)  '  Reuse existing GF-file', this%ReUseGF
    else
       write(*,f11)  '  In-core GF'
    end if
    write(*,f10) '  Electrode TSHS file', trim(this%HSfile)
    write(*,f5)  '  # atoms used in electrode', this%na_used
    write(*,f15) '  Electrode repetition [A1 x A2 x A3]', this%Rep(:)
    write(*,f20) '  Position in geometry', this%idx_a, &
         this%idx_a + TotUsedAtoms(this) - 1
    if ( this%t_dir == 1 ) then
       chars = 'A1'
    else if ( this%t_dir == 2 ) then
       chars = 'A2'
    else if ( this%t_dir == 3 ) then
       chars = 'A3'
    end if
    if ( this%inf_dir == INF_POSITIVE ) then
       chars = 'positive wrt. '//trim(chars)
    else
       chars = 'negative wrt. '//trim(chars)
    end if
    write(*,f10) '  Semi-infinite direction for electrode', trim(chars)
    write(*,f7)  '  Chemical shift', this%mu%mu/eV,'eV'
    write(*,f7)  '  Electronic temperature', this%mu%kT/Kelvin,'K'
    write(*,f1)  '  Bulk values in electrode', this%Bulk
    if ( product(this%Rep) > 1 .and. this%out_of_core ) then
       if ( this%pre_expand == 0 ) then
          chars = 'none'
       else if ( this%pre_expand == 1 ) then
          chars = 'GS'
       else
          chars = 'H, S, GS'
       end if
       write(*,f10)  '  Pre-expansion to reduce computation', trim(chars)
    end if
#ifndef TBTRANS
    if ( this%DM_update == 0 ) then
       write(*,f11)  '  Cross-terms are not updated'
    else if ( this%DM_update == 1 ) then
       write(*,f11)  '  Cross-terms are updated'
    else if ( this%DM_update == 2 ) then
       write(*,f11)  '  Cross-terms and electrode region are updated'
    end if
#endif
    if ( abs(this%mu%mu) > 1.e-10_dp ) then
       write(*,f8)  '  Hamiltonian E-C Ef fractional shift', this%Ef_frac_CT
    end if
#ifndef TBTRANS
    if ( .not. this%kcell_check ) then
       write(*,f11)  '  Will NOT check the kgrid-cell! Ensure sampling!'
    end if
#endif
    write(*,f9)  '  Electrode imaginary Eta', this%Eta/eV,' eV'
#ifndef TBTRANS
    if ( present(plane) ) then
    if ( plane ) then
       write(*,f11) '  Hartree fix plane:'
       write(*,f3)  '    plane origo',this%p%c / Ang, ' Ang'
       write(*,f3)  '    plane normal vector',this%p%n
    end if
    end if
    if ( present(box) ) then
    if ( box ) then
       write(*,f11) '  Hartree potential box:'
       write(*,f3)  '    box origo',this%box%c / Ang, ' Ang'
       write(*,f3)  '    box v1',this%box%v(:,1) / Ang, ' Ang'
       write(*,f3)  '    box v2',this%box%v(:,2) / Ang, ' Ang'
       write(*,f3)  '    box v3',this%box%v(:,3) / Ang, ' Ang'
    end if
    end if
#endif

  end subroutine print_settings

  function Elec_idx(N_Elec,Elecs,El) result(idx)
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec), El
    integer :: idx
    do idx = 1 , N_Elec
       if ( Elecs(idx) == El ) return
    end do
    idx = 0
  end function Elec_idx
  
end module m_ts_electype
