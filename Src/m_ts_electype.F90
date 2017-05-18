module m_ts_electype

  use precision, only : dp

  use class_Sparsity
  use class_dSpData1D
  use class_dSpData2D
  use m_region

  use m_geom_plane, only: geo_plane_delta
  use m_geom_plane, only: in_basal_Elec => voxel_in_plane_delta

  use m_geom_box, only: geo_box_delta
  use m_geom_box, only: in_Elec => voxel_in_box_delta

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
  public :: q_exp, q_exp_all, Elec_kpt

  public :: fdf_nElec, fdf_elec

  public :: read_Elec
  public :: create_sp2sp01
  public :: print_Elec
  public :: print_settings
  public :: init_Elec_sim, check_Elec_sim
  public :: check_connectivity
  public :: delete

  public :: Elec_box2grididx, Elec_frac

  public :: in_basal_Elec, in_Elec

  public :: operator(.eq.)

  public :: copy_DM

  integer, parameter, public :: FILE_LEN = 256
  integer, parameter, public :: NAME_LEN = 32

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
     ! Bloch expansions (repetitions)
     integer :: Bloch(3) = 1
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
     ! The inter-layer distance between succeeding layers along
     ! the semi-infinite direction
     real(dp) :: dINF_layer
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

     ! The accuracy required for the self-energy
     !  == 1e-13 * eV
     real(dp) :: accu = 7.349806700083788e-15_dp

     ! The imaginary part in the electrode
     !  == 0.0001 * eV
     real(dp) :: Eta = 7.3498067e-6_dp

     ! The region of the down-folded region
     type(tRgn) :: o_inD, inDpvt

     ! The basal plane of the electrode
     type(geo_plane_delta) :: p
     ! A box containing all atoms of the electrode in the
     ! simulation box
     type(geo_box_delta) :: box

  end type Elec

  ! Fraction of alignment for considering two vectors having a similar
  ! component
  real(dp), parameter :: cell_unit_align = 1.e-5_dp

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
    
    ! If there are no electrodes, return immediately
    ! This will signal a new read, or default values of the
    ! electrodes
    if ( n == 0 ) return

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
    use parallel, only : IONode

    use fdf
    use intrinsic_missing, only: VNORM
    use m_os, only : file_exist
    use m_ts_io, only : ts_read_TSHS_opt
    use m_ts_io_ctype, only : pline_E_parse

    character(len=*), intent(in) :: prefix,slabel
    type(Elec), intent(inout) :: this
    integer, intent(in) :: N_mu
    type(ts_mu), intent(in), target :: mus(N_mu)
    integer, intent(in), optional :: idx_a
    character(len=*), intent(in), optional :: name_prefix

    logical :: found, mix_bloch, bloch_rep

    ! prepare to read in the data...
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline => null()
    logical :: info(5)
    integer :: i, j
    integer :: cidx_a 
    real(dp) :: rcell(3,3), fmin, fmax, rc

    character(len=FILE_LEN) :: bName, name, ln, tmp

    info(:) = .false.

    bName = trim(prefix)//'.Elec.'//trim(this%name)

    ! Allow the filename to be read in individually
    name = trim(bName)//'.HS'
    if ( fdf_defined(trim(name)) ) then
       this%HSfile = trim(fdf_get(trim(name),''))
       info(1) = .true.
    end if

    do i = 1 , 3
       write(name,'(2a,i0)') trim(bName),'.Bloch.A',i
       this%Bloch(i) = fdf_get(trim(name), 1)
    end do
    name = trim(bName)//'.GF'
    if ( fdf_defined(trim(name)) ) this%GFfile = trim(fdf_get(trim(name),''))

    ! Default to use the chemical potential with the same
    ! name as the electrode
    ! If the block specifies the chemical
    ! potential, it will be the preferred method
    ln = trim(this%name)
    nullify(this%mu)
    do i = 1 , N_mu
       if ( leqi(ln,mus(i)%name) ) then
          this%mu => mus(i)
          info(3) = .true.
          exit
       end if
    end do
    ! If there is only one chemical potential
    ! then, of course, they are associated.
    if ( N_mu == 1 ) then
       this%mu => mus(1)
       info(3) = .true.
    end if

    
    ! Figure out if we should return immediately
    found = fdf_block(trim(bName),bfdf)
    if ( .not. found ) return

    
    cidx_a = 0
    ! Denote that no bloch expansion coefficients has
    ! been set.
    mix_bloch = .false.
    bloch_rep = .false.

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
       if ( leqi(ln,'HS') .or. leqi(ln,'HS-file') .or. &
            leqi(ln,'TSHS') .or. &
            leqi(ln,'TSHS-file') ) then
          if ( fdf_bnnames(pline) < 2 ) call die('HS name not supplied')
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
                  '" for the electrode: "'//trim(this%name)//'". '//&
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

       else if ( leqi(ln,'bloch') ) then
          if ( fdf_bnintegers(pline) < 3 ) &
               call die('Bloch expansion for all directions are not supplied <A1> <A2> <A3>')
          this%Bloch(1) = fdf_bintegers(pline,1)
          this%Bloch(2) = fdf_bintegers(pline,2)
          this%Bloch(3) = fdf_bintegers(pline,3)
          mix_bloch = .true.
          
       else if ( leqi(ln,'bloch-a') .or. leqi(ln,'bloch-a1') ) then
          if ( fdf_bnintegers(pline) < 1 ) &
               call die('Bloch expansion of A1 is not supplied')
          this%Bloch(1) = fdf_bintegers(pline,1)
          mix_bloch = .true.
          
       else if ( leqi(ln,'bloch-b') .or. leqi(ln,'bloch-a2') ) then
          if ( fdf_bnintegers(pline) < 1 ) &
               call die('Bloch expansion of A2 is not supplied')
          this%Bloch(2) = fdf_bintegers(pline,1)
          mix_bloch = .true.
          
       else if ( leqi(ln,'bloch-c') .or. leqi(ln,'bloch-a3') ) then
          if ( fdf_bnintegers(pline) < 1 ) &
               call die('Bloch expansion of A3 is not supplied')
          this%Bloch(3) = fdf_bintegers(pline,1)
          mix_bloch = .true.

       else if ( leqi(ln,'replicate-a') .or. leqi(ln,'rep-a') .or. &
            leqi(ln,'replicate-a1') .or. leqi(ln,'rep-a1') ) then
          if ( fdf_bnintegers(pline) < 1 ) &
               call die('Bloch expansion of A1 is not supplied')
          this%Bloch(1) = fdf_bintegers(pline,1)
          if ( mix_bloch ) then
             call die('A "Bloch" keyword was found previously. &
                  &No mixing of rep/Bloch may be performed.')
          end if
          bloch_rep = .true.

       else if ( leqi(ln,'replicate-b') .or. leqi(ln,'rep-b') .or. &
            leqi(ln,'replicate-a2') .or. leqi(ln,'rep-a2') ) then
          if ( fdf_bnintegers(pline) < 1 ) &
               call die('Bloch expansion of A2 is not supplied')
          this%Bloch(2) = fdf_bintegers(pline,1)
          if ( mix_bloch ) then
             call die('A "Bloch" keyword was found previously. &
                  &No mixing of rep/Bloch may be performed.')
          end if
          bloch_rep = .true.

       else if ( leqi(ln,'replicate-c') .or. leqi(ln,'rep-c') .or. &
            leqi(ln,'replicate-a3') .or. leqi(ln,'rep-a3') ) then
          if ( fdf_bnintegers(pline) < 1 ) &
               call die('Bloch expansion of A3 is not supplied')
          this%Bloch(3) = fdf_bintegers(pline,1)
          if ( mix_bloch ) then
             call die('A "Bloch" keyword was found previously. &
                  &No mixing of rep/Bloch may be performed.')
          end if
          bloch_rep = .true.
          
       else if ( leqi(ln,'replicate') .or. leqi(ln,'rep') ) then
          if ( fdf_bnintegers(pline) < 3 ) &
               call die('Bloch expansion for all directions are not supplied <A1> <A2> <A3>')
          this%Bloch(1) = fdf_bintegers(pline,1)
          this%Bloch(2) = fdf_bintegers(pline,2)
          this%Bloch(3) = fdf_bintegers(pline,3)
          if ( mix_bloch ) then
             call die('A "Bloch" keyword was found previously. &
                  &No mixing of rep/Bloch may be performed.')
          end if
          bloch_rep = .true.

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

       else if ( leqi(ln,'DE') .or. leqi(ln,'DE-file') .or. &
            leqi(ln,'TSDE') .or. leqi(ln,'TSDE-file') ) then
          if ( fdf_bnnames(pline) < 2 ) call die('DE name not supplied')
          this%DEfile = trim(fdf_bnames(pline,2))

#ifdef TBTRANS
       else if ( leqi(ln,'tbt.Accuracy') .or. leqi(ln,'Accuracy') ) then
#else
       else if ( leqi(ln,'Accuracy') ) then
#endif
          call pline_E_parse(pline,1,ln, &
               val = this%accu, before=3)

#ifdef TBTRANS
       else if ( leqi(ln,'tbt.Eta') .or. leqi(ln,'Eta') ) then
#else
       else if ( leqi(ln,'Eta') ) then
#endif
          call pline_E_parse(pline,1,ln, &
               val = this%Eta, before=3)
#ifdef TBTRANS
#ifdef TBT_PHONON
          ! eta value needs to be squared as it is phonon spectrum
          this%Eta = this%Eta ** 2
#endif
#endif

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
               &for electrode: '//trim(this%name))

       end if

    end do
    
    if ( any(this%Bloch(:) < 1) ) then
       call die("Bloch expansion in "//trim(this%name)//" electrode must be >= 1.")
    end if

    if ( .not. all(info(1:4)) ) then
       write(*,*)'You need to supply at least:'
       write(*,*)' - HS'
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

    ! If the user will not use bulk, and haven't set DM-update,
    ! default to 'all'
    if ( .not. info(5) .and. .not. this%Bulk ) then
       this%DM_update = 2 ! set 'all'
    end if
    
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

    ! Calculate the interlayer distance between succeeding layers
    ! along the semi-infinite direction
    fmin =  huge(1._dp)
    fmax = -huge(1._dp)
    call reclat(this%cell,rcell,0)
    do i = 1 , this%na_u
       rc = sum(this%xa(:,i) * rcell(:,this%t_dir))
       fmin = min(fmin,rc)
       fmax = max(fmax,rc)
    end do
    ! Get distance to the cell boundary
    rc = 1._dp - (fmax-fmin)
    ! Calculate the inter-layer distance
    this%dINF_layer = VNORM( rc * this%cell(:,this%t_dir) )

    ! We deallocate xa and lasto as they are not needed
    deallocate(this%xa,this%lasto)

    ! Check that the Bloch expansion is not in the transport-direction
    if ( this%Bloch(this%t_dir) /= 1 ) then
       call die('Bloch expansion in the transport direction &
            &is not allowed.')
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

    ! Print out the error if using the repetition keyword
    if ( IONode .and. bloch_rep ) then
       write(*,'(/a)')'DEPRECATION WARNING:'
       write(*,'(5(a,/))')'Electrode Bloch expansion keyword in electrode '&
            //trim(this%name)//':', &
            '  replicate','has been superseeded by:', &
            '  bloch', 'The replicate keyword may be ignored in future versions.'
    end if

  end function fdf_Elec

  ! Initialize variables for the electrode according
  ! to the simulation variables
  subroutine init_Elec_sim(this,cell,na_u,xa)

    use units, only : Pi
    use intrinsic_missing, only : VNORM, SPC_PROJ, IDX_SPC_PROJ
    use intrinsic_missing, only : VEC_PROJ_SCA

    ! The electrode that needs to be processed
    type(Elec), intent(inout) :: this
    ! The simulation parameters.
    ! Unit-cell of simulation cell
    real(dp), intent(in) :: cell(3,3)
    integer, intent(in) :: na_u
    real(dp), intent(in) :: xa(3,na_u)
    
    real(dp) :: p(3), contrib
    integer :: i, j, ia, na

    ! First figure out the minimum bond-length
    ia = this%idx_a
    na = TotUsedAtoms(this)

    ! Calculate the pivoting table
    do i = 1 , 3

       ! We just want to find the cell vector
       ! which is closests to 1, that ensures a parallel
       ! cell vector
       ! Note that here we do not enforce the direction
       ! of the cell vector.
       ! This is because it _can_ be opposite for directions
       ! without k-point samplings.
       p = SPC_PROJ(cell, this%cell(:,i))
       p = p / VNORM(this%cell(:,i))
       if ( abs(abs(p(1)) - 1._dp) < cell_unit_align ) then
          this%pvt(i) = 1
       else if ( abs(abs(p(2)) - 1._dp) < cell_unit_align ) then
          this%pvt(i) = 2
       else if ( abs(abs(p(3)) - 1._dp) < cell_unit_align ) then
          this%pvt(i) = 3
       else
          ! We simply take the largest one
          this%pvt(i) = IDX_SPC_PROJ(cell,this%cell(:,i), mag=.true.)
       end if
       
    end do
    if ( sum(this%pvt) /= 6 .or. count(this%pvt==2) /= 1 ) then
       print *, this%pvt
       call die('The pivoting table for the electrode unit-cell, &
            &onto the simulation unit-cell is not unique. &
            &Please check your simulation cells.')
    end if
    
    ! Create the basal plane of the electrode
    ! Decide which end of the electrode we use
    ! Calculate planes of the electrodes
    p = this%cell(:,this%t_dir)
    p = p / VNORM(p)

    ! Select the atom farthest from the device region
    ! along the semi-infinite direction
    j = ia
    contrib = VEC_PROJ_SCA(p,xa(:,j))
    if ( this%inf_dir == INF_POSITIVE ) then
       ! We need to utilize the last atom 
       do i = ia + 1, ia + na - 1
          if ( VEC_PROJ_SCA(p,xa(:,i)) > contrib ) then
             j = i
             contrib = VEC_PROJ_SCA(p,xa(:,j))
          end if
       end do
    else
       ! We need to utilize the first atom
       do i = ia + 1, ia + na - 1
          if ( VEC_PROJ_SCA(p,xa(:,i)) < contrib ) then
             j = i
             contrib = VEC_PROJ_SCA(p,xa(:,j))
          end if
       end do
    end if
    this%p%c = xa(:,j)

    ! We add a vector with length of half the minimal bond length
    ! to the vector, to do the averaging 
    ! not on-top of an electrode atom.
    contrib = this%dINF_layer * 0.5_dp
    if ( this%inf_dir == INF_POSITIVE ) then
       this%p%c = this%p%c + p * contrib ! add vector
    else
       this%p%c = this%p%c - p * contrib ! subtract vector
    end if

    ! Normal vector to electrode basal plane.
    ! This coincides with the electrode semi-infinite direction
    this%p%n = p ! normalized semi-infinite direction vector

    ! The distance parameter used when calculating
    ! the +/- voxel placements
    this%p%d = sum( this%p%n(:) * this%p%c(:) )

    ! *** Now we have created the Hartree plane
    !     where the potential is fixed.         ***


    ! *** Calculate the electrode box for N-electrodes
    ! We do this by:
    !  1. find the average Cartesian coordinate.
    !  2. subtract half the unit-cell vector of the
    !     electrode.
    !  3. Create the box from the electrode unit-cell

    do i = 1 , 3

       ! Calculate the Cartesian ith contribution
       contrib = sum(this%cell(i,:)) * 0.5_dp

       ! Calculate average position
       ! and find the lower left corner of the electrode
       p(i) = sum(xa(i,ia:ia+na-1)) / na
       
       ! transfer to the corner of the Hartree box
       this%box%c(i) = p(i) - contrib

    end do

    ! The box is the same as the electrode cell
    this%box%v = this%cell
    
    ! If we do not use all the atoms we need to reduce
    ! the cell vector in the semi-infinite direction
    if ( this%na_u /= this%na_used ) then
       contrib = real(this%na_used,dp) / real(this%na_u,dp)
       i = this%t_dir
       ! Scale vector
       this%box%v(:,i) = this%box%v(:,i) * contrib
    end if

  end subroutine init_Elec_sim

  ! Initialize variables for the electrode according
  ! to the simulation variables
  subroutine check_Elec_sim(this,nspin,s_cell,na_u,xa,xa_EPS,&
       lasto, Gamma3, &
       kcell, kdispl )

    use parallel, only : IONode
    use units, only : Pi, Ang
    use intrinsic_missing, only : VNORM

    use m_ts_io, only: ts_read_TSHS_opt

    ! The electrode that needs to be processed
    type(Elec), intent(inout) :: this
    ! The number of spin components from the siesta calculation
    integer, intent(in) :: nspin
    ! The simulation parameters.
    ! Unit-cell of simulation cell
    real(dp), intent(in) :: s_cell(3,3)
    integer, intent(in) :: na_u
    real(dp), intent(in) :: xa(3,na_u), xa_EPS ! epsilon for space check
    ! The last orbital of each atom
    integer, intent(in) :: lasto(0:na_u)
    ! Whether we check for Gamma point
    ! If it is Gamma we do not have as strict requirements
    logical, intent(in) :: Gamma3(3)
    ! The transiesta kcell, kdispl
    integer, intent(in), optional :: kcell(3,3)
    real(dp), intent(in), optional :: kdispl(3)
    
    ! Local variables
    integer :: this_kcell(3,3)
    real(dp) :: xa_o(3), this_xa_o(3), cell(3,3), this_kdispl(3)
    real(dp) :: max_xa(3), cur_xa(3)
    real(dp), pointer :: this_xa(:,:)
    integer :: i, j, k, ia, na, pvt(3), iaa, idir(3)
    logical :: ldie, er, Gamma
    
    na = TotUsedAtoms(this)
    if ( na == 1 ) then
       call die('One atom electrodes are not allowed')
    end if

    ldie = .false.

    ! Copy
    cell = this%cell
    pvt = this%pvt
    this_xa => this%xa_used
    xa_o(:) = xa(:,this%idx_a)
    this_xa_o(:) = this_xa(:,1)

    max_xa = 0._dp
    er = .false.
    iaa = this%idx_a
    do ia = 1 , this%na_used
       
       do k = 1 , this%Bloch(3)
       idir(3) = k - 1
       do j = 1 , this%Bloch(2)
       idir(2) = j - 1
       do i = 1 , this%Bloch(1)
          idir(1) = i - 1

          ! Calculate repetition vector
          cur_xa(1) = sum(cell(1,:)*idir)
          cur_xa(2) = sum(cell(2,:)*idir)
          cur_xa(3) = sum(cell(3,:)*idir)
          
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

    if ( er .and. IONode ) then
       iaa = this%idx_a

       ! We need to write out the correct formatting of the electrodes.
       ! We shift to the first electrode coordinate, and add the 
       ! first system electrode atom. 
       ! If the user has the %block AtomicCoord... in:
       !    LatticeConstant         1. Ang
       !    AtomicCoordinatesFormat    Ang
       ! then it is direct copy paste ! :)
       xa_o(:) = -this_xa(:,1) + xa(:,iaa)

       i = iaa
       j = iaa + TotUsedAtoms(this) - 1
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
       do ia = 1 , this%na_used
          do k = 0 , this%Bloch(3)-1
          do j = 0 , this%Bloch(2)-1
          do i = 0 , this%Bloch(1)-1
             write(*,'(t2,3(tr1,f20.10))') &
                  (this_xa(1,ia)+xa_o(1)+sum(cell(1,:)*(/i,j,k/)))/Ang, &
                  (this_xa(2,ia)+xa_o(2)+sum(cell(2,:)*(/i,j,k/)))/Ang, &
                  (this_xa(3,ia)+xa_o(3)+sum(cell(3,:)*(/i,j,k/)))/Ang
          end do
          end do
          end do
       end do

    end if

    iaa = this%idx_a
    er = .false.
    do ia = 1 , this%na_used
       do k = 1 , this%Bloch(3)
       do j = 1 , this%Bloch(2)
       do i = 1 , this%Bloch(1)

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

    if ( er .and. IONode) then
       write(*,'(a)') "Number of orbitals per atom in the electrode does not match the system electrode"
       write(*,'(a)') 'Have you changed your basis size?'
       write(*,'(t3,3a20)') "ia system","n_orb_el","n_orb_sys"
       iaa = this%idx_a
       do ia = 1 , this%na_used
          do k = 1 , this%Bloch(3)
          do j = 1 , this%Bloch(2)
          do i = 1 , this%Bloch(1)
                
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

    if ( nspin /= this%nspin ) then
       write(*,'(a)') 'ERROR: Electrode: '//trim(this%name)
       write(*,'(a,i0)') 'Electrode, nspin = ',this%nspin
       write(*,'(a,i0)') 'transiesta, nspin = ',nspin
       ldie = .true.
    end if

    er = .false.
    if ( present(kcell) .and. this%kcell_check ) then
       
       call ts_read_TSHS_opt(this%HSfile, &
            kscell=this_kcell,kdispl=this_kdispl, &
            Gamma=Gamma, Bcast=.true.)

       ! Check that there is actually k-points in the transport direction
       j = this%t_dir
       i = this_kcell(j,j)
       if ( i < 20 .and. IONode ) then
          write(*,'(a)') 'Electrode: '//trim(this%name)//' has very few &
               &k-points in the semi-infinite direction, at least 20 is recommended.'
       else if ( i < 5 .and. IONode ) then
          write(*,'(a)') 'Electrode: '//trim(this%name)//' has exceptionally few &
               &k-points in the semi-infinite direction, at least 5 is required.'
          ldie = .true.
       end if

       ! If the system is not a Gamma calculation, then the file must
       ! not be either (the Bloch expansion will only increase the number of
       ! k-points, hence the above)
       do j = 1 , 3
          k = this%Bloch(j)
          ! The displacements are not allowed non-equivalent.
          er = er .or. ( abs(this_kdispl(j) - kdispl(pvt(j))) > 1.e-7_dp )
          if ( j == this%t_dir ) cycle
          do i = 1 , 3
             if ( i == this%t_dir ) cycle
             if ( j == i ) then
                er = er .or. ( this_kcell(i,j) /= kcell(pvt(i),pvt(j))*k )
             else 
                er = er .or. ( this_kcell(i,j) /= kcell(pvt(i),pvt(j)) )
             end if
          end do
       end do

       if ( er .and. IONode ) then
          
          write(*,'(a)') 'Incompatible k-grids...'
          write(*,'(a)') 'Electrode file k-grid:'
          do j = 1 , 3
             write(*,'(3(i4,tr1),f8.4)') this_kcell(:,j), this_kdispl(j)
          end do
          write(*,'(a)') 'System k-grid:'
          do j = 1 , 3
             write(*,'(3(i4,tr1),f8.4)') kcell(:,j), kdispl(j)
          end do
          write(*,'(a)') 'Electrode file k-grid should probably be:'
          ! Loop the electrode directions
          do j = 1 , 3
             if ( j == this%t_dir ) then
                ! ensure that we retain the semi-infinite
                ! direction k-sampling, and suggest more k-points
                ! if necessary
                if ( this_kcell(j,j) < 20 ) then
                   this_kcell(j,j) = 50
                end if
             else
                this_kcell(:,j) = kcell(:,pvt(j)) * this%Bloch(j)
             end if
             write(*,'(3(i4,tr1),f8.4)') this_kcell(:,j), kdispl(pvt(j))
          end do
          
       end if

    else
       
       call ts_read_TSHS_opt(this%HSfile, Gamma=Gamma, Bcast=.true.)

    end if

    ldie = ldie .or. er

    if ( Gamma ) then
       write(*,'(a)') 'Electrode: '//trim(this%name)//' is a Gamma-only &
            &calculation, this is not feasible.'
       ldie = .true.
    end if

    if ( sum(this%pvt) /= 6 .or. count(this%pvt==2) /= 1 ) then
       write(*,'(tr2,i2)') this%pvt
       ldie = .true.
       write(*,'(a)')'The pivoting table for the electrode unit-cell, &
            &onto the simulation unit-cell is not unique. &
            &Please check your simulation cells.'
    end if
    
    if ( ldie ) then
       call die('Erroneous electrode setup, check out-put')
    end if
    
  end subroutine check_Elec_sim

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

  elemental function Name_(this) result(name)
    type(Elec), intent(in) :: this
    character(len=NAME_LEN) :: Name
    name = this%name
  end function Name_

  elemental function TotUsedAtoms(this) result(val)
    type(Elec), intent(in) :: this
    integer :: val
    val = this%na_used * product(this%Bloch)
  end function TotUsedAtoms

  pure function q_exp_all(this,i,j,k) result(q)
    type(Elec), intent(in) :: this
    integer, intent(in) :: i,j,k
    real(dp) :: q(3)
    
    ! TODO, the current implementation assumes k-symmetry
    ! of the electrode electronic structure.
    ! Using Bloch expansion with non-symmetry will, likely, produce
    ! wrong results.
    ! Luckily this is not a problem currently.
    ! Perhaps one should consider this in tbtrans

    q(1) = 1._dp*(i-1) / real(this%Bloch(1),dp)
    q(2) = 1._dp*(j-1) / real(this%Bloch(2),dp)
    q(3) = 1._dp*(k-1) / real(this%Bloch(3),dp)
    
  end function q_exp_all

  pure function q_exp(this,idx) result(q)
    type(Elec), intent(in) :: this
    integer, intent(in) :: idx
    real(dp) :: q(3)
    integer :: i,j,k,ii
    i =     this%Bloch(1)
    j = i * this%Bloch(2)
    k = j * this%Bloch(3)
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

  subroutine Elec_kpt(this,cell,k1,k2,opt)
    type(Elec), intent(in) :: this
    ! Device unit-cell
    real(dp), intent(in) :: cell(3,3)
    ! Input k-point
    real(dp), intent(in) :: k1(3)
    ! Output k-point
    real(dp), intent(out) :: k2(3)
    ! Whether or not this is an electrode 'output' k-point.
    !  opt ==  1, from device k-point [1/Bohr] to electrode k-point [1/Bohr]
    !  opt ==  2, from device k-point [1/Bohr] to electrode k-point [1/b]
    !  opt == -1, from electrode k-point [1/Bohr] to device k-point [1/Bohr]
    !  opt == -2, from electrode k-point [1/Bohr] to device  k-point [1/b]
    integer, intent(in), optional :: opt
    
    ! Local variables
    integer :: iop
    real(dp) :: tmp(3)

    iop = 1
    if ( present(opt) ) iop = opt

    select case ( iop ) 
    case ( 1 , 2 )
       
       ! Convert system-unit-cell kpoint to reciprocal units
       call kpoint_convert(cell,k1,k2,1)
       ! Scale with Bloch expansion
       tmp(1) = k2(this%pvt(1)) / this%Bloch(1)
       tmp(2) = k2(this%pvt(2)) / this%Bloch(2)
       tmp(3) = k2(this%pvt(3)) / this%Bloch(3)
       ! Remove semi-infinite direction
       tmp(this%t_dir) = 0._dp

       if ( iop == 1 ) then
          ! Convert back to 1 / Bohr
          call kpoint_convert(this%cell,tmp,k2,-1)
       else
          k2 = tmp
       end if

    case ( -2 , -1 )

       ! Convert from electrode to device
       ! Convert system-unit-cell kpoint to reciprocal units
       call kpoint_convert(this%cell,k1,k2,1)
       ! Scale with Bloch expansion
       tmp(this%pvt(1)) = k2(1) * this%Bloch(1)
       tmp(this%pvt(2)) = k2(2) * this%Bloch(2)
       tmp(this%pvt(3)) = k2(3) * this%Bloch(3)

       if ( iop == -1 ) then
          ! Convert back to 1 / Bohr
          call kpoint_convert(cell,tmp,k2,-1)
       else
          k2 = tmp
       end if

    end select

  end subroutine Elec_kpt

  elemental function TotUsedOrbs(this) result(val)
    type(Elec), intent(in) :: this
    integer :: val
    val = this%no_used * product(this%Bloch)
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

  subroutine Elec_box2grididx(this,mesh,dL,imin,imax)
    ! Returns grid indices for the minimum and maximum of
    ! the electrode box
    type(Elec), intent(in) :: this
    ! The global mesh size
    integer, intent(in) :: mesh(3)
    ! The mesh stepping per increment
    real(dp), intent(in) :: dL(3,3)
    ! The minimum/maximum grid indices
    integer, intent(out) :: imin(3), imax(3)

    real(dp) :: LHS(3,3), RHS(3)
    integer :: idx(3), i

    ! Initialize the indices
    imin = huge(1)
    imax = -huge(1)

    call get_idx(0,0,0,imin,imax)
    call get_idx(1,0,0,imin,imax)
    call get_idx(0,1,0,imin,imax)
    call get_idx(0,0,1,imin,imax)
    call get_idx(1,1,0,imin,imax)
    call get_idx(1,0,1,imin,imax)
    call get_idx(0,1,1,imin,imax)
    call get_idx(1,1,1,imin,imax)

    ! This pre-step will move them both simultaneously
    ! as that corresponds to equal shifts and no crossing
    ! of cell-boundaries
    do i = 1 , 3
       do while ( imin(i) <= 0 .and. imax(i) <= 0 )
          imin(i) = imin(i) + mesh(i)
          imax(i) = imax(i) + mesh(i)
       end do
       do while ( mesh(i) < imin(i) .and. mesh(i) < imax(i) )
          imin(i) = imin(i) - mesh(i)
          imax(i) = imax(i) - mesh(i)
       end do
    end do

  contains
    
    subroutine get_idx(ix,iy,iz,imin,imax)
      integer, intent(in) :: ix,iy,iz
      integer, intent(inout) :: imin(3), imax(3)
      
      ! Copy the LHS
      LHS = dL
      ! Create RHS
      RHS = this%box%c + this%box%v(:,1) * ix &
           + this%box%v(:,2) * iy &
           + this%box%v(:,3) * iz
      
      ! Calculate pqosition in the grid
      call dgesv(3,1,LHS,3,idx,RHS,3,i)
      if ( i /= 0 ) then
         call die('ts_voltage: Could not solve linear system.')
      end if
      
      ! Convert to integer position
      idx = nint(RHS)
      
      imin(1) = min(imin(1),idx(1))
      imin(2) = min(imin(2),idx(2))
      imin(3) = min(imin(3),idx(3))
      imax(1) = max(imax(1),idx(1))
      imax(2) = max(imax(2),idx(2))
      imax(3) = max(imax(3),idx(3))
      
    end subroutine get_idx
    
  end subroutine Elec_box2grididx

  ! Returns the fractional coordinates of the
  ! electrode atoms in the unit-cell.
  ! It can do so along any of the three cell
  ! directions.
  subroutine Elec_frac(this,cell,na_u,xa,dir, &
       fmin, fmax)
    type(Elec), intent(in) :: this
    real(dp), intent(in) :: cell(3,3)
    integer, intent(in) :: na_u
    real(dp), intent(in) :: xa(3,na_u)
    ! The cell vector that we will project onto
    integer, intent(in) :: dir
    ! Results (optional)
    real(dp), intent(out), optional :: fmin, fmax

    ! Local variables
    real(dp) :: lfmin, lfmax, rcell(3,3), r
    integer :: i

    call reclat(cell,rcell,0) ! without 2pi
    lfmin =  huge(1._dp)
    lfmax = -huge(1._dp)
    do i = this%idx_a , this%idx_a + TotUsedAtoms(this) - 1
       ! Get supercell in the transport direction
       r = sum(xa(:,i) * rcell(:,dir))
       lfmin = min(lfmin,r)
       lfmax = max(lfmax,r)
    end do

    if ( present(fmin) ) fmin = lfmin
    if ( present(fmax) ) fmax = lfmax

  end subroutine Elec_frac
       

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

    character(len=FILE_LEN) :: fN
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

    ! Check that there is a transfer matrix!
    if ( nnzs(this%sp01) == 0 ) then
       write(*,'(a)') 'Electrode '//trim(this%name)//' has no transfer matrix.'
       write(*,'(a)') 'The self-energy cannot be calculated with a zero transfer matrix!'
       call die('Elec: transfer matrix has 0 elements. The self-&
            &energy cannot be calculated. Please check your electrode &
            &electronic structure.')
    end if

  end subroutine create_sp2sp01

  subroutine delete_(this)
    type(Elec), intent(inout) :: this

    ! Full matrices
    call delete(this%H)
    call delete(this%S)
    call delete(this%sp)

    ! 00 matrices
    call delete(this%H00)
    call delete(this%S00)
    call delete(this%sp00)

    ! 01 matrices
    call delete(this%H01)
    call delete(this%S01)
    call delete(this%sp01)

    if ( associated(this%xa) ) deallocate(this%xa)
    if ( associated(this%lasto) ) deallocate(this%lasto)
    nullify(this%xa,this%lasto)
    if ( associated(this%isc_off) ) deallocate(this%isc_off)
    !if ( associated(this%xa_used) ) deallocate(this%xa_used)
    !if ( associated(this%lasto_used) ) deallocate(this%lasto_used)
    !nullify(this%xa_used,this%lasto_used)

  end subroutine delete_

  
  function check_connectivity(this) result(good)

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
    logical :: good

    real(dp), pointer :: H(:,:)
    real(dp), pointer :: S(:)
    type(OrbitalDistribution), pointer :: fdist
    type(Sparsity) :: sp02

    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer, pointer :: ncol01(:), ptr01(:), col01(:)
    integer, pointer :: ncol02(:), ptr02(:), col02(:)

    integer :: no_l, no_u, i, io, j, iu, ind, ind01, ind02, ia
    integer :: n_nzs
    integer :: istart, iend
    integer :: tm(3)
    ! Print-out values stored...
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
    n_nzs = nnzs(sp02)

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

    ! clean-up
    call delete(sp02)

    ! Secondly, we should also check for cases where the user requests
    ! num-atoms, which is the equivalent of creating a smaller electrode.
    if ( this%na_used /= this%na_u ) then

       call attach(this%sp01,n_col=ncol01,list_ptr=ptr01,list_col=col01)
       ! To assert the connectivity the removed atoms
       ! must not connect with them-selves
       if ( this%inf_dir == INF_NEGATIVE ) then
          ! We select the last atoms!
          ! so opposite here...
          istart = 1 
          iend = this%lasto(this%na_u-this%na_used)
       else
          ! We select the first atoms!
          ! so opposite here...
          istart = this%lasto(this%na_used) + 1
          iend = no_l
       end if
       
       ! loop and assign data elements
       do i = istart , iend
          
          ! Shift out of the buffer region
          io = index_local_to_global(fdist,i)
          ia = iaorb(io,this%lasto)
          
          ! Loop number of entries in the row...
          do j = 1 , ncol01(i)
             
             ! The index in the pointer array is retrieved
             ind01 = ptr01(i) + j
             
             ! Loop in the super-set sparsity pattern
             idx01: do ind = l_ptr(i) + 1 , l_ptr(i) + l_ncol(i)
                
                ! If we have the same column index it must be
                ! the same entry they represent
                iu = ucorb(col01(ind01), no_l)
                if ( col01(ind01) == l_col(ind) .and. &
                     istart <= iu .and. iu <= iend ) then

                   ! increment number of elements extending..
                   n_nzs = n_nzs + 1
                   
                   if ( any(abs(H(ind,:)) > maxH) ) then
                      maxH  = maxval(abs(H(ind,:)))
                   maxS  = S(ind)
                   maxi  = io
                   maxia = ia
                   maxj  = col01(ind01)
                   maxja = iaorb(col01(ind01),this%lasto)
                end if
                exit idx01
             end if

          end do idx01
          
       end do
       
    end do
       

    end if


    ! If both number
    good = n_nzs == 0
    if ( .not. good ) good = maxi == 0

    if ( .not. IONode ) return

    if ( n_nzs == 0 ) then
       write(*,'(t2,a)') trim(this%name)//' principal cell is perfect!'
    else if ( maxi == 0 ) then
       write(*,'(t2,a,i0,a)') trim(this%name)//' principal cell is extending out &
            &with ',n_nzs,' elements, all being zero.'
    else
       write(*,'(t2,a,i0,a)') trim(this%name)//' principal cell is extending out with ',n_nzs,' elements:'
       write(*,'(t5,2(a,i0))') 'Atom ',maxia,' connects with atom ',maxja
       write(*,'(t5,2(a,i0))') 'Orbital ',UCORB(maxi,no_u) , &
            ' connects with orbital ',UCORB(maxj,no_u)
       write(*,'(t5,3(a,i0),a,g10.3,a)') 'Hamiltonian value: |H(',&
            maxi,',',maxj,')|@R=',tm(this%t_dir),' = ',maxH/eV,' eV'
       write(*,'(t5,3(a,i0),a,g10.3)') 'Overlap          :  S(',&
            maxi,',',maxj,')|@R=',tm(this%t_dir),' = ',maxS
    end if

  end function check_connectivity

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
       do k = 0 , this%Bloch(3) - 1
       do j = 0 , this%Bloch(2) - 1
       do i = 0 , this%Bloch(1) - 1
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
       do k = 0 , this%Bloch(3) - 1
       do j = 0 , this%Bloch(2) - 1
       do i = 0 , this%Bloch(1) - 1
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
    integer, parameter :: One3(3) = (/1,1,1/)
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
       call read_ts_dm( this%DEfile, fake_dit, &
            f_DM_2D, f_EDM_2D, Ef, found, &
            Bcast = .true. )
    else
       call read_dm( this%DEfile, fake_dit, f_DM_2D, found, &
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
         this%cell, One3, this%Bloch, &
         product(this%nsc), this%isc_off, &
         na_u,xa,lasto,DM_2D,cell,product(nsc),isc_off, this%idx_a, &
         print = .true., allowed_a = allowed)

    if ( is_TSDE ) then
       call expand_spd2spd_2D(i,this%na_used, &
            this%na_u,this%lasto,this%xa,f_EDM_2D, &
            this%cell, One3, this%Bloch, &
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
    f10 = '('''//trim(prefix)//': '',a,t53,''='',tr1,a)'
    f11 = '('''//trim(prefix)//': '',a)'
    f15 = '('''//trim(prefix)//': '',a,t53,''= '',2(i0,'' x ''),i0)'
    f16 = '('''//trim(prefix)//': '',a,t53,''= A'',2(i0,'', A''),i0)'


    write(*,f11) '>> '//trim(name(this))
    write(*,f16) '  Electrode cell pivoting: E1, E2, E3', this%pvt
    if ( this%out_of_core ) then
       write(*,f10) '  GF file', trim(this%GFfile)
       write(*,f1)  '  Reuse existing GF-file', this%ReUseGF
    else
       write(*,f11)  '  In-core GF'
    end if
    write(*,f10) '  Electrode TSHS file', trim(this%HSfile)
    write(*,f5)  '  # atoms used in electrode', this%na_used
    write(*,f15) '  Electrode Bloch expansion [E1 x E2 x E3]', this%Bloch(:)
    write(*,f20) '  Position in geometry', this%idx_a, &
         this%idx_a + TotUsedAtoms(this) - 1
    if ( this%t_dir == 1 ) then
       chars = 'E1'
    else if ( this%t_dir == 2 ) then
       chars = 'E2'
    else if ( this%t_dir == 3 ) then
       chars = 'E3'
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
    if ( product(this%Bloch) > 1 .and. this%out_of_core ) then
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
#ifdef TBTRANS
#ifdef TBT_PHONON
    write(*,f9)  '  Electrode self-energy imaginary Eta', sqrt(this%Eta)/eV,' eV'
#else
    write(*,f9)  '  Electrode self-energy imaginary Eta', this%Eta/eV,' eV'
#endif
#else
    write(*,f9)  '  Electrode self-energy imaginary Eta', this%Eta/eV,' eV'
#endif
    write(*,f9)  '  Electrode self-energy accuracy', this%accu/eV,' eV'
    write(*,f6)  '  Electrode inter-layer distance (semi-inf)', this%dINF_layer/Ang,' Ang'


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
