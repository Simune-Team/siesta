module m_ts_electype

  use precision, only : dp

  use class_Sparsity
  use class_dSpData1D
  use class_dSpData2D

  implicit none

  private

  public :: assign, read_Elec
  public :: create_sp2sp01
  public :: delete_TSHS
  
  public :: Elec
  public :: HSfile, GFFile, GFTitle
  public :: Atoms, UsedAtoms, TotUsedAtoms
  public :: Orbs, UsedOrbs, TotUsedOrbs
  public :: SCOrbs
  public :: unitcell
  public :: spin, EFermi
  public :: Rep
  public :: RepA1, RepA2, RepA3

  public :: INF_NEGATIVE, INF_POSITIVE

  integer, parameter :: FILE_LEN = 200

  integer, parameter :: INF_NEGATIVE = 0
  integer, parameter :: INF_POSITIVE = 1

  integer, parameter :: HAS_NOTHING = 0
  integer, parameter :: HAS_HS = 1
  integer, parameter :: HAS_HS00_HS01 = 2

  type :: Elec
     character(len=FILE_LEN) :: HSfile, GFfile, GFtitle
     ! These variables are relative to the big system
     integer :: idx_na
     integer :: na_used
     integer :: no_used
     integer :: RepA1 = 1, RepA2 = 1, RepA3 = 1
     real(dp) :: mu
     integer :: inf_dir = INF_NEGATIVE
     integer :: t_dir = 3
     ! Used xa and lasto
     real(dp), pointer :: xa_used(:,:) => null()
     integer, pointer :: lasto_used(:) => null()
     ! ---v--- Below we have the content of the TSHS file
     integer :: nspin, na_u, no_u, no_s
     real(dp) :: ucell(3,3), Ef, Qtot
     real(dp), pointer :: xa(:,:) => null()
     integer, pointer :: lasto(:) => null()
     type(Sparsity) :: sp
     type(dSpData2D) :: H, xij
     type(dSpData1D) :: S
     ! --- --- completed the content of the TSHS file
     ! Below we create the content for the self-energy creation
     type(Sparsity)  :: sp00,               sp01
     type(dSpData2D) :: H00, xij00, xijo00, H01, xij01, xijo01
     type(dSpData1D) :: S00,                S01
  end type Elec

contains

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
    character(len=FILE_LEN) :: GFtitle
    GFtitle = this%GFtitle
  end function GFtitle

  elemental function Spin(this) result(val)
    type(Elec), intent(in) :: this
    integer :: val
    val = this%nspin
  end function Spin

  elemental function EFermi(this) result(val)
    type(Elec), intent(in) :: this
    real(dp) :: val
    val = this%Ef
  end function EFermi

  elemental function Atoms(this) result(val)
    type(Elec), intent(in) :: this
    integer :: val
    val = this%na_u
  end function Atoms

  elemental function UsedAtoms(this) result(val)
    type(Elec), intent(in) :: this
    integer :: val
    val = this%na_used
  end function UsedAtoms
  elemental function TotUsedAtoms(this) result(val)
    type(Elec), intent(in) :: this
    integer :: val
    val = this%na_used * Rep(this)
  end function TotUsedAtoms

  elemental function Rep(this) result(val)
    type(Elec), intent(in) :: this
    integer :: val
    val = RepA1(this)*RepA2(this)*RepA3(this)
  end function Rep
  elemental function RepA1(this) result(val)
    type(Elec), intent(in) :: this
    integer :: val
    val = this%RepA1
  end function RepA1
  elemental function RepA2(this) result(val)
    type(Elec), intent(in) :: this
    integer :: val
    val = this%RepA2
  end function RepA2
  elemental function RepA3(this) result(val)
    type(Elec), intent(in) :: this
    integer :: val
    val = this%RepA3
  end function RepA3

  pure function q_exp(this,i,j,k) result(q)
    type(Elec), intent(in) :: this
    integer, intent(in) :: i,j,k
    real(dp) :: q(3)
    q(1) = 1._dp*(i-1) / real(RepA1(this),dp)
    q(2) = 1._dp*(j-1) / real(RepA2(this),dp)
    q(3) = 1._dp*(k-1) / real(RepA3(this),dp)
  end function q_exp

  elemental function Orbs(this) result(val)
    type(Elec), intent(in) :: this
    integer :: val
    val = this%no_u
  end function Orbs
  elemental function UsedOrbs(this) result(val)
    type(Elec), intent(in) :: this
    integer :: val
    val = this%no_used
  end function UsedOrbs
  elemental function TotUsedOrbs(this) result(val)
    type(Elec), intent(in) :: this
    integer :: val
    val = this%no_used * Rep(this)
  end function TotUsedOrbs
  elemental function SCOrbs(this) result(val)
    type(Elec), intent(in) :: this
    integer :: val
    val = this%no_s
  end function SCOrbs

  function unitcell(this) result(val)
    type(Elec), intent(in) :: this
    real(dp) :: val(3,3)
    val = this%ucell
  end function unitcell


  subroutine read_Elec(this,Bcast)
    use fdf
    use m_ts_io
    use parallel
    use class_OrbitalDistribution
#ifdef MPI
    use mpi_siesta
#endif

    type(Elec), intent(inout) :: this
    logical, intent(in), optional :: Bcast

    character(len=200) :: fN
    integer :: fL, kscell(3,3), istep, ia1
    logical :: onlyS, Gamma_file, TSGamma
    real(dp) :: temp, kdispl(3)
    
    ! Sparsity pattern
    integer, pointer :: iza(:), numh(:), listhptr(:), listh(:), indxuo(:)
    real(dp), pointer :: H(:,:), S(:), xij(:,:), t2D(:,:), t1D(:)
    integer :: n_nzs

    type(OrbitalDistribution) :: fdist

    fN = trim(HSfile(this))
    ! We read in the information
    fL = len_trim(fN)
    if ( leqi(fN(fL-4:fL),'.TSHS') ) then
       call ts_read_tshs(fN, &
            onlyS, Gamma_file, TSGamma, &
            this%ucell, this%na_u, this%no_u, this%no_u, this%no_s, n_nzs, this%nspin,  &
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
    t2D = H
    t1D => val(this%S)
    t1D = S

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

    call print_type(this%sp)

  end subroutine read_Elec


  subroutine create_sp2sp01(this,calc_xijo)
    use parallel, only : Node

    use class_OrbitalDistribution

    use create_Sparsity_SC
    use geom_helper, only : iaorb,cell_c
#ifdef MPI
    use mpi_siesta
#endif

    type(Elec), intent(inout) :: this
    logical, intent(in), optional :: calc_xijo

    logical :: lcalc_xijo

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

    integer :: no_l, i, io, j, ind, ind00, ind01, ia, ja,t
    integer :: tm(3)

    ! Retrieve distribution
    fdist => dist(this%H)

    lcalc_xijo = .false.
    if ( present(calc_xijo) ) lcalc_xijo = calc_xijo

    H   => val(this%H)
    S   => val(this%S)
    xij => val(this%xij)
    tm(:) = TM_ALL
    tm(this%t_dir) = 0
    call crtSparsity_SC(this%sp,this%sp00, &
         TM=tm, ucell=this%ucell, &
         lasto=this%lasto, xa=this%xa, xij=xij)
    tm(this%t_dir) = 1
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
         nrows=io)
    if ( io /= no_l ) call die('Could not do index matching due to &
         &inconsistent sparsity patterns')
    call attach(this%sp01,n_col=ncol01,list_ptr=ptr01,list_col=col01, &
         nrows=io)
    if ( io /= no_l ) call die('Could not do index matching due to &
         &inconsistent sparsity patterns')

    ! loop and assign data elements
    do i = 1 , no_l

       ! Shift out of the buffer region
       io = index_local_to_global(fdist,i,Node)
       ia = iaorb(io,this%lasto)

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

    call print_type(this%sp00)
    call print_type(this%sp01)

  end subroutine create_sp2sp01

  subroutine delete_TSHS(this)
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

  end subroutine delete_TSHS

end module m_ts_electype
