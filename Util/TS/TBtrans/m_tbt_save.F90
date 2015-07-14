! Module for saving data

module m_tbt_save

  use precision, only : dp

  implicit none

  private

  integer, save :: cmp_lvl  = 0
  logical, save :: save_parallel = .false.
  
  ! Optional directory
  character(len=100), save :: save_dir = ' '
  public :: save_dir

  public :: init_save_options
  public :: name_save
#ifdef NCDF_4
  public :: tbt_cdf_precision
  public :: init_cdf_save
  public :: init_cdf_E_check
  public :: cdf_get_E_idx
  public :: cdf_get_kpt_idx
  public :: cdf_save_E
  public :: state_cdf_save
  public :: state_cdf_save_J
  public :: state_cdf_save_kpt
  public :: state_cdf2ascii
#else
  public :: init_save
  public :: step_kpt_save
  public :: state_save
  public :: end_save
#endif

  ! Type to control the energies that is contained
  ! This lets every processor know what the other processors have
  type :: tNodeE
     integer, allocatable :: iE(:)
     real(dp), allocatable :: E(:)
  end type tNodeE
  public :: tNodeE
  public :: MPI_BcastNode
  public :: save_parallel

contains

  subroutine MPI_BcastNode(iE, cE, nE)
#ifdef MPI
    use mpi_siesta, only : MPI_AllGather, MPI_Comm_World
    use mpi_siesta, only : MPI_Integer, MPI_Double_Precision
#endif
    integer, intent(in) :: iE
    complex(dp), intent(in) :: cE
    type(tNodeE), intent(inout) :: nE

#ifdef MPI
    real(dp) :: rE
    integer :: ierr

    rE = real(cE,dp)
    call MPI_AllGather(iE, 1, MPI_Integer, &
         nE%iE(0), 1, MPI_Integer, &
         MPI_Comm_World, ierr)
    call MPI_AllGather(rE, 1, MPI_Double_Precision, &
         nE%E(0), 1, MPI_Double_Precision, &
         MPI_Comm_World, ierr)
#else
    nE%iE(0) = iE
    nE%E(0) = real(cE,dp)
#endif
  end subroutine MPI_BcastNode
    

  subroutine init_save_options()
#ifdef NCDF_4
    use parallel, only : Node
#endif
    use fdf
    use m_io_s, only : dir_exist
    integer :: ldir

    cmp_lvl = fdf_get('CDF.Compress',0)
    cmp_lvl = fdf_get('TBT.CDF.Compress',cmp_lvl)
    if ( cmp_lvl < 0 ) cmp_lvl = 0
    if ( cmp_lvl > 9 ) cmp_lvl = 9
#ifdef NCDF_PARALLEL
    save_parallel = fdf_get('TBT.CDF.MPI',.false.)
    if ( save_parallel ) then
       cmp_lvl = 0
    end if
#endif

    save_dir = fdf_get('TBT.Directory.Save',' ')
    ! Correct with suffix
    ldir = len_trim(save_dir)
    if ( save_dir(ldir-1:ldir) == '/.' ) save_dir = save_dir(1:ldir-1)
    ldir = len_trim(save_dir)
    if ( save_dir(ldir:ldir) /= '/' ) save_dir = trim(save_dir)//'/'

    if ( save_parallel ) then
       if ( .not. dir_exist(save_dir, all = .true. ) ) then
          call die('Directory: '//trim(save_dir)//' not visible &
               &to all processors, or simply does not exist.')
       end if
    else if ( .not. dir_exist(save_dir, Bcast = .true. ) ) then
       call die('Directory: '//trim(save_dir)//' does not exist.')
    end if

#ifdef NCDF_4
    if ( Node == 0 ) then

       if ( len_trim(save_dir) > 0 ) then
          write(*,6)'Storing saved TBT files in',trim(save_dir)
       end if

       if ( cmp_lvl > 0 ) then
          write(*,5) 'Compression level of TBT.nc files',cmp_lvl
       else
          write(*,11)'No compression of TBT.nc files'
       end if

    end if

5   format('tbt_options: ',a,t53,'=',i5,a)
6   format('tbt_options: ',a,t53,'=',tr1,a)
11  format('tbt_options: ',a)

#endif

  end subroutine init_save_options

#ifdef NCDF_4
  subroutine tbt_cdf_precision(name,default,prec)

    use fdf, only : fdf_get, leqi
    use parallel, only : IONode
    use nf_ncdf, only : NF90_FLOAT, NF90_DOUBLE

    character(len=*), intent(in) :: name, default
    integer, intent(out) :: prec

    character(len=20) :: tmp

    ! Default, unless otherwise stated
    prec = NF90_FLOAT

    tmp = fdf_get('TBT.CDF.Precision',default)
    if ( leqi(tmp,'double') ) then
       prec = NF90_DOUBLE
    else if ( leqi(tmp,'single') ) then
       prec = NF90_FLOAT
    else if ( leqi(tmp,'float') ) then
       prec = NF90_FLOAT
    else if ( IONode ) then
       write(*,'(a)')'WARNING: Could not recognize TBT.CDF.Precision, &
            &must be double|single|float, defaults to double.'
    end if
    
    tmp = fdf_get('TBT.CDF.'//trim(name)//'.Precision',tmp)
    if ( leqi(tmp,'double') ) then
       prec = NF90_DOUBLE
    else if ( leqi(tmp,'single') ) then
       prec = NF90_FLOAT
    else if ( leqi(tmp,'float') ) then
       prec = NF90_FLOAT
    end if

  end subroutine tbt_cdf_precision

  subroutine init_cdf_save(fname,TSHS,r,ispin,N_Elec, Elecs, &
       nkpt, kpt, wkpt, NE, &
       a_Dev, a_Buf, sp_dev, &
       save_DATA ) 

    use parallel, only : Node

    use m_io_s, only : file_exist

    use dictionary
    use nf_ncdf, ncdf_parallel => parallel
    use m_ncdf_io, only : cdf_w_Sp
    use m_timestamp, only : datestring
#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD, MPI_Integer, MPI_Logical
    use mpi_siesta, only : MPI_Comm_Self
#endif
    use m_tbt_hs, only : tTSHS
    use m_ts_electype
    use m_region
    use class_OrbitalDistribution
    use class_Sparsity

    ! The file-name
    character(len=*), intent(in) :: fname
    ! The full Hamiltonian and system at present investigation.
    ! Note the H have been shifted to zero energy
    type(tTSHS), intent(in) :: TSHS
    ! The device region that we are checking
    ! This is the device regions pivot-table!
    type(tRgn), intent(in) :: r 
    integer, intent(in) :: ispin
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    integer, intent(in) :: nkpt
    real(dp), intent(in), target :: kpt(3,nkpt), wkpt(nkpt)
    integer, intent(in) :: NE
    type(tRgn), intent(in) :: a_Dev
    ! In case the system has some buffer atoms.
    type(tRgn), intent(in) :: a_Buf
    ! The device sparsity pattern
    type(Sparsity), intent(inout) :: sp_dev
    ! Options read from tbt_options
    type(dict), intent(in) :: save_DATA

    character(len=50) :: tmp
    type(hNCDF) :: ncdf, grp
    type(dict) :: dic
    logical :: exist, sme, isGamma
    integer :: iEl, jEl, i, nnzs_dev
    integer :: prec_DOS, prec_T, prec_J
    type(OrbitalDistribution) :: fdit
    real(dp), allocatable :: r2(:,:)
#ifdef MPI
    integer :: MPIerror
#endif

    ! In case the user thinks the double precision
    ! is too much
    call tbt_cdf_precision('DOS','single',prec_DOS)
    call tbt_cdf_precision('T','single',prec_T)
    call tbt_cdf_precision('Current','single',prec_J)

    isGamma = all(TSHS%nsc(:) == 1)

    ! If compiled with net-cdf we ALWAYS save in this format
    exist = file_exist(fname, Bcast = .true. )

    ! in case it already exists...
    if ( exist ) then

       ! We check the content, and if it is the same,
       ! we allow continuation.
       call ncdf_open(ncdf,fname,mode=NF90_NOWRITE)

       dic = ('no_u'.kv. TSHS%no_u) // ('na_u'.kv. TSHS%na_u ) &
            // ('no_d'.kv. r%n ) // ('nkpt'.kv. nkpt )
       dic = dic // ('na_d'.kv. a_Dev%n)
       if ( a_Buf%n > 0 ) then
          dic = dic // ('na_b'.kv. a_Buf%n)
       end if
       call ncdf_assert(ncdf,sme,dims=dic)
       call delete(dic)
#ifdef MPI
       call MPI_Bcast(sme,1,MPI_Logical,0, &
            MPI_Comm_World,MPIerror)
#endif
       if ( .not. sme ) then
          call die('Dimensions in the '//trim(fname)//' file &
               &does not conform to the current simulation.')
       end if

       ! Check the variables
       dic = ('lasto'.kvp. TSHS%lasto(1:TSHS%na_u) ) // &
            ('pivot'.kvp. r%r )
       dic = dic // ('xa'.kvp. TSHS%xa) // ('a_dev'.kvp.a_Dev%r )
       if ( a_Buf%n > 0 )then
          dic = dic // ('a_buf'.kvp.a_Buf%r )
       end if
       call ncdf_assert(ncdf,sme,vars=dic, d_EPS = 1.e-4_dp )
       call delete(dic,dealloc=.false.) ! we have them pointing...
#ifdef MPI
       call MPI_Bcast(sme,1,MPI_Logical,0, &
            MPI_Comm_World,MPIerror)
#endif
       if ( .not. sme ) then
          call die('pivot, lasto, xa or a_buf in the '//trim(fname)//' &
               &file does not conform to the current simulation.')
       end if

       if ( .not. isGamma ) then
          ! Check the k-points
          allocate(r2(3,nkpt))
          do i = 1 , nkpt
             call kpoint_convert(TSHS%cell,kpt(:,i),r2(:,i),1)
          end do
          dic = ('kpt'.kvp.r2) // ('wkpt'.kvp. wkpt)
          call ncdf_assert(ncdf,sme,vars=dic, d_EPS = 1.e-7_dp )
          if ( .not. sme ) then
             call die('k-points or k-weights are not the same')
          end if
          call delete(dic,dealloc = .false. )
          deallocate(r2)
       end if

       call ncdf_close(ncdf)
       
       ! The file has exactly the same content..
       if ( sme ) then

          ! We just need to set the weights to ensure
          ! unity weight.
          if ( Node == 0 .and. .not. isGamma ) then
             !call ncdf_open(ncdf,fname,mode=NF90_WRITE)
             !call ncdf_put_var(ncdf,'wkpt',wkpt,start=(/1/))
             !call ncdf_close(ncdf)
             write(*,'(a)') 'tbtrans: Continuation run on old TBT.nc file'
             write(*,'(a)') 'tbtrans: *** WARNING ***'
             write(*,'(a)') 'tbtrans: *** You need to make sure all &
                  &energy-points are fully contained in your current &
                  &energy range ***'
          end if

          call die('Currently the '//trim(fname)//' file exists, &
               &we do not currently implement a continuation scheme.')
          
          return

       end if

       ! We complain to the user about it and DIE
       call die('The file content in '//trim(fname)//' &
            &is not consistent with this setup. Please delete the &
            &file.')

    else
       
       if ( Node == 0 ) then
          write(*,'(2a)')'tbtrans: Initializing save file: ',trim(fname)
       end if

    end if

    ! We need to create the file
#ifdef NCDF_PARALLEL
    if ( save_parallel ) then
       call ncdf_create(ncdf,fname, mode=NF90_MPIIO, overwrite=.true., &
            comm = MPI_COMM_WORLD, &
            parallel = .true. )
    else
       call ncdf_create(ncdf,fname, mode=NF90_NETCDF4, overwrite=.true.)
    end if
#else
    call ncdf_create(ncdf,fname, mode=NF90_NETCDF4, overwrite=.true.)
#endif

    ! Save the current system size
    call ncdf_def_dim(ncdf,'no_u',TSHS%no_u)
    call ncdf_def_dim(ncdf,'na_u',TSHS%na_u)
    ! Even for Gamma, it makes files unified
    call ncdf_def_dim(ncdf,'nkpt',NF90_UNLIMITED)
    call ncdf_def_dim(ncdf,'xyz',3)
    call ncdf_def_dim(ncdf,'one',1)
    call ncdf_def_dim(ncdf,'na_d',a_Dev%n)
    call ncdf_def_dim(ncdf,'no_d',r%n)
    call ncdf_def_dim(ncdf,'ne',NF90_UNLIMITED)
    if ( a_Buf%n > 0 ) then
       call ncdf_def_dim(ncdf,'na_b',a_Buf%n) ! number of buffer-atoms
    end if

#ifdef TBT_PHONON
    dic = ('source'.kv.'PHtrans')
#else
    dic = ('source'.kv.'TBtrans')
#endif

    tmp = datestring()
    dic = dic//('date'.kv.tmp(1:10))
    if ( isGamma ) then
       dic = dic // ('Gamma'.kv.'true')
    else
       dic = dic // ('Gamma'.kv.'false')
    end if
    if ( TSHS%nspin > 1 ) then
       if ( ispin == 1 ) then
          dic = dic // ('spin'.kv.'UP')
       else
          dic = dic // ('spin'.kv.'DOWN')
       end if
    end if

    ! We initialize the counter for the current reached
    ! k-point and energy-point
    !if ( isGamma ) then
    !   dic = dic // ('k_idx_cur'.kv.1)
    !   dic = dic // ('k_idx_cur'.kv.0)
    !end if
    
    call ncdf_put_gatt(ncdf, atts = dic )
    call delete(dic)

    ! Create all the variables needed to save the states
    dic = ('info'.kv.'Last orbitals of the equivalent atom')
    call ncdf_def_var(ncdf,'lasto',NF90_INT,(/'na_u'/), &
         atts = dic)
    dic = dic//('info'.kv.'Unit cell')
    dic = dic//('unit'.kv.'Bohr')
    call ncdf_def_var(ncdf,'cell',NF90_DOUBLE,(/'xyz','xyz'/), &
         atts = dic)
    dic = dic//('info'.kv.'Atomic coordinates')
    dic = dic//('unit'.kv.'Bohr')
    call ncdf_def_var(ncdf,'xa',NF90_DOUBLE,(/'xyz ','na_u'/), &
         atts = dic , chunks = (/3, TSHS%na_u/) )
    call delete(dic)

    dic = ('info'.kv.'Device region orbital pivot table')
    call ncdf_def_var(ncdf,'pivot',NF90_INT,(/'no_d'/), &
         atts = dic)

    dic = dic//('info'.kv.'Index of device atoms')
    call ncdf_def_var(ncdf,'a_dev',NF90_INT,(/'na_d'/), &
         atts = dic)

    if ( a_Buf%n > 0 ) then
       dic = dic//('info'.kv.'Index of buffer atoms')
       call ncdf_def_var(ncdf,'a_buf',NF90_INT,(/'na_b'/), &
            atts = dic)
    end if

    if ( 'DOS-Gf' .in. save_DATA ) then
       dic = dic//('info'.kv.'Density of states')//('unit'.kv.'1/Ry')
       call ncdf_def_var(ncdf,'DOS',prec_DOS,(/'no_d','ne  ','nkpt'/), &
            atts = dic , chunks = (/r%n,1,1/) , compress_lvl = cmp_lvl )
    end if

    if ( .not. isGamma ) then

       dic = dic//('info'.kv.'k point')//('unit'.kv.'b')
       call ncdf_def_var(ncdf,'kpt',NF90_DOUBLE,(/'xyz ','nkpt'/), &
            atts = dic)
       call delete(dic)
       dic = ('info'.kv.'k point weights')
       call ncdf_def_var(ncdf,'wkpt',NF90_DOUBLE,(/'nkpt'/), &
            atts = dic , chunks = (/1/) )

    end if

#ifdef TBT_PHONON
    dic = dic//('info'.kv.'Frequency')//('unit'.kv.'Ry')
#else
    dic = dic//('info'.kv.'Energy')//('unit'.kv.'Ry')
#endif
    call ncdf_def_var(ncdf,'E',NF90_DOUBLE,(/'ne'/), &
         atts = dic, chunks = (/1/) )
    call delete(dic)

    call ncdf_put_var(ncdf,'pivot',r%r)
    call ncdf_put_var(ncdf,'cell',TSHS%cell)
    call ncdf_put_var(ncdf,'xa',TSHS%xa)
    call ncdf_put_var(ncdf,'lasto',TSHS%lasto(1:TSHS%na_u))
    call ncdf_put_var(ncdf,'a_dev',a_Dev%r)
    if ( a_Buf%n > 0 )then
       call ncdf_put_var(ncdf,'a_buf',a_Buf%r)
    end if

    ! Save all k-points
    ! Even though they are in an unlimited dimension,
    ! we save them instantly.
    ! This ensures that a continuation can check for 
    ! the same k-points in the same go.
    if ( .not. isGamma ) then
       allocate(r2(3,nkpt))
       do i = 1 , nkpt
          call kpoint_convert(TSHS%cell,kpt(:,i),r2(:,i),1)
       end do
       call ncdf_put_var(ncdf,'kpt',r2)
       call ncdf_put_var(ncdf,'wkpt',wkpt)
       deallocate(r2)
    end if

    if ( initialized(sp_dev) ) then
       ! In case we need to save the device sparsity pattern
       ! Create dimensions
       nnzs_dev = nnzs(sp_dev)
       call ncdf_def_dim(ncdf,'nnzs',nnzs_dev)

       call delete(dic)

       dic = ('info'.kv.'Number of non-zero elements per row')
       call ncdf_def_var(ncdf,'n_col',NF90_INT,(/'no_u'/), &
            compress_lvl=cmp_lvl,atts=dic)

       dic = dic//('info'.kv. &
            'Supercell column indices in the sparse format ')
       call ncdf_def_var(ncdf,'list_col',NF90_INT,(/'nnzs'/), &
            compress_lvl=cmp_lvl,atts=dic, chunks = (/nnzs_dev/) )

#ifdef MPI
       call newDistribution(TSHS%no_u,MPI_Comm_Self,fdit,name='TBT-fake dist')
#else
       call newDistribution(TSHS%no_u,-1           ,fdit,name='TBT-fake dist')
#endif

       call cdf_w_Sp(ncdf,fdit,sp_dev)
       call delete(fdit)
       call delete(dic)

    end if

    do iEl = 1 , N_Elec

       call delete(dic)

       if ( iEl == N_Elec ) then
          ! check if all are calculated
          if ( ('DOS-A-all' .nin. save_DATA) .and. &
               ('T-all'.nin. save_DATA) ) cycle
       end if

       call ncdf_def_grp(ncdf,trim(Elecs(iEl)%name),grp)

       if ( 'DOS-A' .in. save_DATA ) then
          dic = ('info'.kv.'Spectral function density of states')// &
               ('unit'.kv.'1/Ry')
          call ncdf_def_var(grp,'ADOS',prec_DOS,(/'no_d','ne  ','nkpt'/), &
               atts = dic, chunks = (/r%n,1,1/) , compress_lvl=cmp_lvl)
       end if

       ! Save information about electrode
       dic = dic//('info'.kv.'Chemical potential')//('unit'.kv.'Ry')
       call ncdf_def_var(grp,'mu',NF90_DOUBLE,(/'one'/), &
            atts = dic)
       call ncdf_put_var(grp,'mu',Elecs(iEl)%mu%mu)
#ifdef TBT_PHONON
       dic = dic//('info'.kv.'Phonon temperature')
#else
       dic = dic//('info'.kv.'Electronic temperature')
#endif
       call ncdf_def_var(grp,'kT',NF90_DOUBLE,(/'one'/), &
            atts = dic)
       call ncdf_put_var(grp,'kT',Elecs(iEl)%mu%kT)

       call delete(dic)

       if ( 'orb-current' .in. save_DATA ) then
          dic = ('info'.kv.'Orbital current')
          
          call ncdf_def_var(grp,'J',prec_J,(/'nnzs','ne  ','nkpt'/), &
               atts = dic , chunks = (/nnzs_dev/) , compress_lvl=cmp_lvl)
          
       end if
       
       tmp = trim(Elecs(iEl)%name)
       do jEl = 1 , N_Elec
          if ( ('T-all' .nin. save_DATA ) .and. &
               jEl < iEl ) cycle
          if ( ('T-reflect' .nin. save_DATA ) .and. &
               iEl == jEl ) cycle

          if ( iEl /= jEl ) then

             dic = dic//('info'.kv.'Transmission')
             call ncdf_def_var(grp,trim(Elecs(jEl)%name)//'.T',prec_T,(/'ne  ','nkpt'/), &
                  atts = dic )
             
          else

             ! For the same electrode we retain the group
             ! and utilise this for saving the reflection.
             dic = dic//('info'.kv.'Reflection')
             call ncdf_def_var(grp,trim(tmp)//'.R',prec_T,(/'ne  ','nkpt'/), &
                  atts = dic )

             dic = dic//('info'.kv.'Gf transmission')
             call ncdf_def_var(grp,trim(tmp)//'.T',prec_T,(/'ne  ','nkpt'/), &
                  atts = dic )

          end if
          
       end do

       call delete(dic)

    end do

    call ncdf_close(ncdf)

  end subroutine init_cdf_save

  subroutine init_cdf_E_check(fname,E,NE)

    use parallel, only : Node

    use nf_ncdf, ncdf_parallel => parallel
#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD, MPI_Bcast
    use mpi_siesta, only : MPI_Integer
    use mpi_siesta, only : Mpi_double_precision
#endif

    character(len=*), intent(in) :: fname
    real(dp), intent(inout), allocatable :: E(:)
    integer, intent(in) :: NE

    type(hNCDF) :: ncdf
    integer :: cur_NE
#ifdef MPI
    integer :: MPIerror
#endif

    ! Open the netcdf file
    if ( Node == 0 ) then

       call ncdf_open(ncdf,trim(fname), mode=NF90_NOWRITE)

       call ncdf_inq_dim(ncdf,'ne',len=cur_NE)

    end if

#ifdef MPI
    call MPI_BCast(cur_NE,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
#endif

    ! Allocate all previous E
    allocate(E(NE+cur_NE))
    E(:) = huge(1._dp) ! Initialize to something no sane person would sample

    call ncdf_get_var(ncdf,'E',E,count=(/cur_NE/))
#ifdef MPI
    call MPI_Bcast(E(1),cur_NE,MPI_DOUBLE_PRECISION,0, &
         MPI_Comm_World, MPIerror)
#endif

    call ncdf_close(ncdf)
    
  end subroutine init_cdf_E_check

  subroutine cdf_get_E_idx(E,NE,zE,iE,have)
    ! the last NE energy points are, "possibly" dublicate
    real(dp), intent(in) :: E(:)
    ! The current number of searching energy-points
    integer, intent(in) :: NE
    ! The current energy-point
    complex(dp), intent(in) :: zE
    ! The current index in the local pattern
    ! If this is 0 or negative, we know we
    ! are dealing with a fake energy-point.
    ! Then we will immediately return
    integer, intent(inout) :: iE
    logical, intent(out) :: have

    real(dp), parameter :: Eta = 7.349806700083788e-06_dp ! 0.0001 eV
    real(dp) :: rE
    integer :: i

    have = .false.

    ! It is already a fake energy point
    if ( iE <= 0 ) return

    rE = real(zE,dp)
    i = size(E)
    do iE = 1 , i
       if ( abs(rE - E(iE)) < Eta ) then
          ! We have a match
          ! if the index is not in the end NE points
          ! then we know that it is already in the
          ! file.
          have = iE <= i - NE
          return
       end if
    end do

  end subroutine cdf_get_E_idx

  subroutine cdf_get_kpt_idx(fname,bkpt,ikpt)

    use parallel, only : Node
    use nf_ncdf, ncdf_parallel => parallel
#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD, MPI_Bcast
    use mpi_siesta, only : MPI_Integer
#endif

    character(len=*), intent(in) :: fname
    real(dp), intent(in) :: bkpt(3)
    integer, intent(inout) :: ikpt

    type(hNCDF) :: ncdf
    real(dp), allocatable :: rkpt(:,:)
    ! This is limiting us at 1000000 k-points, in each direction
    ! Sigh...
    real(dp), parameter :: Eta = 1.e-6_dp 
    integer :: nk, ik
#ifdef MPI
    integer :: MPIerror
#endif

    ! First we check that the kpoint exists
    ! If it does not, then definetely we do not have the 
    ! energy point.

    if ( Node == 0 ) then
       
       call ncdf_open(ncdf,trim(fname), mode=NF90_WRITE)

       ! Retrieve the attributes 
       call ncdf_inq_dim(ncdf,'nkpt',len=nk)

       ! Initialize to signal that we are going to extend
       ! the position
       ikpt = nk + 1
       
       allocate(rkpt(3,nk))
       call ncdf_get_var(ncdf,'kpt',rkpt)
       
       ! We check the k-point
       ! When they are equal we have processed all energy points
       ! currently in it.
       ! This puts a restriction on the continuation
       ! process, we only allow increasing the energy density
       ! at will, increasing the k-points at will is
       ! only allowed by using the 
       ! user input k-point file and have the SAME k-points
       ! in the beginning (the USER HAS TO DO THIS!!!!!)

       do ik = 1 , nk
          if ( all(abs(rkpt(:,ik) - bkpt(:)) < Eta) ) then
             ikpt = ik
             exit
          end if
       end do
       
       if ( ikpt == nk + 1 ) then

          ! We need to add the k-point to the list
          call ncdf_put_var(ncdf,'kpt',bkpt,start=(/1,nk+1/))

       else

          ! we signal that the k-point 
          ! is already present in the save-file.
          ikpt = - ikpt

       end if

       call ncdf_close(ncdf)

    end if

#ifdef MPI
    call MPI_Bcast(ikpt,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
#endif

  end subroutine cdf_get_kpt_idx

  subroutine cdf_save_E(fname,nE)
    use parallel, only : Node, Nodes
    use nf_ncdf, ncdf_parallel => parallel

    character(len=*), intent(in) :: fname
    type(tNodeE), intent(in) :: nE

    type(hNCDF) :: ncdf
#ifdef MPI
    integer :: iN
#endif

    ! Open the netcdf file
    call ncdf_open(ncdf,fname, mode=NF90_WRITE)

    ! We save the energy
    if ( nE%iE(Node) > 0 ) then
       call ncdf_put_var(ncdf,'E',nE%E(Node),start = (/nE%iE(Node)/) )
    end if
    
#ifdef MPI
    do iN = 1 , Nodes - 1
       if ( nE%iE(iN) <= 0 ) cycle
       call ncdf_put_var(ncdf,'E',nE%E(iN),start = (/nE%iE(iN)/) )
    end do
#endif

    call ncdf_close(ncdf)

  end subroutine cdf_save_E

  subroutine state_cdf_save(fname, ikpt, nE, N_Elec, Elecs, DOS, T, &
       save_DATA)
    
    use parallel, only : Node, Nodes

    use dictionary
    use nf_ncdf, ncdf_parallel => parallel
#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD, MPI_Gather
    use mpi_siesta, only : MPI_Send, MPI_Recv, MPI_DOUBLE_COMPLEX
    use mpi_siesta, only : MPI_Integer, MPI_STATUS_SIZE
    use mpi_siesta, only : Mpi_double_precision
#endif
    use m_ts_electype

    character(len=*), intent(in) :: fname
    integer, intent(in) :: ikpt
    type(tNodeE), intent(in) :: nE
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    real(dp), intent(in) :: DOS(:,:)
    real(dp), intent(in) :: T(N_Elec+1,N_Elec)
    type(dict), intent(in) :: save_DATA

    type(hNCDF) :: ncdf, grp
    integer :: iEl, jEl
    character(len=30) :: tmp, tmp2
#ifdef MPI
    integer :: iN, NDOS, NT
    real(dp), allocatable :: rT(:,:,:)
    real(dp), allocatable :: rDOS(:)
    integer :: MPIerror, status(MPI_STATUS_SIZE)
#endif

#ifdef TBTRANS_TIMING
    call timer('cdf-w-T',1)
#endif

#ifdef MPI
    if ( .not. save_parallel .and. Nodes > 1 ) then
       NDOS = size(DOS,dim=1)
       allocate(rDOS(NDOS))
    end if
    NT = ( N_Elec + 1 ) * N_Elec
#endif

    ! Open the netcdf file
#ifdef NCDF_PARALLEL
    if ( save_parallel ) then
       call ncdf_open(ncdf,fname, mode=NF90_WRITE, &
            comm = MPI_COMM_WORLD )
    else
#endif
       call ncdf_open(ncdf,fname, mode=NF90_WRITE)
#ifdef NCDF_PARALLEL
    end if
#endif

    ! Save the different options given to this routine
    ! We need to save the energy
    if ( nE%iE(Node) > 0 ) then
       call ncdf_put_var(ncdf,'E',nE%E(Node),start = (/nE%iE(Node)/) )
    end if

#ifdef MPI
    if ( Node == 0 .and. .not. save_parallel ) then
       do iN = 1 , Nodes - 1
          if ( nE%iE(iN) <= 0 ) cycle
          call ncdf_put_var(ncdf,'E',nE%E(iN),start = (/nE%iE(iN)/) )
       end do
    end if
#endif

    if ( 'DOS-Gf' .in. save_DATA ) then

       ! We save the DOS
       ! This is the DOS from the Green's function
       call save_DOS(ncdf,'DOS',ikpt,nE,DOS(:,1))

    end if

    if ( 'DOS-A' .in. save_DATA ) then

       ! We save the DOS calculated from the spectral function

       do iEl = 1 , N_Elec
          if ( iEl == N_Elec ) then
             ! check if all are calculated
             if ( ('DOS-A-all' .nin. save_DATA) .and. &
                  ('T-all'.nin. save_DATA) ) cycle
          end if
          
          call ncdf_open_grp(ncdf,trim(Elecs(iEl)%name),grp)
          
          call save_DOS(grp,'ADOS',ikpt,nE,DOS(:,1+iEl))

       end do
       
    end if

#ifdef MPI
    if ( Node == 0 .and. .not. save_parallel ) then
       allocate(rT(N_Elec+1,N_Elec,Nodes-1))
       do iN = 1 , Nodes - 1
          call MPI_Recv(rT(1,1,iN),NT,Mpi_double_precision, &
               iN, iN, Mpi_comm_world,status,MPIerror)
       end do
    else if ( Node /= 0 .and. .not. save_parallel ) then
       call MPI_Send(T(1,1),NT,Mpi_double_precision, &
            0, Node, Mpi_comm_world,MPIerror)
    end if
#endif

    ! Save transmission function
    do iEl = 1 , N_Elec
       if ( iEl == N_Elec .and. ('T-all' .nin. save_DATA) ) cycle

       ! Open group of electrode
       call ncdf_open_grp(ncdf,trim(Elecs(iEl)%name),grp)

       do jEl = 1 , N_Elec
          if ( ('T-all' .nin. save_DATA ) .and. &
               jEl < iEl ) cycle
          if ( ('T-reflect' .nin. save_DATA ) .and. &
               iEl == jEl ) cycle

          if ( jEl == iEl ) then
             tmp  = trim(Elecs(jEl)%name)//'.R'
             tmp2 = trim(Elecs(jEl)%name)//'.T'
          else
             tmp  = trim(Elecs(jEl)%name)//'.T'
          end if

          ! Save data
          if ( nE%iE(Node) > 0 ) then
             call ncdf_put_var(grp,tmp,T(jEl,iEl),start = (/nE%iE(Node),ikpt/) )
             if ( iEl == jEl ) then
                call ncdf_put_var(grp,tmp2,T(N_Elec+1,iEl), &
                     start = (/nE%iE(Node),ikpt/) )
             end if
          end if
       
#ifdef MPI
          if ( Node == 0 .and. .not. save_parallel ) then
             do iN = 1 , Nodes - 1
                if ( nE%iE(iN) > 0 ) then
                   call ncdf_put_var(grp,tmp,rT(jEl,iEl,iN), &
                        start = (/nE%iE(iN),ikpt/) )
                   if ( iEl == jEl ) then
                      call ncdf_put_var(grp,tmp2,rT(N_Elec+1,iEl,iN), &
                           start = (/nE%iE(iN),ikpt/) )
                   end if
                end if
             end do
          end if
#endif

       end do
    end do

    call ncdf_close(ncdf)
       
#ifdef MPI
    if ( allocated(rDOS) ) deallocate(rDOS)
    if ( allocated(rT) ) deallocate(rT)
#endif

#ifdef TBTRANS_TIMING
    call timer('cdf-w-T',2)
#endif

  contains

    subroutine save_DOS(grp,var,ikpt,nE,DOS)
      type(hNCDF), intent(inout) :: grp
      character(len=*), intent(in) :: var
      integer, intent(in) :: ikpt
      type(tNodeE), intent(in) :: nE
      real(dp), intent(in) :: DOS(:)

      if ( nE%iE(Node) > 0 ) then
         call ncdf_put_var(grp,var,DOS,start = (/1,nE%iE(Node),ikpt/) )
      end if
      
#ifdef MPI
      if ( .not. save_parallel ) then
         if ( Node == 0 ) then
            do iN = 1 , Nodes - 1
               if ( nE%iE(iN) <= 0 ) cycle
               call MPI_Recv(rDOS,NDOS,Mpi_double_precision,iN,iN, &
                    Mpi_comm_world,status,MPIerror)
               call ncdf_put_var(grp,var,rDOS,start = (/1,nE%iE(iN),ikpt/) )
            end do
         else if ( nE%iE(Node) > 0 ) then
            call MPI_Send(DOS,NDOS,Mpi_double_precision,0,Node, &
                 Mpi_comm_world,MPIerror)
         end if
      end if
#endif

    end subroutine save_DOS

  end subroutine state_cdf_save

  subroutine state_cdf_save_J(fname, ikpt, nE, El, orb_J, save_DATA)
    
    use parallel, only : Node, Nodes
    use class_dSpData1D

    use dictionary
    use nf_ncdf, ncdf_parallel => parallel
#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD
    use mpi_siesta, only : MPI_Send, MPI_Recv
    use mpi_siesta, only : MPI_STATUS_SIZE
    use mpi_siesta, only : Mpi_double_precision
#endif
    use m_ts_electype

    character(len=*), intent(in) :: fname
    integer, intent(in) :: ikpt
    type(tNodeE), intent(in) :: nE
    type(Elec), intent(in) :: El
    type(dSpData1D), intent(inout) :: orb_J
    type(dict), intent(in) :: save_DATA

    type(hNCDF) :: ncdf, grp
    integer :: nnzs_dev
    real(dp), pointer :: J(:)
#ifdef MPI
    integer :: iN
    integer :: MPIerror, status(MPI_STATUS_SIZE)
#endif

    ! Open the netcdf file
#ifdef NCDF_PARALLEL
    if ( save_parallel ) then
       call ncdf_open(ncdf,fname, mode=NF90_WRITE, &
            comm = MPI_COMM_WORLD )
    else
#endif
       call ncdf_open(ncdf,fname, mode=NF90_WRITE)
#ifdef NCDF_PARALLEL
    end if
#endif

    J => val(orb_J)
    nnzs_dev = size(J)
    ! We save the orbital current
    
    call ncdf_open_grp(ncdf,trim(El%name),grp)

    ! Save the current
    if ( nE%iE(Node) > 0 ) then
       call ncdf_put_var(grp,'J',J,start = (/1,nE%iE(Node),ikpt/) )
    end if

#ifdef MPI
    if ( .not. save_parallel ) then
       if ( Node == 0 ) then
          do iN = 1 , Nodes - 1
             if ( nE%iE(iN) > 0 ) then
                call MPI_Recv(J,nnzs_dev,Mpi_double_precision, &
                     iN, iN, Mpi_comm_world,status,MPIerror)
                call ncdf_put_var(grp,'J',J,start = (/1,nE%iE(iN),ikpt/) )
             end if
          end do
       else if ( nE%iE(Node) > 0 ) then
          call MPI_Send(J(1),nnzs_dev,Mpi_double_precision, &
               0, Node, Mpi_comm_world,MPIerror)
       end if
    end if
#endif
       
    call ncdf_close(ncdf)
       
  end subroutine state_cdf_save_J

  subroutine state_cdf_save_kpt(fname,ikpt)

    use parallel, only : Node

    use nf_ncdf, ncdf_parallel => parallel

    ! We step the k-point index to indicate that we 
    ! have calculated all k-points.
    character(len=*), intent(in) :: fname
    integer, intent(in) :: ikpt

    type(hNCDF) :: ncdf

    if ( Node /= 0 ) return

    call ncdf_open(ncdf,fname, mode=NF90_WRITE)

    call ncdf_put_gatt(ncdf,'k_idx_cur',ikpt)

    call ncdf_close(ncdf)

  end subroutine state_cdf_save_kpt

  ! Routine for reading in the TBT.nc file
  ! and convert it to regular transmission files.
  subroutine state_cdf2ascii(fname,nspin,ispin,N_Elec,Elecs,N_E,rW,save_DATA)

    use parallel, only : Node
    use units, only : eV
#ifdef TBT_PHONON
    use units, only : Kelvin
#endif

    use variable
    use dictionary
    use nf_ncdf, ncdf_parallel => parallel

    use m_interpolate, only : crt_pivot

    use m_timestamp, only : datestring
    use m_ts_electype

    character(len=*), intent(in) :: fname
    integer, intent(in) :: nspin, ispin
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    integer, intent(in) :: N_E
    ! This quantity is the dE weight, with dE in Ry units
    real(dp), intent(in) :: rW(N_E)
    type(dict), intent(in) :: save_DATA

    character(len=250) :: ascii_file, tmp
    type(hNCDF) :: ncdf, grp
    logical :: exist
    integer :: iEl, jEl, i
    integer :: NE, nkpt, no_d
    real(dp), allocatable :: rkpt(:,:), rwkpt(:)
    real(dp), allocatable :: rE(:)
    real(dp), allocatable :: r2(:,:), r3(:,:,:)
#ifdef TBT_PHONON
    real(dp) :: Flow, dT
#else
    real(dp) :: Current, Power, V, dd
#endif
    integer, allocatable :: pvt(:)

    call timer('cdf2ascii',1)

    ! In case we are doing something parallel, 
    ! we simply read in and write them in text based formats
    if ( Node /= 0 ) then
       call timer('cdf2ascii',2)
       return
    end if

    tmp = datestring()
    tmp = tmp(1:10)

    ! Open the netcdf file
    call ncdf_open(ncdf,fname, mode=NF90_NOWRITE)

    ! First we read in all dimensions
    call ncdf_inq_dim(ncdf,'ne',len=NE)
    if ( NE /= N_E ) call die('Error when re-reading the number of &
         &energy-points')
    call ncdf_inq_dim(ncdf,'nkpt',len=nkpt)
    call ncdf_inq_dim(ncdf,'no_d',len=no_d)

    ! Allocate space
    allocate(rE(NE),pvt(NE))
    allocate(rkpt(3,nkpt),rwkpt(nkpt))

    ! Read in common information
    call ncdf_get_var(ncdf,'E',rE)
    ! Convert energy to eV
    rE = rE / eV
    ! Create pivot table
    call crt_pivot(NE,rE,pvt)
    call ncdf_inq_var(ncdf,'kpt',exist=exist)
    if ( exist ) then
       call ncdf_get_var(ncdf,'kpt',rkpt)
       call ncdf_get_var(ncdf,'wkpt',rwkpt)
    else
       ! It MUST be a Gamma calculation
       rkpt = 0._dp
       rwkpt = 1._dp
    end if

    if ( 'DOS-Gf' .in. save_DATA ) then

       allocate(r3(no_d,NE,nkpt))

       ! Get DOS
       call ncdf_get_var(ncdf,'DOS',r3)

       ! Correct unit, from 1/Ry to 1/eV
!$OMP parallel workshare default(shared)
       r3 = r3 * eV
!$OMP end parallel workshare
       
       if ( nkpt > 1 ) then
          call name_save(ispin,nspin,ascii_file,end='DOS')
          call save_DAT(ascii_file,nkpt,rkpt,rwkpt,NE,rE,pvt,no_d,r3,'DOS', &
               '# DOS calculated from the Green function, k-resolved')
       end if
       call name_save(ispin,nspin,ascii_file,end='AVDOS')
       call save_DAT(ascii_file,1,rkpt,rwkpt,NE,rE,pvt,no_d,r3,'DOS', &
            '# DOS calculated from the Green function, k-averaged')

       deallocate(r3)

    end if

#ifdef TBT_PHONON
    if ( Node == 0 .and. N_Elec > 1 ) then
       write(*,'(/,a)')'Heatflow (ensure freqency range covers temperature tails):'
    end if
#else
    if ( Node == 0 .and. N_Elec > 1 ) then
       write(*,'(/,a)')'Currents (ensure entire Fermi function window):'
    end if
#endif

    ! We should now be able to create all the files
    do iEl = 1 , N_Elec

       ! We do not calculate the last electrode
       ! unless requested
       if ( iEl == N_Elec ) then
          ! check if all are calculated
          if ( ('DOS-A-all' .nin. save_DATA) .and. &
               ('T-all'.nin. save_DATA) ) cycle
       end if


       call ncdf_open_grp(ncdf,trim(Elecs(iEl)%name),grp)

       if ( 'DOS-A' .in. save_DATA ) then

          allocate(r3(no_d,NE,nkpt))
          call ncdf_get_var(grp,'ADOS',r3)

          ! Correct unit, from 1/Ry to 1/eV
!$OMP parallel workshare default(shared)
          r3 = r3 * eV
!$OMP end parallel workshare
          
          if ( nkpt > 1 ) then
             call name_save(ispin,nspin,ascii_file,end='ADOS',El1=Elecs(iEl))
             call save_DAT(ascii_file,nkpt,rkpt,rwkpt,NE,rE,pvt,no_d,r3,'DOS',&
                  '# DOS calculated from the spectral function, k-resolved')
          end if
          call name_save(ispin,nspin,ascii_file,end='AVADOS',El1=Elecs(iEl))
          call save_DAT(ascii_file,1,rkpt,rwkpt,NE,rE,pvt,no_d,r3,'DOS', &
               '# DOS calculated from the spectral function, k-averaged')

          deallocate(r3)
          
       end if

       allocate(r2(NE,nkpt))

       do jEl = 1 , N_Elec
          ! Calculating iEl -> jEl is the
          ! same as calculating jEl -> iEl, hence if we
          ! do not wish to assert this is true, we do not
          ! calculate this.
          if ( ('T-all' .nin. save_DATA ) .and. &
               jEl < iEl ) cycle
          if ( ('T-reflect' .nin. save_DATA ) .and. &
               iEl == jEl ) cycle

          if ( iEl == jEl ) then
             call ncdf_get_var(grp,trim(Elecs(jEl)%name)//'.R',r2)
             ! Save the variable to ensure the correct sum in the transmission
             if ( nkpt > 1 ) then
                call name_save(ispin,nspin,ascii_file,end='REFL', El1=Elecs(iEl))
                call save_DAT(ascii_file,nkpt,rkpt,rwkpt,NE,rE,pvt,1,r2,'Reflection',&
                     '# Reflection, k-resolved')
             end if
             call name_save(ispin,nspin,ascii_file,end='AVREFL', El1=Elecs(iEl))
             call save_DAT(ascii_file,1,rkpt,rwkpt,NE,rE,pvt,1,r2,'Reflection',&
                  '# Reflection, k-averaged')

             ! The transmission is now the total incoming wave 
             ! [G-G^\dagger].\Gamma
             call ncdf_get_var(grp,trim(Elecs(jEl)%name)//'.T',r2)
          else
             call ncdf_get_var(grp,trim(Elecs(jEl)%name)//'.T',r2)
          end if
          
          ! Save transmission
          if ( nkpt > 1 ) then
             call name_save(ispin,nspin,ascii_file,end='TRANS', &
                  El1=Elecs(iEl), El2=Elecs(jEl))
             call save_DAT(ascii_file,nkpt,rkpt,rwkpt,NE,rE,pvt,1,r2,'Transmission',&
                  '# Transmission, k-resolved')
          end if
          call name_save(ispin,nspin,ascii_file,end='AVTRANS', &
               El1=Elecs(iEl), El2=Elecs(jEl))
          call save_DAT(ascii_file,1,rkpt,rwkpt,NE,rE,pvt,1,r2,'Transmission',&
               '# Transmission, k-averaged')

          ! The array r2 now contains the k-averaged transmission.
#ifdef TBT_PHONON
          ! Now we calculate the heat-flow
          ! nb function is: nb(E-E1) - nb(E-E2) IMPORTANT
          Flow = 0._dp
!$OMP parallel do default(shared), private(i), &
!$OMP&reduction(+:Flow)
          do i = 1 , NE
             ! We have rE in eV, hence the conversion
             Flow = Flow + r2(i,1) * rW(i) * rE(i) * nb(rE(i)*eV, &
                  Elecs(iEl)%mu%mu, Elecs(iEl)%mu%kT, &
                  Elecs(jEl)%mu%mu, Elecs(jEl)%mu%kT )
          end do
!$OMP end parallel do

          ! rE is already in eV, r2 and nb are unit-less
          ! rW is in Ry => / eV
          !     eV to 1 / s => / hbar[eV s]
          !     1 / s to 1 / fs => / 1e15
          Flow = Flow / eV / 0.658211928_dp
          dT = ( Elecs(iEl)%mu%kT - Elecs(jEl)%mu%kT ) / Kelvin

          if ( Node == 0 ) then
             write(*,'(4a,2(g12.6,a))') trim(Elecs(iEl)%name), &
                  ' -> ',trim(Elecs(jEl)%name),', dT [K] / E-flow [eV/fs]: ', &
                  dT, ' K / ',Flow,' eV/fs'
          end if
#else
          ! Now we calculate the current
          ! nf function is: nF(E-E1) - nF(E-E2) IMPORTANT
          Current = 0._dp
          Power = 0._dp
          if ( iEl == jEl ) then
             ! Do nothing
          else
!$OMP parallel do default(shared), private(i,dd), &
!$OMP&reduction(+:Current,Power)
             do i = 1 , NE
                ! We have rE in eV, hence the conversion
                dd = r2(i,1) * rW(i) * nf2(rE(i)*eV, &
                     Elecs(iEl)%mu%mu, Elecs(iEl)%mu%kT, &
                     Elecs(jEl)%mu%mu, Elecs(jEl)%mu%kT )
                Current = Current + dd
                ! rE is in eV, mu is in Ry
                Power = Power + dd * ( rE(i) - Elecs(iEl)%mu%mu / eV )
             end do
!$OMP end parallel do
          end if

          ! This quantity is e / h * [J] / [eV]: [q] / [J] / [s] * [J] / [eV] = [q] / [eV] / [s]
          ! NOTE: [q] = [J] / [eV] in numeric quantity
          dd = 3.87405e-05_dp

          ! 1. rW is in [Ry], hence we need to convert to [eV] Ry => / eV
          ! 2. Current is in [eV], so dividing by [q] / [eV] / [s] yields [A]
          Current = Current / eV * dd
          ! 1. rW is in [Ry], hence we need to convert to [eV] Ry => / eV
          ! 2. Power is in [eV] ** 2, so dividing by [J] / [eV] / [eV] / [s] yields [W]
          Power = Power / eV * dd

          V = ( Elecs(iEl)%mu%mu - Elecs(jEl)%mu%mu ) / eV

          if ( Node == 0 .and. iEl /= jEl ) then
             write(*,'(4a,2(g12.6,a))') trim(Elecs(iEl)%name), &
                  ' -> ',trim(Elecs(jEl)%name),', V [V] / I [A]: ', &
                  V, ' V / ',Current,' A'
             write(*,'(4a,2(g12.6,a))') trim(Elecs(iEl)%name), &
                  ' -> ',trim(Elecs(jEl)%name),', V [V] / P [W]: ', &
                  V, ' V / ',Power,' W'
          end if
#endif

       end do
       
       deallocate(r2)

    end do

    if ( Node == 0 ) write(*,*) ! new-line

    ! Clean-up
    deallocate(rE,rkpt,rwkpt,pvt)

    call ncdf_close(ncdf)

    call timer('cdf2ascii',2)

  contains
    
    subroutine save_DAT(fname,nkpt,kpt,wkpt,NE,E,ipiv,no_d,DAT,value,header)
      character(len=*), intent(in) :: fname
      integer, intent(in) :: nkpt, NE, no_d, ipiv(NE)
      real(dp), intent(in) :: kpt(3,nkpt), wkpt(nkpt), E(NE)
      real(dp), intent(inout) :: DAT(no_d,NE,nkpt)
      character(len=*), intent(in) :: value, header
      integer :: iu, ik, i
      real(dp) :: rno
      if ( no_d == 1 ) then
         rno = 1._dp
      else
         rno = 1._dp / no_d
      end if
      
      call io_assign(iu)
      open( iu, file=trim(fname), form='formatted', status='unknown' ) 
      
      write(iu,'(a)') trim(header)
      write(iu,'(a)') '# Date: '//trim(tmp)
#ifdef TBT_PHONON
      write(iu,'(a,a9,tr1,a16)')"#","Omega [eV]", value
#else
      write(iu,'(a,a9,tr1,a16)')"#","E [eV]", value
#endif
      do ik = 1 , nkpt 
         if ( nkpt > 1 ) then
            write(iu,'(/,a6,3(f10.6,'', ''),a,f10.6)') &
                 '# kb  = ',kpt(:,ik) ,'w= ',wkpt(ik)
         end if
         do i = 1 , NE
            ! We sum the orbital contributions
            write(iu,'(f10.5,tr1,e16.8)') E(ipiv(i)),sum(DAT(:,ipiv(i),ik)) * rno
         end do
         if ( nkpt > 1 ) then
            if ( ik == 1 ) then
!$OMP parallel workshare default(shared)
               DAT(:,:,1) = DAT(:,:,1) * wkpt(ik)
!$OMP end parallel workshare
            else
!$OMP parallel workshare default(shared)
               DAT(:,:,1) = DAT(:,:,1) + DAT(:,:,ik) * wkpt(ik)
!$OMP end parallel workshare
            end if
         end if
      end do

      call io_close(iu)
      
    end subroutine save_DAT

#ifdef TBT_PHONON
    elemental function nb(E,E1,kT1,E2,kT2)
      real(dp), intent(in) :: E,E1,kT1,E2,kT2
      real(dp) :: nb
      nb = 1._dp/(exp((E-E1)/kT1)-1._dp) - 1._dp/(exp((E-E2)/kT2)-1._dp)
    end function nb
#else
    elemental function nf2(E,E1,kT1,E2,kT2)
      real(dp), intent(in) :: E,E1,kT1,E2,kT2
      real(dp) :: nf2
      nf2 = nf(E,E1,kT1) - nf(E,E2,kT2)
    end function nf2
    elemental function nf(E,Ef,kT)
      real(dp), intent(in) :: E,Ef,kT
      real(dp) :: nf
      nf = 1._dp/(exp((E-Ef)/kT)+1._dp)
    end function nf
#endif

  end subroutine state_cdf2ascii

#endif

  ! Get the file name
  subroutine name_save(ispin,nspin,fname,end,El1,El2)
    use files, only : slabel
    use m_ts_electype
    integer, intent(in) :: ispin, nspin
    character(len=*), intent(out) :: fname
    character(len=*), intent(in), optional :: end ! designator of the file
    type(Elec), intent(in), optional :: El1, El2

    fname = trim(save_dir)//trim(slabel)//'.TBT'

    ! Now figure out the file name
    if ( nspin > 1 ) then
       if( ispin .eq. 1 ) fname = trim(fname)//"_UP"
       if( ispin .eq. 2 ) fname = trim(fname)//"_DN"
    end if

    ! Add the designator
    if ( present(end) ) then
       fname = trim(fname)//'.'//trim(end)
    end if

    if ( present(El1) ) then
       fname = trim(fname)//'_'//trim(El1%name)
       if ( present(El2) ) then
          fname = trim(fname)//'-'//trim(El2%name)
       end if
    end if

  end subroutine name_save

#ifndef NCDF_4
  ! These routines are for creating the output data in
  ! pure ASCII format.
  

  ! This routine prepares the files for saving ASCII format
  ! data.
  ! NOTE that ASCII data will only be created in case
  ! of Netcdf not being compiled in
  subroutine init_save(iounits,ispin,nspin,N_Elec,Elecs,save_DATA)
    
    use parallel, only : Node

    use dictionary
    use m_timestamp, only : datestring
#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD
#endif
    use m_ts_electype

    integer, intent(inout) :: iounits(:)
    integer, intent(in)    :: ispin, nspin
    integer, intent(in)    :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    type(dict), intent(in) :: save_DATA

    character(len=100) :: ascii_file, tmp
    integer :: iu, cu

    integer :: iEl, jEl

    ! Only the IO-Node can prepare the output data.
    if ( Node /= 0 ) return

    tmp = datestring()
    tmp = tmp(1:10)

    cu = 1

    if ( 'DOS-Gf' .in. save_DATA ) then

       call name_save(ispin,nspin,ascii_file,end='DOS')

       call io_assign(iu)
       open( iu, file=trim(ascii_file), form='formatted', status='unknown' ) 
       write(iu,'(a)') '# DOS calculated from the Greens function, k-resolved'
       write(iu,'(a)') '# Date: '//trim(tmp)
       write(iu,'(a,a9,tr1,a16)')'#','E [eV]', 'DOS'
       
       iounits(cu) = iu

       cu = cu + 1
       
    end if

    do iEl = 1 , N_Elec

       ! We do not calculate the last electrode
       ! unless requested
       if ( iEl == N_Elec ) then
          ! check if all are calculated
          if ( ('DOS-A-all' .nin. save_DATA) .and. &
               ('T-all'.nin. save_DATA) ) cycle
       end if
       
       if ( 'DOS-A' .in. save_DATA ) then

          call name_save(ispin,nspin,ascii_file,end='ADOS',El1=Elecs(iEl))

          call io_assign(iu)
          open( iu, file=trim(ascii_file), form='formatted', status='unknown' ) 
          write(iu,'(a)') '# DOS calculated from the spectral function, k-resolved'
          write(iu,'(a)') '# Date: '//trim(tmp)
          write(iu,'(a,a9,tr1,a16)')'#','E [eV]', 'DOS'

          iounits(cu) = iu
          
          cu = cu + 1

       end if

       do jEl = 1 , N_Elec

          ! Calculating iEl -> jEl is the
          ! same as calculating jEl -> iEl, hence if we
          ! do not wish to assert this is true, we do not
          ! calculate this.
          if ( ('T-all' .nin. save_DATA ) .and. &
               jEl < iEl ) cycle
          if ( ('T-reflect' .nin. save_DATA ) .and. &
               iEl == jEl ) cycle


          if ( iEl == jEl ) then

             call name_save(ispin,nspin,ascii_file,end='REFL', &
                  El1=Elecs(iEl))
             
             call io_assign(iu)
             open( iu, file=trim(ascii_file), form='formatted', status='unknown' ) 
             write(iu,'(a)') '# Reflection, k-resolved'
             write(iu,'(a)') '# Date: '//trim(tmp)
             write(iu,'(a,a9,tr1,a16)')'#','E [eV]', 'Reflection'

             iounits(cu) = iu
             
             cu = cu + 1

          end if
             
          call name_save(ispin,nspin,ascii_file,end='TRANS', &
               El1=Elecs(iEl), El2=Elecs(jEl))

          call io_assign(iu)
          open( iu, file=trim(ascii_file), form='formatted', status='unknown' ) 
          write(iu,'(a)') '# Transmission, k-resolved'
          write(iu,'(a)') '# Date: '//trim(tmp)
          write(iu,'(a,a9,tr1,a16)')'#','E [eV]', 'Transmission'

          iounits(cu) = iu
          
          cu = cu + 1

       end do
       
    end do

  end subroutine init_save

  subroutine step_kpt_save(iounits,nkpt,bkpt,wkpt,N_Elec,save_DATA)
    
    use parallel, only : Node

    use dictionary

    integer,    intent(in) :: iounits(:), nkpt
    real(dp),   intent(in) :: bkpt(3), wkpt
    integer,    intent(in) :: N_Elec
    type(dict), intent(in) :: save_DATA
    integer :: cu, iEl, jEl

    if ( Node /= 0 ) return
    if ( nkpt == 1 ) return

    cu = 1

    if ( 'DOS-Gf' .in. save_DATA ) then

       call wrt_k(iounits(cu))
       cu = cu + 1
       
    end if

    do iEl = 1 , N_Elec

       ! We do not calculate the last electrode
       ! unless requested
       if ( iEl == N_Elec ) then
          ! check if all are calculated
          if ( ('DOS-A-all' .nin. save_DATA) .and. &
               ('T-all'.nin. save_DATA) ) cycle
       end if
       
       if ( 'DOS-A' .in. save_DATA ) then

          call wrt_k(iounits(cu))
          cu = cu + 1
          
       end if
       
       do jEl = 1 , N_Elec
          
          ! Calculating iEl -> jEl is the
          ! same as calculating jEl -> iEl, hence if we
          ! do not wish to assert this is true, we do not
          ! calculate this.
          if ( ('T-all' .nin. save_DATA ) .and. &
               jEl < iEl ) cycle
          if ( ('T-reflect' .nin. save_DATA ) .and. &
               iEl == jEl ) cycle

          if ( iEl == jEl ) then

             call wrt_k(iounits(cu))
             cu = cu + 1

          end if
          
          call wrt_k(iounits(cu))
          cu = cu + 1
          
       end do
       
    end do

  contains

    subroutine wrt_k(iu)
      integer, intent(in) :: iu
      write(iu,'(/,a6,3(f10.6,'', ''),a,f10.6)') &
           '# kb  = ',bkpt(:) ,'w= ',wkpt
    end subroutine wrt_k

  end subroutine step_kpt_save

  subroutine end_save(iounits,N_Elec,save_DATA)
    
    use parallel, only : Node

    use dictionary

    integer,    intent(in) :: iounits(:), N_Elec
    type(dict), intent(in) :: save_DATA
    integer :: cu, iEl, jEl

    if ( Node /= 0 ) return

    cu = 1

    if ( 'DOS-Gf' .in. save_DATA ) then

       call io_close(iounits(cu))
       cu = cu + 1
       
    end if

    do iEl = 1 , N_Elec

       ! We do not calculate the last electrode
       ! unless requested
       if ( iEl == N_Elec ) then
          ! check if all are calculated
          if ( ('DOS-A-all' .nin. save_DATA) .and. &
               ('T-all'.nin. save_DATA) ) cycle
       end if
       
       if ( 'DOS-A' .in. save_DATA ) then
          
          call io_close(iounits(cu))
          cu = cu + 1
          
       end if
       
       do jEl = 1 , N_Elec
          
          ! Calculating iEl -> jEl is the
          ! same as calculating jEl -> iEl, hence if we
          ! do not wish to assert this is true, we do not
          ! calculate this.
          if ( ('T-all' .nin. save_DATA ) .and. &
               jEl < iEl ) cycle
          if ( ('T-reflect' .nin. save_DATA ) .and. &
               iEl == jEl ) cycle

          if ( iEl == jEl ) then
             call io_close(iounits(cu))
             cu = cu + 1
          end if
          
          call io_close(iounits(cu))
          cu = cu + 1
          
       end do
       
    end do

  end subroutine end_save

  ! This routine prepares the files for saving ASCII format
  ! data.
  ! NOTE that ASCII data will only be created in case
  ! of Netcdf not being compiled in
  subroutine state_save(iounits,nE,N_Elec,Elecs,DOS, T, &
       save_DATA )
    
    use parallel, only : Node, Nodes
    use units, only : eV
    use m_interpolate, only : crt_pivot

    use dictionary
#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD, MPI_Gather
    use mpi_siesta, only : MPI_Send, MPI_Recv, MPI_DOUBLE_COMPLEX
    use mpi_siesta, only : MPI_Integer, MPI_STATUS_SIZE
    use mpi_siesta, only : Mpi_double_precision
#endif
    use m_ts_electype

    integer, intent(in)    :: iounits(:)
    type(tNodeE), intent(in) :: nE
    integer, intent(in)    :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    real(dp), intent(in)   :: DOS(:,:), T(:,:)
    type(dict), intent(in) :: save_DATA

    integer :: cu
    integer :: iEl, jEl
    real(dp) :: rE
    integer, allocatable :: ipvt(:)
#ifdef MPI
    integer :: NDOS
    real(dp), allocatable :: rDOS(:)
    integer :: MPIerror, status(MPI_STATUS_SIZE)
#endif

    allocate(ipvt(0:Nodes-1))
    ipvt = 1

#ifdef MPI
    if ( Nodes > 1 ) then
       call crt_pivot(Nodes,nE%E,ipvt)
       NDOS = size(DOS,dim=1)
       allocate(rDOS(NDOS))
    end if
#endif

    ! Correct for rank indices
    ipvt(:) = ipvt(:) - 1

    cu = 1

    if ( 'DOS-Gf' .in. save_DATA ) then

       call save_DAT(iounits(cu),DOS(:,1),fact=eV)

       cu = cu + 1

    end if

    do iEl = 1 , N_Elec

       ! We do not calculate the last electrode
       ! unless requested
       if ( iEl == N_Elec ) then
          ! check if all are calculated
          if ( ('DOS-A-all' .nin. save_DATA) .and. &
               ('T-all'.nin. save_DATA) ) cycle
       end if
       
       if ( 'DOS-A' .in. save_DATA ) then

          call save_DAT(iounits(cu),DOS(:,1+iEl),fact=eV)
          
          cu = cu + 1
          
       end if

       do jEl = 1 , N_Elec

          ! Calculating iEl -> jEl is the
          ! same as calculating jEl -> iEl, hence if we
          ! do not wish to assert this is true, we do not
          ! calculate this.
          if ( ('T-all' .nin. save_DATA ) .and. &
               jEl < iEl ) cycle
          if ( ('T-reflect' .nin. save_DATA ) .and. &
               iEl == jEl ) cycle

          call save_DAT(iounits(cu),T(jEl:jEl,iEl))
          
          cu = cu + 1

          if ( jEl == iEl ) then
             ! Note this is reversed according to the 
             ! creation of the arrays.
             ! This is because the reflection is 1->1
             ! and the transmission is the G.\Gamma
             ! flux. Hence we simply reverse the print-outs.
             call save_DAT(iounits(cu),T(N_Elec+1:N_Elec+1,iEl))
             cu = cu + 1
          end if

       end do
       
    end do

    deallocate(ipvt)
#ifdef MPI
    if ( allocated(rDOS) ) deallocate(rDOS)
#endif

  contains
    
    subroutine save_DAT(iu,DATA,fact)
      integer, intent(in) :: iu
      real(dp), intent(in) :: DATA(:)
      real(dp), intent(in), optional :: fact

      integer :: iN, i, ND
      real(dp) :: rnd
      ND = size(DATA)
      if ( ND == 1 ) then
         rnd = 1._dp
      else
         rnd = 1._dp / ND
      end if
      if ( present(fact) ) rnd = rnd * fact

      if ( Node == 0 ) then
         do iN = 0 , Nodes - 1
            i = ipvt(iN) ! sorting E
#ifdef MPI
            if ( nE%iE(i) <= 0 ) cycle ! if the energy point is fake, discard
#endif
            if ( i == 0 ) then ! local node
               write(iu,'(f10.5,tr1,e16.8)') nE%E(i) / eV,sum(DATA(:)) * rnd
            else
#ifdef MPI
               call MPI_Recv(rDOS,ND,Mpi_double_precision,i,i, &
                    Mpi_comm_world,status,MPIerror)
               write(iu,'(f10.5,tr1,e16.8)') nE%E(i) / eV,sum(rDOS(1:ND)) * rnd
#else
               call die('Error')
#endif
            end if
         end do
      else if ( nE%iE(Node) > 0 ) then
#ifdef MPI
         call MPI_Send(DATA(1),ND,Mpi_double_precision,0,Node, &
              Mpi_comm_world,MPIerror)
#endif
      end if

    end subroutine save_DAT
    
  end subroutine state_save

#endif

end module m_tbt_save
