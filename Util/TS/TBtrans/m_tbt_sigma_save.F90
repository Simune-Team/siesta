
! Creation of the hssigma files.

module m_tbt_sigma_save

  use units, only : dp

  use m_tbt_hs, only : tTSHS
  use m_tbt_save, only : tNodeE
  
  implicit none

  private 

  logical, save :: sigma_save      = .false.
  logical, save :: sigma_mean_save = .false.
  logical, save :: sigma_parallel  = .false.
  integer, save :: compress_lvl    = 0

#ifdef NCDF_4
  public :: init_Sigma_options
  public :: init_Sigma_save
  public :: state_Sigma_save
  public :: state_Sigma2mean
#endif

contains

#ifdef NCDF_4

  subroutine init_Sigma_options(save_DATA)

    use dictionary
    use parallel, only : Node
    use fdf

    type(dict), intent(inout) :: save_DATA

    sigma_save   = fdf_get('TBT.Sigma.CDF.Save',.false.)
    if ( sigma_save ) then
       sigma_mean_save = fdf_get('TBT.Sigma.CDF.Save.Mean',.false.)
    end if
    compress_lvl = fdf_get('TBT.CDF.Compress',0)
    compress_lvl = fdf_get('TBT.Sigma.CDF.Compress',compress_lvl)
    if ( compress_lvl < 0 ) compress_lvl = 0
    if ( compress_lvl > 9 ) compress_lvl = 9
#ifdef NCDF_PARALLEL
    sigma_parallel = fdf_get('TBT.Sigma.CDF.MPI',.false.)
    if ( sigma_parallel ) then
       compress_lvl = 0
    end if
#endif

    if ( fdf_get('TBT.Sigma.Only',.false.) ) then
       save_DATA = save_DATA // ('Sigma-only'.kv.1)
    end if

    if ( Node == 0 ) then

       write(*,1)'Saving down-folded self-energies',sigma_save
       write(*,1)'Only calc down-folded self-energies', &
            ('Sigma-only'.in.save_DATA)
       if ( sigma_save ) then
          if ( compress_lvl > 0 ) then
             write(*,5)'Compression level of TBT.Sigma.nc files',compress_lvl
          else
             write(*,11)'No compression level of TBT.Sigma.nc files'
          end if
          write(*,1)'k-average down-folded self-energies',sigma_mean_save
       end if

    end if

1   format('tbt_options: ',a,t53,'=',4x,l1)
5   format('tbt_options: ',a,t53,'=',i5,a)
11  format('tbt_options: ',a)
    
  end subroutine init_Sigma_options

  ! Save the self-energies of the electrodes and
  subroutine init_Sigma_save(fname, TSHS, r, ispin, N_Elec, Elecs, &
       nkpt, kpt, wkpt, NE, &
       r_Buf)

    use parallel, only : Node
    use m_io_s, only : file_exist

    use nf_ncdf
    use m_timestamp, only : datestring
#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD, MPI_Bcast, MPI_Logical
#endif
    use m_ts_electype
    use m_region
    use dictionary

    ! The file name that we save in
    character(len=*), intent(in) :: fname
    ! The full Hamiltonian and system at present investigation.
    ! Note the H have been shifted to zero energy
    type(tTSHS), intent(in) :: TSHS
    ! The device region that we are checking
    ! This is the device regions pivot-table!
    type(tRegion), intent(in) :: r 
    integer, intent(in) :: ispin
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    integer, intent(in) :: nkpt
    real(dp), intent(in) :: kpt(3,nkpt), wkpt(nkpt)
    integer, intent(in) :: NE
    ! Buffer atoms
    type(tRegion), intent(in) :: r_Buf

    type(hNCDF) :: ncdf, grp
    type(dict) :: dic
    logical :: exist, isGamma, same
    character(len=200) :: char
    integer :: i, iEl
    real(dp), allocatable :: r2(:,:)
#ifdef MPI
    integer :: MPIerror
#endif

    if ( .not. sigma_save ) return

    isGamma = all(TSHS%nsc(:) == 1)

    exist = file_exist(fname, Bcast = .true. )

    ! in case it already exists...
    if ( exist ) then

       ! Create a dictionary to check that the sigma file is the
       ! same
       call ncdf_open(ncdf,fname)

       dic = ('no_u'.kv.TSHS%no_u) // ('na_u'.kv.TSHS%na_u) // &
            ('nkpt'.kv.nkpt ) // ('no_d'.kv.r%n) // &
            ('ne'.kv. NE )
       if ( r_Buf%n > 0 ) then
          dic = dic // ('na_b'.kv.r_Buf%n)
       end if
       call ncdf_assert(ncdf,exist,dims=dic)
       call delete(dic)
#ifdef MPI
       call MPI_Bcast(same,1,MPI_Logical,0, &
            MPI_Comm_World,MPIerror)
#endif
       if ( .not. same ) then
          call die('Dimensions in the TBT.nc file does not conform &
               &to the current simulation.')
       end if

       do iEl = 1 , N_Elec
          call ncdf_open_grp(ncdf,Elecs(iEl)%name,grp)
          dic = dic // ('no_e'.kv.Elecs(iEl)%o_inD%n)
          call ncdf_assert(grp,exist,dims=dic)
          if ( .not. exist ) then
             write(*,*) 'Assertion of dimensions in file: '//trim(fname)//' failed.'

             call die('We could not assert the dimensions TBT.Sigma.nc file.')
          end if
       end do
       call delete(dic)

       ! Check the variables
       ! Check the variables
       dic = ('lasto'.kvp. TSHS%lasto(1:TSHS%na_u) ) // &
            ('pivot'.kvp. r%r )
       dic = dic // ('xa'.kvp. TSHS%xa)
       if ( r_Buf%n > 0 )then
          dic = dic // ('a_buf'.kvp.r_Buf%r )
       end if
       call ncdf_assert(ncdf,same,vars=dic, d_EPS = 1.e-4_dp )
       call delete(dic,dealloc=.false.) ! we have them pointing...
#ifdef MPI
       call MPI_Bcast(same,1,MPI_Logical,0, &
            MPI_Comm_World,MPIerror)
#endif
       if ( .not. same ) then
          call die('pivot, lasto, xa or a_buf in the TBT.nc file does &
               &not conform to the current simulation.')
       end if

       if ( .not. isGamma ) then
          ! Check the k-points
          allocate(r2(3,nkpt))
          do i = 1 , nkpt
             call kpoint_convert(TSHS%cell,kpt(:,i),r2(:,i),1)
          end do
          dic = ('kpt'.kvp. r2) // ('wkpt'.kvp. wkpt)
          call ncdf_assert(ncdf,same,vars=dic, d_EPS = 1.e-7_dp )
          if ( .not. same ) then
             call die('k-points or k-weights are not the same')
          end if
          call delete(dic,dealloc = .false. )
          deallocate(r2)
       end if

       call die('Currently the TBT.Sigma.nc file exists, &
            &we do not currently implement a continuation scheme.')

       ! We currently overwrite the Sigma-file
       if ( Node == 0 ) then
          write(*,'(2a)')'tbtrans: Overwriting Sigma file: ',trim(fname)
       end if

    else
       
       if ( Node == 0 ) then
          write(*,'(2a)')'tbtrans: Initializing Sigma file: ',trim(fname)
       end if

    end if

    ! We need to create the file
#ifdef NCDF_PARALLEL
    if ( sigma_parallel ) then
       call ncdf_create(ncdf,fname, mode=NF90_MPIIO, overwrite=.true., &
            comm = MPI_COMM_WORLD, &
            parallel = .true. )
    else
       call ncdf_create(ncdf,fname, mode=NF90_NETCDF4, overwrite=.true., &
            compress_lvl = compress_lvl )
    end if
#else
    call ncdf_create(ncdf,fname, mode=NF90_NETCDF4, overwrite=.true., &
         compress_lvl = compress_lvl )
#endif

    ! Save the current system size
    call ncdf_def_dim(ncdf,'no_u',TSHS%no_u)
    call ncdf_def_dim(ncdf,'na_u',TSHS%na_u)
    call ncdf_def_dim(ncdf,'nkpt',nkpt) ! Even for Gamma, it makes files unified
    call ncdf_def_dim(ncdf,'xyz',3)
    call ncdf_def_dim(ncdf,'one',1)
    call ncdf_def_dim(ncdf,'ne',NE)
    call ncdf_def_dim(ncdf,'no_d',r%n)
    if ( r_Buf%n > 0 ) then
       call ncdf_def_dim(ncdf,'na_b',r_Buf%n) ! number of buffer-atoms
    end if

    dic = ('source'.kv.'TBtrans')

    char = datestring()
    dic = dic//('date'.kv.char(1:10))
    if ( all(TSHS%nsc(:) == 1) ) then
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
         atts = dic)
    call delete(dic)

    dic = ('info'.kv.'Device region orbital pivot table')
    call ncdf_def_var(ncdf,'pivot',NF90_INT,(/'no_d'/), &
         atts = dic)

    if ( r_Buf%n > 0 ) then
       dic = dic//('info'.kv.'Index of buffer atoms')
       call ncdf_def_var(ncdf,'a_buf',NF90_INT,(/'na_b'/), &
            atts = dic)
    end if

    if ( .not. isGamma ) then

       dic = dic//('info'.kv.'k point')//('unit'.kv.'b')
       call ncdf_def_var(ncdf,'kpt',NF90_DOUBLE,(/'xyz ','nkpt'/), &
            atts = dic)
       call delete(dic)
       dic = dic//('info'.kv.'k point weights')
       call ncdf_def_var(ncdf,'wkpt',NF90_DOUBLE,(/'nkpt'/), &
            atts = dic)

    end if

    dic = dic//('info'.kv.'Energy points')//('unit'.kv.'Ry')
    call ncdf_def_var(ncdf,'E',NF90_DOUBLE,(/'ne'/), &
         atts = dic)
    call delete(dic)

    call ncdf_put_var(ncdf,'pivot',r%r)
    call ncdf_put_var(ncdf,'cell',TSHS%cell)
    call ncdf_put_var(ncdf,'xa',TSHS%xa)
    call ncdf_put_var(ncdf,'lasto',TSHS%lasto(1:TSHS%na_u))
    if ( r_Buf%n > 0 ) then
       call ncdf_put_var(ncdf,'a_buf',r_Buf%r)
    end if

    ! Save all k-points
    if ( .not. isGamma ) then
       allocate(r2(3,nkpt))
       do i = 1 , nkpt
          call kpoint_convert(TSHS%cell,kpt(:,i),r2(:,i),1)
       end do
       call ncdf_put_var(ncdf,'kpt',r2)
       call ncdf_put_var(ncdf,'wkpt',wkpt)
       deallocate(r2)
    end if

    do iEl = 1 , N_Elec

       call ncdf_def_grp(ncdf,trim(Elecs(iEl)%name),grp)

       ! Save information about electrode
       dic = dic//('info'.kv.'Chemical potential')//('unit'.kv.'Ry')
       call ncdf_def_var(grp,'mu',NF90_DOUBLE,(/'one'/), &
            atts = dic)
       call ncdf_put_var(grp,'mu',Elecs(iEl)%mu%mu)

       dic = ('info'.kv.'Imaginary part of self-energy')//('unit'.kv.'Ry')
       call ncdf_def_var(grp,'Eta',NF90_DOUBLE,(/'one'/), atts = dic)
       call delete(dic)

       call ncdf_def_dim(grp,'no_e',Elecs(iEl)%o_inD%n)

       dic = ('info'.kv.'Orbital pivot table for self-energy')
       call ncdf_def_var(grp,'pivot',NF90_INT,(/'no_e'/), atts = dic)
       dic = dic//('info'.kv.'Down-folded self-energy')
       dic = dic//('unit'.kv.'Ry')
       ! Chunking greatly reduces IO cost
       i = Elecs(iEl)%o_inD%n
       call ncdf_def_var(grp,'Sigma',NF90_DOUBLE_COMPLEX, &
            (/'no_e','no_e','ne  ','nkpt'/), &
            atts = dic , chunks = (/i,i,1,1/) )
       call delete(dic)

       call ncdf_put_var(grp,'Eta',Elecs(iEl)%Eta)
       call ncdf_put_var(grp,'pivot',Elecs(iEl)%o_inD%r)

    end do

    call ncdf_close(ncdf)

  end subroutine init_Sigma_save

  subroutine state_Sigma_save(fname, ikpt, nE, N_Elec, Elecs)

    use parallel, only : Node, Nodes

    use nf_ncdf
#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD, MPI_Gather
    use mpi_siesta, only : MPI_Send, MPI_Recv, MPI_DOUBLE_COMPLEX
    use mpi_siesta, only : MPI_Integer, MPI_STATUS_SIZE
    use mpi_siesta, only : Mpi_double_precision
#endif
    use m_ts_electype

    ! The file name we save too
    character(len=*), intent(in) :: fname
    integer, intent(in) :: ikpt
    type(tNodeE), intent(in) :: nE
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)

    type(hNCDF) :: ncdf, grp
    integer :: iEl, i
#ifdef MPI
    complex(dp), allocatable :: Sigma(:)
    integer :: MPIerror, iN, status(MPI_STATUS_SIZE)
#endif

    if ( .not. sigma_save ) return

#ifdef NCDF_PARALLEL
    if ( sigma_parallel ) then
       call ncdf_open(ncdf,fname, mode=NF90_WRITE, &
            comm = MPI_COMM_WORLD )
    else
#endif
       call ncdf_open(ncdf,fname, mode=NF90_WRITE)
#ifdef NCDF_PARALLEL
    end if
#endif

    ! Save the energy-point
    if ( parallel_io(ncdf) ) then
       if ( nE%iE(Node) > 0 ) then
          call ncdf_put_var(ncdf,'E',nE%E(Node),start = (/nE%iE(Node)/) )
       end if
    else
       do iN = 0 , Nodes - 1
          if ( nE%iE(iN) <= 0 ) cycle
          call ncdf_put_var(ncdf,'E',nE%E(iN),start = (/nE%iE(iN)/) )
       end do
    end if

#ifdef MPI
    if ( .not. sigma_parallel .and. Nodes > 1 ) then
       i = 0
       do iEl = 1 , N_Elec
          i = max(i,Elecs(iEl)%o_inD%n)
       end do
       allocate(Sigma(i**2))
#endif

    end if

    do iEl = 1 , N_Elec
       
       call ncdf_open_grp(ncdf,trim(Elecs(iEl)%name),grp)

       i = Elecs(iEl)%o_inD%n
       if ( nE%iE(Node) > 0 ) then
          call ncdf_put_var(grp,'Sigma', &
               reshape(Elecs(iEl)%Sigma(1:i*i),(/i,i/)), &
               start = (/1,1,nE%iE(Node),ikpt/) )
       end if

#ifdef MPI
       if ( .not. sigma_parallel .and. Nodes > 1 ) then
          if ( Node == 0 ) then
             do iN = 1 , Nodes - 1
                if ( nE%iE(iN) <= 0 ) cycle
                call MPI_Recv(Sigma,i*i,Mpi_double_complex,iN,iN, &
                     Mpi_comm_world,status,MPIerror)
                call ncdf_put_var(grp,'Sigma',reshape(Sigma(1:i*i),(/i,i/)), &
                     start = (/1,1,nE%iE(iN),ikpt/) )
             end do
          else if ( nE%iE(Node) > 0 ) then
             call MPI_Send(Elecs(iEl)%Sigma(1),i*i,Mpi_double_complex,0,Node, &
                  Mpi_comm_world,MPIerror)
          end if
       end if
#endif

    end do

#ifdef MPI
    if ( allocated(Sigma) ) deallocate(Sigma)
#endif

    call ncdf_close(ncdf)
    
  end subroutine state_Sigma_save

  subroutine state_Sigma2mean(fname,N_Elec,Elecs)

    use parallel, only : Node

    use dictionary
    use nf_ncdf
#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD, MPI_Gather
    use mpi_siesta, only : MPI_Send, MPI_Recv, MPI_DOUBLE_COMPLEX
    use mpi_siesta, only : MPI_Integer, MPI_STATUS_SIZE
    use mpi_siesta, only : Mpi_double_precision
    use mpi_siesta, only : MPI_Barrier
#endif
    use m_ts_electype

    ! The file name we save too
    character(len=*), intent(in) :: fname
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)

    type(dict) :: dic
    type(hNCDF) :: ncdf, grp
    integer :: iEl, iE, ikpt
    integer :: NE, nkpt, no_e
    real(dp), allocatable :: rwkpt(:)
    complex(dp), allocatable :: c2(:,:)
    complex(dp), pointer :: Sigma(:,:)

#ifdef MPI
    integer :: MPIerror
#endif

    ! If we should not save the mean, we return immediately.
    if ( .not. sigma_mean_save ) return

    if ( Node /= 0 ) then
#ifdef MPI
       call MPI_Barrier(Mpi_comm_world,MPIerror)
#endif
       return
    end if

    ! We do this on one processor
    call ncdf_open(ncdf,fname, mode=NF90_WRITE)

    ! We read in the dimensions
    call ncdf_inq_dim(ncdf,'ne',len=NE)
    call ncdf_inq_dim(ncdf,'nkpt',len=nkpt)

    ! Allocate space
    allocate(rwkpt(nkpt))
    call ncdf_get_var(ncdf,'wkpt',rwkpt)

    ! When taking the mean of self-energies
    ! we need the transpose, hence we need half the
    ! contribution from Sigma and Sigma^T
    rwkpt(:) = 0.5_dp * rwkpt(:)

    ! Loop over all electrodes
    do iEl = 1 , N_Elec

       ! We need to extend the netcdf file with the SigmaMean
       ! variable

       call ncdf_open_grp(ncdf,trim(Elecs(iEl)%name),grp)

       ! Get size of Sigma
       call ncdf_inq_dim(grp,'no_e',len=no_e)

       dic = ('info'.kv.'Down-folded self-energy, k-averaged')
       dic = dic//('unit'.kv.'Ry')
       ! Chunking greatly reduces IO cost
       call ncdf_def_var(grp,'SigmaMean',NF90_DOUBLE_COMPLEX, &
            (/'no_e','no_e','ne  '/), chunks = (/no_e,no_e,1/) , &
            atts = dic )
       call delete(dic)

       ! Allocate space for the self-energy mean
       allocate(c2(no_e,no_e))

       ! Point the sigma
       ! This is a hack to ease the processing
       call pass2pnt(no_e,Elecs(iEl)%Sigma(1:no_e**2),Sigma)

       ! loop over all energy points
       do iE = 1 , NE

          do ikpt = 1 , nkpt

             ! Loop over k-points to average
             call ncdf_get_var(grp,'Sigma',Sigma, &
                  start=(/1,1,iE,ikpt/) )
             
             if ( ikpt == 1 ) then
                c2(:,:) = rwkpt(ikpt) * ( &
                     Sigma + transpose(Sigma) &
                     )
             else
                c2(:,:) = c2(:,:) + rwkpt(ikpt) * ( &
                     Sigma + transpose(Sigma) &
                     )
             end if
             
          end do

          call ncdf_put_var(grp,'SigmaMean',c2, start=(/1,1,iE/) )

       end do

       deallocate(c2)

    end do

    deallocate(rwkpt)
    
    call ncdf_close(ncdf)

#ifdef MPI
    call MPI_Barrier(Mpi_comm_world,MPIerror)
#endif

  end subroutine state_Sigma2mean

  subroutine pass2pnt(no,Sigma,new_pnt)
    integer :: no
    complex(dp), target :: Sigma(no,no)
    complex(dp), pointer :: new_pnt(:,:)
    new_pnt => Sigma(:,:)
  end subroutine pass2pnt

#endif

end module m_tbt_sigma_save

  

