! Module for easy writing of all information
! in a single netCDF file.

! We rely on the NetCDF4 library and can
! utilise parallel IO or standard IO

module m_ncdf_io

  use precision, only : sp, dp, grid_p
  use parallel, only : Node, Nodes
#ifdef NCDF_4
  use variable
  use dictionary
  use nf_ncdf, ncdf_parallel => parallel
#endif
#ifdef MPI
  use mpi_siesta
#endif

  implicit none

  private

#ifdef NCDF_4

  public :: cdf_init_file
  public :: cdf_save_settings
  public :: cdf_save_md
  public :: cdf_save_state
  public :: cdf_w_sp
  public :: cdf_w_d1D, cdf_w_d2D

#endif

contains

#ifdef NCDF_4

  subroutine cdf_init_file(fname,is_md)
    use class_Sparsity
    use m_io_s, only : file_exist
    use dictionary
    use files, only : slabel
    use m_gamma, only : Gamma
    use atomlist, only: no_u, no_s, lasto, Qtot
    use siesta_geom, only: na_u
    use sparse_matrices, only: sparse_pattern
    use m_spin, only : nspin
    use siesta_options, only: cdf_comp_lvl, cdf_w_parallel
    use siesta_options, only: sname, isolve
    use siesta_options, only: SOLVE_DIAGON, SOLVE_ORDERN, SOLVE_TRANSI
    use m_timestamp, only: datestring
#ifdef TRANSIESTA
    use m_ts_electype, elec_name => name
    use m_ts_options, only : Volt, N_Elec, Elecs
#endif

    character(len=*), intent(in) :: fname
    logical, intent(in) :: is_md

    ! Local variables
    type(hNCDF) :: ncdf, grp, grp2
    type(dict) :: dic, d
    type(var) :: avar
    character(len=DICT_KEY_LENGTH) :: key
    integer :: n_nzs, tmp, iEl, i
#ifdef TRANSIESTA
    integer, allocatable :: ibuf(:)
#endif
    logical :: exists
#ifdef MPI
    integer :: MPIerror
#endif

    ! We always re-write the file...

    if ( Nodes > 1 .and. cdf_w_parallel ) then
       call ncdf_create(ncdf,fname,&
            mode=NF90_MPIIO, parallel = .true., &
            overwrite = .true. , &
            compress_lvl=cdf_comp_lvl, &
            comm=MPI_Comm_World)
       ! parallel writes are not allowed with compression
       ! Offset positions are not well defined.
    else
       call ncdf_create(ncdf,fname,&
            mode=NF90_NETCDF4, overwrite=.true., &
            compress_lvl=cdf_comp_lvl)
       
    end if

#ifdef MPI
    tmp = nnzs(sparse_pattern)
    call MPI_Reduce(tmp,n_nzs,1,MPI_Integer, &
         MPI_Sum, 0, MPI_Comm_World, MPIerror)
#else
    n_nzs = nnzs(sparse_pattern)
#endif
    
    ! First we create all the dimensions
    ! necessary
    d = ('DIMna_u'.kv.na_u)//('DIMno_u'.kv.no_u)
    d = d//('DIMno_s'.kv.no_s)//('DIMspin'.kv.nspin)
    d = d//('DIMxyz'.kv. 3) // ('DIMone'.kv. 1)
    call ncdf_crt(ncdf,d)
    call delete(d)
    
    ! Create all necessary containers...
    dic = ('info'.kv.'Number of supercells') 
    call ncdf_def_var(ncdf,'nsc',NF90_INT,(/'xyz'/), &
         atts=dic)

    dic = dic//('info'.kv.'Last orbital of equivalent atom')
    call ncdf_def_var(ncdf,'lasto',NF90_INT,(/'na_u'/), &
         atts=dic)
    
    dic = dic//('info'.kv.'Total charge')
    call ncdf_def_var(ncdf,'Qtot',NF90_DOUBLE,(/'one'/), &
         atts=dic)

    dic = dic//('info'.kv.'Fermi level')//('unit'.kv.'Ry')
    call ncdf_def_var(ncdf,'Ef',NF90_DOUBLE,(/'one'/), &
         atts=dic)
    
    dic = dic//('info'.kv.'Atomic coordinates')//('unit'.kv.'Bohr')
    call ncdf_def_var(ncdf,'xa',NF90_DOUBLE,(/'xyz ','na_u'/), &
         atts=dic)
    
    dic = dic//('info'.kv.'Unit cell')
    call ncdf_def_var(ncdf,'cell',NF90_DOUBLE,(/'xyz','xyz'/), &
         atts=dic)
    
    dic = dic//('info'.kv.'Atomic forces')//('unit'.kv.'Ry/Bohr')
    call ncdf_def_var(ncdf,'fa',NF90_DOUBLE,(/'xyz ','na_u'/), &
         atts=dic)
    
    dic = dic//('info'.kv.'Cell stress')//('unit'.kv.'Ry/Bohr**3')
    call ncdf_def_var(ncdf,'stress',NF90_DOUBLE,(/'xyz','xyz'/), &
         atts=dic)
    call delete(dic)
    
    ! Create matrix group
    call ncdf_def_grp(ncdf,'SPARSE',grp)

    call ncdf_def_dim(grp,'nnzs',n_nzs)

    dic = dic//('info'.kv.'Number of non-zero elements per row')
    call ncdf_def_var(grp,'n_col',NF90_INT,(/'no_u'/), &
         compress_lvl=cdf_comp_lvl,atts=dic)
    
    dic = dic//('info'.kv. &
         'Supercell column indices in the sparse format ')
    call ncdf_def_var(grp,'list_col',NF90_INT,(/'nnzs'/), &
         compress_lvl=cdf_comp_lvl,atts=dic)
    
    ! Create the overlap matrix (we know it will not change)
    dic = ('info'.kv.'Overlap matrix')
    call ncdf_def_var(grp,'S',NF90_DOUBLE,(/'nnzs'/), &
         compress_lvl=cdf_comp_lvl,atts=dic)
    
    dic = dic//('info'.kv.'Density matrix')
    call ncdf_def_var(grp,'DM',NF90_DOUBLE,(/'nnzs','spin'/), &
         compress_lvl=cdf_comp_lvl,atts=dic)
       
    dic = dic//('info'.kv.'Energy density matrix')
    dic = dic//('unit'.kv.'Ry')
    call ncdf_def_var(grp,'EDM',NF90_DOUBLE,(/'nnzs','spin'/), &
         compress_lvl=cdf_comp_lvl,atts=dic)
    
    dic = dic//('info'.kv.'Hamiltonian')
    call ncdf_def_var(grp,'H',NF90_DOUBLE,(/'nnzs','spin'/), &
         compress_lvl=cdf_comp_lvl,atts=dic)

    ! Even though I think we could do without, I add the
    ! xij array to the file
    if ( .not. Gamma ) then
       dic = dic//('info'.kv. &
            'Distance between orbital i and j')
       dic = dic//('unit'.kv.'Bohr')
       call ncdf_def_var(grp,'xij',NF90_DOUBLE,(/'xyz ','nnzs'/), &
            compress_lvl=cdf_comp_lvl,atts=dic)
    end if
    
    ! Delete the dictionary
    call delete(dic)

    call ncdf_def_grp(ncdf,'SETTINGS',grp)

    dic = ('info'.kv.'Tolerance for converging the density matrix')
    call ncdf_def_var(grp,'DMTolerance',NF90_DOUBLE,(/'one'/), &
         atts=dic)

    dic = dic//('info'.kv.'Tolerance for converging the Hamiltonian')
    call ncdf_def_var(grp,'HTolerance',NF90_DOUBLE,(/'one'/), &
         atts=dic)

    dic = dic//('info'.kv.'Net charge of the system')
    call ncdf_def_var(grp,'NetCharge',NF90_DOUBLE,(/'one'/), &
         atts=dic)

    dic = dic//('info'.kv.'Mixing weight')
    call ncdf_def_var(grp,'MixingWeight',NF90_DOUBLE,(/'one'/), &
         atts=dic)

    dic = dic//('info'.kv.'Grid used for the Brillouin zone integration')
    call ncdf_def_var(grp,'BZ',NF90_INT,(/'xyz','xyz'/), &
         atts=dic)

    dic = dic//('info'.kv.'Grid displacement used in Brillouin zone') &
         //('unit'.kv.'b**-1')
    call ncdf_def_var(grp,'BZ_displ',NF90_DOUBLE,(/'xyz'/), &
         atts=dic)

    dic = dic//('info'.kv.'Temperature for electrons')// &
         ('unit'.kv.'Ry')
    call ncdf_def_var(grp,'ElectronicTemperature',NF90_DOUBLE,(/'one'/), &
         atts=dic)

    dic = dic//('info'.kv.'Mesh cutoff for real space grid')
    call ncdf_def_var(grp,'MeshCutoff',NF90_DOUBLE,(/'one'/), &
         atts=dic)

    ! As enddef is collective
    ! we force it know
    call ncdf_enddef(ncdf)

    call ncdf_put_var(ncdf,'Qtot',Qtot)
    ! Save lasto
    call ncdf_put_var(ncdf,'lasto',lasto(1:na_u))
    
    ! Create matrix group
    if ( is_MD ) then

       ! ensure redefinition
       call ncdf_redef(ncdf)
       
       call ncdf_def_grp(ncdf,'MD',grp)

       ! Create arrays for containers
       call ncdf_def_dim(grp,'MD',NF90_UNLIMITED)

       dic = ('info'.kv.'Temperature')//('unit'.kv.'K')
       call ncdf_def_var(grp,'T',NF90_DOUBLE,(/'MD'/), &
            compress_lvl=0,atts=dic) ! do not compress unlimited D

       dic = dic//('info'.kv.'Kohn-Sham energy')//('unit'.kv.'Ry')
       call ncdf_def_var(grp,'E_KS',NF90_DOUBLE,(/'MD'/), &
            compress_lvl=0,atts=dic) ! do not compress unlimited D

       dic = dic//('info'.kv.'Total energy')
       call ncdf_def_var(grp,'E_tot',NF90_DOUBLE,(/'MD'/), &
            compress_lvl=0,atts=dic) ! do not compress unlimited D

       dic = dic//('info'.kv.'Kinetic energy')
       call ncdf_def_var(grp,'E_kin',NF90_DOUBLE,(/'MD'/), &
            compress_lvl=0,atts=dic) ! do not compress unlimited D

       dic = dic//('info'.kv.'Fermi level')
       call ncdf_def_var(grp,'Ef',NF90_DOUBLE,(/'MD'/), &
            compress_lvl=0,atts=dic) ! do not compress unlimited D

       dic = dic//('info'.kv.'Atomic coordinates')//('unit'.kv.'Bohr')
       call ncdf_def_var(grp,'xa',NF90_DOUBLE,(/'xyz ','na_u','MD  '/), &
            compress_lvl=0,atts=dic) ! do not compress unlimited D

       dic = dic//('info'.kv.'Atomic forces')//('unit'.kv.'Ry/Bohr')
       call ncdf_def_var(grp,'fa',NF90_DOUBLE,(/'xyz ','na_u','MD  '/), &
            compress_lvl=0,atts=dic) ! do not compress unlimited D

       dic = dic//('info'.kv.'Atomic velocities')//('unit'.kv.'Bohr/fs')
       call ncdf_def_var(grp,'va',NF90_DOUBLE,(/'xyz ','na_u','MD  '/), &
            compress_lvl=0,atts=dic) ! do not compress unlimited D

       if ( parallel_io(grp) ) then
          ! All are accessed independently
          ! The sparse routines force them to be COLLECTIVE
          call ncdf_default(grp,access=NF90_INDEPENDENT)
       end if

       ! end definition
       call ncdf_enddef(ncdf)

    end if

    call delete(dic)

    dic = ('time'.kv.datestring())
    dic = dic//('name'.kv.trim(sname))
    dic = dic//('label'.kv.trim(slabel))

    if ( isolve == SOLVE_DIAGON ) then
       dic = dic//('method'.kv.'diagon')
    else if ( isolve == SOLVE_ORDERN ) then
       dic = dic//('method'.kv.'order-n')
#ifdef TRANSIESTA
    else if ( isolve == SOLVE_TRANSI ) then
       dic = dic//('method'.kv.'transiesta')
#endif
    end if

    call ncdf_put_gatt(ncdf,atts=dic)

    call delete(dic)

#ifdef TRANSIESTA
    if ( isolve == SOLVE_TRANSI ) then
       ! Save all information about the transiesta method
       call ncdf_def_grp(ncdf,'TRANSIESTA',grp)

       
       ! Add all the electrodes
       do iEl = 1 , N_Elec
          
          call ncdf_def_grp(grp,trim(Elecs(iEl)%name),grp2)

          tmp = TotUsedAtoms(Elecs(iEl))
          call ncdf_def_dim(grp2,'na',tmp)

          allocate(ibuf(tmp))
          do i = 1 , tmp
             ibuf(i) = Elecs(iEl)%idx_a + i - 1
          end do
          dic = dic//('info'.kv.'Atoms belonging to electrode')
          call ncdf_def_var(grp2,'a_idx',NF90_INT,(/'na'/), &
               compress_lvl=0,atts=dic) ! do not compress unlimited D
          call ncdf_put_var(grp2,'a_idx',ibuf)
          deallocate(ibuf)
          
          
          dic = dic//('info'.kv.'Chemical potential')//('unit'.kv.'Ry')
          call ncdf_def_var(grp2,'mu',NF90_DOUBLE,(/'one'/), &
               compress_lvl=0,atts=dic) ! do not compress unlimited D

          call ncdf_put_var(grp2,'mu',Elecs(iEl)%mu%mu)

          call delete(dic)
       
       end do

       dic = dic//('info'.kv.'Bias')//('unit'.kv.'Ry')
       call ncdf_def_var(grp,'Volt',NF90_DOUBLE,(/'one'/), &
            compress_lvl=0,atts=dic) ! do not compress unlimited D

       call ncdf_put_var(grp,'Volt',Volt)

    end if

    call delete(dic)
#endif

    ! Close the file
    call ncdf_close(ncdf)
    
  end subroutine cdf_init_file

  subroutine cdf_save_settings(fname)

    use kpoint_grid, only : kscell, kdispl
    use siesta_options, only : cdf_w_parallel
    use siesta_options, only : dDtol, dHtol, charnet, wmix, temp, g2cut

    character(len=*), intent(in) :: fname
    
    type(hNCDF) :: ncdf, grp

    ! We just open it (prepending)
    if ( Nodes > 1 .and. cdf_w_parallel ) then
       call ncdf_open(ncdf,fname, groupname='SETTINGS', &
            mode=ior(NF90_WRITE,NF90_MPIIO), parallel = .true., &
            comm=MPI_Comm_World)
    else
       call ncdf_open(ncdf,fname, groupname='SETTINGS', &
            mode=ior(NF90_WRITE,NF90_NETCDF4))
    end if

    ! Write settings in settings group
    call ncdf_redef(ncdf)
    
    ! Save settings
    call ncdf_put_var(ncdf,'BZ',kscell)
    call ncdf_put_var(ncdf,'BZ_displ',kdispl)
    call ncdf_put_var(ncdf,'DMTolerance',dDtol)
    call ncdf_put_var(ncdf,'HTolerance',dHtol)
    call ncdf_put_var(ncdf,'NetCharge',charnet)
    call ncdf_put_var(ncdf,'MixingWeight',wmix)
    call ncdf_put_var(ncdf,'ElectronicTemperature',Temp)
    call ncdf_put_var(ncdf,'MeshCutoff',g2cut)

    ! I suggest that all MD settings are 
    ! put in the MD-group

    call ncdf_close(ncdf)

  end subroutine cdf_save_settings


  subroutine cdf_save_state(fname,dic_save)
    use m_gamma, only : Gamma
    use m_energies, only: Ef
    use atomlist, only : Qtot
    use siesta_options, only : cdf_w_parallel
    use siesta_geom, only: na_u, ucell, xa, va, nsc
    use sparse_matrices, only: sparse_pattern, block_dist
    use sparse_matrices, only: S_1D, DM_2D, EDM_2D, xij_2D, H_2D
    use m_stress, only : stress
    use m_forces, only: fa

    character(len=*), intent(in) :: fname
    ! Dictionary containing keys that we will save
    type(dict), intent(in) :: dic_save
    type(hNCDF) :: ncdf, grp

    call timer('CDF',1)

    ! We just open it (prepending)
    if ( Nodes > 1 .and. cdf_w_parallel ) then
       call ncdf_open(ncdf,fname, &
            mode=ior(NF90_WRITE,NF90_MPIIO), parallel = .true., &
            comm=MPI_Comm_World)
    else
       call ncdf_open(ncdf,fname,mode=ior(NF90_WRITE,NF90_NETCDF4))
    end if

    ! Add attribute of current Fermi-level
    if ( ('Ef' .in. dic_save) .and. Node == 0 ) &
         call ncdf_put_var(ncdf,'Ef',Ef)
    if ( ('Qtot' .in. dic_save) .and. Node == 0 ) &
         call ncdf_put_var(ncdf,'Qtot',Qtot)
    ! Save nsc, xa, fa, lasto, ucell
    if ( ('nsc' .in. dic_save) .and. Node == 0 ) &
         call ncdf_put_var(ncdf,'nsc',nsc)
    if ( ('cell' .in. dic_save) .and. Node == 0 ) &
         call ncdf_put_var(ncdf,'cell',ucell)
    if ( ('xa' .in. dic_save) .and. Node == 0 ) &
         call ncdf_put_var(ncdf,'xa',xa(:,1:na_u))
    if ( ('va' .in. dic_save) .and. Node == 0 ) &
         call ncdf_put_var(ncdf,'va',va(:,1:na_u))
    if ( ('fa' .in. dic_save) .and. Node == 0 ) &
         call ncdf_put_var(ncdf,'fa',fa(:,1:na_u))
    if ( ('stress' .in. dic_save) .and. Node == 0 ) &
         call ncdf_put_var(ncdf,'stress',stress)

    ! Sparsity format
    call ncdf_open_grp(ncdf,'SPARSE',grp)

    if ( 'sp' .in. dic_save ) &
         call cdf_w_Sp(grp,block_dist,sparse_pattern)
    if ( 'S' .in. dic_save ) &
         call cdf_w_d1D(grp,'S',S_1D)
    if ( .not. Gamma .and. ('xij' .in. dic_save) ) then
       ! Write the xij array, it will not change during SCF
       call cdf_w_d2D(grp,'xij',xij_2D)
    end if
    if ( 'H' .in. dic_save ) &
         call cdf_w_d2D(grp,'H',H_2D)
    if ( 'DM' .in. dic_save ) &
         call cdf_w_d2D(grp,'DM',DM_2D)
    if ( 'EDM' .in. dic_save ) &
         call cdf_w_d2D(grp,'EDM',EDM_2D)

    call ncdf_close(ncdf)
    
    call timer('CDF',2)

  end subroutine cdf_save_state

  subroutine cdf_save_md(fname)
    use dictionary
    use siesta_geom, only: na_u, xa, va
    use m_forces, only: fa
    use m_energies, only: Etot, Ef, Ekinion
    use m_kinetic, only : Tempion

    character(len=*), intent(in) :: fname

    type(hNCDF) :: ncdf
    integer :: MD

    ! open the file...
    call ncdf_open(ncdf,fname,mode=ior(NF90_WRITE,NF90_NETCDF4),groupname='MD')

    ! Inquire the current size of the MD-variable
    MD = 0
    call ncdf_inq_dim(ncdf,'MD',len=MD)
    MD = MD + 1 ! increment

    if ( Node == 0 ) then
       ! Save current state
       call ncdf_put_var(ncdf,'T',Tempion,start=(/MD/))
       call ncdf_put_var(ncdf,'Ef',Ef,start=(/MD/))
       !call ncdf_put_var(ncdf,'E_KS',EKS,start=(/MD/))
       call ncdf_put_var(ncdf,'E_tot',Etot,start=(/MD/))
       call ncdf_put_var(ncdf,'E_kin',Ekinion,start=(/MD/))
       
       call ncdf_put_var(ncdf,'xa',xa,count=(/3,na_u/),start=(/1,1,MD/))
       call ncdf_put_var(ncdf,'fa',fa,count=(/3,na_u/),start=(/1,1,MD/))
       call ncdf_put_var(ncdf,'va',va,count=(/3,na_u/),start=(/1,1,MD/))
    end if

    call ncdf_close(ncdf)

  end subroutine cdf_save_md

  subroutine cdf_w_Sp(ncdf,dit,sp)
    use m_io_s,only : Node_Sp_gncol
    use class_OrbitalDistribution
    use class_Sparsity
  
    type(hNCDF), intent(inout) :: ncdf
    type(OrbitalDistribution), intent(inout) :: dit
    type(Sparsity), intent(inout) :: sp

    ! Local variables
    integer, pointer :: ncol(:), l_col(:)
    integer, allocatable :: buf(:)
    integer :: no_l, no_u, n_nzs, gio, io, ind, gind, max_n
#ifdef MPI
    integer :: BNode, MPIerror, MPIstatus(MPI_STATUS_SIZE), MPIreq
#endif
    integer :: n_nnzs, g_nzs

    ! Write the sparsity to the file...
    call attach(sp,nrows=no_l, nrows_g=no_u, &
         n_col=ncol,list_col=l_col,nnzs=n_nzs)

#ifdef MPI
    call MPI_Reduce(n_nzs,g_nzs,1,MPI_Integer, &
         MPI_Sum, 0, MPI_Comm_World, MPIerror)
#else
    g_nzs = nnzs(sparse_pattern)
#endif

    ! Read dimension nnzs
    call ncdf_inq_dim(ncdf,'nnzs',len=n_nnzs)
    if ( Node == 0 .and. n_nnzs /= g_nzs ) then
       call die('Number of non-zero elements is not equivalent.')
    end if
    
    if ( parallel_io(ncdf) ) then

       print *,'Here:',node,' this isformatted erroneously'

       ! Calculate the offset for the global
       ! partitioning
       gind = global_offset(dit,no_l)

       ! we write it using MPI
       call ncdf_par_access(ncdf,name='n_col',access=NF90_COLLECTIVE)
       call ncdf_put_var(ncdf,'n_col',ncol,start=(/gind+1/), &
            count=(/no_l/))

       gind = global_offset(dit,n_nnzs)

       ! we write it using MPI
       call ncdf_par_access(ncdf,name='list_col',access=NF90_COLLECTIVE)
       call ncdf_put_var(ncdf,'list_col',l_col,start=(/gind+1/), &
            count=(/n_nzs/))

    else

#ifdef MPI
       allocate(buf(no_u))
       call Node_Sp_gncol(0,sp,dit,no_u,buf)
       call ncdf_put_var(ncdf,'n_col',buf)
#else
       call ncdf_put_var(ncdf,'n_col',ncol)
#endif

       if ( Node == 0 ) then
#ifdef MPI
          max_n = maxval(buf)
          deallocate(buf)
#else
          max_n = maxval(ncol)
#endif
          allocate(buf(max_n))
       end if

       ! Write list_col
#ifdef MPI

       MPIreq = MPI_REQUEST_NULL

       ! Loop size
       ind = 0
       gind = 1
       do gio = 1 , no_u
          BNode = node_handling_element(dit,gio)

          if ( Node == BNode ) then
             io = index_global_to_local(dit,gio,Node)
             if ( Node == 0 ) then
                call ncdf_put_var(ncdf,'list_col',l_col(ind+1:ind+ncol(io)), &
                     count=(/ncol(io)/),start=(/gind/))
                gind = gind + ncol(io)
             else
                call MPI_Send( l_col(ind+1) , ncol(io), MPI_Integer, &
                     0, gio, MPI_Comm_World, MPIerror)
             end if
             ind = ind + ncol(io)
          else if ( Node == 0 ) then
             call MPI_Recv( buf(1) , max_n, MPI_Integer, &
                  BNode, gio, MPI_Comm_World, MPIstatus, MPIerror )
             if ( MPIerror /= MPI_Success ) &
                  call die('Error in code: ncdf_write_Sp')
             call MPI_Get_Count(MPIstatus, MPI_Integer, io, MPIerror)
             call ncdf_put_var(ncdf,'list_col',buf(1:io), &
                  count=(/io/),start=(/gind/))
             gind = gind + io
          end if
       end do

       if ( Node == 0 ) then
          deallocate(buf)
       end if

#else
       call ncdf_put_var(ncdf,'list_col',l_col)
#endif

    end if

  end subroutine cdf_w_Sp

  subroutine cdf_w_d1D(ncdf,vname,dSp1D)
    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData1D

    type(hNCDF), intent(inout) :: ncdf
    character(len=*), intent(in) :: vname
    type(dSpData1D), intent(inout) :: dSp1D

    ! Local variables
    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: sp
    integer, pointer :: ncol(:), l_col(:)
    integer :: no_l, no_u, n_nzs, gio, io, ind, gind, max_n
    real(dp), pointer :: a(:)
    real(dp), allocatable :: buf(:)

#ifdef MPI
    integer :: BNode, MPIerror, MPIstatus(MPI_STATUS_SIZE), MPIreq
#endif

    dit => dist(dSp1D)
    sp => spar(dSp1D)
    ! Write the sparsity to the file...
    call attach(sp,nrows=no_l, nrows_g=no_u, &
         n_col=ncol,list_col=l_col,nnzs=n_nzs)

    a => val(dSp1D)
    
    if ( parallel_io(ncdf) ) then

       call ncdf_par_access(ncdf,name=vname,access=NF90_COLLECTIVE)

       ! Calculate the processors offset
       gind = global_offset(dit,n_nzs)

       ! we write it using MPI
       call ncdf_put_var(ncdf,vname,a,start=(/gind+1/), &
            count=(/n_nzs/))

    else

#ifdef MPI
       io = maxval(ncol)
       call MPI_Reduce(io,max_n,1,MPI_Integer, &
            MPI_Max, 0, MPI_Comm_World, MPIerror)
       if ( Node == 0 ) then
          allocate(buf(max_n))
       end if

       MPIreq = MPI_REQUEST_NULL

       ind = 0
       gind = 1
       do gio = 1 , no_u
          BNode = node_handling_element(dit,gio)
          
          if ( Node == BNode ) then
             io = index_global_to_local(dit,gio,Node)
             if ( Node == 0 ) then
                call ncdf_put_var(ncdf,vname,a(ind+1:ind+ncol(io)), &
                     count=(/ncol(io)/),start=(/gind/))
                gind = gind + ncol(io)
             else
                call MPI_Send( a(ind+1) , ncol(io), MPI_Double_Precision, &
                     0, gio, MPI_Comm_World, MPIerror)
             end if
             ind = ind + ncol(io)
          else if ( Node == 0 ) then
             call MPI_Recv( buf(1) , max_n, MPI_Double_Precision, &
                  BNode, gio, MPI_Comm_World, MPIstatus, MPIerror )
             if ( MPIerror /= MPI_Success ) &
                  call die('Error in code: ncdf_write_d1D')
             call MPI_Get_Count(MPIstatus, MPI_Double_Precision, io, MPIerror)
             call ncdf_put_var(ncdf,vname,buf(1:io), &
                  count=(/io/),start=(/gind/))
             gind = gind + io
          end if
       end do

#else
       call ncdf_put_var(ncdf,trim(vname),a)
#endif

    end if
    
  end subroutine cdf_w_d1D

  subroutine cdf_w_d2D(ncdf,vname,dSp2D)
    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData2D

    type(hNCDF), intent(inout) :: ncdf
    character(len=*), intent(in) :: vname
    type(dSpData2D), intent(inout) :: dSp2D

    ! Local variables
    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: sp
    integer, pointer :: ncol(:), l_col(:)
    integer :: no_l, no_u, n_nzs, gio, io, ind, gind, max_n
    integer :: id2, dim2, sp_dim
    real(dp), pointer :: a(:,:)
    real(dp), allocatable :: buf(:)

#ifdef MPI
    integer :: BNode, MPIerror, MPIstatus(MPI_STATUS_SIZE), MPIreq
#endif

    dit => dist(dSp2D)
    sp => spar(dSp2D)
    ! Write the sparsity to the file...
    call attach(sp,nrows=no_l, nrows_g=no_u, &
         n_col=ncol,list_col=l_col,nnzs=n_nzs)

    a => val(dSp2D)

    if ( size(a,dim=2) == n_nzs ) then
       sp_dim = 2
       dim2 = size(a,dim=1)
    else
       sp_dim = 1
       dim2 = size(a,dim=2)
    end if

    if ( parallel_io(ncdf) ) then

       call ncdf_par_access(ncdf,name=vname,access=NF90_COLLECTIVE)

       ! we write it using MPI
       gind = global_offset(dit,n_nzs)

       if ( sp_dim == 1 ) then
          do io = 1 , dim2
             call ncdf_put_var(ncdf,trim(vname),a(:,io), &
                  start=(/gind+1,io/), count=(/n_nzs/))
          end do
       else
          call ncdf_put_var(ncdf,trim(vname),a,start=(/1,gind+1/), &
               count=(/dim2,n_nzs/))
       end if

    else

#ifdef MPI
       io = maxval(ncol)
       call MPI_Reduce(io,max_n,1,MPI_Integer, &
            MPI_Max, 0, MPI_Comm_World, MPIerror)
       max_n = max_n * dim2
       if ( Node == 0 ) then
          allocate(buf(max_n))
       end if

       MPIreq = MPI_REQUEST_NULL

    if ( sp_dim == 1 ) then

       do id2 = 1 , dim2
          ind = 0
          gind = 1
          do gio = 1 , no_u
             BNode = node_handling_element(dit,gio)
             
             if ( Node == BNode ) then
                io = index_global_to_local(dit,gio,Node)
                if ( Node == 0 ) then
                   call ncdf_put_var(ncdf,trim(vname),a(ind+1:ind+ncol(io),id2), &
                        count=(/ncol(io)/),start=(/gind,id2/))
                   gind = gind + ncol(io)
                else
                   call MPI_Send( a(ind+1,id2) , ncol(io), MPI_Double_Precision, &
                        0, gio, MPI_Comm_World, MPIerror)
                   
                end if
                ind = ind + ncol(io)
             else if ( Node == 0 ) then
                call MPI_Recv( buf(1) , max_n, MPI_Double_Precision, &
                     BNode, gio, MPI_Comm_World, MPIstatus, MPIerror )
                if ( MPIerror /= MPI_Success ) &
                     call die('Error in code (1): ncdf_write_d2D')
                call MPI_Get_Count(MPIstatus, MPI_Double_Precision, io, MPIerror)
                call ncdf_put_var(ncdf,trim(vname),buf(1:io), &
                     count=(/io/),start=(/gind,id2/))
                gind = gind + io
             end if
          end do ! gio
       end do ! id2

    else

       ind = 0
       gind = 1
       do gio = 1 , no_u
          BNode = node_handling_element(dit,gio)
          
          if ( Node == BNode ) then
             io = index_global_to_local(dit,gio,Node)
             if ( Node == 0 ) then
                call ncdf_put_var(ncdf,trim(vname),a(1:dim2,ind+1:ind+ncol(io)), &
                     count=(/dim2,ncol(io)/),start=(/1,gind/))
                gind = gind + ncol(io)
             else
                call MPI_Send( a(1,ind+1) , dim2*ncol(io), MPI_Double_Precision, &
                     0, gio, MPI_Comm_World, MPIerror)
             end if
             ind = ind + ncol(io)
          else if ( Node == 0 ) then
             call MPI_Recv( buf(1) , max_n, MPI_Double_Precision, &
                  BNode, gio, MPI_Comm_World, MPIstatus, MPIerror )
             if ( MPIerror /= MPI_Success ) &
                  call die('Error in code: ncdf_write_d2D')
             call MPI_Get_Count(MPIstatus, MPI_Double_Precision, io, MPIerror)
             io = io / dim2
             call ncdf_put_var(ncdf,trim(vname),reshape(buf(1:io*dim2),(/dim2,io/)), &
                  count=(/dim2,io/),start=(/1,gind/))
             gind = gind + io
          end if
       end do

    end if
#else
       call ncdf_put_var(ncdf,trim(vname),a)
#endif

    end if
    
  end subroutine cdf_w_d2D
  
#else
  subroutine dummy()
    !pass
  end subroutine dummy
#endif

end module m_ncdf_io
