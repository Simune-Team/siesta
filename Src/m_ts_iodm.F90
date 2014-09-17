!     
!     This file is part of the SIESTA package.
!     
!     Copyright (c) Fundacion General Universidad Autonoma de Madrid:
!     E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
!     and J.M.Soler, 1996- .
!     
!     Use of this software constitutes agreement with the full conditions
!     given in the SIESTA license, as signed by all legitimate users.
!    

! Fully created by Nick Papior Andersen to conform with the io_s
! library.
module m_ts_iodm
  
  use precision, only: dp
  use parallel, only : Node, Nodes
#ifdef MPI
  use mpi_siesta
#endif
  
  implicit none
  
  private
  public :: ts_init_dm
  public :: write_ts_dm, read_ts_dm

  ! We do not allow formatted writes now...

contains
  
  subroutine read_ts_dm ( slabel, nspin, dit, sp_def, DM, EDM, Ef, found )

    use class_Sparsity
    use class_OrbitalDistribution
    use class_dSpData2D
    use m_io_s

! **********************
! * INPUT variables    *
! **********************
    character(len=*), intent(in) :: slabel
    integer, intent(in) :: nspin
    type(OrbitalDistribution), intent(in) :: dit
    type(Sparsity), intent(inout) :: sp_def
    type(dSpData2D), intent(inout) :: DM, EDM
    real(dp), intent(inout) :: Ef
    logical, intent(out) :: found

! ************************
! * LOCAL variables      *
! ************************
    type(Sparsity) :: sp
    logical :: exists
    integer :: iu, io, n_nzs, no_u
    integer, dimension(:), pointer :: oncol, ol_col
    integer, dimension(:), pointer :: ncol, l_col
#ifdef MPI
    integer :: MPIerror
#endif

    external :: io_assign, io_close

    if ( Node == 0 ) then
       inquire(file=trim(slabel)//'.TSDE', exist=exists)
    end if

#ifdef MPI
    call MPI_Bcast(exists,1,MPI_Logical,0, &
         MPI_Comm_World,MPIerror)
#endif
    
    found = .false.
    if ( .not. exists ) return

    call attach(sp_def,nrows_g=no_u)

    if ( Node == 0 ) then
       write(*,'(/,a)') 'ts_iodm: Reading Density Matrix from files'
       call io_assign(iu)
       open( iu, file=trim(slabel)//'.TSDE', &
            form='unformatted', status='old' )
       rewind(iu)
       read(iu) io, n_nzs
    end if

#ifdef MPI
    call MPI_Bcast(io,1,MPI_integer,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(n_nzs,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif

    ! This checks for equality of spin and number of rows.
    call chkdim( 'iodm', 'no_u', no_u, io, 0 )
    call chkdim( 'iodm', 'nspin',  nspin,  n_nzs, 0 )
    
    ! Read in the sparsity pattern (distributed)
    call io_read_Sp(iu, no_u, sp, 'temp-TSIO', dit)

    ! Check number of sparsity dimensions
    call attach(sp_def,nnzs=n_nzs)
    call attach(sp    ,nnzs=io)
    call chkdim( 'iodm', 'nnzs',  n_nzs,  io, 0 )

    ! Read DM
    call io_read_d2D(iu,sp,DM,nspin,'ts-IODM',dit=dit)

    ! Read EDM
    call io_read_d2D(iu,sp,EDM,nspin,'ts-IODM',dit=dit)

    ! Clean-up
    call delete(sp)

    ! Read Ef and close
    if ( Node == 0 ) then

       read(iu) Ef

       call io_close(iu)

    end if

#ifdef MPI
    call MPI_BCast(Ef,1,MPI_Double_Precision,0, &
         MPI_Comm_World, MPIerror)
#endif

    found = .true.
    
  end subroutine read_ts_dm

  subroutine write_ts_dm(slabel,nspin, DM, EDM, Ef )
    
    use class_Sparsity
    use class_OrbitalDistribution
    use class_dSpData2D
    use m_io_s

! **********************
! * INPUT variables    *
! **********************
    character(len=*), intent(in) :: slabel
    integer, intent(in) :: nspin
    type(dSpData2D), intent(inout) :: DM, EDM
    real(dp), intent(in) :: Ef
    
! ************************
! * LOCAL variables      *
! ************************
    type(Sparsity), pointer :: sp
    type(OrbitalDistribution), pointer :: dit
    
    external :: io_assign, io_close
    
    integer :: no_u
    integer :: iu

    ! Gather sparse pattern
    dit => dist(DM)
    sp => spar(DM)
    call attach(sp,nrows_g=no_u)
    
    if ( Node == 0 ) then

       ! Open file
       call io_assign( iu )
       open( iu, file=trim(slabel)//'.TSDE', &
            form='unformatted', status='unknown' )
       rewind(iu)
       
       write(iu) no_u, nspin

    end if

    ! Write sparsity pattern...
    call io_write_Sp(iu,sp,dit=dit)

    ! Write density matrix
    call io_write_d2D(iu,DM)

    ! Write energy density matrix
    call io_write_d2D(iu,EDM)

    ! Write Ef and close
    if ( Node == 0 ) then
       
       write(iu) Ef
       
       call io_close(iu)
       
    end if
    
  end subroutine write_ts_dm
  
  subroutine ts_init_dm(found)
    
    use parallel, only: IOnode 
    use m_ts_global_vars, only : TSinit, TSrun
    
    logical, intent(in) :: found
    
    if( .not. found )  then ! not a TS continuation run
       TSinit = .true.  ! start to converge a diagon run
       TSrun = .false.
       
       if(IONode) then
          write(*,'(a)') 'TRANSIESTA: No TS-DensityMatrix file found'
          write(*,'(a)') 'TRANSIESTA: Initialization runs using diagon'
       end if
    else
       TSinit = .false.
       TSrun = .true.
       
       if ( IONode ) then
          write(*,'(a,/)') 'TRANSIESTA: Continuation run'
          write(*,'(a)') '                     ************************'
          write(*,'(a)') '                     *   TRANSIESTA BEGIN   *'
          write(*,'(a)') '                     ************************'
       end if
    end if !found TSDM-file
    
  end subroutine ts_init_dm
  
end module m_ts_iodm
