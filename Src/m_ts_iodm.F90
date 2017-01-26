!     
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!    

! Fully created by Nick Papior Andersen to conform with the io_s
! library.
module m_ts_iodm
  
  use precision, only: dp
  use parallel, only : Node
  
  use class_Sparsity
  use class_OrbitalDistribution
  use class_dSpData2D

  use m_os, only : file_exist
  use m_io_s

  implicit none
  
  private
  public :: write_ts_dm, read_ts_dm

contains
  
  subroutine read_ts_dm( file, dit, DM, EDM, Ef, found, Bcast)

#ifdef MPI
  use mpi_siesta
#endif

! **********************
! * INPUT variables    *
! **********************
    character(len=*), intent(in) :: file
    type(OrbitalDistribution), intent(in) :: dit
    type(dSpData2D), intent(inout) :: DM, EDM
    real(dp), intent(inout) :: Ef
    logical, intent(out) :: found
    logical, intent(in), optional :: Bcast

! ************************
! * LOCAL variables      *
! ************************
    type(Sparsity) :: sp
    character(len=500) :: fn
    logical :: lBcast
    integer :: iu, two(2), no_u, nspin
    integer, allocatable, target :: gncol(:)
#ifdef MPI
    integer :: MPIerror
#endif

    external :: io_assign, io_close

    lBcast = .false.
    if ( present(Bcast) ) lBcast = Bcast

    found = file_exist(file, Bcast = Bcast )

    if ( .not. found ) then

       ! Clean-output
       call delete(DM)
       call delete(EDM)

       return

    end if

    ! Make name for the readed TSDE
    fn = 'IO-TSDE: '//trim(file)

    if ( Node == 0 ) then
       call io_assign(iu)
       open( iu, file=file, form='unformatted', status='old' )
       rewind(iu)
       read(iu) two
    end if
    
#ifdef MPI
    call MPI_Bcast(two,2,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif

    no_u = two(1)
    nspin = two(2)

    allocate(gncol(no_u))
    gncol(1) = 1
    
    ! Read in the sparsity pattern (distributed)
    if ( lBcast ) then
       call io_read_Sp(iu, no_u, sp, trim(fn), gncol=gncol, Bcast=Bcast)
    else
       call io_read_Sp(iu, no_u, sp, trim(fn), dit=dit, gncol=gncol)
    end if

    ! Read DM
    if ( lBcast ) then
       call io_read_d2D(iu,sp,DM ,nspin, trim(fn), gncol=gncol, Bcast=Bcast)
    else
       call io_read_d2D(iu,sp,DM ,nspin, trim(fn), dit=dit, gncol=gncol)
    end if

    ! Read EDM
    if ( lBcast ) then
       call io_read_d2D(iu,sp,EDM,nspin, trim(fn), gncol=gncol, Bcast=Bcast)
    else
       call io_read_d2D(iu,sp,EDM,nspin, trim(fn), dit=dit, gncol=gncol)
    end if

    ! Clean-up
    call delete(sp)

    deallocate(gncol)

    ! Read Ef and close
    if ( Node == 0 ) then

       read(iu) Ef

       call io_close(iu)

    end if

#ifdef MPI
    call MPI_BCast(Ef,1,MPI_Double_Precision, &
         0,MPI_Comm_World, MPIerror)
#endif

  end subroutine read_ts_dm

  subroutine write_ts_dm(file, DM, EDM, Ef )
    
! **********************
! * INPUT variables    *
! **********************
    character(len=*), intent(in) :: file
    type(dSpData2D), intent(inout) :: DM, EDM
    real(dp), intent(in) :: Ef
    
! ************************
! * LOCAL variables      *
! ************************
    type(Sparsity), pointer :: sp
    type(OrbitalDistribution), pointer :: dit
    integer, allocatable, target :: gncol(:)
    integer :: no_u, nspin
    integer :: iu

    external :: io_assign, io_close

    ! Gather sparse pattern
    dit => dist(DM)
    sp  => spar(DM)
    call attach(sp,nrows_g=no_u)
    ! Retrieve number of spin-components
    nspin = size(DM, 2)
    
    if ( Node == 0 ) then

       ! Open file
       call io_assign( iu )
       open( iu, file=file, form='unformatted', status='unknown' )
       rewind(iu)
       
       write(iu) no_u, nspin

    end if

    allocate(gncol(no_u))
    gncol(1) = -1

    ! Write sparsity pattern...
    call io_write_Sp(iu,sp, dit=dit, gncol=gncol)

    ! Write density matrix
    call io_write_d2D(iu, DM , gncol=gncol)

    ! Write energy density matrix
    call io_write_d2D(iu, EDM, gncol=gncol)

    deallocate(gncol)

    ! Write Ef and close
    if ( Node == 0 ) then
       
       write(iu) Ef
       
       call io_close(iu)
       
    end if
    
  end subroutine write_ts_dm
  
end module m_ts_iodm
