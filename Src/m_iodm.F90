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
! This module has been rewritten to conform to Nick Papior Andersens IO routines
! It has also been coded by Nick Papior Andersen
module m_iodm

  use parallel, only : Node

  use class_Sparsity
  use class_OrbitalDistribution
  use class_dSpData2D

  use m_io_s

  implicit none
  
  private
  public :: write_dm, read_dm

contains
        
  subroutine read_dm( file, nspin, dit, no_u, DM, found , Bcast )

#ifdef MPI
    use mpi_siesta
#endif

! **********************
! * INPUT variables    *
! **********************
    character(len=*), intent(in) :: file
    integer, intent(in) :: nspin
    type(OrbitalDistribution), intent(inout) :: dit
    integer, intent(in) :: no_u
    type(dSpData2D), intent(inout) :: DM
    logical, intent(out) :: found
    logical, intent(in), optional :: Bcast

! ************************
! * LOCAL variables      *
! ************************
    type(Sparsity) :: sp
    character(len=500) :: fn
    logical :: exists, lBcast
    integer :: iu, two(2)
    integer, allocatable, target :: gncol(:)
#ifdef MPI
    integer :: MPIerror
#endif

    external :: io_assign, io_close

    if ( Node == 0 ) then
       inquire(file=file, exist=exists)
    end if

    lBcast = .false.
    if ( present(Bcast) ) lBcast = Bcast

#ifdef MPI
    call MPI_Bcast(exists,1,MPI_Logical,0, &
         MPI_Comm_World,MPIerror)
#endif

    fn = 'IO-DM: '//trim(file)
    
    found = .false.
    if ( .not. exists ) return

    if ( Node == 0 ) then
       call io_assign(iu)
       open( iu, file=file, form='unformatted', status='old' )
       rewind(iu)
       read(iu) two
    end if

#ifdef MPI
    call MPI_Bcast(two,2,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif

    ! If the total number of orbitals does not match, bail out
    if ( no_u /= two(1) .or. nspin /= two(2) ) then
       if ( Node == 0 ) then
          write(*,"(a,i6,/,a)") &
               "WARNING: Wrong number of orbitals in DM file: ",two(1), &
               "WARNING: Falling back to atomic initialization of DM."
          write(0,"(a,i6,/,a)") &
               "WARNING: Wrong number of orbitals in DM file: ",two(1), &
               "WARNING: Falling back to atomic initialization of DM."
          call io_close(iu)
       endif

       found = .false.

       return
    end if

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

    ! Clean-up (sp is not fully deleted, it just only resides in DM)
    call delete(sp)

    ! All have this allocated (Node == 0 have just a larger
    ! one...)
    deallocate(gncol)

    ! Close
    if ( Node == 0 ) then

       call io_close(iu)

    end if

    found = .true.
    
  end subroutine read_dm
  
  subroutine write_dm(file, nspin, DM )
    
! **********************
! * INPUT variables    *
! **********************
    character(len=*), intent(in) :: file
    integer, intent(in) :: nspin
    type(dSpData2D), intent(inout) :: DM
    
! ************************
! * LOCAL variables      *
! ************************
    type(Sparsity), pointer :: sp
    type(OrbitalDistribution), pointer :: dit
    integer, allocatable, target :: gncol(:)
    integer :: no_u
    integer :: iu

    external :: io_assign, io_close
    
    ! Gather sparse pattern
    dit => dist(DM)
    sp => spar(DM)
    call attach(sp,nrows_g=no_u)
    
    if ( Node == 0 ) then

       ! Open file
       call io_assign( iu )
       open( iu, file=file, &
            form='unformatted', status='unknown' )
       rewind(iu)
       
       write(iu) no_u, nspin

    end if

    allocate(gncol(no_u))
    gncol(1) = -1

    ! Write sparsity pattern...
    call io_write_Sp(iu,sp,dit=dit,gncol=gncol)

    ! Write density matrix
    call io_write_d2D(iu,DM,gncol=gncol)

    deallocate(gncol)

    ! Close
    if ( Node == 0 ) then
       
       call io_close(iu)
       
    end if
    
  end subroutine write_dm

end module m_iodm
