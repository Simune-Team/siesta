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
        
  subroutine read_dm ( slabel, nspin, dit, sp_def, DM, found )

#ifdef MPI
    use mpi_siesta
#endif

! **********************
! * INPUT variables    *
! **********************
    character(len=*), intent(in) :: slabel
    integer, intent(in) :: nspin
    type(OrbitalDistribution), intent(in) :: dit
    type(Sparsity), intent(inout) :: sp_def
    type(dSpData2D), intent(inout) :: DM
    logical, intent(out) :: found

! ************************
! * LOCAL variables      *
! ************************
    type(Sparsity) :: sp
    logical :: exists
    integer :: iu, no_u, two(2)
#ifdef MPI
    integer :: MPIerror
#endif

    external :: io_assign, io_close

    if ( Node == 0 ) then
       inquire(file=trim(slabel)//'.DM', exist=exists)
    end if

#ifdef MPI
    call MPI_Bcast(exists,1,MPI_Logical,0, &
         MPI_Comm_World,MPIerror)
#endif
    
    found = .false.
    if ( .not. exists ) return

    call attach(sp_def,nrows_g=no_u)

    if ( Node == 0 ) then
       write(*,'(/,a)') 'iodm: Reading Density Matrix from files'
       call io_assign(iu)
       open( iu, file=trim(slabel)//'.DM', &
            form='unformatted', status='old' )
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

    ! Read in the sparsity pattern (distributed)
    call io_read_Sp(iu, no_u, sp, 'temp-IO', dit)

    ! Read DM
    call io_read_d2D(iu,sp,DM,nspin,'iodm',dit=dit)

    ! Clean-up (sp is not fully deleted, it just only resides in DM)
    call delete(sp)

    ! Close
    if ( Node == 0 ) then

       call io_close(iu)

    end if

    found = .true.
    
  end subroutine read_dm
  
  subroutine write_dm(slabel, nspin, DM )
    
! **********************
! * INPUT variables    *
! **********************
    character(len=*), intent(in) :: slabel
    integer, intent(in) :: nspin
    type(dSpData2D), intent(inout) :: DM
    
! ************************
! * LOCAL variables      *
! ************************
    type(Sparsity), pointer :: sp
    type(OrbitalDistribution), pointer :: dit
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
       open( iu, file=trim(slabel)//'.DM', &
            form='unformatted', status='unknown' )
       rewind(iu)
       
       write(iu) no_u, nspin

    end if

    ! Write sparsity pattern...
    call io_write_Sp(iu,sp,dit=dit)

    ! Write density matrix
    call io_write_d2D(iu,DM)

    ! Close
    if ( Node == 0 ) then
       
       call io_close(iu)
       
    end if
    
  end subroutine write_dm

end module m_iodm
