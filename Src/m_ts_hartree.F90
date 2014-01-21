!
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
!
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
! This code segment has been improved or fully created by:
! Nick Papior Andersen, 2014, nickpapior@gmail.com
!
module m_ts_hartree

! Module for fixing the Hartree potential so that the potential fluctuations
! does not go wild.
! This is necessary to get a stable SCF solution
!
! Created and copyrighted by: Nick Papior Andersen, 2014
! The use of this program is allowed for not-for-profit research only.
! Copy or disemination of all or part of this package is not
! permitted without prior and explicit authorization by the author.
  
  use precision, only : dp

  implicit none
  
  private

  ! The idea is to have sub routines in this module to do
  ! various Hartree potential fixes
  public :: ts_init_hartree_fix
  public :: ts_hartree_fix

  ! We construct the square at which we fix the potential

  ! The lower-left corner and vectors spanning the square
  real(dp) :: ll_c(3), v1(3), v2(3)
  ! The normal-vector
  real(dp) :: n(3)
  ! an auxillary length to ease computations
  real(dp) :: d

contains

  subroutine ts_init_hartree_fix(ucell,na_u,xa,meshG,nsm)
    use mesh, only : cmesh
! ***********************
! * INPUT variables     *
! ***********************
    real(dp),      intent(in) :: ucell(3,3)
    integer,       intent(in) :: na_u
    real(dp),      intent(in) :: xa(3,na_u)
    integer,       intent(in) :: meshG(3), nsm

    ! currently we don't do anything...

  end subroutine ts_init_hartree_fix

  ! Fix the potential
  subroutine ts_hartree_fix( Vscf )
    use precision, only : dp, grid_p
    use sys, only : die
    use parallel, only : Node, Nodes
#ifdef MPI
    use mpi_siesta, only : MPI_AllReduce, MPI_Sum
    use mpi_siesta, only : MPI_Comm_World, MPI_integer
    use mpi_siesta, only : MPI_double_precision
#endif
    use m_ts_tdir
    use m_ts_mesh, only : meshl, offset_i
    
    real(grid_p), intent(inout) :: Vscf(:,:)

! Internal variables
    integer :: i1, i2, i3, imesh, ntemp
    integer :: nlp
    integer, target  :: i10, i20, i30
    integer, pointer :: iT
#ifdef MPI
    integer :: MPIerror, npl
#endif
    real(dp) :: Vav, Vtot, temp

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE TSVHfix' )
#endif

    ! Set up counter
    if ( ts_tdir == 1 ) then
       iT => i10
    else if ( ts_tdir == 2 ) then
       iT => i20
    else
       iT => i30
    end if

    Vtot = 0._dp
    nlp  = 0

    ! Test whether we should do anything (note that iT => [i10|i20|i30]):
    i10 = 0
    i20 = offset_i(2) - 1
    i30 = offset_i(3) - 1
    if ( iT <= 0 ) then
       imesh = 0
       i30 = offset_i(3) - 1
       do i3 = 0,meshl(3)-1
          i30 = i30 + 1
          i20 = offset_i(2) - 1
          do i2 = 0,meshl(2)-1
             i20 = i20 + 1
             do i10 = 0,meshl(1)-1
                imesh = imesh + 1
                if (iT.eq.0) then
                   nlp = nlp + 1
                   Vtot = Vtot + Vscf(imesh,1)
                end if
             end do
          end do
       end do
    end if

#ifdef MPI
    call MPI_AllReduce(Vtot,temp,1,MPI_double_precision,MPI_Sum, &
         MPI_Comm_World,MPIerror)
    Vtot = temp
    call MPI_AllReduce(nlp,ntemp,1,MPI_integer,MPI_Sum, &
         MPI_Comm_World,MPIerror)
    nlp = ntemp
#endif

    Vav = Vtot/real(nlp,dp)

    imesh = 0
    do i30 = 1 , meshl(3)
       do i20 = 1 , meshl(2)
          do i10 = 1 , meshl(1)
             imesh = imesh + 1
             Vscf(imesh,1) = Vscf(imesh,1) - Vav
          end do
       end do
    end do

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS TSVHfix' )
#endif

  end subroutine ts_hartree_fix

end module m_ts_hartree
