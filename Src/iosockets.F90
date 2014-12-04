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
module iosockets
! Handles communications between siesta-as-a-subroutine and a driver 
! program, transfering coordinates and forces trough POSIX sockets.
! The routines that handle the other side of the communication are
! in module fsiesta.
! modified from iopipes J.M.Soler and A.Garcia. Nov.2003
! M. Ceriotti 2014

use precision, only: dp
use parallel, only: IOnode
use fdf
use m_fdf_global, only: fdf_global_get
use f90sockets, ONLY : open_socket, writebuffer, readbuffer
use sys, only: die, bye
use  m_mpi_utils, only: broadcast
#ifdef MPI
      use mpi_siesta
#endif

  implicit none

  external :: io_assign, io_close

PUBLIC :: coordsFromSocket, forcesToSocket

PRIVATE  ! Nothing is declared public beyond this point

  character(len=*), parameter :: siesta_xunit = 'Bohr'
  character(len=*), parameter :: siesta_eunit = 'Ry'
  character(len=*), parameter :: siesta_funit = 'Ry/Bohr'
  character(len=*), parameter :: siesta_sunit = 'Ry/Bohr**3'

  INTEGER, SAVE :: socket, nat
  REAL*8, ALLOCATABLE,SAVE :: combuf(:), cellv

  integer,           save :: iuc, iuf

#ifdef MPI
  integer :: MPIerror
#endif

CONTAINS

subroutine coordsFromSocket( na, xa, cell )
  ! Reads coordinates from socket
  implicit none
  integer,  intent(in)  :: na         ! Number of atoms
  real(dp), intent(out) :: xa(3,na)   ! Atomic coordinates (bohr)
  real(dp), intent(out) :: cell(3,3)  ! Lattice vectors (bohr)

  logical, save     :: firstTime = .true.

  INTEGER, PARAMETER :: MSGLEN=12
  LOGICAL :: isinit=.false., hasdata=.false., firststep=.true., isunix
  CHARACTER*12 :: header
  CHARACTER*1024 :: parbuffer, host
  INTEGER ::  inet, port, i
  REAL *8 :: cellh(3,3), cellih(3,3), mtxbuf(9)

  ! Open the socket -- this will be shared from the receive and the send side of communication
  if (firstTime) then
    write(*,*) " Opening socket fot two-way communication. "
    call fdf_global_get(host, "Master.address", "localhost")
    call fdf_global_get(port, "Master.port", 31415)
    call fdf_global_get(isunix, "Master.unix", .false.)

    if (isunix) inet=0
    host = TRIM(host)//achar(0)

    if (IOnode) call open_socket(socket, inet, port, host) 

    firstTime = .false.
  end if ! (firstTime .and. IOnode)
  header=""

! Read coordinates from pipe
  if (IOnode) then
    write(*,*) " Receiving coordinates from socket. "
    do
      call readbuffer(socket, header, MSGLEN)
      if (trim(header)/='STATUS') exit
      call writebuffer(socket,"READY       ",MSGLEN)      
    enddo
  endif

#ifdef MPI !broadcasts the message we just received
  call MPI_Bcast(header, 12, MPI_Character, 0,  &
                 MPI_Comm_World, MPIerror)
#endif

  if (trim(header) == "POSDATA") then 
    !receives the positions & the cell data
    ! first reads cell and the number of atoms
    if (ionode) call readbuffer(socket, mtxbuf, 9)
    cellh = RESHAPE(mtxbuf,(/3,3/))
    if (ionode) call readbuffer(socket, mtxbuf, 9)
    cellih = RESHAPE(mtxbuf,(/3,3/))
    if (ionode) call readbuffer(socket, nat)
  
#ifdef MPI
    call MPI_Bcast(cellh(1,1),9, MPI_Double_Precision,0,  &
            MPI_Comm_World, MPIerror)
    call MPI_Bcast(cellih(1,1),9, MPI_Double_Precision,0,  &
            MPI_Comm_World, MPIerror)
    call MPI_Bcast(nat,1, MPI_Integer,0,  &
            MPI_Comm_World, MPIerror)
#endif
    if (.not.allocated(combuf)) then
      allocate(combuf(3*nat))
    end if
    if (ionode) call readbuffer(socket, combuf, nat*3)

#ifdef MPI
    call MPI_Bcast(combuf,3*nat, MPI_Double_Precision,0,  &
                 MPI_Comm_World, MPIerror)
#endif

   print '(/,3a,/,(3f12.6))', &
      'coordsFromSocket: cell (',trim('Bohr'),') =', cellh
   print '(  3a,/,(3f12.6))', &
      'coordsFromSocket: xa (',trim('Bohr'),') =', xa
   ! Convert coordinate units         
    cell = cellh * fdf_convfac( 'Bohr', siesta_xunit )
    xa   = RESHAPE(combuf, (/ 3 , nat /) )  &
           * fdf_convfac( 'Bohr', siesta_xunit )
   ! Computes the volume of the cell to be used to scale the stress tensor
    cellv = ABS(cell(1,1) * cell(2,2) * cell(3,3)) ! the cell is upper triangular, so super easy!
    write(*,*) "CELL VOLUME", cellv  
 else
   call die('coordsFromSocket: Unexpected message from server')
 end if
   
end subroutine coordsFromSocket


subroutine forcesToSocket( na, energy, forces, stress )
! Writes stress and forces to socket
  implicit none
  integer,  intent(in) :: na            ! Number of atoms
  real(dp), intent(in) :: energy        ! Total energy (Ry)
  real(dp), intent(in) :: forces(3,na)  ! Atomic forces (Ry/bohr)
  real(dp), intent(in) :: stress(3,3)   ! Stress tensor (Ry/bohr^3)

  integer           :: i, ia
  real*8   :: pot, vir(3,3)
  INTEGER, PARAMETER :: MSGLEN=12
  CHARACTER*12 :: header

  if (IOnode) write(*,*) " Sending forces to socket. "

! Convert physical units
  pot = energy * fdf_convfac( siesta_eunit, 'Ry' ) * 0.5 ! i-PI wants Hartree!
  ! we want to return the potential energy virial term, *without* the volume scaling
  ! also i-PI and Siesta use different sign conventions for the virial, so we need to cope with that.
  vir =  -1.0 * stress * cellv * fdf_convfac( siesta_eunit, 'Ry' ) * 0.5 
  combuf = RESHAPE(forces, (/ 3 * nat /) ) * fdf_convfac( siesta_funit, 'Ry/Bohr' ) * 0.5

! Write forces to socket
  if (IOnode) then
    do
      call readbuffer(socket, header, MSGLEN)
      if (trim(header)/='STATUS') exit
      call writebuffer(socket,"HAVEDATA    ",MSGLEN)
    enddo
    if (trim(header)/='GETFORCE') &
        call die('forcesToSocket: Error in socket communication!') 
    call writebuffer(socket,"FORCEREADY  ",MSGLEN)
    call writebuffer(socket,pot)
    call writebuffer(socket,nat)
    call writebuffer(socket,combuf,3*nat)
    call writebuffer(socket,reshape(vir,(/9/)),9)
          
    ! i-pi can also receive an arbitrary string, that will be printed out to the "extra" 
    ! trajectory file. this is useful if you want to return additional information, e.g.
    ! atomic charges, wannier centres, etc. one must return the number of characters, then
    ! the string. here we just send back zero characters.
    nat = 0
    call writebuffer(socket,nat)
  end if ! IOnode

end subroutine forcesToSocket

end module iosockets
