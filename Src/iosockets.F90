!!@LICENSE

!*****************************************************************************
! module iosockets
!*****************************************************************************
! Handles communications between siesta-as-a-subroutine and a driver 
! program, transfering coordinates and forces trough POSIX sockets.
! The routines that handle the other side of the communication are
! in file fsiesta_sockets.f90 and/or in the i-PI code.
! J.M.Soler, A.Garcia, and M.Ceriotti, 2014
!*****************************************************************************
!
!   Modules used:
! use precision,    only: dp
! use parallel,     only: IOnode
! use fdf
! use m_fdf_global, only: fdf_global_get
! use f90sockets,   only: open_socket, writebuffer, readbuffer
! use sys,          only: die, bye
! use m_mpi_utils,  only: broadcast
! use cellSubs,     only: volcel
! use mpi_siesta
!
!   Public derived types:
! none
!
!   Public variables and arrays:
! none
!
!   Public procedures:
! coordsFromSocket, &! receive coords and cell vectors from master program
! forcesToSocket     ! send energy, forces, and stress to master program
!
!*****************************************************************************
! subroutine coordsFromSocket( na, xa, cell )
!   integer,  intent(in) :: na            ! Number of atoms
!   real(dp), intent(out):: xa(3,na)      ! Atomic coordinates (bohr)
!   real(dp), intent(out):: cell(3,3)     ! Lattice vectors (bohr)
!*****************************************************************************
! subroutine forcesToSocket( na, energy, forces, stress )
!   integer,  intent(in) :: na            ! Number of atoms
!   real(dp), intent(in) :: energy        ! Total energy (Ry)
!   real(dp), intent(in) :: forces(3,na)  ! Atomic forces (Ry/bohr)
!   real(dp), intent(in) :: stress(3,3)   ! Stress tensor (Ry/bohr^3)
!*****************************************************************************

module iosockets

! Used modules
  use precision,    only: dp
  use parallel,     only: IOnode
  use fdf
  use m_fdf_global, only: fdf_global_get
  use f90sockets,   only: open_socket, writebuffer, readbuffer
  use sys,          only: die, bye
  use m_mpi_utils,  only: broadcast
  use cellSubs,     only: volcel
#ifdef MPI
  use mpi_siesta
#endif

  implicit none

! Public procedures
PUBLIC :: &
  coordsFromSocket, &! receive coords and cell vectors from master program
  forcesToSocket     ! send energy, forces, and stress to master program

PRIVATE  ! Nothing is declared public beyond this point

! Private module parameters
  integer,          parameter :: unit_len = 12
  character(len=*), parameter :: siesta_xunit = 'bohr'
  character(len=*), parameter :: siesta_eunit = 'ry'
  character(len=*), parameter ::    ipi_xunit = 'bohr'
  character(len=*), parameter ::    ipi_eunit = 'hartree'

! Private module variables
  integer, save :: socket
  real(dp),save :: cellv
  character(len=unit_len),save:: master='unknown'
  character(len=32),      save:: master_eunit='unknown', &
                                 master_xunit='unknown'
#ifdef MPI
  integer :: MPIerror
#endif

CONTAINS

!*****************************************************************************

subroutine coordsFromSocket( na, xa, cell )
  ! Reads coordinates from socket
  implicit none
  integer,  intent(in)  :: na         ! Number of atoms
  real(dp), intent(out) :: xa(3,na)   ! Atomic coordinates (bohr)
  real(dp), intent(out) :: cell(3,3)  ! Lattice vectors (bohr)

! Local parameters
  character(len=*),parameter:: myName='coordsFromSocket '
  INTEGER, PARAMETER :: MSGLEN=12

! Local variables and arrays
  logical, save  :: firstTime = .true.
  LOGICAL        :: isunix
  CHARACTER*12   :: code, header
  CHARACTER*1024 :: parbuffer, host
  INTEGER        :: inet, port, n
  REAL*8         :: aux(9), c(9)
  REAL*8,ALLOCATABLE:: x(:)

! Open the socket, shared by the receive and send sides of communication
  if (firstTime) then
    call fdf_global_get(master, "Master.code", "fsiesta")
    call fdf_global_get(host,   "Master.address", "localhost")
    call fdf_global_get(port,   "Master.port", 31415)
    call fdf_global_get(isunix, "Master.unix", .false.)

    inet = 1
    if (isunix) inet=0
    host = TRIM(host)//achar(0)

    if (IOnode) then
      print*,myName//"Opening socket fot two-way communication."
      call open_socket(socket, inet, port, host)
    endif

    firstTime = .false.
  end if ! (firstTime .and. IOnode)

! Read header from socket
  header=""
  if (IOnode) then
    print'(/,a)',myName//"Receiving coordinates from socket. "
    do
      call readbuffer(socket, header, MSGLEN)
      if (trim(master)=='ipi') then
        if (trim(header)/='STATUS') exit 
        call writebuffer(socket,"READY       ",MSGLEN)
      elseif (trim(master)=='fsiesta') then
        if (trim(header)/='wait') exit
        call writebuffer(socket,"ready       ",MSGLEN)
      endif ! trim(master)
    enddo
  endif ! IOnode
#ifdef MPI
  call MPI_Bcast(header, 12, MPI_Character, 0, MPI_Comm_World, MPIerror)
#endif

! Read cell vectors from socket, as a single buffer vector
  if (IOnode) then
    if (master=="i-pi" .and. trim(header)=="POSDATA") then 
      call readbuffer(socket, c, 9)
      call readbuffer(socket, aux, 9)    ! not used in siesta
      master_xunit = ipi_xunit
      master_eunit = ipi_eunit
    elseif (master=="fsiesta" .and. trim(header)=="begin_coords") then
      call readbuffer(socket, master_xunit, unit_len)
      call readbuffer(socket, master_eunit, unit_len)
      call readbuffer(socket, c, 9)
    else
      call die(myName//'Unexpected message from master')
    end if
  endif

! Broadcast and copy cell from buffer
#ifdef MPI
  call MPI_Bcast(master_xunit, unit_len, MPI_Character,0,  &
                 MPI_Comm_World, MPIerror)
  call MPI_Bcast(master_eunit, unit_len, MPI_Character,0,  &
                 MPI_Comm_World, MPIerror)
  call MPI_Bcast(c,9, MPI_Double_Precision,0, MPI_Comm_World, MPIerror)
#endif
  cell = RESHAPE( c, (/3,3/) )

! Read and check number of atoms
  if (ionode) then
    call readbuffer(socket, n)
    if (n/=na) call die(myName//'ERROR: unexpected number of atoms')
  endif

! Read atomic coordinates
  allocate(x(3*na))
  if (ionode) call readbuffer(socket, x, 3*na)
#ifdef MPI
  call MPI_Bcast(x,3*na, MPI_Double_Precision,0, MPI_Comm_World, MPIerror)
#endif
  xa = RESHAPE( x, (/3,na/) )
  deallocate(x)

! Print coordinates and cell vectors received
  print '(/,4a,/,(3f12.6))', myName,'cell (',trim(master_xunit),') =', cell
  print '(  4a,/,(3f12.6))', myName,'coords (',trim(master_xunit),') =', xa

! Convert physical units
  cell = cell * fdf_convfac( master_xunit, siesta_xunit )
  xa   = xa   * fdf_convfac( master_xunit, siesta_xunit )

! Compute and save the cell volume, to be used by forcesToSocket
  cellv = volcel(cell)
   
end subroutine coordsFromSocket

!*****************************************************************************

subroutine forcesToSocket( na, energy, forces, stress )
! Writes stress and forces to socket
  implicit none
  integer,  intent(in) :: na            ! Number of atoms
  real(dp), intent(in) :: energy        ! Total energy (Ry)
  real(dp), intent(in) :: forces(3,na)  ! Atomic forces (Ry/bohr)
  real(dp), intent(in) :: stress(3,3)   ! Stress tensor (Ry/bohr^3)

! Local parameters
  character(len=*),parameter:: myName='forcesToSocket '
  INTEGER, PARAMETER :: MSGLEN=12

! Local variables and arrays
  CHARACTER*12      :: header
  real*8            :: e, s(9), vir(9)
  real*8,allocatable:: f(:)

! Copy input to local variables
  allocate(f(3*na))
  e = energy
  f = reshape( forces, (/3*na/) )
  s = reshape( stress, (/9/) )
  
! Convert physical units
  e = e * fdf_convfac( siesta_eunit, master_eunit )
  f = f * fdf_convfac( siesta_eunit, master_eunit ) &
        / fdf_convfac( siesta_xunit, master_xunit )
  s = s * fdf_convfac( siesta_eunit, master_eunit ) &
        / fdf_convfac( siesta_xunit, master_xunit )**3

! Find virial tensor for i-pi, in master's units
  vir = -s * cellv * fdf_convfac( siesta_xunit, master_xunit )**3

! Write forces to socket
  if (IOnode) then
    print'(/,a)',myName//"Sending forces to socket. "
    if (trim(master)=='i-pi') then
      do
        call readbuffer(socket, header, MSGLEN)
        if (trim(header)/='STATUS') exit
        call writebuffer(socket,"HAVEDATA    ",MSGLEN)
      enddo
      if (trim(header)/='GETFORCE') &
          call die(myName//'ERROR in socket communication!') 
      call writebuffer(socket,"FORCEREADY  ",MSGLEN)
      call writebuffer(socket,e)
      call writebuffer(socket,na)
      call writebuffer(socket,f,3*na)
      call writebuffer(socket,vir,9)       
      ! i-pi can also receive an arbitrary string, that will be printed 
      ! out to the "extra" trajectory file. This is useful if you want 
      ! to return additional information, e.g. atomic charges, wannier 
      ! centres, etc. one must return the number of characters, then
      ! the string. here we just send back zero characters.
      call writebuffer(socket,0)
!      call writebuffer(socket,'')  ! would this not be clearer?
    elseif (trim(master)=='fsiesta') then
      do
        call readbuffer(socket, header, MSGLEN)
        if (trim(header)/='wait') exit
      enddo
      call writebuffer(socket,'begin_forces',12)
      call writebuffer(socket,e)
      call writebuffer(socket,s,9)
      call writebuffer(socket,na)
      call writebuffer(socket,f,3*na)
      call writebuffer(socket,'end_forces',10)
    else
      call die(myName//'ERROR: unknown master')
    endif ! trim(master)
  end if ! IOnode

! Print energy, forces, and stress tensor sent to master
  print '(/,a,f12.6)',      myName// &
    'energy ('//trim(master_eunit)//') =', e
  print '(  a,/,(3f12.6))', myName// &
    'forces ('//trim(master_eunit)//'/'//trim(master_xunit)//') =', f
  print '(  a,/,(3f12.6))', myName// &
    'stress ('//trim(master_eunit)//'/'//trim(master_xunit)//'^3) =', s
  deallocate(f)

end subroutine forcesToSocket

end module iosockets
