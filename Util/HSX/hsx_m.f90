module hsx_m

!
! Experimental module to process HSX files
!
implicit none

integer, parameter, private :: dp = selected_real_kind(14,100)
integer, parameter, private :: sp = selected_real_kind(6,30)

public :: read_hsx_file

type, public :: hsx_t
  integer :: nspecies, na_u, no_u, no_s, nspin, nh
  logical :: gamma
  integer, pointer :: no(:)=>null(), nquant(:,:)=>null(), lquant(:,:)=>null(), zeta(:,:)=>null()
  integer, pointer :: iaorb(:)=>null(), iphorb(:)=>null()
  character(len=20), pointer :: label(:)=>null()

  integer, pointer :: numh(:)=>null() , listhptr(:)=>null(), listh(:)=>null()  
  integer, pointer :: indxuo(:)=>null()
  real(dp), pointer :: hamilt(:,:)=>null()
  real(dp), pointer :: Sover(:)=>null(), xij(:,:)=>null()
  integer, pointer :: isa(:)=>null()
  real(dp), pointer :: zval(:)=>null()

end type

private

CONTAINS

subroutine read_hsx_file(hsx,fname)
type(hsx_t), intent(out)  :: hsx
character(len=*), intent(in) :: fname

!
! Reads HSX file "fname" and stores the info in the hsx data structure
!
  integer, allocatable  :: ibuff(:)
  real(sp), allocatable  :: hbuff(:)
  real(sp), allocatable  :: buff3(:,:)

  integer numx, ind, no_u, nnz, na_u, nspecies, nspin, nh, i
  integer :: im, is, hs_u, ia, io, iostat, k, naoatx, no_s
  real(dp) :: aux1, aux2
  logical  :: debug = .false.

  call get_unit_number(hs_u)
  print *, "Using unit: ", hs_u

  open(hs_u,file=trim(fname),status='old',form='unformatted')

  read(hs_u,iostat=iostat) hsx%no_u, hsx%no_s, hsx%nspin, hsx%nh
  if (iostat /= 0) STOP "nnao, no_s..."

  no_u = hsx%no_u
  no_s = hsx%no_s

  read(hs_u,iostat=iostat) hsx%gamma
  if (iostat /= 0) STOP "gamma"
  IF (DEBUG) PRINT *, "GAMMA=", hsx%gamma
  if (.not. hsx%gamma) then
     allocate(hsx%indxuo(no_s))
     read(hs_u) (hsx%indxuo(i),i=1,hsx%no_s)
  else
     allocate(hsx%indxuo(hsx%no_u))
     do i=1,hsx%no_u
        hsx%indxuo(i) = i
     enddo
  endif

  nh  = hsx%nh
  nspin = hsx%nspin
  print *, "nh: ", nh
  allocate (hsx%numh(no_u), hsx%listhptr(no_u), hsx%listh(nh))

       allocate (hsx%xij(3,nh),stat=iostat)
       allocate (hsx%hamilt(nh,nspin),stat=iostat)
       allocate (hsx%Sover(nh),stat=iostat)

  read(hs_u,iostat=iostat) (hsx%numh(io), io=1,no_u)      
  if (iostat /= 0) STOP "numh"

  numx = maxval(hsx%numh(1:no_u))
  allocate(ibuff(numx), hbuff(numx), buff3(3,numx))

  nnz = sum(hsx%numh(1:hsx%no_u))
   if (nnz > nh) STOP "nh overflow in HS"
  ! Create listhptr 
  hsx%listhptr(1)=0
  do io=2,hsx%no_u
     hsx%listhptr(io)=hsx%listhptr(io-1)+hsx%numh(io-1)
  enddo

  do io=1,hsx%no_u
     read(hs_u,iostat=iostat) (ibuff(im), im=1,hsx%numh(io))
     if (iostat /= 0) STOP "listh"
     do im=1,hsx%numh(io)
        hsx%listh(hsx%listhptr(io)+im) = ibuff(im)
     enddo
  enddo

  do is=1,hsx%nspin
     do io=1,hsx%no_u
        read(hs_u,iostat=iostat) (hbuff(im), im=1,hsx%numh(io))
        if (iostat /= 0) STOP "Hamilt"
        do im=1,hsx%numh(io)
           hsx%hamilt(hsx%listhptr(io)+im,is) = hbuff(im)
           if (debug) print *, "Hamilt ", io, im, hbuff(im)
        enddo
     enddo
  enddo
  !
  !       Read overlap matrix
  !
  do io=1,hsx%no_u
     read(hs_u,iostat=iostat) (hbuff(im), im=1,hsx%numh(io))
     if (iostat /= 0) STOP "Overlap matrix read error"
     do im=1,hsx%numh(io)
        hsx%Sover(hsx%listhptr(io)+im) = hbuff(im)
        if (debug) print *, "S ", io, im, hbuff(im)
     enddo
  enddo

  read(hs_u,iostat=iostat) aux1, aux2                   ! qtot, temp
  if (debug) print *, "QTOT, Temp: ", aux1, aux2
  if (iostat /= 0) then
     if (debug) print *, "iostat:", iostat
     STOP "qtot, temp"
  endif

  !
  !        Always read xijk
  !
  do io=1,hsx%no_u
     read(hs_u,iostat=iostat) ((buff3(k,im), k=1,3), im=1,hsx%numh(io))
     if (iostat /= 0) STOP "xij(k)"
     do im=1,hsx%numh(io)
        ind = hsx%listhptr(io)+im
        if (debug) print *, "xijk ", buff3(:,im)
        hsx%xij(1:3,ind) = buff3(1:3,im)
     enddo
  enddo
  !
  !        Read auxiliary info
  !
  read(hs_u) hsx%nspecies
  nspecies = hsx%nspecies
  print *, "nspecies: ", nspecies
  allocate(hsx%label(nspecies), hsx%zval(nspecies), hsx%no(nspecies))
  read(hs_u) (hsx%label(is),hsx%zval(is),hsx%no(is), is=1,nspecies)
  naoatx = maxval(hsx%no(1:nspecies))
  allocate (hsx%nquant(nspecies,naoatx), hsx%lquant(nspecies,naoatx), &
       hsx%zeta(nspecies,naoatx))
  do is=1, nspecies
     do io=1, hsx%no(is)
        read(hs_u) hsx%nquant(is,io), hsx%lquant(is,io), hsx%zeta(is,io)
     enddo
  enddo
  read(hs_u) hsx%na_u
  na_u = hsx%na_u
  allocate(hsx%isa(na_u))
  allocate(hsx%iaorb(no_u), hsx%iphorb(no_u))
  read(hs_u) (hsx%isa(ia), ia=1,na_u)
  read(hs_u) (hsx%iaorb(io), hsx%iphorb(io), io=1,no_u)

  close(hs_u)
  deallocate(ibuff, hbuff, buff3)

end subroutine read_hsx_file

!----------------------------------------------------------
subroutine get_unit_number(lun)
integer, intent(out) :: lun

logical :: used
integer :: iostat

do lun= 10, 99
   inquire(unit=lun, opened=used, iostat=iostat)
   if (iostat .ne. 0) used = .true.
   if (.not. used) return
enddo
STOP "Cannot get unit"
end subroutine get_unit_number

end module hsx_m
