module io_hs

implicit none
public :: read_hs_file

CONTAINS
subroutine read_hs_file()
  use main_vars
  use precision, only: sp

  integer, allocatable  :: ibuff(:)
  real(sp), allocatable  :: hbuff(:)
  real(sp), allocatable  :: buff3(:,:)

  integer numx

  write(6,"(1x,a,'.HS ...')",advance='no') trim(sflnm)
  open(hs_u,file=trim(sflnm)//'.HSX',status='old',form='unformatted')

  read(hs_u,iostat=iostat) nnao, no_s, nspin, nh
  print *, "nnao, no_s, nspin, nh:",  nnao, no_s, nspin, nh
  if (iostat /= 0) STOP "nnao, no_s..."
  if (nnao /= nao) STOP "norbs inconsistency"
  no_u = nao

  read(hs_u,iostat=iostat) gamma
  if (iostat /= 0) STOP "gamma"
  IF (DEBUG) PRINT *, "GAMMA=", gamma
  if (.not. gamma) then
     allocate(indxuo(no_s))
     read(hs_u) (indxuo(i),i=1,no_s)
  else
     allocate(indxuo(no_u))
     do i=1,no_u
        indxuo(i) = i
     enddo
  endif

  if (debug) print *, "HS read: nh, nsp, nnao: ", nh, nspin, nnao
  if (nnao.ne.nao) STOP " nnao .ne. nao in HS"
  if (wfs_x.and.(nspin.ne.nsp)) STOP " nspin .ne. nsp in HS"
  nsp=nspin
  allocate (numh(nao), listhptr(nao), listh(nh))

!!$        allocate (hamilt(nh,nspin))
!!$        allocate (Sover(nh))
!!$        allocate (xij(3,nh))

  read(hs_u,iostat=iostat) (numh(io), io=1,no_u)         ! numhg
  if (iostat /= 0) STOP "numh(io)"
  do io=1,no_u
     if (debug) print *, "numhg ", io, numh(io)
  enddo

  numx = maxval(numh(:))
  allocate(ibuff(numx), hbuff(numx), buff3(3,numx))


  ! Create listhptr 
  listhptr(1)=0
  do io=2,no_u
     listhptr(io)=listhptr(io-1)+numh(io-1)
  enddo
  if (listhptr(no_u)+numh(no_u).gt.nh) STOP "nh overflow in HS"

  allocate (nsr(no_u,no_u))
  nsr(:,:)=0  
  do io=1,no_u
     read(hs_u,iostat=iostat) (ibuff(im), im=1,numh(io))
     if (iostat /= 0) STOP "listh"
     do im=1,numh(io)
        listh(listhptr(io)+im) = ibuff(im)
        ii=listh(listhptr(io)+im)
        io2=ii-(ii-1)/no_u*no_u
        if (io2 /= indxuo(ii)) then
           print *, "mismatch: ii, io2, indxuo(ii):", ii, io2, indxuo(ii)
           stop
        endif
        nsr(io,io2)=nsr(io,io2)+1
        if (debug) print *, "listh ", io, im, listh(listhptr(io)+im)
     enddo
  enddo

  nsrmx = maxval(nsr)
  print *, "Maximum supercell multiplicity: ", nsrmx

  allocate (sr(no_u,no_u,nsrmx))
  allocate (hr(no_u,no_u,nsrmx,nspin))

  do is=1,nspin
     nsr(:,:)=0   ! Wasteful to repeat this (and below), but the
     ! structure of the HS file is crazy
     do io=1,no_u
        read(hs_u,iostat=iostat) (hbuff(im), im=1,numh(io))
        if (iostat /= 0) STOP "Hamilt"
        do im=1,numh(io)
           ii=listh(listhptr(io)+im)
           io2=ii-(ii-1)/no_u*no_u
           nsr(io,io2)=nsr(io,io2)+1
           if (nsr(io,io2).gt.nsrmx)  &
                stop "* ERROR Cells in supercell limit exceeded."
           hr(io,io2,nsr(io,io2),is) = hbuff(im)
!!!             hamilt(listhptr(io)+im,is) = hr(io,io2,nsr(io,io2),is)
           if (debug) print *, "Hamilt ", io, im, hr(io,io2,nsr(io,io2),is)
        enddo
     enddo
  enddo
  !
  !       Read overlap matrix
  !
  nsr(:,:)=0
  do io=1,no_u
     read(hs_u,iostat=iostat) (hbuff(im), im=1,numh(io))
     if (iostat /= 0) STOP "Overlap matrix read error"
     do im=1,numh(io)
        ii=listh(listhptr(io)+im)
        io2=ii-(ii-1)/no_u*no_u
        nsr(io,io2)=nsr(io,io2)+1
        sr(io,io2,nsr(io,io2)) = hbuff(im)
!!!            Sover(listhptr(io)+im) = sr(io,io2,nsr(io,io2)) 
        if (debug) print *, "S ", io, im, sr(io,io2,nsr(io,io2))
     enddo
  enddo

  read(hs_u,iostat=iostat) aux1, aux2                   ! qtot, temp
  if (debug) print *, "QTOT, Temp: ", aux1, aux2
  if (iostat /= 0) then
     if (debug) print *, "iostat:", iostat
     STOP "qtot, temp"
  endif

  mxsr=maxval(nsr(:no_u,:no_u))
  allocate (rn(3,no_u,no_u,mxsr), dt(no_u,no_u,mxsr))
  !
  !        Always read xijk
  !
  nsr(:,:)=0
  do io=1,no_u
     read(hs_u,iostat=iostat) ((buff3(k,im), k=1,3), im=1,numh(io))
     if (iostat /= 0) STOP "xij(k)"
     do im=1,numh(io)
        ii=listh(listhptr(io)+im)
        io2=ii-(ii-1)/no_u*no_u
        nsr(io,io2)=nsr(io,io2)+1
        rn(:,io,io2,nsr(io,io2))= buff3(:,im)

        if (debug) print *, "xijk ", buff3(:,im)
!!!            xij(1:3,listhptr(io)+im) = (/ ax, ay, az /)
!!!        dt(io,io2,nsr(io,io2))= sqrt(ax*ax+ay*ay+az*az) / Ang
        dt(io,io2,nsr(io,io2))=                 &
               sqrt(dot_product(buff3(:,im),buff3(:,im))) / Ang
     enddo
  enddo
  !
  !        Read auxiliary info
  !
  read(hs_u) nspecies
  allocate(label(nspecies), zval(nspecies), no(nspecies))
  read(hs_u) (label(is),zval(is),no(is), is=1,nspecies)
  naoatx = maxval(no(1:nspecies))
  allocate (nquant(nspecies,naoatx), lquant(nspecies,naoatx), &
       zeta(nspecies,naoatx))
  do is=1, nspecies
     do io=1, no(is)
        read(hs_u) nquant(is,io), lquant(is,io), zeta(is,io)
     enddo
  enddo
  read(hs_u) na_u
  allocate(isa(na_u))
  allocate(iaorb(no_u), iphorb(no_u))
  read(hs_u) (isa(ia), ia=1,na_u)
  read(hs_u) (iaorb(io), iphorb(io), io=1,no_u)

  close(hs_u)
  deallocate(ibuff, hbuff, buff3)

end subroutine read_hs_file
end module io_hs
