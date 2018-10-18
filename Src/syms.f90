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
! written by Matthieu Verstraete around 12 June 2014
! 

module m_syms

  !use spglib_f08 ! this is also included in spglib and 
  !  could be a nice way to avoid recoding the interfaces for everything.
  !  however, it does imply making spglib with exactly the same fortran compiler as siesta
  implicit none

  type, public :: syms_type
    integer :: nsymop
    integer :: space_group
    character (len=11) :: symbol
    integer, allocatable :: symops(:,:,:)
    !integer, allocatable :: symops_cart(:,:,:)
    double precision, allocatable :: symops_cart(:,:,:)
    integer, allocatable :: symops_recip(:,:,:)
    double precision, allocatable :: trans(:,:)
    double precision, allocatable :: trans_cart(:,:)
  end type syms_type

  ! for the moment this single instance of the syms is exported by the module.
  ! if you want to play with several, declare 1 per subroutine,
  ! and use it as the intent(out) argument to the init_syms subroutine
  type(syms_type), public :: syms_t

  public :: init_syms, delete_syms

!these functions are defined in spglib_f.c: interfaces copied from Xavier
!Andrade's module in octopus and updated to spglib 1.5.2 on 10 June 2014
  interface
    subroutine spg_get_multiplicity(max_size, cell, xa, isa, na, symprec)
      integer, intent(out) :: max_size
      real(8), intent(in) :: cell
      real(8), intent(in) :: xa
      integer, intent(in) :: isa
      integer, intent(in) :: na
      real(8), intent(in) :: symprec
    end subroutine spg_get_multiplicity

    subroutine spg_get_symmetry(nsymop, symop, trans, max_size, cell, xa, isa, na, symprec)
      integer, intent(out) :: nsymop
      integer, intent(out) :: symop
      real(8), intent(out) :: trans
      integer, intent(in)  :: max_size
      real(8), intent(in)  :: cell
      real(8), intent(in)  :: xa
      integer, intent(in)  :: isa
      integer, intent(in)  :: na
      real(8), intent(in)  :: symprec
    end subroutine spg_get_symmetry

! this is no longer in spglib 1.5.2
!    subroutine show_symmetry(cell, xa, isa, na, symprec)
!      real(8), intent(in) :: cell
!      real(8), intent(in) :: xa
!      integer, intent(in) :: isa
!      integer, intent(in) :: na
!      real(8), intent(in) :: symprec
!    end subroutine show_symmetry

    subroutine spg_get_international(spgrp, symbol, cell, xa, isa, na, symprec)
      integer, intent(out) :: spgrp
      character*11, intent(out) :: symbol
      real(8), intent(in) :: cell
      real(8), intent(in) :: xa
      integer, intent(in) :: isa
      integer, intent(in) :: na
      real(8), intent(in) :: symprec
    end subroutine spg_get_international

  end interface

contains

  subroutine init_syms( syms_this, cell, na, isa, xa )
! *****************************************************************
! Arguments:
! real*8  cell(3,3)    : input lattice vectors (Bohr)
! integer na           : input number of atoms
! integer isa(na)      : input species indexes
! real*8  xa(3,na)     : input atomic cartesian coordinates (Bohr)
! *****************************************************************
    implicit         none
    type(syms_type), intent(out) :: syms_this
    integer, intent(in)          :: na, isa(na)
    double precision, intent(inout) :: cell(3,3), xa(3,na)

! local vars
    integer :: max_size
    integer :: ii, info
    double precision :: symprec
    double precision :: cellinv(3,3)
    double precision :: celltransp(3,3)
    double precision :: dsymop(3,3)
    double precision :: dsymopinv(3,3)
    double precision :: zerovec(3) = (/0.0d0, 0.0d0,0.0d0/)

    double precision, allocatable :: xreda(:,:)
    integer, allocatable :: symops_tmp(:,:,:)
    double precision, allocatable :: trans_tmp(:,:)

    symprec= 1.e-6 ! this is hard coded for the moment - could be input var for init_syms

    allocate (xreda(3,na))
    call cart2red(cell, na, xa, xreda)
! not necessary - no effect on spgrp recognition by spglib
!    do ii = 1, na
!      xreda(:,ii) = wrapvec_zero_one(xreda(:,ii))
!    end do
!DEBUG
print *,  'xred = '
print '(3E20.10)', xreda
!END DEBUG

    ! refines positions with highest symmetry found within symprec
    ! REMOVED 13/6/2014!!! - this returns the conventional cell!!!
    !call spg_refine_cell(cell, xreda, isa, na, symprec)
    !print *, 'xred refined'
    !print '(3E20.10)', xreda
    !print *, 'cell refined'
    !print '(3E20.10)', cell

    celltransp = transpose(cell)
print *, 'transposed cell array: 1 vector on each line?'
print '(3E20.10)', celltransp(:,1)
print '(3E20.10)', celltransp(:,2)
print '(3E20.10)', celltransp(:,3)

    call spg_get_international(syms_this%space_group, syms_this%symbol, celltransp(1,1), zerovec(1), isa(1), 1, symprec)
print *, 'for empty lattice space group is ', syms_this%space_group, syms_this%symbol

    call spg_get_international(syms_this%space_group, syms_this%symbol, celltransp(1,1), xreda(1,1), isa(1), na, symprec)
print *, 'for full lattice space group is ', syms_this%space_group, syms_this%symbol

    call spg_get_multiplicity(max_size, celltransp(1,1), xreda(1,1), isa(1), na, symprec)

    ! spglib returns row-major not column-major matrices!!! --DAS
    ! should we transpose these after reading in???
    allocate(symops_tmp(1:3, 1:3, 1:max_size))
    allocate(trans_tmp(1:3, 1:max_size))

    call spg_get_symmetry(syms_this%nsymop, symops_tmp(1, 1, 1), trans_tmp(1, 1), &
      max_size, celltransp(1, 1), xreda(1, 1), isa(1), na, symprec)

    allocate(syms_this%symops(1:3, 1:3, 1:syms_this%nsymop))
    do ii = 1, max_size
! transpose needed or not? only visible on triclinic or hexagonal etc..
! Checked for 1 case with P6_3 c m (#185) comparing to abinit - a transpose is needed!!!
      syms_this%symops(:,:,ii) = transpose(symops_tmp(:,:,ii))
      !syms_this%symops(:,:,ii) = symops_tmp(:,:,ii)
    end do
    deallocate(symops_tmp)

! get symops in cartesian coordinates - should still be integer -1 0 1 elements
    allocate(syms_this%symops_cart(1:3, 1:3, 1:syms_this%nsymop))
    call INVER(cell,cellinv,3,3,info)
    if (info .ne. 0) stop 'subroutine init_syms: error in inverse matrix of cell'

    do ii = 1, syms_this%nsymop
      syms_this%symops_cart(:, :, ii) = nint( matmul(cell, &
&       matmul(dble(syms_this%symops(:, :, ii)), cellinv)) )

! it appears these matrices can be non integer for hexagonal systems... 
!  for the moment, make them dble
!      if (  any(abs(dble(syms_this%symops_cart(:, :, ii)) &
!&       - matmul(cell, matmul(dble(syms_this%symops(:, :, ii)), cellinv))) > 1.e-10)  ) then
!         print *, ii, syms_this%symops(:, :, ii)
!         print *, matmul(cell, matmul(dble(syms_this%symops(:, :, ii)), cellinv))
!         stop 'error : cartesian symop element is not -1 0 +1'
!      end if

    end do

! get symops in reciprocal space = normal^-1 ^T
    allocate(syms_this%symops_recip(1:3, 1:3, 1:syms_this%nsymop))
    do ii = 1, syms_this%nsymop
      dsymop = dble(syms_this%symops(:, :, ii)) 
      call INVER(dsymop,dsymopinv,3,3,info)
      if (info .ne. 0) stop 'subroutine init_syms: error in inverse matrix of symop'
      dsymopinv = transpose (dsymopinv)
      syms_this%symops_recip(:,:,ii) = int(dsymopinv)
! DEBUG - comment out later
  if (any(abs(dble(syms_this%symops_recip(:, :, ii)) - dsymopinv(:,:)) > 1.e-10)) then
    stop 'error : recip symop element is not -1 0 +1'
  end if
! END DEBUG
    end do

    allocate(syms_this%trans(1:3, 1:syms_this%nsymop))
    syms_this%trans = trans_tmp(:,1:syms_this%nsymop)
    deallocate(trans_tmp)
    
    allocate(syms_this%trans_cart(1:3, 1:syms_this%nsymop))
    syms_this%trans_cart = matmul(cell, syms_this%trans)

  end subroutine init_syms

  subroutine delete_syms( syms_this )
! *****************************************************************
! Arguments:
! syms_type syms_t     : object with symops etc...
! *****************************************************************
    implicit         none
    type(syms_type), intent(inout) :: syms_this

    if(allocated(syms_this%symops)) deallocate (syms_this%symops)
    if(allocated(syms_this%symops_cart)) deallocate (syms_this%symops_cart)
    if(allocated(syms_this%trans)) deallocate (syms_this%trans)
    if(allocated(syms_this%trans_cart)) deallocate (syms_this%trans_cart)

  end subroutine delete_syms


  subroutine print_syms( syms_this )
! *****************************************************************
! Arguments:
! syms_type syms_t     : object with symops etc...
! *****************************************************************
    implicit         none
    type(syms_type), intent(in) :: syms_this


! local
    integer :: ii

    write (*,'(a)') " output of crystal symmetries "
    write (*,'(a,I6)') " nsymop = ", syms_this%nsymop
    write (*,'(2a)') " symbol = ", syms_this%symbol
    write (*,'(a)') " symops in reduced coordinates = "
    do ii = 1, syms_this%nsymop
      write (*,'(3(3I10,2x))') syms_this%symops(:,:,ii)
    end do
    write (*,'(a)') " trans in reduced coordinates = "
    do ii = 1, syms_this%nsymop
      write (*, '(3(E20.10,2x))') syms_this%trans(:,ii)
    end do
  end subroutine print_syms


  subroutine red2cart(cell, na, xred, xa)
! *****************************************************************
! Arguments:
! cell : unit cell vectors
! na : number of atoms
! xa : cartesian position vectors
! xred : reduced position vectors
! *****************************************************************
    implicit         none
    double precision,intent(in) :: cell(3,3)
    integer,intent(in) :: na
    double precision,intent(in) :: xred(3,na)
    double precision,intent(out) :: xa(3,na)

    xa = matmul(cell, xred)

  end subroutine red2cart


  subroutine cart2red(cell, na, xa, xred)
! *****************************************************************
! Arguments:
! cell : unit cell vectors
! na : number of atoms
! xa : cartesian position vectors
! xred : reduced position vectors
! *****************************************************************
    implicit         none
    double precision,intent(in) :: cell(3,3)
    integer,intent(in) :: na
    double precision,intent(in) :: xa(3,na)
    double precision,intent(out) :: xred(3,na)

    double precision :: cellinv(3,3)
    integer :: info

    call INVER(cell,cellinv,3,3,info)
    if (info .ne. 0) stop 'subroutine cart2red: error in inverse matrix of cell'

    xred = matmul(cellinv, xa)

  end subroutine cart2red

  function wrapvec_zero_one(num) result (red)
    implicit none
    double precision, intent(in) :: num(3)
    double precision :: red(3)

    integer :: ii

    do ii=1, 3
      red(ii) = wrap2zero_one(num(ii))
    end do

  end function wrapvec_zero_one

  function wrap2zero_one(num) result (red)
    implicit none
    double precision, intent(in) :: num
    double precision :: red

    if (num>0.0d0) then
      red=mod((num+1.d-13),1.0d0)-1.d-13
    else
      red=-mod(-(num-1.0d0+1.d-13),1.0d0)+1.0d0-1.d-13
    end if
    if(abs(red)<1.d-13)red=0.0d0

  end function wrap2zero_one
   
  subroutine mk_irred_k_grid(nkgrid, shiftk, nkpt, kpt_list)
    implicit none
    integer, intent(in) :: nkgrid(3)
    double precision, intent(in) :: shiftk(3)
    integer, intent(out) :: nkpt
    double precision, allocatable, intent(out) :: kpt_list(:,:)


    

  end subroutine mk_irred_k_grid

end module m_syms
