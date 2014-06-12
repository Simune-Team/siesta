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
! 

module m_syms

  implicit none

  type, public :: syms_type
    integer :: nsymop
    integer :: space_group
    character (len=11) :: symbol
    integer, allocatable :: symops(:,:,:)
    double precision, allocatable :: trans(:,:)
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
    double precision, intent(in) :: cell(3,3), xa(3,na)

! local vars
    integer :: max_size
    integer :: ii, info
    double precision :: symprec
    double precision :: cellinv(3,3)

    double precision, allocatable :: xreda(:,:)
    integer, allocatable :: symops_tmp(:,:,:)
    double precision, allocatable :: trans_tmp(:,:)

    symprec= 1.e-6 ! this is hard coded for the moment - could be input var for init_syms

    call INVER(cell,cellinv,3,3,info)
    if (info .ne. 0) stop 'subroutine constr: error in inverse matrix of cell'

    allocate (xreda(3,na))
    ! might have to check if the cellinv needs to be transposed here...
    xreda = matmul(cellinv, xa)
    call spg_get_international(syms_this%space_group, syms_this%symbol, cell(1,1), xreda(1,1), isa(1), na, symprec)

    call spg_get_multiplicity(max_size, cell(1,1), xreda(1,1), isa(1), na, symprec)

    ! spglib returns row-major not column-major matrices!!! --DAS
    ! should we transpose these after reading in???
    allocate(symops_tmp(1:3, 1:3, 1:max_size))
    allocate(trans_tmp(1:3, 1:max_size))

    call spg_get_symmetry(syms_this%nsymop, symops_tmp(1, 1, 1), trans_tmp(1, 1), &
      max_size, cell(1, 1), xreda(1, 1), isa(1), na, symprec)

    allocate(syms_this%symops(1:3, 1:3, 1:syms_this%nsymop))
    do ii = 1, max_size
      syms_this%symops(:,:,ii) = transpose(symops_tmp(:,:,ii)) ! transpose needed or not? only visible on triclinic or hexagonal etc..
    end do
    deallocate(symops_tmp)

    allocate(syms_this%trans(1:3, 1:syms_this%nsymop))
    syms_this%trans = trans_tmp(:,1:syms_this%nsymop)
    deallocate(trans_tmp)
    
  end subroutine init_syms

  subroutine delete_syms( syms_this )
! *****************************************************************
! Arguments:
! syms_type syms_t     : object with symops etc...
! real*8  cell(3,3)    : input lattice vectors (Bohr)
! integer na           : input number of atoms
! integer isa(na)      : input species indexes
! real*8  xa(3,na)     : input atomic cartesian coordinates (Bohr)
! *****************************************************************
    implicit         none
    type(syms_type), intent(inout) :: syms_this

    if(allocated(syms_this%symops)) deallocate (syms_this%symops)
    if(allocated(syms_this%trans)) deallocate (syms_this%trans)

  end subroutine delete_syms

end module m_syms
