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
  subroutine constr( cell, na, isa, amass, xa, stress, fa, ntcon )
! *****************************************************************
! User-written routine to implement specific geometric constraints,
! by orthogonalizing the forces and stress to undesired changes.
! Arguments:
! real*8  cell(3,3)    : input lattice vectors (Bohr)
! integer na           : input number of atoms
! integer isa(na)      : input species indexes
! real*8  amass(na)    : input atomic masses
! real*8  xa(3,na)     : input atomic cartesian coordinates (Bohr)
! real*8  stress( 3,3) : input/output stress tensor (Ry/Bohr**3)
! real*8  fa(3,na)     : input/output atomic forces (Ry/Bohr)
! integer ntcon        : total number of positions constr. imposed
! *****************************************************************
    use m_syms
    implicit         none
    integer          :: na, isa(na), ntcon
    double precision :: amass(na), cell(3,3), fa(3,na), &
      stress(3,3), xa(3,na)

! local vars
    integer :: ii


    ! call init routine with main syms_t object from m_syms
    call init_syms( syms_t, cell, na, isa, xa )

print *, "test outputs "
print *, " nsymop = ", syms_t%nsymop
print *, " symbol = ", syms_t%symbol
print *, " symops = "
do ii = 1, syms_t%nsymop
  print '(3(3I10,2x))', syms_t%symops(:,:,ii)
end do
print *, " trans = "
do ii = 1, syms_t%nsymop
  print '(3(E20.10,2x))', syms_t%trans(:,ii)
end do

! UP TO HERE!!!

! Write here your problem-specific code.



    ! clean up
    call delete_syms(syms_t)

  end subroutine constr
