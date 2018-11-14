! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
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

! Write here your problem-specific code.

  end subroutine constr
