c $Id: constr.f,v 1.5 1999/02/23 12:05:22 wdpgaara Exp $

      subroutine constr( cell, na, isa, amass, xa, stress, fa )
c *****************************************************************
c User-written routine to implement specific geometric constraints,
c by orthogonalizing the forces and stress to undesired changes.
c Arguments:
c real*8  cell(3,3)    : input lattice vectors (Bohr)
c integer na           : input number of atoms
c integer isa(na)      : input species indexes
c real*8  amass(na)    : input atomic masses
c real*8  xa(3,na)     : input atomic cartesian coordinates (Bohr)
c real*8  stress( 3,3) : input/output stress tensor (Ry/Bohr**3)
c real*8  fa(3,na)     : input/output atomic forces (Ry/Bohr)
c *****************************************************************
      implicit         none
      integer          na, isa(na)
      double precision amass(na), cell(3,3), fa(3,na),
     .                 stress(3,3), xa(3,na)

c Write here your problem-specific code.

      end

