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
    integer :: ii, iatom, jatom, isym

    integer, allocatable :: equiv_atom(:,:) ! atom found by application of operation isym on atom ia
    integer, allocatable :: invequiv_atom(:,:) 
    logical :: foundsymatom
    double precision :: xredb(3)
    double precision,allocatable :: xreda(:,:)

    integer :: force_sym_flag(na)
    double precision :: avgforce(3)
    double precision,allocatable :: fa_sym(:,:)

    integer :: pos_sym_flag(na)
    double precision :: avgpos(3)
    double precision,allocatable :: xreda_sym(:,:)



    ! call init routine with main syms_t object from m_syms
    call init_syms( syms_t, cell, na, isa, xa )

    call print_syms(syms_t)
! UP TO HERE!!!

! Write here your problem-specific code.

! find equivalent atoms - this can be encapsulated later
    allocate (equiv_atom(na, syms_t%nsymop))
    allocate (invequiv_atom(na, syms_t%nsymop))
    allocate (xreda(3, na))
    equiv_atom = 0
    invequiv_atom = 0
    call cart2red(cell, na, xa, xreda)
    do iatom = 1, na
      do isym = 1, syms_t%nsymop
        foundsymatom = .false.
        xredb = matmul (syms_t%symops(:,:,isym), xreda(:,iatom)) + syms_t%trans(:,isym)
        do jatom = 1, na
          if (isa(jatom) /= isa(iatom)) cycle
! the tolerance here should correspond to symprec in syms.f90
          if (sum(abs(  wrapvec_zero_one(xredb(:)-xreda(:,jatom))  )) > 1.e-6) cycle
          foundsymatom = .true.
          equiv_atom (iatom, isym) = jatom
          invequiv_atom (jatom, isym) = iatom
        end do
        if (.not. foundsymatom) then
          print *, 'error: no symmetric for atom ', iatom, ' under space group sym ', isym 
          stop
        end if
      end do
    end do
    deallocate(xreda)
! end encapsulation  - output is in equiv_atom invequiv_atom
    
! symmetrize positions encapsulate from here 
    pos_sym_flag = 0
    allocate (xreda(3, na))
    allocate (xreda_sym(3,na))
    call cart2red(cell, na, xa, xreda)
    do iatom = 1, na
      ! skip already symmetrized positions
      if (pos_sym_flag(iatom) /= 0) cycle

      ! average all positions of atoms equiv to iatom
      avgpos = 0.0d0
      do isym = 1, syms_t%nsymop
        avgpos = avgpos + wrapvec_zero_one(matmul (syms_t%symops(:,:,isym),& 
                 xreda(:,invequiv_atom(iatom,isym))) + syms_t%trans(:,isym))
      end do
      avgpos = avgpos / syms_t%nsymop

      ! copy symmetrized position to all equiv positions and flag them as done
      do isym = 1, syms_t%nsymop
        jatom = equiv_atom(iatom,isym)
        if (pos_sym_flag(jatom) /= 0) cycle
        xreda_sym(:,jatom) = wrapvec_zero_one(matmul (syms_t%symops(:,:,isym), avgpos) + syms_t%trans(:,isym))
!DEBUG
print '(a,3E30.20)', 'change in xreda due to symmetrization = ', xreda(:,jatom)-xreda_sym(:,jatom)
!END DEBUG
        pos_sym_flag(jatom) = 1
      end do
    end do
    xreda = xreda_sym
    deallocate (xreda_sym)
    call red2cart(cell, na, xreda, xa)
    deallocate(xreda)
! end encapsulation 

! symmetrize forces (could be any input vector, and could do this for random
! tensors too...) encapsulate from here?
! remember: in reciprocal space the inverse transpose matrices should be used...
    force_sym_flag = 0
    allocate (fa_sym(3,na))
    do iatom = 1, na
      ! skip already symmetrized forces
      if (force_sym_flag(iatom) /= 0) cycle

      ! average all forces on atoms equiv to iatom
      avgforce = 0.0d0
      do isym = 1, syms_t%nsymop
        avgforce = avgforce + matmul (syms_t%symops_cart(:,:,isym),& 
                 fa(:,invequiv_atom(iatom,isym)))
      end do
      avgforce = avgforce / syms_t%nsymop
      ! copy symmetrized force to all equiv positions and flag them as done
      do isym = 1, syms_t%nsymop
        jatom = equiv_atom(iatom,isym)
        if (force_sym_flag(jatom) /= 0) cycle
        fa_sym(:,jatom) = matmul (syms_t%symops_cart(:,:,isym), avgforce)
!DEBUG
print '(a,3E30.20)', 'change in fa due to symmetrization = ', fa(:,jatom)-fa_sym(:,jatom)
!END DEBUG
        force_sym_flag(jatom) = 1
      end do
    end do
    fa = fa_sym
    deallocate (fa_sym)
! end encapsulation
    

    ! clean up
    call delete_syms(syms_t)

  end subroutine constr
