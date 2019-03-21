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
! written by Matthieu Verstraete 2018
! content is minimal - the aim of this module is to isolate siesta specific 
! objects and procedures from the spglib interface

module m_syms

#ifdef SIESTA__HAS__SPGLIB
use spglib_f03

! for the moment this single instance of the syms is exported by the module.
! if you want to play with several, declare 1 per subroutine,
! and use it as the intent(out) argument to the init_syms subroutine
type(SpglibDataset_ext) :: syms_global

private

! TODO: could be incorporated into second object in m_syms, with siesta specific info
integer, allocatable :: syms_map_atoms(:,:)
integer, allocatable :: syms_inv_map_atoms(:,:) 

double precision :: syms_prec= 1.e-6 ! this is hard coded for the moment - could be input var for init_syms

public :: syms_global, syms_prec
public :: syms_equivpos, syms_pos, syms_cell, syms_forces, syms_stresses
public :: syms_print, syms_init, syms_delete

! helper functions - might be made public
private :: syms_wrapvec_zero_one, syms_wrap2zero_one
private :: syms_wrapvec_pmhalf, syms_wrap2pmhalf
private :: syms_red2cart, syms_cart2red


contains 

subroutine syms_print(syms_this)
  implicit none
  type(SpglibDataset_ext), intent(in) :: syms_this
  call spgf_print( syms_this ) 
  ! could add some output of sym maps as well, and reference symprec
end subroutine syms_print


subroutine syms_init(syms_this, cell, na, isa, xa)
  implicit none
  type(SpglibDataset_ext), intent(out) :: syms_this
  integer, intent(in)                  :: na, isa(na)
  double precision, intent(inout)      :: cell(3,3), xa(3,na)
  
  double precision :: xred(3,na)

  call syms_cart2red(cell, na, xa, xred)
  call spgf_init (syms_this, cell, na, isa, xred, syms_prec )
  call syms_equivpos(syms_this)
end subroutine syms_init


subroutine syms_delete(syms_this)
  implicit none

  type(SpglibDataset_ext), intent(inout) :: syms_this

  call spgf_delete(syms_this)

  if (allocated(syms_map_atoms)) deallocate(syms_map_atoms)
  if (allocated(syms_inv_map_atoms)) deallocate(syms_inv_map_atoms)

end subroutine syms_delete


!
! find equivalent atoms - this subroutine modifies global variables in the present module!
!
subroutine syms_equivpos(syms_this)
  implicit none
  type(SpglibDataset_ext), intent(in) :: syms_this

! local vars
  integer :: iatom, jatom, isym

! atom found by application of operation isym on atom ia
  logical :: foundsymatom
  double precision :: xredb(3)
  character (len=80) :: message

  ! should check that the syms_this is properly inited

  allocate (syms_map_atoms(syms_this%num_atom, syms_this%n_operations))
  allocate (syms_inv_map_atoms(syms_this%num_atom, syms_this%n_operations))

  syms_map_atoms = 0
  syms_inv_map_atoms = 0
  do iatom = 1, syms_this%num_atom
    do isym = 1, syms_this%n_operations
      foundsymatom = .false.
      xredb = matmul (syms_this%rotations(:,:,isym), syms_this%xred(:,iatom)) + syms_this%translations(:,isym)
      do jatom = 1, syms_this%num_atom
        if (syms_this%types(jatom) /= syms_this%types(iatom)) cycle
! the tolerance here should correspond to symprec in syms.f90
        if ( any(abs(  syms_wrapvec_pmhalf(xredb(:)-syms_this%xred(:,jatom))  ) > syms_prec) ) then
          cycle
        end if
        foundsymatom = .true.
        syms_map_atoms     (iatom, isym) = jatom
        syms_inv_map_atoms (jatom, isym) = iatom
      end do
      if (.not. foundsymatom) then
        write (message, '(a,i6,a,i6)') 'error: no symmetric for atom ', iatom, ' under space group sym ', isym 
        call die (message) 
      end if
    end do
  end do

end subroutine syms_equivpos    


! symmetrize positions according to space group we found
! NB: refine_cell in spglib is not usable as it returns the conventional cell!
subroutine syms_pos( syms_this, na, xa, cell )
  implicit none
  type(SpglibDataset_ext), intent(inout) :: syms_this
  integer,intent(in)             :: na
  double precision,intent(inout) :: xa(3,na)
  double precision,intent(inout) :: cell(3,3)

  integer :: iatom, isym, jatom, info
  integer :: pos_sym_flag(na)
  double precision :: avgpos(3)
  double precision, allocatable :: xreda_sym(:,:)

!DEBUG
 do jatom = 1, na
   print '(a,3E25.15)', 'input xreda = ', syms_this%xred(:,jatom)
 end do
!END DEBUG

! check na==syms_this%num_atom
  pos_sym_flag = 0
  allocate (xreda_sym(3,na))
  do iatom = 1, na
    ! skip already symmetrized positions
    if (pos_sym_flag(iatom) /= 0) cycle

    ! average all positions of atoms equiv to iatom
    avgpos = 0.0d0
    do isym = 1, syms_this%n_operations
! wrap to pmhalf to accumulate comparable xred values
      avgpos = avgpos + syms_wrapvec_pmhalf(matmul (syms_this%rotations(:,:,isym),& 
               syms_this%xred(:,syms_inv_map_atoms(iatom,isym))) + syms_this%translations(:,isym))
    end do
    avgpos = avgpos / syms_this%n_operations

    ! copy symmetrized positions to all equiv positions and flag them as done
    do isym = 1, syms_this%n_operations
      jatom = syms_map_atoms(iatom,isym)
      if (pos_sym_flag(jatom) /= 0) cycle
      ! use closest position by translation to the original input xred
      xreda_sym(:,jatom) = syms_this%xred(:,jatom) + syms_wrapvec_pmhalf( &
&       matmul (syms_this%rotations(:,:,isym), avgpos) + syms_this%translations(:,isym) &
&       - syms_this%xred(:,jatom) )
!DEBUG
print '(a,3E25.15)', 'change in xreda due to symmetrization = ', xreda_sym(:,jatom)-syms_this%xred(:,jatom)
!END DEBUG
      pos_sym_flag(jatom) = 1
    end do
  end do

!DEBUG
 do jatom = 1, na
   print '(a,3E25.15)', 'symmetrized xreda = ', xreda_sym(:,jatom)
 end do
!END DEBUG

  ! also update the copy inside the object
  syms_this%xred = xreda_sym
  deallocate (xreda_sym)

  ! overwrite input xa with symmetrized version
  call syms_red2cart(cell, na, syms_this%xred, xa)

end subroutine syms_pos


! symmetrize cell according to space group we found
! NB: refine_cell in spglib is not usable as it returns the conventional cell!
subroutine syms_cell( syms_this, na, xa, cell )
  implicit none
  type(SpglibDataset_ext), intent(inout) :: syms_this
  integer,intent(in)             :: na
  double precision,intent(inout) :: xa(3,na)
  double precision,intent(inout) :: cell(3,3)

  integer :: iatom, isym, jatom, info
  double precision :: avgcell(3,3)
  double precision :: trialcell(3,3)
  double precision :: itrialcell(3,3)
  double precision :: revertcell(3,3)
  double precision, allocatable :: xreda_sym(:,:)

  ! Now unit cell:
  avgcell = 0.0d0
  do isym = 1, syms_this%n_operations
    trialcell = matmul(cell, syms_this%rotations(:,:,isym))
    ! for each of the new vectors, use translations as best as possible to approach the old ones
    call INVER(trialcell,itrialcell,3,3,info)
    if (info .ne. 0) then
      call die ("inver failed in syms_cell")
    end if
    revertcell = dble(nint(matmul(itrialcell, cell)))
    trialcell = matmul(trialcell, revertcell)
    avgcell = avgcell + trialcell
  end do
  avgcell = avgcell / syms_this%n_operations
  

!DEBUG
print *, 'avgcell'
print '(3E20.10)', avgcell(:,1)
print '(3E20.10)', avgcell(:,2)
print '(3E20.10)', avgcell(:,3)
print '(a)', 'change in cell due to symmetrization = '
print '(3E20.10)', avgcell(:,1)-cell(:,1)
print '(3E20.10)', avgcell(:,2)-cell(:,2)
print '(3E20.10)', avgcell(:,3)-cell(:,3)
!END DEBUG

  cell = avgcell

  ! overwrite input with symmetrized version
  call syms_red2cart(cell, na, syms_this%xred, xa)
  syms_this%celltransp = transpose(cell)

end subroutine syms_cell



! symmetrize forces (could be any input vector, and could do this for random
! tensors too...) 
! remember: in reciprocal space the inverse transpose matrices should be used...
subroutine syms_forces(syms_this, na, fa)
  type(SpglibDataset_ext), intent(in) :: syms_this
  integer, intent(in) :: na
  double precision, intent(inout) :: fa(3,na)

  integer :: iatom, isym, jatom
  integer :: force_sym_flag(na)
  double precision :: maxdiff
  double precision :: avgforce(3)
  double precision, allocatable :: fa_sym(:,:)

  force_sym_flag = 0
  allocate (fa_sym(3,na))
  do iatom = 1, na
    ! skip already symmetrized forces
    if (force_sym_flag(iatom) /= 0) cycle

    ! average all forces on atoms equiv to iatom
    avgforce = 0.0d0
    do isym = 1, syms_this%n_operations
      avgforce = avgforce + matmul (syms_this%symops_cart(:,:,isym),& 
               fa(:,syms_inv_map_atoms(iatom,isym)))
    end do
    avgforce = avgforce / syms_this%n_operations
    ! copy symmetrized force to all equiv positions and flag them as done
    do isym = 1, syms_this%n_operations
      jatom = syms_map_atoms(iatom,isym)
      if (force_sym_flag(jatom) /= 0) cycle
      fa_sym(:,jatom) = matmul (syms_this%symops_cart(:,:,isym), avgforce)
!DEBUG
print '(a,3E30.20)', 'change in fa due to symmetrization = ', fa(:,jatom)-fa_sym(:,jatom)
!END DEBUG
      force_sym_flag(jatom) = 1
    end do
  end do
  maxdiff = maxval(abs(fa(:,:)-fa_sym(:,:)))
  print '(a,E20.10)', 'Max change in any component of fa due to symmetrization = ', maxdiff
  maxdiff = maxval(abs((fa(:,:)-fa_sym(:,:))/(1.d-16+fa(:,:))))
  print '(a,E20.10)', 'Max relative change = ', maxdiff
  fa = fa_sym
  deallocate (fa_sym)

end subroutine syms_forces

! stress
subroutine syms_stresses (syms_this, stress)
  type(SpglibDataset_ext), intent(in) :: syms_this
  double precision, intent(inout) :: stress(3,3)

  double precision :: stress_sym(3,3), dsymop(3,3), dsymopi(3,3)
  double precision :: maxdiff

  stress_sym = 0.0d0
  do isym = 1, syms_this%n_operations
    dsymop = syms_this%symops_recip_cart(:, :, isym)
    dsymopi = transpose(syms_this%symops_recip_cart(:,:,isym))
! TODO: make trials with hex and trigonal space groups to check the order of recip and normal symops
    stress_sym = stress_sym + matmul(dsymopi, matmul(stress, dsymop))
  end do
  stress_sym = stress_sym/syms_this%n_operations
  print '(a)', 'old  stress = '
  print '(3E20.10)',stress(:,1)
  print '(3E20.10)',stress(:,2)
  print '(3E20.10)',stress(:,3)
  print '(a)', 'stress with symmetrization = '
  print '(3E20.10)', stress_sym(:,1)
  print '(3E20.10)', stress_sym(:,2)
  print '(3E20.10)', stress_sym(:,3)
  print '(a)', 'change in components of stress due to symmetrization = '
  print '(3E20.10)', stress_sym(:,1)-stress(:,1)
  print '(3E20.10)', stress_sym(:,2)-stress(:,2)
  print '(3E20.10)', stress_sym(:,3)-stress(:,3)
  maxdiff = maxval(abs(stress - stress_sym))
  print '(a,E20.10)', 'Max change in any component of stress due to symmetrization = ', maxdiff

  stress = stress_sym

end subroutine syms_stresses
    

subroutine syms_red2cart(cell, num_atom, xred, xcart)
! *****************************************************************
! Arguments:
! cell : unit cell vectors
! num_atom : number of atoms
! xcart : cartesian position vectors
! xred : reduced position vectors
! *****************************************************************
    implicit         none
    double precision,intent(in) :: cell(3,3)
    integer,intent(in) :: num_atom
    double precision,intent(in) :: xred(3,num_atom)
    double precision,intent(out) :: xcart(3,num_atom)

    ! the cell input is now the habitual one with cell(:,1) the first vector
    xcart = matmul(cell, xred)

  end subroutine syms_red2cart


  subroutine syms_cart2red(cell, num_atom, xcart, xred)
! *****************************************************************
! Arguments:
! cell : unit cell vectors
! num_atom : number of atoms
! xcart : cartesian position vectors
! xred : reduced position vectors
! *****************************************************************
    implicit         none
    double precision,intent(in) :: cell(3,3)
    integer,intent(in) :: num_atom
    double precision,intent(in) :: xcart(3,num_atom)
    double precision,intent(out) :: xred(3,num_atom)

    double precision :: cellinv(3,3)
    integer :: info

    call INVER(cell,cellinv,3,3,info)
! TODO: this level of routine should not stop, but throw an error message
    if (info .ne. 0) then 
        call die ('subroutine syms_cart2red: error in inverse matrix of cell') 
    end if

    ! the cell input is the real cell with cell(:,1) the first vector
    ! and INVER does not transpose its output
    xred = matmul(cellinv,xcart)

  end subroutine syms_cart2red

  function syms_wrapvec_zero_one(num) result (red)
    implicit none
    double precision, intent(in) :: num(3)
    double precision :: red(3)

    integer :: ii

    do ii=1, 3
      red(ii) = syms_wrap2zero_one(num(ii))
    end do

  end function syms_wrapvec_zero_one

  function syms_wrap2zero_one(num) result (red)
    implicit none
    double precision, intent(in) :: num
    double precision :: red

    if (num>0.0d0) then
      red=mod((num+1.d-13),1.0d0)-1.d-13
    else
      red=-mod(-(num-1.0d0+1.d-13),1.0d0)+1.0d0-1.d-13
    end if
    if(abs(red)<1.d-13)red=0.0d0

  end function syms_wrap2zero_one

  function syms_wrapvec_pmhalf(num) result (red)
    implicit none
    double precision, intent(in) :: num(3)
    double precision :: red(3)

    integer :: ii

    do ii=1, 3
      red(ii) = syms_wrap2pmhalf(num(ii))
    end do

  end function syms_wrapvec_pmhalf

  function syms_wrap2pmhalf(num) result (red)
    implicit none
    double precision, intent(in) :: num
    double precision :: red

    red = syms_wrap2zero_one(num+0.5d0)-0.5d0

  end function syms_wrap2pmhalf
   
#endif

end module m_syms
