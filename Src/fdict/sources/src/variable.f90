! @LICENSE@, see README.md
! Generic purpose variable as in any scripting language
! It has the power to transform into any variable at any time
module variable
  use iso_var_str
  implicit none
  private
  integer, parameter :: ih = selected_int_kind(4)
  integer, parameter :: is = selected_int_kind(9)
  integer, parameter :: il = selected_int_kind(18)
  integer, parameter :: sp = selected_real_kind(p=6)
  integer, parameter :: dp = selected_real_kind(p=15)
  ! To create a constant transfer data-type of the
  ! pointer methods
  character(len=1) :: local_enc_type(1)
  type :: var
     character(len=4) :: t = '    '
     ! The encoding placement of all data
     character(len=1), dimension(:), allocatable :: enc
  end type var
  public :: var
! public :: size
  interface which
     module procedure which_
  end interface which
  public :: which
  interface delete
     module procedure delete_
  end interface delete
  public :: delete
  interface nullify
     module procedure nullify_
  end interface nullify
  public :: nullify
  interface print
     module procedure print_
  end interface print
  public :: print
  ! Specific routines for passing types to variables
  interface associate_type
     module procedure associate_type_
  end interface associate_type
  public :: associate_type
  interface enc
     module procedure enc_
  end interface enc
  public :: enc
  interface size_enc
     module procedure size_enc_
  end interface size_enc
  public :: size_enc
interface assign
module procedure assign_get_char0
module procedure assign_set_char0
module procedure assign_var
module procedure assign_get_V0
module procedure assign_set_V0
module procedure assign_get_s0
module procedure assign_set_s0
module procedure assign_get_s1
module procedure assign_set_s1
module procedure assign_get_s2
module procedure assign_set_s2
module procedure assign_get_d0
module procedure assign_set_d0
module procedure assign_get_d1
module procedure assign_set_d1
module procedure assign_get_d2
module procedure assign_set_d2
module procedure assign_get_c0
module procedure assign_set_c0
module procedure assign_get_c1
module procedure assign_set_c1
module procedure assign_get_c2
module procedure assign_set_c2
module procedure assign_get_z0
module procedure assign_set_z0
module procedure assign_get_z1
module procedure assign_set_z1
module procedure assign_get_z2
module procedure assign_set_z2
module procedure assign_get_b0
module procedure assign_set_b0
module procedure assign_get_b1
module procedure assign_set_b1
module procedure assign_get_b2
module procedure assign_set_b2
module procedure assign_get_h0
module procedure assign_set_h0
module procedure assign_get_h1
module procedure assign_set_h1
module procedure assign_get_h2
module procedure assign_set_h2
module procedure assign_get_i0
module procedure assign_set_i0
module procedure assign_get_i1
module procedure assign_set_i1
module procedure assign_get_i2
module procedure assign_set_i2
module procedure assign_get_l0
module procedure assign_set_l0
module procedure assign_get_l1
module procedure assign_set_l1
module procedure assign_get_l2
module procedure assign_set_l2
end interface assign
public :: assign
interface associate
module procedure associate_var
module procedure associate_get_V0
module procedure associate_set_V0
module procedure associate_get_s0
module procedure associate_set_s0
module procedure associate_get_s1
module procedure associate_set_s1
module procedure associate_get_s2
module procedure associate_set_s2
module procedure associate_get_d0
module procedure associate_set_d0
module procedure associate_get_d1
module procedure associate_set_d1
module procedure associate_get_d2
module procedure associate_set_d2
module procedure associate_get_c0
module procedure associate_set_c0
module procedure associate_get_c1
module procedure associate_set_c1
module procedure associate_get_c2
module procedure associate_set_c2
module procedure associate_get_z0
module procedure associate_set_z0
module procedure associate_get_z1
module procedure associate_set_z1
module procedure associate_get_z2
module procedure associate_set_z2
module procedure associate_get_b0
module procedure associate_set_b0
module procedure associate_get_b1
module procedure associate_set_b1
module procedure associate_get_b2
module procedure associate_set_b2
module procedure associate_get_h0
module procedure associate_set_h0
module procedure associate_get_h1
module procedure associate_set_h1
module procedure associate_get_h2
module procedure associate_set_h2
module procedure associate_get_i0
module procedure associate_set_i0
module procedure associate_get_i1
module procedure associate_set_i1
module procedure associate_get_i2
module procedure associate_set_i2
module procedure associate_get_l0
module procedure associate_set_l0
module procedure associate_get_l1
module procedure associate_set_l1
module procedure associate_get_l2
module procedure associate_set_l2
end interface associate
public :: associate
interface associatd
module procedure associatd_l_V0
module procedure associatd_r_V0
module procedure associatd_l_s0
module procedure associatd_r_s0
module procedure associatd_l_s1
module procedure associatd_r_s1
module procedure associatd_l_s2
module procedure associatd_r_s2
module procedure associatd_l_d0
module procedure associatd_r_d0
module procedure associatd_l_d1
module procedure associatd_r_d1
module procedure associatd_l_d2
module procedure associatd_r_d2
module procedure associatd_l_c0
module procedure associatd_r_c0
module procedure associatd_l_c1
module procedure associatd_r_c1
module procedure associatd_l_c2
module procedure associatd_r_c2
module procedure associatd_l_z0
module procedure associatd_r_z0
module procedure associatd_l_z1
module procedure associatd_r_z1
module procedure associatd_l_z2
module procedure associatd_r_z2
module procedure associatd_l_b0
module procedure associatd_r_b0
module procedure associatd_l_b1
module procedure associatd_r_b1
module procedure associatd_l_b2
module procedure associatd_r_b2
module procedure associatd_l_h0
module procedure associatd_r_h0
module procedure associatd_l_h1
module procedure associatd_r_h1
module procedure associatd_l_h2
module procedure associatd_r_h2
module procedure associatd_l_i0
module procedure associatd_r_i0
module procedure associatd_l_i1
module procedure associatd_r_i1
module procedure associatd_l_i2
module procedure associatd_r_i2
module procedure associatd_l_l0
module procedure associatd_r_l0
module procedure associatd_l_l1
module procedure associatd_r_l1
module procedure associatd_l_l2
module procedure associatd_r_l2
end interface associatd
public :: associatd
contains
  subroutine print_(this)
    type(var), intent(in) :: this
    write(*,'(t2,a)') this%t
  end subroutine print_
  elemental function which_(this) result(t)
    type(var), intent(in) :: this
    character(len=4) :: t
    t = this%t
  end function which_
  subroutine delete_(this,dealloc)
    type(var), intent(inout) :: this
    logical, intent(in), optional :: dealloc
    logical :: ldealloc
type :: ptV0
 type(var_str), pointer :: p => null()
end type ptV0
type(ptV0) :: pV0
type :: pts0
 real(sp), pointer :: p => null()
end type pts0
type(pts0) :: ps0
type :: pts1
 real(sp), pointer :: p(:) => null()
end type pts1
type(pts1) :: ps1
type :: pts2
 real(sp), pointer :: p(:,:) => null()
end type pts2
type(pts2) :: ps2
type :: ptd0
 real(dp), pointer :: p => null()
end type ptd0
type(ptd0) :: pd0
type :: ptd1
 real(dp), pointer :: p(:) => null()
end type ptd1
type(ptd1) :: pd1
type :: ptd2
 real(dp), pointer :: p(:,:) => null()
end type ptd2
type(ptd2) :: pd2
type :: ptc0
 complex(sp), pointer :: p => null()
end type ptc0
type(ptc0) :: pc0
type :: ptc1
 complex(sp), pointer :: p(:) => null()
end type ptc1
type(ptc1) :: pc1
type :: ptc2
 complex(sp), pointer :: p(:,:) => null()
end type ptc2
type(ptc2) :: pc2
type :: ptz0
 complex(dp), pointer :: p => null()
end type ptz0
type(ptz0) :: pz0
type :: ptz1
 complex(dp), pointer :: p(:) => null()
end type ptz1
type(ptz1) :: pz1
type :: ptz2
 complex(dp), pointer :: p(:,:) => null()
end type ptz2
type(ptz2) :: pz2
type :: ptb0
 logical, pointer :: p => null()
end type ptb0
type(ptb0) :: pb0
type :: ptb1
 logical, pointer :: p(:) => null()
end type ptb1
type(ptb1) :: pb1
type :: ptb2
 logical, pointer :: p(:,:) => null()
end type ptb2
type(ptb2) :: pb2
type :: pth0
 integer(ih), pointer :: p => null()
end type pth0
type(pth0) :: ph0
type :: pth1
 integer(ih), pointer :: p(:) => null()
end type pth1
type(pth1) :: ph1
type :: pth2
 integer(ih), pointer :: p(:,:) => null()
end type pth2
type(pth2) :: ph2
type :: pti0
 integer(is), pointer :: p => null()
end type pti0
type(pti0) :: pi0
type :: pti1
 integer(is), pointer :: p(:) => null()
end type pti1
type(pti1) :: pi1
type :: pti2
 integer(is), pointer :: p(:,:) => null()
end type pti2
type(pti2) :: pi2
type :: ptl0
 integer(il), pointer :: p => null()
end type ptl0
type(ptl0) :: pl0
type :: ptl1
 integer(il), pointer :: p(:) => null()
end type ptl1
type(ptl1) :: pl1
type :: ptl2
 integer(il), pointer :: p(:,:) => null()
end type ptl2
type(ptl2) :: pl2
    ldealloc = .true.
    if ( present(dealloc) ) ldealloc = dealloc
    if ( ldealloc ) then
if (this%t == 'V0') then
  pV0 = transfer(this%enc,pV0)
  pV0%p = ''
  deallocate(pV0%p)
end if
if (this%t == 's0') then
  ps0 = transfer(this%enc,ps0)
  deallocate(ps0%p)
end if
if (this%t == 's1') then
  ps1 = transfer(this%enc,ps1)
  deallocate(ps1%p)
end if
if (this%t == 's2') then
  ps2 = transfer(this%enc,ps2)
  deallocate(ps2%p)
end if
if (this%t == 'd0') then
  pd0 = transfer(this%enc,pd0)
  deallocate(pd0%p)
end if
if (this%t == 'd1') then
  pd1 = transfer(this%enc,pd1)
  deallocate(pd1%p)
end if
if (this%t == 'd2') then
  pd2 = transfer(this%enc,pd2)
  deallocate(pd2%p)
end if
if (this%t == 'c0') then
  pc0 = transfer(this%enc,pc0)
  deallocate(pc0%p)
end if
if (this%t == 'c1') then
  pc1 = transfer(this%enc,pc1)
  deallocate(pc1%p)
end if
if (this%t == 'c2') then
  pc2 = transfer(this%enc,pc2)
  deallocate(pc2%p)
end if
if (this%t == 'z0') then
  pz0 = transfer(this%enc,pz0)
  deallocate(pz0%p)
end if
if (this%t == 'z1') then
  pz1 = transfer(this%enc,pz1)
  deallocate(pz1%p)
end if
if (this%t == 'z2') then
  pz2 = transfer(this%enc,pz2)
  deallocate(pz2%p)
end if
if (this%t == 'b0') then
  pb0 = transfer(this%enc,pb0)
  deallocate(pb0%p)
end if
if (this%t == 'b1') then
  pb1 = transfer(this%enc,pb1)
  deallocate(pb1%p)
end if
if (this%t == 'b2') then
  pb2 = transfer(this%enc,pb2)
  deallocate(pb2%p)
end if
if (this%t == 'h0') then
  ph0 = transfer(this%enc,ph0)
  deallocate(ph0%p)
end if
if (this%t == 'h1') then
  ph1 = transfer(this%enc,ph1)
  deallocate(ph1%p)
end if
if (this%t == 'h2') then
  ph2 = transfer(this%enc,ph2)
  deallocate(ph2%p)
end if
if (this%t == 'i0') then
  pi0 = transfer(this%enc,pi0)
  deallocate(pi0%p)
end if
if (this%t == 'i1') then
  pi1 = transfer(this%enc,pi1)
  deallocate(pi1%p)
end if
if (this%t == 'i2') then
  pi2 = transfer(this%enc,pi2)
  deallocate(pi2%p)
end if
if (this%t == 'l0') then
  pl0 = transfer(this%enc,pl0)
  deallocate(pl0%p)
end if
if (this%t == 'l1') then
  pl1 = transfer(this%enc,pl1)
  deallocate(pl1%p)
end if
if (this%t == 'l2') then
  pl2 = transfer(this%enc,pl2)
  deallocate(pl2%p)
end if
if (this%t == 'USER') then
print '(a)','var: Cannot deallocate UT, proceed:'
print '(a)','     1) retrieve type, 2) deallocate, 3) call nullify(var)'
end if
    end if
    call nullify(this)
  end subroutine delete_
  elemental subroutine nullify_(this)
    type(var), intent(inout) :: this
    this%t = '  '
    if ( allocated(this%enc) ) deallocate(this%enc)
  end subroutine nullify_
  ! Returns the bare encoding of this variable
  ! This can ease the process of assigning
  ! user-types to a variable.
  ! An encoding might be 2, or 10000000 bytes big.
  ! Therefore we use a subroutine to determine
  ! the size of the returning encoding characters.
  ! If the size of the returning enc is not
  ! big enough it will be reset to ' '
  subroutine enc_(this,enc)
    type(var), intent(in) :: this
    character(len=1), intent(out) :: enc(:)
    integer :: i
    if ( this%t == '  ' ) then
       enc = ' '
    else
       ! We do have an encoding
       i = size(this%enc)
       if ( i > size(enc) ) then
          enc = ' '
       else
          enc(1:i) = this%enc
       end if
    end if
  end subroutine enc_
  function size_enc_(this) result(len)
    type(var), intent(in) :: this
    integer :: len
    if ( this%t == '  ' ) then
       len = 0
    else
       len = size(this%enc)
    end if
  end function size_enc_
  ! We allow the user to pass an encoding field.
  ! As this is the same as passing a char
  ! we MUST use a specific routine for this.
  ! One _could_, in principle, add an optional
  ! logical flag for the assign_var_char0, however
  ! one cannot assign a type by passing a reference
  ! and hence we ONLY allow associate_type
  ! This also means that any de-allocation of variables
  ! containing an external type will only de-reference it.
  ! A bit counter-intuitive, yet the variable type needs
  ! all information about the type to successfully de-allocate it.
  ! It is ALSO very important that the user
  ! passed the full-encoding WITHOUT padding of ' '.
  ! We cannot know for sure whether the encoding actually terminates
  ! in a bit corresponding to char(' ')!
  subroutine associate_type_(this,enc,dealloc)
    type(var), intent(inout) :: this
    character(len=1), intent(in) :: enc(:)
    logical, intent(in), optional :: dealloc
    logical :: ldealloc
    ldealloc = .false.
    if(present(dealloc))ldealloc = dealloc
    if (.not. ldealloc) then
       ! if we do not deallocate, nullify
       call nullify(this)
    else
       call delete(this)
    end if
    this%t = 'USER'
    allocate(this%enc(size(enc)))
    this%enc = enc
  end subroutine associate_type_
  subroutine assign_var(this,rhs,dealloc)
    type(var), intent(inout) :: this
    type(var), intent(in) :: rhs
    logical, intent(in), optional :: dealloc
    logical :: ldealloc
type :: ptV0
 type(var_str), pointer :: p => null()
end type ptV0
type(ptV0) :: pV0_1, pV0_2
type :: pts0
 real(sp), pointer :: p => null()
end type pts0
type(pts0) :: ps0_1, ps0_2
type :: pts1
 real(sp), pointer :: p(:) => null()
end type pts1
type(pts1) :: ps1_1, ps1_2
type :: pts2
 real(sp), pointer :: p(:,:) => null()
end type pts2
type(pts2) :: ps2_1, ps2_2
type :: ptd0
 real(dp), pointer :: p => null()
end type ptd0
type(ptd0) :: pd0_1, pd0_2
type :: ptd1
 real(dp), pointer :: p(:) => null()
end type ptd1
type(ptd1) :: pd1_1, pd1_2
type :: ptd2
 real(dp), pointer :: p(:,:) => null()
end type ptd2
type(ptd2) :: pd2_1, pd2_2
type :: ptc0
 complex(sp), pointer :: p => null()
end type ptc0
type(ptc0) :: pc0_1, pc0_2
type :: ptc1
 complex(sp), pointer :: p(:) => null()
end type ptc1
type(ptc1) :: pc1_1, pc1_2
type :: ptc2
 complex(sp), pointer :: p(:,:) => null()
end type ptc2
type(ptc2) :: pc2_1, pc2_2
type :: ptz0
 complex(dp), pointer :: p => null()
end type ptz0
type(ptz0) :: pz0_1, pz0_2
type :: ptz1
 complex(dp), pointer :: p(:) => null()
end type ptz1
type(ptz1) :: pz1_1, pz1_2
type :: ptz2
 complex(dp), pointer :: p(:,:) => null()
end type ptz2
type(ptz2) :: pz2_1, pz2_2
type :: ptb0
 logical, pointer :: p => null()
end type ptb0
type(ptb0) :: pb0_1, pb0_2
type :: ptb1
 logical, pointer :: p(:) => null()
end type ptb1
type(ptb1) :: pb1_1, pb1_2
type :: ptb2
 logical, pointer :: p(:,:) => null()
end type ptb2
type(ptb2) :: pb2_1, pb2_2
type :: pth0
 integer(ih), pointer :: p => null()
end type pth0
type(pth0) :: ph0_1, ph0_2
type :: pth1
 integer(ih), pointer :: p(:) => null()
end type pth1
type(pth1) :: ph1_1, ph1_2
type :: pth2
 integer(ih), pointer :: p(:,:) => null()
end type pth2
type(pth2) :: ph2_1, ph2_2
type :: pti0
 integer(is), pointer :: p => null()
end type pti0
type(pti0) :: pi0_1, pi0_2
type :: pti1
 integer(is), pointer :: p(:) => null()
end type pti1
type(pti1) :: pi1_1, pi1_2
type :: pti2
 integer(is), pointer :: p(:,:) => null()
end type pti2
type(pti2) :: pi2_1, pi2_2
type :: ptl0
 integer(il), pointer :: p => null()
end type ptl0
type(ptl0) :: pl0_1, pl0_2
type :: ptl1
 integer(il), pointer :: p(:) => null()
end type ptl1
type(ptl1) :: pl1_1, pl1_2
type :: ptl2
 integer(il), pointer :: p(:,:) => null()
end type ptl2
type(ptl2) :: pl2_1, pl2_2
    ! collect deallocation option (default as =)
    ! ASSIGNMENT in fortran is per default destructive
    ldealloc = .true.
    if(present(dealloc))ldealloc = dealloc
    if (.not. ldealloc) then
       ! if we do not deallocate, nullify
       call nullify(this)
    else
       call delete(this)
    end if
    this%t = rhs%t
    ! First allocate the LHS
if ( this%t == 'V0' ) then
pV0_2 = transfer(rhs%enc,pV0_2)
allocate(pV0_1%p)
endif
if ( this%t == 's0' ) then
ps0_2 = transfer(rhs%enc,ps0_2)
allocate(ps0_1%p)
elseif ( this%t == 's1' ) then
ps1_2 = transfer(rhs%enc,ps1_2)
allocate(ps1_1%p(size(ps1_2%p)))
elseif ( this%t == 's2' ) then
ps2_2 = transfer(rhs%enc,ps2_2)
allocate(ps2_1%p(size(ps2_2%p,1),size(ps2_2%p,2)))
endif
if ( this%t == 'd0' ) then
pd0_2 = transfer(rhs%enc,pd0_2)
allocate(pd0_1%p)
elseif ( this%t == 'd1' ) then
pd1_2 = transfer(rhs%enc,pd1_2)
allocate(pd1_1%p(size(pd1_2%p)))
elseif ( this%t == 'd2' ) then
pd2_2 = transfer(rhs%enc,pd2_2)
allocate(pd2_1%p(size(pd2_2%p,1),size(pd2_2%p,2)))
endif
if ( this%t == 'c0' ) then
pc0_2 = transfer(rhs%enc,pc0_2)
allocate(pc0_1%p)
elseif ( this%t == 'c1' ) then
pc1_2 = transfer(rhs%enc,pc1_2)
allocate(pc1_1%p(size(pc1_2%p)))
elseif ( this%t == 'c2' ) then
pc2_2 = transfer(rhs%enc,pc2_2)
allocate(pc2_1%p(size(pc2_2%p,1),size(pc2_2%p,2)))
endif
if ( this%t == 'z0' ) then
pz0_2 = transfer(rhs%enc,pz0_2)
allocate(pz0_1%p)
elseif ( this%t == 'z1' ) then
pz1_2 = transfer(rhs%enc,pz1_2)
allocate(pz1_1%p(size(pz1_2%p)))
elseif ( this%t == 'z2' ) then
pz2_2 = transfer(rhs%enc,pz2_2)
allocate(pz2_1%p(size(pz2_2%p,1),size(pz2_2%p,2)))
endif
if ( this%t == 'b0' ) then
pb0_2 = transfer(rhs%enc,pb0_2)
allocate(pb0_1%p)
elseif ( this%t == 'b1' ) then
pb1_2 = transfer(rhs%enc,pb1_2)
allocate(pb1_1%p(size(pb1_2%p)))
elseif ( this%t == 'b2' ) then
pb2_2 = transfer(rhs%enc,pb2_2)
allocate(pb2_1%p(size(pb2_2%p,1),size(pb2_2%p,2)))
endif
if ( this%t == 'h0' ) then
ph0_2 = transfer(rhs%enc,ph0_2)
allocate(ph0_1%p)
elseif ( this%t == 'h1' ) then
ph1_2 = transfer(rhs%enc,ph1_2)
allocate(ph1_1%p(size(ph1_2%p)))
elseif ( this%t == 'h2' ) then
ph2_2 = transfer(rhs%enc,ph2_2)
allocate(ph2_1%p(size(ph2_2%p,1),size(ph2_2%p,2)))
endif
if ( this%t == 'i0' ) then
pi0_2 = transfer(rhs%enc,pi0_2)
allocate(pi0_1%p)
elseif ( this%t == 'i1' ) then
pi1_2 = transfer(rhs%enc,pi1_2)
allocate(pi1_1%p(size(pi1_2%p)))
elseif ( this%t == 'i2' ) then
pi2_2 = transfer(rhs%enc,pi2_2)
allocate(pi2_1%p(size(pi2_2%p,1),size(pi2_2%p,2)))
endif
if ( this%t == 'l0' ) then
pl0_2 = transfer(rhs%enc,pl0_2)
allocate(pl0_1%p)
elseif ( this%t == 'l1' ) then
pl1_2 = transfer(rhs%enc,pl1_2)
allocate(pl1_1%p(size(pl1_2%p)))
elseif ( this%t == 'l2' ) then
pl2_2 = transfer(rhs%enc,pl2_2)
allocate(pl2_1%p(size(pl2_2%p,1),size(pl2_2%p,2)))
endif
    ! copy over RHS and Save encoding
if ( this%t == 'V0' ) then
pV0_1%p = pV0_2%p
allocate(this%enc(size(transfer(pV0_1, local_enc_type))))
this%enc = transfer(pV0_1, local_enc_type)
endif
if ( this%t == 's0' ) then
ps0_1%p = ps0_2%p
allocate(this%enc(size(transfer(ps0_1, local_enc_type))))
this%enc = transfer(ps0_1, local_enc_type)
elseif ( this%t == 's1' ) then
ps1_1%p = ps1_2%p
allocate(this%enc(size(transfer(ps1_1, local_enc_type))))
this%enc = transfer(ps1_1, local_enc_type)
elseif ( this%t == 's2' ) then
ps2_1%p = ps2_2%p
allocate(this%enc(size(transfer(ps2_1, local_enc_type))))
this%enc = transfer(ps2_1, local_enc_type)
endif
if ( this%t == 'd0' ) then
pd0_1%p = pd0_2%p
allocate(this%enc(size(transfer(pd0_1, local_enc_type))))
this%enc = transfer(pd0_1, local_enc_type)
elseif ( this%t == 'd1' ) then
pd1_1%p = pd1_2%p
allocate(this%enc(size(transfer(pd1_1, local_enc_type))))
this%enc = transfer(pd1_1, local_enc_type)
elseif ( this%t == 'd2' ) then
pd2_1%p = pd2_2%p
allocate(this%enc(size(transfer(pd2_1, local_enc_type))))
this%enc = transfer(pd2_1, local_enc_type)
endif
if ( this%t == 'c0' ) then
pc0_1%p = pc0_2%p
allocate(this%enc(size(transfer(pc0_1, local_enc_type))))
this%enc = transfer(pc0_1, local_enc_type)
elseif ( this%t == 'c1' ) then
pc1_1%p = pc1_2%p
allocate(this%enc(size(transfer(pc1_1, local_enc_type))))
this%enc = transfer(pc1_1, local_enc_type)
elseif ( this%t == 'c2' ) then
pc2_1%p = pc2_2%p
allocate(this%enc(size(transfer(pc2_1, local_enc_type))))
this%enc = transfer(pc2_1, local_enc_type)
endif
if ( this%t == 'z0' ) then
pz0_1%p = pz0_2%p
allocate(this%enc(size(transfer(pz0_1, local_enc_type))))
this%enc = transfer(pz0_1, local_enc_type)
elseif ( this%t == 'z1' ) then
pz1_1%p = pz1_2%p
allocate(this%enc(size(transfer(pz1_1, local_enc_type))))
this%enc = transfer(pz1_1, local_enc_type)
elseif ( this%t == 'z2' ) then
pz2_1%p = pz2_2%p
allocate(this%enc(size(transfer(pz2_1, local_enc_type))))
this%enc = transfer(pz2_1, local_enc_type)
endif
if ( this%t == 'b0' ) then
pb0_1%p = pb0_2%p
allocate(this%enc(size(transfer(pb0_1, local_enc_type))))
this%enc = transfer(pb0_1, local_enc_type)
elseif ( this%t == 'b1' ) then
pb1_1%p = pb1_2%p
allocate(this%enc(size(transfer(pb1_1, local_enc_type))))
this%enc = transfer(pb1_1, local_enc_type)
elseif ( this%t == 'b2' ) then
pb2_1%p = pb2_2%p
allocate(this%enc(size(transfer(pb2_1, local_enc_type))))
this%enc = transfer(pb2_1, local_enc_type)
endif
if ( this%t == 'h0' ) then
ph0_1%p = ph0_2%p
allocate(this%enc(size(transfer(ph0_1, local_enc_type))))
this%enc = transfer(ph0_1, local_enc_type)
elseif ( this%t == 'h1' ) then
ph1_1%p = ph1_2%p
allocate(this%enc(size(transfer(ph1_1, local_enc_type))))
this%enc = transfer(ph1_1, local_enc_type)
elseif ( this%t == 'h2' ) then
ph2_1%p = ph2_2%p
allocate(this%enc(size(transfer(ph2_1, local_enc_type))))
this%enc = transfer(ph2_1, local_enc_type)
endif
if ( this%t == 'i0' ) then
pi0_1%p = pi0_2%p
allocate(this%enc(size(transfer(pi0_1, local_enc_type))))
this%enc = transfer(pi0_1, local_enc_type)
elseif ( this%t == 'i1' ) then
pi1_1%p = pi1_2%p
allocate(this%enc(size(transfer(pi1_1, local_enc_type))))
this%enc = transfer(pi1_1, local_enc_type)
elseif ( this%t == 'i2' ) then
pi2_1%p = pi2_2%p
allocate(this%enc(size(transfer(pi2_1, local_enc_type))))
this%enc = transfer(pi2_1, local_enc_type)
endif
if ( this%t == 'l0' ) then
pl0_1%p = pl0_2%p
allocate(this%enc(size(transfer(pl0_1, local_enc_type))))
this%enc = transfer(pl0_1, local_enc_type)
elseif ( this%t == 'l1' ) then
pl1_1%p = pl1_2%p
allocate(this%enc(size(transfer(pl1_1, local_enc_type))))
this%enc = transfer(pl1_1, local_enc_type)
elseif ( this%t == 'l2' ) then
pl2_1%p = pl2_2%p
allocate(this%enc(size(transfer(pl2_1, local_enc_type))))
this%enc = transfer(pl2_1, local_enc_type)
endif
if ( this%t == 'USER' ) then
print '(a)','var: Cannot assign a UT, USE call associate(..)'
end if
  end subroutine assign_var
  subroutine associate_var(this,rhs,dealloc,success)
    type(var), intent(inout) :: this
    type(var), intent(in) :: rhs
    logical, intent(in), optional :: dealloc
    logical, intent(out), optional :: success
    logical :: ldealloc
    ! collect deallocation option (default as =)
    ! ASSOCIATION in fortran is per default non-destructive
    ldealloc = .false.
    if ( present(success) ) success = .true.
    if ( present(dealloc) ) ldealloc = dealloc
    if (.not. ldealloc) then
       ! if we do not deallocate, nullify
       call nullify(this)
    else
       call delete(this)
    end if
    ! Association is done by copying the encoding
    this%t = rhs%t
    allocate(this%enc(size(rhs%enc)))
    this%enc = rhs%enc
  end subroutine associate_var
  pure function associatd_var(this,rhs) result(ret)
    type(var), intent(in) :: this
    type(var), intent(in) :: rhs
    logical :: ret
type :: ptV0
 type(var_str), pointer :: p => null()
end type ptV0
type(ptV0) :: pV0_1, pV0_2
type :: pts0
 real(sp), pointer :: p => null()
end type pts0
type(pts0) :: ps0_1, ps0_2
type :: pts1
 real(sp), pointer :: p(:) => null()
end type pts1
type(pts1) :: ps1_1, ps1_2
type :: pts2
 real(sp), pointer :: p(:,:) => null()
end type pts2
type(pts2) :: ps2_1, ps2_2
type :: ptd0
 real(dp), pointer :: p => null()
end type ptd0
type(ptd0) :: pd0_1, pd0_2
type :: ptd1
 real(dp), pointer :: p(:) => null()
end type ptd1
type(ptd1) :: pd1_1, pd1_2
type :: ptd2
 real(dp), pointer :: p(:,:) => null()
end type ptd2
type(ptd2) :: pd2_1, pd2_2
type :: ptc0
 complex(sp), pointer :: p => null()
end type ptc0
type(ptc0) :: pc0_1, pc0_2
type :: ptc1
 complex(sp), pointer :: p(:) => null()
end type ptc1
type(ptc1) :: pc1_1, pc1_2
type :: ptc2
 complex(sp), pointer :: p(:,:) => null()
end type ptc2
type(ptc2) :: pc2_1, pc2_2
type :: ptz0
 complex(dp), pointer :: p => null()
end type ptz0
type(ptz0) :: pz0_1, pz0_2
type :: ptz1
 complex(dp), pointer :: p(:) => null()
end type ptz1
type(ptz1) :: pz1_1, pz1_2
type :: ptz2
 complex(dp), pointer :: p(:,:) => null()
end type ptz2
type(ptz2) :: pz2_1, pz2_2
type :: ptb0
 logical, pointer :: p => null()
end type ptb0
type(ptb0) :: pb0_1, pb0_2
type :: ptb1
 logical, pointer :: p(:) => null()
end type ptb1
type(ptb1) :: pb1_1, pb1_2
type :: ptb2
 logical, pointer :: p(:,:) => null()
end type ptb2
type(ptb2) :: pb2_1, pb2_2
type :: pth0
 integer(ih), pointer :: p => null()
end type pth0
type(pth0) :: ph0_1, ph0_2
type :: pth1
 integer(ih), pointer :: p(:) => null()
end type pth1
type(pth1) :: ph1_1, ph1_2
type :: pth2
 integer(ih), pointer :: p(:,:) => null()
end type pth2
type(pth2) :: ph2_1, ph2_2
type :: pti0
 integer(is), pointer :: p => null()
end type pti0
type(pti0) :: pi0_1, pi0_2
type :: pti1
 integer(is), pointer :: p(:) => null()
end type pti1
type(pti1) :: pi1_1, pi1_2
type :: pti2
 integer(is), pointer :: p(:,:) => null()
end type pti2
type(pti2) :: pi2_1, pi2_2
type :: ptl0
 integer(il), pointer :: p => null()
end type ptl0
type(ptl0) :: pl0_1, pl0_2
type :: ptl1
 integer(il), pointer :: p(:) => null()
end type ptl1
type(ptl1) :: pl1_1, pl1_2
type :: ptl2
 integer(il), pointer :: p(:,:) => null()
end type ptl2
type(ptl2) :: pl2_1, pl2_2
    ret = this%t==rhs%t
    if ( .not. ret ) return
if ( this%t == 'V0' ) then
pV0_1 = transfer(this%enc,pV0_1)
pV0_2 = transfer(rhs%enc,pV0_2)
ret = associated(pV0_1%p,pV0_2%p)
endif
if ( this%t == 's0' ) then
ps0_1 = transfer(this%enc,ps0_1)
ps0_2 = transfer(rhs%enc,ps0_2)
ret = associated(ps0_1%p,ps0_2%p)
elseif ( this%t == 's1' ) then
ps1_1 = transfer(this%enc,ps1_1)
ps1_2 = transfer(rhs%enc,ps1_2)
ret = associated(ps1_1%p,ps1_2%p)
elseif ( this%t == 's2' ) then
ps2_1 = transfer(this%enc,ps2_1)
ps2_2 = transfer(rhs%enc,ps2_2)
ret = associated(ps2_1%p,ps2_2%p)
endif
if ( this%t == 'd0' ) then
pd0_1 = transfer(this%enc,pd0_1)
pd0_2 = transfer(rhs%enc,pd0_2)
ret = associated(pd0_1%p,pd0_2%p)
elseif ( this%t == 'd1' ) then
pd1_1 = transfer(this%enc,pd1_1)
pd1_2 = transfer(rhs%enc,pd1_2)
ret = associated(pd1_1%p,pd1_2%p)
elseif ( this%t == 'd2' ) then
pd2_1 = transfer(this%enc,pd2_1)
pd2_2 = transfer(rhs%enc,pd2_2)
ret = associated(pd2_1%p,pd2_2%p)
endif
if ( this%t == 'c0' ) then
pc0_1 = transfer(this%enc,pc0_1)
pc0_2 = transfer(rhs%enc,pc0_2)
ret = associated(pc0_1%p,pc0_2%p)
elseif ( this%t == 'c1' ) then
pc1_1 = transfer(this%enc,pc1_1)
pc1_2 = transfer(rhs%enc,pc1_2)
ret = associated(pc1_1%p,pc1_2%p)
elseif ( this%t == 'c2' ) then
pc2_1 = transfer(this%enc,pc2_1)
pc2_2 = transfer(rhs%enc,pc2_2)
ret = associated(pc2_1%p,pc2_2%p)
endif
if ( this%t == 'z0' ) then
pz0_1 = transfer(this%enc,pz0_1)
pz0_2 = transfer(rhs%enc,pz0_2)
ret = associated(pz0_1%p,pz0_2%p)
elseif ( this%t == 'z1' ) then
pz1_1 = transfer(this%enc,pz1_1)
pz1_2 = transfer(rhs%enc,pz1_2)
ret = associated(pz1_1%p,pz1_2%p)
elseif ( this%t == 'z2' ) then
pz2_1 = transfer(this%enc,pz2_1)
pz2_2 = transfer(rhs%enc,pz2_2)
ret = associated(pz2_1%p,pz2_2%p)
endif
if ( this%t == 'b0' ) then
pb0_1 = transfer(this%enc,pb0_1)
pb0_2 = transfer(rhs%enc,pb0_2)
ret = associated(pb0_1%p,pb0_2%p)
elseif ( this%t == 'b1' ) then
pb1_1 = transfer(this%enc,pb1_1)
pb1_2 = transfer(rhs%enc,pb1_2)
ret = associated(pb1_1%p,pb1_2%p)
elseif ( this%t == 'b2' ) then
pb2_1 = transfer(this%enc,pb2_1)
pb2_2 = transfer(rhs%enc,pb2_2)
ret = associated(pb2_1%p,pb2_2%p)
endif
if ( this%t == 'h0' ) then
ph0_1 = transfer(this%enc,ph0_1)
ph0_2 = transfer(rhs%enc,ph0_2)
ret = associated(ph0_1%p,ph0_2%p)
elseif ( this%t == 'h1' ) then
ph1_1 = transfer(this%enc,ph1_1)
ph1_2 = transfer(rhs%enc,ph1_2)
ret = associated(ph1_1%p,ph1_2%p)
elseif ( this%t == 'h2' ) then
ph2_1 = transfer(this%enc,ph2_1)
ph2_2 = transfer(rhs%enc,ph2_2)
ret = associated(ph2_1%p,ph2_2%p)
endif
if ( this%t == 'i0' ) then
pi0_1 = transfer(this%enc,pi0_1)
pi0_2 = transfer(rhs%enc,pi0_2)
ret = associated(pi0_1%p,pi0_2%p)
elseif ( this%t == 'i1' ) then
pi1_1 = transfer(this%enc,pi1_1)
pi1_2 = transfer(rhs%enc,pi1_2)
ret = associated(pi1_1%p,pi1_2%p)
elseif ( this%t == 'i2' ) then
pi2_1 = transfer(this%enc,pi2_1)
pi2_2 = transfer(rhs%enc,pi2_2)
ret = associated(pi2_1%p,pi2_2%p)
endif
if ( this%t == 'l0' ) then
pl0_1 = transfer(this%enc,pl0_1)
pl0_2 = transfer(rhs%enc,pl0_2)
ret = associated(pl0_1%p,pl0_2%p)
elseif ( this%t == 'l1' ) then
pl1_1 = transfer(this%enc,pl1_1)
pl1_2 = transfer(rhs%enc,pl1_2)
ret = associated(pl1_1%p,pl1_2%p)
elseif ( this%t == 'l2' ) then
pl2_1 = transfer(this%enc,pl2_1)
pl2_2 = transfer(rhs%enc,pl2_2)
ret = associated(pl2_1%p,pl2_2%p)
endif
if ( this%t == 'USER' ) then
ret = all(this%enc == rhs%enc)
end if
  end function associatd_var
  subroutine assign_set_char0(this,rhs,dealloc)
    type(var), intent(inout) :: this
    character(len=*), intent(in) :: rhs
    logical, intent(in), optional :: dealloc
    type(var_str) :: str
    str = rhs
    call assign(this,str,dealloc=dealloc)
    str = "" ! deallocation
  end subroutine assign_set_char0
  subroutine assign_get_char0(lhs,this,success)
    character(len=*), intent(out) :: lhs
    type(var), intent(inout) :: this
    logical, intent(out), optional :: success
    type(var_str) :: str
    logical :: lsuccess
    call assign(str,this,success=lsuccess)
    if ( present(success) ) success = lsuccess
    if ( lsuccess ) lhs = str
    str = "" ! deallocation
  end subroutine assign_get_char0
subroutine assign_set_V0(this,rhs,dealloc)
  type(var), intent(inout) :: this
  type(var_str), intent(in) :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    type(var_str), pointer :: p => null()
  end type
  type(pt) :: p
  ! ASSIGNMENT in fortran is per default destructive
  ldealloc = .true.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  ! With pointer transfer we need to deallocate
  ! else bounds might change...
  this%t = "V0"
  allocate(p%p) ! allocate space
  p%p = rhs ! copy data over
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
  ! We already have shipped it
  nullify(p%p)
end subroutine assign_set_V0
subroutine assign_get_V0(lhs,this,success)
  type(var_str), intent(out) :: lhs
  type(var), intent(in) :: this
  logical, intent(out), optional :: success
  logical :: lsuccess
  type :: pt
    type(var_str), pointer :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "V0"
  if (present(success)) success = lsuccess
  if (.not. lsuccess) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs = p%p
end subroutine assign_get_V0
subroutine associate_get_V0(lhs,this,dealloc,success)
  type(var_str), pointer :: lhs
  type(var), intent(in) :: this
  logical, intent(in), optional :: dealloc
  logical, intent(out), optional :: success
  logical :: ldealloc, lsuccess
  type :: pt
    type(var_str), pointer :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "V0"
  if (present(success)) success = lsuccess
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  ! there is one problem, say if lhs is not nullified...
  if (ldealloc.and.associated(lhs)) then
     deallocate(lhs)
     nullify(lhs)
  end if
  if (.not. lsuccess ) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs => p%p
end subroutine associate_get_V0
subroutine associate_set_V0(this,rhs,dealloc)
  type(var), intent(inout) :: this
  type(var_str), intent(in), target :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    type(var_str), pointer :: p => null()
  end type
  type(pt) :: p
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  this%t = "V0"
  p%p => rhs
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
end subroutine associate_set_V0
pure function associatd_l_V0(lhs,this) result(ret)
  type(var_str), pointer :: lhs
  type(var), intent(in) :: this
  logical :: ret
  type :: pt
    type(var_str), pointer :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "V0"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(lhs,p%p)
  endif
end function associatd_l_V0
pure function associatd_r_V0(this,rhs) result(ret)
  type(var), intent(in) :: this
  type(var_str), pointer :: rhs
  logical :: ret
  type :: pt
    type(var_str), pointer :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "V0"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(p%p,rhs)
  endif
end function associatd_r_V0
! All boolean functions
subroutine assign_set_s0(this,rhs,dealloc)
  type(var), intent(inout) :: this
  real(sp), intent(in) :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    real(sp), pointer :: p => null()
  end type
  type(pt) :: p
  ! ASSIGNMENT in fortran is per default destructive
  ldealloc = .true.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  ! With pointer transfer we need to deallocate
  ! else bounds might change...
  this%t = "s0"
  allocate(p%p) ! allocate space
  p%p = rhs ! copy data over
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
  ! We already have shipped it
  nullify(p%p)
end subroutine assign_set_s0
subroutine assign_get_s0(lhs,this,success)
  real(sp), intent(out) :: lhs
  type(var), intent(in) :: this
  logical, intent(out), optional :: success
  logical :: lsuccess
  type :: pt
    real(sp), pointer :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "s0"
  if (present(success)) success = lsuccess
  if (.not. lsuccess) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs = p%p
end subroutine assign_get_s0
subroutine associate_get_s0(lhs,this,dealloc,success)
  real(sp), pointer :: lhs
  type(var), intent(in) :: this
  logical, intent(in), optional :: dealloc
  logical, intent(out), optional :: success
  logical :: ldealloc, lsuccess
  type :: pt
    real(sp), pointer :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "s0"
  if (present(success)) success = lsuccess
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  ! there is one problem, say if lhs is not nullified...
  if (ldealloc.and.associated(lhs)) then
     deallocate(lhs)
     nullify(lhs)
  end if
  if (.not. lsuccess ) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs => p%p
end subroutine associate_get_s0
subroutine associate_set_s0(this,rhs,dealloc)
  type(var), intent(inout) :: this
  real(sp), intent(in), target :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    real(sp), pointer :: p => null()
  end type
  type(pt) :: p
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  this%t = "s0"
  p%p => rhs
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
end subroutine associate_set_s0
pure function associatd_l_s0(lhs,this) result(ret)
  real(sp), pointer :: lhs
  type(var), intent(in) :: this
  logical :: ret
  type :: pt
    real(sp), pointer :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "s0"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(lhs,p%p)
  endif
end function associatd_l_s0
pure function associatd_r_s0(this,rhs) result(ret)
  type(var), intent(in) :: this
  real(sp), pointer :: rhs
  logical :: ret
  type :: pt
    real(sp), pointer :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "s0"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(p%p,rhs)
  endif
end function associatd_r_s0
! All boolean functions
subroutine assign_set_s1(this,rhs,dealloc)
  type(var), intent(inout) :: this
  real(sp), intent(in), dimension(:) :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    real(sp), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  ! ASSIGNMENT in fortran is per default destructive
  ldealloc = .true.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  ! With pointer transfer we need to deallocate
  ! else bounds might change...
  this%t = "s1"
  allocate(p%p(size(rhs))) ! allocate space
  p%p = rhs ! copy data over
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
  ! We already have shipped it
  nullify(p%p)
end subroutine assign_set_s1
subroutine assign_get_s1(lhs,this,success)
  real(sp), intent(out), dimension(:) :: lhs
  type(var), intent(in) :: this
  logical, intent(out), optional :: success
  logical :: lsuccess
  type :: pt
    real(sp), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "s1"
  if (lsuccess) then
    p = transfer(this%enc,p) ! retrieve pointer encoding
    lsuccess = all(shape(p%p)==shape(lhs)) !&
     ! .and. all((lbound(p%p) == lbound(lhs))) &
     ! .and. all((ubound(p%p) == ubound(lhs)))
  end if
  if (present(success)) success = lsuccess
  if (.not. lsuccess) return
  lhs = p%p
end subroutine assign_get_s1
subroutine associate_get_s1(lhs,this,dealloc,success)
  real(sp), pointer , dimension(:) :: lhs
  type(var), intent(in) :: this
  logical, intent(in), optional :: dealloc
  logical, intent(out), optional :: success
  logical :: ldealloc, lsuccess
  type :: pt
    real(sp), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "s1"
  if (present(success)) success = lsuccess
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  ! there is one problem, say if lhs is not nullified...
  if (ldealloc.and.associated(lhs)) then
     deallocate(lhs)
     nullify(lhs)
  end if
  if (.not. lsuccess ) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs => p%p
end subroutine associate_get_s1
subroutine associate_set_s1(this,rhs,dealloc)
  type(var), intent(inout) :: this
  real(sp), intent(in), dimension(:), target :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    real(sp), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  this%t = "s1"
  p%p => rhs
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
end subroutine associate_set_s1
pure function associatd_l_s1(lhs,this) result(ret)
  real(sp), pointer , dimension(:) :: lhs
  type(var), intent(in) :: this
  logical :: ret
  type :: pt
    real(sp), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "s1"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(lhs,p%p)
  endif
end function associatd_l_s1
pure function associatd_r_s1(this,rhs) result(ret)
  type(var), intent(in) :: this
  real(sp), pointer , dimension(:) :: rhs
  logical :: ret
  type :: pt
    real(sp), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "s1"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(p%p,rhs)
  endif
end function associatd_r_s1
! All boolean functions
subroutine assign_set_s2(this,rhs,dealloc)
  type(var), intent(inout) :: this
  real(sp), intent(in), dimension(:,:) :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    real(sp), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  ! ASSIGNMENT in fortran is per default destructive
  ldealloc = .true.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  ! With pointer transfer we need to deallocate
  ! else bounds might change...
  this%t = "s2"
  allocate(p%p(size(rhs,1),size(rhs,2))) ! allocate space
  p%p = rhs ! copy data over
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
  ! We already have shipped it
  nullify(p%p)
end subroutine assign_set_s2
subroutine assign_get_s2(lhs,this,success)
  real(sp), intent(out), dimension(:,:) :: lhs
  type(var), intent(in) :: this
  logical, intent(out), optional :: success
  logical :: lsuccess
  type :: pt
    real(sp), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "s2"
  if (lsuccess) then
    p = transfer(this%enc,p) ! retrieve pointer encoding
    lsuccess = all(shape(p%p)==shape(lhs)) !&
     ! .and. all((lbound(p%p) == lbound(lhs))) &
     ! .and. all((ubound(p%p) == ubound(lhs)))
  end if
  if (present(success)) success = lsuccess
  if (.not. lsuccess) return
  lhs = p%p
end subroutine assign_get_s2
subroutine associate_get_s2(lhs,this,dealloc,success)
  real(sp), pointer , dimension(:,:) :: lhs
  type(var), intent(in) :: this
  logical, intent(in), optional :: dealloc
  logical, intent(out), optional :: success
  logical :: ldealloc, lsuccess
  type :: pt
    real(sp), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "s2"
  if (present(success)) success = lsuccess
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  ! there is one problem, say if lhs is not nullified...
  if (ldealloc.and.associated(lhs)) then
     deallocate(lhs)
     nullify(lhs)
  end if
  if (.not. lsuccess ) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs => p%p
end subroutine associate_get_s2
subroutine associate_set_s2(this,rhs,dealloc)
  type(var), intent(inout) :: this
  real(sp), intent(in), dimension(:,:), target :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    real(sp), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  this%t = "s2"
  p%p => rhs
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
end subroutine associate_set_s2
pure function associatd_l_s2(lhs,this) result(ret)
  real(sp), pointer , dimension(:,:) :: lhs
  type(var), intent(in) :: this
  logical :: ret
  type :: pt
    real(sp), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "s2"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(lhs,p%p)
  endif
end function associatd_l_s2
pure function associatd_r_s2(this,rhs) result(ret)
  type(var), intent(in) :: this
  real(sp), pointer , dimension(:,:) :: rhs
  logical :: ret
  type :: pt
    real(sp), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "s2"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(p%p,rhs)
  endif
end function associatd_r_s2
! All boolean functions
subroutine assign_set_d0(this,rhs,dealloc)
  type(var), intent(inout) :: this
  real(dp), intent(in) :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    real(dp), pointer :: p => null()
  end type
  type(pt) :: p
  ! ASSIGNMENT in fortran is per default destructive
  ldealloc = .true.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  ! With pointer transfer we need to deallocate
  ! else bounds might change...
  this%t = "d0"
  allocate(p%p) ! allocate space
  p%p = rhs ! copy data over
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
  ! We already have shipped it
  nullify(p%p)
end subroutine assign_set_d0
subroutine assign_get_d0(lhs,this,success)
  real(dp), intent(out) :: lhs
  type(var), intent(in) :: this
  logical, intent(out), optional :: success
  logical :: lsuccess
  type :: pt
    real(dp), pointer :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "d0"
  if (present(success)) success = lsuccess
  if (.not. lsuccess) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs = p%p
end subroutine assign_get_d0
subroutine associate_get_d0(lhs,this,dealloc,success)
  real(dp), pointer :: lhs
  type(var), intent(in) :: this
  logical, intent(in), optional :: dealloc
  logical, intent(out), optional :: success
  logical :: ldealloc, lsuccess
  type :: pt
    real(dp), pointer :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "d0"
  if (present(success)) success = lsuccess
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  ! there is one problem, say if lhs is not nullified...
  if (ldealloc.and.associated(lhs)) then
     deallocate(lhs)
     nullify(lhs)
  end if
  if (.not. lsuccess ) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs => p%p
end subroutine associate_get_d0
subroutine associate_set_d0(this,rhs,dealloc)
  type(var), intent(inout) :: this
  real(dp), intent(in), target :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    real(dp), pointer :: p => null()
  end type
  type(pt) :: p
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  this%t = "d0"
  p%p => rhs
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
end subroutine associate_set_d0
pure function associatd_l_d0(lhs,this) result(ret)
  real(dp), pointer :: lhs
  type(var), intent(in) :: this
  logical :: ret
  type :: pt
    real(dp), pointer :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "d0"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(lhs,p%p)
  endif
end function associatd_l_d0
pure function associatd_r_d0(this,rhs) result(ret)
  type(var), intent(in) :: this
  real(dp), pointer :: rhs
  logical :: ret
  type :: pt
    real(dp), pointer :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "d0"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(p%p,rhs)
  endif
end function associatd_r_d0
! All boolean functions
subroutine assign_set_d1(this,rhs,dealloc)
  type(var), intent(inout) :: this
  real(dp), intent(in), dimension(:) :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    real(dp), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  ! ASSIGNMENT in fortran is per default destructive
  ldealloc = .true.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  ! With pointer transfer we need to deallocate
  ! else bounds might change...
  this%t = "d1"
  allocate(p%p(size(rhs))) ! allocate space
  p%p = rhs ! copy data over
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
  ! We already have shipped it
  nullify(p%p)
end subroutine assign_set_d1
subroutine assign_get_d1(lhs,this,success)
  real(dp), intent(out), dimension(:) :: lhs
  type(var), intent(in) :: this
  logical, intent(out), optional :: success
  logical :: lsuccess
  type :: pt
    real(dp), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "d1"
  if (lsuccess) then
    p = transfer(this%enc,p) ! retrieve pointer encoding
    lsuccess = all(shape(p%p)==shape(lhs)) !&
     ! .and. all((lbound(p%p) == lbound(lhs))) &
     ! .and. all((ubound(p%p) == ubound(lhs)))
  end if
  if (present(success)) success = lsuccess
  if (.not. lsuccess) return
  lhs = p%p
end subroutine assign_get_d1
subroutine associate_get_d1(lhs,this,dealloc,success)
  real(dp), pointer , dimension(:) :: lhs
  type(var), intent(in) :: this
  logical, intent(in), optional :: dealloc
  logical, intent(out), optional :: success
  logical :: ldealloc, lsuccess
  type :: pt
    real(dp), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "d1"
  if (present(success)) success = lsuccess
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  ! there is one problem, say if lhs is not nullified...
  if (ldealloc.and.associated(lhs)) then
     deallocate(lhs)
     nullify(lhs)
  end if
  if (.not. lsuccess ) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs => p%p
end subroutine associate_get_d1
subroutine associate_set_d1(this,rhs,dealloc)
  type(var), intent(inout) :: this
  real(dp), intent(in), dimension(:), target :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    real(dp), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  this%t = "d1"
  p%p => rhs
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
end subroutine associate_set_d1
pure function associatd_l_d1(lhs,this) result(ret)
  real(dp), pointer , dimension(:) :: lhs
  type(var), intent(in) :: this
  logical :: ret
  type :: pt
    real(dp), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "d1"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(lhs,p%p)
  endif
end function associatd_l_d1
pure function associatd_r_d1(this,rhs) result(ret)
  type(var), intent(in) :: this
  real(dp), pointer , dimension(:) :: rhs
  logical :: ret
  type :: pt
    real(dp), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "d1"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(p%p,rhs)
  endif
end function associatd_r_d1
! All boolean functions
subroutine assign_set_d2(this,rhs,dealloc)
  type(var), intent(inout) :: this
  real(dp), intent(in), dimension(:,:) :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    real(dp), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  ! ASSIGNMENT in fortran is per default destructive
  ldealloc = .true.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  ! With pointer transfer we need to deallocate
  ! else bounds might change...
  this%t = "d2"
  allocate(p%p(size(rhs,1),size(rhs,2))) ! allocate space
  p%p = rhs ! copy data over
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
  ! We already have shipped it
  nullify(p%p)
end subroutine assign_set_d2
subroutine assign_get_d2(lhs,this,success)
  real(dp), intent(out), dimension(:,:) :: lhs
  type(var), intent(in) :: this
  logical, intent(out), optional :: success
  logical :: lsuccess
  type :: pt
    real(dp), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "d2"
  if (lsuccess) then
    p = transfer(this%enc,p) ! retrieve pointer encoding
    lsuccess = all(shape(p%p)==shape(lhs)) !&
     ! .and. all((lbound(p%p) == lbound(lhs))) &
     ! .and. all((ubound(p%p) == ubound(lhs)))
  end if
  if (present(success)) success = lsuccess
  if (.not. lsuccess) return
  lhs = p%p
end subroutine assign_get_d2
subroutine associate_get_d2(lhs,this,dealloc,success)
  real(dp), pointer , dimension(:,:) :: lhs
  type(var), intent(in) :: this
  logical, intent(in), optional :: dealloc
  logical, intent(out), optional :: success
  logical :: ldealloc, lsuccess
  type :: pt
    real(dp), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "d2"
  if (present(success)) success = lsuccess
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  ! there is one problem, say if lhs is not nullified...
  if (ldealloc.and.associated(lhs)) then
     deallocate(lhs)
     nullify(lhs)
  end if
  if (.not. lsuccess ) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs => p%p
end subroutine associate_get_d2
subroutine associate_set_d2(this,rhs,dealloc)
  type(var), intent(inout) :: this
  real(dp), intent(in), dimension(:,:), target :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    real(dp), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  this%t = "d2"
  p%p => rhs
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
end subroutine associate_set_d2
pure function associatd_l_d2(lhs,this) result(ret)
  real(dp), pointer , dimension(:,:) :: lhs
  type(var), intent(in) :: this
  logical :: ret
  type :: pt
    real(dp), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "d2"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(lhs,p%p)
  endif
end function associatd_l_d2
pure function associatd_r_d2(this,rhs) result(ret)
  type(var), intent(in) :: this
  real(dp), pointer , dimension(:,:) :: rhs
  logical :: ret
  type :: pt
    real(dp), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "d2"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(p%p,rhs)
  endif
end function associatd_r_d2
! All boolean functions
subroutine assign_set_c0(this,rhs,dealloc)
  type(var), intent(inout) :: this
  complex(sp), intent(in) :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    complex(sp), pointer :: p => null()
  end type
  type(pt) :: p
  ! ASSIGNMENT in fortran is per default destructive
  ldealloc = .true.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  ! With pointer transfer we need to deallocate
  ! else bounds might change...
  this%t = "c0"
  allocate(p%p) ! allocate space
  p%p = rhs ! copy data over
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
  ! We already have shipped it
  nullify(p%p)
end subroutine assign_set_c0
subroutine assign_get_c0(lhs,this,success)
  complex(sp), intent(out) :: lhs
  type(var), intent(in) :: this
  logical, intent(out), optional :: success
  logical :: lsuccess
  type :: pt
    complex(sp), pointer :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "c0"
  if (present(success)) success = lsuccess
  if (.not. lsuccess) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs = p%p
end subroutine assign_get_c0
subroutine associate_get_c0(lhs,this,dealloc,success)
  complex(sp), pointer :: lhs
  type(var), intent(in) :: this
  logical, intent(in), optional :: dealloc
  logical, intent(out), optional :: success
  logical :: ldealloc, lsuccess
  type :: pt
    complex(sp), pointer :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "c0"
  if (present(success)) success = lsuccess
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  ! there is one problem, say if lhs is not nullified...
  if (ldealloc.and.associated(lhs)) then
     deallocate(lhs)
     nullify(lhs)
  end if
  if (.not. lsuccess ) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs => p%p
end subroutine associate_get_c0
subroutine associate_set_c0(this,rhs,dealloc)
  type(var), intent(inout) :: this
  complex(sp), intent(in), target :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    complex(sp), pointer :: p => null()
  end type
  type(pt) :: p
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  this%t = "c0"
  p%p => rhs
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
end subroutine associate_set_c0
pure function associatd_l_c0(lhs,this) result(ret)
  complex(sp), pointer :: lhs
  type(var), intent(in) :: this
  logical :: ret
  type :: pt
    complex(sp), pointer :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "c0"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(lhs,p%p)
  endif
end function associatd_l_c0
pure function associatd_r_c0(this,rhs) result(ret)
  type(var), intent(in) :: this
  complex(sp), pointer :: rhs
  logical :: ret
  type :: pt
    complex(sp), pointer :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "c0"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(p%p,rhs)
  endif
end function associatd_r_c0
! All boolean functions
subroutine assign_set_c1(this,rhs,dealloc)
  type(var), intent(inout) :: this
  complex(sp), intent(in), dimension(:) :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    complex(sp), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  ! ASSIGNMENT in fortran is per default destructive
  ldealloc = .true.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  ! With pointer transfer we need to deallocate
  ! else bounds might change...
  this%t = "c1"
  allocate(p%p(size(rhs))) ! allocate space
  p%p = rhs ! copy data over
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
  ! We already have shipped it
  nullify(p%p)
end subroutine assign_set_c1
subroutine assign_get_c1(lhs,this,success)
  complex(sp), intent(out), dimension(:) :: lhs
  type(var), intent(in) :: this
  logical, intent(out), optional :: success
  logical :: lsuccess
  type :: pt
    complex(sp), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "c1"
  if (lsuccess) then
    p = transfer(this%enc,p) ! retrieve pointer encoding
    lsuccess = all(shape(p%p)==shape(lhs)) !&
     ! .and. all((lbound(p%p) == lbound(lhs))) &
     ! .and. all((ubound(p%p) == ubound(lhs)))
  end if
  if (present(success)) success = lsuccess
  if (.not. lsuccess) return
  lhs = p%p
end subroutine assign_get_c1
subroutine associate_get_c1(lhs,this,dealloc,success)
  complex(sp), pointer , dimension(:) :: lhs
  type(var), intent(in) :: this
  logical, intent(in), optional :: dealloc
  logical, intent(out), optional :: success
  logical :: ldealloc, lsuccess
  type :: pt
    complex(sp), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "c1"
  if (present(success)) success = lsuccess
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  ! there is one problem, say if lhs is not nullified...
  if (ldealloc.and.associated(lhs)) then
     deallocate(lhs)
     nullify(lhs)
  end if
  if (.not. lsuccess ) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs => p%p
end subroutine associate_get_c1
subroutine associate_set_c1(this,rhs,dealloc)
  type(var), intent(inout) :: this
  complex(sp), intent(in), dimension(:), target :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    complex(sp), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  this%t = "c1"
  p%p => rhs
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
end subroutine associate_set_c1
pure function associatd_l_c1(lhs,this) result(ret)
  complex(sp), pointer , dimension(:) :: lhs
  type(var), intent(in) :: this
  logical :: ret
  type :: pt
    complex(sp), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "c1"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(lhs,p%p)
  endif
end function associatd_l_c1
pure function associatd_r_c1(this,rhs) result(ret)
  type(var), intent(in) :: this
  complex(sp), pointer , dimension(:) :: rhs
  logical :: ret
  type :: pt
    complex(sp), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "c1"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(p%p,rhs)
  endif
end function associatd_r_c1
! All boolean functions
subroutine assign_set_c2(this,rhs,dealloc)
  type(var), intent(inout) :: this
  complex(sp), intent(in), dimension(:,:) :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    complex(sp), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  ! ASSIGNMENT in fortran is per default destructive
  ldealloc = .true.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  ! With pointer transfer we need to deallocate
  ! else bounds might change...
  this%t = "c2"
  allocate(p%p(size(rhs,1),size(rhs,2))) ! allocate space
  p%p = rhs ! copy data over
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
  ! We already have shipped it
  nullify(p%p)
end subroutine assign_set_c2
subroutine assign_get_c2(lhs,this,success)
  complex(sp), intent(out), dimension(:,:) :: lhs
  type(var), intent(in) :: this
  logical, intent(out), optional :: success
  logical :: lsuccess
  type :: pt
    complex(sp), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "c2"
  if (lsuccess) then
    p = transfer(this%enc,p) ! retrieve pointer encoding
    lsuccess = all(shape(p%p)==shape(lhs)) !&
     ! .and. all((lbound(p%p) == lbound(lhs))) &
     ! .and. all((ubound(p%p) == ubound(lhs)))
  end if
  if (present(success)) success = lsuccess
  if (.not. lsuccess) return
  lhs = p%p
end subroutine assign_get_c2
subroutine associate_get_c2(lhs,this,dealloc,success)
  complex(sp), pointer , dimension(:,:) :: lhs
  type(var), intent(in) :: this
  logical, intent(in), optional :: dealloc
  logical, intent(out), optional :: success
  logical :: ldealloc, lsuccess
  type :: pt
    complex(sp), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "c2"
  if (present(success)) success = lsuccess
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  ! there is one problem, say if lhs is not nullified...
  if (ldealloc.and.associated(lhs)) then
     deallocate(lhs)
     nullify(lhs)
  end if
  if (.not. lsuccess ) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs => p%p
end subroutine associate_get_c2
subroutine associate_set_c2(this,rhs,dealloc)
  type(var), intent(inout) :: this
  complex(sp), intent(in), dimension(:,:), target :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    complex(sp), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  this%t = "c2"
  p%p => rhs
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
end subroutine associate_set_c2
pure function associatd_l_c2(lhs,this) result(ret)
  complex(sp), pointer , dimension(:,:) :: lhs
  type(var), intent(in) :: this
  logical :: ret
  type :: pt
    complex(sp), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "c2"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(lhs,p%p)
  endif
end function associatd_l_c2
pure function associatd_r_c2(this,rhs) result(ret)
  type(var), intent(in) :: this
  complex(sp), pointer , dimension(:,:) :: rhs
  logical :: ret
  type :: pt
    complex(sp), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "c2"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(p%p,rhs)
  endif
end function associatd_r_c2
! All boolean functions
subroutine assign_set_z0(this,rhs,dealloc)
  type(var), intent(inout) :: this
  complex(dp), intent(in) :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    complex(dp), pointer :: p => null()
  end type
  type(pt) :: p
  ! ASSIGNMENT in fortran is per default destructive
  ldealloc = .true.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  ! With pointer transfer we need to deallocate
  ! else bounds might change...
  this%t = "z0"
  allocate(p%p) ! allocate space
  p%p = rhs ! copy data over
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
  ! We already have shipped it
  nullify(p%p)
end subroutine assign_set_z0
subroutine assign_get_z0(lhs,this,success)
  complex(dp), intent(out) :: lhs
  type(var), intent(in) :: this
  logical, intent(out), optional :: success
  logical :: lsuccess
  type :: pt
    complex(dp), pointer :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "z0"
  if (present(success)) success = lsuccess
  if (.not. lsuccess) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs = p%p
end subroutine assign_get_z0
subroutine associate_get_z0(lhs,this,dealloc,success)
  complex(dp), pointer :: lhs
  type(var), intent(in) :: this
  logical, intent(in), optional :: dealloc
  logical, intent(out), optional :: success
  logical :: ldealloc, lsuccess
  type :: pt
    complex(dp), pointer :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "z0"
  if (present(success)) success = lsuccess
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  ! there is one problem, say if lhs is not nullified...
  if (ldealloc.and.associated(lhs)) then
     deallocate(lhs)
     nullify(lhs)
  end if
  if (.not. lsuccess ) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs => p%p
end subroutine associate_get_z0
subroutine associate_set_z0(this,rhs,dealloc)
  type(var), intent(inout) :: this
  complex(dp), intent(in), target :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    complex(dp), pointer :: p => null()
  end type
  type(pt) :: p
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  this%t = "z0"
  p%p => rhs
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
end subroutine associate_set_z0
pure function associatd_l_z0(lhs,this) result(ret)
  complex(dp), pointer :: lhs
  type(var), intent(in) :: this
  logical :: ret
  type :: pt
    complex(dp), pointer :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "z0"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(lhs,p%p)
  endif
end function associatd_l_z0
pure function associatd_r_z0(this,rhs) result(ret)
  type(var), intent(in) :: this
  complex(dp), pointer :: rhs
  logical :: ret
  type :: pt
    complex(dp), pointer :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "z0"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(p%p,rhs)
  endif
end function associatd_r_z0
! All boolean functions
subroutine assign_set_z1(this,rhs,dealloc)
  type(var), intent(inout) :: this
  complex(dp), intent(in), dimension(:) :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    complex(dp), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  ! ASSIGNMENT in fortran is per default destructive
  ldealloc = .true.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  ! With pointer transfer we need to deallocate
  ! else bounds might change...
  this%t = "z1"
  allocate(p%p(size(rhs))) ! allocate space
  p%p = rhs ! copy data over
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
  ! We already have shipped it
  nullify(p%p)
end subroutine assign_set_z1
subroutine assign_get_z1(lhs,this,success)
  complex(dp), intent(out), dimension(:) :: lhs
  type(var), intent(in) :: this
  logical, intent(out), optional :: success
  logical :: lsuccess
  type :: pt
    complex(dp), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "z1"
  if (lsuccess) then
    p = transfer(this%enc,p) ! retrieve pointer encoding
    lsuccess = all(shape(p%p)==shape(lhs)) !&
     ! .and. all((lbound(p%p) == lbound(lhs))) &
     ! .and. all((ubound(p%p) == ubound(lhs)))
  end if
  if (present(success)) success = lsuccess
  if (.not. lsuccess) return
  lhs = p%p
end subroutine assign_get_z1
subroutine associate_get_z1(lhs,this,dealloc,success)
  complex(dp), pointer , dimension(:) :: lhs
  type(var), intent(in) :: this
  logical, intent(in), optional :: dealloc
  logical, intent(out), optional :: success
  logical :: ldealloc, lsuccess
  type :: pt
    complex(dp), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "z1"
  if (present(success)) success = lsuccess
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  ! there is one problem, say if lhs is not nullified...
  if (ldealloc.and.associated(lhs)) then
     deallocate(lhs)
     nullify(lhs)
  end if
  if (.not. lsuccess ) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs => p%p
end subroutine associate_get_z1
subroutine associate_set_z1(this,rhs,dealloc)
  type(var), intent(inout) :: this
  complex(dp), intent(in), dimension(:), target :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    complex(dp), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  this%t = "z1"
  p%p => rhs
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
end subroutine associate_set_z1
pure function associatd_l_z1(lhs,this) result(ret)
  complex(dp), pointer , dimension(:) :: lhs
  type(var), intent(in) :: this
  logical :: ret
  type :: pt
    complex(dp), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "z1"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(lhs,p%p)
  endif
end function associatd_l_z1
pure function associatd_r_z1(this,rhs) result(ret)
  type(var), intent(in) :: this
  complex(dp), pointer , dimension(:) :: rhs
  logical :: ret
  type :: pt
    complex(dp), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "z1"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(p%p,rhs)
  endif
end function associatd_r_z1
! All boolean functions
subroutine assign_set_z2(this,rhs,dealloc)
  type(var), intent(inout) :: this
  complex(dp), intent(in), dimension(:,:) :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    complex(dp), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  ! ASSIGNMENT in fortran is per default destructive
  ldealloc = .true.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  ! With pointer transfer we need to deallocate
  ! else bounds might change...
  this%t = "z2"
  allocate(p%p(size(rhs,1),size(rhs,2))) ! allocate space
  p%p = rhs ! copy data over
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
  ! We already have shipped it
  nullify(p%p)
end subroutine assign_set_z2
subroutine assign_get_z2(lhs,this,success)
  complex(dp), intent(out), dimension(:,:) :: lhs
  type(var), intent(in) :: this
  logical, intent(out), optional :: success
  logical :: lsuccess
  type :: pt
    complex(dp), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "z2"
  if (lsuccess) then
    p = transfer(this%enc,p) ! retrieve pointer encoding
    lsuccess = all(shape(p%p)==shape(lhs)) !&
     ! .and. all((lbound(p%p) == lbound(lhs))) &
     ! .and. all((ubound(p%p) == ubound(lhs)))
  end if
  if (present(success)) success = lsuccess
  if (.not. lsuccess) return
  lhs = p%p
end subroutine assign_get_z2
subroutine associate_get_z2(lhs,this,dealloc,success)
  complex(dp), pointer , dimension(:,:) :: lhs
  type(var), intent(in) :: this
  logical, intent(in), optional :: dealloc
  logical, intent(out), optional :: success
  logical :: ldealloc, lsuccess
  type :: pt
    complex(dp), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "z2"
  if (present(success)) success = lsuccess
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  ! there is one problem, say if lhs is not nullified...
  if (ldealloc.and.associated(lhs)) then
     deallocate(lhs)
     nullify(lhs)
  end if
  if (.not. lsuccess ) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs => p%p
end subroutine associate_get_z2
subroutine associate_set_z2(this,rhs,dealloc)
  type(var), intent(inout) :: this
  complex(dp), intent(in), dimension(:,:), target :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    complex(dp), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  this%t = "z2"
  p%p => rhs
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
end subroutine associate_set_z2
pure function associatd_l_z2(lhs,this) result(ret)
  complex(dp), pointer , dimension(:,:) :: lhs
  type(var), intent(in) :: this
  logical :: ret
  type :: pt
    complex(dp), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "z2"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(lhs,p%p)
  endif
end function associatd_l_z2
pure function associatd_r_z2(this,rhs) result(ret)
  type(var), intent(in) :: this
  complex(dp), pointer , dimension(:,:) :: rhs
  logical :: ret
  type :: pt
    complex(dp), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "z2"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(p%p,rhs)
  endif
end function associatd_r_z2
! All boolean functions
subroutine assign_set_b0(this,rhs,dealloc)
  type(var), intent(inout) :: this
  logical, intent(in) :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    logical, pointer :: p => null()
  end type
  type(pt) :: p
  ! ASSIGNMENT in fortran is per default destructive
  ldealloc = .true.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  ! With pointer transfer we need to deallocate
  ! else bounds might change...
  this%t = "b0"
  allocate(p%p) ! allocate space
  p%p = rhs ! copy data over
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
  ! We already have shipped it
  nullify(p%p)
end subroutine assign_set_b0
subroutine assign_get_b0(lhs,this,success)
  logical, intent(out) :: lhs
  type(var), intent(in) :: this
  logical, intent(out), optional :: success
  logical :: lsuccess
  type :: pt
    logical, pointer :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "b0"
  if (present(success)) success = lsuccess
  if (.not. lsuccess) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs = p%p
end subroutine assign_get_b0
subroutine associate_get_b0(lhs,this,dealloc,success)
  logical, pointer :: lhs
  type(var), intent(in) :: this
  logical, intent(in), optional :: dealloc
  logical, intent(out), optional :: success
  logical :: ldealloc, lsuccess
  type :: pt
    logical, pointer :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "b0"
  if (present(success)) success = lsuccess
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  ! there is one problem, say if lhs is not nullified...
  if (ldealloc.and.associated(lhs)) then
     deallocate(lhs)
     nullify(lhs)
  end if
  if (.not. lsuccess ) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs => p%p
end subroutine associate_get_b0
subroutine associate_set_b0(this,rhs,dealloc)
  type(var), intent(inout) :: this
  logical, intent(in), target :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    logical, pointer :: p => null()
  end type
  type(pt) :: p
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  this%t = "b0"
  p%p => rhs
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
end subroutine associate_set_b0
pure function associatd_l_b0(lhs,this) result(ret)
  logical, pointer :: lhs
  type(var), intent(in) :: this
  logical :: ret
  type :: pt
    logical, pointer :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "b0"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(lhs,p%p)
  endif
end function associatd_l_b0
pure function associatd_r_b0(this,rhs) result(ret)
  type(var), intent(in) :: this
  logical, pointer :: rhs
  logical :: ret
  type :: pt
    logical, pointer :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "b0"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(p%p,rhs)
  endif
end function associatd_r_b0
! All boolean functions
subroutine assign_set_b1(this,rhs,dealloc)
  type(var), intent(inout) :: this
  logical, intent(in), dimension(:) :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    logical, pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  ! ASSIGNMENT in fortran is per default destructive
  ldealloc = .true.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  ! With pointer transfer we need to deallocate
  ! else bounds might change...
  this%t = "b1"
  allocate(p%p(size(rhs))) ! allocate space
  p%p = rhs ! copy data over
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
  ! We already have shipped it
  nullify(p%p)
end subroutine assign_set_b1
subroutine assign_get_b1(lhs,this,success)
  logical, intent(out), dimension(:) :: lhs
  type(var), intent(in) :: this
  logical, intent(out), optional :: success
  logical :: lsuccess
  type :: pt
    logical, pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "b1"
  if (lsuccess) then
    p = transfer(this%enc,p) ! retrieve pointer encoding
    lsuccess = all(shape(p%p)==shape(lhs)) !&
     ! .and. all((lbound(p%p) == lbound(lhs))) &
     ! .and. all((ubound(p%p) == ubound(lhs)))
  end if
  if (present(success)) success = lsuccess
  if (.not. lsuccess) return
  lhs = p%p
end subroutine assign_get_b1
subroutine associate_get_b1(lhs,this,dealloc,success)
  logical, pointer , dimension(:) :: lhs
  type(var), intent(in) :: this
  logical, intent(in), optional :: dealloc
  logical, intent(out), optional :: success
  logical :: ldealloc, lsuccess
  type :: pt
    logical, pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "b1"
  if (present(success)) success = lsuccess
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  ! there is one problem, say if lhs is not nullified...
  if (ldealloc.and.associated(lhs)) then
     deallocate(lhs)
     nullify(lhs)
  end if
  if (.not. lsuccess ) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs => p%p
end subroutine associate_get_b1
subroutine associate_set_b1(this,rhs,dealloc)
  type(var), intent(inout) :: this
  logical, intent(in), dimension(:), target :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    logical, pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  this%t = "b1"
  p%p => rhs
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
end subroutine associate_set_b1
pure function associatd_l_b1(lhs,this) result(ret)
  logical, pointer , dimension(:) :: lhs
  type(var), intent(in) :: this
  logical :: ret
  type :: pt
    logical, pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "b1"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(lhs,p%p)
  endif
end function associatd_l_b1
pure function associatd_r_b1(this,rhs) result(ret)
  type(var), intent(in) :: this
  logical, pointer , dimension(:) :: rhs
  logical :: ret
  type :: pt
    logical, pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "b1"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(p%p,rhs)
  endif
end function associatd_r_b1
! All boolean functions
subroutine assign_set_b2(this,rhs,dealloc)
  type(var), intent(inout) :: this
  logical, intent(in), dimension(:,:) :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    logical, pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  ! ASSIGNMENT in fortran is per default destructive
  ldealloc = .true.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  ! With pointer transfer we need to deallocate
  ! else bounds might change...
  this%t = "b2"
  allocate(p%p(size(rhs,1),size(rhs,2))) ! allocate space
  p%p = rhs ! copy data over
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
  ! We already have shipped it
  nullify(p%p)
end subroutine assign_set_b2
subroutine assign_get_b2(lhs,this,success)
  logical, intent(out), dimension(:,:) :: lhs
  type(var), intent(in) :: this
  logical, intent(out), optional :: success
  logical :: lsuccess
  type :: pt
    logical, pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "b2"
  if (lsuccess) then
    p = transfer(this%enc,p) ! retrieve pointer encoding
    lsuccess = all(shape(p%p)==shape(lhs)) !&
     ! .and. all((lbound(p%p) == lbound(lhs))) &
     ! .and. all((ubound(p%p) == ubound(lhs)))
  end if
  if (present(success)) success = lsuccess
  if (.not. lsuccess) return
  lhs = p%p
end subroutine assign_get_b2
subroutine associate_get_b2(lhs,this,dealloc,success)
  logical, pointer , dimension(:,:) :: lhs
  type(var), intent(in) :: this
  logical, intent(in), optional :: dealloc
  logical, intent(out), optional :: success
  logical :: ldealloc, lsuccess
  type :: pt
    logical, pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "b2"
  if (present(success)) success = lsuccess
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  ! there is one problem, say if lhs is not nullified...
  if (ldealloc.and.associated(lhs)) then
     deallocate(lhs)
     nullify(lhs)
  end if
  if (.not. lsuccess ) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs => p%p
end subroutine associate_get_b2
subroutine associate_set_b2(this,rhs,dealloc)
  type(var), intent(inout) :: this
  logical, intent(in), dimension(:,:), target :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    logical, pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  this%t = "b2"
  p%p => rhs
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
end subroutine associate_set_b2
pure function associatd_l_b2(lhs,this) result(ret)
  logical, pointer , dimension(:,:) :: lhs
  type(var), intent(in) :: this
  logical :: ret
  type :: pt
    logical, pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "b2"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(lhs,p%p)
  endif
end function associatd_l_b2
pure function associatd_r_b2(this,rhs) result(ret)
  type(var), intent(in) :: this
  logical, pointer , dimension(:,:) :: rhs
  logical :: ret
  type :: pt
    logical, pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "b2"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(p%p,rhs)
  endif
end function associatd_r_b2
! All boolean functions
subroutine assign_set_h0(this,rhs,dealloc)
  type(var), intent(inout) :: this
  integer(ih), intent(in) :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    integer(ih), pointer :: p => null()
  end type
  type(pt) :: p
  ! ASSIGNMENT in fortran is per default destructive
  ldealloc = .true.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  ! With pointer transfer we need to deallocate
  ! else bounds might change...
  this%t = "h0"
  allocate(p%p) ! allocate space
  p%p = rhs ! copy data over
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
  ! We already have shipped it
  nullify(p%p)
end subroutine assign_set_h0
subroutine assign_get_h0(lhs,this,success)
  integer(ih), intent(out) :: lhs
  type(var), intent(in) :: this
  logical, intent(out), optional :: success
  logical :: lsuccess
  type :: pt
    integer(ih), pointer :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "h0"
  if (present(success)) success = lsuccess
  if (.not. lsuccess) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs = p%p
end subroutine assign_get_h0
subroutine associate_get_h0(lhs,this,dealloc,success)
  integer(ih), pointer :: lhs
  type(var), intent(in) :: this
  logical, intent(in), optional :: dealloc
  logical, intent(out), optional :: success
  logical :: ldealloc, lsuccess
  type :: pt
    integer(ih), pointer :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "h0"
  if (present(success)) success = lsuccess
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  ! there is one problem, say if lhs is not nullified...
  if (ldealloc.and.associated(lhs)) then
     deallocate(lhs)
     nullify(lhs)
  end if
  if (.not. lsuccess ) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs => p%p
end subroutine associate_get_h0
subroutine associate_set_h0(this,rhs,dealloc)
  type(var), intent(inout) :: this
  integer(ih), intent(in), target :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    integer(ih), pointer :: p => null()
  end type
  type(pt) :: p
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  this%t = "h0"
  p%p => rhs
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
end subroutine associate_set_h0
pure function associatd_l_h0(lhs,this) result(ret)
  integer(ih), pointer :: lhs
  type(var), intent(in) :: this
  logical :: ret
  type :: pt
    integer(ih), pointer :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "h0"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(lhs,p%p)
  endif
end function associatd_l_h0
pure function associatd_r_h0(this,rhs) result(ret)
  type(var), intent(in) :: this
  integer(ih), pointer :: rhs
  logical :: ret
  type :: pt
    integer(ih), pointer :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "h0"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(p%p,rhs)
  endif
end function associatd_r_h0
! All boolean functions
subroutine assign_set_h1(this,rhs,dealloc)
  type(var), intent(inout) :: this
  integer(ih), intent(in), dimension(:) :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    integer(ih), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  ! ASSIGNMENT in fortran is per default destructive
  ldealloc = .true.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  ! With pointer transfer we need to deallocate
  ! else bounds might change...
  this%t = "h1"
  allocate(p%p(size(rhs))) ! allocate space
  p%p = rhs ! copy data over
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
  ! We already have shipped it
  nullify(p%p)
end subroutine assign_set_h1
subroutine assign_get_h1(lhs,this,success)
  integer(ih), intent(out), dimension(:) :: lhs
  type(var), intent(in) :: this
  logical, intent(out), optional :: success
  logical :: lsuccess
  type :: pt
    integer(ih), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "h1"
  if (lsuccess) then
    p = transfer(this%enc,p) ! retrieve pointer encoding
    lsuccess = all(shape(p%p)==shape(lhs)) !&
     ! .and. all((lbound(p%p) == lbound(lhs))) &
     ! .and. all((ubound(p%p) == ubound(lhs)))
  end if
  if (present(success)) success = lsuccess
  if (.not. lsuccess) return
  lhs = p%p
end subroutine assign_get_h1
subroutine associate_get_h1(lhs,this,dealloc,success)
  integer(ih), pointer , dimension(:) :: lhs
  type(var), intent(in) :: this
  logical, intent(in), optional :: dealloc
  logical, intent(out), optional :: success
  logical :: ldealloc, lsuccess
  type :: pt
    integer(ih), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "h1"
  if (present(success)) success = lsuccess
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  ! there is one problem, say if lhs is not nullified...
  if (ldealloc.and.associated(lhs)) then
     deallocate(lhs)
     nullify(lhs)
  end if
  if (.not. lsuccess ) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs => p%p
end subroutine associate_get_h1
subroutine associate_set_h1(this,rhs,dealloc)
  type(var), intent(inout) :: this
  integer(ih), intent(in), dimension(:), target :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    integer(ih), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  this%t = "h1"
  p%p => rhs
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
end subroutine associate_set_h1
pure function associatd_l_h1(lhs,this) result(ret)
  integer(ih), pointer , dimension(:) :: lhs
  type(var), intent(in) :: this
  logical :: ret
  type :: pt
    integer(ih), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "h1"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(lhs,p%p)
  endif
end function associatd_l_h1
pure function associatd_r_h1(this,rhs) result(ret)
  type(var), intent(in) :: this
  integer(ih), pointer , dimension(:) :: rhs
  logical :: ret
  type :: pt
    integer(ih), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "h1"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(p%p,rhs)
  endif
end function associatd_r_h1
! All boolean functions
subroutine assign_set_h2(this,rhs,dealloc)
  type(var), intent(inout) :: this
  integer(ih), intent(in), dimension(:,:) :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    integer(ih), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  ! ASSIGNMENT in fortran is per default destructive
  ldealloc = .true.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  ! With pointer transfer we need to deallocate
  ! else bounds might change...
  this%t = "h2"
  allocate(p%p(size(rhs,1),size(rhs,2))) ! allocate space
  p%p = rhs ! copy data over
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
  ! We already have shipped it
  nullify(p%p)
end subroutine assign_set_h2
subroutine assign_get_h2(lhs,this,success)
  integer(ih), intent(out), dimension(:,:) :: lhs
  type(var), intent(in) :: this
  logical, intent(out), optional :: success
  logical :: lsuccess
  type :: pt
    integer(ih), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "h2"
  if (lsuccess) then
    p = transfer(this%enc,p) ! retrieve pointer encoding
    lsuccess = all(shape(p%p)==shape(lhs)) !&
     ! .and. all((lbound(p%p) == lbound(lhs))) &
     ! .and. all((ubound(p%p) == ubound(lhs)))
  end if
  if (present(success)) success = lsuccess
  if (.not. lsuccess) return
  lhs = p%p
end subroutine assign_get_h2
subroutine associate_get_h2(lhs,this,dealloc,success)
  integer(ih), pointer , dimension(:,:) :: lhs
  type(var), intent(in) :: this
  logical, intent(in), optional :: dealloc
  logical, intent(out), optional :: success
  logical :: ldealloc, lsuccess
  type :: pt
    integer(ih), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "h2"
  if (present(success)) success = lsuccess
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  ! there is one problem, say if lhs is not nullified...
  if (ldealloc.and.associated(lhs)) then
     deallocate(lhs)
     nullify(lhs)
  end if
  if (.not. lsuccess ) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs => p%p
end subroutine associate_get_h2
subroutine associate_set_h2(this,rhs,dealloc)
  type(var), intent(inout) :: this
  integer(ih), intent(in), dimension(:,:), target :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    integer(ih), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  this%t = "h2"
  p%p => rhs
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
end subroutine associate_set_h2
pure function associatd_l_h2(lhs,this) result(ret)
  integer(ih), pointer , dimension(:,:) :: lhs
  type(var), intent(in) :: this
  logical :: ret
  type :: pt
    integer(ih), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "h2"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(lhs,p%p)
  endif
end function associatd_l_h2
pure function associatd_r_h2(this,rhs) result(ret)
  type(var), intent(in) :: this
  integer(ih), pointer , dimension(:,:) :: rhs
  logical :: ret
  type :: pt
    integer(ih), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "h2"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(p%p,rhs)
  endif
end function associatd_r_h2
! All boolean functions
subroutine assign_set_i0(this,rhs,dealloc)
  type(var), intent(inout) :: this
  integer(is), intent(in) :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    integer(is), pointer :: p => null()
  end type
  type(pt) :: p
  ! ASSIGNMENT in fortran is per default destructive
  ldealloc = .true.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  ! With pointer transfer we need to deallocate
  ! else bounds might change...
  this%t = "i0"
  allocate(p%p) ! allocate space
  p%p = rhs ! copy data over
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
  ! We already have shipped it
  nullify(p%p)
end subroutine assign_set_i0
subroutine assign_get_i0(lhs,this,success)
  integer(is), intent(out) :: lhs
  type(var), intent(in) :: this
  logical, intent(out), optional :: success
  logical :: lsuccess
  type :: pt
    integer(is), pointer :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "i0"
  if (present(success)) success = lsuccess
  if (.not. lsuccess) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs = p%p
end subroutine assign_get_i0
subroutine associate_get_i0(lhs,this,dealloc,success)
  integer(is), pointer :: lhs
  type(var), intent(in) :: this
  logical, intent(in), optional :: dealloc
  logical, intent(out), optional :: success
  logical :: ldealloc, lsuccess
  type :: pt
    integer(is), pointer :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "i0"
  if (present(success)) success = lsuccess
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  ! there is one problem, say if lhs is not nullified...
  if (ldealloc.and.associated(lhs)) then
     deallocate(lhs)
     nullify(lhs)
  end if
  if (.not. lsuccess ) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs => p%p
end subroutine associate_get_i0
subroutine associate_set_i0(this,rhs,dealloc)
  type(var), intent(inout) :: this
  integer(is), intent(in), target :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    integer(is), pointer :: p => null()
  end type
  type(pt) :: p
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  this%t = "i0"
  p%p => rhs
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
end subroutine associate_set_i0
pure function associatd_l_i0(lhs,this) result(ret)
  integer(is), pointer :: lhs
  type(var), intent(in) :: this
  logical :: ret
  type :: pt
    integer(is), pointer :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "i0"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(lhs,p%p)
  endif
end function associatd_l_i0
pure function associatd_r_i0(this,rhs) result(ret)
  type(var), intent(in) :: this
  integer(is), pointer :: rhs
  logical :: ret
  type :: pt
    integer(is), pointer :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "i0"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(p%p,rhs)
  endif
end function associatd_r_i0
! All boolean functions
subroutine assign_set_i1(this,rhs,dealloc)
  type(var), intent(inout) :: this
  integer(is), intent(in), dimension(:) :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    integer(is), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  ! ASSIGNMENT in fortran is per default destructive
  ldealloc = .true.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  ! With pointer transfer we need to deallocate
  ! else bounds might change...
  this%t = "i1"
  allocate(p%p(size(rhs))) ! allocate space
  p%p = rhs ! copy data over
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
  ! We already have shipped it
  nullify(p%p)
end subroutine assign_set_i1
subroutine assign_get_i1(lhs,this,success)
  integer(is), intent(out), dimension(:) :: lhs
  type(var), intent(in) :: this
  logical, intent(out), optional :: success
  logical :: lsuccess
  type :: pt
    integer(is), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "i1"
  if (lsuccess) then
    p = transfer(this%enc,p) ! retrieve pointer encoding
    lsuccess = all(shape(p%p)==shape(lhs)) !&
     ! .and. all((lbound(p%p) == lbound(lhs))) &
     ! .and. all((ubound(p%p) == ubound(lhs)))
  end if
  if (present(success)) success = lsuccess
  if (.not. lsuccess) return
  lhs = p%p
end subroutine assign_get_i1
subroutine associate_get_i1(lhs,this,dealloc,success)
  integer(is), pointer , dimension(:) :: lhs
  type(var), intent(in) :: this
  logical, intent(in), optional :: dealloc
  logical, intent(out), optional :: success
  logical :: ldealloc, lsuccess
  type :: pt
    integer(is), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "i1"
  if (present(success)) success = lsuccess
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  ! there is one problem, say if lhs is not nullified...
  if (ldealloc.and.associated(lhs)) then
     deallocate(lhs)
     nullify(lhs)
  end if
  if (.not. lsuccess ) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs => p%p
end subroutine associate_get_i1
subroutine associate_set_i1(this,rhs,dealloc)
  type(var), intent(inout) :: this
  integer(is), intent(in), dimension(:), target :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    integer(is), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  this%t = "i1"
  p%p => rhs
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
end subroutine associate_set_i1
pure function associatd_l_i1(lhs,this) result(ret)
  integer(is), pointer , dimension(:) :: lhs
  type(var), intent(in) :: this
  logical :: ret
  type :: pt
    integer(is), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "i1"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(lhs,p%p)
  endif
end function associatd_l_i1
pure function associatd_r_i1(this,rhs) result(ret)
  type(var), intent(in) :: this
  integer(is), pointer , dimension(:) :: rhs
  logical :: ret
  type :: pt
    integer(is), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "i1"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(p%p,rhs)
  endif
end function associatd_r_i1
! All boolean functions
subroutine assign_set_i2(this,rhs,dealloc)
  type(var), intent(inout) :: this
  integer(is), intent(in), dimension(:,:) :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    integer(is), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  ! ASSIGNMENT in fortran is per default destructive
  ldealloc = .true.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  ! With pointer transfer we need to deallocate
  ! else bounds might change...
  this%t = "i2"
  allocate(p%p(size(rhs,1),size(rhs,2))) ! allocate space
  p%p = rhs ! copy data over
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
  ! We already have shipped it
  nullify(p%p)
end subroutine assign_set_i2
subroutine assign_get_i2(lhs,this,success)
  integer(is), intent(out), dimension(:,:) :: lhs
  type(var), intent(in) :: this
  logical, intent(out), optional :: success
  logical :: lsuccess
  type :: pt
    integer(is), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "i2"
  if (lsuccess) then
    p = transfer(this%enc,p) ! retrieve pointer encoding
    lsuccess = all(shape(p%p)==shape(lhs)) !&
     ! .and. all((lbound(p%p) == lbound(lhs))) &
     ! .and. all((ubound(p%p) == ubound(lhs)))
  end if
  if (present(success)) success = lsuccess
  if (.not. lsuccess) return
  lhs = p%p
end subroutine assign_get_i2
subroutine associate_get_i2(lhs,this,dealloc,success)
  integer(is), pointer , dimension(:,:) :: lhs
  type(var), intent(in) :: this
  logical, intent(in), optional :: dealloc
  logical, intent(out), optional :: success
  logical :: ldealloc, lsuccess
  type :: pt
    integer(is), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "i2"
  if (present(success)) success = lsuccess
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  ! there is one problem, say if lhs is not nullified...
  if (ldealloc.and.associated(lhs)) then
     deallocate(lhs)
     nullify(lhs)
  end if
  if (.not. lsuccess ) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs => p%p
end subroutine associate_get_i2
subroutine associate_set_i2(this,rhs,dealloc)
  type(var), intent(inout) :: this
  integer(is), intent(in), dimension(:,:), target :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    integer(is), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  this%t = "i2"
  p%p => rhs
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
end subroutine associate_set_i2
pure function associatd_l_i2(lhs,this) result(ret)
  integer(is), pointer , dimension(:,:) :: lhs
  type(var), intent(in) :: this
  logical :: ret
  type :: pt
    integer(is), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "i2"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(lhs,p%p)
  endif
end function associatd_l_i2
pure function associatd_r_i2(this,rhs) result(ret)
  type(var), intent(in) :: this
  integer(is), pointer , dimension(:,:) :: rhs
  logical :: ret
  type :: pt
    integer(is), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "i2"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(p%p,rhs)
  endif
end function associatd_r_i2
! All boolean functions
subroutine assign_set_l0(this,rhs,dealloc)
  type(var), intent(inout) :: this
  integer(il), intent(in) :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    integer(il), pointer :: p => null()
  end type
  type(pt) :: p
  ! ASSIGNMENT in fortran is per default destructive
  ldealloc = .true.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  ! With pointer transfer we need to deallocate
  ! else bounds might change...
  this%t = "l0"
  allocate(p%p) ! allocate space
  p%p = rhs ! copy data over
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
  ! We already have shipped it
  nullify(p%p)
end subroutine assign_set_l0
subroutine assign_get_l0(lhs,this,success)
  integer(il), intent(out) :: lhs
  type(var), intent(in) :: this
  logical, intent(out), optional :: success
  logical :: lsuccess
  type :: pt
    integer(il), pointer :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "l0"
  if (present(success)) success = lsuccess
  if (.not. lsuccess) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs = p%p
end subroutine assign_get_l0
subroutine associate_get_l0(lhs,this,dealloc,success)
  integer(il), pointer :: lhs
  type(var), intent(in) :: this
  logical, intent(in), optional :: dealloc
  logical, intent(out), optional :: success
  logical :: ldealloc, lsuccess
  type :: pt
    integer(il), pointer :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "l0"
  if (present(success)) success = lsuccess
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  ! there is one problem, say if lhs is not nullified...
  if (ldealloc.and.associated(lhs)) then
     deallocate(lhs)
     nullify(lhs)
  end if
  if (.not. lsuccess ) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs => p%p
end subroutine associate_get_l0
subroutine associate_set_l0(this,rhs,dealloc)
  type(var), intent(inout) :: this
  integer(il), intent(in), target :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    integer(il), pointer :: p => null()
  end type
  type(pt) :: p
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  this%t = "l0"
  p%p => rhs
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
end subroutine associate_set_l0
pure function associatd_l_l0(lhs,this) result(ret)
  integer(il), pointer :: lhs
  type(var), intent(in) :: this
  logical :: ret
  type :: pt
    integer(il), pointer :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "l0"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(lhs,p%p)
  endif
end function associatd_l_l0
pure function associatd_r_l0(this,rhs) result(ret)
  type(var), intent(in) :: this
  integer(il), pointer :: rhs
  logical :: ret
  type :: pt
    integer(il), pointer :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "l0"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(p%p,rhs)
  endif
end function associatd_r_l0
! All boolean functions
subroutine assign_set_l1(this,rhs,dealloc)
  type(var), intent(inout) :: this
  integer(il), intent(in), dimension(:) :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    integer(il), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  ! ASSIGNMENT in fortran is per default destructive
  ldealloc = .true.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  ! With pointer transfer we need to deallocate
  ! else bounds might change...
  this%t = "l1"
  allocate(p%p(size(rhs))) ! allocate space
  p%p = rhs ! copy data over
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
  ! We already have shipped it
  nullify(p%p)
end subroutine assign_set_l1
subroutine assign_get_l1(lhs,this,success)
  integer(il), intent(out), dimension(:) :: lhs
  type(var), intent(in) :: this
  logical, intent(out), optional :: success
  logical :: lsuccess
  type :: pt
    integer(il), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "l1"
  if (lsuccess) then
    p = transfer(this%enc,p) ! retrieve pointer encoding
    lsuccess = all(shape(p%p)==shape(lhs)) !&
     ! .and. all((lbound(p%p) == lbound(lhs))) &
     ! .and. all((ubound(p%p) == ubound(lhs)))
  end if
  if (present(success)) success = lsuccess
  if (.not. lsuccess) return
  lhs = p%p
end subroutine assign_get_l1
subroutine associate_get_l1(lhs,this,dealloc,success)
  integer(il), pointer , dimension(:) :: lhs
  type(var), intent(in) :: this
  logical, intent(in), optional :: dealloc
  logical, intent(out), optional :: success
  logical :: ldealloc, lsuccess
  type :: pt
    integer(il), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "l1"
  if (present(success)) success = lsuccess
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  ! there is one problem, say if lhs is not nullified...
  if (ldealloc.and.associated(lhs)) then
     deallocate(lhs)
     nullify(lhs)
  end if
  if (.not. lsuccess ) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs => p%p
end subroutine associate_get_l1
subroutine associate_set_l1(this,rhs,dealloc)
  type(var), intent(inout) :: this
  integer(il), intent(in), dimension(:), target :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    integer(il), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  this%t = "l1"
  p%p => rhs
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
end subroutine associate_set_l1
pure function associatd_l_l1(lhs,this) result(ret)
  integer(il), pointer , dimension(:) :: lhs
  type(var), intent(in) :: this
  logical :: ret
  type :: pt
    integer(il), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "l1"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(lhs,p%p)
  endif
end function associatd_l_l1
pure function associatd_r_l1(this,rhs) result(ret)
  type(var), intent(in) :: this
  integer(il), pointer , dimension(:) :: rhs
  logical :: ret
  type :: pt
    integer(il), pointer , dimension(:) :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "l1"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(p%p,rhs)
  endif
end function associatd_r_l1
! All boolean functions
subroutine assign_set_l2(this,rhs,dealloc)
  type(var), intent(inout) :: this
  integer(il), intent(in), dimension(:,:) :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    integer(il), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  ! ASSIGNMENT in fortran is per default destructive
  ldealloc = .true.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  ! With pointer transfer we need to deallocate
  ! else bounds might change...
  this%t = "l2"
  allocate(p%p(size(rhs,1),size(rhs,2))) ! allocate space
  p%p = rhs ! copy data over
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
  ! We already have shipped it
  nullify(p%p)
end subroutine assign_set_l2
subroutine assign_get_l2(lhs,this,success)
  integer(il), intent(out), dimension(:,:) :: lhs
  type(var), intent(in) :: this
  logical, intent(out), optional :: success
  logical :: lsuccess
  type :: pt
    integer(il), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "l2"
  if (lsuccess) then
    p = transfer(this%enc,p) ! retrieve pointer encoding
    lsuccess = all(shape(p%p)==shape(lhs)) !&
     ! .and. all((lbound(p%p) == lbound(lhs))) &
     ! .and. all((ubound(p%p) == ubound(lhs)))
  end if
  if (present(success)) success = lsuccess
  if (.not. lsuccess) return
  lhs = p%p
end subroutine assign_get_l2
subroutine associate_get_l2(lhs,this,dealloc,success)
  integer(il), pointer , dimension(:,:) :: lhs
  type(var), intent(in) :: this
  logical, intent(in), optional :: dealloc
  logical, intent(out), optional :: success
  logical :: ldealloc, lsuccess
  type :: pt
    integer(il), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  lsuccess = this%t == "l2"
  if (present(success)) success = lsuccess
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  ! there is one problem, say if lhs is not nullified...
  if (ldealloc.and.associated(lhs)) then
     deallocate(lhs)
     nullify(lhs)
  end if
  if (.not. lsuccess ) return
  p = transfer(this%enc,p) ! retrieve pointer encoding
  lhs => p%p
end subroutine associate_get_l2
subroutine associate_set_l2(this,rhs,dealloc)
  type(var), intent(inout) :: this
  integer(il), intent(in), dimension(:,:), target :: rhs
  logical, intent(in), optional :: dealloc
  logical :: ldealloc
  type :: pt
    integer(il), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  ! ASSOCIATION in fortran is per default non-destructive
  ldealloc = .false.
  if(present(dealloc))ldealloc = dealloc
  if (ldealloc) then
     call delete(this)
  else
     call nullify(this)
  end if
  this%t = "l2"
  p%p => rhs
  allocate(this%enc(size(transfer(p, local_enc_type)))) ! allocate encoding
  this%enc = transfer(p, local_enc_type) ! transfer pointer type to the encoding
end subroutine associate_set_l2
pure function associatd_l_l2(lhs,this) result(ret)
  integer(il), pointer , dimension(:,:) :: lhs
  type(var), intent(in) :: this
  logical :: ret
  type :: pt
    integer(il), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "l2"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(lhs,p%p)
  endif
end function associatd_l_l2
pure function associatd_r_l2(this,rhs) result(ret)
  type(var), intent(in) :: this
  integer(il), pointer , dimension(:,:) :: rhs
  logical :: ret
  type :: pt
    integer(il), pointer , dimension(:,:) :: p => null()
  end type
  type(pt) :: p
  ret = this%t == "l2"
  if (ret) then
     p = transfer(this%enc,p)
     ret = associated(p%p,rhs)
  endif
end function associatd_r_l2
! All boolean functions
end module variable
