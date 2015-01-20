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

! This code has been re-coded by:
!  Nick Papior Andersen to handle other than standard constrainst.

! Using CG/Broyden optimization *should* work as it is done per 
! force magnitude, and we correct for that.

module m_fixed

  use precision

  implicit none

  public :: init_fixed, resetFixedPointers
  public :: print_fixed
  public :: fixed

  ! Allows to check whether certain
  ! atoms are fixed
  public :: is_fixed, is_constr

  private

  integer, parameter :: TYPE_LEN = 10
  type tFix
     ! Number of atoms belonging to this fix
     integer :: n = 0
     ! the type of fixation
     character(TYPE_LEN) :: type
     ! Atoms in this fixation
     integer, allocatable :: a(:)
     ! direction of fixation
     ! Note that this is a single value for all
     ! If more atoms belong to one fixation 
     ! we fix them all to move along that direction
     real(dp) :: fix(3) = 0._dp
  end type tFix

  ! Containers for fixations
  integer, save :: N_fix = 0
  type(tFix), allocatable, save :: fixs(:)

  ! Stress fixation
  real(dp), save :: xs(6) = 0._dp

  ! Use constr routine
  logical, save :: use_constr = .false.

contains

  subroutine resetFixedPointers( )

    integer :: i

    if ( N_fix == 0 ) return

    do i = 1 , N_fix
       deallocate(fixs(i)%a)
    enddo
    deallocate(fixs)
    N_fix = 0

  end subroutine resetFixedPointers

  ! **********************************************************************
  ! Reads and imposes required constraints to atomic displacements by
  ! making zero the forces in those directions. Constraints are specified
  ! by the FDF data block GeometryConstraints (see example below).
  ! Only types position and routine implemented in this version.
  ! Written by J.M.Soler. Feb., 1998
  ! Modified by P. Ordejon to output the total number of constraints
  !    imposed.  June, 2003
  ! Modified by Nick Papior Andersen for additional constraint types and
  !    printing of explicite constraints.
  !    January, 2015
  ! *********** INPUT ****************************************************
  ! real*8  cell(3,3)    : Lattice vectors
  ! real*8  stress( 3,3) : Stress tensor
  ! integer na           : Number of atoms
  ! integer isa(na)      : Species indexes
  ! real*8  amass(na)    : Atomic masses
  ! real*8  xa(3,na)     : Atomic cartesian coordinates
  ! real*8  fa(3,na)     : Atomic forces
  ! *********** OUTPUT ***************************************************
  ! real*8  cstress( 3,3) : Constrained stress tensor
  ! real*8  cfa(3,na)     : Constrained atomic forces
  ! integer ntcon         : Total number of position constraints imposed
  ! *********** UNITS ****************************************************
  ! Units are arbitrary but cell and xa must be in the same units
  ! *********** BEHAVIOUR ************************************************
  ! cstress may be the same physical array as stress, and cfa the same 
  ! as fa, i.e. it is allowed:
  !     call fixed( cell, stress, na, isa, amass, xa, fa, stress, fa, ntcon )
  ! NOTE: This is NOT good practice, and might be flagged by some compilers
  ! *********** USAGE ****************************************************
  ! Example: consider a diatomic molecule (atoms 1 and 2) above a surface, 
  ! represented by a slab of 5 atomic layers, with 10 atoms per layer.
  ! To fix the cell height, the slab's botom layer (last 10 atoms),
  ! the molecule's interatomic distance, its height above the surface
  ! (but not its azimutal orientation and lateral position), and the
  ! relative height of the two atoms:
  !
  !   %block GeometryConstraints
  !   cellside   c 
  !   cellangle  alpha  beta  gamma
  !   position  from -1 to -10
  !   rigid  1  2
  !   center 1  2   0.0  0.0  1.0
  !   routine constr
  !   stress 1  2  3
  !   %endblock GeometryConstraints
  !
  ! where constr is the following user-written subroutine:
  !
  !      subroutine constr( cell, na, isa, amass, xa, stress, fa, ntcon )
    !c real*8  cell(3,3)    : input lattice vectors (Bohr)
    !c integer na           : input number of atoms
    !c integer isa(na)      : input species indexes
    !c real*8  amass(na)    : input atomic masses
    !c real*8  xa(3,na)     : input atomic cartesian coordinates (Bohr)
    !c real*8  stress( 3,3) : input/output stress tensor (Ry/Bohr**3)
    !c real*8  fa(3,na)     : input/output atomic forces (Ry/Bohr)
    !c integer ntcon        : total number of position constraints, accumulative
    !      integer na, isa(na), ntcon
    !      double precision amass(na), cell(3,3), fa(3,na),
    !     .                 stress(3,3), xa(3,na), fz
    !      fz =  fa(3,1) + fa(3,2) 
    !      fa(3,1) = fz*amass(1)/(amass(1)+amass(2))
    !      fa(3,2) = fz*amass(2)/(amass(1)+amass(2))
    !      ntcon = ntcon+1
    !      end
    ! **********************************************************************
  
  subroutine fixed( cell, stress, na, isa, amass, xa, fa, cstress, cfa, ntcon , &
       magnitude_usage )

    use intrinsic_missing, only : VNORM

    integer,  intent(in)  :: na, isa(na)
    real(dp), intent(in)  :: amass(na), cell(3,3), fa(3,na), stress(3,3), xa(3,na)
    integer, intent(out)  :: ntcon
    real(dp), intent(out) :: cfa(3,na), cstress(3,3)

    ! Whether we use a routine to find the gradient
    ! of the magnitudes, in which case we do not
    ! want to scale with mass for certain optimizations.
    logical, intent(in), optional :: magnitude_usage

    integer :: ix, jx, ia
    integer :: if, N, i
    character(len=TYPE_LEN) :: namec
    real(dp) :: fxc, ca(3), am, cf(3)
    logical :: lmag_use

#ifdef DEBUG
    call write_debug( '    PRE fixed' )
#endif
    
    ! Copy stress and forces to output arrays 
    do ix = 1,3
       do jx = 1,3
          cstress(jx,ix) = stress(jx,ix)
       end do
    end do
    do ia = 1,na
       do ix = 1,3
          cfa(ix,ia) = fa(ix,ia)
       end do
    end do

    ! If we should call the routine
    if ( use_constr ) then
       call constr( cell, na, isa, amass, xa, cstress, cfa, ntcon )
    end if

    ! apply the stress constraint
    cstress(1,1) = cstress(1,1) - xs(1) * cstress(1,1)
    cstress(2,2) = cstress(2,2) - xs(2) * cstress(2,2)
    cstress(3,3) = cstress(3,3) - xs(3) * cstress(3,3)
    cstress(3,2) = cstress(3,2) - xs(4) * cstress(3,2)
    cstress(2,3) = cstress(2,3) - xs(4) * cstress(2,3)
    cstress(3,1) = cstress(3,1) - xs(5) * cstress(3,1)
    cstress(1,3) = cstress(1,3) - xs(5) * cstress(1,3)
    cstress(1,2) = cstress(1,2) - xs(6) * cstress(1,2)
    cstress(2,1) = cstress(2,1) - xs(6) * cstress(2,1)

    lmag_use = .false.
    if ( present(magnitude_usage) ) lmag_use = magnitude_usage

    ! The number of constraints
    ntcon = 0

    if ( N_fix > 0 ) then
    do if = 1 , N_fix

       ! apply the designated constraints

       N = fixs(if)%N

       namec = fixs(if)%type
       
       ! Select type of constraint
       if ( namec == 'pos' ) then

          ! this is the easy one, all atoms should not move
          do i = 1 , N

             ia = fixs(if)%a(i)
             cfa(:,ia) = 0._dp

          end do

          ! we remove N * 3 degrees of freedom
          ntcon = ntcon + 3 * N

       else if ( namec == 'pos-dir' ) then
          
          do i = 1 , N

             ia = fixs(if)%a(i)

             ! Calculate the force projected onto the fixation vector
             fxc = sum( fixs(if)%fix(:) * cfa(:,ia) )
             cfa(:,ia) = cfa(:,ia) - fxc * fixs(if)%fix(:)

          end do

          ! Only one direction is fixed, hence all atoms
          ! remove one degree of freedom
          ntcon = ntcon + N

       else if ( namec == 'center' ) then
          
          ! Maintain the same center of the molecule
          cf(:) = 0._dp

          ! we center the molecule by constraining the 
          ! average acceleration to 0 (for MD)
          ! or by subtracting the average force for magnitude used forces
          do i = 1 , N

             ia = fixs(if)%a(i)

             if ( lmag_use ) then
                cf(:) = cf(:) + cfa(:,ia)
             else
                cf(:) = cf(:) + cfa(:,ia) / amass(ia)
             end if

          end do

          ! This is the average force/acceleration
          cf(:) = cf(:) / real(N,dp)

          ! fix for the correct direction, project onto the fixation vector
          do i = 1 , N

             ia = fixs(if)%a(i)

             if ( lmag_use ) then
                ! subtract the average force so that the net-force is zero.
                cfa(:,ia) = cfa(:,ia) - cf(:)
             else
                cfa(:,ia) = cfa(:,ia) - cf(:) * amass(ia)
             end if

          end do

          ! in principle we do not constrain anything, we just scale 
          ! the forces so that the center of the system will be the same.

          if ( .not. lmag_use ) then
             call die('Center of motion does not currently work')
          end if

       else if ( namec == 'com' ) then
          
          ! Calculate center of mass of the molecule
          ca(:) = 0._dp
          am = 0._dp

          do i = 1 , N

             ia = fixs(if)%a(i)

             ! center of mass
             ca(:) = ca(:) + xa(:,ia) * amass(ia)
             am = am + amass(ia)

          end do

          ! get the center of mass
          ca(:) = ca(:) / am

          ! correct the forces so that the center of mass is the same

          call die('Center-of-mass does not currently work')
          
       else if ( namec == 'mol' ) then
          
          ! Calculate total force on the molecule
          ! this is done using the center-of-force method
          cf(:) = 0._dp
          am    = 0._dp
          
          do i = 1 , N
             
             ia = fixs(if)%a(i)

             ! sum the force
             cf(:) = cf(:) + cfa(:,ia)
             am    = am + amass(ia)

          end do

          if ( lmag_use ) then
             ! use the average mass.
             cf(:) = cf(:) / real(N,dp)
          else
             cf(:) = cf(:) / am
          end if

          do i = 1 , N

             ia = fixs(if)%a(i)
             
             if ( lmag_use ) then
                cfa(:,ia) = cf(:)
             else
                cfa(:,ia) = cf(:) * amass(ia)
             end if

          end do

          ! we constrain (N - 1) * 3 degrees of freedom (only the relative positions are
          ! constrained)
          ntcon = ntcon + ( N - 1 ) * 3

       else if ( namec == 'mol-dir' ) then
          
          ! Calculate total force on the molecule
          ! this is done using the center-of-force method
          cf(:) = 0._dp
          am    = 0._dp
          
          ! not so straight forward constraint
          do i = 1 , N
             
             ia = fixs(if)%a(i)

             ! sum the force
             cf(:) = cf(:) + cfa(:,ia)
             am = am + amass(ia)

          end do

          if ( lmag_use ) then
             ! use the average force
             cf(:) = cf(:) / real(N,dp)
          else
             ! use the average mass
             cf(:) = cf(:) / am
          end if

          ! only set the molecular force in the
          ! specified direction, project onto the fixation vector
          fxc = sum(fixs(if)%fix(:) * cf(:))
          cf(:) = cf(:) - fxc * fixs(if)%fix(:)

          do i = 1 , N

             ia = fixs(if)%a(i)
             
             if ( lmag_use ) then
                cfa(:,ia) = cf(:)
             else
                cfa(:,ia) = cf(:) * amass(ia)
             end if

          end do

          ! we constrain N atoms to not move in one direction
          ! hence we remove N degrees of freedom
          ntcon = ntcon + N
         
       else if ( namec == 'mol-max' ) then
          
          ! Calculate total force on the molecule
          ! this is done using the center-of-force method
          cf(:) = 0._dp
          am    = 0._dp
          
          do i = 1 , N
             
             ia = fixs(if)%a(i)

             if ( VNORM(cfa(:,ia)) > VNORM(cf) ) then
                cf(:) = cfa(:,ia)
                am    = amass(ia)
             end if

          end do

          if ( lmag_use ) then
             ! do nothing, we have the max force
          else
             cf(:) = cf(:) / am
          end if

          do i = 1 , N

             ia = fixs(if)%a(i)
             
             if ( lmag_use ) then
                cfa(:,ia) = cf(:)
             else
                cfa(:,ia) = cf(:) * amass(ia)
             end if

          end do

          ! we constrain (N - 1) * 3 degrees of freedom (only the relative positions are
          ! constrained)
          ntcon = ntcon + ( N - 1 ) * 3

       else if ( namec == 'mol-max-dir' ) then
          
          ! Calculate total force on the molecule
          ! this is done using the center-of-force method
          cf(:) = 0._dp
          am    = 0._dp
          
          ! not so straight forward constraint
          do i = 1 , N
             
             ia = fixs(if)%a(i)

             if ( VNORM(cfa(:,ia)) > VNORM(cf) ) then
                cf(:) = cfa(:,ia)
                am    = amass(ia)
             end if
          end do

          if ( lmag_use ) then
             ! do nothing, we have the max force
          else
             cf(:) = cf(:) / am
          end if

          ! only set the molecular force in the
          ! specified direction, project onto the fixation vector
          fxc = sum(fixs(if)%fix(:) * cf(:))
          cf(:) = cf(:) - fxc * fixs(if)%fix(:)

          do i = 1 , N

             ia = fixs(if)%a(i)
             
             if ( lmag_use ) then
                cfa(:,ia) = cf(:)
             else
                cfa(:,ia) = cf(:) * amass(ia)
             end if

          end do

          ! we constrain N atoms to not move in one direction
          ! hence we remove N degrees of freedom
          ntcon = ntcon + N
         
       end if
       
    end do
    end if

#ifdef DEBUG
    call write_debug( '    POS fixed' )
#endif

  end subroutine fixed

  subroutine print_fixed( )

    use parallel, only : Node
    use m_region

    integer :: if
    character(len=70) :: name

    type(tRgn) :: r


    if ( Node /= 0 ) return

    write(*,*) ! newline

    if ( use_constr ) then
       write(*,'(a)') 'siesta: Constraints using custom constr routine'
    end if

    if ( any(xs /= 0._dp) ) then
       write(*,'(a)') 'siesta: Constraint (stress):'
       write(*,'(3(2(tr2,e11.4),/))') xs(:)
    end if

    if ( N_fix == 0 ) return

    if ( N_fix > 1 ) then
       write(*,'(a)') 'siesta: Constraints applied in the following order:'
    end if

    do if = 1 , N_fix 

       call rgn_list(r,fixs(if)%n,fixs(if)%a)
       r%name = trim(fixs(if)%type)
       name = 'siesta: Constraint'
       if ( index(r%name,'dir') > 0 ) then
          ! the name should include the fixation vector
          write(name,'(a,2(tr1,e10.4,'',''),tr1,e10.4,a)') 'siesta: Constraint v=[',fixs(if)%fix,']'
       end if
       call rgn_print(r, name = trim(name), seq_max = 12)
       
    end do
    call rgn_delete(r)

    write(*,*) ! newline

  end subroutine print_fixed
  
  subroutine init_fixed( na )

    use fdf
    use fdf_extra
    use intrinsic_missing, only : VNORM
    use parallel, only : IONode

    use m_region

    integer,  intent(in)  :: na

    ! Internal variables
    character(len=TYPE_LEN) :: namec

    integer :: ifix, i, ix, N

    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline

    type(tRgn) :: rr

    ! Initialise stress constraints to unconstrained state
    xs(:) = 0._dp

    ! No fixation atoms
    N_fix = 0

    ! Look for constraints data block
    if ( .not. fdf_block('GeometryConstraints',bfdf) ) return

#ifdef DEBUG
    call write_debug( '    PRE init_fixed' )
#endif

    ! First read in number of constraints
    do while ( fdf_bline(bfdf,pline) )

       ! If no names exists, we have no constraint line
       if ( fdf_bnnames(pline) == 0 ) cycle

       namec = fdf_bnames(pline,1)

       if ( leqi(namec,'position') .or. leqi(namec,'atom') .or. &
            leqi(namec,'molecule') .or. leqi(namec,'molecule-max') .or. &
            leqi(namec,'center') .or. &
            leqi(namec,'center-of-mass') ) then

          ! All these names belong to atomic constraints
          N_fix = N_fix + 1
          
       end if

    end do

    ! Allocate space
    allocate(fixs(N_fix))

    call fdf_brewind(bfdf)

    ! First read in number of constraints
    ifix = 0
    do while ( fdf_bline(bfdf,pline) )

       ! If no names exists, we have no constraint line
       if ( fdf_bnnames(pline) == 0 ) cycle

       namec = fdf_bnames(pline,1)

       ! Take those not specific to atomic positions
       if ( leqi(namec,'stress') ) then
          
          N = fdf_bnvalues(pline)

          do i = 1, N

             ix = nint(fdf_bvalues(pline,i))

             if ( ix < 1 .or. 6 < ix ) then
                call die('fixed: Stress restriction not &
                     &with expected input [1:6]')
             end if

             xs(ix) = 1._dp

          end do

       else if ( leqi(namec,'routine') ) then

          use_constr = .true.


       ! ****** Now we only look at atomic specifications for constraints ******
          
       else if ( leqi(namec,'position') .or. leqi(namec,'atom') ) then
          
          ifix = ifix + 1

          ! Create a list of atoms from this line
          call fdf_brange(pline,rr,1,na)

          fixs(ifix)%n = rr%n
          allocate(fixs(ifix)%a(rr%n))
          fixs(ifix)%a(:) = rr%r(:)

          ! Clean-up
          call rgn_delete(rr)

          ! Now read in the constraints if available
          
          ! Here are two variants
          ! If no reals exists on this line we will force
          ! the atoms to not move (i.e. fa = 0)
          N = fdf_bnreals(pline)
          if ( N == 0 ) then

             ! We fix them
             fixs(ifix)%type = 'pos'

          else if ( N == 3 ) then

             ! Fix the direction specified by the 3 reals
             fixs(ifix)%type = 'pos-dir'
             
             fixs(ifix)%fix(1) = fdf_breals(pline,1)
             fixs(ifix)%fix(2) = fdf_breals(pline,2)
             fixs(ifix)%fix(3) = fdf_breals(pline,3)

             fixs(ifix)%fix(:) = fixs(ifix)%fix(:) / &
                  VNORM( fixs(ifix)%fix(:) )

          else

             call die('You *must* specify 0 or 3 real values (not integers) &
                  &to do a constraint on atomic positions.')

          end if
             
       else if ( leqi(namec,'molecule') ) then

          ! We restrict the entire molecule to move "together"

          ifix = ifix + 1

          ! Create a list of atoms from this line
          call fdf_brange(pline,rr,1,na)

          fixs(ifix)%n = rr%n
          allocate(fixs(ifix)%a(rr%n))
          fixs(ifix)%a(:) = rr%r(:)

          ! Clean-up
          call rgn_delete(rr)

          ! Here are two variants
          ! If no reals exists on this line we will force
          ! the atoms to move relative
          N = fdf_bnreals(pline)
          if ( N == 0 ) then

             ! We just force the molecule to move together
             ! i.e. the relative positions will be the same.
             fixs(ifix)%type = 'mol'

          else if ( N == 3 ) then

             ! Fix the direction specified by the 3 reals
             ! Once we have calculated the relative displacements
             ! we fix the movement by these reals
             fixs(ifix)%type = 'mol-dir'
             
             fixs(ifix)%fix(1) = fdf_breals(pline,1)
             fixs(ifix)%fix(2) = fdf_breals(pline,2)
             fixs(ifix)%fix(3) = fdf_breals(pline,3)

             fixs(ifix)%fix(:) = fixs(ifix)%fix(:) / &
                  VNORM( fixs(ifix)%fix(:) )

          else

             call die('You *must* specify 0 or 3 real values (not integers) &
                  &to do a constraint on a molecule.')

          end if

       else if ( leqi(namec,'molecule-max') ) then

          ! We restrict the entire molecule to move "together"

          ifix = ifix + 1

          ! Create a list of atoms from this line
          call fdf_brange(pline,rr,1,na)

          fixs(ifix)%n = rr%n
          allocate(fixs(ifix)%a(rr%n))
          fixs(ifix)%a(:) = rr%r(:)

          ! Clean-up
          call rgn_delete(rr)

          ! Here are two variants
          ! If no reals exists on this line we will force
          ! the atoms to move relative
          N = fdf_bnreals(pline)
          if ( N == 0 ) then

             ! We just force the molecule to move together
             ! i.e. the relative positions will be the same.
             fixs(ifix)%type = 'mol-max'

          else if ( N == 3 ) then

             ! Fix the direction specified by the 3 reals
             ! Once we have calculated the relative displacements
             ! we fix the movement by these reals
             fixs(ifix)%type = 'mol-max-dir'
             
             fixs(ifix)%fix(1) = fdf_breals(pline,1)
             fixs(ifix)%fix(2) = fdf_breals(pline,2)
             fixs(ifix)%fix(3) = fdf_breals(pline,3)

             fixs(ifix)%fix(:) = fixs(ifix)%fix(:) / &
                  VNORM( fixs(ifix)%fix(:) )

          else

             call die('You *must* specify 0 or 3 real values (not integers) &
                  &to do a constraint on a molecule.')

          end if

       else if ( leqi(namec,'center-of-mass') ) then

          ! We restrict the center-of-mass for the specified region
          ! of atoms

          ifix = ifix + 1

          ! Create a list of atoms from this line
          call fdf_brange(pline,rr,1,na)

          fixs(ifix)%n = rr%n
          allocate(fixs(ifix)%a(rr%n))
          fixs(ifix)%a(:) = rr%r(:)

          ! Clean-up
          call rgn_delete(rr)

          fixs(ifix)%type = 'com'

       else if ( leqi(namec,'center') ) then

          ! We restrict the entire molecule to have the same center
          ! of coordinates

          ifix = ifix + 1

          ! Create a list of atoms from this line
          call fdf_brange(pline,rr,1,na)

          fixs(ifix)%n = rr%n
          allocate(fixs(ifix)%a(rr%n))
          fixs(ifix)%a(:) = rr%r(:)

          ! Clean-up
          call rgn_delete(rr)

          fixs(ifix)%type = 'center'

       else if ( IONode ) then
          
          write(*,'(2a)') 'siesta: Unrecognized geometry &
               &constraint: ',trim(namec)

       end if

    end do

#ifdef DEBUG
    call write_debug( '    POS init_fixed' )
#endif

  end subroutine init_fixed

  function is_fixed(ia) result(fixed)
    integer, intent(in) :: ia
    logical :: fixed, fi(3)
    integer :: if, i, N

    fixed = .false.
    fi(:) = .false.

    do if = 1 , N_fix
       
       N = fixs(if)%n

       do i = 1 , N

          if ( ia == fixs(if)%a(i) ) then

             if ( fixs(if)%type == 'pos' ) then
                fi(:) = .true.
             else if ( fixs(if)%fix(1) == 1._dp ) then
                fi(1) = .true.
             else if ( fixs(if)%fix(2) == 1._dp ) then
                fi(2) = .true.
             else if ( fixs(if)%fix(3) == 1._dp ) then
                fi(3) = .true.
             end if

             exit 

          end if

       end do

       fixed = all(fi)
       if ( fixed ) return

    end do

  end function is_fixed

  function is_constr(ia,type) result(constr)
    integer, intent(in) :: ia
    character(len=*), intent(in), optional :: type
    logical :: constr
    integer :: if, i, N

    constr = .false.

    do if = 1 , N_fix
       
       N = fixs(if)%n

       do i = 1 , N

          if ( ia == fixs(if)%a(i) ) then
             if ( present(type) ) then
                constr = fixs(if)%type == type
             else
                constr = .true.
             end if
             if ( constr ) return
          end if

       end do

    end do

  end function is_constr

end module m_fixed
