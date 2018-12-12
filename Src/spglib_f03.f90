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
! written by Matthieu Verstraete 2014-2018
! Bulk of the interfaces are copied from spglib_f08 in the spglib distribution, 
! by Keith Refson. I here avoid mapping C and f90 datatypes, only int char float
! reimplemented fortran initialization/destruction/printing routines

module spglib_f03

  use iso_c_binding, only:  c_char, c_int, c_double, c_ptr, c_f_pointer
  !use spglib_f08 ! this is also included in spglib and 
  !  could be a nice way to avoid recoding the interfaces for everything.
  !  however, it does imply making spglib with exactly the same fortran compiler as siesta
  implicit none

  ! this type now contains complementary derived information, with the main stuff stored in spglib's own datastructure above
  type :: SpglibDataset_ext
    ! saves stuff from which it was built as well
    real(c_double) :: celltransp(3,3)
    real(c_double), allocatable :: xred (:,:)
    real(c_double), allocatable :: spins_vect (:,:)
    real(c_double), allocatable :: spins_coll (:)
    integer, allocatable :: types (:)
    integer :: num_atom

    integer :: spacegroup_number
    integer :: hall_number
    integer :: bravais_number
    character(len=11) :: bravais_symbol
    character(len=11) :: international_symbol
    character(len=17) :: hall_symbol
    real(c_double)  :: transformation_matrix(3,3)
    real(c_double)  :: origin_shift(3)
    integer :: n_operations
    integer, allocatable :: rotations(:,:,:)
    real(c_double), allocatable :: translations(:,:)
    integer :: n_atoms
    integer, allocatable :: wyckoffs(:)
    integer, allocatable :: equivalent_atoms(:) !Beware mapping refers to positions starting at 0

    ! The following are extensions to Keith Refson's datastructure
    !integer, allocatable :: symops_cart(:,:,:)
    double precision, allocatable :: symops_cart(:,:,:)
    double precision, allocatable :: symops_recip_cart(:,:,:)
    integer, allocatable :: symops_recip(:,:,:)
    double precision, allocatable :: trans_cart(:,:)
  end type SpglibDataset_ext


  public :: SpglibDataset_ext
  public :: spgf_init, spgf_delete, spgf_print
  public :: spgf_mk_irred_k_grid

  private :: matr3inv

!
! spglib C For interfaces copied from spglib_f08, and some more added
!
   interface 
   
   function spg_get_symmetry( rotation, translation, max_size, lattice, &
                              & position, types, num_atom, symprec) bind(c)           
      import c_int, c_double
      integer(c_int), intent(inout) :: rotation(3,3,*)
      real(c_double), intent(inout) :: translation(3,*)
      integer(c_int), intent(in), value :: max_size
      real(c_double), intent(in) :: lattice(3,3), position(3,*)
      integer(c_int), intent(in) :: types(*)
      integer(c_int), intent(in), value :: num_atom
      real(c_double), intent(in), value :: symprec
      integer(c_int) :: spg_get_symmetry
   end function spg_get_symmetry
   
      
   function spgat_get_symmetry( rotation, translation, max_size, lattice, &
               & position, types, num_atom, symprec, angle_tolerance) bind(c)           
      import c_int, c_double
      integer(c_int), intent(inout) :: rotation(3,3,*)
      real(c_double), intent(inout) :: translation(3,*)
      integer(c_int), intent(in), value :: max_size
      real(c_double), intent(in) :: lattice(3,3), position(3,*)
      integer(c_int), intent(in) :: types(*)
      integer(c_int), intent(in), value :: num_atom
      real(c_double), intent(in), value :: symprec, angle_tolerance
      integer(c_int) :: spgat_get_symmetry
   end function spgat_get_symmetry

   
   function spg_get_symmetry_with_collinear_spin( rotation, translation, &
      & equivalent_atoms, max_size, lattice, position, types, spins, num_atom, symprec) bind(c)
      import c_int, c_double      
      integer(c_int), intent(inout) :: rotation(3,3,*)
      real(c_double), intent(inout) :: translation(3,*)
      integer(c_int), intent(inout) :: equivalent_atoms(*)
      integer(c_int), intent(in), value :: max_size
      real(c_double), intent(in) :: lattice(3,3), position(3,*)
      integer(c_int), intent(in) :: types(*)
      real(c_double), intent(in) :: spins(*)
      integer(c_int), intent(in), value :: num_atom
      real(c_double), intent(in), value :: symprec
      integer(c_int) :: spg_get_symmetry_with_collinear_spin
   end function spg_get_symmetry_with_collinear_spin

   
   function spgat_get_symmetry_with_collinear_spin( rotation, translation, &
      & equivalent_atoms, max_size, lattice, position, types, spins, num_atom, symprec, angle_tolerance) bind(c)
      import c_int, c_double      
      integer(c_int), intent(inout) :: rotation(3,3,*)
      real(c_double), intent(inout) :: translation(3,*)
      integer(c_int), intent(inout) :: equivalent_atoms(*)
      integer(c_int), intent(in), value :: max_size
      real(c_double), intent(in) :: lattice(3,3), position(3,*)
      integer(c_int), intent(in) :: types(*)
      real(c_double), intent(in) :: spins(*)
      integer(c_int), intent(in), value :: num_atom
      real(c_double), intent(in), value :: symprec, angle_tolerance
      integer(c_int) :: spgat_get_symmetry_with_collinear_spin
   end function spgat_get_symmetry_with_collinear_spin

   

   function spg_get_multiplicity( lattice, position, types, num_atom, symprec) bind(c)
      import c_int, c_double   
      real(c_double), intent(in) :: lattice(3,3), position(3,*)
      integer(c_int), intent(in) :: types(*)
      integer(c_int), intent(in), value :: num_atom
      real(c_double), intent(in), value :: symprec
      integer(c_int) :: spg_get_multiplicity
   end function spg_get_multiplicity
   

   function spgat_get_multiplicity( lattice, position, types, num_atom, symprec, angle_tolerance) bind(c)
      import c_int, c_double   
      real(c_double), intent(in) :: lattice(3,3), position(3,*)
      integer(c_int), intent(in) :: types(*)
      integer(c_int), intent(in), value :: num_atom
      real(c_double), intent(in), value :: symprec, angle_tolerance
      integer(c_int) :: spgat_get_multiplicity
   end function spgat_get_multiplicity

   
   function spg_get_smallest_lattice( smallest_lattice, lattice, symprec) bind(c)
      import c_int, c_double   
      real(c_double), intent(inout) :: smallest_lattice(3,3)
      real(c_double), intent(in) :: lattice(3,3)
      real(c_double), intent(in), value :: symprec
      integer(c_int) :: spg_get_smallest_lattice
   end function spg_get_smallest_lattice


   function spg_find_primitive(lattice, position, types, num_atom, symprec) bind(c)
      import c_int, c_double   
      real(c_double), intent(inout) :: lattice(3,3), position(3,*)
      integer(c_int), intent(inout) :: types(*)
      integer(c_int), intent(in), value :: num_atom
      real(c_double), intent(in), value :: symprec
      integer(c_int) :: spg_find_primitive
   end function spg_find_primitive

 
   function spgat_find_primitive(lattice, position, types, num_atom, symprec, angle_tolerance) bind(c)
      import c_int, c_double   
      real(c_double), intent(inout) :: lattice(3,3), position(3,*)
      integer(c_int), intent(inout) :: types(*)
      integer(c_int), intent(in), value :: num_atom
      real(c_double), intent(in), value :: symprec, angle_tolerance
      integer(c_int) :: spgat_find_primitive
   end function spgat_find_primitive

     
   
   function spg_get_international( symbol, lattice, position, types, num_atom, symprec) bind(c)
      import c_char, c_int, c_double
      character(kind=c_char), intent(out) :: symbol(11)
      real(c_double), intent(in) :: lattice(3,3), position(3, *)
      integer(c_int), intent(in) :: types(*)
      integer(c_int), intent(in), value :: num_atom
      real(c_double), intent(in), value :: symprec   
      integer(c_int) :: spg_get_international ! the number corresponding to 'symbol'. 0 on failure
   end function spg_get_international
   
   
   function spgat_get_international( symbol, lattice, position, types, num_atom, symprec, angle_tolerance) bind(c)
      import c_char, c_int, c_double
      character(kind=c_char), intent(out) :: symbol(11)
      real(c_double), intent(in) :: lattice(3,3), position(3, *)
      integer(c_int), intent(in) :: types(*)
      integer(c_int), intent(in), value :: num_atom
      real(c_double), intent(in), value :: symprec, angle_tolerance
      integer(c_int) :: spgat_get_international ! the number corresponding to 'symbol'. 0 on failure
   end function spgat_get_international
   
   
   

   function spg_get_schoenflies( symbol, lattice, position, types, num_atom, symprec) bind(c)
      import c_char, c_int, c_double
      character(kind=c_char), intent(out) :: symbol(10)
      real(c_double), intent(in) :: lattice(3,3), position(3, *)
      integer(c_int), intent(in) :: types(*)
      integer(c_int), intent(in), value :: num_atom
      real(c_double), intent(in), value :: symprec   
      integer(c_int) :: spg_get_schoenflies ! the number corresponding to 'symbol'. 0 on failure
   end function spg_get_schoenflies

   function spgat_get_schoenflies( symbol, lattice, position, types, num_atom, symprec, angle_tolerance) bind(c)
      import c_char, c_int, c_double
      character(kind=c_char), intent(out) :: symbol(10)
      real(c_double), intent(in) :: lattice(3,3), position(3, *)
      integer(c_int), intent(in) :: types(*)
      integer(c_int), intent(in), value :: num_atom
      real(c_double), intent(in), value :: symprec, angle_tolerance
      integer(c_int) :: spgat_get_schoenflies ! the number corresponding to 'symbol'. 0 on failure
   end function spgat_get_schoenflies



   function spg_get_pointgroup( symbol, trans_mat, rotations, num_rotations) bind(c)
      import c_char, c_int, c_double   
      character(kind=c_char) :: symbol(6)
      integer(c_int), intent(inout) :: trans_mat(3,3)
      integer(c_int), intent(in) :: rotations(3,3,*)
      integer(c_int), intent(in), value :: num_rotations
      integer(c_int) :: spg_get_pointgroup
   end function spg_get_pointgroup
   
   
   function spg_refine_cell( lattice, position, types, num_atom, symprec) bind(c)
      import c_int, c_double   
      real(c_double), intent(inout) :: lattice(3,3), position(3,*)
      integer(c_int), intent(inout) :: types(*)
      integer(c_int), intent(in), value :: num_atom
      real(c_double), intent(in), value :: symprec
      integer(c_int) :: spg_refine_cell
   end function spg_refine_cell
      
   function spgat_refine_cell( lattice, position, types, num_atom, symprec, angle_tolerance) bind(c)
      import c_int, c_double   
      real(c_double), intent(inout) :: lattice(3,3), position(3,*)
      integer(c_int), intent(inout) :: types(*)
      integer(c_int), intent(in), value :: num_atom
      real(c_double), intent(in), value :: symprec, angle_tolerance
      integer(c_int) :: spgat_refine_cell
   end function spgat_refine_cell
   
   
   function spg_get_ir_kpoints( map, kpoints, num_kpoints, lattice, position, &
                           & types, num_atom, is_time_reversal, symprec) bind(c)
!   Beware the map refers to positions starting at 0
      import c_int, c_double                           
      integer(c_int), intent(inout) :: map(*)
      real(c_double), intent(in) :: kpoints(3,*)
      integer(c_int), intent(in), value :: num_kpoints
      real(c_double), intent(in) :: lattice(3,3), position(3,*)
      integer(c_int), intent(in) :: types(*)
      integer(c_int), intent(in), value :: num_atom
      integer(c_int), intent(in), value :: is_time_reversal
      real(c_double), intent(in), value :: symprec
      integer(c_int) :: spg_get_ir_kpoints
   end function spg_get_ir_kpoints
   

   function spg_get_ir_reciprocal_mesh(grid_point, map, mesh, &
      & is_shift, is_time_reversal, lattice, position, types, num_atom, symprec) bind(c)
      import c_int, c_double
!   Beware the map refers to positions starting at 0      
      integer(c_int), intent(out) :: grid_point(3, *), map(*) ! size is product(mesh)
      integer(c_int), intent(in) :: mesh(3), is_shift(3)
      integer(c_int), intent(in), value :: is_time_reversal
      real(c_double), intent(in) :: lattice(3,3), position(3, *)
      integer(c_int), intent(in) :: types(*)
      integer(c_int), intent(in), value :: num_atom
      real(c_double), intent(in), value :: symprec   
      integer(c_int) :: spg_get_ir_reciprocal_mesh ! the number of points in the reduced mesh
   end function spg_get_ir_reciprocal_mesh


   function spg_get_stabilized_reciprocal_mesh(grid_point, map, mesh, is_shift, &
      & is_time_reversal, lattice, num_rot, rotations, num_q, qpoints) bind(c)
      import c_int, c_double      
!   Beware the map refers to positions starting at 0     
      integer(c_int), intent(inout) :: grid_point(3,*), map(*)
      integer(c_int), intent(in) :: mesh(3)
      integer(c_int), intent(in) :: is_shift(3)
      integer(c_int), intent(in), value :: is_time_reversal
      real(c_double), intent(in) :: lattice(3,3)
      integer(c_int), intent(in), value :: num_rot
      integer(c_int), intent(in) :: rotations(3,3,*)
      integer(c_int), intent(in), value :: num_q
      real(c_double), intent(in) :: qpoints(3,*)
      integer(c_int) :: spg_get_stabilized_reciprocal_mesh
   end function spg_get_stabilized_reciprocal_mesh
   
   
   function spg_get_triplets_reciprocal_mesh_at_q(weights, grid_points, third_q, &
      & grid_point, mesh, is_time_reversal, lattice, num_rot, rotations) bind(c)
      import c_int, c_double      
      integer(c_int), intent(inout) :: weights(*)
      integer(c_int), intent(inout) :: grid_points(3,*)
      integer(c_int), intent(inout) :: third_q(*)
      integer(c_int), intent(in), value :: grid_point
      integer(c_int), intent(in) :: mesh(3)
      integer(c_int), intent(in), value :: is_time_reversal
      real(c_double), intent(in) :: lattice(3,3)
      integer(c_int), intent(in), value :: num_rot
      integer(c_int), intent(in) :: rotations(3,3,*)
      integer(c_int) :: spg_get_triplets_reciprocal_mesh_at_q
   end function spg_get_triplets_reciprocal_mesh_at_q

   
   
   function spg_extract_triplets_reciprocal_mesh_at_q(triplets_at_q, &
         & weight_triplets_at_q, fixed_grid_number, num_triplets, triplets, &
         & mesh, is_time_reversal, lattice, num_rot, rotations) bind(c)
      import c_int, c_double         
      integer(c_int), intent(inout) :: triplets_at_q(3,*)
      integer(c_int), intent(inout) :: weight_triplets_at_q(*)
      integer(c_int), intent(in), value :: fixed_grid_number
      integer(c_int), intent(in), value :: num_triplets
      integer(c_int), intent(in) :: triplets(3,*)
      integer(c_int), intent(in) :: mesh(3)
      integer(c_int), intent(in), value :: is_time_reversal
      real(c_double), intent(in) :: lattice(3,3)
      integer(c_int), intent(in), value :: num_rot
      integer(c_int), intent(in) :: rotations(3,3,*)
      integer(c_int) :: spg_extract_triplets_reciprocal_mesh_at_q
   end function spg_extract_triplets_reciprocal_mesh_at_q
         
   function spg_get_hall_number_from_symmetry(rotation, translation, num_operations, &
      & symprec) bind(c)
      import c_int, c_double         
      integer(c_int), intent(in) :: rotation(3,3,*)
      real(c_double), intent(in) :: translation(3,*)
      integer(c_int), intent(in), value :: num_operations
      real(c_double), intent(in), value :: symprec   
      integer(c_int) :: spg_get_hall_number_from_symmetry
   end function spg_get_hall_number_from_symmetry

   end interface 
      
   
!   public :: SpglibDataset, spg_get_dataset, &
   public ::  &
      & spg_get_symmetry, spgat_get_symmetry, &
      & spg_get_symmetry_with_collinear_spin, spgat_get_symmetry_with_collinear_spin, &
      & spg_get_multiplicity, spgat_get_multiplicity, spg_get_smallest_lattice, &
      & spg_find_primitive, spgat_find_primitive, &
      & spg_get_international, spgat_get_international, &
      & spg_get_schoenflies, spgat_get_schoenflies, &
      & spg_get_pointgroup, spg_refine_cell, spgat_refine_cell, &
      & spg_get_ir_kpoints, spg_get_ir_reciprocal_mesh, spg_get_stabilized_reciprocal_mesh, &
      & spg_get_triplets_reciprocal_mesh_at_q, spg_extract_triplets_reciprocal_mesh_at_q, &
      & spg_get_hall_number_from_symmetry
   
  private

contains

  subroutine spgf_init (syms_this, cell, num_atom, types, xred, symprec )
! *****************************************************************
! Arguments:
! real*8  cell(3,3)    : input lattice vectors (Bohr)
! integer num_atom           : input number of atoms
! integer types(num_atom)      : input species indexes
! real*8  xred(3,num_atom)     : input atomic cartesian coordinates (Bohr)
! *****************************************************************
    implicit         none
    type(SpglibDataset_ext), intent(out) :: syms_this
    integer, intent(in)          :: num_atom, types(num_atom)
    double precision, intent(inout) :: cell(3,3), xred(3,num_atom)
    double precision, intent(in) :: symprec

! local vars
    integer :: max_size
    integer :: ii, isym, jj
    double precision :: cellinv(3,3)
    double precision :: dsymop(3,3)
    double precision :: dsymopinv(3,3)
    double precision :: zerovec(3) = (/0.0d0, 0.0d0, 0.0d0/)


    syms_this%num_atom = num_atom

    allocate (syms_this%types(num_atom))
    syms_this%types = types

    allocate (syms_this%spins_coll(num_atom))
    syms_this%spins_coll = 1.0d0

    allocate (syms_this%xred(3,num_atom))
    syms_this%xred = xred

    ! refines positions with highest symmetry found within symprec
    ! REMOVED 13/6/2014!!! - this returns the conventional cell!!!
    !call spg_refine_cell(cell, syms_this%xred, types, num_atom, symprec)
    !print *, 'xred refined'
    !print '(3E20.10)', syms_this%xred
    !print *, 'cell refined'
    !print '(3E20.10)', cell

    syms_this%celltransp = transpose(cell)

    syms_this%bravais_number = spg_get_international(syms_this%bravais_symbol, syms_this%celltransp(1,1), &
           & zerovec(1), types(1), 1, symprec)

    syms_this%spacegroup_number = spg_get_international(syms_this%international_symbol, syms_this%celltransp(1,1), &
           & syms_this%xred(1,1), types(1), num_atom, symprec)

    max_size = spg_get_multiplicity(syms_this%celltransp(1,1), syms_this%xred(1,1), types(1), num_atom, symprec)

    allocate(syms_this%rotations(1:3, 1:3, 1:max_size))
    allocate(syms_this%translations(1:3, 1:max_size))
    allocate(syms_this%equivalent_atoms(1:syms_this%num_atom))

    syms_this%n_operations =  spg_get_symmetry_with_collinear_spin(syms_this%rotations(1, 1, 1), &
       & syms_this%translations(1, 1), syms_this%equivalent_atoms(1), max_size, &
       & syms_this%celltransp(1, 1), syms_this%xred(1, 1), types(1), syms_this%spins_coll(1), &
       & num_atom, symprec)

    ! complete object with Hall number - how do I extract the symbol???
    syms_this%hall_number = spg_get_hall_number_from_symmetry(syms_this%rotations(1, 1, 1),&
       & syms_this%translations(1, 1), syms_this%n_operations, symprec)

    ! spglib returns row-major not column-major matrices!!! --DAS
    ! we transpose these after reading in: later use in fortran is x'(:) = Rotation(:,:) * x(:) + trans(:)
    do isym = 1, syms_this%n_operations
      syms_this%rotations(:,:,isym) = transpose(syms_this%rotations(:,:,isym))
    end do

! transpose needed or not? only visible on triclinic or hexagonal etc..
! Checked for 1 case with P6_3 c m (#185) comparing to abinit - a transpose is needed!!!

! get symops in cartesian coordinates - should still be integer -1 0 1 elements but store in floats
    allocate(syms_this%symops_cart(1:3, 1:3, 1:syms_this%n_operations))

    !call matr3inv(syms_this%celltransp, cellinv) ! NB: returns inverse transpose, which is what we want
    call matr3inv(cell, cellinv) ! NB: returns inverse transpose
    cellinv = transpose(cellinv)

    ! symcart = rprim * symred * gprim^T
    do ii = 1, syms_this%n_operations
      !dsymop = matmul(syms_this%celltransp, matmul(dble(syms_this%rotations(:, :, ii)), cellinv))
      dsymop = matmul(cell, matmul(dble(syms_this%rotations(:, :, ii)), cellinv))
      !syms_this%symops_cart(:, :, ii) = nint(dsymop)
      syms_this%symops_cart(:, :, ii) = dsymop

! it appears these matrices can be non integer for hexagonal systems...  feels wrong to me, but still:
!  for the moment, make them dble
!DEBUG
      if (  any(abs(dble(syms_this%symops_cart(:, :, ii)) - dsymop) > 1.e-10)  ) then
         print *, 'isym ', ii
         print *, 'rot        ', syms_this%rotations(:, :, ii)
         print *, 'rotcart    ', syms_this%symops_cart(:, :, ii)
         print *, 'rotcart_dp ', dsymop
         print *, 'error : cartesian symop element is not -1 0 +1'
         !stop
      end if
!END DEBUG
    end do

! get symops in reciprocal space = real space op^-1 ^T
    allocate(syms_this%symops_recip(1:3, 1:3, 1:syms_this%n_operations))
    allocate(syms_this%symops_recip_cart(1:3, 1:3, 1:syms_this%n_operations))
    do ii = 1, syms_this%n_operations
      dsymop = dble(syms_this%rotations(:, :, ii)) 
      call matr3inv(dsymop,dsymopinv) ! still inverse transpose here
      syms_this%symops_recip(:,:,ii) = nint(dsymopinv)
! DEBUG - perhaps comment out later
      if (any(abs(dble(syms_this%symops_recip(:, :, ii)) - dsymopinv(:,:)) > 1.e-10)) then
        stop 'error : recip symop element is not integer (normally -1 0 +1)'
      end if
! END DEBUG
      dsymop = matmul(transpose(cellinv), matmul(dble(syms_this%symops_recip(:, :, ii)), syms_this%celltransp))
      syms_this%symops_recip_cart(:, :, ii) = nint(dsymop)
! DEBUG - perhaps comment out later: one check is that the group of operations is closed.
!      do jj = 1, ii-1
!         print '(3(3E20.10, 3x))', matmul(syms_this%symops_recip_cart(:, :, ii), syms_this%symops_recip_cart(:, :, jj))
!      end do
! END DEBUG
    end do

    allocate(syms_this%trans_cart(1:3, 1:syms_this%n_operations))
    syms_this%trans_cart = matmul(syms_this%celltransp, syms_this%translations)

! TODO: add wyckoff positions ?

  end subroutine spgf_init



  subroutine spgf_delete( syms_this )
! *****************************************************************
! Arguments:
! SpglibDataset_ext syms_t     : object with symops etc...
! *****************************************************************
    implicit         none
    type(SpglibDataset_ext), intent(inout) :: syms_this

    if(allocated(syms_this%symops_cart))      deallocate (syms_this%symops_cart)
    if(allocated(syms_this%trans_cart))       deallocate (syms_this%trans_cart)

    if(allocated(syms_this%rotations))        deallocate (syms_this%rotations)
    if(allocated(syms_this%translations))     deallocate (syms_this%translations)
    if(allocated(syms_this%wyckoffs))         deallocate (syms_this%wyckoffs)
    if(allocated(syms_this%equivalent_atoms)) deallocate (syms_this%equivalent_atoms)

  end subroutine spgf_delete




  subroutine spgf_print( syms_this )
! *****************************************************************
! Arguments:
! SpglibDataset_ext syms_t     : object with symops etc...
! *****************************************************************
    implicit         none
    type(SpglibDataset_ext), intent(in) :: syms_this

! local
    integer :: ii
    
    write (*,*)
    write (*,'(a)') "----------------------------------------------------------------------"
    write (*,'(a)') " Output of crystal symmetries "
    write (*,*)
    write (*,'(a,I6)') " Int Sp Group number = ", syms_this%spacegroup_number
    write (*,'(2a)')   " Int Sp Group symbol = ", trim(syms_this%international_symbol)
    write (*,'(a,I6,2a)') ' Bravais lattice is ', syms_this%bravais_number, "  ",&
          & trim(syms_this%bravais_symbol)
    write (*,'(a,I6)') " Hall number = ", syms_this%hall_number
    !write (*,'(2a)') " Hall symbol = ", trim(syms_this%hall_symbol)

    write (*,*)
    write (*,'(a,I6)') " number of symops = ", syms_this%n_operations
    write (*,'(a)') " symops and translations in reduced coordinates = "
    do ii = 1, syms_this%n_operations
      write (*,'(3I6)') syms_this%rotations(:,1,ii)
      write (*,'(3I6)') syms_this%rotations(:,2,ii)
      write (*,'(3I6)') syms_this%rotations(:,3,ii)
      write (*, '(3(E20.10,2x))') syms_this%translations(:,ii)
      write (*,*)
    end do
    write (*,'(a)') " symops and translations in cartesian coordinates = "
    do ii = 1, syms_this%n_operations
      write (*,'(3E20.10)') syms_this%symops_cart(:,1,ii)
      write (*,'(3E20.10)') syms_this%symops_cart(:,2,ii)
      write (*,'(3E20.10)') syms_this%symops_cart(:,3,ii)
      write (*, '(3(E20.10,2x))') syms_this%trans_cart(:,ii)
      write (*,*)
    end do
    write (*,'(a)') "----------------------------------------------------------------------"
    write (*,*)
  end subroutine spgf_print

  subroutine spgf_mk_irred_k_grid(syms_this, symprec, is_time_reversal, &
&            nkgrid, shiftk, nkpt, kpt_list, irrk_map)
    implicit none
    type(SpglibDataset_ext), intent(in) :: syms_this
    double precision, intent(in) :: symprec
    integer, intent(in) :: is_time_reversal
    integer, intent(in) :: nkgrid(3)
    double precision, intent(in) :: shiftk(3)
    integer, intent(out) :: nkpt
    double precision, allocatable, intent(out) :: kpt_list(:,:)
    integer, allocatable, intent(out) :: irrk_map(:)

    integer :: nkmax, ik, idir
    integer, allocatable :: grid_point(:,:)


!    nkmax = product(nkgrid)
!    allocate (irrk_map(nkmax))
!    allocate (grid_point(3, nkmax))
!
!    nkpt = spg_get_ir_reciprocal_mesh(grid_point, irrk_map, nkgrid, &
!      & shiftk, is_time_reversal, syms_this%celltransp, &
!      & syms_this%xred, syms_this%types, syms_this%num_atom, symprec)
!
!    allocate (kpt_list(3,nkpt))
!    do ik = 1, nkpt
!      do idir = 1, 3
!        kpt_list(idir,ik) = dble(grid_point(idir,ik))/dble(nkgrid(idir))
!      end do
!    end do
!    deallocate (grid_point)

  end subroutine spgf_mk_irred_k_grid


!
! COPYRIGHT
!  Copyright (C) 1998-2018 ABINIT group (RC, XG, GMR, MG, JWZ)
!  This file is distributed under the terms of the
!  GNU General Public License, see ~abinit/COPYING
!  or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine matr3inv(aa, ait)

 implicit none

!Arguments ------------------------------------
!arrays
 double precision, intent(in) :: aa(3,3)
 double precision, intent(out) :: ait(3,3)

!Local variables-------------------------------
!scalars
 double precision :: dd,det,t1,t2,t3

! *************************************************************************

 t1 = aa(2,2) * aa(3,3) - aa(3,2) * aa(2,3)
 t2 = aa(3,2) * aa(1,3) - aa(1,2) * aa(3,3)
 t3 = aa(1,2) * aa(2,3) - aa(2,2) * aa(1,3)
 det  = aa(1,1) * t1 + aa(2,1) * t2 + aa(3,1) * t3

!Make sure matrix is not singular
 if (abs(det)>1.d-16) then
   dd=1.0d0/det
 else
   write(*, '(2a,2x)' ) 'Attempting to invert real(8) 3x3 array'
   write(*, '(9es16.8)' ) aa(:,:)
   write(*, '(a,es16.8,a)' ) '   ==> determinant=',det,' is zero.'
   stop
 end if

 ait(1,1) = t1 * dd
 ait(2,1) = t2 * dd
 ait(3,1) = t3 * dd
 ait(1,2) = (aa(3,1)*aa(2,3)-aa(2,1)*aa(3,3)) * dd
 ait(2,2) = (aa(1,1)*aa(3,3)-aa(3,1)*aa(1,3)) * dd
 ait(3,2) = (aa(2,1)*aa(1,3)-aa(1,1)*aa(2,3)) * dd
 ait(1,3) = (aa(2,1)*aa(3,2)-aa(3,1)*aa(2,2)) * dd
 ait(2,3) = (aa(3,1)*aa(1,2)-aa(1,1)*aa(3,2)) * dd
 ait(3,3) = (aa(1,1)*aa(2,2)-aa(2,1)*aa(1,2)) * dd

end subroutine matr3inv


end module spglib_f03

