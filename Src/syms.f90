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

  use iso_c_binding, only:  c_char, c_int, c_double, c_ptr, c_f_pointer
  !use spglib_f08 ! this is also included in spglib and 
  !  could be a nice way to avoid recoding the interfaces for everything.
  !  however, it does imply making spglib with exactly the same fortran compiler as siesta
  implicit none

! Copied from spglib_f08 for info
!      integer :: spacegroup_number
!      integer :: hall_number
!      character(len=11) :: international_symbol
!      character(len=17) :: hall_symbol
!      real(c_double)  :: transformation_matrix(3,3)
!      real(c_double)  :: origin_shift(3)
!      integer :: n_operations
!      integer, allocatable :: rotations(:,:,:)
!      real(c_double), allocatable :: translations(:,:)
!      integer :: n_atoms
!      integer, allocatable :: wyckoffs(:)
!      integer, allocatable :: equivalent_atoms(:) !Beware mapping refers to positions starting at 0
! end copy

  ! this type now contains complementary derived information, with the main stuff stored in spglib's own datastructure above
  type, public :: syms_type
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

    !integer, allocatable :: symops_cart(:,:,:)
    double precision, allocatable :: symops_cart(:,:,:)
    integer, allocatable :: symops_recip(:,:,:)
    double precision, allocatable :: trans_cart(:,:)
  end type syms_type

  ! for the moment this single instance of the syms is exported by the module.
  ! if you want to play with several, declare 1 per subroutine,
  ! and use it as the intent(out) argument to the init_syms subroutine
  type(syms_type), public :: syms_global

  public :: init_syms, delete_syms, print_syms
  public :: wrapvec_zero_one, wrap2zero_one
  public :: red2cart, cart2red

! spglib C For interfaces copied from spglib_f08
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
      & max_size, lattice, position, types, spins, num_atom, symprec) bind(c)
      import c_int, c_double      
      integer(c_int), intent(inout) :: rotation(3,3,*)
      real(c_double), intent(inout) :: translation(3,*)
      integer(c_int), intent(in), value :: max_size
      real(c_double), intent(in) :: lattice(3,3), position(3,*)
      integer(c_int), intent(in) :: types(*)
      real(c_double), intent(in) :: spins(*)
      integer(c_int), intent(in), value :: num_atom
      real(c_double), intent(in), value :: symprec
      integer(c_int) :: spg_get_symmetry_with_collinear_spin
   end function spg_get_symmetry_with_collinear_spin

   
   function spgat_get_symmetry_with_collinear_spin( rotation, translation, &
      & max_size, lattice, position, types, spins, num_atom, symprec, angle_tolerance) bind(c)
      import c_int, c_double      
      integer(c_int), intent(inout) :: rotation(3,3,*)
      real(c_double), intent(inout) :: translation(3,*)
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
   public :: &
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
    ! not necessary to wrap positions into unit cell: no effect on spgrp recognition by spglib

    ! refines positions with highest symmetry found within symprec
    ! REMOVED 13/6/2014!!! - this returns the conventional cell!!!
    !call spg_refine_cell(cell, xreda, isa, na, symprec)
    !print *, 'xred refined'
    !print '(3E20.10)', xreda
    !print *, 'cell refined'
    !print '(3E20.10)', cell

    celltransp = transpose(cell)

    syms_this%bravais_number = spg_get_international(syms_this%bravais_symbol, celltransp(1,1), &
           & zerovec(1), isa(1), 1, symprec)
! for information

    syms_this%spacegroup_number = spg_get_international(syms_this%international_symbol, celltransp(1,1), &
           & xreda(1,1), isa(1), na, symprec)



    max_size = spg_get_multiplicity(celltransp(1,1), xreda(1,1), isa(1), na, symprec)

    ! spglib returns row-major not column-major matrices!!! --DAS
    ! should we transpose these after reading in???
    allocate(syms_this%rotations(1:3, 1:3, 1:max_size))
    allocate(syms_this%translations(1:3, 1:max_size))

    syms_this%n_operations =  spg_get_symmetry(syms_this%rotations(1, 1, 1), syms_this%translations(1, 1), &
      max_size, celltransp(1, 1), xreda(1, 1), isa(1), na, symprec)

    ! complete object with Hall number - how do I extract the symbol???
    syms_this%hall_number = spg_get_hall_number_from_symmetry(syms_this%rotations(1, 1, 1),&
       & syms_this%translations(1, 1), syms_this%n_operations, symprec)

! transpose needed or not? only visible on triclinic or hexagonal etc..
! Checked for 1 case with P6_3 c m (#185) comparing to abinit - a transpose is needed!!!

! get symops in cartesian coordinates - should still be integer -1 0 1 elements but store in floats
    allocate(syms_this%symops_cart(1:3, 1:3, 1:syms_this%n_operations))
    call INVER(cell,cellinv,3,3,info)
    if (info .ne. 0) stop 'subroutine init_syms: error in inverse matrix of cell'

    do ii = 1, syms_this%n_operations
      syms_this%symops_cart(:, :, ii) = nint( matmul(cell, &
&       matmul(dble(syms_this%rotations(:, :, ii)), cellinv)) )

! it appears these matrices can be non integer for hexagonal systems...  feels wrong to me, but still:
!  for the moment, make them dble
!      if (  any(abs(dble(syms_this%symops_cart(:, :, ii)) &
!&       - matmul(cell, matmul(dble(syms_this%symops(:, :, ii)), cellinv))) > 1.e-10)  ) then
!         print *, ii, syms_this%symops(:, :, ii)
!         print *, matmul(cell, matmul(dble(syms_this%symops(:, :, ii)), cellinv))
!         stop 'error : cartesian symop element is not -1 0 +1'
!      end if

    end do

! get symops in reciprocal space = real space op^-1 ^T
    allocate(syms_this%symops_recip(1:3, 1:3, 1:syms_this%n_operations))
    do ii = 1, syms_this%n_operations
      dsymop = dble(syms_this%rotations(:, :, ii)) 
      call INVER(dsymop,dsymopinv,3,3,info)
      if (info .ne. 0) stop 'subroutine init_syms: error in inverse matrix of symop'
      dsymopinv = transpose (dsymopinv)
      syms_this%symops_recip(:,:,ii) = int(dsymopinv)
! DEBUG - perhaps comment out later
      if (any(abs(dble(syms_this%symops_recip(:, :, ii)) - dsymopinv(:,:)) > 1.e-10)) then
        stop 'error : recip symop element is not integer (normally -1 0 +1)'
      end if
! END DEBUG
    end do

    allocate(syms_this%trans_cart(1:3, 1:syms_this%n_operations))
    syms_this%trans_cart = matmul(cell, syms_this%translations)

  end subroutine init_syms



  subroutine delete_syms( syms_this )
! *****************************************************************
! Arguments:
! syms_type syms_t     : object with symops etc...
! *****************************************************************
    implicit         none
    type(syms_type), intent(inout) :: syms_this

    if(allocated(syms_this%symops_cart))          deallocate (syms_this%symops_cart)
    if(allocated(syms_this%trans_cart))           deallocate (syms_this%trans_cart)

    if(allocated(syms_this%rotations))        deallocate (syms_this%rotations)
    if(allocated(syms_this%translations))     deallocate (syms_this%translations)
    if(allocated(syms_this%wyckoffs))         deallocate (syms_this%wyckoffs)
    if(allocated(syms_this%equivalent_atoms)) deallocate (syms_this%equivalent_atoms)

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
    
    write (*,*)
    write (*,'(a)') "----------------------------------------------------------------------"
    write (*,'(a)') " Output of crystal symmetries "
    write (*,*)
    write (*,'(a,I6)') " Int Sp Group number = ", syms_this%spacegroup_number
    write (*,'(2a)') " Int Sp Group symbol = ", trim(syms_this%international_symbol)
    write (*,'(a,I6,2a)') ' for empty lattice space group is ', syms_this%bravais_number, "  ",&
          & trim(syms_this%bravais_symbol)
    write (*,'(a,I6)') " Hall number = ", syms_this%hall_number
    !write (*,'(2a)') " Hall symbol = ", trim(syms_this%hall_symbol)

    write (*,*)
    write (*,'(a,I6)') " number of symops = ", syms_this%n_operations
    write (*,'(a)') " symops in reduced coordinates = "
    do ii = 1, syms_this%n_operations
      write (*,'(3(3I10,2x))') syms_this%rotations(:,:,ii)
    end do
    write (*,'(a)') " trans in reduced coordinates = "
    do ii = 1, syms_this%n_operations
      write (*, '(3(E20.10,2x))') syms_this%translations(:,ii)
    end do
    write (*,'(a)') "----------------------------------------------------------------------"
    write (*,*)
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
