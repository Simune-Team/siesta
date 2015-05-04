
! ---------------------------------------------------------------------
! MODULE m_planewavematrix
! In this module we define all the variables and subroutines required
! to compute the matrix elements of a plane wave in a basis of
! numerical atomic orbitals
! Implemented by J. Junquera, June-July 2013
! ---------------------------------------------------------------------

module m_planewavematrix

  CONTAINS

  subroutine planewavematrix( isigneikr, kptpw )

! ------------------------------ INPUT -----------------------------------------
!   integer     : isigneikr    ! Sign of the exponent of the planewave
!   real(dp)    : kptpw(3)     ! Wavevector of the planewave
!                                (Units in Bohr^-1) 
! ------------------------------ OUTPUT ----------------------------------------
!   complex(dp) : delkmat(:)   ! Matrix element of the plane wave in sparse
!                                format (same as in the Hamiltonian or Overlap
!                                matrices). 
!                                Note that the output is not part of the 
!                                variables of the subroutine, but is defined
!                                as a variable of this module.
! ------------------------------ BEHAVIOUR -------------------------------------
!   This subroutine calls dhscf to compute the matrix elements of a
!   plane wave, with wave vector kptpw and sign isigneikr.
!   The matrix elements are computed in a subroutine called delk,
!   essentially a copy of vmat.
!   delk is only called if isigneikr = +1 or -1. 
! ------------------------------------------------------------------------------

!   Used module variables
!   Most of the following parameters are only used to call dhscf.
    use precision,        only: dp       ! Real double precision type
    use m_spin,           only: nspin    ! Number of spin components 
    use siesta_geom,      only: na_u     ! Number of atoms in the unit cell
    use siesta_geom,      only: na_s     ! Number of atoms in the supercell
    use siesta_geom,      only: isa      ! Species index of each atom
    use siesta_geom,      only: xa       ! Atomic coordinates (in Bohr)
    use atomlist,         only: no_s     ! Number of orbitals in the supercell
    use atomlist,         only: iaorb    ! Atom to which each orbital belongs
    use atomlist,         only: iphorb   ! Orbital index within atom
    use atomlist,         only: no_l     ! Number of orbitals in the local node
    use atomlist,         only: no_u     ! Number of orbitals in the unit cell
    use atomlist,         only: indxua   ! Index of equivalent atom in unit cell
    use atomlist,         only: Datm     ! Occupations of basis orbitals
                                         !   in free atom
    use m_ntm,            only: ntm      ! Number of integration mesh divisions
                                         !   of each cell vector
    use sparse_matrices,  only: maxnh    ! Maximum number of orbitals
                                         !   interacting
    use sparse_matrices,  only: numh     ! Number of nonzero elements of each
                                         !   row of the hamiltonian matrix 
    use sparse_matrices,  only: listh    ! Nonzero hamiltonian-matrix elements
    use sparse_matrices,  only: listhptr ! Pointer to start of each row of the
                                         !   hamiltonian matrix
    use sparse_matrices,  only: Dscf     ! Density matrix in sparse form
    use sparse_matrices,  only: H        ! Hamiltonian matrix in sparse form
    use m_energies,       only: Enaatm   ! Integral of Vna * rhoatm
    use m_energies,       only: Enascf   ! Integral of Vna * rhoscf
    use m_energies,       only: Uatm     ! Harris hartree electron energy
    use m_energies,       only: Uscf     ! SCF hartree electron energy
    use m_energies,       only: DUscf    ! Electrostatic energy of 
                                         !   (rhoscf- rhoatm)
    use m_energies,       only: DUext    ! Interaction energy with 
                                         !   external  electric field
    use m_energies,       only: Exc      ! Exchange-correlation energy
    use m_energies,       only: Dxc      ! Integral((epsxc-Vxc)*Rho)
    use m_dipol,          only: dipol    ! Electric dipole
    use m_stress,         only: stress   ! Stress tensor
    use files,            only: filesOut_t ! derived type for output file names

!   Used module procedures
    use m_dhscf,          only: delk_wrapper    ! Mesh subroutine
    use alloc,            only: re_alloc ! Re-allocation routine
    use m_planewavematrixvar, only: wavevector ! Wave vector of the plane wave
 
    implicit none

    integer,  intent(in) :: isigneikr    ! Sign of the exponent of the planewave
    real(dp), intent(in) :: kptpw(3)     ! Wavevector of the planewave
                                         !   (Units in Bohr^-1) 

!   Local variables
    real(dp)             :: stressl(3,3) ! Local node part of stress
    real(dp), pointer    :: fal(:,:)     ! Local node part of atomic forces

    type(filesOut_t)     :: filesOut  ! blank output file names

!   Initialize the wave vector
    wavevector = kptpw

!   Initialize local forces and stresses 
!   Those are not used afterwards, so their values are not so important here
    nullify(fal)
    call re_alloc( fal, 1, 3, 1, na_u, 'fal', 'planewavematrix' )

    fal(1:3,1:na_u)  = 0.0_dp
    stressl(1:3,1:3) = 0.0_dp  
    call delk_wrapper(isigneikr, no_s, maxnh,  &
                            numh, listhptr, listh,  &
                            no_l,  no_u, iaorb, iphorb, isa )
 
  endsubroutine planewavematrix

endmodule m_planewavematrix



