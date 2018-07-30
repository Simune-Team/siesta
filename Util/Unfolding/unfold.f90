! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

program unfold

! Reads the .fdf, .ion, and .WFSX files of a SIESTA calculation and generates
! unfolded bands bands.
! Ref: "Band unfolding made simple", S.Garcia-Mayo and J.M.Soler, Draft to be published
! S.Garcia-Mayo and J.M.Soler, Aug.2018

  use precision,    only: dp
  use chemical,     only: species_label
  use atm_types,    only: species_info, nspecies, species
  use atmfuncs,     only: lofio, mofio, nofis, rcut, rphiatm, zetafio
  use atom_options, only: get_atom_options
  use atom,         only: setup_atom_tables
  use basis_specs,  only: read_basis_specs
  use basis_types,  only: nsp, basis_specs_transfer
  use basis_types,  only: write_basis_specs
  use m_getopts
  use basis_io,     only: read_ion_ascii
  use sys,          only: die
  use fdf

  implicit none

  ! Internal parameters
  integer, parameter :: nr=1024     ! number of radial points for basis orbitals

  ! Internal variables
  integer  :: iostat, iorb, ir, isp, maxorb, norb
  real(dp) :: dr, grad, gradv(3), r, rc
  real(dp),allocatable:: phi(:,:,:)
  type(species_info), pointer :: spp
  character(len=50):: filein

  ! Read and process input file
  filein = "stdin"
  call fdf_init(filein,'unfold_fdf.log')
  call get_atom_options()
  call read_xc_info()
  call read_basis_specs()
  call basis_specs_transfer()
  call setup_atom_tables(nsp)

  ! Allocate local arrays
  nspecies = nsp
  allocate(species(nsp))
  maxorb = 0               ! max. number of orbitals in any atom
  do isp=1,nsp
    spp => species(isp)
    spp%label = species_label(isp)
    spp%read_from_file = .true.
    call read_ion_ascii(spp)
    maxorb = max(maxorb,nofis(isp))
  enddo
  allocate(phi(nsp,maxorb,nr))

  ! Find atomic orbitals in real space
  do isp = 1,nsp
    do iorb = 1,nofis(isp)
      write(6,*) "# Orbital (#, l, z, m, rc):", &
        lofio(isp,iorb), zetafio(isp,iorb), mofio(isp,iorb), rcut(isp,iorb)
      rc = rcut(isp,iorb)
      dr = rc/nr
      do ir=0,nr
        r = dr*ir
        call rphiatm(isp,iorb,r,phi(isp,iorb,ir),grad)
      enddo
    enddo ! iorb
  enddo ! isp

end program unfold


