! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

program unfold

! Reads the .fdf, .ion, .psf and .WFSX files of a SIESTA calculation and generates
! unfolded bands.
! Ref: "Band unfolding made simple", S.Garcia-Mayo and J.M.Soler, Draft to be published
! S.Garcia-Mayo and J.M.Soler, Aug.2018

  use precision,    only: dp
  use chemical,     only: species_label
  use atm_types,    only: species_info, nspecies, species
  use atmfuncs,     only: lofio, mofio, nofis, rcut, rphiatm, zetafio
  use atom_options, only: get_atom_options
  use atom,         only: setup_atom_tables
  use basis_specs,  only: read_basis_specs
  use basis_types,  only: nsp, basis_specs_transfer, write_basis_specs
  use basis_io,     only: read_ion_ascii
  use sys,          only: die
  use fdf,          only: block_fdf, fdf_bintegers, fdf_bline, fdf_block, fdf_bmatch, &
                          fdf_bnames, fdf_bvalues, fdf_convfac, fdf_init, parsed_line
  use m_get_kpoints_scale, &
                    only: get_kpoints_scale
  use m_radfft,     only: radfft
  use spher_harm,   only: rlylm

  implicit none

  ! Internal parameters
  integer, parameter :: nr=1024         ! number of radial points for basis orbitals
  integer, parameter :: maxlines = 100  ! max number of unfolded band lines

  ! Internal variables
  integer  :: ierr, iostat, iorb, ir, isp, & 
              maxorb, nlines, ne(maxlines), norb, nq(maxlines)
  real(dp) :: dr, emax(maxlines), emin(maxlines), grad, gradv(3), &
              latConst, r, rc, qcell(3,3), qmax(3,maxlines)
  real(dp),allocatable:: phir(:,:,:), phik(:,:,:)
  logical  :: found
  character(len=50):: eunit, filein
  type(species_info),pointer :: spp
  type(block_fdf)            :: bfdf
  type(parsed_line), pointer :: pline

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
  allocate(phir(nsp,maxorb,nr))
  allocate(phik(nsp,maxorb,nr))

  ! Find atomic orbitals in real space
  do isp = 1,nsp
    do iorb = 1,nofis(isp)
      write(6,*) "# Orbital (#, l, z, m, rc):", &
        lofio(isp,iorb), zetafio(isp,iorb), mofio(isp,iorb), rcut(isp,iorb)
      rc = rcut(isp,iorb)
      dr = rc/nr
      do ir=0,nr
        r = dr*ir
        call rphiatm(isp,iorb,r,phir(isp,iorb,ir),grad)
      enddo
    enddo ! iorb
  enddo ! isp

  ! Fourier transform orbitals
  do isp = 1,nsp
    write(6,*) 'isp = ', isp 
    do iorb = 1,nofis(isp)
      write(6,*) 'iorb = ', iorb
      call radfft(lofio(isp,iorb),nr,rcut(isp,iorb),phir(isp,iorb,:),phik(isp,iorb,:))
!      write(6,*) phik(isp,iorb,:)
    enddo
  enddo

  ! Read UnfoldedBandLines
  call get_kpoints_scale('BandLinesScale',qcell,ierr)
  if (ierr/=0) &
    call die('unfold: ERROR calling get_kpoints_scale')
  found = fdf_block('UnfoldedBandLines',bfdf)
  if (.not.found) &
    call die('unfold ERROR: fdf block UnfoldedBandLines not found')
  nlines = 0
  do while( fdf_bline(bfdf,pline) )
    if ( fdf_bmatch(pline,'iivvvvvs') ) then
      nlines = nlines+1
      if (nlines>maxlines) &
        call die('unfold ERROR: parameter maxlines too small')
      nq(nlines) = fdf_bintegers(pline,1)
      ne(nlines) = fdf_bintegers(pline,2)
      qmax(:,nlines) = qcell(:,1)*fdf_bvalues(pline,1,after=2) &
                     + qcell(:,2)*fdf_bvalues(pline,2,after=2) &
                     + qcell(:,3)*fdf_bvalues(pline,3,after=2)
      eunit        = fdf_bnames(pline,1)          ! energy unit
      emin(nlines) = fdf_bvalues(pline,1,after=5)*fdf_convfac(eunit,'ry')
      emax(nlines) = fdf_bvalues(pline,2,after=5)*fdf_convfac(eunit,'ry')
      print*,'iline,nq,ne,qmax,emin,emax=', &
        nlines,nq(nlines),ne(nlines),qmax(:,nlines),emin(nlines),emax(nlines)
    else
      call die('unfold ERROR: wrong format in fdf block UnfoldedBandLines')
    endif
  enddo

end program unfold


