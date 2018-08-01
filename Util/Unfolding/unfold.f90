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

  use precision,    only: dp, sp
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
!  use files,        only: slabel

  implicit none

  ! Internal parameters
  integer, parameter :: nr=1024         ! number of radial points for basis orbitals
  integer, parameter :: maxlines = 100  ! max number of unfolded band lines

  ! Internal variables
  integer  :: ierr, iik, iispin, ik, indwf, iostat, iorb, ir, isp, & 
              ispin, iu, iw, j, lmax, maxorb, nk, nlines, ne(maxlines), &
              norb, nq(maxlines), nspin, nuotot, nwflist
  integer, allocatable:: iaorb(:), iphorb(:), cnfigfio(:)
  real(sp), allocatable :: psi(:,:)
  real(dp) :: dr, emax(maxlines), emin(maxlines), energy, grad, gradv(3), &
              latConst, modk, r, rc, qcell(3,3), qmax(3,maxlines)
  real(dp),allocatable:: k(:,:), phir(:,:,:), phik(:,:,:), ylm(:,:), &
                         gylm(:,:)
  logical  :: found, gamma
  character(len=50):: eunit, filein, fname
  character(len=20), allocatable :: symfio(:), labelfis(:)
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
    do iorb = 1,nofis(isp)
      call radfft(lofio(isp,iorb),nr,rcut(isp,iorb),phir(isp,iorb,:),phik(isp,iorb,:))
    enddo
  enddo
  print*,'Orbitals fourier-transformed'

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

  ! Read wavefunctions at k points of the BZ
  fname = 'si_bulk.fullBZ.WFSX'
!  fname = trim(slabel) // '.WFSX'                  NOT WORKING
!  call io_assing(iu)
  open(iu, file=fname, form='unformatted', status='old' )
  print*,'unfold, read wfsx'  ! ERROR HERE (last mod 31.07 20:24)
  rewind (iu)
  read(iu) nk, gamma
  allocate(k(nk,3))
  read(iu) nspin
  read(iu) nuotot
  if (gamma) then
    allocate(psi(1,nuotot))
  else
    allocate(psi(2,nuotot))
  endif
  allocate(iaorb(nuotot),labelfis(nuotot),iphorb(nuotot),cnfigfio(nuotot), &
           symfio(nuotot))
  read(iu) (iaorb(j),labelfis(j),iphorb(j),cnfigfio(j),symfio(j),j=1,nuotot)
  do iik = 1,nk
    do iispin = 1,nspin
      read(iu) ik,k(iik,1),k(iik,2),k(iik,3)   ! save all kpoints coord
      if (ik .ne. iik) &
        call die('unfold: ERROR in index of k-point')
      read(iu) ispin
      if (ispin .ne. iispin) &
        call die('unfold: ERROR in index of spin')
      read(iu) nwflist
      do iw = 1,nwflist
        read(iu) indwf
        read(iu) energy
        read(iu) (psi(1:,j), j=1,nuotot)
      enddo
    enddo
  enddo
  close (iu)

  
  ! Spherical harmonics: computed at k(nk,3) as we want to interpolate DOS
  ! rlylm( lmax, k(3), rly(0:) grly(1:,0:)) --- we give k normalized
  do iik = 1,nk
    modk = k(iik,1)**2 + k(iik,2)**2 + k(iik,3)**2
    k(iik,1) = k(iik,1)/modk
    k(iik,2) = k(iik,2)/modk
    k(iik,3) = k(iik,3)/modk
  enddo
  PRINT*,'START SPH HARM'
  lmax = 0
  do isp = 1,nsp
    lmax = max(lmax,lofio(isp,nofis(isp)))
  enddo
  allocate(ylm(nk,lmax*lmax),gylm(3,lmax*lmax))
  PRINT*,'CHECKPOINT'
  do iik = 1,nk
    call rlylm(lmax,k(iik,:),ylm(iik,:),gylm(:,:))
  enddo
  print*, 'lmax=',lmax, 'ylm=',ylm(1,2,3) 
  print*,'unfold, end spher_harm'


end program unfold

