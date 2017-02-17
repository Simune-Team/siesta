! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine initatom(ns)

! Routine to initialize the Pseudopotentials and Atomic Orbitals.
! Substantially modified by Alberto Garcia (2000)
!
! The PAO and KB information can optionally be read from ASCII files
! (those produced by a standard run of Siesta or Base, but with the extension 
! renamed to '.input' instead of '.dump'), or from NetCDF files (if NetCDF
! is available). Note that there must be files for *all* species.
!
! This behavior is controlled by the FDF logicals 'user-basis' and
! 'user-basis-netcdf'.
!
! The old 'USER' basis type has been removed.
!
! This routine also outputs information about the basis specification
! determined by the routines in the 'basis_specs' modules. 
!
      use fdf
      use precision
      use basis_types, only: basis_specs_transfer, nsp
      use basis_types, only: deallocate_spec_arrays, initialize
      use basis_types, only: iz, lmxkb, nkbl, 
     &           erefkb, lmxo, nzeta, rco, 
     &           lambda, filtercut,
     &           atm_label, polorb, semic, nsemic,
     &           cnfigmx, charge, smass, basistype,
     &           rinn, vcte, qcoe, qyuk, qwid, split_norm
      use basis_types, only: write_basis_specs
      use basis_types, only: basis_def_t, basis_parameters
      use basis_specs, only: read_basis_specs

      use basis_io, only: read_basis_ascii, read_basis_netcdf
      use basis_io, only: dump_basis_ascii, dump_basis_netcdf
      use basis_io, only: dump_basis_xml

      use old_atmfuncs, only: nsmax, allocate_old_arrays, nkblsave
      use old_atmfuncs, only: clear_tables, deallocate_old_arrays
      use atom, only: atom_main, prinput
      use electrostatic, only: elec_corr_setup
      use atmparams, only: lmaxd, nkbmx, nsemx, nzetmx
      use atom_options, only: get_atom_options
      use ldau_specs, only: read_ldau_specs
      use ldau_specs, only: ldau_proj_gen

      use pseudopotential, only: pseudo_read
    
      use chemical

! CC RC  Added for the offSpOrb
      use m_spin,   only: spin
      use alloc,    only: de_alloc
      use parallel, only: IONode
! CC RC  Added for the offSpOrb


      implicit none
      integer,         intent(out) :: ns   ! Number of species
!     Internal variables ...................................................
      integer                      :: is
      logical                      :: user_basis, user_basis_netcdf
      logical :: req_init_setup
      type(basis_def_t),   pointer :: basp

      external atm_transfer

      call get_atom_options()


!     Reading input for the pseudopotentials and atomic orbitals
      write(6,'(/2a)') 
     &    'initatom: Reading input for the pseudopotentials ',
     &    'and atomic orbitals ----------'

      user_basis = fdf_boolean('user-basis',.false.)
      user_basis_netcdf = fdf_boolean('user-basis-netcdf',.false.)

      ! Create list of options NOT compatible with psf/vps file
      ! reads.
      req_init_setup = fdf_defined('LDAU.proj')
      ! Add any other dependencies here...

      ! Check that the user can perform a legal action
      req_init_setup = req_init_setup .and.
     & ( user_basis_netcdf .or. user_basis )
      if ( req_init_setup ) then
         call die('Reading PAOs and KBs from NetCDF/ascii files '//
     &'is not possible with LDAU.Proj')
      end if
      
      if (user_basis_netcdf) then
        write(6,'(/a)') 'Reading PAOs and KBs from NetCDF files...'
        call read_basis_netcdf(ns)
        call elec_corr_setup()
      else if (user_basis) then
       if ( spin%SO ) then  
          write(6,'(a)') ' initatom: Spin configuration = spin-orbit'
          call read_chemical_types()
          nsp = number_of_species()
          
          allocate(basis_parameters(nsp))
          do is = 1 , nsp
             call initialize(basis_parameters(is))
          enddo
          
          do is = 1 , nsp
             basp => basis_parameters(is)
             
             basp%label = species_label(is)
             call pseudo_read(basp%label,basp%pseudopotential)
          end do
       end if
       write(6,'(/a)') 'Reading PAOs and KBs from ascii files...'
       call read_basis_ascii(ns)
       call elec_corr_setup()
      else
        if ( IONode .and. spin%deb_offSO ) write(spin%iout_SO,'(a)') 
     &     '    initatom: Calling read_basis_specs... '
!       New routines in basis_specs and basis_types.
        call read_basis_specs()
        if ( IONode .and. spin%deb_offSO ) write(spin%iout_SO,'(a)') 
     &     '    initatom: Calling basis_specs_transfer... '
        call basis_specs_transfer()

!       Get the parameters for the generation of the LDA+U projectors
        call read_ldau_specs()

        nsmax = nsp             !! For old_atmfuncs
        if ( IONode .and. spin%deb_offSO ) write(spin%iout_SO,'(a)') 
     &     '    initatom: Calling allocate_old_arrays... '
        call allocate_old_arrays()
        if ( IONode .and. spin%deb_offSO ) write(spin%iout_SO,'(a)') 
     &     '    initatom: Calling clear_tables... '
        call clear_tables()

        do is = 1,nsp
         if ( IONode .and. spin%deb_offSO ) write(spin%iout_SO,'(a)') 
     &     '    initatom: Calling write_basis_specs... '
          call write_basis_specs(6,is)
          basp=>basis_parameters(is)
         if ( IONode .and. spin%deb_offSO ) write(spin%iout_SO,'(a,i3)')
     &     '    initatom: Calling ATOM_MAIN for specie ', is
          call ATOM_MAIN( iz(is), lmxkb(is), nkbl(0:lmaxd,is),
     &                    erefkb(1:nkbmx,0:lmaxd,is), lmxo(is),
     &                    nzeta(0:lmaxd,1:nsemx,is),
     &                    rco(1:nzetmx,0:lmaxd,1:nsemx,is),
     &                    lambda(1:nzetmx,0:lmaxd,1:nsemx,is),
     &                    atm_label(is), polorb(0:lmaxd,1:nsemx,is),
     &                    semic(is), nsemic(0:lmaxd,is),
     &                    cnfigmx(0:lmaxd,is), charge(is), smass(is),
     &                    basistype(is), is, rinn(0:lmaxd,1:nsemx,is),
     &                    vcte(0:lmaxd,1:nsemx,is),
     &                    qcoe(0:lmaxd,1:nsemx,is),
     &                    qyuk(0:lmaxd,1:nsemx,is),
     &                    qwid(0:lmaxd,1:nsemx,is),
     &                    split_norm(0:lmaxd,1:nsemx,is), 
     &                    filtercut(0:lmaxd,1:nsemx,is), basp)
!         Generate the projectors for the LDA+U simulations (if requested)
          call ldau_proj_gen(is)
        enddo 

        call prinput(nsp)

!       Create the new data structures for atmfuncs.
         if ( IONode .and. spin%deb_offSO .or. spin%deb_P ) 
     &    write(spin%iout_SO,'(a)') 
     &     '    initatom: Calling atm_transfer... '
        call atm_transfer()

         if ( IONode .and. spin%deb_offSO .or. spin%deb_P) 
     &    write(spin%iout_SO,'(a)') 
     &     '    initatom: Calling deallocate_old_arrays... '
        call deallocate_old_arrays()
         if ( IONode .and. spin%deb_offSO .or. spin%deb_P) 
     &    write(spin%iout_SO,'(a)') 
     &     '    initatom: Leaving deallocate_old_arrays... '
        call elec_corr_setup()
        ns = nsp               ! Set number of species for main program

      endif

      if ( IONode .and. spin%deb_offSO .or. spin%deb_P) 
     &     write(spin%iout_SO,'(a)') 
     &     '    initatom: Calling dump_basis_ascii... '
      call dump_basis_ascii()

! CC RC  Added for the offSpOrb
      if ( spin%SO ) then
       write(6,*) 
     & ' WARNING: NETCDF is not supported by the SO implementation'
      else
       call dump_basis_netcdf()
      endif
! CC RC  Added for the offSpOrb

      if ( IONode .and. spin%deb_offSO .or. spin%deb_P ) 
     &     write(spin%iout_SO,'(a)') 
     &     '    initatom: Calling dump_basis_xml... '
      call dump_basis_xml()

! CC RC  Added for the offSpOrb
! It is deallocated here instead of do it in deallocate_old_arrays
! Because the variable is needed to dump xml and ascii files
! Check if it'd be better to add an if condition to split the offSO
! from other type of calculation
      call de_alloc( nkblsave,    'nkblsave',    'initatom' )
! CC RC  Added for the offSpOrb

      if (.not. user_basis .and. .not. user_basis_netcdf) then
        call deallocate_spec_arrays()
      endif

      end subroutine initatom

