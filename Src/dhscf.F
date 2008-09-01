! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      subroutine dhscf_init(nspin, norb, iaorb, iphorb,
     &                      nuo, nuotot, nua, na,
     &                      isa, xa, indxua, ucell,
     &                      mscell, G2max, ntm,
     &                      maxnd, numd, listdptr, listd, datm,
     &                      Fal, stressl)
C     Modules
      use precision
      use parallel,       only : Node, Nodes
      use atmfuncs,       only : rcut, rcore
      use fdf
      use sys,            only : die
      use mesh,           only : xdsp, nsm, nsp
      use parsing
      use xcmod,          only : nXCfunc, XCauth
      use alloc,          only : re_alloc, de_alloc
      use siesta_options, only : harrisfun
      use moreMeshSubs,   only : setMeshDistr, distMeshData
      use m_dhscf,        only : rhopcc, rhoatm, Vna, spiral, rmax,
     &                           IsDiag, cell, ntopl, dvol, ntml, ntpl,
     &                           field, isefld, npcc, uharrs, scell,
     &                           shape, qspiral, g2mesh, nm, nml,
     &                           bcell, nmpl
      use iogrid_netcdf, only: set_box_limits
      use cellxc_mod,     only : setGGA
#ifdef MPI
      use mpi_siesta
#endif
      implicit none
      integer, intent(in)     :: nspin, norb, iaorb(norb), iphorb(norb),
     &                           nuo, nuotot, nua, na, isa(na),
     &                           indxua(na), mscell(3,3), maxnd,
     &                           numd(nuo), listdptr(nuo), listd(maxnd)
      real(dp), intent(in)    :: xa(3,na), ucell(3,3), G2max, datm(norb)
      integer, intent(inout)  :: ntm(3)
      real(dp), intent(inout) :: Fal(3,nua), stressl(3,3)

      real(dp), parameter     :: tiny  = 1.e-12_dp 
      integer                 :: io, ia, iphi, is, n, i, j
      integer                 :: nsc(3), nbcell, nsd
      real(dp)                :: DStres(3,3), volume
      logical                 :: samesh, samexa
      real(dp), external      :: volcel, ddot
      integer, parameter      :: vzero(3) = (/ 0, 0, 0 /)
      real(grid_p)            :: dummy_Drho(1,1), dummy_Vaux(1),
     &                           dummy_Vscf(1)
      logical, save           :: frstme = .true.   ! Keeps state
      real(grid_p),   pointer :: Vscf(:,:), rhoatm_seq(:)
      integer,        pointer :: numphi(:), numphi_par(:)

!------------------------------------------------------------------------- BEGIN
#ifdef DEBUG
      call write_debug( '    PRE dhscf_init' )
#endif
C ----------------------------------------------------------------------
C     General initialisation
C ----------------------------------------------------------------------
C     Start time counter
      call timer( 'DHSCF', 1 )

      nullify(Vscf,rhoatm_seq)
      if (frstme) then
        nullify( xdsp, rhopcc, Vna, rhoatm )
C       nsm lives in module mesh in meshsubs
        nsm = fdf_integer( 'MeshSubDivisions', 2 )
        nsm = max( nsm, 1 )

C       Set mesh sub-division variables & perform one off allocation
        nsp = nsm*nsm*nsm
        call re_alloc( xdsp, 1, 3, 1, nsp, 'xdsp', 'dhscf_init' )

C       Read spin-spiral wavevector (if defined)
        call readsp( qspiral, spiral )
        if (spiral .and. nspin.lt.4)
     $    call die('dhscf: ERROR: spiral defined but nspin < 4')
      endif   ! First time

      if (harrisfun) then
        do n = 1,nXCfunc
          if (.not.(leqi(XCauth(n),'PZ').or.leqi(XCauth(n),'CA'))) then
            if (Node .eq. 0) then
              write(6,100)
              write(6,101)
              write(6,'(a)')
     &         'dhscf:   Harris forces for a XC functional different of'
              write(6,'(a)')
     &         'dhscf:   Perdew-Zunger has not been implemented yet'
              write(6,102)
            endif
            call die()
          endif
        enddo
      endif

100   format(/,'dhscf: ',71(1h*))
101   format('dhscf:                  INPUT ERROR')
102   format('dhscf: ',71(1h*))

C ----------------------------------------------------------------------
C     Orbital initialisation : part 1
C ----------------------------------------------------------------------

C     Find the maximum orbital radius
      rmax = 0.0_dp
      do io = 1, norb
        ia   = iaorb(io)    ! Atomic index of each orbital
        iphi = iphorb(io)   ! Orbital index of each  orbital in its atom
        is   = isa(ia)      ! Species index of each atom
        rmax = max( rmax, rcut(is,iphi) )
      enddo

C ----------------------------------------------------------------------
C     Check whether mesh initialisation is needed
C ----------------------------------------------------------------------
      call SameMeshAndAtoms( na, xa, ucell, rmax, G2max, G2mesh,
     &                       samesh, samexa )

C ----------------------------------------------------------------------
C     Initialization due to change of mesh
C ----------------------------------------------------------------------
      if ( .not.samesh ) then
C       Start time counter for mesh initialization
        call timer( 'DHSCF1', 1 )

C ----------------------------------------------------------------------
C       Unit cell handling
C ----------------------------------------------------------------------
C       Find diagonal unit cell and supercell
        call digcel( ucell, mscell, cell, scell, nsc, IsDiag )
        if (.not.IsDiag) then
          if (Node.eq.0) then
            write(6,'(/,a,3(/,a,3f12.6,a,i6))')
     &      'DHSCF: WARNING: New shape of unit cell and supercell:',
     &     ('DHSCF:',(cell(i,j),i=1,3),'   x',nsc(j),j=1,3)
          endif
        endif

C       Find the system shape
        call shaper( cell, nua, isa, xa, shape, nbcell, bcell )

C       Find system volume
        volume = volcel( cell )

C ----------------------------------------------------------------------
C       Mesh initialization 
C ----------------------------------------------------------------------
#ifdef DEBUG
        call write_debug( '      PRE InitMesh' )
#endif
        call InitMesh( na, cell, norb, iaorb, iphorb, isa, rmax, 
     &                 G2max, G2mesh, nsc, nmpl, nm,
     &                 nml, ntm, ntml, ntpl, ntopl, dvol )

        ! Setup box descriptors for each processor,
        ! held in module iogrid_netcdf
        call set_box_limits(ntm,nsm)

C       Stop time counter for mesh initialization
        call timer( 'DHSCF1', 2 )
      endif
C ----------------------------------------------------------------------
C     End of mesh initialization 
C ----------------------------------------------------------------------


      if (.not.samexa) then
C ----------------------------------------------------------------------
C       Initialize atomic orbitals, density and potential
C ----------------------------------------------------------------------
C       Start time counter for atomic initializations
        call timer( 'DHSCF2', 1 )

C       Initialise quantities relating to the atom-mesh positioning
        call InitAtomMesh( na, xa)

C       Check if we need extencils in cellxc
        call setGGA( )

C       Compute the number of orbitals on the mesh and recompute the
C       partions for every processor in order to have a similar load
C       in each of them.
#ifdef MPI
        call setMeshDistr( 1, nsm, nsp, nml, nmpl, ntml, ntpl )
#endif
        nullify( numphi )
        call re_alloc( numphi, 1, nmpl, 'numphi', 'dhscf_init' )
#ifdef DEBUG
        call write_debug( '      PRE distriPhiOnMesh' )
#endif
        call distriPhiOnMesh( nm, nmpl, norb, iaorb, iphorb,
     &                        isa, numphi )
#ifdef MPI
        call setMeshDistr( 2, nsm, nsp, nml, nmpl, ntml, ntpl )
#endif
C       Find if there are partial-core-corrections for any atom
        npcc = 0
        do ia = 1,na
          if (rcore(isa(ia)) .gt. tiny) npcc = 1
        enddo

C       Find partial-core-correction energy density
C       Vscf and Vaux are not used here
        call re_alloc( rhopcc, 1, ntpl*npcc+1, 'rhopcc', 'dhscf_init' )
        if (npcc .eq. 1) then
          call PartialCoreOnMesh( na, isa, ntpl, rhopcc, indxua,
     &      nsd, dvol, volume, dummy_Vscf, dummy_Vaux, Fal, stressl,
     &      .false., .false. )
        endif

C       Find neutral-atom potential
C       Drho is not used here
        call re_alloc( Vna, 1, ntpl, 'Vna', 'dhscf_init' )
        call NeutralAtomOnMesh( na, isa, ntpl, Vna, indxua, dvol, 
     &                          volume, dummy_DRho, Fal, stressl, 
     &                          .false., .false. )

C       Calculate orbital values on mesh
#ifdef DEBUG
        call write_debug( '      PRE PhiOnMesh' )
#endif
#ifdef MPI
        call re_alloc( numphi_par, 1, nmpl, 'numphi_par', 'dhscf_init' )
        call distMeshData( 1, numphi, 2, numphi_par, 0 )
#else
        numphi_par => numphi
#endif
        call PhiOnMesh( nmpl, norb, iaorb, iphorb, isa, numphi_par )
#ifdef MPI
        call de_alloc( numphi_par, 'numphi_par', 'dhscf_init' )
#endif
        call de_alloc( numphi, 'numphi', 'dhscf_init' )
C ----------------------------------------------------------------------
C       Create sparse indexing for Dscf as needed for local mesh
C ----------------------------------------------------------------------
        if (Nodes.gt.1) then
          call CreateLocalDscfPointers( nmpl, nuotot, numd, listdptr, 
     &                                  listd )
        endif

C ----------------------------------------------------------------------
C       Calculate terms relating to the neutral atoms on the mesh
C ----------------------------------------------------------------------
C       Find Harris (sum of atomic) electron density
        call re_alloc( rhoatm, 1, ntpl, 'rhoatm', 'dhscf_init' )
        call rhooda( norb, nmpl, datm, rhoatm, iaorb, iphorb, isa )

C       Routine Poison should use the uniform data distribution
#ifdef MPI
        call setMeshDistr( 1, nsm, nsp, nml, nmpl, ntml, ntpl )
#endif
C       Allocate Temporal array
        call re_alloc( Vscf, 1, ntpl, 1, nspin, 'Vscf', 'dhscf_init' )
        call re_alloc( rhoatm_seq, 1, ntpl, 'rhoatm_seq', 'dhscf_init' )

C       Vscf is filled here but not used later
C       Uharrs is computed (and saved)
C       DStres is computed but not used later
        call distMeshData( 2, rhoatm, 1, rhoatm_seq, +1 )
#ifdef DEBUG
        call write_debug( '      PRE poison' )
#endif
        call poison( cell, ntml(1), ntml(2), ntml(3), ntm, rhoatm_seq,
     &               Uharrs, Vscf, DStres, nsm )
        call de_alloc( Vscf, 'Vscf', 'dhscf_init' )
        call de_alloc( rhoatm_seq, 'rhoatm_seq', 'dhscf_init' )

#ifdef MPI
        call setMeshDistr( 2, nsm, nsp, nml, nmpl, ntml, ntpl )
#endif

C       Stop time counter for atomic initializations
        call timer( 'DHSCF2', 2 )

C ----------------------------------------------------------------------
C     End of atomic initializations
C ----------------------------------------------------------------------
      endif

      if (frstme) then
C       Vscf is not used here because we pass "0,0,0" as mesh
        call efield( cell, nua, isa, xa, vzero, nsm,
     &               dummy_Vscf, field )
        IsEfld = .false.
        if (sqrt(ddot(3,field,1,field,1)) .gt. tiny) IsEfld = .true.
      endif
      frstme = .false.

      call timer( 'DHSCF', 2 )
#ifdef DEBUG
      call write_debug( '    POS dhscf_init' )
#endif
!--------------------------------------------------------------------------- END
      end subroutine dhscf_init

      subroutine dhscf( nspin, norb, iaorb, iphorb, nuo,
     &                  nuotot, nua, na, isa, xa, indxua,
     &                  ntm, ifa, istr, iHmat, filrho, fildrh, 
     &                  filevh, filevt, filepsch, filetoch, maxnd, numd,
     &                  listdptr, listd, Dscf, datm, maxnh, Hmat,
     &                  Enaatm, Enascf, Uatm, Uscf, DUscf, DUext, 
     &                  Exc, Dxc, dipol, stress, Fal, stressl )

C
C Calculates the self-consistent field contributions to Hamiltonian
C matrix elements, total energy and atomic forces.
C Coded by J.M. Soler, August 1996. July 1997.
C Modified by J.D. Gale, February 2000.
C
C ----------------------------------------------------------------------
C Input :
C ----------------------------------------------------------------------
C integer nspin         : Number of different spin polarisations
C                         nspin=1 => Unpolarized, nspin=2 => polarized
C                         nspin=4 => Noncollinear spin
C integer norb          : Total number of basis orbitals in supercell
C integer iaorb(norb)   : Atom to which each orbital belongs
C integer iphorb(norb)  : Orbital index (within atom) of each orbital
C integer nuo           : Number of orbitals in a unit cell in this node
C integer nuotot        : Number of orbitals in a unit cell
C integer nua           : Number of atoms in unit cell
C integer na            : Number of atoms in supercell
C integer isa(na)       : Species index of all atoms in supercell
C real*8  xa(3,na)      : Atomic positions of all atoms in supercell
C integer indxua        : Index of equivalent atom in unit cell
C integer ifa           : Switch which fixes whether the SCF contrib.
C                         to atomic forces is calculated and added to fa
C integer istr          : Switch which fixes whether the SCF contrib.
C                         to stress is calculated and added to stress
C integer iHmat         : Switch which fixes whether the Hmat matrix
C                         elements are calculated or not.
C character*(*) filrho  : Name of file to saving the electron density
C                         (If blank => not saved)
C character*(*) fildrh  : Name of file to saving Delta Rho (Rho - Rho_atoms)
C                         (If blank => not saved)
C character*(*) filevh  : Name of file to save electrostatic potential
C                         (If blank => not saved)
C character*(*) filevt  : Name of file to save total potential
C                         (If blank => not saved)
C character*(*) filepsch: Name of file to save the soft diffuse ionic charge
C                         (If blank => not saved)
C character*(*) filetoch: Name of file to save the total ionic+electronic charge
C                         (If blank => not saved)
C integer maxnd             : First dimension of listd and Dscf
C integer numd(nuo)         : Number of nonzero density-matrix
C                             elements for each matrix row
C integer listdptr(nuo)     : Pointer to start of rows of density-matrix
C integer listd(maxnd)      : Nonzero-density-matrix-element column
C                             indexes for each matrix row
C real*8  Dscf(maxnd,nspin): SCF density-matrix elements
C real*8  datm(norb)        : Harris density-matrix diagonal elements
C                             (atomic occupation charges of orbitals)
C integer maxnh             : First dimension of listh and Hmat
C real*8  Hmat(maxnh,nspin) : Hamiltonian matrix in sparse form,
C                             to which are added the matrix elements
C                                 <ORB_I | DeltaV | ORB_J>, where
C                             DeltaV = Vna + Vxc(SCF) + 
C                                      Vhartree(RhoSCF-RhoHarris)
C ----------------------------------------------------------------------
C Input/output :
C ----------------------------------------------------------------------
C integer ntm(3) : Number of mesh divisions of each cell
C                  vector, including subgrid.
C ----------------------------------------------------------------------
C Output :
C ----------------------------------------------------------------------
C real*8  Enaatm : Integral of Vna * rhoatm
C real*8  Enascf : Integral of Vna * rhoscf
C real*8  Uatm   : Harris hartree electron-interaction energy
C real*8  Uscf   : SCF hartree electron-interaction energy
C real*8  DUscf  : Electrostatic (Hartree) energy of
C                    (rhoscf - rhoatm) density
C real*8  DUext  : Interaction energy with external electric field
C real*8  Exc    : SCF exchange-correlation energy
C real*8  Dxc    : SCF double-counting correction to Exc
C                    Dxc = integral of ( (epsxc - Vxc) * Rho )
C                    All energies in Rydbergs
C real*8  dipol(3): Electric dipole (in a.u.)
C                   only when the system is a molecule
C ----------------------------------------------------------------------
C Input/output :
C ----------------------------------------------------------------------
C real*8  Fal(3,nua) : Atomic forces, to which the SCF contribution
C                        is added by this routine when ifa=1.
C                        the SCF contribution is minus the derivative
C                        of ( Enascf - Enaatm + DUscf + Exc ) with
C                        respect to atomic positions, in  Ry/Bohr.
C                        contributions local to this node
C real*8 stressl(3,3): Stress tensor, to which the SCF contribution
C                      is added by this routine when ifa=1.
C                      the SCF contribution is minus the derivative of
C                         ( Enascf - Enaatm + DUscf + Exc ) / volume
C                      with respect to the strain tensor, in Ry.
C                        contributions local to this node
C ----------------------------------------------------------------------
C Units :
C ----------------------------------------------------------------------
C Energies in Rydbergs
C Distances in Bohr
C ----------------------------------------------------------------------
C     Modules
      use precision
      use parallel,      only  : Node, Nodes
      use atmfuncs,      only  : rcut, rcore
      use fdf
      use sys,           only  : die
      use mesh,          only  : xdsp, nsm, nsp
      use parsing
      use m_iorho,       only  : write_rho
      use m_forhar,      only  : forhar
      use xcmod,         only  : nXCfunc, XCauth
      use alloc,         only  : re_alloc, de_alloc
      use siesta_options, only : harrisfun
      use m_dhscf,        only : rhopcc, rhoatm, Vna, uharrs, bcell,
     &                           ntml, dvol, npcc, field, isefld,
     &                           spiral, ntpl, cell, nmpl, shape, nml,
     &                           qspiral
      use moreMeshSubs,   only : setMeshDistr, distMeshData
      use cellxc_mod,     only : GGA
#ifdef MPI
      use mpi_siesta
#endif
      use iogrid_netcdf, only: write_grid_netcdf
      use iogrid_netcdf, only: read_grid_netcdf
      use siesta_options, only: read_charge_cdf
      use siesta_options, only: read_deformation_charge_cdf
      implicit none
      integer
     &  maxnd, maxnh, nua, na, norb, nspin, nuo, nuotot,
     &  iaorb(norb), ifa, iHmat, indxua(na), iphorb(norb), isa(na), 
     &  istr, listd(*), listdptr(nuo), ntm(3), numd(nuo)

      real(dp)
     &  datm(norb), dipol(3), Dscf(maxnd,nspin),
     &  DUext, DUscf, Dxc, Enaatm, Enascf, Exc,
     &  Hmat(maxnh,nspin), Uatm, Uscf, xa(3,na),
     &  stress(3,3), Fal(3,nua), stressl(3,3)

      character
     &  filevh*(*), filevt*(*), fildrh*(*), filrho*(*), 
     &  filepsch*(*), filetoch*(*)

C ----------------------------------------------------------------------
C Routines called internally:
C ----------------------------------------------------------------------
C        cellxc(...)    : Finds total exch-corr energy and potential
C        cross(a,b,c)   : Finds the cross product of two vectors
C        dfscf(...)     : Finds SCF contribution to atomic forces
C        dipole(...)    : Finds electric dipole moment
C        efield(...)    : Adds potential of an external electric field
C        write_rho(...)     : Saves electron density on a file
C        poison(...)    : Solves Poisson equation
C        reord(...)     : Reorders electron density and potential arrays
C        rhooda(...)    : Finds Harris electron density in the mesh
C        rhoofd(...)    : Finds SCF electron density in the mesh
C        rhoofdsp(...)  : Finds SCF electron density in the mesh for
C                         spiral arrangement of spins
C        timer(...)     : Finds CPU times
C        vmat(...)      : Finds matrix elements of SCF potential
C        vmatsp(...)    : Finds matrix elements of SCF potential for
C                         spiral arrangement of spins
C real*8 volcel( cell ) : Returns volume of unit cell
C ----------------------------------------------------------------------
C Internal variables and arrays:
C ----------------------------------------------------------------------
C real*8  bcell(3,3)    : Bulk lattice vectors
C real*8  cell(3,3)     : Auxiliary lattice vectors (same as ucell)
C real*8  const         : Auxiliary variable (constant within a loop)
C real*8  DEc           : Auxiliary variable to call cellxc
C real*8  DEx           : Auxiliary variable to call cellxc
C real*8  dvol          : Mesh-cell volume
C real*8  Ec            : Correlation energy
C real*8  Ex            : Exchange energy
C real*8  field(3)      : External electric field
C integer i             : General-purpose index
C integer ia            : Atom index
C integer io            : Orbital index
C integer ip            : Point index
C integer is            : Species index
C logical IsDiag        : Is supercell diagonal?
C logical IsEfld        : Is there an external electric field?
C integer ispin         : Spin index
C integer j             : General-purpose index
C integer nbcell        : Number of independent bulk lattice vectors
C integer nm(3)         : Number of Mesh divisions of each cell vector
C integer npcc          : Partial core corrections? (0=no, 1=yes)
C integer nsd           : Number of diagonal spin values (1 or 2)
C integer ntpl          : Number of mesh Total Points in unit cell
C                         (including subpoints) locally
C integer ntopl         : Total number of nonzero orbital points locally
C real*4  rhoatm(ntpl)  : Harris electron density
C real*4  rhopcc(ntpl)  : Partial-core-correction density for xc
C real*4  DRho(ntpl)    : Selfconsistent electron density difference
C real*8  rhotot        : Total density at one point
C real*8  rmax          : Maximum orbital radius
C logical samesh        : Same mesh of last call?
C logical samexa        : Same atomic positions of last call?
C real*8  scell(3,3)    : Supercell vectors
C character shape*10    : Name of system shape
C real*4  Vaux(ntpl)    : Auxiliary potential array
C real*4  Vna(ntpl)     : Sum of neutral-atom potentials
C real*8  volume        : Unit cell volume
C real*4  Vscf(ntpl)    : Hartree potential of selfconsistent density
C real*8  x0(3)         : Center of molecule
C logical harrisfun     : Harris functional or Kohn-Sham?
C ----------------------------------------------------------------------

C
C Internal variables
C
      integer  :: i, ia, ip, ispin, nsd,
     &            ntpl_1, ntpl_2, ntpl_3

      real(dp) :: b1Xb2(3), const, DEc, DEx, DStres(3,3),
     $            Ec, Ex, rhotot, x0(3),
     $            volume

#ifdef MPI
      integer  MPIerror
#endif

!
!     Dummy arrays for cellxc call
!
      real(grid_p) ::  aux3(3,1)
      real(grid_p), dimension(1,1,1)  :: dummy_DVxcdn

      real(dp), external :: volcel, ddot

      external
     &  cellxc, cross,
     &  dfscf,  dipole,
     &  efield, poison, 
     &  reord, rhooda, rhoofd, rhoofdsp, 
     &  timer, vmat, vmatsp, 
     &  readsp

      ! Work arrays
      real(grid_p), pointer :: Vscf(:,:), Vscf_seq(:), Vscf_gga(:,:),
     &                         DRho(:,:), DRho_seq(:), DRho_gga(:,:),
     &                         Vaux(:), Chlocal(:), Totchar(:),
     &                         IOBuff(:,:)
      real(dp)              :: sbuffer(7), rbuffer(7)
!------------------------------------------------------------------------- BEGIN
C ----------------------------------------------------------------------
C Start of SCF iteration part
C ----------------------------------------------------------------------
      nullify( Vscf, Vscf_seq, Vscf_gga, DRho, DRho_seq, DRho_gga,
     &         Vaux, Chlocal, Totchar, IOBuff )
      call timer( 'DHSCF', 1 )
      call timer( 'DHSCF3', 1 )

      volume = volcel(cell)
      ntpl_2 = ntpl
C Allocate local memory -------------------------------------------------
      call re_alloc( DRho, 1, ntpl, 1, nspin, 'DRho','dhscf' )

C-------------------------------------------------------------------------


C Find number of diagonal spin values
      nsd = min( nspin, 2 )

C ----------------------------------------------------------------------
C Find SCF electron density at mesh points. Store it in array DRho
C ----------------------------------------------------------------------
#ifdef DEBUG
        call write_debug( '    PRE rhoofd' )
#endif
!
!     We are now in distribution "2". We need to expand the capabilities
!     of the reading routine to handle arbitrary distributions.
!
      if (read_charge_cdf .or. read_deformation_charge_cdf) then
         read_charge_cdf = .false.
         read_deformation_charge_cdf = .false.
         call die("Reading of Charge for SCF not implemented yet")
      endif
      if (read_charge_cdf) then
         call read_grid_netcdf(ntm(1:3),nspin,ntpl,DRho,"Rho")
         do ispin = 1,nspin
            call reord( DRho(1,ispin), DRho(1,ispin), nml, nsm, -1)
         enddo
         read_charge_cdf = .false.
      else if (read_deformation_charge_cdf) then
         call read_grid_netcdf(ntm(1:3),nspin,ntpl,DRho,"DeltaRho")
         do ispin = 1,nspin
            call reord( DRho(1,ispin), DRho(1,ispin), nml, nsm, -1)
            DRho(:,ispin) = DRho(:,ispin) + rhoatm(:)
         enddo
         read_deformation_charge_cdf = .false.
      else
         if (spiral) then
            call rhoofdsp(norb, nmpl, maxnd, numd, listdptr, listd,
     &                nspin, Dscf, DRho, nuo, nuotot, iaorb, iphorb,
     &                isa, qspiral)
         else
            call rhoofd( norb, nmpl, maxnd, numd, listdptr, listd,
     &               nspin,
     &               Dscf, DRho, nuo, nuotot, iaorb, iphorb, isa)
         endif
      endif

C ----------------------------------------------------------------------
C Save electron density
C ----------------------------------------------------------------------
      if (filrho .ne. ' ') then
#ifdef MPI
        call setMeshDistr( 1, nsm, nsp, nml, nmpl, ntml, ntpl )
#endif
        ntpl_1 = ntpl

        call re_alloc( IOBuff, 1, ntpl, 1, nspin, name='IOBuff',
     &               routine='dhscf' )
        do ispin = 1,nspin
          call distMeshData( 2, DRho(:,ispin),
     &                       1, IOBuff(:,ispin), +1 )
        enddo

        call write_rho( filrho, cell, ntm, nsm, ntpl, nspin, IOBuff )
        call write_grid_netcdf( cell, ntm, nspin, ntpl, IOBuff, "Rho")

        call de_alloc( IOBuff, name='IOBuff', routine='dhscf' )
#ifdef MPI
        call setMeshDistr( 2, nsm, nsp, nml, nmpl, ntml, ntpl )
#endif
      endif
C ----------------------------------------------------------------------
C Save the diffuse ionic charge and/or the total (ionic+electronic) charge
C ----------------------------------------------------------------------
      if (filepsch .ne. ' ' .or. filetoch .ne. ' ') then
        call re_alloc( Chlocal, 1, ntpl, name='Chlocal',
     &                 routine='dhscf' )

C Find diffuse ionic charge on mesh
        call LocalChargeOnMesh( na, isa, ntpl, Chlocal, indxua )
      
C Save diffuse ionic charge 
        if (filepsch .ne. ' ') then
#ifdef MPI
          call setMeshDistr( 1, nsm, nsp, nml, nmpl, ntml, ntpl )
#endif
          call re_alloc( IOBuff, 1, ntpl, 1, 1, name='IOBuff',
     &                   routine='dhscf' )
          call distMeshData( 2, Chlocal(:),
     &                       1, IOBuff(:,1), +1 )

          call write_rho( filepsch, cell, ntm, nsm, ntpl, 1, IOBuff)
          call write_grid_netcdf( cell, ntm, 1, ntpl,
     $                            IOBuff, "Chlocal")
        call de_alloc( IOBuff, name='IOBuff', routine='dhscf' )
#ifdef MPI
        call setMeshDistr( 2, nsm, nsp, nml, nmpl, ntml, ntpl )
#endif
        endif

C Save total (ionic+electronic) charge 
        if ( filetoch .ne. ' ') then
           call re_alloc( Totchar, 1, ntpl, name='Totchar',
     &          routine='dhscf' )
        
          Totchar(:) = 0.0_grid_p
          do ispin = 1,nsd
            Totchar(1:ntpl) = Totchar(1:ntpl) + 
     .                        Chlocal(1:ntpl)/nsd + DRho(1:ntpl,ispin)
          enddo
#ifdef MPI
          call setMeshDistr( 1, nsm, nsp, nml, nmpl, ntml, ntpl )
#endif
          call re_alloc( IOBuff, 1, ntpl, 1, 1, name='IOBuff',
     &                   routine='dhscf' )
          call distMeshData( 2, Totchar(:),
     &                       1, IOBuff(:,1), +1 )
          call write_rho( filetoch, cell, ntm, nsm, ntpl, 1, IOBuff )
          call write_grid_netcdf( cell, ntm, 1, ntpl,
     $                            IOBuff, "TotalCharge")
          call de_alloc( IOBuff, name='IOBuff', routine='dhscf' )
#ifdef MPI
          call setMeshDistr( 2, nsm, nsp, nml, nmpl, ntml, ntpl )
#endif
          call de_alloc( Totchar, name='Totchar', routine='dhscf' )
        endif 

        call de_alloc( Chlocal, name='Chlocal', routine='dhscf' )
      endif

C Find difference between selfconsistent and atomic densities

      do ispin = 1,nsd
        do ip = 1,ntpl
          DRho(ip,ispin) = DRho(ip,ispin) -
     &                         rhoatm(ip)/nsd
        enddo
      enddo

C ----------------------------------------------------------------------
C Save electron density difference
C ----------------------------------------------------------------------
      if (fildrh .ne. ' ') then
#ifdef MPI
        call setMeshDistr( 1, nsm, nsp, nml, nmpl, ntml, ntpl )
#endif
        call re_alloc( IOBuff, 1, ntpl, 1, nspin, name='IOBuff',
     &               routine='dhscf' )
        do ispin = 1,nspin
          call distMeshData( 2, DRho(:,ispin),
     &                       1, IOBuff(:,ispin), +1 )
        enddo

        call write_rho( fildrh, cell, ntm, nsm, ntpl, nspin, IOBuff )
        call write_grid_netcdf( cell, ntm, nspin, ntpl,
     $                          IOBuff, "DeltaRho")
        call de_alloc( IOBuff, name='IOBuff', routine='dhscf' )
#ifdef MPI
        call setMeshDistr( 2, nsm, nsp, nml, nmpl, ntml, ntpl )
#endif
      endif

C Transform spin density into sum and difference
      if (nsd .eq. 2) then
        do ip = 1,ntpl
          rhotot     = DRho(ip,1) + DRho(ip,2)
          DRho(ip,2) = DRho(ip,2) - DRho(ip,1)
          DRho(ip,1) = rhotot
        enddo
      endif
C ----------------------------------------------------------------------
C Calculate the dipole moment
C ----------------------------------------------------------------------
      dipol(1:3) = 0.0_dp
      if (shape .ne. 'bulk') then

C Find center of system
        x0(1:3) = 0.0_dp
        do ia = 1,nua
          x0(1:3) = x0(1:3) + xa(1:3,ia) / nua
        enddo

C Find dipole
        !! AG** Is this correct? Or do we need an extra transformation?
        call reord( DRho, DRho, nml, nsm, +1)
        call dipole( cell, ntm, ntml(1), ntml(2), ntml(3), nsm,
     &    DRho, x0, dipol )
        call reord( DRho, DRho, nml, nsm, -1)

C Orthogonalize dipole to bulk directions
        if (shape .eq. 'chain') then
          const = ddot(3,dipol,1,bcell,1) / ddot(3,bcell,1,bcell,1)
          dipol(1:3) = dipol(1:3) - const * bcell(1:3,1)
        elseif (shape .eq. 'slab') then
          call cross( bcell(1,1), bcell(1,2), b1Xb2 )
          const = ddot(3,dipol,1,b1Xb2,1) / ddot(3,b1Xb2,1,b1Xb2,1)
          dipol(1:3) = const * b1Xb2(1:3)
        endif
      endif
C ----------------------------------------------------------------------
C     Set the uniform data distribution
C ----------------------------------------------------------------------
#ifdef MPI
      call setMeshDistr( 1, nsm, nsp, nml, nmpl, ntml, ntpl )
#endif
      ntpl_1 = ntpl

      call re_alloc( DRho_seq, 1, ntpl, name='DRho_seq',
     &               routine='dhscf' )
      call distMeshData( 2, DRho(:,1),
     &                   1, DRho_seq, +1 )

C ----------------------------------------------------------------------
C     Find Hartree potential of DRho = rhoscf-rhoatm. Store it in Vscf
C ----------------------------------------------------------------------
C     Solve Poisson's equation
      call re_alloc( Vscf_seq, 1, ntpl, name='Vscf_seq',
     &               routine='dhscf' )
#ifdef DEBUG
        call write_debug( '    PRE Poison' )
#endif
      call poison( cell, ntml(1), ntml(2), ntml(3), ntm, DRho_seq,
     &             DUscf, Vscf_seq, DStres, nsm )
      call de_alloc( DRho_seq,   name='DRho_seq', routine='dhscf' )

C Add contribution to stress from electrostatic energy of rhoscf-rhoatm
      if (istr .eq. 1) then
        stressl(1:3,1:3) = stressl(1:3,1:3) + DStres(1:3,1:3)
      endif

C ----------------------------------------------------------------------
C     Move data to the nested disection distribution
C ----------------------------------------------------------------------
#ifdef MPI
      call setMeshDistr( 2, nsm, nsp, nml, nmpl, ntml, ntpl )
#endif
      call re_alloc( Vaux, 1, ntpl, name='Vaux', routine='dhscf' )

      call distMeshData( 1, Vscf_seq,
     &                   2, Vaux, -1 )
      call de_alloc( Vscf_seq, name='Vscf_seq', routine='dhscf' )

C ----------------------------------------------------------------------
C     Find electrostatic (Hartree) energy of full SCF electron density
C     using the original data distribution 
C ----------------------------------------------------------------------
      Uatm = Uharrs
      Uscf = Uatm + DUscf
      do ip = 1,ntpl
        Uscf = Uscf + dvol * Vaux(ip) * rhoatm(ip)
      enddo

C ----------------------------------------------------------------------
C Add neutral-atom potential to Vscf
C ----------------------------------------------------------------------
      Enaatm = 0.0_dp
      Enascf = 0.0_dp
      do ip = 1,ntpl
        Enaatm         = Enaatm + dvol*Vna(ip)*rhoatm(ip)
        Enascf         = Enascf + dvol*Vna(ip)*
     &                   (rhoatm(ip) + DRho(ip,1))
        Vaux(ip)       = Vaux(ip) + Vna(ip)
      enddo

C ----------------------------------------------------------------------
C Add potential from external electric field (if present)
C ----------------------------------------------------------------------
      !! ** AG: Is the Efield treatment correct?
      if (IsEfld) then
        call reord( Vaux, Vaux, nml, nsm, +1)
        call efield( cell, nua, isa, xa, ntm, nsm, Vaux, field )
        call reord( Vaux, Vaux, nml, nsm, -1)
        DUext = - ddot(3,field,1,dipol,1)
      endif

C ----------------------------------------------------------------------
C Save electrostatic potential
C ----------------------------------------------------------------------

      if (filevh .ne. ' ') then
#ifdef MPI
        call setMeshDistr( 1, nsm, nsp, nml, nmpl, ntml, ntpl )
#endif
        call re_alloc( IOBuff, 1, ntpl, 1, 1, name='IOBuff',
     &               routine='dhscf' )
        call distMeshData( 2, Vaux(:),
     &                     1, IOBuff(:,1), +1 )

        call write_rho( filevh, cell, ntm, nsm, ntpl, 1, IOBuff)
        call write_grid_netcdf( cell, ntm, 1, ntpl,
     $                          IOBuff, "ElectrostaticPotential")
        call de_alloc( IOBuff, name='IOBuff', routine='dhscf' )
#ifdef MPI
        call setMeshDistr( 2, nsm, nsp, nml, nmpl, ntml, ntpl )
#endif
      endif

C Get back spin density from sum and difference
      if (nsd .eq. 2) then
        do ip = 1,ntpl
          rhotot         = DRho(ip,1)
          DRho(ip,1) = 0.5_dp*(rhotot - DRho(ip,2))
          DRho(ip,2) = 0.5_dp*(rhotot + DRho(ip,2))
        enddo
      endif

C ----------------------------------------------------------------------
C Exchange-correlation energy
C ----------------------------------------------------------------------
      do ispin = 1,nsd
        DRho(1:ntpl,ispin) = DRho(1:ntpl,ispin) +
     &                       rhoatm(1:ntpl) / nsd
        if (npcc .eq. 1) then
          DRho(1:ntpl,ispin) = DRho(1:ntpl,ispin) +
     &                         rhopcc(1:ntpl)/nsd
        endif
      enddo

#ifdef MPI
      call setMeshDistr( 3, nsm, nsp, nml, nmpl, ntml, ntpl )
      ntpl_3 = ntpl
      call re_alloc( Vscf_gga, 1, ntpl, 1, nspin, name='Vscf_gga',
     &               routine='dhscf' )
      call re_alloc( DRho_gga, 1, ntpl, 1, nspin, name='DRho_gga',
     &               routine='dhscf' )

      do ispin = 1,nsd
        call distMeshData( 2, DRho(:,ispin),
     &                     3, DRho_gga(:,ispin), +1 )
      enddo
#else
      DRho_gga => DRho
      call re_alloc( Vscf_gga, 1, ntpl, 1, nspin, name='Vscf_gga',
     &               routine='dhscf' )      
#endif

#ifdef DEBUG
        call write_debug( '    PRE cellxc' )
#endif
      call cellxc( 0, 0, cell, ntml, ntml, ntpl, 0, aux3, nspin,
     &             DRho_gga, Ex, Ec, DEx, DEc, Vscf_gga,
     &             dummy_DVxcdn, stressl )
#ifdef MPI
      call setMeshDistr( 2, nsm, nsp, nml, nmpl, ntml, ntpl )

      call re_alloc( Vscf, 1, ntpl, 1, nspin, name='Vscf',
     &               routine='dhscf' )

      do ispin = 1,nsd
        call distMeshData( 3, Vscf_gga(:,ispin),
     &                     2, Vscf(:,ispin), -1 )
      enddo
      call de_alloc( DRho_gga,   name='DRho_gga', routine='dhscf' )
      call de_alloc( Vscf_gga,   name='Vscf_gga', routine='dhscf' )
#else
      Vscf => Vscf_gga
#endif
      Exc = Ex + Ec
      Dxc = DEx + DEc
      do ispin = 1,nsd
        do ip = 1,ntpl
          DRho(ip,ispin) = DRho(ip,ispin) - rhoatm(ip) / nsd
          if (npcc .eq. 1) then
            DRho(ip,ispin) = DRho(ip,ispin) - rhopcc(ip)/nsd
          endif
          Vscf(ip,ispin) = Vscf(ip,ispin) + Vaux(ip)
        enddo
      enddo

C ----------------------------------------------------------------------
C     Save total potential
C ----------------------------------------------------------------------
      if (filevt .ne. ' ') then
#ifdef MPI
        call setMeshDistr( 1, nsm, nsp, nml, nmpl, ntml, ntpl )
#endif
        call re_alloc( IOBuff, 1, ntpl, 1, nspin, name='IOBuff',
     &               routine='dhscf' )
        do ispin = 1,nspin
          call distMeshData( 2, Vscf(:,ispin),
     &                       1, IOBuff(:,ispin), +1 )
        enddo

        call write_rho( filevt, cell, ntm, nsm, ntpl, nspin, IOBuff )
        call write_grid_netcdf( cell, ntm, nspin, ntpl,
     $                          IOBuff, "TotalPotential")
        call de_alloc( IOBuff, name='IOBuff', routine='dhscf' )
#ifdef MPI
        call setMeshDistr( 2, nsm, nsp, nml, nmpl, ntml, ntpl )
#endif
      endif

C ----------------------------------------------------------------------
C     Find SCF contribution to hamiltonian matrix elements
C ----------------------------------------------------------------------
#ifdef DEBUG
        call write_debug( '    PRE vmat' )
#endif
      if (iHmat .eq. 1) then
        if (spiral) then
          call vmatsp( norb, nmpl, dvol, nspin, Vscf, maxnd,
     &                 numd, listdptr, listd, Hmat, nuo,
     &                 nuotot, iaorb, iphorb, isa, qspiral )
        else
          call vmat( norb, nmpl, dvol, nspin, Vscf, maxnd,
     &               numd, listdptr, listd, Hmat, nuo,
     &               nuotot, iaorb, iphorb, isa )
        endif
      endif
#ifdef MPI
C Global reduction of Uscf/DUscf/Uatm/Enaatm/Enascf
      sbuffer(1) = Uscf
      sbuffer(2) = DUscf
      sbuffer(3) = Uatm
      sbuffer(4) = Enaatm
      sbuffer(5) = Enascf
      sbuffer(6) = Exc
      sbuffer(7) = Dxc

#ifdef DEBUG
        call write_debug( '    PRE reduce Energies' )
#endif
      call MPI_AllReduce( sbuffer, rbuffer, 7, MPI_double_precision,
     &                     MPI_Sum, MPI_Comm_World, MPIerror )

      Uscf   = rbuffer(1)
      DUscf  = rbuffer(2)
      Uatm   = rbuffer(3)
      Enaatm = rbuffer(4)
      Enascf = rbuffer(5)
      Exc    = rbuffer(6)
      Dxc    = rbuffer(7)
#endif
C Add contribution to stress from the derivative of the Jacobian of ---
C r->r' (strained r) in the integral of Vna*(rhoscf-rhoatm)
      if (istr .eq. 1) then
        do i = 1,3
          stress(i,i) = stress(i,i) + ( Enascf - Enaatm ) / volume
        enddo
      endif

C Stop time counter for SCF iteration part
      call timer( 'DHSCF3', 2 )

C ----------------------------------------------------------------------
C End of SCF iteration part
C ----------------------------------------------------------------------

      if (ifa.eq.1 .or. istr.eq.1) then
C ----------------------------------------------------------------------
C Forces and stress : SCF contribution
C ----------------------------------------------------------------------
C Start time counter for force calculation part
        call timer( 'DHSCF4', 1 )

C Find contribution of partial-core-correction
        if (npcc .eq. 1) then
          call PartialCoreOnMesh( na, isa, ntpl, rhopcc, indxua,
     &      nsd, dvol, volume, Vscf, Vaux, Fal, stressl, ifa.ne.0,
     &      istr.ne.0 )
        endif

        if ( harrisfun)
     &    call forhar( ntpl, nspin, nml, ntml, npcc, cell, 
     &                 rhoatm, rhopcc, Vna, DRho, Vscf, Vaux )

C Transform spin density into sum and difference
        if (nsd .eq. 2) then
          do ip = 1,ntpl
            rhotot     = DRho(ip,1) + DRho(ip,2)
            DRho(ip,2) = DRho(ip,2) - DRho(ip,1)
            DRho(ip,1) = rhotot
          enddo
        endif

C Find contribution of neutral-atom potential
        call NeutralAtomOnMesh( na, isa, ntpl, Vna, indxua, dvol, 
     &                          volume, DRho, Fal, stressl, 
     &                          ifa.ne.0, istr.ne.0 )

C Vaux is (minus) the potential which multiplies rhoatm
*       if (nsd .eq. 2) then
*         Vaux(1:ntpl) = 0.5_grid_p * Vaux(1:ntpl)
*       endif
        call dfscf( ifa, istr, na, norb, nuo, nuotot, nmpl, nspin,
     &              indxua, isa, iaorb, iphorb,
     &              maxnd, numd, listdptr, listd, Dscf, datm,
     &              Vscf, Vaux, dvol, volume, Fal, stressl )

C Stop time counter for force calculation part
        call timer( 'DHSCF4', 2 )
C ----------------------------------------------------------------------
C End of force and stress calculation
C ----------------------------------------------------------------------
      endif


C Stop time counter
      call timer( 'DHSCF', 2 )

C ----------------------------------------------------------------------
C Free locally allocated memory
C ----------------------------------------------------------------------
      call de_alloc( Vaux, name='Vaux', routine='dhscf' )
      call de_alloc( Vscf, name='Vscf', routine='dhscf' )
      call de_alloc( DRho, name='DRho', routine='dhscf' )
!--------------------------------------------------------------------------- END
      end subroutine dhscf
