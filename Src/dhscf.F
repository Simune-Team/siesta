      subroutine dhscf( nspin, norb, iaorb, iphorb, nuo,
     .                  nuotot, nua, na, isa, xa, indxua, ucell, mscell, 
     .                  G2max, ntm, ifa, istr, iHmat, filrho, fildrh, 
     .                  filevh, filevt, filepsch, filetoch, maxnd, numd,
     .                  listdptr, listd, 
     .                  Dscf, datm, maxnh, numh, listhptr, listh, Hmat,
     .                  Enaatm, Enascf, Uatm, Uscf, DUscf, DUext, 
     .                  Exc, Dxc, dipol, Fa, stress, fal, stressl )

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
C integer nuo           : Number of orbitals in a unit cell
C integer nua           : Number of atoms in unit cell
C integer na            : Number of atoms in supercell
C integer isa(na)       : Species index of all atoms in supercell
C real*8  xa(3,na)      : Atomic positions of all atoms in supercell
C integer indxua        : Index of equivalent atom in unit cell
C real*8  ucell(3,3)    : Unit cell lattice vectors: ucell(ixyz,ivect)
C integer mscell(3,3)   : Supercell vectors in units of ucell:
C                         Supercell(ix,iv) = Sum_jv( ucell(ix,jv) *
C                                                   mscell(jv,iv) )
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
C integer numh(nuo)         : Number of nonzero Hamiltonian-matrix
C                             elements for each matrix row
C integer listhptr(nuo)     : Pointer to start of rows of H-matrix
C integer listh(maxnh)      : Nonzero-Hamiltonian-matrix-element column
C                             indexes for each matrix row
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
C real*8  G2max      : Effective planewave cutoff in Ry. determines
C                        the mesh density and the precision of integrals
C                          on input : value required
C                          on output: value used, which may be larger
C real*8  Fa(3,nua)  : Atomic forces, to which the SCF contribution
C                        is added by this routine when ifa=1.
C                        the SCF contribution is minus the derivative
C                        of ( Enascf - Enaatm + DUscf + Exc ) with
C                        respect to atomic positions, in Ry/Bohr
C real*8 stress(3,3) : Stress tensor, to which the SCF contribution
C                      is added by this routine when ifa=1.
C                      the SCF contribution is minus the derivative of
C                         ( Enascf - Enaatm + DUscf + Exc ) / volume
C                      with respect to the strain tensor, in Ry.
C real*8  fal(3,nua) : Atomic forces, to which the SCF contribution
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

C
C  Modules
C
      use precision
      use parallel
      use atmfuncs, only: rcut, rcore
      use fdf
      use sys
      use mesh
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

      integer
     .  maxnd, maxnh, nua, na, norb, nspin, nuo, nuotot,
     .  iaorb(norb), ifa, iHmat, indxua(na), iphorb(norb), isa(na), 
     .  istr, listd(*), listh(maxnh), listdptr(nuo), listhptr(nuo),
     .  mscell(3,3), ntm(3), numd(nuo), numh(nuo)

      real*8
     .  datm(norb), dipol(3), Dscf(maxnd,nspin),
     .  DUext, DUscf, Dxc,
     .  Enaatm, Enascf, Exc, Fa(3,nua), G2max, G2mesh,
     .  Hmat(maxnh,nspin),
     .  stress(3,3), Uatm, ucell(3,3), Uscf, xa(3,na),
     .  fal(3,nua), stressl(3,3)

      character
     .  filevh*(*), filevt*(*), fildrh*(*), filrho*(*), 
     .  filepsch*(*), filetoch*(*)

C ----------------------------------------------------------------------
C Routines called internally:
C ----------------------------------------------------------------------
C        cellxc(...)    : Finds total exch-corr energy and potential
C        cross(a,b,c)   : Finds the cross product of two vectors
C        dfscf(...)     : Finds SCF contribution to atomic forces
C        dipole(...)    : Finds electric dipole moment
C real*8 dot(v1,v2,n)   : Finds the dot product of vectors V1 and V2
C                         of length N.
C        efield(...)    : Adds potential of an external electric field
C        iorho(...)     : Saves electron density on a file
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
C real*8  DFa(3)        : Auxiliary variable to call dfscf
C real*8  DFaatm(3)     : Auxiliary variable to call dfscf
C real*8  DFascf(3)     : Auxiliary variable to call dfscf
C real*8  dvol          : Mesh-cell volume
C real*8  Ec            : Correlation energy
C real*8  Ex            : Exchange energy
C real*8  field(3)      : External electric field
C logical frstme        : First time this routine is called?
C logical found         : File found?
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
C integer nzero(3)      : Auxiliary array with zeros, to call efield
C real*4  rhoatm(ntpl)  : Harris electron density
C real*4  rhopcc(ntpl)  : Partial-core-correction density for xc
C real*4  DRho(ntpl)    : Selfconsistent electron density difference
C real*4  DVxcdn(ntpl,nspin,nspin): Derivatice of exchange-correlation
C                         potential respect charge density (see forhar.f)
C                         DVxcdn(ip,spin1,spin2) = D Vxc(spin2) / D n(spin1)
C real*8  rhotot        : Total density at one point
C real*8  rmax          : Maximum orbital radius
C logical samesh        : Same mesh of last call?
C logical samexa        : Same atomic positions of last call?
C real*8  scell(3,3)    : Supercell vectors
C character shape*10    : Name of system shape
C real*8  tiny          : A small constant
C real*8  Uharrs        : Hartree energy of Harris electron density
C real*4  Vaux(ntpl)    : Auxiliary potential array
C real*4  Vna(ntpl)     : Sum of neutral-atom potentials
C real*8  volume        : Unit cell volume
C real*4  Vscf(ntpl)    : Hartree potential of selfconsistent density
C character XCauth*10   : Initials of xc version authors
C character XCfunc*10   : Name of xc functional
C real*8  x0(3)         : Center of molecule
C logical harrisfun     : Harris functional or Kohn-Sham?
C ----------------------------------------------------------------------

C
C Internal variables
C
      integer
     .  i, ia, io, ip, iphi, is, ispin, j, ifor, iu, 
     .  nbcell, nm(3), npcc, nsd, nzero(3), 
     .  Node, Nodes, ntopl, nmpl, nsc(3), ntpl, 
     .  nml(3), ntml(3)

#ifdef MPI
      integer
     .  MPIerror
#endif

      real, dimension(:), allocatable, save ::
     .  Vaux

      real, dimension(:), allocatable, save ::
     .  rhopcc, rhoatm, Vna

      real, dimension(:,:), allocatable, save ::
     .  DRho, Vscf

      real, dimension(:,:,:), allocatable, save ::
     .  DVxcdn

      real, dimension(:), allocatable::
     .  Chlocal, Totchar
!
!     Dummy arrays for cellxc call
!
      real aux1(1), aux3(3,1)

      real*8
     .  bcell(3,3), b1Xb2(3), cell(3,3), const,
     .  DEc, DEx, DFa(3), DFaatm(3), DFascf(3),
     .  dot, DStres(3,3), dvol, Ec, Ex, field(3),
     .  rhotot, rmax, scell(3,3), tiny,
     .  Uharrs, volume, volcel, x0(3), qspiral(3)

#ifdef MPI
      real*8
     .  Eloc
#endif

      logical
     .  found, frstme, IsDiag, IsEfld, samesh, samexa,
     .  harrisfun, leqi, spiral

      character
     .  shape*10, XCauth*10, XCfunc*10

      external
     .  cellxc, cross,
     .  dfscf,  dipole, dot,  
     .  efield, iorho, poison, 
     .  reord, rhooda, rhoofd, rhoofdsp, 
     .  timer, vmat, vmatsp, volcel,
     .  memory, readsp

      save
     .  bcell, cell, dvol, field, frstme, G2mesh,
     .  IsDiag, IsEfld, nm, nml, ntml,
     .  npcc, nzero, nmpl, ntpl, ntopl,
     .  rmax, samesh, samexa, scell, shape, 
     .  tiny, Uharrs, volume, XCauth, XCfunc, harrisfun

      data
     .  G2mesh / 0.0d0 /
     .  nzero / 3*0 /,
     .  tiny  / 1.D-12 /,
     .  frstme /.true./

C ----------------------------------------------------------------------

C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
#else
      Node = 0
      Nodes = 1
#endif

C ----------------------------------------------------------------------
C General initialisation
C ----------------------------------------------------------------------

C Start time counter
      call timer( 'DHSCF', 1 )

      if (frstme) then
        if (Node.eq.0) then
          XCfunc = fdf_string('xc.functional','LDA')
          XCauth = fdf_string('xc.authors','PZ')
          harrisfun = fdf_boolean('Harris_functional',.false.)
        endif
#ifdef MPI
        call MPI_Bcast(XCfunc,10,MPI_character,0,MPI_Comm_World,
     *    MPIerror)
        call MPI_Bcast(XCauth,10,MPI_character,0,MPI_Comm_World,
     *    MPIerror)
        call MPI_Bcast(harrisfun,1,MPI_logical,0,MPI_Comm_World,
     *    MPIerror)
#endif
      endif

      if( .not. harrisfun) then
        ifor = 0
      else if (harrisfun) then
        ifor = 1
        if (.not. ( leqi(XCauth,'PZ') .or. leqi(XCauth,'CA') ) ) then
           if (Node .eq. 0) then
             write(6,100)
             write(6,101)
             write(6,'(a)')
     .        'dhscf:   Harris forces for a XC functional different of'      
             write(6,'(a)')
     .        'dhscf:   Perdew-Zunger has not been implemented yet'
             write(6,102)
           endif
           call die
        endif
      endif

C Read spin-spiral wavevector (if defined)
      call readsp( qspiral, spiral )
      if (spiral .and. nspin.lt.4) then
        if (Node.eq.0) then
        write(6,*) 'dhscf: ERROR: spiral is only defined for nspin = 4'
        endif
        stop
      endif

C ----------------------------------------------------------------------
C Orbital initialisation : part 1
C ----------------------------------------------------------------------

C Find the maximum orbital radius
      rmax = 0.0d0
      do io = 1,norb
        ia = iaorb(io)
        iphi = iphorb(io)
        is = isa(ia)
        rmax = max( rmax, rcut(is,iphi) )
      enddo

C ----------------------------------------------------------------------
C Check whether mesh initialisation is needed
C ----------------------------------------------------------------------
      call SameMeshAndAtoms( na, xa, ucell, rmax, G2max, G2mesh,
     .  samesh, samexa )

C ----------------------------------------------------------------------
C Initialization due to change of mesh
C ----------------------------------------------------------------------
      if ( .not.samesh ) then

C Start time counter for mesh initialization
        call timer( 'DHSCF1', 1 )


C ----------------------------------------------------------------------
C Unit cell handling
C ----------------------------------------------------------------------
C Find diagonal unit cell and supercell
        call digcel( ucell, mscell, cell, scell, nsc, IsDiag )
        if (.not.IsDiag) then
          if (Node.eq.0) then
            write(6,'(/,a,3(/,a,3f12.6,a,i6))')
     .      'DHSCF: WARNING: New shape of unit cell and supercell:',
     .     ('DHSCF:',(cell(i,j),i=1,3),'   x',nsc(j),j=1,3)
          endif
        endif

C Find the system shape
        call shaper( cell, nua, isa, xa, shape, nbcell, bcell )

C Find system volume
        volume = volcel( cell )

C ----------------------------------------------------------------------
C Mesh initialization 
C ----------------------------------------------------------------------
        call InitMesh( na, xa, cell, norb, iaorb, iphorb, isa, rmax, 
     .                 G2max, G2mesh, Node, Nodes, nsc, nmpl, nm,
     .                 nml, ntm, ntml, ntpl, ntopl, dvol)

C ----------------------------------------------------------------------
C Allocate memory that depends on the mesh
C ----------------------------------------------------------------------

C Allocate memory related to ntpl
        if (allocated(rhoatm)) then
          call memory('D','S',size(rhoatm),'dhscf')
          deallocate(rhoatm)
        endif
        if (allocated(Vna)) then
          call memory('D','S',size(Vna),'dhscf')
          deallocate(Vna)
        endif
        allocate(rhoatm(ntpl))
        call memory('A','S',ntpl,'dhscf')
        allocate(Vna(ntpl))
        call memory('A','S',ntpl,'dhscf')

C Find if there are partial-core-corrections for any atom
        npcc = 0
        do ia = 1,na
          if (rcore(isa(ia)) .gt. tiny) npcc = 1
        enddo

C Allocate memory for PCC density
        if (allocated(rhopcc)) then
          call memory('D','S',size(rhopcc),'dhscf')
          deallocate(rhopcc)
        endif
        allocate(rhopcc(ntpl*npcc+1))
        call memory('A','S',ntpl*npcc,'dhscf')

C Stop time counter for mesh initialization
        call timer( 'DHSCF1', 2 )

      endif
C ----------------------------------------------------------------------
C End of mesh initialization 
C ----------------------------------------------------------------------

C Allocate local memory
      allocate(DRho(ntpl,nspin))
      call memory('A','S',ntpl*nspin,'dhscf')
      allocate(Vscf(ntpl,nspin))
      call memory('A','S',ntpl*nspin,'dhscf')
      allocate(Vaux(ntpl))
      call memory('A','S',ntpl,'dhscf')
      allocate(DVxcdn(ntpl,nspin,nspin))
      call memory('A','S',1*nspin*nspin,'dhscf')

C Initialization of DVxcdn
      DVxcdn(1:ntpl,1:nspin,1:nspin) = 0.0



      if (.not.samexa) then
C ----------------------------------------------------------------------
C Initialize atomic orbitals, density and potential
C ----------------------------------------------------------------------

C Start time counter for atomic initializations
        call timer( 'DHSCF2', 1 )

C Initialise quantities relating to the atom-mesh positioning
        call InitAtomMesh( na, xa)

C Find partial-core-correction energy density
        if (npcc .eq. 1) then
          call PartialCoreOnMesh( na, isa, ntpl, rhopcc, indxua,
     .      nsd, dvol, volume, Vscf, Vaux, fal, stressl, .false.,
     .      .false. )
        endif

C Find neutral-atom potential
        call NeutralAtomOnMesh( na, isa, ntpl, Vna, indxua, dvol, 
     .                          volume, DRho, fal, stressl, 
     .                          .false., .false. )

C Calculate orbital values on mesh
        call PhiOnMesh( nmpl, ntopl, norb, iaorb, iphorb, isa, Node)

C ----------------------------------------------------------------------
C Create sparse indexing for Dscf as needed for local mesh
C ----------------------------------------------------------------------

        if (nuo .ne. nuotot) then
          call CreateLocalDscfPointers( nmpl, nuotot, numd, listdptr, 
     .      listd, Node, Nodes)
        endif

C ----------------------------------------------------------------------
C Calculate terms relating to the neutral atoms on the mesh
C ----------------------------------------------------------------------

C Find Harris (sum of atomic) electron density
        call rhooda( norb, nmpl, datm, rhoatm, iaorb, iphorb, isa )

C Solve Poisson's equation to find Hartree energy of rhoatm
        call reord( rhoatm, rhoatm, nml, nsm, +1)
        call poison( cell, ntml(1), ntml(2), ntml(3), ntm, rhoatm,
     .               Uharrs, Vscf, DStres, nsm )
        call reord( rhoatm, rhoatm, nml, nsm, -1)

C Stop time counter for atomic initializations
        call timer( 'DHSCF2', 2 )

C ----------------------------------------------------------------------
C End of atomic initializations
C ----------------------------------------------------------------------
      endif

C ----------------------------------------------------------------------
C Start of SCF iteration part
C ----------------------------------------------------------------------
      call timer( 'DHSCF3', 1 )

C Find number of diagonal spin values
      nsd = min( nspin, 2 )

C ----------------------------------------------------------------------
C Find SCF electron density at mesh points. Store it in array DRho
C ----------------------------------------------------------------------

      if (spiral) then
        call rhoofdsp(norb, nmpl, maxnd, numd, listdptr, listd,
     .                nspin, Dscf, DRho, nuo, nuotot, iaorb, iphorb,
     .                isa, qspiral)
      else
        call rhoofd(norb, nmpl, maxnd, numd, listdptr, listd,
     .              nspin, Dscf, DRho, nuo, nuotot, iaorb, iphorb, isa)
      endif

C ----------------------------------------------------------------------
C Save electron density
C ----------------------------------------------------------------------
      if (filrho .ne. ' ') then
        do ispin = 1,nspin
          call reord( DRho(1,ispin), DRho(1,ispin), nml, nsm, +1)
        enddo
        call iorho( 'write', filrho, cell, ntm, nsm, ntpl, nspin,
     .              DRho, found )
        do ispin = 1,nspin
          call reord( DRho(1,ispin), DRho(1,ispin), nml, nsm, -1)
        enddo
      endif

C ----------------------------------------------------------------------
C Save the diffuse ionic charge and/or the total (ionic+electronic) charge
C ----------------------------------------------------------------------

      if (filepsch .ne. ' ' .or. filetoch .ne. ' ') then

        allocate(Chlocal(ntpl))
        call memory('A','S',ntpl,'dhscf')

C Find diffuse ionic charge on mesh
        call LocalChargeOnMesh( na, isa, ntpl, Chlocal, indxua )
      
C Save diffuse ionic charge 
        if (filepsch .ne. ' ') then
          call reord( Chlocal, Chlocal, nml, nsm, +1)
          call iorho( 'write', filepsch, cell, ntm, nsm, ntpl, 1, 
     .      Chlocal, found )
          call reord( Chlocal, Chlocal, nml, nsm, -1)
        endif

C Save total (ionic+electronic) charge 
        if ( filetoch .ne. ' ') then
          allocate(Totchar(ntpl))
          call memory('A','S',ntpl,'dhscf')
        
          Totchar(:) = 0.0
          do ispin = 1,nsd
            Totchar(1:ntpl) = Totchar(1:ntpl) + 
     .                        Chlocal(1:ntpl)/nsd + DRho(1:ntpl,ispin)
          enddo
          call reord( Totchar, Totchar, nml, nsm, +1)
          call iorho( 'write', filetoch, cell, ntm, nsm, ntpl, 1, 
     .      Totchar, found )
          call reord( Totchar, Totchar, nml, nsm, -1)

          call memory('D','S',ntpl,'dhscf')
          deallocate(Totchar)
        endif 

        call memory('D','S',ntpl,'dhscf')
        deallocate(Chlocal)
      endif

C Find difference between selfconsistent and atomic densities
      do ispin = 1,nsd
        DRho(1:ntpl,ispin) = DRho(1:ntpl,ispin) - rhoatm(1:ntpl)/nsd
      enddo

C ----------------------------------------------------------------------
C Save electron density difference
C ----------------------------------------------------------------------
      if (fildrh .ne. ' ') then
        do ispin = 1,nspin
          call reord( DRho(1,ispin), DRho(1,ispin), nml, nsm, +1)
        enddo
        call iorho( 'write', fildrh, cell, ntm, nsm, ntpl, nspin,
     .              DRho, found )
        do ispin = 1,nspin
          call reord( DRho(1,ispin), DRho(1,ispin), nml, nsm, -1)
        enddo
      endif

C Transform spin density into sum and difference
      if (nsd .eq. 2) then
        do ip = 1,ntpl
          rhotot = DRho(ip,1) + DRho(ip,2)
          DRho(ip,2) = DRho(ip,2) - DRho(ip,1)
          DRho(ip,1) = rhotot
        enddo
      endif

C ----------------------------------------------------------------------
C Calculate the dipole moment
C ----------------------------------------------------------------------
      dipol(1:3) = 0.0d0
      if (shape .ne. 'bulk') then

C Find center of system
        x0(1:3) = 0.0d0
        do ia = 1,nua
          x0(1:3) = x0(1:3) + xa(1:3,ia) / nua
        enddo

C Find dipole
        call reord( DRho, DRho, nml, nsm, +1)
        call dipole( cell, ntm, ntml(1), ntml(2), ntml(3), nsm,
     .    DRho, x0, dipol )
        call reord( DRho, DRho, nml, nsm, -1)

C Orthogonalize dipole to bulk directions
        if (shape .eq. 'chain') then
          const = dot( dipol, bcell, 3 ) / dot( bcell, bcell, 3 )
          dipol(1:3) = dipol(1:3) - const * bcell(1:3,1)
        elseif (shape .eq. 'slab') then
          call cross( bcell(1,1), bcell(1,2), b1Xb2 )
          const = dot( dipol, b1Xb2, 3 ) / dot( b1Xb2, b1Xb2, 3 )
          dipol(1:3) = const * b1Xb2(1:3)
        endif
      endif


C ----------------------------------------------------------------------
C Find Hartree potential of DRho = rhoscf-rhoatm. Store it in Vscf
C ----------------------------------------------------------------------

C Initialise Vscf
      Vscf(1:ntpl,1) = 0.0d0

C Solve Poisson's equation
      call reord( DRho, DRho, nml, nsm, +1)
      call poison( cell, ntml(1), ntml(2), ntml(3), ntm, DRho,
     .             DUscf, Vscf, DStres, nsm )
      call reord( DRho, DRho, nml, nsm, -1)
      call reord( Vscf, Vscf, nml, nsm, -1)

C Add contribution to stress from electrostatic energy of rhoscf-rhoatm
      if (istr .eq. 1) then
        stressl(1:3,1:3) = stressl(1:3,1:3) + DStres(1:3,1:3)
      endif

C ----------------------------------------------------------------------
C Find electrostatic (Hartree) energy of full SCF electron density
C ----------------------------------------------------------------------
      Uatm = Uharrs
      Uscf = Uatm + DUscf
      do ip = 1,ntpl
        Uscf = Uscf + dvol * Vscf(ip,1) * rhoatm(ip)
      enddo

#ifdef MPI
C Global reduction of Uscf/DUscf/Uatm
#ifdef NODAT
      call MPI_AllReduce(Uscf,Eloc,1,MPI_double_precision,MPI_Sum,
     .  MPI_Comm_World,MPIerror)
#else
      call MPI_AllReduce(Uscf,Eloc,1,DAT_double,MPI_Sum,
     .  MPI_Comm_World,MPIerror)
#endif
      Uscf = Eloc
#ifdef NODAT
      call MPI_AllReduce(DUscf,Eloc,1,MPI_double_precision,MPI_Sum,
     .  MPI_Comm_World,MPIerror)
#else
      call MPI_AllReduce(DUscf,Eloc,1,DAT_double,MPI_Sum,
     .  MPI_Comm_World,MPIerror)
#endif
      DUscf = Eloc
#ifdef NODAT
      call MPI_AllReduce(Uatm,Eloc,1,MPI_double_precision,MPI_Sum,
     .  MPI_Comm_World,MPIerror)
#else
      call MPI_AllReduce(Uatm,Eloc,1,DAT_double,MPI_Sum,
     .  MPI_Comm_World,MPIerror)
#endif
      Uatm = Eloc
#endif

C ----------------------------------------------------------------------
C Add neutral-atom potential to Vscf
C ----------------------------------------------------------------------

      Enaatm = 0.0d0
      Enascf = 0.0d0
      do ip = 1,ntpl
        Enaatm = Enaatm + dvol * Vna(ip) * rhoatm(ip)
        Enascf = Enascf + dvol * Vna(ip) * ( rhoatm(ip) + DRho(ip,1) )
        Vscf(ip,1) = Vscf(ip,1) + Vna(ip)
      enddo

#ifdef MPI
C Global reduction of Enaatm and Enascf
      Eloc = Enaatm
#ifdef NODAT
      call MPI_AllReduce(Eloc,Enaatm,1,MPI_double_precision,MPI_Sum,
     .  MPI_Comm_World,MPIerror)
#else
      call MPI_AllReduce(Eloc,Enaatm,1,DAT_double,MPI_Sum,
     .  MPI_Comm_World,MPIerror)
#endif
      Eloc = Enascf
#ifdef NODAT
      call MPI_AllReduce(Eloc,Enascf,1,MPI_double_precision,MPI_Sum,
     .  MPI_Comm_World,MPIerror)
#else
      call MPI_AllReduce(Eloc,Enascf,1,DAT_double,MPI_Sum,
     .  MPI_Comm_World,MPIerror)
#endif
#endif

C ----------------------------------------------------------------------
C Add potential from external electric field (if present)
C ----------------------------------------------------------------------

      if (frstme) then
        call efield( cell, nua, isa, xa, nzero, nsm, Vscf, field )
        IsEfld = .false.
        if (sqrt(dot(field,field,3)) .gt. tiny) IsEfld = .true.
      endif
      if (IsEfld) then
        call reord( Vscf, Vscf, nml, nsm, +1)
        call efield( cell, nua, isa, xa, ntm, nsm, Vscf, field )
        call reord( Vscf, Vscf, nml, nsm, -1)
        DUext = - dot( field, dipol, 3 )
      endif

C ----------------------------------------------------------------------
C Save electrostatic potential
C ----------------------------------------------------------------------

      if (filevh .ne. ' ') then
        call reord( Vscf, Vscf, nml, nsm, +1)
        call iorho( 'write', filevh, cell, ntm, nsm, ntpl, 1, Vscf, 
     .    found )
        call reord( Vscf, Vscf, nml, nsm, -1)
      endif

C Add contribution to stress from the derivative of the Jacobian of ---
C r->r' (strained r) in the integral of Vna*(rhoscf-rhoatm)
        if (istr .eq. 1) then
          do i = 1,3
            stress(i,i) = stress(i,i) + ( Enascf - Enaatm ) / volume
          enddo
        endif

C Get back spin density from sum and difference
      if (nsd .eq. 2) then
        do ip = 1,ntpl
          rhotot = DRho(ip,1)
          DRho(ip,1) = 0.5d0 * ( rhotot - DRho(ip,2) )
          DRho(ip,2) = 0.5d0 * ( rhotot + DRho(ip,2) )
        enddo
      endif

C ----------------------------------------------------------------------
C Exchange-correlation energy
C ----------------------------------------------------------------------

      Vaux(1:ntpl) = Vscf(1:ntpl,1)
      do ispin = 1,nsd
        DRho(1:ntpl,ispin) = DRho(1:ntpl,ispin) + rhoatm(1:ntpl) / nsd
        if (npcc .eq. 1) 
     .    DRho(1:ntpl,ispin) = DRho(1:ntpl,ispin) + rhopcc(1:ntpl)/nsd
      enddo

      do ispin = 1,nspin
        call reord(DRho(1,ispin),DRho(1,ispin),nml,nsm,+1)
      enddo

      call cellxc( XCfunc, XCauth, 0, 0, cell, ntml, ntml, ntpl, 0, 
     .             aux3, nspin, DRho, Ex, Ec, DEx, DEc, Vscf, 
     .             DVxcdn, stressl, 1, aux1 )

      do ispin = 1,nspin
        call reord(DRho(1,ispin),DRho(1,ispin),nml,nsm,-1)
        call reord(Vscf(1,ispin),Vscf(1,ispin),nml,nsm,-1)
      enddo

      Exc = Ex + Ec
      Dxc = DEx + DEc
      do ispin = 1,nsd
        do ip = 1,ntpl
          DRho(ip,ispin) = DRho(ip,ispin) - rhoatm(ip) / nsd
          if (npcc .eq. 1) 
     .      DRho(ip,ispin) = DRho(ip,ispin) - rhopcc(ip) / nsd
          Vscf(ip,ispin) = Vscf(ip,ispin) + Vaux(ip)
        enddo
      enddo

C ----------------------------------------------------------------------
C Save total potential
C ----------------------------------------------------------------------

      if (filevt .ne. ' ') then
        do ispin = 1,nspin
          call reord( Vscf(1,ispin), Vscf(1,ispin), nml, nsm, +1)
        enddo
        call iorho( 'write', filevt, cell, ntm, nsm, ntpl, nspin,
     .              Vscf, found )
        do ispin = 1,nspin
          call reord( Vscf(1,ispin), Vscf(1,ispin), nml, nsm, -1)
        enddo
      endif

C ----------------------------------------------------------------------
C Find SCF contribution to hamiltonian matrix elements
C ----------------------------------------------------------------------

      if (iHmat .eq. 1) then
        if (spiral) then
          call vmatsp( norb, nmpl, dvol, nspin, Vscf, maxnd,
     .               numd, listdptr, listd, Hmat, nuo,
     .               nuotot, iaorb, iphorb, isa, qspiral )
        else
          call vmat( norb, nmpl, dvol, nspin, Vscf, maxnd,
     .               numd, listdptr, listd, Hmat, nuo,
     .               nuotot, iaorb, iphorb, isa )
        endif
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
     .      nsd, dvol, volume, Vscf, Vaux, fal, stressl, ifa.ne.0,
     .      istr.ne.0 )
        endif

        if ( ifor .eq. 1)
     .    call forhar( ntpl, nspin, nml, ntml, ntm, npcc,
     .                 XCfunc, XCauth, cell, rhoatm, rhopcc,
     .                 Vna, DRho, Vscf, Vaux )

C Transform spin density into sum and difference
        if (nsd .eq. 2) then
          do ip = 1,ntpl
            rhotot = DRho(ip,1) + DRho(ip,2)
            DRho(ip,2) = DRho(ip,2) - DRho(ip,1)
            DRho(ip,1) = rhotot
          enddo
        endif

C Find contribution of neutral-atom potential
        call NeutralAtomOnMesh( na, isa, ntpl, Vna, indxua, dvol, 
     .                          volume, DRho, fal, stressl, 
     .                          ifa.ne.0, istr.ne.0 )

C Vaux is (minus) the potential which multiplies rhoatm
*       if (nsd .eq. 2) then
*         Vaux(1:ntpl) = 0.5d0 * Vaux(1:ntpl)
*       endif

*       call dfscf( norb, nuo, nuotot, ifa, istr, nmpl, maxnd,
*    .              numd, listdptr, listd, Dscf, datm, dvol, volume, 
*    .              Vscf, Vaux, DFascf, DFaatm, iaorb, indxua, iphorb, 
*    .              isa, nspin, fal, stressl )

        call dfscf( ifa, istr, na, norb, nuo, nuotot, nmpl, nspin,
     .              indxua, isa, iaorb, iphorb,
     .              maxnd, numd, listdptr, listd, Dscf, datm,
     .              Vscf, Vaux, dvol, volume, fal, stressl )

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
      if (allocated(DRho)) then
        call memory('D','S',size(DRho),'dhscf')
        deallocate(DRho)
      endif
      if (allocated(Vscf)) then
        call memory('D','S',size(Vscf),'dhscf')
        deallocate(Vscf)
      endif
      if (allocated(Vaux)) then
        call memory('D','S',size(Vaux),'dhscf')
        deallocate(Vaux)
      endif
      if (allocated(DVxcdn)) then
        call memory('D','S',size(DVxcdn),'dhscf')
        deallocate(DVxcdn)
      endif

      frstme = .false.

100   format(/,'dhscf: ',71(1h*))
101   format('dhscf:                  INPUT ERROR')
102   format('dhscf: ',71(1h*))

      end
