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
!***********************************************************************
! Contained modules:
!   module mesh    : Stores variables and arrays associated to the mesh
!   module meshphi : Stores variables and arrays connected with the
!                    calculation of basis functions phi on the mesh
! ----------------------------------------------------------------------
! Contained external subroutines:
! InitMesh          : initialises the mesh
! SameMeshAndAtoms  : checks whether the mesh or atoms have changed
! InitAtomMesh      : initialises quantities relating to the atoms
!                     and the mesh
! PartialCoreOnMesh : calculates the partial core density on mesh
! NeutralAtomOnMesh : calculates the neutral atom potential on mesh
! PhiOnMesh         : calculates the orbital values on the mesh
!***********************************************************************
!
!      module mesh
! ----------------------------------------------------------------------
! Stores the following variables and arrays connected with the mesh:
! ----------------------------------------------------------------------
! real*8  cmesh(3,3)    : Mesh-cell vectors
! real*8  dxa(3,na)     : Atom position within mesh-cell
! integer idop(mop)     : Extended-mesh-index displacement of points
!                       : within a sphere of radius rmax
! integer indexp(nep)   : Translation from extended to normal mesh index
! integer ipa(na)       : Mesh cell in which atom is
! integer mop           : Maximum number of non-zero orbital points
! integer nmeshg(3)     : Total number of mesh points in each direction
! integer ne(3)         : Number of mesh-Extension intervals 
!                       : in each direction
! integer nem(3)        : Extended-mesh divisions in each direction
! integer nmsc(3)       : Mesh divisions of each supercell vector
! integer nmuc(3)       : Mesh divisions of each unit cell vector
! integer nsm           : Number of mesh sub-divisions in each direction
! integer nsp           : Number of sub-points of each mesh point
! integer nusc          : Number of unit cells in each supercell dir.
! real*8  rcmesh(3,3)   : Reciprocal mesh-cell vectors
!                       : (WITHOUT 2*pi factor)
! real*8  xdop(3,mop)   : Vector to mesh points within rmax
! real*8  xdsp(3,nsp)   : Vector to mesh sub-points
! real*8  myBox(2,3)    : My processor's box of mesh points:
!                           myBox(1,ivec) : first point in box along ivec
!                           myBox(2,ivec) : last point in box alon ivec
! real*8  myExtBox(2,3) : My processor's extended box of mesh points
!                         In this module, myBox and myExtBox refer to
!                         super points, each containing nsm**3 subpoints.
!                         This is in contrast with their use in cellXC.
!
!***********************************************************************
!
! module meshphi
! ----------------------------------------------------------------------
! Stores variables and arrays related to basis orbitals on the mesh
! ----------------------------------------------------------------------
! logical DirectPhi     : If true the phi is calculated on the fly
! integer endpht(0:nmpl): Last position occupied by points in lstpht
! integer lstpht(ntopl) : List of non-zero orbitals at point
! integer listp2(ntopl) : Maps orbital-mesh point to iop
! integer nphi          : Length of phi array second dimension
! real    phi(nsp,ntopl): Basis orbitals at mesh points (sparse)
!
!***********************************************************************
!
! subroutine InitMesh
! ----------------------------------------------------------------------
! Initialises the mesh
! ----------------------------------------------------------------------
! Input :
! ----------------------------------------------------------------------
! integer na            : Number of atoms in supercell
! real*8  cell(3,3)     : Auxillary lattice vectors
! integer norb          : Total number of basis orbitals in supercell
! integer iaorb(norb)   : Atom to which each orbital belongs
! integer iphorb(norb)  : Orbital index (within atom) of each orbital
! integer isa(na)       : Species index of all atoms in supercell
! real*8  rmax          : Maximum orbital radius
! integer nsc(3)        : Number of unit-cells in each supercell direct.
! ----------------------------------------------------------------------
! Input and output :
! ----------------------------------------------------------------------
! real*8  G2max         : Effective planewave cutoff (Ry) determines
!                       : mesh density and the precision of integrals
!                       : On input : Value required
!                       : On output: Value used, which may be larger
! ----------------------------------------------------------------------
! Output :
! ----------------------------------------------------------------------
! real*8  G2mesh        : Effective planewave cutoff of mesh used
! integer nmpl          : Number of mesh points in unit cell locally
! integer nm(3)         : Number of Mesh divisions of each cell vector
! integer nml(3)        : Local form of nm
! integer ntm(3)        : Total number of mesh points
! integer ntml(3)       : Total number of mesh points stored locally
! integer ntpl          : Number of mesh Total Points in unit cell
!                       : (including subpoints) locally
! integer ntopl         : Total number of nonzero orbital points locally
! real*8  dvol          : Mesh-cell volume
! ----------------------------------------------------------------------
! Internal variables and arrays:
! ----------------------------------------------------------------------
! real*8  dx(3)         : Vector from atom to mesh sub-point
! real*8  dxp(3)        : Vector from atom to mesh point
! integer i             : General-purpose index
! integer ia            : Looping variable for number of atoms
! integer i1,i2,i3      : Mesh indexes in each mesh direction
! integer is            : Species index
! integer isp           : Sub-Point index
! integer ity           : Orbital-type index
! integer j             : General-porpose index
! integer j1,j2,j3      : Mesh indexes in each mesh direction
! real*8  k0(3)         : Zero-vector argument for routine chkgmx
! integer nep           : Number of extended-mesh points
! integer nmp           : Number of mesh points in unit cell
! integer nty           : Number of orbital types
! integer noty(nty)     : Number of orbitals of each type
! real*8  pldist        : Distance between mesh planes
! real*8  r             : Distance between atom and mesh point
! real*8  rcty(nty)     : Radius of each orbital type
! logical same          : Indicates whether orbital types are the same
! real*8  vecmod        : Vector modulus
! real*8  volume        : Unit cell volume
! logical within        : Is a mesh point within orbital range?
! ----------------------------------------------------------------------
! Units :
! ----------------------------------------------------------------------
! Energies in Rydbergs
! Distances in Bohr
!
!***********************************************************************
!
! subroutine SameMeshAndAtoms
! ----------------------------------------------------------------------
! Checks whether anything has changed that requires the mesh to be
! reinitialised or quantities relating to the atoms relative to
! the mesh.
! ----------------------------------------------------------------------
! Input :
! ----------------------------------------------------------------------
! integer na           : Number of atoms in the supercell
! real*8  xa(3,na)     : Coordinates of the atoms in the supercell
! real*8  ucell(3,3)   : Current unit cell vectors
! real*8  rmax         : Maximum orbital radius
! real*8  G2max        : Requested mesh cut-off
! real*8  G2mesh       : Mesh cut-off from last call
! ----------------------------------------------------------------------
! Output :
! ----------------------------------------------------------------------
! logical samesh       : If .true. then the mesh must be reinitialised
! logical samexa       : If .true. then atom related quantities need
!                      : to be recalculated
!
!***********************************************************************
!
! subroutine InitAtomMesh
! ----------------------------------------------------------------------
! Initialises the atom information relating to the mesh
! ----------------------------------------------------------------------
! Input :
! ----------------------------------------------------------------------
! integer na            : Number of atoms in supercell
! real*8  xa(3,na)      : Atomic positions of all atoms in supercell
! ----------------------------------------------------------------------
! Output :
! ----------------------------------------------------------------------
! all output quantities are in the mesh module
!
!***********************************************************************
!
! subroutine PartialCoreOnMesh
! ----------------------------------------------------------------------
! Finds the partial-core-correction energy density on the mesh
! ----------------------------------------------------------------------
! Input :
! ----------------------------------------------------------------------
! integer na            : Number of atoms in supercell
! integer isa(na)       : Species index of all atoms in supercell
! integer ntpl          : Number of mesh Total Points in unit cell
!                       : (including subpoints) locally
! integer indxua        : Index of equivalent atom in unit cell
! integer nsd           : Number of diagonal spin values (1 or 2)
! real*8  dvol          : Mesh-cell volume
! real*8  volume        : Unit cell volume
! real*4  Vscf(ntpl)    : Hartree potential of SCF density
! real*4  Vaux(ntpl)    : Auxiliary potential array
! logical Forces        : Are the forces to be calculated?
! logical Stress        : Are the stresses to be calculated?
! ----------------------------------------------------------------------
! Output : (non-gradient case)
! ----------------------------------------------------------------------
! real*4  rhopcc(ntpl)  : Partial-core-correction density for xc
! ----------------------------------------------------------------------
! Output : (gradient case)
! ----------------------------------------------------------------------
! real*8  fal(3,:)     : Local copy of atomic forces
! real*8  stressl(3,3)  : Local copy of stress tensor
!
!***********************************************************************
!
! subroutine NeutralAtomOnMesh
! ----------------------------------------------------------------------
! Finds the potential due to the neutral atoms on the mesh
! ----------------------------------------------------------------------
! Input :
! ----------------------------------------------------------------------
! integer na            : Number of atoms in supercell
! integer isa(na)       : Species index of all atoms in supercell
! integer ntpl          : Number of mesh Total Points in unit cell
!                       : (including subpoints) locally
! integer indxua        : Index of equivalent atom in unit cell
! real*8  dvol          : Mesh-cell volume
! real*8  volume        : Unit cell volume
! real    drho(ntpl)    : SCF density difference
! logical Forces        : Should the forces be calculated?
! logical Stress        : Should the stress be calculated?
! ----------------------------------------------------------------------
! Input/Output : (Output for non-gradient call, must be kept otherwise)
! ----------------------------------------------------------------------
! real    vna(ntpl)     : Sum of neutral-atom potentials
! ----------------------------------------------------------------------
! Output : (gradient call)
! ----------------------------------------------------------------------
! real*8  fal(3,:)     : Local copy of atomic forces
! real*8  stressl(3,3)  : Local copy of stress tensor
!
!***********************************************************************
!
! subroutine LocalChargeOnMesh
! ----------------------------------------------------------------------
! Finds the diffuse ionic charge, whose electrostatic potential is equal
! to Vlocal on the mesh
! ----------------------------------------------------------------------
! Input :
! ----------------------------------------------------------------------
! integer na            : Number of atoms in supercell
! integer isa(na)       : Species index of all atoms in supercell
! integer ntpl          : Number of mesh Total Points in unit cell
!                       : (including subpoints) locally
! integer indxua        : Index of equivalent atom in unit cell
! ----------------------------------------------------------------------
! Output : 
! ----------------------------------------------------------------------
! real    Chlocal(ntpl)     : Sum of diffuse ionic charge densities
! 
!***********************************************************************
!
! subroutine PhiOnMesh
! ----------------------------------------------------------------------
! Calculates the values of the orbitals at the mesh points
! ----------------------------------------------------------------------
! Input :
! ----------------------------------------------------------------------
! integer iaorb(norb)   : Atom to which each orbital belongs
! integer iphorb(norb)  : Orbital index (within atom)
! integer isa(na)       : Species index of all atoms in supercell
! integer nmpl          : Number of mesh points in unit cell locally
! integer norb          : Total number of basis orbitals in supercell
! ----------------------------------------------------------------------
! Output :
! ----------------------------------------------------------------------
! All output quantities are in the module meshphi
! 
!***********************************************************************

      module mesh

      use precision, only : dp

      implicit none

      integer, pointer, save ::
     .  idop(:)   =>null(), ! Mesh-index span of points 
                            ! within an atomic sphere
     .  indexp(:) =>null(), ! Translation between mesh indexes in 
                            ! myExtBox and myBox
     .  ipa(:)    =>null()  ! Index of atoms' mesh cells within myExtBox

      real(dp), pointer, save :: 
     .  dxa(:,:)  =>null(), ! Atomic coordinates within their mesh cells
     .  xdop(:,:) =>null(), ! Coordinates of (super)points 
                            ! within an atomic sphere
     .  xdsp(:,:) =>null()  ! Coordinates of subpoints
                            ! relative to their superpoint

      integer,save :: myBox(2,3)    ! My processor's box of mesh points:
                                    !   myBox(1,:) : lower bounds
                                    !   myBox(2,:) : upper bounds
      integer,save :: myExtBox(2,3) ! My extended box of mesh points
      integer,save :: mop           ! Accumulated num. of orbital points
      integer,save :: ne(3)         ! (Half) width of mesh extension 
                                    ! 'wings' (not including subpoints)
      integer,save :: nem(3)        ! Widths of my processor's extended
                                    ! box of mesh (super)points
      integer,save :: nmsc(3)       ! Number of (super)points along each
                                    ! supercell vector
      integer,save :: nmeshg(3)     ! Number of mesh points in each unit
                                    ! cell vector, including subpoints
      integer,save :: nsm           ! Number of subpoints per superpoint
                                    ! along each mesh direction
      integer,save :: nsp           ! Total number of subpoints per
                                    ! superpoint: nsp=nsm**3
      integer,save :: nmuc(3)       ! Mesh points in each unit cell vector
      integer,save :: nusc(3)       ! Unit cells in each supercell dir.

      real(dp),save:: cmesh(3,3)    ! Mesh cell vectors
      real(dp),save:: rcmesh(3,3)   ! Reciprocal mesh cell vectors

      end module mesh

!***********************************************************************

      module meshphi

      use precision, only : grid_p

      implicit none

      logical, save :: DirectPhi  ! (Re)calculate on the fly basis
                                  ! orbitals at mesh points?
      integer, save :: nphi

      integer,pointer,save:: endpht(:)=>null()
      integer,pointer,save:: lstpht(:)=>null()
      integer,pointer,save:: listp2(:)=>null()

      real(grid_p),pointer,save:: 
     .  phi(:,:)=>null() ! Basis orbitals at mesh points (sparse format)

      end module meshphi

!***********************************************************************

      subroutine InitMesh( na, cell, norb, iaorb, iphorb, isa,
     .                     rmax, G2max, G2mesh, nsc, nmpl, 
     .                     nm, nml, ntm, ntml, ntpl, ntopl, dvol)

! Modules
      use precision,     only : dp
      use parallel,      only : Node, Nodes
      use parallel,      only : ProcessorY
      use parallelsubs,  only : HowManyMeshPerNode, GlobalToLocalMesh
      use atmfuncs,      only : rcut
      use mesh
      use meshphi,       only : endpht
      use siesta_cml
      use cellsubs,      only : reclat  ! Finds reciprocal unit cell
      use cellsubs,      only : volcel  ! Finds unit cell volume
      use m_chkgmx,      only : chkgmx  ! Checks planewave cutoff of a mesh
      use fft1d,         only : nfft    ! Finds allowed value for 1-D FFT
      use alloc,    only : de_alloc     ! Deallocation routine
      use alloc,    only : re_alloc     ! (Re)allocation routine
      use siestaXC, only : myMeshBox    ! Returns my processor mesh box
      use siestaXC, only : setMeshDistr ! Sets a distribution of mesh
                                        ! points over parallel processors
      implicit none

! Passed arguments
      integer, intent(in):: na           ! Number of atoms in supercell
      real(dp),intent(in):: cell(3,3)    ! Unit cell vectors
      integer, intent(in):: norb         ! Total number of basis 
                                         ! orbitals in supercell
      integer, intent(in):: iaorb(norb)  ! Atom to which orbitals belong
      integer, intent(in):: iphorb(norb) ! Orbital index (within atom) 
                                         ! of each orbital
      integer, intent(in):: isa(na)      ! Species index of each atom
      real(dp),intent(in):: rmax         ! Maximum orbital radius
      integer, intent(in):: nsc(3)       ! Number of unit-cells in each
                                         ! supercell direction
      real(dp),intent(inout):: G2max     ! Effective planewave cutoff
                                         ! On input : Value required
                                         ! On output: Value used, which
                                         ! may be larger
      real(dp),intent(out):: G2mesh      ! Planewave cutoff of mesh used
                                         ! (same as G2max on output)
      integer, intent(out):: nmpl        ! Number of mesh points stored
                                         ! by my processor
      integer, intent(out):: nm(3)       ! Number of mesh divisions of
                                         ! each unit cell vector
      integer, intent(out):: nml(3)      ! Sizes of my processor's box
                                         ! of mesh points
      integer, intent(out):: ntm(3)      ! Mesh divisions of each unit
                                         ! cell vector, incl. subpoints
      integer, intent(out):: ntml(3)     ! Sizes of my processor's mesh
                                         ! box, including subpoints
      integer, intent(out):: ntpl        ! Mesh points stored by my 
                                         ! processor, incl. subpoints
      integer, intent(out):: ntopl       ! Total number of nonzero 
                                         ! orbital points stored by my
                                         ! processor
      real(dp),intent(out):: dvol        ! Volume per mesh (super)point

! External functions used
      real(dp),external:: dismin ! Distance of a point to a unit cell

! Local parameters, variables, and arrays
      character(len=*),parameter:: myName = 'InitMesh '
      logical :: same, within
      integer,save:: JDGdistr=-1
      integer :: i, ia, i1, i2, i3, io, iphi, is, isp, ity, 
     .           j, j1, j2, j3, k, k1, k2, k3, 
     .           ncells, nep, ntop, ntp, nty
      real(dp):: boxWidth(3), dx(3), dxp(3), extWidth(3),
     .           pldist, r, rcell(3,3), vecmod, volume
      real(dp),         save :: k0(3)=0
      integer, pointer, save :: iphty(:)=>null(), isty(:)=>null(),
     .                          noty(:)=>null()
      real(dp),pointer, save :: rcty(:)=>null()

! ----------------------------------------------------------------------
! Orbital type initialization 
! ----------------------------------------------------------------------

! Find number of orbital types
      call re_alloc( isty,  1,norb, myName//'isty'  )
      call re_alloc( iphty, 1,norb, myName//'iphty' )
      nty = 0
      do io = 1,norb
        ia = iaorb(io)
        is = isa(ia)
        iphi = iphorb(io)
        ity = 0
        same = .false.
        do while (ity .lt. nty .and. .not.same)
          ity = ity + 1
          same = (is.eq.isty(ity) .and. iphty(ity).eq.iphi)
        enddo
        if (.not.same) then
          nty = nty + 1
          isty(nty) = is
          iphty(nty) = iphi
        endif
      enddo

! Allocate local arrays dependent on nty
      call re_alloc( iphty, 1,nty, myName//'iphty', shrink=.true. )
      call re_alloc( isty,  1,nty, myName//'isty',  shrink=.true. )
      call re_alloc( noty,  1,nty, myName//'noty' )
      call re_alloc( rcty,  1,nty, myName//'rcty' )

! Find number of orbital types, number of orbitals of each type,
! and cutoff radius of each type.
      noty(1:nty) = 0
      nty = 0
      io_loop: do io = 1,norb
        ia = iaorb(io)
        is = isa(ia)
        iphi = iphorb(io)
        do ity = 1,nty
          if (is.eq.isty(ity) .and. iphty(ity).eq.iphi) then
            noty(ity) = noty(ity) + 1
            cycle io_loop
          endif
        enddo
        nty = nty + 1
        isty(nty) = is
        iphty(nty) = iphi
        noty(nty) = 1
        rcty(nty) = rcut( is, iphi )
      enddo io_loop

      call de_alloc( iphty, myName//'iphty' )
      call de_alloc( isty,  myName//'isty'  )

! ----------------------------------------------------------------------
! Mesh initialization 
! ----------------------------------------------------------------------

! Find reciprocal cell vectors (multiplied by 2*pi)
      call reclat( cell, rcell, 1 )

! Find number of mesh intervals for each cell vector.
! The reciprocal vectors of the mesh unit cell (cell/ntm)
! are rcell*ntm, and must be larger than 2*G2max
      do i = 1,3
        vecmod = sqrt(dot_product(rcell(:,i),rcell(:,i)))
        ntm(i) = 2 * sqrt(G2max) / vecmod + 1
      enddo

! Impose that mesh cut-off is large enough
      impose_cutoff: do

        ! Require that ntm is suitable for FFTs and a multiple of nsm
        do i = 1,3
          impose_fft: do 
            ! nfft increases ntm to the next integer suitable for FFTs
            call nfft( ntm(i) )
            if ( mod( ntm(i), nsm )==0 ) exit impose_fft
            ntm(i) = ntm(i) + 1
          end do impose_fft
          nm(i) = ntm(i) / nsm
          nmsc(i) = nm(i) * nsc(i)
        enddo ! i

        ! Check that effective cut-off is large enough as for non-right
        ! angled unit cells this is not guaranteed to be the case. 
        ! If cut-off needs to be larger, increase ntm and try again.
        ! chkgmx decreases G2mesh to the right cutoff value of the mesh
        G2mesh = 1.0d6
        call chkgmx( k0, rcell, ntm, G2mesh )
        if (G2mesh >= G2max) then
          exit impose_cutoff
        else
          ntm(1:3) = ntm(1:3) + 1
        end if

      end do impose_cutoff

! Save some numbers
      nmeshg(:) = ntm(:) ! Total number of mesh points, incl. subpoints
      nmuc(:) = nm(:)    ! Mesh points in each unit cell vector
      nusc(:) = nsc(:)   ! Unit cells in each supercell direction

! Find number of mesh points in unit cell.
      ntp = ntm(1) * ntm(2) * ntm(3)

! Find local number of Mesh points of each kind
      call HowManyMeshPerNode(nm,Node,Nodes,nmpl,nml)
      ntpl = nmpl*nsp
      ntml(1:3) = nml(1:3)*nsm

! Find volume of unit cell and of mesh cell
      volume = volcel( cell )
      dvol = volume / ntp

! Output current mesh dimensions and cut-off
      if (Node.eq.0) then
        write(6,'(/,a,3(i6,a),i12)') 'InitMesh: MESH =',
     .    ntm(1),' x',ntm(2),' x',ntm(3),' =', ntp
        write(6,'(a,2f10.3,a)')
     .   'InitMesh: Mesh cutoff (required, used) =',G2max,G2mesh,' Ry'
      endif
      if (cml_p) then
        call cmlStartPropertyList(mainXML)
        call cmlAddProperty(xf=mainXML, value=ntm,
     .       dictref='siesta:ntm', title='Mesh',
     .       units='cmlUnits:countable')
        call cmlAddProperty(xf=mainXML, value=G2max,
     .       units='siestaUnits:Ry',
     .       dictref='siesta:g2max', title='Requested Cut-Off')     
        call cmlAddProperty(xf=mainXML, value=G2mesh,
     .       units='siestaUnits:Ry',
     .       dictref='siesta:g2mesh', title='Actual Cut-Off')
        call cmlEndPropertyList(mainXML)
      endif
      G2max = G2mesh

! Find mesh-cell vectors
      do i = 1,3
        do j = 1,3
          cmesh(j,i) = cell(j,i) / nm(i)
        enddo
      enddo

! Find reciprocal mesh-cell vectors (not multiplied by 2*pi)
      call reclat( cmesh, rcmesh, 0 )

! Find number ne(:) of extended-mesh intervals for each cell vector.
      do i = 1,3
        ! pldist is the distance between mesh planes
        pldist = 1.0_dp / sqrt(dot_product(rcmesh(:,i),rcmesh(:,i)))
        ! Find number of planes spanned by rmax
        ne(i) = rmax / pldist
      enddo ! i

! For an atom at x=0, ne=rmax/pldist is the last mesh point within rmax.
! But an atom at mesh cell ix=0 can be almost at x=dx. Therefore we have
! to go up to ne+1. And subpoints in mesh cell ix=-ne-1 can be almost at
! x=-ne*dx. Therefore we have to go up to -ne-1 to the left. Thus, we 
! just increase ne and forget about these two effects from now on.
      ne(:) = ne(:) + 1

! Set uniform distribution of mesh points and find my processor mesh box.
! This is the distribution presently used throughout this JDG subroutine
! but it could be changed in the future.
      call setMeshDistr( distrID=JDGdistr, nMesh=ntm, 
     .                   nNodesX=1, nNodesY=ProcessorY, nBlock=nsm )
      call myMeshBox( ntm, JDGdistr, myBox )
      myBox = myBox / nsm   ! Here myBox will refer to 'super' points

! Add 'wings extensions' to myBox of mesh points, containing all points 
! within the orbital spheres that may intersect myBox. Thus, wings must
! be one diameter wide.
      myExtBox(1,:) = myBox(1,:) - 2*ne(:)
      myExtBox(2,:) = myBox(2,:) + 2*ne(:)

! Find number of extended-box points.
      nem(:) = myExtBox(2,:) - myExtBox(1,:) + 1
      nep = nem(1) * nem(2) * nem(3)

! Allocate local memory relating to nep
      call re_alloc( indexp, 1,nep, myName//'indexp' )

! Find relationship between extended and normal box points
      boxWidth(:) =    myBox(2,:)    - myBox(1,:) + 1
      extWidth(:) = myExtBox(2,:) - myExtBox(1,:) + 1
      do i3 = myExtBox(1,3),myExtBox(2,3)
        do i2 = myExtBox(1,2),myExtBox(2,2)
          do i1 = myExtBox(1,1),myExtBox(2,1)

            ! Find periodic indexes within first supercell
            j1 = modulo( i1, nmsc(1) )
            j2 = modulo( i2, nmsc(2) )
            j3 = modulo( i3, nmsc(3) )

            ! Find indexes relative to box origins
            j1 = j1 -    myBox(1,1)
            j2 = j2 -    myBox(1,2)
            j3 = j3 -    myBox(1,3)
            k1 = i1 - myExtBox(1,1)
            k2 = i2 - myExtBox(1,2)
            k3 = i3 - myExtBox(1,3)

            ! Find combined mesh indexes.
            j = 1 + j1 + boxWidth(1)*j2 + boxWidth(1)*boxWidth(2)*j3
            k = 1 + k1 + extWidth(1)*k2 + extWidth(1)*extWidth(2)*k3

            ! Find myExtBox -> myBox index translation
            if (j1>=0 .and. j1<boxWidth(1) .and.
     .          j2>=0 .and. j2<boxWidth(2) .and.
     .          j3>=0 .and. j3<boxWidth(3)) then
              indexp(k) = j   ! Point is within myBox
            else
              indexp(k) = 0   ! Point is outside myBox
            endif ! (j>=myBox(1) and j<=myBox(2))

          enddo ! i1
        enddo ! i2
      enddo ! i3

! Find sub-points
      isp = 0
      do i3 = 0, nsm-1
        do i2 = 0, nsm-1
          do i1 = 0, nsm-1
            isp = isp + 1
            do i = 1,3
              xdsp(i,isp) = ( cmesh(i,1) * i1 +
     .                        cmesh(i,2) * i2 +
     .                        cmesh(i,3) * i3 ) / nsm
            enddo
          enddo
        enddo
      enddo

! Find number of points within rmax (orbital points)
      mop = 0
      ntop = 0
      do i3 = -ne(3), ne(3)    ! Loop over neighbor (super) points
        do i2 = -ne(2), ne(2)
          do i1 = -ne(1), ne(1)
            dxp(:) = cmesh(:,1) * i1 +  ! (Super) point coordinates
     .               cmesh(:,2) * i2 +
     .               cmesh(:,3) * i3 
            ! Find if any subpoint can be within rmax of an atom that
            ! is in the origin's mesh cell
            within = .false.
            do isp = 1,nsp  ! Loop over sub-points
              dx(1:3) = dxp(1:3) + xdsp(1:3,isp)  ! Subpoint coords.
              ! Find distance from point to the mesh cell at the origin
              ! (since the atom might be at any point of the mesh cell)
              r = dismin( cmesh, dx )
              if ( r .lt. rmax ) within = .true.
            enddo ! isp
            if ( within ) then
              mop = mop + 1      ! Total number of points within rmax
              ! Find added num. of points within rmax of all atom types
              do ity = 1,nty
                if ( r .lt. rcty(ity) ) ntop = ntop + noty(ity)
              enddo ! ity
            endif ! (within)
          enddo ! i1
        enddo ! i2
      enddo ! i3

! Calculate ntopl
      ncells = nsc(1) * nsc(2) * nsc(3)
      ntop = ntop / ncells
      ntopl = ((ntop/Nodes) + 1)

! Allocate local memory that depends on mop
      call re_alloc( idop,      1,mop, myName//'idop' )
      call re_alloc( xdop, 1,3, 1,mop, myName//'xdop' )

! Find points within rmax (orbital points)
      mop = 0
      do i3 = -ne(3), ne(3)
        do i2 = -ne(2), ne(2)
          do i1 = -ne(1), ne(1)
            dxp(:) = cmesh(:,1) * i1 +
     .               cmesh(:,2) * i2 +
     .               cmesh(:,3) * i3 
            within = .false.
            do isp = 1,nsp
              dx(1:3) = dxp(1:3) + xdsp(1:3,isp) 
              r = dismin( cmesh, dx )
              if ( r .lt. rmax ) within = .true.
            enddo
            if ( within ) then
              mop = mop + 1
              ! Store index-distance and vector-distance to point.
              idop(mop) = i1 + nem(1) * i2 + nem(1) * nem(2) * i3
              xdop(1:3,mop) = dxp(1:3)
            endif ! (within)
          enddo ! i1
        enddo ! i2
      enddo ! i3

! Deallocate orbital related memory
      call de_alloc( noty, myName//'noty' )
      call de_alloc( rcty, myName//'rcty' )

! Allocate memory related to nmpl
      call re_alloc( endpht, 0,nmpl, myName//'endpht' )

      end subroutine InitMesh

!***********************************************************************

      subroutine SameMeshAndAtoms( na, xa, ucell, rmax, G2max, G2mesh,
     .                             samesh, samexa)

! Used modules
      use precision, only : dp
      use alloc,     only : de_alloc     ! Deallocation routine
      use alloc,     only : re_alloc     ! (Re)allocation routine

      implicit none

! Passed arguments
      integer,   intent(in)  :: na
      real(dp),  intent(in)  :: G2max, G2mesh
      real(dp),  intent(in)  :: rmax
      real(dp),  intent(in)  :: ucell(3,3)
      real(dp),  intent(in)  :: xa(3,na)
      logical,   intent(out) :: samesh, samexa

! Local variables
! Saved
      character(len=*),parameter:: myName = 'SameMeshAndAtoms '
      integer, save :: lastna = 0  ! Number of atoms in last call
      real(dp), save :: tiny       = 1.0e-12_dp,
     .  lastc(3,3) = 0.743978657912656e50_dp, ! Unit cell in last call
     .  lastra     = 0.0_dp                   ! Rmax in last call

! non-saved
      integer:: i, ia, j
      real(dp),pointer,save::
     .  lastxa(:,:)=>null()  ! Position of atoms from last call

! If number of atoms has changed then reallocate lastxa and
! initialise it to an unlikely number
      if (na.ne.lastna) then
        call re_alloc( lastxa, 1,3, 1,na, myName//'lastxa',
     .                 shrink=.true. )
        lastxa(1:3,1:na) = 0.0_dp
        lastxa(1,1) = 0.987654321273567e50_dp
      endif

! Find if mesh has to be changed due to unit cell
      samesh = .true.
      do i = 1,3
        do j = 1,3
          if ( ucell(j,i) .ne. lastc(j,i) ) samesh = .false.
          lastc(j,i) = ucell(j,i)
        enddo
      enddo

! Find if mesh has to be changed due to unit cell
      if ( G2max .gt. G2mesh * (1.0_dp + tiny) ) samesh = .false.

! Find if mesh has to be changed due to rmax
      if (rmax .ne. lastra) samesh = .false.
      lastra = rmax

! Find if atoms have moved having checked the number of atoms first
      samexa = (na.eq.lastna)
      if (samexa) then
        do ia = 1,na
          do i = 1,3
            if ( xa(i,ia) .ne. lastxa(i,ia) ) samexa = .false.
          enddo
        enddo
      endif

! Copy the number of atoms and coordinates to save for next call
      lastna = na
      lastxa(1:3,1:na) = xa(1:3,1:na)

! If cell has changed then it is necessary to reinitialise coordinates
      if (.not.samesh) samexa = .false.

      end subroutine SameMeshAndAtoms

!***********************************************************************

      subroutine InitAtomMesh( na, xa)

      use precision, only: dp
      use mesh, only: ipa      ! atoms' mesh points within myExtBox
      use mesh, only: dxa      ! atom's position relative to mesh point
      use mesh, only: cmesh    ! Mesh cell vectors
      use mesh, only: myBox    ! My processor's box of mesh points
      use mesh, only: myExtBox ! My extended box of mesh points
      use mesh, only: ne       ! Points in rmax along each lat. vector
      use mesh, only: nem      ! Widths of myExtBox
      use mesh, only: nmuc     ! Mesh points in each unit cell vector
      use mesh, only: nmsc     ! Mesh points in each supercell vector
      use mesh, only: nusc     ! Unit cells in each supercell direction
      use mesh, only: rcmesh   ! Reciprocal mesh cell vectors
      use alloc,only: de_alloc ! Deallocation routine
      use alloc,only: re_alloc ! (Re)allocation routine
      use parallel, only: node, nodes

      implicit none

! Passed arguments
      integer,  intent(in) :: na        ! Number of atoms in supercell
      real(dp), intent(in) :: xa(3,na)  ! Atomic coordinates

! Internal variables
      character(len=*),parameter:: myName = 'InitAtomMesh '
      integer, save :: lastna = 0 ! Number of atoms on last call
      integer :: ia, ixa(3), j1, j2, j3, jsc(3), jxa(3)
      real(dp):: cxa(3)

! (Re)allocate atom-mesh arrays
      if (na.ne.lastna) then
        call re_alloc( ipa,      1,na, myName//'ipa' )
        call re_alloc( dxa, 1,3, 1,na, myName//'dxa' )
      endif
      lastna = na

! Find atomic positions relative to mesh
      do ia = 1,na

        ! Find atomic coordinates in mesh basis vectors
        cxa(:) = matmul( xa(:,ia), rcmesh(:,:) )

        ! Find indexes of mesh cell in which atom is
        ixa(:) = floor(cxa(:))

        ! Find atom position within mesh cell
        cxa(:) = cxa(:) - ixa(:)                  ! Mesh coordinates
        dxa(:,ia) = matmul( cmesh(:,:), cxa(:) )  ! Cartesian coords

        ! Find atom's mesh index within myExtBox, but only if the
        ! (periodic) atomic sphere (of width ne) can intersect myBox.
        ! Otherwise, make ipa=0. nem is the width of myExtBox.
        ! nmsc is the width of the supercell.
        ! Notice that, even if several periodic images of the same 
        ! atom may fit within myExtBox, only one must be considered,
        ! since the points within its sphere, once folded into myBox,
        ! will be equivalent to those of any other periodic image.
        ipa(ia) = 0
        sc: do j3 = -1,1   ! Loop on neighbor supercells
            do j2 = -1,1
            do j1 = -1,1
              jsc(:) = (/j1,j2,j3/)
              jxa(:) = ixa(:) + jsc(:)*nmsc(:)
              if ( all(jxa(:)>=myBox(1,:)-ne(:)) .and. 
     .             all(jxa(:)<=myBox(2,:)+ne(:)) ) then
                jxa(:) = jxa(:) - myExtBox(1,:)
                ipa(ia) = 1 + jxa(1) + nem(1)*jxa(2)
     .                      + nem(1)*nem(2)*jxa(3)
                
                exit sc
              end if
            end do    ! j1
          end do      ! j2
        end do sc     ! j3
          
      enddo ! ia

      end subroutine InitAtomMesh

!***********************************************************************

      subroutine PartialCoreOnMesh( na, isa, ntpl, rhopcc, indxua, 
     .                              nsd, dvol, volume, Vscf, Vaux, 
     .                              fal, stressl, Forces, Stress )

! Used modules
      use precision,     only : dp, grid_p
      use atmfuncs, only: rcore, chcore_sub
      use mesh,     only: dxa, idop, indexp, ipa, mop, nsp, xdop, xdsp

      implicit none

! Passed arguments
      integer,  intent(in)    :: na, isa(na), ntpl
      integer,  intent(in)    :: indxua(*)
      integer,  intent(in)    :: nsd
      real(dp), intent(in)    :: dvol, volume
      real(grid_p),     intent(in)    :: Vscf(ntpl,*)
      real(grid_p),     intent(in)    :: Vaux(ntpl)
      logical,  intent(in)    :: Forces
      logical,  intent(in)    :: Stress
      real(grid_p),     intent(out)   :: rhopcc(ntpl)
      real(dp), intent(inout) :: fal(3,*)
      real(dp), intent(inout) :: stressl(3,3)

! Internal variables
      real(dp),parameter:: tiny = 1.0e-12_dp
      integer :: i, ia, iop, ip, ip0, is, isp, ispin, iua
      real(dp):: dfa(3), dx(3), grrho(3), r, ra, rhop, vxc
      logical :: Gradients

! Find out whether this is a gradient run based on the arguments passed
      Gradients = ( Forces .or. Stress )

! Initialise array to zero
      if (.not.Gradients) rhopcc(1:ntpl) = 0.0_dp

      do ia = 1,na
        if (ipa(ia)==0) cycle ! Atomic sphere does not intersec myBox
        iua = indxua(ia)
        is = isa(ia)
        ra = rcore( is )
        if (ra .gt. tiny) then

! Loop over mesh points inside rmax
          do iop = 1,mop
            ip0 = indexp( ipa(ia) + idop(iop) )
            if (ip0 .gt. 0) then

              ! Loop over sub-points
              do isp = 1,nsp
                dx(1:3) = xdop(1:3,iop) + xdsp(1:3,isp) - dxa(1:3,ia)
                r = sqrt(dot_product(dx,dx))
!                if (r .lt. ra .and. r .gt. tiny) then
                if (r .lt. ra) then
                  ip = isp + nsp * (ip0 - 1)
                  call chcore_sub( is, dx, rhop, grrho )
                  if (Gradients) then
                    ! Calculate gradients of PCC
                    do ispin = 1,nsd
                      vxc = Vscf(ip,ispin) - Vaux(ip)
                      do i = 1,3
                        dfa(i) = dvol * vxc * grrho(i) / nsd
                        if (Forces) fal(i,iua) = fal(i,iua) + dfa(i)
                        if (Stress) stressl(1:3,i) = stressl(1:3,i) + 
     .                    dx(1:3)*dfa(i)/volume
                      enddo
                    enddo
                  else
                    ! Calculate density due to PCC
                    rhopcc(ip) = rhopcc(ip) + rhop
                  endif
                endif
              enddo

            endif
          enddo

        endif
      enddo

      end subroutine PartialCoreOnMesh

!***********************************************************************

      subroutine NeutralAtomOnMesh( na, isa, ntpl, vna,
     .                              indxua, dvol, volume,
     .                              drho, fal, stressl,
     .                              Forces, Stress )

! Used modules
      use precision,     only : dp, grid_p
      use atmfuncs, only: rcut, phiatm
      use mesh,     only: dxa, idop, indexp, ipa, mop, nsp, xdop, xdsp

      implicit none

! Passed arguments
      integer,  intent(in)    :: na, isa(na), ntpl
      integer,  intent(in)    :: indxua(*)
      real(dp), intent(in)    :: dvol, volume
      real(grid_p),     intent(in)    :: drho(*)
      logical,  intent(in)    :: Forces, Stress
      real(grid_p),     intent(inout)   :: vna(ntpl)
      real(dp), intent(inout) :: fal(3,*)
      real(dp), intent(inout) :: stressl(3,3)

! Internal variables
      integer :: i, ia, iop, ip, ip0, is, isp, iua
      real(dp):: dx(3), grva(3), r, ra, va
      logical :: Gradients

! Check whether forces and/or stress has been requested based on arguments
      Gradients = (Forces.or.Stress)

! Initialise array to zero if not computing gradients
! If we are computing gradients, the current value of vna
! must be kept. This is particularly important for grid-cell sampling...

      if (.not.Gradients) vna(1:ntpl) = 0.0_dp

      do ia = 1,na
        if (ipa(ia)==0) cycle ! Atomic sphere does not intersec myBox
        iua = indxua(ia)
        is = isa(ia)
        ra = rcut( is, 0 )

        ! Loop over mesh points inside rmax
        do iop = 1,mop
          ip0 = indexp( ipa(ia) + idop(iop) )
          if (ip0 .gt. 0) then

            ! Loop over sub-points
            do isp = 1,nsp
              dx(1:3) = xdop(1:3,iop) + xdsp(1:3,isp) - dxa(1:3,ia)
              r = sqrt(dot_product(dx,dx))
              if (r .lt. ra) then
                ip = isp + nsp * (ip0 - 1)
                call phiatm( is, 0, dx, va, grva )
                if (Gradients) then
                  do i = 1,3
                    grva(i) = dvol * grva(i) * drho(ip)
                    if (Forces) fal(i,iua) = fal(i,iua) + grva(i)
                    if (Stress) stressl(1:3,i) = stressl(1:3,i) +
     .                dx(1:3) * grva(i) / volume
                  enddo
                else
                  vna(ip) = vna(ip) + va
                endif
              endif
            enddo

          endif
        enddo

      enddo

      end subroutine NeutralAtomOnMesh

!***********************************************************************

      subroutine LocalChargeOnMesh( na, isa, ntpl, Chlocal, indxua )

! Used modules
      use precision,     only : dp, grid_p
      use atmfuncs, only: rcut, psch
      use atm_types, only : species
      use radial, only : rad_func
      use mesh,     only: dxa, idop, indexp, ipa, mop, nsp, xdop, xdsp
 
      implicit none
 
! Passed arguments
      integer, intent(in)  :: na, isa(na), ntpl
      integer, intent(in)  :: indxua(*)
      real(grid_p),    intent(out) :: Chlocal(ntpl)
 
! Internal variables
      integer :: ia, iop, ip, ip0, is, isp, iua
      real(dp):: dx(3), grpsch(3), r, ra, pschloc
      type(rad_func), pointer:: func   
 
! Initialise array to zero 
      Chlocal(1:ntpl) = 0.0_grid_p

! Loop on atoms
      do ia = 1,na
        if (ipa(ia)==0) cycle ! Atomic sphere does not intersec myBox
        iua = indxua(ia)
        is = isa(ia)
        func => species(is)%chlocal
        ra = func%cutoff
 
        ! Loop over mesh points inside rmax
        do iop = 1,mop
          ip0 = indexp( ipa(ia) + idop(iop) )
          if (ip0 .gt. 0) then
 
            ! Loop over sub-points
            do isp = 1,nsp
              dx(1:3) = xdop(1:3,iop) + xdsp(1:3,isp) - dxa(1:3,ia)
              r = sqrt(dot_product(dx,dx))
              if (r .lt. ra) then
                ip = isp + nsp * (ip0 - 1)
                call psch( is, dx, pschloc, grpsch )
                Chlocal(ip) = Chlocal(ip) + pschloc
              endif
            enddo
 
          endif
        enddo
 
      enddo
 
      end subroutine LocalChargeOnMesh

!***********************************************************************

      subroutine PhiOnMesh( nmpl, norb, iaorb, iphorb, isa)

! Modules
      use precision,     only : dp, grid_p
      use atmfuncs, only : rcut, phiatm
      use mesh,     only : dxa, idop, indexp, ipa, mop, nsp, xdop, xdsp
      use meshphi,  only : directphi, nphi, phi, lstpht, listp2, endpht
      use parallel, only : Node
      use fdf,      only : fdf_boolean
      use alloc,    only : de_alloc ! Deallocation routine
      use alloc,    only : re_alloc ! (Re)allocation routine
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

! Passed arguments
      integer,intent(in):: nmpl         ! Number of mesh (super) points 
                                        ! stored by my processor
      integer,intent(in):: norb         ! Total number of basis orbitals
                                        ! in supercell
      integer,intent(in):: iaorb(norb)  ! Atom to which each orbital belongs
      integer,intent(in):: iphorb(norb) ! Orbital index (within atom)
      integer,intent(in):: isa(*)       ! Species index of all atoms in
                                        ! supercell

! Local variables
      character(len=*),parameter:: myName = 'PhiOnMesh '
      logical :: within
      integer :: ia, io, iop, ip, ip0, iphi, is, isp,
     .           n, nlist, nliste, nmax
      real(dp):: dxsp(3,nsp), grphi(3), phip, r2o, r2sp(nsp)
      logical, save :: firsttime = .true.
      integer, pointer, save :: numphi(:)=>null()
#ifdef MPI
      integer:: MPIerror
#endif

! On first call set the logical DirectPhi
      if (firsttime) then
        if (Node.eq.0) then
          DirectPhi  = fdf_boolean('DirectPhi',.false.)
        endif
#ifdef MPI
        call MPI_Bcast(DirectPhi,1,MPI_logical,0,MPI_Comm_World,
     *    MPIerror)
#endif
        firsttime = .false.
      endif

! Allocate local scratch array
      call re_alloc( numphi, 1,nmpl, myName//'numphi' )
      numphi(1:nmpl) = 0

! Find number of atomic orbitals at mesh points
      do io = 1,norb
        ia = iaorb(io)
        if (ipa(ia)==0) cycle ! Atomic sphere does not intersec myBox
        iphi = iphorb(io)
        is = isa(ia)
        r2o = rcut(is,iphi)**2

        ! Loop over mesh points inside rmax
        do iop = 1,mop
          ip0 = indexp( ipa(ia) + idop(iop) )
          if (ip0 .gt. 0) then

            ! Loop over sub-points to find if point is within range
            within = .false.
            do isp = 1,nsp
              dxsp(1:3,isp) = xdop(1:3,iop)+xdsp(1:3,isp)-dxa(1:3,ia)
              r2sp(isp) = sum( dxsp(1:3,isp)**2 )
              if (r2sp(isp) .lt. r2o) within = .true.
            enddo

            ! If within range, add one to number of point orbitals
            if (within) numphi(ip0) = numphi(ip0) + 1

          endif
        enddo
      enddo

! Initialise pointer array
      endpht(0) = 0
      do ip = 1,nmpl
        endpht(ip) = endpht(ip-1) + numphi(ip)
      enddo

! Allocate phi if this is not a direct calculation
      nlist = endpht(nmpl)
      if (DirectPhi) then
        nphi = 1
      else
        nphi = nlist
      endif

! Add an extra margin of error to nlist to minimise reallocations
      nliste = 1.01*nlist

! Adjust dimensions of phi if necessary
      if (associated(phi)) then
        nmax = size(phi,2)
      else
        nmax = 0
      end if
      if (nphi.gt.nmax) then
        if (DirectPhi) then
          call re_alloc( phi, 1,nsp, 1,nphi,   myName//'phi' )
        else                                                                
          call re_alloc( phi, 1,nsp, 1,nliste, myName//'phi' )
        endif 
      endif                                                              

! Adjust dimensions of list arrays if necessary
      if (associated(lstpht) .and. associated(listp2)) then
        nmax = min( size(lstpht), size(listp2) )
      else
        nmax = 0
      end if
      if (nlist.gt.nmax) then
        call re_alloc( lstpht, 1,nliste, myName//'lstpht' )
        call re_alloc( listp2, 1,nliste, myName//'listp2' )
      endif

! Find indexes and values of atomic orbitals at mesh points
      numphi = 0
      do io = 1,norb
        ia = iaorb(io)
        if (ipa(ia)==0) cycle ! Atomic sphere does not intersec myBox
        iphi = iphorb(io)
        is = isa(ia)
        r2o = rcut(is,iphi)**2
        do iop = 1,mop
          ip0 = indexp( ipa(ia) + idop(iop) )
          if (ip0 .gt. 0) then
            within = .false.
            do isp = 1,nsp
              dxsp(1:3,isp) = xdop(1:3,iop)+xdsp(1:3,isp)-dxa(1:3,ia)
              r2sp(isp) = sum( dxsp(1:3,isp)**2 )
              if (r2sp(isp) .lt. r2o) within = .true.
            enddo
            if (within) then
              numphi(ip0) = numphi(ip0) + 1
              n = endpht(ip0-1) + numphi(ip0)
              lstpht(n) = io
              listp2(n) = iop
              if (.not.DirectPhi) then
                do isp = 1,nsp
                  if (r2sp(isp) .lt. r2o) then
                    call phiatm( is, iphi, dxsp(1,isp), phip, grphi )
                    phi(isp,n) = phip
                  else
                    phi(isp,n) = 0.0_grid_p
                  endif
                enddo
              endif
            endif
          endif
        enddo
      enddo

! Free local memory
      call de_alloc( numphi, myName//'numphi' )

      end subroutine PhiOnMesh

