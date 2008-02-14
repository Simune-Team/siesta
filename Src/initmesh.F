      subroutine InitMesh( na, cell, norb, iaorb, iphorb, isa,
     &                     rmax, G2max, G2mesh, nsc, nmpl, 
     &                     nm, nml, ntm, ntml, ntpl, ntopl, dvol)
C
C Initialises the mesh
C
C ----------------------------------------------------------------------
C Input :
C ----------------------------------------------------------------------
C integer na            : Number of atoms in supercell
C real*8  cell(3,3)     : Auxillary lattice vectors
C integer norb          : Total number of basis orbitals in supercell
C integer iaorb(norb)   : Atom to which each orbital belongs
C integer iphorb(norb)  : Orbital index (within atom) of each orbital
C integer isa(na)       : Species index of all atoms in supercell
C real*8  rmax          : Maximum orbital radius
C integer nsc(3)        : Number of unit-cells in each supercell direct.
C ----------------------------------------------------------------------
C Input and output :
C ----------------------------------------------------------------------
C real*8  G2max         : Effective planewave cutoff (Ry) determines
C                       : mesh density and the precision of integrals
C                       : On input : Value required
C                       : On output: Value used, which may be larger
C ----------------------------------------------------------------------
C Output :
C ----------------------------------------------------------------------
C real*8  G2mesh        : Effective planewave cutoff of mesh used
C integer nmpl          : Number of mesh points in unit cell locally
C integer nm(3)         : Number of Mesh divisions of each cell vector
C integer nml(3)        : Local form of nm
C integer ntm(3)        : Total number of mesh points
C integer ntml(3)       : Total number of mesh points stored locally
C integer ntpl          : Number of mesh Total Points in unit cell
C                       : (including subpoints) locally
C integer ntopl         : Total number of nonzero orbital points locally
C real*8  dvol          : Mesh-cell volume
C ----------------------------------------------------------------------
C Internal variables and arrays:
C ----------------------------------------------------------------------
C real*8  dx(3)         : Vector from atom to mesh sub-point
C real*8  dxp(3)        : Vector from atom to mesh point
C integer i             : General-purpose index
C integer ia            : Looping variable for number of atoms
C integer i1,i2,i3      : Mesh indexes in each mesh direction
C integer is            : Species index
C integer isp           : Sub-Point index
C integer ity           : Orbital-type index
C integer j             : General-porpose index
C integer j1,j2,j3      : Mesh indexes in each mesh direction
C real*8  k0(3)         : Zero-vector argument for routine chkgmx
C integer nep           : Number of extended-mesh points
C integer nmp           : Number of mesh points in unit cell
C integer nty           : Number of orbital types
C integer noty(nty)     : Number of orbitals of each type
C real*8  pldist        : Distance between mesh planes
C real*8  r             : Distance between atom and mesh point
C real*8  rcty(nty)     : Radius of each orbital type
C logical same          : Indicates whether orbital types are the same
C real*8  vecmod        : Vector modulus
C real*8  volume        : Unit cell volume
C logical within        : Is a mesh point within orbital range?
C ----------------------------------------------------------------------
C Units :
C ----------------------------------------------------------------------
C
C Energies in Rydbergs
C Distances in Bohr
C

C  Modules
      use precision,     only : dp
      use parallel,      only : Node, Nodes
      use parallelsubs,  only : HowManyMeshPerNode, GlobalToLocalMesh
      use atmfuncs,      only : rcut
      use moreMeshSubs,  only : initMeshDistr, setMeshDistr,
     &                          resetMeshDistr
      use mesh
      use siesta_cml
#ifdef MPI
      use mpi_siesta
#endif
      implicit none
C     Passed arguments
      integer     :: na, ntpl, ntopl, nm(3), nml(3), nmpl,
     &               norb, nsc(3), ntm(3), ntml(3), iaorb(norb),
     &               iphorb(norb), isa(na)
      real(dp)    :: cell(3,3), dvol, G2max, G2mesh, rmax
C     Local variables
      integer                    :: i, ia, i1, i2, i3, indi, io, iphi,
     &                              is, isp, ity, j, j1, j2, j3,
     &                              ncells, nep, ntop, ntp, nty, odistr
      integer,  allocatable, save :: iphty(:), isty(:), noty(:)
      real(dp)                    :: dx(3), dxp(3), pldist, r,
     &                               rcell(3,3), vecmod, volume
      real(dp),              save :: k0(3)
      real(dp), allocatable, save :: rcty(:)

      logical                    :: same, within
C     Functions
      real(dp)                   :: dismin, volcel

      external                      chkgmx, dismin, nfft, reclat,
     &                              timer, volcel, memory
      data                          k0 / 3*0.0_dp /
#ifdef MPI
      integer             :: MPIerror
#endif

!------------------------------------------------------------------------- BEGIN
C ----------------------------------------------------------------------
C     Orbital type initialization 
C ----------------------------------------------------------------------
C     Find number of orbital types
      allocate(isty(norb))
      call memory('A','I',norb,'InitMesh')
      allocate(iphty(norb))
      call memory('A','I',norb,'InitMesh')
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
      call memory('D','I',size(isty),'InitMesh')
      deallocate(isty)
      call memory('D','I',size(iphty),'InitMesh')
      deallocate(iphty)

C     Allocate local arrays dependent on nty
      allocate(iphty(nty))
      call memory('A','I',nty,'InitMesh')
      allocate(isty(nty))
      call memory('A','I',nty,'InitMesh')
      allocate(noty(nty))
      call memory('A','I',nty,'InitMesh')
      allocate(rcty(nty))
      call memory('A','D',nty,'InitMesh')

C     Find number of orbital types, number of orbitals of each type,
C     and cutoff radius of each type.
      noty(1:nty) = 0
      nty = 0
      do io = 1,norb
        ia = iaorb(io)
        is = isa(ia)
        iphi = iphorb(io)
        do ity = 1,nty
          if (is.eq.isty(ity) .and. iphty(ity).eq.iphi) then
            noty(ity) = noty(ity) + 1
            goto 160
          endif
        enddo
        nty = nty + 1
        isty(nty) = is
        iphty(nty) = iphi
        noty(nty) = 1
        rcty(nty) = rcut( is, iphi )
  160   continue
      enddo

      if (allocated(iphty)) then
        call memory('D','I',size(iphty),'InitMesh')
        deallocate(iphty)
      endif
      if (allocated(isty)) then
        call memory('D','I',size(isty),'InitMesh')
        deallocate(isty)
      endif

C ----------------------------------------------------------------------
C     Mesh initialization 
C ----------------------------------------------------------------------
C     Find reciprocal cell vectors (multiplied by 2*pi)
      call reclat( cell, rcell, 1 )

C     Find number of mesh intervals for each cell vector.
C     Loop over cell vectors
      do i = 1,3
C       The reciprocal vectors of the mesh unit cell (cell/ntm)
C       are rcell*ntm, and must be larger than 2*G2max
        vecmod = sqrt(dot_product(rcell(:,i),rcell(:,i)))
        ntm(i) = 2 * sqrt(G2max) / vecmod + 1
      enddo

C     Return here until mesh cut-off is large enough
   50 continue
      do i = 1,3
C       NFFT selects appropriate number of points for fft
   55   call nfft( ntm(i) )
C       Require that ntm(i) to be a multiple of nsm
        if ( mod( ntm(i), nsm ) .ne. 0 ) then
          ntm(i) = ntm(i) + 1
          goto 55
        endif
        nm(i) = ntm(i) / nsm
        nmsc(i) = nm(i) * nsc(i)
      enddo

C     Check that effective cut-off is large enough as for non-right angled
C     unit cells this is not guaranteed to be the case. If cut-off needs
C     to be larger, increase ntm and try again.
      G2mesh = 1.0d6
      call chkgmx( k0, rcell, ntm, G2mesh )
      if (G2mesh .lt. G2max) then
        ntm(1:3) = ntm(1:3) + 1
        goto 50
      endif

C     Store number of mesh points in module array
      nmeshg(1:3) = ntm(1:3)

C     Find number of mesh points in unit cell.
      ntp = ntm(1) * ntm(2) * ntm(3)

C     Reset previous mesh distributions
      call resetMeshDistr( )

C     Create the first mesh distribution
      call initMeshDistr( nm=nm, oDistr=oDistr )

C     Find and sets local number of Mesh points of each kind
      call setMeshDistr( oDistr, nsm, nsp, nml, nmpl, ntml, ntpl )

C     Find volume of unit cell and of mesh cell
      volume = volcel( cell )
      dvol = volume / ntp

C     Output current mesh dimensions and cut-off
      if (Node.eq.0) then
        write(6,'(/,a,3(i6,a),i12)') 'InitMesh: MESH =',
     &        ntm(1),' x',ntm(2),' x',ntm(3),' =', ntp
        write(6,'(a,2f10.3,a)')
     &        'InitMesh: Mesh cutoff (required, used) =',
     &        G2max, G2mesh, ' Ry'
      endif
      if (cml_p) then
        call cmlStartPropertyList(mainXML)
        call cmlAddProperty( xf=mainXML, property=ntm,
     &                       dictref='siesta:ntm', title='Mesh')
        call cmlAddProperty( xf=mainXML, property=G2max,
     &                       units='siestaUnits:Ry',
     &                       dictref='siesta:g2max',
     &                       title='Requested Cut-Off')
        call cmlAddProperty( xf=mainXML, property=G2mesh,
     &                       units='siestaUnits:Ry',
     &                       dictref='siesta:g2mesh',
     &                       title='Actual Cut-Off' )
        call cmlEndPropertyList(mainXML)
      endif
      G2max = G2mesh

C     Find mesh-cell vectors
      do i = 1,3
        do j = 1,3
          cmesh(j,i) = cell(j,i) / nm(i)
        enddo
      enddo

C     Find reciprocal mesh-cell vectors (not multiplied by 2*pi)
      call reclat( cmesh, rcmesh, 0 )

C     Find number of extended-mesh intervals for each cell vector.
C     Loop over mesh directions
      do i = 1,3
C       pldist is the distance between mesh planes
        pldist = 1.0_dp / sqrt(dot_product(rcmesh(:,i),rcmesh(:,i)))
C       Find number of planes spanned by rmax
        ne(i) = rmax / pldist + real(nsm-1,kind=dp) / real(nsm,kind=dp)
        ne(i) = max( ne(i), nsm-2 )
C       Add ne(i) points to the left and ne(I)+1 points to the
C       right, to cover the spilling rmax from an atom at any
C       possible place within the unit cell.
        nem(i) = nmsc(i) + 2 * ne(i) + 1
      enddo

C     Find total number of extended-mesh points.
      nep = nem(1) * nem(2) * nem(3)

C     Find sub-points
      isp = 0
      do i3 = 0, nsm-1
        do i2 = 0, nsm-1
          do i1 = 0, nsm-1
            isp = isp + 1
            do i = 1,3
              xdsp(i,isp) = ( cmesh(i,1) * i1 +
     &                        cmesh(i,2) * i2 +
     &                        cmesh(i,3) * i3 ) / nsm
            enddo
          enddo
        enddo
      enddo

C     Find number of points only within rmax (orbital points)
      mop = 0
      ntop = 0

C     Loop over possible points within rmax
      do i3 = -ne(3), ne(3)+1
        do i2 = -ne(2), ne(2)+1
          do i1 = -ne(1), ne(1)+1
C          Find point coordinates
            do i = 1,3
              dxp(i) = cmesh(i,1) * i1 +
     &                 cmesh(i,2) * i2 +
     &                 cmesh(i,3) * i3 
            enddo

C           Loop over sub-points
            within = .false.
            do isp = 1,nsp
C             Find point coordinates
              dx(1:3) = dxp(1:3) + xdsp(1:3,isp) 
C             Find distance from point to mesh cell
              r = dismin( cmesh, dx )
              if ( r .lt. rmax ) within = .true.
            enddo
            if ( within ) then
C             mop is the number of mesh points within rmax
              mop = mop + 1
              do ity = 1,nty
                if ( r .lt. rcty(ity) ) ntop = ntop + noty(ity)
              enddo
            endif
          enddo
        enddo
      enddo

C     Calculate ntopl
      ncells = nsc(1) * nsc(2) * nsc(3)
      ntop = ntop / ncells
      ntopl = ((ntop/Nodes) + 1)

C     Allocate local memory that depends on mop
      if (allocated(idop)) then
        call memory('D','I',size(idop),'InitMesh')
        deallocate(idop)
      endif
      if (allocated(xdop)) then
        call memory('D','D',size(xdop),'InitMesh')
        deallocate(xdop)
      endif
      allocate(idop(mop))
      call memory('A','I',mop,'InitMesh')
      allocate(xdop(3,mop))
      call memory('A','D',3*mop,'InitMesh')

C     Find points within rmax (orbital points)
      mop = 0
C     Loop over possible points within rmax
      do i3 = -ne(3), ne(3)+1
        do i2 = -ne(2), ne(2)+1
          do i1 = -ne(1), ne(1)+1
C           Find point coordinates
            do i = 1,3
              dxp(i) = cmesh(i,1) * i1 +
     &                 cmesh(i,2) * i2 +
     &                 cmesh(i,3) * i3 
            enddo
C           Loop over sub-points
            within = .false.
            do isp = 1,nsp
C             Find point coordinates
              dx(1:3) = dxp(1:3) + xdsp(1:3,isp) 
C             Find distance from point to mesh cell
              r = dismin( cmesh, dx )
              if ( r .lt. rmax ) within = .true.
            enddo
            if ( within ) then
              mop = mop + 1
C             Store index-distance and vector-distance to point.
              idop(mop) = i1 + nem(1) * i2 + nem(1) * nem(2) * i3
              xdop(1:3,mop) = dxp(1:3)
            endif
          enddo
        enddo
      enddo

C     Deallocate orbital related memory 
      if (allocated(noty)) then
        call memory('D','I',size(noty),'InitMesh')
        deallocate(noty)
      endif
      if (allocated(rcty)) then
        call memory('D','D',size(rcty),'InitMesh')
        deallocate(rcty)
      endif

!--------------------------------------------------------------------------- END
      end subroutine InitMesh
