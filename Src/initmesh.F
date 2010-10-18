      subroutine InitMesh( na, cell, norb, iaorb, iphorb, isa,
     &                     rmax, G2max, G2mesh, nsc, nmpl, 
     &                     nm, nml, ntm, ntml, ntpl, dvol)

!  Modules
      use precision,     only : dp
      use parallel,      only : Node, Nodes
      use parallelsubs,  only : HowManyMeshPerNode, GlobalToLocalMesh
      use moreMeshSubs,  only : initMeshDistr, setMeshDistr
      use moreMeshSubs,  only : UNIFORM
      use alloc,         only : re_alloc, de_alloc
      use mesh,          only : ne, mop, nmeshg, nmsc, nsm, nsp
      use mesh,          only : cmesh, rcmesh, xdsp, nmuc, nusc
      use siesta_cml
      use cellsubs,      only : reclat  ! Finds reciprocal unit cell
      use cellsubs,      only : volcel  ! Finds unit cell volume
      use m_chkgmx,      only : chkgmx  ! Checks planewave cutoff of a mesh
      use fft1d,         only : nfft    ! Finds allowed value for 1-D FFT
#ifdef MPI
      use mpi_siesta
#endif

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
!!OLD      integer, intent(out):: ntopl       ! Total number of nonzero
                                         ! orbital points stored by my
                                         ! processor
      real(dp),intent(out):: dvol        ! Volume per mesh (super)point

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
C integer j             : General-porpose index
C integer j1,j2,j3      : Mesh indexes in each mesh direction
C real*8  k0(3)         : Zero-vector argument for routine chkgmx
C integer nep           : Number of extended-mesh points
C integer nmp           : Number of mesh points in unit cell
C real*8  pldist        : Distance between mesh planes
C real*8  r             : Distance between atom and mesh point
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

C     Local variables
      integer                 :: i, ia, i1, i2, i3, io, iphi, is, isp,
     &                           ity, j, nep, j1, j2, j3, indi
      integer*8               :: ntp
      real(dp)                :: dx(3), dxp(3), pldist, r,
     &                           rcell(3,3), vecmod, volume
      real(dp), parameter     :: k0(3) = (/ 0.0, 0.0, 0.0 /)
      logical                 :: within
C     Functions
      real(dp)                :: dismin
      external                :: dismin
!--------------------------------------------------------------------- BEGIN
#ifdef DEBUG
      call write_debug( '      PRE InitMesh' )
#endif

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
!     nmeshg in module array
      nmeshg(:) = ntm(:) ! Total number of mesh points, incl. subpoints
      nmuc(:) = nm(:)    ! Mesh points in each unit cell vector
      nusc(:) = nsc(:)   ! Unit cells in each supercell direction

! Find number of mesh points in unit cell.
! Note needed integer*8 !
!!AG** Check proper kind name

      ntp = product(INT(ntm,kind=8))

C     Create the first mesh distribution
      call initMeshDistr( oDistr=UNIFORM, nm=nm )

C     Find and sets local number of Mesh points of each kind
!     Sets nml, nmpl, ntml and ntpl
      call setMeshDistr( UNIFORM, nsm, nsp, nml, nmpl, ntml, ntpl )

C     Find volume of unit cell and of mesh cell
      volume = volcel( cell )
      dvol = volume / ntp

C     Output current mesh dimensions and cut-off
      if (Node.eq.0) then
        write(6,'(/,a,3(i6,a),i12)') 'InitMesh: MESH =',
     &        ntm(1),' x',ntm(2),' x',ntm(3),' =', ntp
        write(6,'(a,3(i6,a),i12)') 'InitMesh: (bp) =',
     &        nm(1),' x',nm(2),' x',nm(3),' =', product(INT(nm,8))
        write(6,'(a,2f10.3,a)')
     &        'InitMesh: Mesh cutoff (required, used) =',
     &        G2max, G2mesh, ' Ry'
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

C     Find mesh-cell vectors
      do i = 1,3
        do j = 1,3
          cmesh(j,i) = cell(j,i) / nm(i)
        enddo
      enddo

C     Find reciprocal mesh-cell vectors (not multiplied by 2*pi)
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

! Find number of points within rmax (orbital points)
      mop = 0
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
            endif ! (within)
          enddo ! i1
        enddo ! i2
      enddo ! i3

C     Create extended mesh arrays for the first data distribution
      call setupExtMesh( UNIFORM, rmax )

#ifdef DEBUG
      call write_debug( '      POS InitMesh' )
#endif
!----------------------------------------------------------------------- END
      end subroutine InitMesh

