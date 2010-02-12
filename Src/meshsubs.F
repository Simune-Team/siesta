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
C Includes the following subroutines connected to the mesh :
C
C InitMesh          = initialises the mesh
C SameMeshAndAtoms  = checks whether the mesh or atoms have changed
C InitAtomMesh      = initialises quantities relating to the atoms
C                     and the mesh
C PartialCoreOnMesh = calculates the partial core density on mesh
C NeutralAtomOnMesh = calculates the neutral atom potential on mesh
C PhiOnMesh         = calculates the orbital values on the mesh
C
      module mesh
C
C Stores quantities that are connected with the mesh
C
      use precision, only : dp
      implicit none
C ----------------------------------------------------------------------
C Mesh variables
C ----------------------------------------------------------------------
C real*8  cmesh(3,3)    : Mesh-cell vectors
C real*8  dxa(3,na)     : Atom position within mesh-cell
C integer idop(mop)     : Extended-mesh-index displacement of points
C                       : within a sphere of radius rmax
C integer indexp(nep)   : Translation from extended to normal mesh index
C integer ipa(na)       : Mesh cell in which atom is
C integer mop           : Maximum number of non-zero orbital points
c integer nmeshg(3)     : Total number of mesh points in each direction
C integer ne(3)         : Number of mesh-Extension intervals 
C                       : in each direction
C integer nem(3)        : Extended-mesh divisions in each direction
C integer nmsc(3)       : Mesh divisions of each supercell vector
C integer meshLim(2,3)  : Upper an lower limits of the mesh in every process
C integer nsm           : Number of mesh sub-divisions in each direction
C integer nsp           : Number of sub-points of each mesh point
C real*8  rcmesh(3,3)   : Reciprocal mesh-cell vectors
C                       : (WITHOUT 2*pi factor)
C real*8  xdop(3,mop)   : Vector to mesh points within rmax
C real*8  xdsp(3,nsp)   : Vector to mesh sub-points
C ----------------------------------------------------------------------
      integer,  pointer, save :: idop(:)
      integer,  pointer, save :: ipa(:) => null()
      real(dp), pointer, save :: dxa(:,:) => null()
      real(dp), pointer, save :: xdop(:,:)
      real(dp), pointer, save :: xdsp(:,:) ! Allocated in dhscf
      integer,           save :: mop
      integer,           save :: ne(3)
      integer,           save :: nem(3)
      integer,           save :: nmsc(3)
      integer,           save :: meshLim(2,3)  ! My processor's box of mesh points:
                                               !   myBox(1,:) : lower bounds
                                               !   myBox(2,:) : upper bounds
      integer,           save :: nmeshg(3)
      integer,           save :: nsm
      integer,           save :: nsp
      real(dp),          save :: cmesh(3,3)
      real(dp),          save :: rcmesh(3,3)
      integer,  pointer, save :: indexp(:) => null()
      end module mesh

      module meshphi
      use precision, only : grid_p
C
C Stores quantities that are connected with phi on the mesh
C
      implicit none

C ----------------------------------------------------------------------
C Phi-mesh variables
C ----------------------------------------------------------------------
C logical DirectPhi     : If true the phi is calculated on the fly
C integer endpht(0:nmpl): Last position occupied by points in lstpht
C integer lstpht(ntopl) : List of non-zero orbitals at point
C integer listp2(ntopl) : Maps orbital-mesh point to iop
C integer nphi          : Length of phi array second dimension
C real    phi(nsp,ntopl): Basis orbitals at mesh points (sparse)
C ----------------------------------------------------------------------
      logical,                   save :: DirectPhi
      integer,                   save :: nphi
      integer,          pointer, save :: endpht(:)
      integer,          pointer, save :: lstpht(:)
      integer,          pointer, save :: listp2(:)
      real(grid_p),     pointer, save :: phi(:,:)
      end module meshphi

      subroutine resetMeshPhi( )
C
C Deallocate arrays of module meshphi if necessary
C
      use meshphi
      use alloc,         only : de_alloc
      implicit none
      if (associated(endpht))
     &  call de_alloc( endpht, 'endpht',  'PhiOnMesh' )

      if (associated(lstpht))
     &  call de_alloc( lstpht, 'lstpht',  'PhiOnMesh' )


      if (associated(listp2))
     &  call de_alloc( listp2, 'listp2',  'PhiOnMesh' )

      if (associated(phi))
     &  call de_alloc( phi, 'phi',  'PhiOnMesh' )

      end subroutine resetMeshPhi

      subroutine SameMeshAndAtoms(na, xa, ucell, rmax, G2max, G2mesh,
     &                            samesh, samexa)
C
C Checks whether anything has changed that requires the mesh to be
C reinitialised or quantities relating to the atoms relative to
C the mesh.
C
C ----------------------------------------------------------------------
C Input :
C ----------------------------------------------------------------------
C integer na           : Number of atoms in the supercell
C real*8  xa(3,na)     : Coordinates of the atoms in the supercell
C real*8  ucell(3,3)   : Current unit cell vectors
C real*8  rmax         : Maximum orbital radius
C real*8  G2max        : Requested mesh cut-off
C real*8  G2mesh       : Mesh cut-off from last call
C ----------------------------------------------------------------------
C Output :
C ----------------------------------------------------------------------
C logical samesh       : If .true. then the mesh must be reinitialised
C logical samexa       : If .true. then atom related quantities need
C                      : to be recalculated
C ----------------------------------------------------------------------
C Internal :
C ----------------------------------------------------------------------
C integer lastna       : Number of atoms from last call
C real*8  lastxa(3,na) : Position of atoms from last call
C real*8  lastra       : Rmax from last call
C real*8  lastc(3,3)   : Unit cell from last call
C ----------------------------------------------------------------------
      use precision,     only : dp
      use alloc,         only : re_alloc, de_alloc
      implicit none
C     Passed arguments
      integer,    intent(in)  :: na
      real(dp),   intent(in)  :: G2max, G2mesh
      real(dp),   intent(in)  :: rmax
      real(dp),   intent(in)  :: ucell(3,3)
      real(dp),   intent(in)  :: xa(3,na)
      logical,    intent(out) :: samesh, samexa
C     Local variables
C     Saved
      integer,           save :: lastna     = 0
      real(dp),          save :: tiny       = 1.0e-12_dp,
     &                           lastc(3,3) = 0.743978657912656e50_dp,
     &                           lastra     = 0.0_dp
C     non-saved
      integer                 :: i, ia, j
      real(dp), pointer, save :: lastxa(:,:)

!------------------------------------------------------------------------- BEGIN
C     If number of atoms has changed then deallocate lastxa
      if (na.ne.lastna) then
        call re_alloc(lastxa, 1, 3, 1, na, 'lastxa', 'SameMeshAndAtoms')
        lastxa(1:3,1:na) = 0.0_dp
        lastxa(1,1) = 0.987654321273567e50_dp
      endif

C     Find if mesh has to be changed due to unit cell
      samesh = .true.
      do i = 1,3
        do j = 1,3
          if ( ucell(j,i) .ne. lastc(j,i) ) samesh = .false.
          lastc(j,i) = ucell(j,i)
        enddo
      enddo

C     Find if mesh has to be changed due to unit cell
      if ( G2max .gt. G2mesh * (1.0_dp + tiny) ) samesh = .false.

C     Find if mesh has to be changed due to rmax
      if (rmax .ne. lastra) samesh = .false.
      lastra = rmax

C     Find if atoms have moved having checked the number of atoms first
      samexa = (na.eq.lastna)
      if (samexa) then
        do ia = 1,na
          do i = 1,3
            if ( xa(i,ia) .ne. lastxa(i,ia) ) samexa = .false.
          enddo
        enddo
      endif

C     Copy the number of atoms and coordinates to save for next call
      lastna = na
      lastxa(1:3,1:na) = xa(1:3,1:na)

C     If cell has changed then it is necessary to reinitialise coordinates
      if (.not.samesh) samexa = .false.

!--------------------------------------------------------------------------- END
      end subroutine SameMeshAndAtoms

      subroutine PartialCoreOnMesh( na, isa, ntpl, rhopcc, indxua, 
     &                              nsd, dvol, volume, Vscf, Vaux, 
     &                              fal, stressl, Forces, Stress )
C
C Finds the partial-core-correction energy density on the mesh
C
C ----------------------------------------------------------------------
C Input :
C ----------------------------------------------------------------------
C integer na            : Number of atoms in supercell
C integer isa(na)       : Species index of all atoms in supercell
C integer ntpl          : Number of mesh Total Points in unit cell
C                       : (including subpoints) locally
C integer indxua        : Index of equivalent atom in unit cell
C integer nsd           : Number of diagonal spin values (1 or 2)
C real*8  dvol          : Mesh-cell volume
C real*8  volume        : Unit cell volume
C real*4  Vscf(ntpl)    : Hartree potential of SCF density
C real*4  Vaux(ntpl)    : Auxiliary potential array
C logical Forces        : Are the forces to be calculated?
C logical Stress        : Are the stresses to be calculated?
C ----------------------------------------------------------------------
C Output : (non-gradient case)
C ----------------------------------------------------------------------
C real*4  rhopcc(ntpl)  : Partial-core-correction density for xc
C ----------------------------------------------------------------------
C Output : (gradient case)
C ----------------------------------------------------------------------
C real*8  fal(3,:)     : Local copy of atomic forces
C real*8  stressl(3,3)  : Local copy of stress tensor
C ----------------------------------------------------------------------
      use precision, only: dp, grid_p
      use atmfuncs,  only: rcore, chcore_sub
      use mesh,      only: dxa, idop, indexp, ipa, mop, nsp, xdop, xdsp
      implicit none
C     Passed arguments
      integer,      intent(in)  :: na, isa(na), ntpl
      integer,      intent(in)  :: indxua(*)
      integer,      intent(in)  :: nsd
      real(dp),     intent(in)  :: dvol, volume
      real(grid_p), intent(in)  :: Vscf(ntpl,*)
      real(grid_p), intent(in)  :: Vaux(ntpl)
      logical,      intent(in)  :: Forces
      logical,      intent(in)  :: Stress
      real(grid_p), intent(out) :: rhopcc(ntpl)
      real(dp),   intent(inout) :: fal(3,*)
      real(dp), intent(inout)   :: stressl(3,3)
C     Internal variables
      integer                   :: i, ia, iop, ip, ip0, is, isp,
     &                             ispin, iua
      real(dp)                  :: dfa(3), dx(3), grrho(3), r, ra,
     &                             rhop, vxc
      real(dp),            save :: tiny
      logical                   :: Gradients
      data tiny / 1.0e-12_dp /

!------------------------------------------------------------------------- BEGIN
#ifdef DEBUG
      call write_debug( '      PRE PartialCoreOnMesh' )
#endif
C     Find out whether this is a gradient run based on the arguments passed
      Gradients = ( Forces .or. Stress )

C     Initialise array to zero
      if (.not.Gradients) rhopcc(1:ntpl) = 0.0_dp

      do ia = 1,na
        if (ipa(ia)==0) cycle ! Atomic sphere does not intersec my Box
        iua = indxua(ia)
        is = isa(ia)
        ra = rcore( is )
        if (ra .gt. tiny) then
C         Loop over mesh points inside rmax
          do iop = 1,mop
            ip0 = indexp( ipa(ia) + idop(iop) )
            if (ip0 .gt. 0) then
C             Loop over sub-points
              do isp = 1,nsp
                dx(1:3) = xdop(1:3,iop) + xdsp(1:3,isp) - dxa(1:3,ia)
                r = sqrt(dot_product(dx,dx))
C               if (r .lt. ra .and. r .gt. tiny) then
                if (r .lt. ra) then
                  ip = isp + nsp * (ip0 - 1)
                  call chcore_sub( is, dx, rhop, grrho )
                  if (Gradients) then
C                   Calculate gradients of PCC
                    do ispin = 1,nsd
                      vxc = Vscf(ip,ispin) - Vaux(ip)
                      do i = 1,3
                        dfa(i) = dvol * vxc * grrho(i) / nsd
                        if (Forces) fal(i,iua) = fal(i,iua) + dfa(i)
                        if (Stress) stressl(1:3,i) = stressl(1:3,i) + 
     &                    dx(1:3)*dfa(i)/volume
                      enddo
                    enddo
                  else
C                   Calculate density due to PCC
                    rhopcc(ip) = rhopcc(ip) + rhop
                  endif
                endif
              enddo
            endif
          enddo
        endif
      enddo
#ifdef DEBUG
      call write_debug( '      POS PartialCoreOnMesh' )
#endif
!--------------------------------------------------------------------------- END
      end subroutine PartialCoreOnMesh

      subroutine NeutralAtomOnMesh( na, isa, ntpl, vna,
     &                              indxua, dvol, volume,
     &                              drho, fal, stressl,
     &                              Forces, Stress )
C
C Finds the potential due to the neutral atoms on the mesh
C
C ----------------------------------------------------------------------
C Input :
C ----------------------------------------------------------------------
C integer na            : Number of atoms in supercell
C integer isa(na)       : Species index of all atoms in supercell
C integer ntpl          : Number of mesh Total Points in unit cell
C                       : (including subpoints) locally
C integer indxua        : Index of equivalent atom in unit cell
C real*8  dvol          : Mesh-cell volume
C real*8  volume        : Unit cell volume
C real    drho(ntpl)    : SCF density difference
C logical Forces        : Should the forces be calculated?
C logical Stress        : Should the stress be calculated?
C ----------------------------------------------------------------------
C Input/Output : (Output for non-gradient call, must be kept otherwise)
C ----------------------------------------------------------------------
C real    vna(ntpl)     : Sum of neutral-atom potentials
C ----------------------------------------------------------------------
C Output : (gradient call)
C ----------------------------------------------------------------------
C real*8  fal(3,:)     : Local copy of atomic forces
C real*8  stressl(3,3)  : Local copy of stress tensor
C ----------------------------------------------------------------------
      use precision, only: dp, grid_p
      use atmfuncs,  only: rcut, phiatm
      use mesh,      only: dxa, idop, indexp, ipa, mop, nsp, xdop, xdsp
      implicit none
C    Passed arguments
      integer,         intent(in) :: na, isa(na), ntpl
      integer,         intent(in) :: indxua(*)
      real(dp),        intent(in) :: dvol, volume
      real(grid_p),    intent(in) :: drho(*)
      logical,         intent(in) :: Forces, Stress
      real(grid_p), intent(inout) :: vna(ntpl)
      real(dp),     intent(inout) :: fal(3,*)
      real(dp),     intent(inout) :: stressl(3,3)
C     Internal variables
      integer                     :: i, ia, iop, ip, ip0, is, isp, iua
      real(dp)                    :: dx(3), grva(3), r, ra, va
      logical                     :: Gradients

#ifdef DEBUG
      call write_debug( '      PRE NeutralAtomOnMesh' )
#endif
C     Check whether forces and/or stress has been
C     requested based on arguments
      Gradients = (Forces.or.Stress)

C Initialise array to zero if not computing gradients
C If we are computing gradients, the current value of vna
C must be kept. This is particularly important for grid-cell sampling...

      if (.not.Gradients) vna(1:ntpl) = 0.0_dp

      do ia = 1,na
        if (ipa(ia)==0) cycle ! Atomic sphere does not intersec my Box
        iua = indxua(ia)
        is = isa(ia)
        ra = rcut( is, 0 )

C       Loop over mesh points inside rmax
        do iop = 1,mop
          ip0 = indexp( ipa(ia) + idop(iop) )
          if (ip0 .gt. 0) then
C           Loop over sub-points
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
     &                dx(1:3) * grva(i) / volume
                  enddo
                else
                  vna(ip) = vna(ip) + va
                endif
              endif
            enddo
          endif
        enddo
      enddo
#ifdef DEBUG
      call write_debug( '      POS NeutralAtomOnMesh' )
#endif
!--------------------------------------------------------------------------- END
      end subroutine NeutralAtomOnMesh

      subroutine LocalChargeOnMesh( na, isa, ntpl, Chlocal, indxua )
C
C Finds the diffuse ionic charge, whose electrostatic potential is equal
C to Vlocal on the mesh
C
C ----------------------------------------------------------------------
C Input :
C ----------------------------------------------------------------------
C integer na            : Number of atoms in supercell
C integer isa(na)       : Species index of all atoms in supercell
C integer ntpl          : Number of mesh Total Points in unit cell
C                       : (including subpoints) locally
C integer indxua        : Index of equivalent atom in unit cell
C ----------------------------------------------------------------------
C Output : 
C ----------------------------------------------------------------------
C real    Chlocal(ntpl)     : Sum of diffuse ionic charge densities
C ----------------------------------------------------------------------
      use precision, only: dp, grid_p
      use atmfuncs,  only: rcut, psch
      use atm_types, only: species
      use radial,    only: rad_func
      use mesh,      only: dxa, idop, indexp, ipa, mop, nsp, xdop, xdsp
      implicit none
C     Passed arguments
      integer, intent(in)       :: na, isa(na), ntpl
      integer, intent(in)       :: indxua(*)
      real(grid_p), intent(out) :: Chlocal(ntpl)
C     Internal variables
      integer                   :: ia, iop, ip, ip0, is, isp, iua
      real(dp)                  :: dx(3), grpsch(3), r, ra, pschloc
      type(rad_func), pointer  :: func   
 
!------------------------------------------------------------------------- BEGIN
#ifdef DEBUG
      call write_debug( '      PRE LocalChargeOnMesh' )
#endif
C     Initialise array to zero 
      Chlocal(1:ntpl) = 0.0_grid_p
 
      do ia = 1,na
        if (ipa(ia)==0) cycle ! Atomic sphere does not intersec my Box
        iua = indxua(ia)
        is = isa(ia)
        func => species(is)%chlocal
        ra = func%cutoff
C       Loop over mesh points inside rmax
        do iop = 1,mop
          ip0 = indexp( ipa(ia) + idop(iop) )
          if (ip0 .gt. 0) then
C           Loop over sub-points
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
#ifdef DEBUG
      call write_debug( '      POS LocalChargeOnMesh' )
#endif
!--------------------------------------------------------------------------- END
      end subroutine LocalChargeOnMesh
