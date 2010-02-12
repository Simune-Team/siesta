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


C ==================================================================
C Computes the number of orbitals that intersects every mesh point
C and calls initMeshDistr in order to compute a new data distribution
C based in the number of orbital intersections.
C ==================================================================
C SUBROUTINE distriPhiOnMesh( nm, nmpl, norb,
C     &                       iaorb, iphorb, isa )
C
C INPUT:
C integer nm(3)         : Number of Mesh divisions of each cell vector
C integer nmpl          : Local number of mesh points
C integer norb          : Total number of basis orbitals in supercell
C integer iaorb(norb)   : Atom to which each orbital belongs
C integer iphorb(norb)  : Orbital index (within atom) of each orbital
C integer isa(na)       : Species index of all atoms in supercell
C
C OUTPUT:
C The output values are in the module moreMeshSubs (New limits of the mesh
C and the information needed to transfer data between data
C distributions)
C
C BEHAVIOR:
C The amount of computation in every point of the mesh depends on the
C number of orbitals that intersecs in that mesh point. This can lead
C to unbalanced distributions in parallel executions. We want to
C compute a new data distribution, where the load is balanced among
C the several processes.
C
C In order to compute the number of orbitals that intersects in every
C mesh point, we loop in all orbitals of every atom and we increase the
C value of numphi in every point of the mesh that is inside of the
C radius of the current orbital.
C
C Finally, we call initMeshDistr that computes the new data distributions
C that will be used inside setup_hamiltonian. We are going to have 3
C different data distributions:
C  - Uniform:   Same amount of mesh in every process. Used in Poison
C  - Quadratic: The load in every point of the mesh is proporcional
C               to the square of the number of orbitals that intersects
C               in that point. Used in rhoofd and vmat.
C  - Lineal:    We try to distribute those points of the mesh that have
C               at least one orbital intersection. Used in cellxc.
C
C ==================================================================
      subroutine distriPhiOnMesh( nm, nmpl, norb, iaorb,
     &                            iphorb, isa, numphi )

      use mesh,         only: mop, ipa, idop, nsp, xdop, xdsp, dxa,
     &                        indexp, meshLim
      use precision,    only: dp, grid_p
      use atmfuncs,     only: rcut
      use parallel,     only: Nodes, node
#ifdef MPI
      use moreMeshSubs, only: initMeshExtencil
#endif
      use moreMeshSubs, only: initMeshDistr
      use cellxc_mod,   only: GGA, setGGA
      use alloc,        only: re_alloc, de_alloc
      implicit none
C     Passed arguments
      integer, intent(in)  :: nm(3), nmpl, norb, iaorb(norb),
     &                        iphorb(norb), isa(*)
      integer, intent(out) :: numphi(nmpl)
C     Local variables
      integer              :: ii, io, ia, iphi, is, iop, ip0, isp
      real(dp)             :: r2o, dxsp(3,nsp), r2sp(nsp)
      logical              :: within
      integer,     pointer :: wload(:)
!------------------------------------------------------------------------- BEGIN
#ifdef DEBUG
      call write_debug( '      PRE distriPhiOnMesh' )
#endif
      call timer( 'REMESH', 1 )

      do ii= 1, nmpl
        numphi(ii) = 0
      enddo

C     Find number of atomic orbitals at mesh points
      do io = 1,norb
        ia   = iaorb(io)
        if (ipa(ia)==0) cycle ! Atomic sphere does not intersec my Box
        iphi = iphorb(io)
        is   = isa(ia)
        r2o  = rcut(is,iphi)**2
C       Loop over mesh points inside rmax
        do iop = 1, mop
          ip0 = indexp( ipa(ia) + idop(iop) )
          if (ip0 .gt. 0) then
C           Loop over sub-points to find if point is within range
            within = .false.
            do isp = 1,nsp
              dxsp(1:3,isp) = xdop(1:3,iop)+xdsp(1:3,isp)-dxa(1:3,ia)
              r2sp(isp)     = sum( dxsp(1:3,isp)**2 )
              if (r2sp(isp) .lt. r2o) within = .true.
            enddo
C           If within range, add one to number of point orbitals
            if (within) numphi(ip0) = numphi(ip0) + 1
          endif
        enddo
      enddo

#ifdef MPI
      if (nodes.gt.1) then
        nullify( wload )
        call re_alloc( wload, 1, nmpl, 'wload', 'distriPhiOnMesh' )
        do ip0= 1, nmpl
          wload(ip0) = numphi(ip0)*numphi(ip0)
        enddo
        call initMeshDistr( 1, 2, nm, wload )

        do ip0= 1, nmpl
          if (numphi(ip0).ne.0) wload(ip0) = 1
        enddo
        call initMeshDistr( 1, 3, nm, wload )
        if (GGA) call initMeshExtencil( 3, nm )

C       Free local memory
        call de_alloc( wload, 'wload', 'distriPhiOnMesh' )

      endif
      call timer( 'REMESH', 2 )
#endif
#ifdef DEBUG
      call write_debug( '      POS distriPhiOnMesh' )
#endif
!--------------------------------------------------------------------------- END
      end subroutine distriPhiOnMesh
