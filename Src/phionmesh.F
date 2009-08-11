      subroutine PhiOnMesh( nmpl, norb, iaorb, iphorb, isa, numphi )
C
C Calculates the values of the orbitals at the mesh points
C
C ----------------------------------------------------------------------
C Input :
C ----------------------------------------------------------------------
C integer nmpl          : Number of mesh points in unit cell locally
C integer norb          : Total number of basis orbitals in supercell
C integer iaorb(norb)   : Atom to which each orbital belongs
C integer iphorb(norb)  : Orbital index (within atom)
C integer isa(na)       : Species index of all atoms in supercell
C ----------------------------------------------------------------------
C Output :
C ----------------------------------------------------------------------
C All output quantities are in the module meshphi
C ----------------------------------------------------------------------

C Modules
      use precision, only : dp, grid_p
      use atmfuncs,  only : rcut, phiatm
      use mesh,      only : dxa, idop, indexp, ipa, mop, nsp, xdop, xdsp
      use meshphi,   only : directphi, nphi, phi, lstpht, listp2, endpht
      use parallel,  only : Node
      use fdf,       only : fdf_boolean
      use alloc,     only : re_alloc
#ifdef MPI
      use mpi_siesta
#endif
      implicit none
C Passed arguments
      integer,    intent(in) :: nmpl, norb, iaorb(norb), iphorb(norb),
     &                          isa(*)
      integer, intent(inout) :: numphi(nmpl)
C     Local variables
      integer                :: ia, io, iop, ip, ip0, iphi, is, isp, n,
     &                          nlist, nliste

#ifdef MPI
      integer                :: MPIerror
#endif
      logical                :: within
      logical,          save :: firsttime = .true.
      real(dp)               :: dxsp(3,nsp), grphi(3), phip, r2o,
     &                          r2sp(nsp)

!------------------------------------------------------------------------- BEGIN
#if defined(MPI_TRACE) && defined(MPI)
      call mpitrace_eventandcounters( 1000, 15 )
#endif
C     On first call set the logical DirectPhi
      if (firsttime) then
        if (Node.eq.0) then
          DirectPhi  = fdf_boolean( 'DirectPhi', .false. )
        endif
#ifdef MPI
        call MPI_Bcast( DirectPhi, 1, MPI_logical, 0, MPI_Comm_World,
     &                  MPIerror )
#endif
        firsttime = .false.
      else
        call resetMeshPhi( )
      endif

C     Allocate memory related to nmpl
      nullify( endpht, phi, lstpht, listp2 )
      call re_alloc( endpht, 0, nmpl, 'endpht', 'PhiOnMesh' )

C     Initialise pointer array
      endpht(0) = 0
      do ip = 1,nmpl
        endpht(ip) = endpht(ip-1) + numphi(ip)
      enddo

C     Allocate phi if this is not a direct calculation
      nlist = endpht(nmpl)
      if (DirectPhi) then
        nphi = 1
      else
        nphi = nlist
      endif

C     Add an extra margin of error to nlist to minimise reallocations
      nliste = 1.01*nlist

C     Adjust dimensions of phi if necessary
      if (DirectPhi) then
        call re_alloc( phi, 1, nsp, 1, nphi, 'phi', 'PhiOnMesh' )
      else
        call re_alloc( phi, 1, nsp, 1, nliste, 'phi', 'PhiOnMesh' )
      endif

C     Adjust dimensions of list arrays if necessary
      call re_alloc( lstpht, 1, nliste, 'lstpht', 'PhiOnMesh' )
      call re_alloc( listp2, 1, nliste, 'listp2', 'PhiOnMesh' )

C     Find indexes and values of atomic orbitals at mesh points
      numphi = 0
      do io = 1,norb
        ia = iaorb(io)
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
                    call phiatm( is, iphi, dxsp(1,isp),
     &                           phip, grphi )
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
#if defined(MPI_TRACE) && defined(MPI)
      call mpitrace_eventandcounters( 1000, 0 )
#endif
!--------------------------------------------------------------------------- END
      end subroutine PhiOnMesh
